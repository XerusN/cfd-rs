use hashbrown::HashMap;
use log::warn;
use std::{
    cell::{RefCell, RefMut},
    clone,
    ops::{Add, DerefMut, Div, Mul, Sub},
};

use cfd_rs_utils::mesh::{computational_mesh::Computational2DMesh, indices::CellIndex};
use nalgebra::DVector;
use nalgebra_sparse::{CooMatrix, CsrMatrix};
use nalgebra_sparse_linalg::iteratives;

use super::{
    base::Field,
    boundary::{BoundaryCondition, FieldsBoundaryConditions},
    case::{Case, GradRequirements, VariableFields},
    config::{CaseConfig, Schemes},
    discretizations::DifferentialOperator,
    error::CfdError,
};

/// Implementation of the creation of calculation graph for matrix creation (OpenFoam style)

#[derive(Clone, Debug, PartialEq)]
pub enum Op {
    Add(Box<(Op, Op)>),
    Sub(Box<(Op, Op)>),
    Discretize(DifferentialOperator),
    MulScalar(f64, Box<Op>),
    DivScalar(f64, Box<Op>),
    Scalar(f64),
    Field(Field),
}

impl Op {
    pub fn collect_differential_operators(&self, collector: &mut Vec<DifferentialOperator>) {
        match self {
            Op::Add(pair) | Op::Sub(pair) => {
                Self::collect_differential_operators(&pair.0, collector);
                Self::collect_differential_operators(&pair.1, collector);
            }
            Op::MulScalar(_, inner) | Op::DivScalar(_, inner) => {
                Self::collect_differential_operators(inner, collector);
            }
            Op::Discretize(dop) => collector.push(dop.clone()),
            Op::Scalar(_) => {}
        }
    }
}

impl Add for Op {
    type Output = Op;

    fn add(self, rhs: Self) -> Self::Output {
        Op::Add(Box::new((self, rhs)))
    }
}

impl Sub for Op {
    type Output = Op;

    fn sub(self, rhs: Self) -> Self::Output {
        Op::Sub(Box::new((self, rhs)))
    }
}

impl Mul<f64> for Op {
    type Output = Op;

    fn mul(self, rhs: f64) -> Self::Output {
        Op::MulScalar(rhs, Box::new(self))
    }
}

impl Mul<Op> for f64 {
    type Output = Op;

    fn mul(self, rhs: Op) -> Self::Output {
        Op::MulScalar(self, Box::new(rhs))
    }
}

impl Div<f64> for Op {
    type Output = Op;

    fn div(self, rhs: f64) -> Self::Output {
        Op::DivScalar(rhs, Box::new(self))
    }
}

#[derive(Clone, Debug, PartialEq)]
pub enum IntegrationCategory {
    Implicit,
    Explicit,
}

/// Will help to implement unit checking
#[derive(Clone, Debug, PartialEq, Hash, Eq)]
pub struct Variable {
    name: String,
    dim: Dimension,
}

/// For now only support of 1D fields
#[derive(Clone, Debug, PartialEq, Hash, Eq)]
pub enum Dimension {
    Scalar,
    Vector2,
    XComponent,
    YComponent,
}

impl Variable {
    pub fn new(name: String, dim: Dimension) -> Self {
        Variable { name, dim }
    }

    pub fn name(&self) -> &str {
        &self.name
    }
}

/// Implicit formulations are allowed on one variable only
#[derive(Clone, Debug, PartialEq)]
pub struct Equation {
    lhs: Op,
    rhs: Op,
    unknown: Variable,
    fields_required: Vec<Variable>,
}

impl Equation {
    pub fn lhs(&self) -> &Op {
        &self.lhs
    }

    pub fn rhs(&self) -> &Op {
        &self.rhs
    }

    pub fn unknown(&self) -> &Variable {
        &self.unknown
    }

    pub fn collect_differential_operators(&self) -> Vec<DifferentialOperator> {
        let mut collector = vec![];
        self.lhs.collect_differential_operators(&mut collector);
        self.rhs.collect_differential_operators(&mut collector);
        collector
    }

    pub fn new(lhs: Op, rhs: Op, schemes: &Schemes) -> Result<Equation, CfdError> {
        let mut collector = vec![];
        lhs.collect_differential_operators(&mut collector);
        rhs.collect_differential_operators(&mut collector);

        let mut unknown_var = None;

        for diff_op in collector {
            match diff_op {
                DifferentialOperator::TimeDerivative(var) => match &unknown_var {
                    None => unknown_var = Some(var),
                    Some(current) => {
                        if current != &var {
                            return Err(CfdError::EquationInconsistentUnknown {
                                lhs,
                                rhs,
                                var1: current.clone(),
                                var2: var,
                            });
                        }
                    }
                },

                DifferentialOperator::Convection {
                    var, integration, ..
                }
                | DifferentialOperator::Divergence(var, integration)
                | DifferentialOperator::Laplacian(var, integration) => {
                    if let IntegrationCategory::Implicit = integration {
                        match &unknown_var {
                            None => unknown_var = Some(var),
                            Some(current) => {
                                if current != &var {
                                    return Err(CfdError::EquationInconsistentUnknown {
                                        lhs,
                                        rhs,
                                        var1: current.clone(),
                                        var2: var,
                                    });
                                }
                            }
                        }
                    }
                }
            }
        }
        
        let mut variable_requirements = HashMap::new();
        let mut integration = IntegrationCategory::Explicit;
        
        let mut collector = vec![];
        (lhs.clone() - rhs.clone()).collect_differential_operators(&mut collector);
        
        for diff_operator in collector {
            let var = diff_operator.variable();
            let (_, required_grad) = diff_operator.required_grads(schemes);
            variable_requirements
                .entry(var.clone())
                .and_modify(|current: &mut GradRequirements| {
                    current.update_requirements(&required_grad)
                })
                .or_insert(required_grad);
            if diff_operator.integration() == Some(&IntegrationCategory::Implicit) {
                integration = IntegrationCategory::Implicit;
            }
        }

        match unknown_var {
            None => Err(CfdError::EquationNoUnknown { lhs: lhs, rhs: rhs }),
            Some(var) => Ok(Equation {
                lhs,
                rhs,
                unknown: var,
                fields_required: variable_requirements
                    .keys()
                    .map(|var| var.clone())
                    .collect(),
            }),
        }
    }
    
    pub fn apply_op(
        &self,
        op: &Op,
        solver: &mut EquationSolver,
        fields: &Vec<RefMut<Field>>,
        mesh: &Computational2DMesh,
        config: &CaseConfig,
        coeff: f64,
    ) {
        match op {
            Op::Add(op) => {
                self.apply_op(&op.as_ref().0, solver, fields, mesh, config, coeff);
                self.apply_op(&op.as_ref().1, solver, fields, mesh, config, coeff);
            }
            Op::Sub(op) => {
                self.apply_op(&op.as_ref().0, solver, fields, mesh, config, coeff);
                self.apply_op(&op.as_ref().1, solver, fields, mesh, config, -coeff);
            }
            Op::MulScalar(scalar, op) => {
                self.apply_op(op.as_ref(), solver, fields, mesh, config, coeff * scalar);
            }
            Op::DivScalar(scalar, op) => {
                self.apply_op(op.as_ref(), solver, fields, mesh, config, coeff / scalar);
            }
            Op::Discretize(d_op) => {
                d_op.discretize(self, solver, fields, mesh, config, coeff);
            }
            Op::Scalar(scalar) => {
                todo!();
                //solver.add_scalar(*scalar * coeff);
            }
            Op::Field(field) => {
                todo!();
            }
        }
    }
    
    pub fn solve(&self, solver: &mut EquationSolver, fields: &mut VariableFields, mesh: &Computational2DMesh, config: &CaseConfig) {
        
    }

}

pub fn 

pub struct EquationSolver {
    matrix: CsrMatrix<f64>,
    rhs: DVector<f64>,
}

impl EquationSolver {
    pub fn new(mesh: &Computational2DMesh) -> Self {
        let mut matrix = CooMatrix::new(mesh.num_cells(), mesh.num_cells());

        for i in 0..mesh.num_cells() {
            matrix.push(i, i, 0.);
            for neighbor in mesh.neighboring_cells_id(CellIndex(i)) {
                matrix.push(i, neighbor.0, 0.)
            }
        }

        let matrix = CsrMatrix::from(&matrix);

        let rhs = DVector::zeros(mesh.num_cells());
        
        EquationSolver { matrix, rhs}
    }
    
    pub fn clear(&mut self) {
        for v in self.matrix.values_mut() {
            *v = 0.;
        }
        for v in self.rhs.iter_mut() {
            *v = 0.;
        }
    }

    pub fn matrix(&self) -> &CsrMatrix<f64> {
        &self.matrix
    }

    pub fn matrix_mut(&mut self) -> &mut CsrMatrix<f64> {
        &mut self.matrix
    }

    pub fn rhs(&self) -> &DVector<f64> {
        &self.rhs
    }

    pub fn rhs_mut(&mut self) -> &mut DVector<f64> {
        &mut self.rhs
    }

    pub fn solver_borrow_mut(&mut self) -> (&mut CsrMatrix<f64>, &mut DVector<f64>) {
        (&mut self.matrix, &mut self.rhs)
    }
    
    fn add_scalar(&mut self, scalar: f64) {
        for v in self.rhs.iter_mut() {
            *v -= scalar;
        }
    }
    
    
}

impl System {

    

    fn add_scalar(&mut self, scalar: f64) {
        for v in self.rhs.iter_mut() {
            *v -= scalar;
        }
    }

    pub fn solve<T: Case>(case: &mut T, name: &str) -> usize {
        let (systems, variable_fields, mesh, config) = case.equation_solver_borrow();

        let system = systems
            .map
            .get_mut(name)
            .expect(&format!("This equation is not defined: {name:?}"));

        let fields_grad_required: Vec<&(RefCell<Field>, GradRequirements)> = system
            .fields_required
            .iter()
            .map(|var| {
                variable_fields.map.get(var).expect(&format!(
                    "A field ({var:?}) needed for equation {name:?} is missing"
                ))
            })
            .collect();

        let grad_requirements = fields_grad_required.iter().map(|tuple| &tuple.1).collect();
        let fields_required: Vec<RefMut<Field>> = fields_grad_required
            .iter()
            .map(|tuple| tuple.0.borrow_mut())
            .collect();

        solve(system, fields_required, grad_requirements, mesh, config)
    }
}

fn solve(
    system: &mut System,
    mut fields: Vec<RefMut<Field>>,
    grad_requirements: Vec<&GradRequirements>,
    mesh: &Computational2DMesh,
    config: &CaseConfig,
) -> usize {
    system.clear();

    let lhs = system.equation.lhs.clone();
    let rhs = system.equation.rhs.clone();

    let eq = lhs - rhs;

    for (i, field) in fields.iter_mut().enumerate() {
        field.update_grads(
            grad_requirements[i],
            mesh,
            &config.schemes.gradients,
            config
                .bc
                .map
                .get(&system.fields_required()[i])
                .expect("Boundary Condition missing for field"),
        );
    }
    
    match system.equation().unknown().dim {
        Dimension::Scalar => {
            let dim = Dimension::Scalar;
            system.apply_op(&eq, &fields, mesh, &config.bc, &config.schemes, 1., dim);

            let index = system
                .fields_required()
                .iter()
                .position(|var_i| system.equation().unknown() == var_i)
                .expect(&format!(
                    "Missing variable {:?} in fields for equation {:?}",
                    system.equation().unknown(),
                    system.equation()
                ));
            

            let field = fields[index].deref_mut();

            let field = match field {
                Field::Scalar(value) => value,
                _ => panic!("Unknown should be scalar"),
            };
            
            warn!("Hard-coded tol and max_iter for solve");
            let result = iteratives::biconjugate_gradient::solve_with_initial_guess(
                system.matrix(),
                &system.rhs,
                field.values_mut(),
                1000,
                1e-3,
            );

            if !result {
                panic!("Did not converge when solving {:?}", system.equation())
            }
        },
        Dimension::Vector2 => {
            
            system.apply_op(&eq, &fields, mesh, &config.bc, &config.schemes, 1.);

            let index = system
                .fields_required()
                .iter()
                .position(|var_i| system.equation().unknown() == var_i)
                .expect(&format!(
                    "Missing variable {:?} in fields for equation {:?}",
                    system.equation().unknown(),
                    system.equation()
                ));
            

            let field = fields[index].deref_mut();

            let field = match field {
                Field::Vector2(value) => value,
                _ => panic!("Unknown should be scalar"),
            };
            
            warn!("Hard-coded tol and max_iter for solve");
            let result = iteratives::biconjugate_gradient::solve_with_initial_guess(
                system.matrix(),
                &system.rhs,
                field.values_mut(),
                1000,
                1e-3,
            );

            if !result {
                panic!("Did not converge when solving {:?}", system.equation())
            }
        },
        _ => (),
    }
    
    

    0
}
