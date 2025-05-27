use hashbrown::HashMap;
use log::warn;
use std::{
    cell::{RefCell, RefMut},
    clone,
    ops::{Add, Deref, DerefMut, Div, Mul, Sub},
};

use cfd_rs_utils::mesh::{computational_mesh::Computational2DMesh, indices::CellIndex};
use nalgebra::DVector;
use nalgebra_sparse::{csr::CsrRowMut, CooMatrix, CsrMatrix};
use nalgebra_sparse_linalg::iteratives;

use super::{
    base::Field,
    boundary::{BoundaryCondition, FieldsBoundaryConditions},
    case::{Case, GradRequirements, VariableFields},
    config::{CaseConfig, Schemes},
    discretizations::{find_var_in_fields, DifferentialOperator},
    error::CfdError,
};

/// Implementation of the creation of calculation graph for matrix creation (OpenFoam style)

#[derive(Clone, Debug, PartialEq)]
pub enum Op {
    Add(Box<(Op, Op)>),
    Sub(Box<(Op, Op)>),
    FieldOperator(FieldOperator),
    MulScalar(f64, Box<Op>),
    DivScalar(f64, Box<Op>),
    Scalar(f64),
}

impl Op {
    pub fn collect_field_operators(&self, collector: &mut Vec<FieldOperator>) {
        match self {
            Op::Add(pair) | Op::Sub(pair) => {
                Self::collect_field_operators(&pair.0, collector);
                Self::collect_field_operators(&pair.1, collector);
            }
            Op::MulScalar(_, inner) | Op::DivScalar(_, inner) => {
                Self::collect_field_operators(inner, collector);
            }
            Op::FieldOperator(f_op) => collector.push(f_op.clone()),
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

#[derive(Clone, Debug, PartialEq, Hash, Eq)]
pub enum Dimension {
    Scalar,
    Vector2,
}

#[derive(Clone, Debug, PartialEq, Hash, Eq)]
pub enum Component {
    X,
    Y,
}

impl Variable {
    pub fn new(name: String, dim: Dimension) -> Self {
        Variable { name, dim }
    }

    pub fn name(&self) -> &str {
        &self.name
    }

    pub fn dim(&self) -> &Dimension {
        &self.dim
    }
}

#[derive(Clone, Debug, PartialEq)]
pub enum FieldOperator {
    DifferentialOperator(DifferentialOperator),
    Field(Variable, IntegrationCategory),
    Gradient(Variable),
}

impl FieldOperator {
    pub fn variable(&self) -> &Variable {
        match self {
            FieldOperator::DifferentialOperator(diff_op) => diff_op.variable(),
            FieldOperator::Gradient(var) => &var,
            FieldOperator::Field(var, _) => &var,
        }
    }

    pub fn required_grads(&self, schemes: &Schemes) -> (&Variable, GradRequirements) {
        match self {
            FieldOperator::DifferentialOperator(diff_op) => diff_op.required_grads(schemes),
            FieldOperator::Gradient(var) => (&var, GradRequirements::new(true, false)),
            FieldOperator::Field(var, _) => (&var, GradRequirements::new(false, false)),
        }
    }
}

/// Implicit formulations are allowed on one variable only
#[derive(Clone, Debug, PartialEq)]
pub struct Equation {
    lhs: Op,
    rhs: Op,
    unknown: Variable,
    variables_requirements: HashMap<Variable, GradRequirements>,
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

    pub fn fields_required(&self) -> &HashMap<Variable, GradRequirements> {
        &self.variables_requirements
    }

    pub fn collect_field_operators(&self) -> Vec<FieldOperator> {
        let mut collector = vec![];
        self.lhs.collect_field_operators(&mut collector);
        self.rhs.collect_field_operators(&mut collector);
        collector
    }

    pub fn new(lhs: Op, rhs: Op, schemes: &Schemes) -> Result<Equation, CfdError> {
        let mut unknown_var = None;

        let mut variables_requirements = HashMap::new();

        let mut collector = vec![];
        (lhs.clone() - rhs.clone()).collect_field_operators(&mut collector);

        for f_op in &collector {
            match f_op {
                FieldOperator::DifferentialOperator(diff_op) => match diff_op {
                    DifferentialOperator::TimeDerivative(var) => match &unknown_var {
                        None => unknown_var = Some(var),
                        Some(current) => {
                            if current != &var {
                                return Err(CfdError::EquationInconsistentUnknown {
                                    lhs,
                                    rhs,
                                    var1: current.clone().clone(),
                                    var2: var.clone(),
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
                                            var1: current.clone().clone(),
                                            var2: var.clone(),
                                        });
                                    }
                                }
                            }
                        }
                    }
                },
                FieldOperator::Gradient(_) => (),
                FieldOperator::Field(var, integration) => {
                    if let IntegrationCategory::Implicit = integration {
                        match &unknown_var {
                            None => unknown_var = Some(var),
                            Some(current) => {
                                if current != &var {
                                    return Err(CfdError::EquationInconsistentUnknown {
                                        lhs,
                                        rhs,
                                        var1: current.clone().clone(),
                                        var2: var.clone(),
                                    });
                                }
                            }
                        }
                    }
                }
            }
        }

        for f_op in &collector {
            let var = f_op.variable();
            let (_, required_grad) = f_op.required_grads(schemes);
            variables_requirements
                .entry(var.clone())
                .and_modify(|current: &mut GradRequirements| {
                    current.update_requirements(&required_grad)
                })
                .or_insert(required_grad);
        }

        match unknown_var {
            None => Err(CfdError::EquationNoUnknown { lhs: lhs, rhs: rhs }),
            Some(var) => Ok(Equation {
                lhs,
                rhs,
                unknown: var.clone(),
                variables_requirements,
            }),
        }
    }

    pub fn solve<T: Case>(case: &mut T, name: &str) -> usize {
        let time_step = case.time_step();
        let (solver, variable_fields, equations, mesh, config) = case.equation_solver_borrow();

        let equation = equations
            .map
            .get(name)
            .expect(&format!("This equation is not defined: {name:?}"));

        solve(solver, equation, variable_fields, mesh, config, time_step)
    }
}

#[derive(Clone, Debug, PartialEq)]
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

        EquationSolver { matrix, rhs }
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

    fn add_gradient(&mut self, var: &Variable, component: &Component, fields: &VariableFields, coeff: f64) {
        let field = find_var_in_fields(var, fields);
        let field = field.borrow();
        let field = match *field {
            Field::Scalar(ref values) => values,
            Field::Vector2(_) => panic!("Can't add the gradient of a vector field"),
        };
        
        match component {
            Component::X => for (i, v) in self.rhs.iter_mut().enumerate() {
                *v -= field.grads_cell()[i].x;
            },
            Component::Y => for (i, v) in self.rhs.iter_mut().enumerate() {
                *v -= field.grads_cell()[i].y;
            },
        }
    }

    pub fn apply_op(
        &mut self,
        op: &Op,
        component: &Component,
        fields: &VariableFields,
        mesh: &Computational2DMesh,
        config: &CaseConfig,
        time_step: f64,
        coeff: f64,
    ) {
        match op {
            Op::Add(op) => {
                self.apply_op(
                    &op.as_ref().0,
                    component,
                    fields,
                    mesh,
                    config,
                    time_step,
                    coeff,
                );
                self.apply_op(
                    &op.as_ref().1,
                    component,
                    fields,
                    mesh,
                    config,
                    time_step,
                    coeff,
                );
            }
            Op::Sub(op) => {
                self.apply_op(
                    &op.as_ref().0,
                    component,
                    fields,
                    mesh,
                    config,
                    time_step,
                    coeff,
                );
                self.apply_op(
                    &op.as_ref().1,
                    component,
                    fields,
                    mesh,
                    config,
                    time_step,
                    -coeff,
                );
            }
            Op::MulScalar(scalar, op) => {
                self.apply_op(
                    op.as_ref(),
                    component,
                    fields,
                    mesh,
                    config,
                    time_step,
                    coeff * scalar,
                );
            }
            Op::DivScalar(scalar, op) => {
                self.apply_op(
                    op.as_ref(),
                    component,
                    fields,
                    mesh,
                    config,
                    time_step,
                    coeff / scalar,
                );
            }
            Op::FieldOperator(f_op) => match f_op {
                FieldOperator::DifferentialOperator(diff_op) => {
                    diff_op.discretize(component, self, fields, mesh, config, time_step, coeff)
                }
                FieldOperator::Field(var, integration) => todo!(),
                FieldOperator::Gradient(var) => self.add_gradient(var, component, fields, coeff),
            },
            Op::Scalar(scalar) => {
                self.add_scalar(*scalar * coeff);
            }
        }
    }
}

fn solve(
    solver: &mut EquationSolver,
    equation: &Equation,
    fields: &mut VariableFields,
    mesh: &Computational2DMesh,
    config: &CaseConfig,
    time_step: f64,
) -> usize {
    solver.clear();

    let lhs = equation.lhs.clone();
    let rhs = equation.rhs.clone();

    let eq = lhs - rhs;

    for (var, field) in &fields.map {
        field.0.borrow_mut().update_grads(
            &field.1,
            mesh,
            &config.schemes.gradients,
            config
                .bc
                .map
                .get(var)
                .expect("Boundary Condition missing for field"),
        );
    }

    match equation.unknown().dim {
        Dimension::Scalar => {
            let component = Component::X;
            solver.apply_op(&eq, &component, &fields, mesh, config, time_step, 1.);

            let field_cell = &fields
                .map
                .get_mut(equation.unknown())
                .expect("Missing field for equation")
                .0;

            let mut field = field_cell.borrow_mut();

            let field = match &mut *field {
                Field::Scalar(ref mut scalar_field) => scalar_field,
                _ => panic!("Unknown should be scalar"),
            };

            // To change
            warn!("Matrix cloned for amg");
            warn!("Hard-coded tol and max_iter for solve");
            let result = iteratives::amg::solve_with_initial_guess(
                solver.matrix().clone(),
                &solver.rhs,
                field.values_mut(),
                10000,
                1e-1,
                0.25,
            );

            if !result {
                panic!("Did not converge when solving {:?}", equation)
            }
        }
        Dimension::Vector2 => {
            for component in [Component::X, Component::Y] {
                solver.apply_op(&eq, &component, &fields, mesh, config, time_step, 1.);

                let field_cell = &fields
                    .map
                    .get_mut(equation.unknown())
                    .expect("Missing field for equation")
                    .0;

                let mut field = field_cell.borrow_mut();

                let field = match &mut *field {
                    Field::Vector2(ref mut scalar_field) => match component {
                        Component::X => &mut scalar_field.x,
                        Component::Y => &mut scalar_field.y,
                    },
                    _ => panic!("Unknown should be vector"),
                };

                // To change
                warn!("Matrix cloned for amg");
                warn!("Hard-coded tol and max_iter for solve");
                let result = iteratives::amg::solve_with_initial_guess(
                    solver.matrix().clone(),
                    &solver.rhs,
                    field.values_mut(),
                    10000,
                    1e-1,
                    0.8,
                );

                if !result {
                    panic!("Did not converge when solving {:?}", equation)
                }
            }
        }
    }

    0
}
