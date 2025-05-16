use std::{
    cell::{RefCell, RefMut},
    collections::HashMap,
    ops::{Add, Div, Mul, Sub},
};

use cfd_rs_utils::mesh::{computational_mesh::Computational2DMesh, indices::CellIndex};
use nalgebra::DVector;
use nalgebra_sparse::{CooMatrix, CsrMatrix};

use super::{
    base::Field,
    case::{Case, GradRequirements},
    config::Schemes,
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
    //Vector2,
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

    pub fn new(lhs: Op, rhs: Op) -> Result<Equation, CfdError> {
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

        match unknown_var {
            None => Err(CfdError::EquationNoUnknown { lhs, rhs }),
            Some(var) => Ok(Equation {
                lhs,
                rhs,
                unknown: var,
            }),
        }
    }

    pub fn into_system<T: Case>(self, case: &T) -> (System, HashMap<Variable, GradRequirements>) {
        System::new(self, case)
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct System {
    equation: Equation,
    matrix: CsrMatrix<f64>,
    rhs: DVector<f64>,
    int_cat: IntegrationCategory,
    fields_required: Vec<Variable>,
}

impl System {
    pub fn new<T: Case>(
        equation: Equation,
        case: &T,
    ) -> (System, HashMap<Variable, GradRequirements>) {
        let mut matrix = CooMatrix::new(case.mesh().num_cells(), case.mesh().num_cells());

        for i in 0..case.mesh().num_cells() {
            matrix.push(i, i, 0.);
            for neighbor in case.mesh().neighboring_cells_id(CellIndex(i)) {
                matrix.push(i, neighbor.0, 0.)
            }
        }

        let matrix = CsrMatrix::from(&matrix);

        let rhs = DVector::zeros(case.mesh().num_cells());

        let mut variable_requirements = HashMap::new();
        let mut integration = IntegrationCategory::Explicit;

        for diff_operator in equation.collect_differential_operators() {
            let var = diff_operator.variable();
            let (_, required_grad) = diff_operator.required_grads(case.schemes());
            variable_requirements
                .entry(var.clone())
                .and_modify(|current: &mut GradRequirements| {
                    current.update_requirements(required_grad.clone())
                })
                .or_insert(required_grad);
            if diff_operator.integration() == Some(&IntegrationCategory::Implicit) {
                integration = IntegrationCategory::Implicit;
            }
        }

        (
            System {
                equation,
                matrix,
                rhs,
                int_cat: integration,
                fields_required: variable_requirements
                    .keys()
                    .map(|var| var.clone())
                    .collect(),
            },
            variable_requirements,
        )
    }

    pub fn clear(&mut self) {
        for v in self.matrix.values_mut() {
            *v = 0.;
        }
        for v in self.rhs.iter_mut() {
            *v = 0.;
        }
    }

    pub fn equation(&self) -> &Equation {
        &self.equation
    }

    pub fn equation_mut(&mut self) -> &mut Equation {
        &mut self.equation
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

    pub fn integration_category(&self) -> &IntegrationCategory {
        &self.int_cat
    }

    pub fn fields_required(&self) -> &[Variable] {
        &self.fields_required
    }

    pub fn apply_op(
        &mut self,
        op: &Op,
        fields: &Vec<RefMut<Field>>,
        mesh: &Computational2DMesh,
        schemes: &Schemes,
        coeff: f64,
    ) {
        match op {
            Op::Add(op) => {
                self.apply_op(&op.as_ref().0, fields, mesh, schemes, coeff);
                self.apply_op(&op.as_ref().1, fields, mesh, schemes, coeff);
            }
            Op::Sub(op) => {
                self.apply_op(&op.as_ref().0, fields, mesh, schemes, coeff);
                self.apply_op(&op.as_ref().1, fields, mesh, schemes, -coeff);
            }
            Op::MulScalar(scalar, op) => {
                self.apply_op(op.as_ref(), fields, mesh, schemes, coeff * scalar);
            }
            Op::DivScalar(scalar, op) => {
                self.apply_op(op.as_ref(), fields, mesh, schemes, coeff / scalar);
            }
            Op::Discretize(d_op) => {
                d_op.discretize(self, fields, schemes, coeff);
            }
            Op::Scalar(scalar) => {
                self.add_scalar(*scalar * coeff);
            }
        }
    }

    fn add_scalar(&mut self, scalar: f64) {
        for v in self.rhs.iter_mut() {
            *v -= scalar;
        }
    }

    pub fn solve<T: Case>(case: &mut T, name: &str) -> usize {
        let (systems, variable_fields, mesh, config) = case.equation_solver_borrow();

        let schemes = &config.schemes;

        let system = systems
            .map
            .get_mut(name)
            .expect(&format!("This equation is not defined: {name:?}"));

        let fields_grad_required: Vec<&(RefCell<Field>, GradRequirements)> = system
            .fields_required
            .iter()
            .map(|var| {
                variable_fields
                    .map
                    .get(var)
                    .expect(&format!("A field ({var:?}) needed for equation {name:?}"))
            })
            .collect();

        let fields_required: Vec<RefMut<Field>> = fields_grad_required
            .iter()
            .map(|tuple| tuple.0.borrow_mut())
            .collect();

        solve(system, fields_required, mesh, schemes)
    }
}

fn solve(
    system: &mut System,
    fields: Vec<RefMut<Field>>,
    mesh: &Computational2DMesh,
    schemes: &Schemes,
) -> usize {
    system.clear();

    let lhs = system.equation.lhs.clone();
    let rhs = system.equation.rhs.clone();

    let eq = lhs - rhs;

    system.apply_op(&eq, &fields, mesh, schemes, 1.);

    // Add linear solver

    0
}
