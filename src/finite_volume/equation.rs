use std::{
    collections::HashMap,
    ops::{Add, Div, Mul, Sub},
};

use cfd_rs_utils::mesh::indices::CellIndex;
use nalgebra::DVector;
use nalgebra_sparse::{CooMatrix, CsrMatrix};

use super::{
    case::{Case, GradRequirements}, config::Schemes, discretizations::{laplacian::LaplacianScheme, time_schemes::TimeIntegration}, error::CfdError
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

#[derive(Clone, Debug, PartialEq)]
pub enum DifferentialOperator {
    Laplacian(Variable, IntegrationCategory),
    Convection {
        var: Variable,
        speed: Variable,
        integration: IntegrationCategory,
    },
    Divergence(Variable, IntegrationCategory),
    TimeDerivative(Variable),
}

impl DifferentialOperator {
    pub fn required_grads(&self, schemes: &Schemes) -> (&Variable, GradRequirements) {
        match self {
            Self::Laplacian(var, _) => (var, schemes.laplacian.required_grads()),
            Self::Convection{var, ..} => (var, schemes.convection.required_grads()),
            Self::Divergence(var, _) => (var, schemes.divergence.required_grads()),
            Self::TimeDerivative(var) => (var, schemes.transient.required_grads()),
        }
    }
    
    pub fn variable(&self) -> &Variable {
        match self {
            Self::Laplacian(var, _) => &var,
            Self::Convection{var, ..} => &var,
            Self::Divergence(var, _) => &var,
            Self::TimeDerivative(var) => &var,
        }
    }
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
        
        for diff_operator in equation.collect_differential_operators() {
            let var = diff_operator.variable();
            let (_, required_grad)= diff_operator.required_grads(case.schemes());
            variable_requirements.entry(var.clone()).and_modify(|current: &mut GradRequirements| current.update_requirements(required_grad.clone())).or_insert(required_grad);
        }
        
        (
            System {
                equation,
                matrix,
                rhs,
            },
            variable_requirements,
        )
        
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

    pub fn solve<T: Case>(case: &mut T, name: &str) {
        todo!()
    }
}
