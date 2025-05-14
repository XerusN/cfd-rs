use std::ops::{Add, Sub, Mul, Div};

use nalgebra_sparse::CsrMatrix;
use nalgebra::DVector;

use super::case::Case;

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
pub enum DifferentialOperator {
    Laplacian(Variable),
    Convection {var: Variable, speed: Variable},
    Divergence(Variable),
    TimeDerivative(Variable),
}

#[derive(Clone, Debug, PartialEq, Hash, Eq)]
pub struct Variable {
    name: String,
    grads_cell_required: bool,
    grads_face_required: bool,
}

impl Variable {
    pub fn new(name: String, grads_cell_required: bool, grads_face_required: bool) -> Self {
        Variable{name, grads_cell_required, grads_face_required}
    }
    
    pub fn name(&self) -> &str {
        &self.name
    }
    
    pub fn update_grads_requirements(&mut self, other: &Variable) {
        self.grads_cell_required = self.grads_cell_required | other.grads_cell_required;
        self.grads_face_required = self.grads_face_required | other.grads_face_required;
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct Equation {
    lhs: Op,
    rhs: Op,
    unknown: Variable,
}

impl Equation {
    
    pub fn new(lhs: Op, rhs: Op, unknown: Variable) -> Equation {
        Equation { lhs, rhs, unknown }
    }

    pub fn into_system<T: Case>(self, case: &T) -> (System, Vec<Variable>) {
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
    pub fn new<T: Case>(equation: Equation, case: &T) -> (System, Vec<Variable>) {
        todo!()
        
        //System { equation, matrix: (), rhs: () }
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



