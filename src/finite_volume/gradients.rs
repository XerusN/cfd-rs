use nalgebra::Vector2;
use super::{case::Case, equation::*};

#[derive(Clone, Debug, PartialEq)]
pub enum GradientScheme {
    GreenGaussCompact,
    GreenGaussExtended,
    LeastSquare,
}

#[derive(Clone, Debug, PartialEq)]
pub struct FaceGrad {
    pub var: Variable,
    pub value: Vec<Vector2<f64>>,
}

impl FaceGrad {
    /// Does not compute gradient
    pub fn new<T: Case>(case: &T, var: Variable) -> FaceGrad {
        // Generate a warning if no VolGrad is created for the same variable
        todo!()
    }
    
    pub fn interpolate<T: Case>(&mut self, case: &T) {
        todo!()
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct CellGrad {
    pub var: Variable,
    pub value: Vec<Vector2<f64>>,
}

impl CellGrad {
    /// Does not compute gradient
    pub fn new<T: Case>(case: &T, var: Variable) -> CellGrad {
        todo!()
    }
    
    pub fn update<T: Case>(&mut self, case: &T) {
        todo!()
    }
}