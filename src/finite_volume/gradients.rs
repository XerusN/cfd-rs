use cfd_rs_utils::mesh::computational_mesh::Computational2DMesh;
use nalgebra::Vector2;
use super::{case::Case, equation::*};

use crate::finite_volume::interpolations::GradientInterpConfig;

use super::base::ScalarVariable;

#[derive(Clone, Debug, PartialEq)]
pub enum GradientConfig {
    GreenGaussCompact,
    GreenGaussExtended,
    LeastSquare,
}

#[derive(Clone, Debug, PartialEq)]
pub struct VolGrad {
    pub var: Variable,
    pub value: Vec<Vector2<f64>>,
}

impl VolGrad {
    /// Does not compute gradient
    pub fn new<T: Case>(case: &T, var: Variable) -> VolGrad {
        todo!()
    }
    
    pub fn update<T: Case>(&mut self, case: &T) {
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
        // Generate a warning if no VolGrad is created for the same variable
        todo!()
    }
    
    pub fn interpolate<T: Case>(&mut self, case: &T) {
        todo!()
    }
}