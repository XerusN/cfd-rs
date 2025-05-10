use cfd_rs_utils::mesh::computational_mesh::Computational2DMesh;

use crate::finite_volume::interpolations::GradientInterpConfig;

use super::base::ScalarVariable;

#[derive(Clone, Debug, PartialEq)]
pub enum GradientConfig {
    GreenGaussCompact,
    GreenGaussExtended,
    LeastSquare,
}

impl GradientConfig {
    
    pub fn compute_gradient(&self) -> fn(&Computational2DMesh, &mut ScalarVariable) {
        match *self {
            Self::GreenGaussCompact => green_gauss_compact,
            _ => panic!("{self:?} not implemented"),
        }
    }
}

pub fn green_gauss_compact(mesh: &Computational2DMesh, values: &mut ScalarVariable) {
    todo!()
}