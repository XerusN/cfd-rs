use cfd_rs_utils::mesh::{computational_mesh::*, indices::FaceIndex};

use super::base::ScalarVariable;

#[derive(Clone, Debug, Default, PartialEq)]
pub enum InterpolationConfig {
    #[default]
    Linear,
    RhieChow,
}

impl InterpolationConfig {
    pub fn interpolate(&self) -> fn() {
        match *self {
            Self::Linear => todo!(),
            _ => panic!("{self:?} not implemented"),
        }
    }
}

#[derive(Clone, Debug, Default, PartialEq)]
pub enum GradientInterpConfig {
    Averaged,
    #[default]
    AveragedCorrected,
}

impl GradientInterpConfig {
    pub fn gradient_interp(&self) -> fn(&Computational2DMesh, &mut ScalarVariable) {
        match *self {
            Self::AveragedCorrected => averaged_corrected,
            _ => panic!("{self:?} not implemented"),
        }
    }
}

pub fn averaged_corrected(mesh: &Computational2DMesh, value: &mut ScalarVariable) {
    todo!()
}

#[derive(Clone, Debug, Default, PartialEq)]
pub enum DecompositionConfig {
    #[default]
    MinimumCorrection,
    OrthogonalCorrection,
    OverRelaxed,
}

impl DecompositionConfig {
    pub fn decompose(&self) -> fn() {
        match *self {
            Self::OverRelaxed => todo!(),
            _ => panic!("{self:?} not implemented"),
        }
    }
}


// pub fn linear_interp(face_index: FaceIndex, mut field: CellFaceField, mesh: Computational2DMesh) {
//     todo!()
// }