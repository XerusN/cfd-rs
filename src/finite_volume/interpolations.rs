use cfd_rs_utils::mesh::{computational_mesh::*, indices::FaceIndex};

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

#[derive(Clone, Debug, Default, PartialEq)]
pub enum DecompositionConfig {
    #[default]
    MinimumCorrection,
    OrthogonalCorrection,
    OverRelaxed,
}

// pub fn linear_interp(face_index: FaceIndex, mut field: CellFaceField, mesh: Computational2DMesh) {
//     todo!()
// }