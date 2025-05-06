use cfd_rs_utils::mesh::{computational_mesh::*, indices::FaceIndex};

#[derive(Clone, Debug, Default, PartialEq)]
pub enum InterpolationConfig {
    #[default]
    Linear,
    RhieChow,
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

pub fn linear_interp(face_index: FaceIndex, mut field: CellFaceField, mesh: Computational2DMesh) {
    todo!()
}