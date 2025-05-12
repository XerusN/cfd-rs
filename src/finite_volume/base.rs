use nalgebra::Vector2;

use super::gradients::{CellGrad, FaceGrad};

#[derive(Debug, PartialEq, Clone)]
pub enum Field {
    CellScalarField{value: Vec<f64>, cell_grad: CellGrad, face_grad: FaceGrad},
    FaceField(Vec<f64>),
    Constant(f64),
}