use nalgebra::Vector2;

use super::gradients::{CellGrad, FaceGrad};

#[derive(Debug, PartialEq, Clone)]
pub enum Field {
    CellScalarField(CellScalarField),
    CellVectorField(Vector2<CellScalarField>),
}

/// The grads will only be allocated if necessary
#[derive(Debug, PartialEq, Clone)]
pub struct CellScalarField {
    pub value: Vec<f64>,
    pub cell_grad: Vec<Vector2<f64>>,
    pub face_grad: Vec<Vector2<f64>>
}