use nalgebra::Vector2;

/// For now will be fully allocated every time,
/// Can be optimised with Option for gradients
#[derive(Clone, Debug, PartialEq, Default)]
pub struct ScalarVariable {
    pub value: Vec<f64>,
    pub grad_cell: Vector2<Vec<f64>>,
    pub grad_faces: Vector2<Vec<f64>>,
}