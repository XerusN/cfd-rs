use nalgebra::Vector2;

#[derive(Clone, Debug, PartialEq, Default)]
pub struct ScalarVariable {
    pub value: Vec<f64>,
    pub grad_cell: Option<Vector2<Vec<f64>>>,
    pub grad_faces: Option<Vector2<Vec<f64>>>,
}