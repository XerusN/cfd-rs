use nalgebra::Vector2;
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, PartialEq, Deserialize, Serialize)]
pub struct Pairs {
    centers: Vec<Vector2<f64>>,
    length: Vec<f64>,
    normals: Vec<Vec<Vector2<f64>>>,
    neighboring_nodes: Vec<[usize; 2]>,
    neighboring_cells: Vec<[usize; 2]>,
}