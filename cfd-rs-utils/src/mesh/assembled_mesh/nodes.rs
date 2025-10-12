use nalgebra::Vector2;
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, PartialEq, Deserialize, Serialize)]
pub struct Nodes {
    centers: Vec<Vector2<f64>>,
    volumes: Vec<f64>,
    areas: Vec<Vec<Vector2<f64>>>,
    normals: Vec<Vec<Vector2<f64>>>,
    neighboring_nodes: Vec<Vec<usize>>,
    neighboring_cells: Vec<Vec<usize>>,
    neighboring_pairs: Vec<Vec<usize>>,
}