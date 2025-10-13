use nalgebra::{Point2, Vector2};
use serde::{Deserialize, Serialize};

use crate::mesh::assembled_mesh::core::Patch;

#[derive(Debug, Clone, PartialEq, Deserialize, Serialize)]
pub struct Pairs {
    centers: Vec<Point2<f64>>,
    lengths: Vec<f64>,
    normals: Vec<Vector2<f64>>,
    nodes: Vec<[usize; 2]>,
    neighboring_cells: Vec<[Patch; 2]>,
}

impl Pairs {
    pub fn new(
        centers: Vec<Point2<f64>>,
        lengths: Vec<f64>,
        normals: Vec<Vector2<f64>>,
        nodes: Vec<[usize; 2]>,
        neighboring_cells: Vec<[Patch; 2]>,
    ) -> Self {
        Pairs {
            centers,
            lengths,
            normals,
            nodes,
            neighboring_cells,
        }
    }
}
