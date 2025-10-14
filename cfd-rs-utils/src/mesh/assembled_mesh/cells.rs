use nalgebra::{Point2, Vector2};
use serde::{Deserialize, Serialize};

use crate::mesh::assembled_mesh::core::Patch;

#[derive(Debug, Clone, PartialEq, Deserialize, Serialize)]
pub struct Cells {
    centers: Vec<Point2<f64>>,
    volumes: Vec<f64>,
    areas: Vec<Vec<f64>>,
    normals: Vec<Vec<Vector2<f64>>>,
    neighboring_nodes: Vec<Vec<usize>>,
    neighboring_cells: Vec<Vec<Patch>>,
    neighboring_pairs: Vec<Vec<usize>>,
}

impl Cells {
    pub fn new(
        centers: Vec<Point2<f64>>,
        volumes: Vec<f64>,
        areas: Vec<Vec<f64>>,
        normals: Vec<Vec<Vector2<f64>>>,
        neighboring_nodes: Vec<Vec<usize>>,
        neighboring_cells: Vec<Vec<Patch>>,
        neighboring_pairs: Vec<Vec<usize>>,
    ) -> Self {
        Cells {
            centers,
            volumes,
            areas,
            normals,
            neighboring_nodes,
            neighboring_cells,
            neighboring_pairs,
        }
    }
}
