use nalgebra::{Point2, Vector2};
use serde::{Deserialize, Serialize};

use crate::mesh::assembled_mesh::core::Patch;

#[derive(Debug, Clone, PartialEq, Deserialize, Serialize)]
pub struct Nodes {
    centers: Vec<Point2<f64>>,
    volumes: Vec<f64>,
    areas: Vec<Vec<f64>>,
    normals: Vec<Vec<Vector2<f64>>>,
    neighboring_nodes: Vec<Vec<usize>>,
    neighboring_cells: Vec<Vec<Patch>>,
    neighboring_pairs: Vec<Vec<usize>>,
}

impl Nodes {
    pub fn new(
        centers: Vec<Point2<f64>>,
        volumes: Vec<f64>,
        areas: Vec<Vec<f64>>,
        normals: Vec<Vec<Vector2<f64>>>,
        neighboring_nodes: Vec<Vec<usize>>,
        neighboring_cells: Vec<Vec<Patch>>,
        neighboring_pairs: Vec<Vec<usize>>,
    ) -> Self {
        Nodes {
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
