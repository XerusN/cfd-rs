use nalgebra::{Point2, Vector2};
use serde::{Deserialize, Serialize};

use crate::mesh::assembled_mesh::core::Patch;

#[derive(Debug, Clone, PartialEq, Deserialize, Serialize)]
pub struct Pairs {
    pub n: usize,
    centers: Vec<Point2<f64>>,
    lengths: Vec<f64>,
    normals: Vec<Vector2<f64>>,
    nodes: Vec<[usize; 2]>,
    neighboring_cells: Vec<[Patch; 2]>,
}

impl Pairs {
    pub fn centers(&self) -> &[Point2<f64>] {
        &self.centers
    }
    
    pub fn lengths(&self) -> &[f64] {
        &self.lengths
    }
    
    pub fn normals(&self) -> &[Vector2<f64>] {
        &self.normals
    }
    
    pub fn nodes(&self) -> &[[usize; 2]] {
        &self.nodes
    }
    
    pub fn neighboring_cells(&self) -> &[[Patch; 2]] {
        &self.neighboring_cells
    }
}

impl Pairs {
    pub fn new(
        centers: Vec<Point2<f64>>,
        lengths: Vec<f64>,
        normals: Vec<Vector2<f64>>,
        nodes: Vec<[usize; 2]>,
        neighboring_cells: Vec<[Patch; 2]>,
    ) -> Self {
        let n = centers.len();
        assert_eq!(n, lengths.len());
        assert_eq!(n, normals.len());
        assert_eq!(n, nodes.len());
        assert_eq!(n, neighboring_cells.len());
        Pairs {
            n,
            centers,
            lengths,
            normals,
            nodes,
            neighboring_cells,
        }
    }
}
