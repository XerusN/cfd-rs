use nalgebra::{Point2, Vector2};
use serde::{Deserialize, Serialize};

use crate::mesh::assembled_mesh::core::Patch;

#[derive(Debug, Clone, PartialEq, Deserialize, Serialize)]
pub struct Pairs {
    pub n: usize,
    centers: Vec<Point2<f64>>,
    lengths: Vec<f64>,
    nodes: Vec<[usize; 2]>,
    vectors: Vec<Vector2<f64>>,
    neighboring_cells: Vec<[Patch; 2]>,
    cells_normals: Vec<Vector2<f64>>,
    cells_areas: Vec<f64>,
    nodes_normals: Vec<Vector2<f64>>,
    nodes_areas: Vec<f64>,
    on_bnd: Vec<bool>,
}

impl Pairs {
    pub fn centers(&self) -> &[Point2<f64>] {
        &self.centers
    }

    pub fn lengths(&self) -> &[f64] {
        &self.lengths
    }

    pub fn nodes(&self) -> &[[usize; 2]] {
        &self.nodes
    }

    pub fn vectors(&self) -> &[Vector2<f64>] {
        &self.vectors
    }

    pub fn neighboring_cells(&self) -> &[[Patch; 2]] {
        &self.neighboring_cells
    }

    pub fn cells_normals(&self) -> &[Vector2<f64>] {
        &self.cells_normals
    }

    pub fn cells_areas(&self) -> &[f64] {
        &self.cells_areas
    }

    pub fn nodes_normals(&self) -> &[Vector2<f64>] {
        &self.nodes_normals
    }

    pub fn nodes_areas(&self) -> &[f64] {
        &self.nodes_areas
    }

    pub fn on_bnd(&self) -> &[bool] {
        &self.on_bnd
    }
}

impl Pairs {
    pub fn new(
        centers: Vec<Point2<f64>>,
        lengths: Vec<f64>,
        nodes: Vec<[usize; 2]>,
        vectors: Vec<Vector2<f64>>,
        neighboring_cells: Vec<[Patch; 2]>,
        cells_normals: Vec<Vector2<f64>>,
        cells_areas: Vec<f64>,
        nodes_normals: Vec<Vector2<f64>>,
        nodes_areas: Vec<f64>,
        on_bnd: Vec<bool>,
    ) -> Self {
        let n = centers.len();
        assert_eq!(n, lengths.len());
        assert_eq!(n, nodes.len());
        assert_eq!(n, vectors.len());
        assert_eq!(n, neighboring_cells.len());
        assert_eq!(n, cells_normals.len());
        assert_eq!(n, cells_areas.len());
        assert_eq!(n, nodes_normals.len());
        assert_eq!(n, nodes_areas.len());
        assert_eq!(n, on_bnd.len());
        Pairs {
            n,
            centers,
            lengths,
            nodes,
            vectors,
            neighboring_cells,
            cells_normals,
            cells_areas,
            nodes_normals,
            nodes_areas,
            on_bnd,
        }
    }
}
