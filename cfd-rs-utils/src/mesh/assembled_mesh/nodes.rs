use nalgebra::{Point2, Vector2};
use serde::{Deserialize, Serialize};

use crate::mesh::assembled_mesh::core::Patch;

#[derive(Debug, Clone, PartialEq, Deserialize, Serialize)]
pub struct Nodes {
    pub n: usize,
    centers: Vec<Point2<f64>>,
    volumes: Vec<f64>,
    areas: Vec<Vec<f64>>,
    normals: Vec<Vec<Vector2<f64>>>,
    neighboring_nodes: Vec<Vec<usize>>,
    neighboring_cells: Vec<Vec<Patch>>,
    neighboring_pairs: Vec<Vec<usize>>,
    cv_nodes: Vec<Vec<Point2<f64>>>,
}

impl Nodes {
    pub fn centers(&self) -> &[Point2<f64>] {
        &self.centers
    }
    
    pub fn volumes(&self) -> &[f64] {
        &self.volumes
    }
    
    pub fn areas(&self) -> &[Vec<f64>] {
        &self.areas
    }
    
    pub fn normals(&self) -> &[Vec<Vector2<f64>>] {
        &self.normals
    }
    
    pub fn cv_nodes(&self) -> &[Vec<Point2<f64>>] {
        &self.cv_nodes
    }
    
    pub fn neighboring_nodes(&self) -> &[Vec<usize>] {
        &self.neighboring_nodes
    }
    
    pub fn neighboring_pairs(&self) -> &[Vec<usize>] {
        &self.neighboring_pairs
    }
    
    pub fn neighboring_cells(&self) -> &[Vec<Patch>] {
        &self.neighboring_cells
    }
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
        cv_nodes: Vec<Vec<Point2<f64>>>,
    ) -> Self {
        let n = centers.len();
        assert_eq!(n, volumes.len());
        assert_eq!(n, areas.len());
        assert_eq!(n, normals.len());
        assert_eq!(n, neighboring_nodes.len());
        assert_eq!(n, neighboring_cells.len());
        assert_eq!(n, neighboring_pairs.len());
        assert_eq!(n, cv_nodes.len());
        Nodes {
            n,
            centers,
            volumes,
            areas,
            normals,
            neighboring_nodes,
            neighboring_cells,
            neighboring_pairs,
            cv_nodes,
        }
    }
}
