use nalgebra::Vector2;
use serde::{de::DeserializeOwned, Serialize, Deserialize};
use std::fmt::Debug;

#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub enum Patch {
    Cell(usize),
    Boundary(usize),
}

pub trait MeshCore: Debug + Clone + PartialEq + DeserializeOwned + Serialize {
    
    fn n_nodes(&self) -> usize;
    
    fn n_faces(&self) -> usize;
    
    fn n_cells(&self) -> usize;
    
    fn node(&self, i_node: usize) -> Vector2<f64>;
    
    fn face_to_nodes(&self, i_face: usize) -> Vec<usize>;
    
    fn face_to_neighbors(&self, i_face: usize) -> [Patch; 2];
    
    fn cell_to_faces(&self, i_cell: usize) -> Vec<usize>;
    
    fn node_to_faces(&self, i_node: usize) -> Vec<usize> {
        
        let mut node_to_faces = vec![];
        
        for i_face in 0..self.n_faces() {
            let face_to_nodes = self.face_to_nodes(i_face);
            for node in face_to_nodes {
                if node == i_node {
                    node_to_faces.push(i_face);
                }
            }
        }
        
        node_to_faces
    }
    
    fn node_to_cells(&self, i_node: usize) -> Vec<usize> {
        
        let mut node_to_cells = vec![];
        
        let node_to_faces = self.node_to_faces(i_node);
        
        for i_cell in 0..self.n_cells() {
            let cell_to_faces = self.cell_to_faces(i_cell);
            for face in cell_to_faces {
                if node_to_faces.contains(&face) {
                    node_to_cells.push(i_cell);
                }
            }
        }
        
        node_to_cells
    }
    
    fn cell_to_neighbors(&self, i_cell: usize) -> Vec<Patch> {
        
        let mut cell_to_neighbors = vec![];
        
        for i_face in self.cell_to_faces(i_cell) {
            let mut adjoint = None;
            let neighbors = self.face_to_neighbors(i_face);
            
            if Patch::Cell(i_cell) == neighbors[0] {
                adjoint = Some(neighbors[1]);
            }
            if Patch::Cell(i_cell) == neighbors[1] {
                adjoint = Some(neighbors[0]);
            }
            
            cell_to_neighbors.push(match adjoint {
                None => panic!("Invalid connectivities in MeshCore"),
                Some(neighbor) => neighbor,
            });
        }
        
        cell_to_neighbors
    }
    
    
    
}