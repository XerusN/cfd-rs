use nalgebra::Point2;
use serde::{de::DeserializeOwned, Deserialize, Serialize};
use std::fmt::Debug;

#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub enum Patch {
    Cell(usize),
    Boundary(usize),
}

/// API for 2D for now
pub trait MeshCore: Debug + Clone + PartialEq + DeserializeOwned + Serialize {
    fn n_nodes(&self) -> usize;

    fn n_faces(&self) -> usize;

    fn n_cells(&self) -> usize;

    fn node(&self, i_node: usize) -> Point2<f64>;

    fn face_to_nodes(&self, i_face: usize) -> [usize; 2];

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
    
    fn node_to_nodes(&self, i_node: usize) -> Vec<usize> {
        let mut node_to_nodes = vec![];
        
        for i_face in self.node_to_faces(i_node) {
            let mut adjoint = None;
            let nodes = self.face_to_nodes(i_face);
            
            if nodes[0] == i_node {
                adjoint = Some(nodes[1])
            } else if nodes[1] == i_node {
                adjoint = Some(nodes[0])
            }
            
            node_to_nodes.push(match adjoint {
                None => panic!("Invalid connectivities in MeshCore"),
                Some(neighbor) => neighbor,
            });
        }
        
        node_to_nodes
    }

    fn cell_to_neighbors(&self, i_cell: usize) -> Vec<Patch> {
        let mut cell_to_neighbors = vec![];

        for i_face in self.cell_to_faces(i_cell) {
            let mut adjoint = None;
            let neighbors = self.face_to_neighbors(i_face);

            if Patch::Cell(i_cell) == neighbors[0] {
                adjoint = Some(neighbors[1]);
            } else if Patch::Cell(i_cell) == neighbors[1] {
                adjoint = Some(neighbors[0]);
            }

            cell_to_neighbors.push(match adjoint {
                None => panic!("Invalid connectivities in MeshCore"),
                Some(neighbor) => neighbor,
            });
        }

        cell_to_neighbors
    }

    fn cell_to_nodes(&self, i_cell: usize) -> Vec<usize> {
        let mut cell_to_nodes = vec![];

        for i_face in self.cell_to_faces(i_cell) {
            for i_node in self.face_to_nodes(i_face) {
                if !cell_to_nodes.contains(&i_node) {
                    cell_to_nodes.push(i_node);
                }
            }
        }

        cell_to_nodes
    }
}
