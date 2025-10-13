use nalgebra::Vector2;
use serde::{Deserialize, Serialize};

pub use cells::Cells;
pub use core::MeshCore;
pub use nodes::Nodes;
pub use pairs::Pairs;

use crate::geometry::area;

mod cells;
mod core;
mod nodes;
mod pairs;

#[derive(Debug, Clone, PartialEq, Deserialize, Serialize)]
#[serde(bound = "T: MeshCore")]
pub struct Mesh<T: MeshCore> {
    core: T,
    pub nodes: Nodes,
    pub cells: Cells,
    pub pairs: Pairs,
}

impl<T: MeshCore> Mesh<T> {}

impl<T: MeshCore> From<T> for Mesh<T> {
    fn from(core: T) -> Self {
        let n_nodes = core.n_nodes();
        let n_cells = core.n_cells();
        // Made for 2D!!
        let n_pairs = core.n_faces();

        let mut nodes_centers = Vec::with_capacity(n_nodes);
        let mut nodes_volumes = Vec::with_capacity(n_nodes);
        let mut nodes_areas = Vec::with_capacity(n_nodes);
        let mut nodes_normals = Vec::with_capacity(n_nodes);
        let mut nodes_neighboring_nodes = Vec::with_capacity(n_nodes);
        let mut nodes_neighboring_cells = Vec::with_capacity(n_nodes);
        let mut nodes_neighboring_pairs = Vec::with_capacity(n_nodes);

        let mut cells_centers = Vec::with_capacity(n_cells);
        let mut cells_volumes = Vec::with_capacity(n_cells);
        let mut cells_areas = Vec::with_capacity(n_cells);
        let mut cells_normals = Vec::with_capacity(n_cells);
        let mut cells_neighboring_nodes = Vec::with_capacity(n_cells);
        let mut cells_neighboring_cells = Vec::with_capacity(n_cells);
        let mut cells_neighboring_pairs = Vec::with_capacity(n_cells);

        let mut pairs_centers = Vec::with_capacity(n_pairs);
        let mut pairs_lengths = Vec::with_capacity(n_pairs);
        let mut pairs_normals = Vec::with_capacity(n_pairs);
        let mut pairs_nodes = Vec::with_capacity(n_pairs);
        let mut pairs_neighboring_cells = Vec::with_capacity(n_pairs);
        
        // ---------------------------
        
        // Connectivity
        for i_node in 0..n_nodes {
            nodes_centers.push(core.node(i_node));
            nodes_neighboring_nodes.push(core.node_to_nodes(i_node));
            // need to check for consistency
            nodes_neighboring_pairs.push(core.node_to_faces(i_node));
            nodes_neighboring_cells.push(core.node_to_cells(i_node));
        }
        for i_cell in 0..n_cells {
            cells_centers.push(core.node(i_cell));
            cells_neighboring_nodes.push(core.node_to_nodes(i_cell));
            // Need to check for consitency face/neighbor
            cells_neighboring_pairs.push(core.node_to_faces(i_cell));
            cells_neighboring_cells.push(core.node_to_cells(i_cell));
        }
        for i_pair in 0..n_pairs {
            pairs_nodes.push(core.face_to_nodes(i_pair));
            pairs_neighboring_cells.push(core.face_to_neighbors(i_pair));
        }
        
        // Geometry
        for i_node in 0..n_nodes {
            let mut cv_nodes = Vec::with_capacity(cells_neighboring_nodes[i_cell].len());
            for i in 0..nodes_neighboring_pairs[i_node].len() {
                cv_nodes.push(pairs_centers[nodes_neighboring_pairs[i_node][i]]);
                cv_nodes.push(cells_centers[nodes_neighboring_cells[i_node][i]]);
            }
            // For 2D only
            let mut areas = Vec::with_capacity(cv_nodes.len()/2);
            let mut normals = Vec::with_capacity(cv_nodes.len()/2);
            for i_node in 0..cv_nodes.len()/2 {
                let vector;
                if i_node == 0 {
                    vector = cv_nodes[i_node] - cv_nodes[cv_nodes.len() - 1];
                } else {
                    vector = cv_nodes[i_node*3] - cv_nodes[i_node*3 - 1];
                }
                let area_1 = vector.magnitude();
                let vector = vector.normalize();
                let normal_1 = Vector2::new(vector.y, - vector.x).normalize();
                let vector = cv_nodes[(i_node*3 + 1) % cv_nodes.len()] - cv_nodes[i_node*3];
                let area_2 = vector.magnitude();
                let vector = vector.normalize();
                let normal_2 = Vector2::new(vector.y, - vector.x).normalize();
                areas.push(area_1 + area_2);
                normals.push(normal_1.lerp(&normal_2, area_2/(area_1 + area_2)));
            }
            cells_areas.push(areas);
            cells_normals.push(normals);
            cells_volumes.push(area(&cv_nodes, &nodes_centers[i_node]));
        }
        
        for i_cell in 0..n_cells {
            let mut nodes = Vec::with_capacity(cells_neighboring_nodes[i_cell].len());
            for i_node in &cells_neighboring_nodes[i_cell] {
                nodes.push(nodes_centers[*i_node]);
            }
            // For 2D only
            let mut areas = Vec::with_capacity(nodes.len());
            let mut normals = Vec::with_capacity(nodes.len());
            for i_node in 0..nodes.len() {
                let vector;
                if i_node == 0 {
                    vector = nodes[i_node] - nodes[nodes.len() - 1];
                } else {
                    vector = nodes[i_node] - nodes[i_node - 1];
                }
                areas.push(vector.magnitude());
                normals.push(Vector2::new(vector.y, - vector.x).normalize());
            }
            cells_areas.push(areas);
            cells_normals.push(normals);
            cells_volumes.push(area(&nodes, &cells_centers[i_cell]));
        }
        for i_pair in 0..n_pairs {
            let i_nodes = pairs_nodes[i_pair];
            let nodes = [core.node(i_nodes[0]), core.node(i_nodes[1])];
            pairs_centers.push(nodes[0].lerp(&nodes[1], 0.5));
            pairs_lengths.push((nodes[1] - nodes[0]).magnitude());
            pairs_normals.push((nodes[1] - nodes[0]).normalize());
        }
        
        // ---------------------------
        
        let nodes = Nodes::new(
            nodes_centers,
            nodes_volumes,
            nodes_areas,
            nodes_normals,
            nodes_neighboring_nodes,
            nodes_neighboring_cells,
            nodes_neighboring_pairs,
        );

        let cells = Cells::new(
            cells_centers,
            cells_volumes,
            cells_areas,
            cells_normals,
            cells_neighboring_nodes,
            cells_neighboring_cells,
            cells_neighboring_pairs,
        );

        let pairs = Pairs::new(
            pairs_centers,
            pairs_lengths,
            pairs_normals,
            pairs_nodes,
            pairs_neighboring_cells,
        );

        Mesh {
            core,
            nodes,
            cells,
            pairs,
        }
    }
}
