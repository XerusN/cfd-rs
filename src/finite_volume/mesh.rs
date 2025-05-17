use cfd_rs_mesh::triangle::advancing_front::{self, advancing_front};
use cfd_rs_utils::{boundary::Boundary, control::OutputControl, mesh::{computational_mesh::Computational2DMesh, indices::{ParentIndex, VertexIndex}, Modifiable2DMesh, Parent}};
use nalgebra::Point2;

use super::config::GeometryConfig;

fn square_4_bc() -> Modifiable2DMesh {
    let parents = vec![Parent::Boundary(Boundary(0)), Parent::Boundary(Boundary(1)), Parent::Boundary(Boundary(2)), Parent::Boundary(Boundary(3))];
    let vertices = vec![
        Point2::new(0.0, 0.0),
        Point2::new(1.0, 0.0),
        Point2::new(1.0, 1.0),
        Point2::new(0.0, 1.0),
    ];

    let edge_to_vertices_and_parent = vec![
        (VertexIndex(0), VertexIndex(1), ParentIndex(0)),
        (VertexIndex(1), VertexIndex(2), ParentIndex(1)),
        (VertexIndex(2), VertexIndex(3), ParentIndex(2)),
        (VertexIndex(3), VertexIndex(0), ParentIndex(3)),
    ];

    let mesh;

    unsafe {
        mesh = Modifiable2DMesh::new_from_boundary(vertices, edge_to_vertices_and_parent, parents);
    }

    mesh
}

/// Needs a lot of rework
pub fn mesh(geometry: GeometryConfig) -> Computational2DMesh {
    let mut mesh = square_4_bc();
    advancing_front(&mut mesh, 0.01, OutputControl::Final);
    Computational2DMesh::new_from_he(mesh.0)
}