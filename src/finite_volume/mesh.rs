use cfd_rs_mesh::triangle::advancing_front::advancing_front;
use cfd_rs_utils::{
    control::OutputControl,
    mesh::{
        computational_mesh::{BoundaryPatch, Computational2DMesh},
        indices::{BoundaryPatchIndex, ParentIndex, VertexIndex},
        Modifiable2DMesh, Parent,
    },
};
use nalgebra::Point2;

use super::config::GeometryConfig;

fn square_4_bc() -> Modifiable2DMesh {
    let parents = vec![Parent::Boundary];
    let vertices = vec![
        Point2::new(0.0, 0.0),
        Point2::new(1.0, 0.0),
        Point2::new(1.0, 1.0),
        Point2::new(0.0, 1.0),
    ];

    let edge_to_vertices_and_parent = vec![
        (
            VertexIndex(0),
            VertexIndex(1),
            (ParentIndex(0), Some(BoundaryPatchIndex(0))),
        ),
        (
            VertexIndex(1),
            VertexIndex(2),
            (ParentIndex(0), Some(BoundaryPatchIndex(1))),
        ),
        (
            VertexIndex(2),
            VertexIndex(3),
            (ParentIndex(0), Some(BoundaryPatchIndex(2))),
        ),
        (
            VertexIndex(3),
            VertexIndex(0),
            (ParentIndex(0), Some(BoundaryPatchIndex(3))),
        ),
    ];

    let boundaries = vec![
        BoundaryPatch::new("Bot".to_string()),
        BoundaryPatch::new("Right".to_string()),
        BoundaryPatch::new("Top".to_string()),
        BoundaryPatch::new("Left".to_string()),
    ];

    let mesh;

    unsafe {
        mesh = Modifiable2DMesh::new_from_boundary(
            vertices,
            edge_to_vertices_and_parent,
            parents,
            boundaries,
        );
    }

    mesh
}

/// Needs a lot of rework
pub fn mesh(geometry: &GeometryConfig) -> Computational2DMesh {
    let mut mesh = square_4_bc();
    advancing_front(&mut mesh, 0.07, OutputControl::None).expect("Error in meshing");
    Computational2DMesh::new_from_he(mesh.0)
}
