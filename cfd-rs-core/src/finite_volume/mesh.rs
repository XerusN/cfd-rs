use cfd_rs_mesh::triangle::advancing_front::advancing_front;
use cfd_rs_utils::{
    control::OutputControl,
    errors::MeshError,
    mesh::{
        assembled_mesh::{Mesh, MeshCore}, computational_mesh::{
            manual_meshes::{quad_square, straight_line},
            BoundaryPatch, Computational2DMesh,
        }, indices::{BoundaryPatchIndex, ParentIndex, VertexIndex}, Modifiable2DMesh, Parent
    },
};
use nalgebra::Point2;
use std::io;

use crate::finite_volume::config::MeshingConfig;

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

/// Needs a rework
pub fn mesh(geometry: &GeometryConfig) -> Mesh<Computational2DMesh> {
    match &geometry.import_path {
        None => meshing_algos(&geometry.meshing).expect("Error in meshing"),
        Some(path) => match Mesh::deserialize_file(&path) {
            Err(err) => match err.kind() {
                io::ErrorKind::NotFound => {
                    let mesh = meshing_algos(&geometry.meshing).expect("Error in meshing");
                    mesh.serialize_file(path).unwrap();
                    mesh
                }
                _ => Err(err).unwrap(),
            },
            Ok(mesh) => mesh,
        },
    }
}

fn meshing_algos(meshing_config: &MeshingConfig) -> Result<Mesh<Computational2DMesh>, MeshError> {
    match meshing_config {
        MeshingConfig::AdvancingFront { element_size } => {
            let mut mesh = square_4_bc();
            advancing_front(&mut mesh, *element_size, OutputControl::None)?;
            let mesh = Computational2DMesh::new_from_he(mesh.0);
            Ok(Mesh::from(mesh))
        }
        MeshingConfig::Cartesian { length, n_elements } => Ok(Mesh::from(quad_square(length, n_elements))),
        MeshingConfig::Line { length, n_elements } => Ok(Mesh::from(straight_line(*length, *n_elements))),
    }
}
