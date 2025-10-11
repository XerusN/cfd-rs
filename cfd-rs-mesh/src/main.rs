use cfd_rs_utils::{
    control::*,
    mesh::{
        computational_mesh::{BoundaryPatch, Computational2DMesh},
        indices::*,
        Modifiable2DMesh, *,
    },
};
use nalgebra::Point2;
use std::f64::consts::PI;
use triangle::advancing_front::*;

pub mod triangle;

fn simple_mesh() -> Modifiable2DMesh {
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
            (ParentIndex(0), Some(BoundaryPatchIndex(0))),
        ),
        (
            VertexIndex(2),
            VertexIndex(3),
            (ParentIndex(0), Some(BoundaryPatchIndex(0))),
        ),
        (
            VertexIndex(3),
            VertexIndex(0),
            (ParentIndex(0), Some(BoundaryPatchIndex(0))),
        ),
    ];

    let boundaries = vec![BoundaryPatch::new("Test".to_string())];

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

fn circle_mesh() -> (Modifiable2DMesh, f64) {
    let parents = vec![Parent::Boundary];
    let n = 60;
    let mut vertices = vec![];
    let mut edge_to_vertices_and_parent = vec![];
    for i in 0..n {
        let angle = 2. * PI * i as f64 / n as f64;
        vertices.push(Point2::new(angle.cos(), angle.sin()));
        edge_to_vertices_and_parent.push((
            VertexIndex(i),
            VertexIndex((i + 1) % n),
            (ParentIndex(0), Some(BoundaryPatchIndex(0))),
        ));
    }

    let element_size =
        ((vertices[0].x - vertices[1].x).powi(2) + (vertices[0].y - vertices[1].y).powi(2)).sqrt();

    let mesh;

    let boundaries = vec![BoundaryPatch::new("Test".to_string())];

    unsafe {
        mesh = Modifiable2DMesh::new_from_boundary(
            vertices,
            edge_to_vertices_and_parent,
            parents,
            boundaries,
        );
    }

    (mesh, element_size)
}

fn main() {
    let mut mesh = simple_mesh();
    let element_size = 0.5;
    // let (mut mesh, element_size) = circle_mesh();
    advancing_front(&mut mesh, element_size, OutputControl::Final).unwrap();
    println!("{:?}", mesh);
    let mesh = Computational2DMesh::new_from_he(mesh.0);
}
