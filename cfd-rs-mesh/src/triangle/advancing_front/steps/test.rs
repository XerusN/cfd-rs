use super::*;
use cfd_rs_utils::mesh::{computational_mesh::BoundaryPatch, Parent};

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

#[test]
fn refine_boundary_test1() {
    let mut mesh = simple_mesh();
    let base_len = mesh.0.he_len();
    refine_boundary(&mut mesh, 0.09);
    assert!(mesh.0.he_len() > base_len * 5);

    let mut mesh = simple_mesh();
    mesh.0
        .export_vtk("./output/test_0.vtk")
        .expect("Error in export");
    refine_boundary(&mut mesh, 0.5);
    mesh.0
        .export_vtk("./output/test_1.vtk")
        .expect("Error in export");
    println!("{:?}", mesh);
}

#[test]
fn ideal_node_test() {
    let mut mesh = simple_mesh();
    let element_size = 0.5;
    refine_boundary(&mut mesh, element_size);
    let mut cell_id = None;
    for (i, parent) in mesh.0.parents().iter().enumerate() {
        if let &Parent::Cell = parent {
            cell_id = Some(i)
        }
    }
    let &base_edge = mesh
        .0
        .he_from_parent(match cell_id {
            None => panic!("No cell in the simple mesh after boundary refinement"),
            Some(i) => ParentIndex(i),
        })
        .get(0)
        .expect("No Half edge in first found cell");

    let node = ideal_node(&mesh, base_edge, element_size);

    match node {
        None => panic!("No ideal node created"),
        Some(node) => {
            if (node.x >= 1.) | (node.x <= 0.) | (node.y >= 1.) | (node.y <= 0.) {
                panic!("Ideal node out of bound")
            }
        }
    }
}
