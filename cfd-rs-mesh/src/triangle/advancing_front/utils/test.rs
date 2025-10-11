use super::*;
use cfd_rs_utils::mesh::{computational_mesh::BoundaryPatch, Modifiable2DMesh, Parent};

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
// No tests for points yet
fn space_normalization_test() {
    let mesh = simple_mesh();

    for i in 0..mesh.0.he_len() {
        let he = HalfEdgeIndex(i);

        let space = NormalizedSpace::new(&mesh.0, he);

        // Vec tests
        let he_vec = mesh.0.he_vector(he);
        let he_vec_norm = he_vec.to_normalized_space(&space);
        let prev = mesh.0.he_to_prev_he()[he];
        let prev_vec = mesh.0.he_vector(prev);
        let prev_vec_norm = prev_vec.to_normalized_space(&space);
        println!("prev {:?} {:?}", prev, prev_vec);

        println!("prev {:?} {:?}", he_vec_norm, prev_vec_norm);
        // Checks that the angle is conserved (invalid if stretching is implemented)
        if he_vec.angle(&prev_vec) != he_vec_norm.0.angle(&prev_vec_norm.0) {
            panic!(
                "Angle not conserved by normalization: {:?} and normalized {:?}",
                he_vec.angle(&prev_vec),
                he_vec_norm.0.angle(&prev_vec_norm.0)
            );
        };

        let vec = Vector2::new(0.5, 1.8);
        let vec_norm = vec.to_normalized_space(&space);
        // Checks that the angle is conserved (invalid if stretching is implemented)
        println!("vec {:?} {:?}", he_vec_norm, vec_norm);
        if he_vec.angle(&vec) != he_vec_norm.0.angle(&vec_norm.0) {
            panic!(
                "Angle not conserved by normalization: {:?} and normalized {:?}",
                he_vec.angle(&vec),
                he_vec_norm.0.angle(&vec_norm.0)
            );
        };

        // Checks that the vector is the same after normalization and inverse
        if vec != Vector2::from_normalized_space(&vec_norm, &space) {
            panic!(
                "Not equal after normalization and inverse: {:?} and {:?}",
                vec,
                Vector2::from_normalized_space(&vec_norm, &space)
            );
        };

        //Points test
        let he_vertices = mesh.0.vertices_from_he(he);
        let he_vertices = [
            mesh.0.vertices(he_vertices[0]),
            mesh.0.vertices(he_vertices[1]),
        ];
        //Middle point of the base edge
        let point = Point2::new(
            he_vertices[0].x + he_vec.x / 2.,
            he_vertices[0].y + he_vec.y / 2.,
        );
        let point_norm = point.to_normalized_space(&space);
        assert_eq!(point_norm.0, Point2::origin());
        assert_eq!(point, Point2::from_normalized_space(&point_norm, &space));
        //Random point
        let point = Point2::new(0.654, -7.6);
        let point_norm = point.to_normalized_space(&space);
        assert_eq!(point, Point2::from_normalized_space(&point_norm, &space));
    }
}

#[test]
fn test_edge_intersect_triangle() {
    let triangle = vec![
        Point2::new(0.0, 0.0),
        Point2::new(1.0, 0.0),
        Point2::new(0.5, 1.0),
    ];

    // Case 1: Edge intersects the triangle
    let edge_intersecting = vec![Point2::new(0.25, -0.5), Point2::new(0.25, 0.5)];
    assert!(edge_intersect_triangle(&triangle, &edge_intersecting));

    // Case 2: Edge does not intersect the triangle
    let edge_not_intersecting = vec![Point2::new(-1.0, -1.0), Point2::new(-0.5, -0.5)];
    assert!(!edge_intersect_triangle(&triangle, &edge_not_intersecting));

    // Case 3: Edge is entirely inside the triangle
    let edge_inside = vec![Point2::new(0.25, 0.25), Point2::new(0.5, 0.25)];
    assert!(edge_intersect_triangle(&triangle, &edge_inside));

    // Case 4: Edge is on one side of the triangle
    let edge_on_side = vec![Point2::new(0.0, 0.0), Point2::new(1.0, 0.0)];
    assert!(edge_intersect_triangle(&triangle, &edge_on_side));

    // Case 5: Edge is exactly on a vertex of the triangle and colinear
    let edge_on_vertex = vec![Point2::new(0.0, 0.0), Point2::new(0.5, 0.0)];
    assert!(edge_intersect_triangle(&triangle, &edge_on_vertex));

    // Case 6: Edge is exactly on a vertex of the triangle
    let edge_on_vertex = vec![Point2::new(0.0, 0.0), Point2::new(0.5, -1.0)];
    assert!(!edge_intersect_triangle(&triangle, &edge_on_vertex));
}
