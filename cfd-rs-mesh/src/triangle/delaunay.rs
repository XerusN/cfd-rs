use cfd_rs_utils::{
    errors::MeshError,
    mesh::{indices::*, *},
};
use nalgebra::Point2;

/// Take as input the points of a triangle (in trigonometric order) and return if the point is in the circumcircle.
pub fn point_in_circumcircle(triangle: [Point2<f64>; 3], point: &Point2<f64>) -> bool {
    let ax_ = triangle[0].x - point.x;
    let ay_ = triangle[0].y - point.y;
    let bx_ = triangle[1].x - point.x;
    let by_ = triangle[1].y - point.y;
    let cx_ = triangle[2].x - point.x;
    let cy_ = triangle[2].y - point.y;

    ((ax_ * ax_ + ay_ * ay_) * (bx_ * cy_ - cx_ * by_)
        - (bx_ * bx_ + by_ * by_) * (ax_ * cy_ - cx_ * ay_)
        + (cx_ * cx_ + cy_ * cy_) * (ax_ * by_ - bx_ * ay_))
        > 0.
}

pub fn delaunay_check_and_swap(
    mesh: &mut Modifiable2DMesh,
    stack: &mut Vec<ParentIndex>,
) -> Result<(), MeshError> {
    let &stack_parent = match stack.last() {
        None => return Ok(()),
        Some(parent) => parent,
    };

    if mesh.0.vertices_from_parent(stack_parent).len() != 3 {
        return Err(MeshError::ParentNotTriangle {
            parent: stack_parent,
        });
    }

    for (parent, _) in mesh.0.neighbors_from_parent(stack_parent) {
        if mesh.0.vertices_from_parent(parent).len() != 3 {
            continue;
        }

        let mut common_he = None;
        'outer: for he_0 in mesh.0.he_from_parent(stack_parent) {
            for he_1 in mesh.0.he_from_parent(parent) {
                if he_0 == mesh.0.he_to_twin()[he_1] {
                    common_he = Some(he_0);
                    break 'outer;
                }
            }
        }

        let common_he = match common_he {
            None => {
                panic!(
                    "No common edge for neighbors : {:?}",
                    MeshError::NoCommonEdge {
                        parent_0: stack_parent,
                        parent_1: parent
                    }
                )
            }
            Some(he) => he,
        };

        let triangle = [
            mesh.0.vertices(mesh.0.vertices_from_he(common_he)[0]),
            mesh.0.vertices(mesh.0.vertices_from_he(common_he)[1]),
            mesh.0
                .vertices(mesh.0.vertices_from_he(mesh.0.he_to_next_he()[common_he])[1]),
        ];
        let point = mesh.0.vertices(
            mesh.0
                .vertices_from_he(mesh.0.he_to_next_he()[mesh.0.he_to_twin()[common_he]])[1],
        );

        if point_in_circumcircle(triangle, &point) {
            stack.push(parent);
            return mesh.swap_edge((stack_parent, parent));
        }
    }

    stack.pop();

    Ok(())
}
