use core::f64;
use std::f64::consts::PI;

use super::utils::*;
use cfd_rs_utils::{
    errors::MeshError,
    mesh::{indices::*, Base2DMesh, Modifiable2DMesh},
};
use nalgebra::{Point2, Vector2};

#[cfg(test)]
mod test;

/// Refines the boundary up to the level specified by the element size.
pub fn refine_boundary(mesh: &mut Modifiable2DMesh, mesh_size: f64) {
    let mesh_size_squared = mesh_size * mesh_size;
    let mut i = 0;
    while i < mesh.0.he_len() {
        while mesh.0.he_vector(HalfEdgeIndex(i)).norm_squared() > mesh_size_squared {
            mesh.split_edge(HalfEdgeIndex(i), 0.5)
                .expect("Error when refining boundary");
        }
        i += 1;
    }
}

/// Selects the edge to be used.
/// For now most simple choice (only size based).
pub fn select_base_edge(mesh: &Modifiable2DMesh, front: &[ParentIndex]) -> HalfEdgeIndex {
    let edges_sizes_squared = mesh
        .0
        .he_from_parent(front[0])
        .iter()
        .map(|he| (*he, mesh.0.he_vector(*he).norm_squared()))
        .collect::<Vec<(HalfEdgeIndex, f64)>>();
    let mut current_min = 0;
    for (i, value) in edges_sizes_squared.iter().enumerate() {
        if value.1 < edges_sizes_squared[current_min].1 {
            current_min = i
        }
    }
    edges_sizes_squared[current_min].0
}

/// Will get proper definition later, for now the mesh size is considered to be constant so the base element size is constant too.
pub fn element_size(mesh_size: f64) -> f64 {
    mesh_size
}

/// Returns the newly created triangles
pub fn new_element(
    mesh: &mut Modifiable2DMesh,
    front: &mut Vec<ParentIndex>,
    base_edge: HalfEdgeIndex,
    element_size: f64,
) -> Result<Vec<ParentIndex>, MeshError> {
    let ideal_nod;
    let mut valid_points = vec![];
    let front_parent = mesh.0.he_to_parent()[base_edge].0;
    'ideal_node: {
        ideal_nod = ideal_node(mesh, base_edge, element_size);
        // println!("{:?}", ideal_nod);
        if ideal_nod.is_none() {
            break 'ideal_node;
        }

        let considered_point = ConsideredPoint::NewPoint(ideal_nod.expect("uh"));

        if !node_validity_check(
            mesh,
            front_parent,
            element_size,
            considered_point,
            base_edge,
        ) {
            // println!("Node not valid");
            break 'ideal_node;
        }
        if !element_validity_check(mesh, front_parent, base_edge, considered_point) {
            //println!("Element not valid");
            break 'ideal_node;
        }
        let tan_a = node_suitability_check(mesh, base_edge, considered_point);
        if tan_a <= 0. {
            break 'ideal_node;
        }
        if tan_a < f64::tan(30.) {
            valid_points.push((considered_point, tan_a));
            break 'ideal_node;
        }

        return add_element(mesh, front, base_edge, considered_point);
    }

    // Implement a choice based on the suitability criterion to enhance quality
    let existing_candidates = find_existing_candidates(mesh, front_parent, base_edge, element_size);

    for point in existing_candidates {
        // if let ConsideredPoint::OldPoint(id) = point {
        //     println!("{:?} {:?}", point, mesh.0.vertices(id));
        // }
        'point: {
            if !node_validity_check(mesh, front_parent, element_size, point, base_edge) {
                // println!("Node not valid");
                break 'point;
            }
            if !element_validity_check(mesh, front_parent, base_edge, point) {
                //println!("Element not valid");
                break 'point;
            }
            let tan_a = node_suitability_check(mesh, base_edge, point);
            if tan_a <= 0. {
                break 'point;
            }
            if tan_a < f64::tan(30.) {
                valid_points.push((point, tan_a));
                break 'point;
            }

            return add_element(mesh, front, base_edge, point);
        }
    }

    let trie_vertices = create_trie_vertices(mesh, base_edge, element_size);
    for point in trie_vertices {
        'point: {
            // if let ConsideredPoint::NewPoint(coord) = point {
            //     println!("{:?} {:?}", point, coord);
            // }
            if !node_validity_check(mesh, front_parent, element_size, point, base_edge) {
                // println!("Node not valid");
                break 'point;
            }
            if !element_validity_check(mesh, front_parent, base_edge, point) {
                // println!("Element not valid");
                break 'point;
            }
            let tan_a = node_suitability_check(mesh, base_edge, point);
            if tan_a <= 0. {
                break 'point;
            }
            if tan_a < f64::tan(30.) {
                valid_points.push((point, tan_a));
                break 'point;
            }

            return add_element(mesh, front, base_edge, point);
        }
    }

    if valid_points.is_empty() {
        Err(MeshError::NoElementCreatable(base_edge))
    } else {
        let mut max = (0, valid_points[0].1);
        for (i, valid) in valid_points.iter().enumerate() {
            if valid.1 > max.1 {
                max = (i, valid.1);
            }
        }
        add_element(mesh, front, base_edge, valid_points[max.0].0)
    }
}

/// Based on J. Frysketig, 1994.
/// Returns a value in the normalized space.
pub fn ideal_node(
    mesh: &Modifiable2DMesh,
    base_edge_id: HalfEdgeIndex,
    element_size: f64,
) -> Option<Point2<f64>> {
    let space_normalization = NormalizedSpace::new(&mesh.0, base_edge_id);
    let base_edge = mesh
        .0
        .he_vector(base_edge_id)
        .to_normalized_space(&space_normalization);
    let prev_edge = (
        mesh.0.he_to_prev_he()[base_edge_id],
        mesh.0
            .he_vector(mesh.0.he_to_prev_he()[base_edge_id])
            .to_normalized_space(&space_normalization),
    );
    let next_edge = (
        mesh.0.he_to_next_he()[base_edge_id],
        mesh.0
            .he_vector(mesh.0.he_to_next_he()[base_edge_id])
            .to_normalized_space(&space_normalization),
    );

    let mut alpha = (-prev_edge.1 .0).angle(&base_edge.0).abs();
    let mut beta = (-next_edge.1 .0).angle(&base_edge.0).abs();
    //println!("alpha, beta : {:?} {:?} {:?}", alpha, beta, PI);

    if alpha > beta {
        std::mem::swap(&mut alpha, &mut beta);
    }

    if alpha < PI * 80. / 180. {
        return None;
    }
    if (alpha < PI * 91. / 180.) && (alpha > PI * 89. / 180.) {
        let node = SpaceNormalized(Point2::new(
            element_size * (PI / 4.).cos(),
            element_size * (PI / 4.).sin(),
        ));
        //println!("ideal node normalized {:?}", node);
        //println!("ideal node {:?}", Point2::from_normalized_space(&node, &space_normalization));
        return Some(Point2::from_normalized_space(&node, &space_normalization));
    }
    let alpha_deg = alpha * 180. / PI;
    let beta_deg = beta * 180. / PI;

    let phi_a = alpha_deg / (alpha_deg / 60.).round() * PI / 180.;
    let phi_b = beta_deg / (beta_deg / 60.).round() * PI / 180.;
    let node = SpaceNormalized(
        Point2::new(phi_a.cos() - phi_b.cos(), phi_a.sin() - phi_b.sin()) * element_size / 2.,
    );
    //println!("ideal node normalized {:?}", node);
    //println!("ideal node {:?}", Point2::from_normalized_space(&node, &space_normalization));
    Some(Point2::from_normalized_space(&node, &space_normalization))
}

pub fn node_validity_check(
    mesh: &Modifiable2DMesh,
    front_parent: ParentIndex,
    element_size: f64,
    considered_point: ConsideredPoint,
    base_edge: HalfEdgeIndex,
) -> bool {
    let proximity_admissibility = 0.67;
    if let ConsideredPoint::OldPoint(old_point) = considered_point {
        if (mesh
            .0
            .vertices_from_he(mesh.0.he_to_next_he()[mesh.0.he_to_twin()[base_edge]])[1]
            == old_point)
            & (mesh
                .0
                .vertices_from_he(mesh.0.he_to_prev_he()[mesh.0.he_to_twin()[base_edge]])[0]
                == old_point)
        {
            return false;
        }
        return true;
    }
    let base_point = considered_point.coordinates(&mesh.0);

    let min = proximity_admissibility * proximity_admissibility * element_size * element_size;

    for node in mesh.0.vertices_from_parent(front_parent) {
        let current_point = mesh.0.vertices(node);
        if nalgebra::distance_squared(&base_point, &current_point) < min {
            return false;
        }
    }

    true
}

pub fn element_validity_check(
    mesh: &Modifiable2DMesh,
    front_parent: ParentIndex,
    base_edge: HalfEdgeIndex,
    considered_point: ConsideredPoint,
) -> bool {
    !triangle_intersect_front(&mesh.0, front_parent, base_edge, considered_point)
}

pub fn node_suitability_check(
    mesh: &Modifiable2DMesh,
    // front: &[ParentIndex],
    base_edge: HalfEdgeIndex,
    considered_point: ConsideredPoint,
) -> f64 {
    let base_edge_vert = mesh.0.vertices_from_he(base_edge);
    let base_edge_vert = (
        mesh.0.vertices(base_edge_vert[0]),
        mesh.0.vertices(base_edge_vert[1]),
    );
    let point = considered_point.coordinates(&mesh.0);

    let mut tan_a = tan_alpha_suitability(&[point, base_edge_vert.0, base_edge_vert.1]);
    tan_a = tan_a.min(tan_alpha_suitability(&[
        base_edge_vert.0,
        base_edge_vert.1,
        point,
    ]));
    tan_a.min(tan_alpha_suitability(&[
        base_edge_vert.1,
        point,
        base_edge_vert.0,
    ]))
}

pub fn tan_alpha_suitability(points: &[Point2<f64>]) -> f64 {
    let point = points[2];
    let base_edge_vec = Vector2::new(points[1].x - points[0].x, points[1].y - points[0].y);

    let t = (base_edge_vec.dot(&Vector2::new(point.x - points[0].x, point.y - points[0].y))
        / base_edge_vec.norm_squared())
    .abs();
    let theta = (t).max(t - 1.);
    Vector2::new(
        point.x - (points[0].x + t * base_edge_vec.x),
        point.y - (points[0].y + t * base_edge_vec.y),
    )
    .norm()
        / (theta * base_edge_vec.norm())
}

pub fn find_existing_candidates(
    mesh: &Modifiable2DMesh,
    front_parent: ParentIndex,
    base_edge: HalfEdgeIndex,
    element_size: f64,
) -> Vec<ConsideredPoint> {
    let rel_radius = 1.33;
    let base_edge_vert = mesh.0.vertices_from_he(base_edge);
    let mid_base_edge = Point2::new(
        mesh.0.vertices(base_edge_vert[0]).x + mesh.0.vertices(base_edge_vert[1]).x,
        mesh.0.vertices(base_edge_vert[0]).y + mesh.0.vertices(base_edge_vert[1]).y,
    ) / 2.;

    let mut candidates = vec![];
    for point in mesh.0.vertices_from_parent(front_parent) {
        if (base_edge_vert[0] != point) & (base_edge_vert[1] != point) {
            let vert = mesh.0.vertices(point);
            // Maybe *4 too small (*4*1.33*1.33 in the book but not same lhs)
            if Vector2::new(vert.x - mid_base_edge.x, vert.y - mid_base_edge.y).norm_squared()
                < element_size * element_size * rel_radius * rel_radius
            {
                candidates.push(ConsideredPoint::OldPoint(point))
            }
        }
    }

    candidates
}

pub fn create_trie_vertices(
    mesh: &Modifiable2DMesh,
    base_edge: HalfEdgeIndex,
    element_size: f64,
) -> Vec<ConsideredPoint> {
    let space_normalization = NormalizedSpace::new(&mesh.0, base_edge);
    let base_edge_vec = mesh
        .0
        .he_vector(base_edge)
        .to_normalized_space(&space_normalization);

    let trie_norm_space = [
        SpaceNormalized(Point2::new(0., element_size)),
        SpaceNormalized(Point2::new(0., element_size * 0.5)),
        SpaceNormalized(Point2::new(0., element_size * 0.5 + 1. / 6.)),
        SpaceNormalized(Point2::new(0., element_size * 0.5 + 1. / 3.)),
        SpaceNormalized(Point2::new(base_edge_vec.0.x / 3., element_size * 0.6)),
        SpaceNormalized(Point2::new(base_edge_vec.0.x / 3., element_size * 0.9)),
        SpaceNormalized(Point2::new(-base_edge_vec.0.x / 3., element_size * 0.6)),
        SpaceNormalized(Point2::new(-base_edge_vec.0.x / 3., element_size * 0.9)),
    ];

    trie_norm_space
        .iter()
        .map(|vert_norm| {
            ConsideredPoint::NewPoint(Point2::from_normalized_space(
                vert_norm,
                &space_normalization,
            ))
        })
        .collect()
}

/// Once the quality and validity is checked modifies the data structure to add the new element.
/// Returns the newly created triangles
///
/// # Warning
///
/// Since the mesh structure is changed any value held before that might represent something else now
pub fn add_element(
    mesh: &mut Modifiable2DMesh,
    front: &mut Vec<ParentIndex>,
    base_edge: HalfEdgeIndex,
    considered_point: ConsideredPoint,
) -> Result<Vec<ParentIndex>, MeshError> {
    //println!("considered point: {:?} {:?}", considered_point, considered_point.coordinates(&mesh.0));
    let point_id = match considered_point {
        ConsideredPoint::NewPoint(point) => {
            let new_parent;
            unsafe { new_parent = mesh.notching(base_edge, point)? };
            // If useless parent added it will get removed
            front.push(new_parent);
            return Ok(clean_front(&mesh.0, front));
        }
        ConsideredPoint::OldPoint(point_id) => point_id,
    };

    // Double check this part
    if mesh.0.vertices_from_he(mesh.0.he_to_next_he()[base_edge])[1] == point_id {
        let front_parent = mesh.0.he_to_parent()[base_edge].0;
        let next_point = mesh.0.vertices_from_he(mesh.0.he_to_next_he()[base_edge])[1];
        let new_parent;
        unsafe {
            new_parent = mesh.trimming(
                (mesh.0.vertices_from_he(base_edge)[0], next_point),
                front_parent,
                None,
            )?;
        }
        front.push(new_parent);
        return Ok(clean_front(&mesh.0, front));
    }
    if mesh.0.vertices_from_he(mesh.0.he_to_prev_he()[base_edge])[0] == point_id {
        let front_parent = mesh.0.he_to_parent()[base_edge].0;
        let next_point = mesh.0.vertices_from_he(mesh.0.he_to_prev_he()[base_edge])[0];
        let new_parent;
        unsafe {
            new_parent = mesh.trimming(
                (mesh.0.vertices_from_he(base_edge)[1], next_point),
                front_parent,
                None,
            )?;
        }
        front.push(new_parent);
        return Ok(clean_front(&mesh.0, front));
    }

    // Only remains the case where both edge have to be created toward an old point, right?
    let front_parent = mesh.0.he_to_parent()[base_edge].0;
    let base_points = mesh.0.vertices_from_he(base_edge);

    // Relies on non-enforced behaviour of the trimming function! (old parent = parent from new he from vertices.0 to vertices.1)
    let new_parent1;
    let new_parent2;
    unsafe {
        new_parent1 = mesh.trimming((point_id, base_points[0]), front_parent, None)?;
    }
    let front_parent = mesh.0.he_to_parent()[mesh.0.he_to_next_he()[base_edge]].0;
    unsafe {
        new_parent2 = mesh.trimming((point_id, base_points[1]), front_parent, None)?;
    }
    front.push(new_parent1);
    front.push(new_parent2);

    Ok(clean_front(&mesh.0, front))
}

pub fn clean_front(mesh: &Base2DMesh, front: &mut Vec<ParentIndex>) -> Vec<ParentIndex> {
    let stack = front
        .clone()
        .into_iter()
        .filter(|&parent| (mesh.vertices_from_parent(parent).len() <= 3))
        .collect();
    front.retain(|&parent| mesh.vertices_from_parent(parent).len() > 3);
    stack
}
