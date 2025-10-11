use cfd_rs_utils::mesh::{indices::*, Base2DMesh};
use nalgebra::{Point2, Vector2};

#[cfg(test)]
mod test;

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum ConsideredPoint {
    NewPoint(Point2<f64>),
    OldPoint(VertexIndex),
}

impl ConsideredPoint {
    pub fn coordinates(&self, mesh: &Base2DMesh) -> Point2<f64> {
        match self {
            ConsideredPoint::NewPoint(point) => *point,
            ConsideredPoint::OldPoint(point_id) => mesh.vertices(*point_id),
        }
    }
}

/// Space Normalization for isotropic mesh control
#[derive(Clone, Debug, PartialEq)]
pub struct NormalizedSpace {
    pub u: Vector2<f64>,
    pub center: Point2<f64>,
}

impl NormalizedSpace {
    pub fn new(mesh: &Base2DMesh, base_edge: HalfEdgeIndex) -> Self {
        let base_edge_vert = mesh.vertices_from_he(base_edge);
        let base_edge_vert = (
            mesh.vertices(base_edge_vert[0]),
            mesh.vertices(base_edge_vert[1]),
        );
        let center = Point2::new(
            base_edge_vert.0.x + base_edge_vert.1.x,
            base_edge_vert.0.y + base_edge_vert.1.y,
        ) / 2.;
        let u = mesh.he_vector(base_edge).normalize();
        NormalizedSpace { u, center }
    }
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct SpaceNormalized<T>(pub T);

pub trait SpaceNormalize<S> {
    fn to_normalized_space(&self, ns: &NormalizedSpace) -> SpaceNormalized<S>;

    fn from_normalized_space(coordinates: &SpaceNormalized<S>, ns: &NormalizedSpace) -> Self;
}

impl SpaceNormalize<Point2<f64>> for Point2<f64> {
    fn to_normalized_space(&self, ns: &NormalizedSpace) -> SpaceNormalized<Point2<f64>> {
        //center at the mid point
        let centered = Point2::new(self.x - ns.center.x, self.y - ns.center.y);
        // Rotate in the direction of the base edge
        let rotated = Point2::new(
            centered.x * ns.u.x + centered.y * ns.u.y,
            -centered.x * ns.u.y + centered.y * ns.u.x,
        );
        SpaceNormalized(rotated)
    }

    fn from_normalized_space(
        coordinates: &SpaceNormalized<Point2<f64>>,
        ns: &NormalizedSpace,
    ) -> Self {
        // Rotate back
        let derotated = Point2::new(
            coordinates.0.x * ns.u.x - coordinates.0.y * ns.u.y,
            coordinates.0.x * ns.u.y + coordinates.0.y * ns.u.x,
        );
        //center at the mid point
        Point2::new(derotated.x + ns.center.x, derotated.y + ns.center.y)
    }
}

impl SpaceNormalize<Vector2<f64>> for Vector2<f64> {
    fn to_normalized_space(&self, ns: &NormalizedSpace) -> SpaceNormalized<Vector2<f64>> {
        //center at the mid point
        // let centered = Vector2::new(self.x - ns.center.x, self.y - ns.center.y);
        let centered = self;
        // Rotate in the direction of the base edge
        let rotated = Vector2::new(
            centered.x * ns.u.x + centered.y * ns.u.y,
            -centered.x * ns.u.y + centered.y * ns.u.x,
        );
        SpaceNormalized(rotated)
    }

    fn from_normalized_space(
        coordinates: &SpaceNormalized<Vector2<f64>>,
        ns: &NormalizedSpace,
    ) -> Self {
        // Rotate back
        let derotated = Vector2::new(
            coordinates.0.x * ns.u.x - coordinates.0.y * ns.u.y,
            coordinates.0.x * ns.u.y + coordinates.0.y * ns.u.x,
        );
        //center at the mid point
        // Vector2::new(derotated.x + ns.center.x, derotated.y + ns.center.y)
        Vector2::new(derotated.x, derotated.y)
    }
}

/// Checks if a triangle intersects the front
pub fn triangle_intersect_front(
    mesh: &Base2DMesh,
    front_parent: ParentIndex,
    base_edge: HalfEdgeIndex,
    considered_point: ConsideredPoint,
) -> bool {
    let point = considered_point.coordinates(mesh);
    //println!("point: {:?}", point);
    // Might cause issue with ordering
    let egde1 = mesh.vertices_from_he(base_edge);
    let triangle = [mesh.vertices(egde1[0]), mesh.vertices(egde1[1]), point];
    //println!("Triangle: {:?}", triangle);
    // Helps to have some assumptions later (points in trigonometric order)
    if half_plane(&[triangle[0], triangle[1]], &triangle[2]) <= 0. {
        return true;
    }

    'half_edge: for he in mesh.he_from_parent(front_parent) {
        let edge = mesh.vertices_from_he(he);
        let edge = [mesh.vertices(edge[0]), mesh.vertices(edge[1])];
        for i in 0..3 {
            if (edge[0] == triangle[i]) && (edge[1] == triangle[(i + 1) % 3]) {
                continue 'half_edge;
            }
        }
        //println!("Tested edge: {:?}", edge);
        if edge_intersect_triangle(&triangle, &edge) {
            //println!("Intersecting edge!");
            //println!();
            return true;
        }
    }

    false
}

pub fn edge_intersect_triangle(triangle: &[Point2<f64>], edge: &[Point2<f64>]) -> bool {
    // 1. If both point are on the same negative half-edge => no intersection
    for i in 0..3 {
        if half_plane(&[triangle[i], triangle[(i + 1) % 3]], &edge[0]) < 0.
            && half_plane(&[triangle[i], triangle[(i + 1) % 3]], &edge[1]) < 0.
        {
            return false;
        }
    }

    for i in 0..3 {
        for j in 0..2 {
            if edge[j] == triangle[i] {
                if (edge[(j + 1) % 2] == triangle[(i + 1) % 3])
                    | (edge[(j + 1) % 2] == triangle[(i + 2) % 3])
                {
                    return true;
                }
                return !((half_plane(&[triangle[i], triangle[(i + 1) % 3]], &edge[(j + 1) % 2])
                    < 0.)
                    | (half_plane(&[triangle[(i + 2) % 3], triangle[i]], &edge[(j + 1) % 2])
                        < 0.));
            }
        }
    }

    // 2. If one or two of the points is in the triangle => intersection
    if point_in_triangle(triangle, &edge[0]) {
        return true;
    }
    if point_in_triangle(triangle, &edge[1]) {
        return true;
    }

    // 3. Intersection directly (rare)
    let mut t0: f64 = 0.;
    let mut t1: f64 = 1.;
    for i in 0..3 {
        let b0 = half_plane(&[triangle[i], triangle[(i + 1) % 3]], &edge[0]);
        let b1 = half_plane(&[triangle[i], triangle[(i + 1) % 3]], &edge[1]);
        if b0 * b1 < 0. {
            if b0 < 0. {
                t0 = t0.max(b0 / (b0 - b1));
            } else {
                t1 = t1.min(b0 / (b0 - b1));
            }
            if t0 > t1 {
                return false;
            }
        }
    }

    true
}

/// Determines if a point is on the positive or negative half-plane of an edge
pub fn half_plane(edge: &[Point2<f64>], point: &Point2<f64>) -> f64 {
    let normal = Vector2::new(edge[0].y - edge[1].y, edge[1].x - edge[0].x).normalize();
    normal.dot(&Vector2::new(point.x - edge[0].x, point.y - edge[0].y))
    // (edge[1].x - edge[0].x) * (point.y - edge[0].y)
    //     - (edge[1].y - edge[0].y) * (point.x - edge[0].x)
}

/// Tells if a point is in a triangle, the triangle points must be in the trigonometric order.
pub fn point_in_triangle(triangle: &[Point2<f64>], point: &Point2<f64>) -> bool {
    // A point is in a triangle if and only if he is on the 3 positive half-planes
    for i in 0..3 {
        if half_plane(&[triangle[i], triangle[(i + 1) % 3]], point) <= 0. {
            return false;
        }
    }

    true
}
