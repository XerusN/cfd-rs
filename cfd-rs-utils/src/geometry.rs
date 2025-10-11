use nalgebra::{Point2, Vector2};

/// Using Heron's formula
pub fn triangle_area(points: &[Point2<f64>]) -> f64 {
    assert_eq!(points.len(), 3);
    let (a, b, c) = (
        line_length(&[points[0], points[1]]),
        line_length(&[points[1], points[2]]),
        line_length(&[points[2], points[0]]),
    );
    let s = 0.5 * (a + b + c);
    (s * (s - a) * (s - b) * (s - c)).sqrt()
}

pub fn line_length(points: &[Point2<f64>; 2]) -> f64 {
    Vector2::new(points[1].x - points[0].x, points[1].y - points[0].y).norm()
}

/// Rotated counter-clockwise and normalized
pub fn line_normal(points: &[Point2<f64>; 2]) -> Vector2<f64> {
    Vector2::new(points[0].y - points[1].y, points[1].x - points[0].x).normalize()
}

pub fn triangle_centroid(points: &[Point2<f64>]) -> Point2<f64> {
    assert_eq!(points.len(), 3);
    Point2::new(
        (points[0].x + points[1].x + points[2].x) / 3.,
        (points[0].y + points[1].y + points[2].y) / 3.,
    )
}

/// Does not need an ordered set of points
pub fn centroid_and_area(points: &[Point2<f64>]) -> (Point2<f64>, f64) {
    if points.len() < 3 {
        panic!("Can't create cell with less than 3 points");
    } else if points.len() == 3 {
        return (triangle_centroid(points), triangle_area(points));
    }
    
    let geometric_center = geometric_center(points);
    
    let mut centroid = Point2::new(0., 0.);
    let mut area = 0.;
    for i in 0..points.len() {
        let triangle = [points[i], points[(i+1)%points.len()], geometric_center];
        let sub_area = triangle_area(&triangle);
        let sub_centroid = triangle_centroid(&triangle);
        area += sub_area;
        centroid.x += sub_centroid.x*sub_area;
        centroid.y += sub_centroid.y*sub_area;
    }
    
    (centroid/area, area)
}

pub fn geometric_center(points: &[Point2<f64>]) -> Point2<f64> {
    let mut geometric_center = Point2::new(0., 0.);
    for point in points {
        geometric_center.x += point.x;
        geometric_center.y += point.y;
    }
    geometric_center /= points.len() as f64;
    geometric_center
}

pub fn geometric_weighting_factor(
    cell_centroids: &[Point2<f64>; 2],
    face_center: &Point2<f64>,
) -> f64 {
    line_length(&[cell_centroids[1], *face_center])
        / line_length(&[cell_centroids[0], cell_centroids[1]])
}

/// Give a set of point in trigonometric order.
/// Unexpected behaviour if two points are at the same place.
/// Home-made algorithm
pub fn order_points(points: &mut Vec<Point2<f64>>) {
    let geomtric_center = geometric_center(&points);
    
    todo!()
}