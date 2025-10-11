use nalgebra::{Point2, Vector2};

use crate::mesh::indices::{BoundaryPatchIndex, CellIndex, FaceIndex, VertexIndex};

use super::{BoundaryPatch, Cell, Computational2DMesh, Face, Patch};

pub fn straight_line(l: f64, n: usize) -> Computational2DMesh {
    let t = l/n as f64;

    let mut vertices = vec![];
    let mut cells = vec![];
    let mut faces = vec![];
    let boundaries = vec![
        BoundaryPatch::new("left".to_string()),
        BoundaryPatch::new("right".to_string()),
        BoundaryPatch::new("empty".to_string()),
    ];

    vertices.push(Point2::new(0., 0.));
    vertices.push(Point2::new(0., t));
    faces.push(Face::new(
        [VertexIndex(0), VertexIndex(1)],
        (
            Patch::Boundary(BoundaryPatchIndex(0)),
            Patch::Cell(CellIndex(0)),
        ),
        &vertices,
    ));

    for i in 0..n {
        vertices.push(Point2::new((i as f64 + 1.) * l / n as f64, 0.));
        vertices.push(Point2::new((i as f64 + 1.) * l / n as f64, t));
        faces.push(Face::new(
            [VertexIndex(i*2), VertexIndex(i*2+2)],
            (
                Patch::Cell(CellIndex(i)),
                Patch::Boundary(BoundaryPatchIndex(2)),
            ),
            &vertices,
        ));
        faces.push(Face::new(
            [VertexIndex(i*2+1), VertexIndex(i*2+3)],
            (
                Patch::Boundary(BoundaryPatchIndex(2)),
                Patch::Cell(CellIndex(i)),
            ),
            &vertices,
        ));
        if i < n-1 {
            faces.push(Face::new(
                [VertexIndex(i*2+2), VertexIndex(i*2+3)],
                (
                    Patch::Cell(CellIndex(i)),
                    Patch::Cell(CellIndex(i+1)),
                ),
                &vertices,
            ));
        } else {
            faces.push(Face::new(
                [VertexIndex(i*2+2), VertexIndex(i*2+3)],
                (
                    Patch::Cell(CellIndex(i)),
                    Patch::Boundary(BoundaryPatchIndex(1)),
                ),
                &vertices,
            ));
        }
        let last_face = faces.len() - 1;
        cells.push(Cell::new(vec![FaceIndex(last_face), FaceIndex(last_face-1), FaceIndex(last_face-2), FaceIndex(last_face-3)], &vertices, &faces, CellIndex(cells.len())));
    }

    Computational2DMesh {
        cells,
        boundaries,
        faces,
        vertices,
    }
}

pub fn quad_square(l: &Vector2<f64>, n: &Vector2<usize>) -> Computational2DMesh {
    
    let delta = Vector2::new(l.x/(n.x-1) as f64, l.y/(n.y-1) as f64);

    let mut vertices = vec![];
    let mut cells = vec![];
    let mut faces = vec![];
    let boundaries = vec![
        BoundaryPatch::new("left".to_string()),
        BoundaryPatch::new("bot".to_string()),
        BoundaryPatch::new("right".to_string()),
        BoundaryPatch::new("top".to_string()),
    ];
    
    // First line
    vertices.push(Point2::new(0., 0.));
    for i in 1..n.x {
        vertices.push(Point2::new(i as f64*delta.x, 0.));
        faces.push(Face::new(
            [VertexIndex(i-1), VertexIndex(i)],
            (
                Patch::Cell(CellIndex(i-1)),
                Patch::Boundary(BoundaryPatchIndex(1)),
            ),
            &vertices,
        ));
    }
    let mut first_row_face = faces.len();
    let mut old_second_row_face = 0;
    for j in 1..n.y {
        vertices.push(Point2::new(0., delta.y*j as f64));
        
        faces.push(Face::new(
            [VertexIndex((j-1)*n.x), VertexIndex(j*n.x)],
            (
                Patch::Boundary(BoundaryPatchIndex(0)),
                Patch::Cell(CellIndex((j-1)*(n.x-1))),
            ),
            &vertices,
        ));
        
        for i in 1..n.x {
            vertices.push(Point2::new(i as f64*delta.x, delta.y*j as f64));
            if j < n.y-1 {
                faces.push(Face::new(
                    [VertexIndex(j*n.x + i - 1), VertexIndex(j*n.x + i)],
                    (
                        Patch::Cell(CellIndex(j*(n.x-1) + i-1)),
                        Patch::Cell(CellIndex((j-1)*(n.x-1) + i-1)),
                    ),
                    &vertices,
                ));
            } else {
                faces.push(Face::new(
                    [VertexIndex(j*n.x + i - 1), VertexIndex(j*n.x + i)],
                    (
                        Patch::Boundary(BoundaryPatchIndex(3)),
                        Patch::Cell(CellIndex((j-1)*(n.x-1) + i-1)),
                    ),
                    &vertices,
                ));
            }
            
            if i < n.x-1 {
                faces.push(Face::new(
                    [VertexIndex(j*n.x + i), VertexIndex((j-1)*n.x + i)],
                    (
                        Patch::Cell(CellIndex((j-1)*(n.x-1) + i)),
                        Patch::Cell(CellIndex((j-1)*(n.x-1) + i-1)),
                    ),
                    &vertices,
                ));
            } else {
                faces.push(Face::new(
                    [VertexIndex(j*n.x + i), VertexIndex((j-1)*n.x + i)],
                    (
                        Patch::Boundary(BoundaryPatchIndex(2)),
                        Patch::Cell(CellIndex((j-1)*(n.x-1) + i-1)),
                    ),
                    &vertices,
                ));
            }
            
            // println!("i {i} j {j} faces {}", faces.len());
            
            if j == 1 {
                cells.push(Cell::new(vec![FaceIndex(first_row_face + 2 + 2*(i-1)), FaceIndex(first_row_face + 1 + 2*(i-1)), FaceIndex(first_row_face + 2*(i-1)), FaceIndex(old_second_row_face + (i-1))], &vertices, &faces, CellIndex(cells.len())));
            } else {
                cells.push(Cell::new(vec![FaceIndex(first_row_face + 2 + 2*(i-1)), FaceIndex(first_row_face + 1 + 2*(i-1)), FaceIndex(first_row_face + 2*(i-1)), FaceIndex(old_second_row_face + 2*(i-1))], &vertices, &faces, CellIndex(cells.len())));
            }
            
        }
        
        old_second_row_face = first_row_face + 1;
        first_row_face = faces.len();
    }
    
    Computational2DMesh {
        cells,
        boundaries,
        faces,
        vertices,
    }
    
    
}
