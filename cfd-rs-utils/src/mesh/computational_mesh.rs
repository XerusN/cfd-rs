use std::{
    fs::File,
    io::{self, Write},
    path::PathBuf,
    usize,
};

use nalgebra::{point, Point2, Vector2};
use serde::{Deserialize, Serialize};

use crate::geometry::*;

use super::{
    indices::{BoundaryPatchIndex, CellIndex, FaceIndex, HalfEdgeIndex, ParentIndex, VertexIndex},
    Base2DMesh, Parent,
};

pub mod manual_meshes;

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub enum Patch {
    Cell(CellIndex),
    Boundary(BoundaryPatchIndex),
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Face {
    vertices: [VertexIndex; 2],
    area: f64,
    normal: Vector2<f64>,
    patches: (Patch, Patch),
}

impl Face {
    /// Patches needs to be given in an order such that the first parent sees the face defined in the trigonometric order
    pub fn new(
        vertices: [VertexIndex; 2],
        patches: (Patch, Patch),
        vertices_glob: &[Point2<f64>],
    ) -> Self {
        let area = line_length(&[vertices_glob[vertices[0].0], vertices_glob[vertices[1].0]]);
        let normal = line_normal(&[vertices_glob[vertices[0].0], vertices_glob[vertices[1].0]]);

        Face {
            vertices,
            area,
            normal,
            patches,
        }
    }

    pub fn vertices(&self) -> &[VertexIndex; 2] {
        &self.vertices
    }

    /// From the first patch to the second
    pub fn normal(&self) -> &Vector2<f64> {
        &self.normal
    }

    pub fn area(&self) -> f64 {
        self.area
    }

    pub fn patches(&self) -> &(Patch, Patch) {
        &self.patches
    }

    pub fn middle_point(&self, vertices_glob: &[Point2<f64>]) -> Point2<f64> {
        vertices_glob[self.vertices[0].0].lerp(&vertices_glob[self.vertices[1].0], 0.5)
    }

    /// Points outward
    pub fn normal_from_cell(&self, cell: CellIndex) -> Option<Vector2<f64>> {
        if let Patch::Cell(id) = self.patches.0 {
            if id == cell {
                return Some(self.normal().clone());
            }
        }
        if let Patch::Cell(id) = self.patches.1 {
            if id == cell {
                return Some((-self.normal()).clone());
            }
        }
        None
    }

    pub fn geometric_weighting_factor(
        &self,
        vertices_glob: &[Point2<f64>],
        cells_glob: &[Cell],
    ) -> (&Patch, &Patch, f64) {
        let r_f = self.middle_point(vertices_glob);
        let r_a = match self.patches.0 {
            Patch::Cell(id) => cells_glob[id.0].centroid,
            _ => {
                let r_b = match self.patches.1 {
                    Patch::Cell(id) => cells_glob[id.0].centroid,
                    _ => panic!("Face with two boundaries as neighbors"),
                };
                return (&self.patches.0, &self.patches.1, (r_b - r_f).magnitude());
            }
        };
        let r_b = match self.patches.1 {
            Patch::Cell(id) => cells_glob[id.0].centroid,
            _ => return (&self.patches.0, &self.patches.1, (r_a - r_f).magnitude()),
        };

        (
            &self.patches.0,
            &self.patches.1,
            (r_b - r_f).magnitude() / (r_b - r_a).magnitude(),
        )
    }
}

#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct Cell {
    volume: f64,
    centroid: Point2<f64>,
    faces: Vec<FaceIndex>,
    vertices: Vec<VertexIndex>,
}

impl Cell {
    /// faces needs to be given in trigonometric order
    pub fn new(faces: Vec<FaceIndex>, vertices_glob: &[Point2<f64>], faces_glob: &[Face], cell_id: CellIndex) -> Self {
        let mut vertices = vec![];
        for (i, face) in faces
            .iter()
            .map(|f_id| &faces_glob[f_id.0])
            .collect::<Vec<&Face>>().iter().enumerate()
        {
            'push_vert: {
                if i == 0 {
                    match face.patches().0 {
                        Patch::Cell(id) if id == cell_id => {
                            vertices.push(face.vertices()[0]);
                            vertices.push(face.vertices()[1]);
                            break 'push_vert;
                        },
                        _ => (),
                    }
                    match face.patches().1 {
                        Patch::Cell(id) if id == cell_id => {
                            vertices.push(face.vertices()[1]);
                            vertices.push(face.vertices()[0]);
                            break 'push_vert;
                        },
                        _ => (),
                    }
                } else if i == faces.len() - 1 {
                    assert!(vertices.contains(&face.vertices[0]), "Last face contains vertices not defined in other faces");
                    assert!(vertices.contains(&face.vertices[1]), "Last face contains vertices not defined in other faces");
                } else {
                    for (j, vert) in face.vertices().iter().enumerate() {
                        if *vert == vertices[vertices.len()-1] {
                            vertices.push(face.vertices[(j+1)%2]);
                            break;
                        }
                    }
                }
            }
            assert!(vertices.len() > 0, "The first face does not contain the cell as one of the patches (Cell {cell_id})");
        }
        assert_eq!(vertices.len(), faces.len(), "Faces not defined in trigonometric order");
        
        let mut points = vec!();
        for vertex in &vertices {
            points.push(vertices_glob[vertex.0]);
        }
        
        let (centroid, volume) = centroid_and_area(&points);
        
        Cell {
            volume,
            centroid,
            faces,
            vertices,
        }
    }

    pub fn volume(&self) -> f64 {
        self.volume
    }

    pub fn centroid(&self) -> &Point2<f64> {
        &self.centroid
    }

    pub fn faces_id(&self) -> &[FaceIndex] {
        &self.faces
    }

    pub fn faces<'a>(&self, faces_glob: &'a [Face]) -> Vec<&'a Face> {
        self.faces
            .iter()
            .map(|f_id| &faces_glob[f_id.0])
            .collect::<Vec<&Face>>()
    }

    pub fn neighboring_cells_id(&self, cells_glob: &[Cell], faces_glob: &[Face]) -> Vec<CellIndex> {
        let faces = self.faces(faces_glob);

        faces
            .iter()
            .filter_map(|&face| match face.patches().0 {
                Patch::Cell(id) if &cells_glob[id.0] == self => match face.patches().1 {
                    Patch::Cell(id) => Some(id),
                    _ => None,
                },
                Patch::Cell(id) => Some(id),
                _ => None,
            })
            .collect()
    }

    pub fn neighboring_patches<'a>(
        &self,
        cells_glob: &[Cell],
        faces_glob: &'a [Face],
    ) -> Vec<&'a Patch> {
        let faces = self.faces(faces_glob);

        faces
            .iter()
            .map(|&face| match face.patches().0 {
                Patch::Cell(id) if &cells_glob[id.0] == self => &face.patches().1,
                _ => &face.patches.0,
            })
            .collect()
    }

    pub fn neighboring_patches_and_faces<'a>(
        &self,
        cells_glob: &[Cell],
        faces_glob: &'a [Face],
    ) -> Vec<(&'a Patch, &'a Face, FaceIndex)> {
        let faces = self.faces(faces_glob);
        let faces_id = self.faces_id();

        faces
            .iter()
            .enumerate()
            .map(|(id, &face)| {
                (
                    match face.patches().0 {
                        Patch::Cell(id) if &cells_glob[id.0] == self => &face.patches().1,
                        _ => &face.patches.0,
                    },
                    face,
                    faces_id[id],
                )
            })
            .collect()
    }

    pub fn vertices_id(&self) -> &[VertexIndex] {
        &self.vertices
    }

    pub fn vertices<'a>(&self, vertices_glob: &'a [Point2<f64>]) -> Vec<&'a Point2<f64>> {
        self.vertices
            .iter()
            .map(|f_id| &vertices_glob[f_id.0])
            .collect::<Vec<&Point2<f64>>>()
    }
}

#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct BoundaryPatch {
    name: String,
    // Difficult to implement, to be considered
    //faces: Vec<FaceIndex>,
}

impl BoundaryPatch {
    pub fn new(name: String) -> Self {
        BoundaryPatch { name }
    }

    pub fn boundary(&self) -> &str {
        &self.name
    }

    // pub fn faces_id(&self) -> &[FaceIndex] {
    //     &self.faces
    // }

    // pub fn faces<'a>(&self, faces_glob: &'a [Face]) -> Vec<&'a Face> {
    //     self.faces
    //         .iter()
    //         .map(|f_id| &faces_glob[f_id.0])
    //         .collect::<Vec<&Face>>()
    // }
}

/// Contains all the topological and geometric data needed by a mesh
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct Computational2DMesh {
    cells: Vec<Cell>,
    boundaries: Vec<BoundaryPatch>,
    faces: Vec<Face>,
    vertices: Vec<Point2<f64>>,
}

impl Computational2DMesh {
    pub fn num_cells(&self) -> usize {
        self.cells.len()
    }

    pub fn num_faces(&self) -> usize {
        self.faces.len()
    }

    pub fn num_vertices(&self) -> usize {
        self.vertices.len()
    }

    pub fn vertices(&self) -> &[Point2<f64>] {
        &self.vertices
    }

    pub fn cells(&self) -> &[Cell] {
        &self.cells
    }

    pub fn faces(&self) -> &[Face] {
        &self.faces
    }

    pub fn neighboring_cells_id(&self, cell: CellIndex) -> Vec<CellIndex> {
        self.cells[cell].neighboring_cells_id(&self.cells, &self.faces)
    }

    pub fn neighboring_patches(&self, cell: CellIndex) -> Vec<&Patch> {
        self.cells[cell].neighboring_patches(&self.cells, &self.faces)
    }

    pub fn neighboring_patches_and_faces(
        &self,
        cell: CellIndex,
    ) -> Vec<(&Patch, &Face, FaceIndex)> {
        self.cells[cell].neighboring_patches_and_faces(&self.cells, &self.faces)
    }

    pub fn normals_from_cell(&self, cell: CellIndex) -> Vec<Vector2<f64>> {
        self.cells[cell]
            .faces(&self.faces)
            .iter()
            .map(|&face| match face.patches.0 {
                Patch::Cell(id) => {
                    if id == cell {
                        face.normal().clone()
                    } else {
                        -face.normal().clone()
                    }
                }
                _ => -face.normal().clone(),
            })
            .collect()
    }

    pub fn normals_from_cell_with_faces_id(
        &self,
        cell: CellIndex,
    ) -> (&[FaceIndex], Vec<Vector2<f64>>) {
        (self.cells[cell].faces_id(), self.normals_from_cell(cell))
    }

    pub fn normal_vectors_from_cell_with_faces_id(
        &self,
        cell: CellIndex,
    ) -> (&[FaceIndex], Vec<Vector2<f64>>) {
        (self.cells[cell].faces_id(), self.normals_from_cell(cell))
    }

    /// Returns the distance to the boundary if there is one
    pub fn geometric_weighting_factor(&self, face_id: FaceIndex) -> (&Patch, &Patch, f64) {
        self.faces[face_id.0].geometric_weighting_factor(&self.vertices, &self.cells)
    }

    pub fn middle_point_from_face(&self, face: FaceIndex) -> Point2<f64> {
        self.faces[face.0].middle_point(&self.vertices)
    }

    pub fn new_from_he(mesh: Base2DMesh) -> Self {
        let mut vertices = Vec::with_capacity(mesh.vertices_len());
        for i in 0..mesh.vertices_len() {
            vertices.push(mesh.vertices(VertexIndex(i)));
        }

        let mut parent_to_patch = vec![];
        let mut cell_id = 0;

        for parent in mesh.parents() {
            match *parent {
                Parent::None => panic!(),
                Parent::Cell => {
                    parent_to_patch.push(Patch::Cell(CellIndex(cell_id)));
                    cell_id += 1;
                }
                Parent::Boundary => {
                    // value stored here must not be used
                    parent_to_patch.push(Patch::Boundary(BoundaryPatchIndex(usize::MAX)));
                    // boundary_id += 1;
                }
            }
        }

        let mut he_to_face = vec![];
        let mut faces = vec![];

        let mut i = 0;
        for (he, &HalfEdgeIndex(twin)) in mesh.he_to_twin().iter().enumerate() {
            match he.cmp(&twin) {
                std::cmp::Ordering::Less => {
                    he_to_face.push(FaceIndex(i));
                    let patches = (
                        if let Parent::Boundary = mesh.parents()[mesh.he_to_parent()[he].0 .0] {
                            Patch::Boundary(
                                mesh.he_to_parent()[he].1.expect("Boundary has no index"),
                            )
                        } else {
                            parent_to_patch[mesh.he_to_parent()[he].0 .0].clone()
                        },
                        if let Parent::Boundary = mesh.parents()[mesh.he_to_parent()[twin].0 .0] {
                            Patch::Boundary(
                                mesh.he_to_parent()[twin].1.expect("Boundary has no index"),
                            )
                        } else {
                            parent_to_patch[mesh.he_to_parent()[twin].0 .0].clone()
                        },
                    );
                    faces.push(Face::new(
                        mesh.vertices_from_he(HalfEdgeIndex(he)),
                        patches,
                        &vertices,
                    ));
                    i += 1;
                }
                std::cmp::Ordering::Greater => he_to_face.push(he_to_face[twin]),
                std::cmp::Ordering::Equal => panic!(),
            };
        }

        let mut cells = vec![];

        for (parent, patch) in parent_to_patch.iter().enumerate() {
            match *patch {
                Patch::Boundary(_) => {
                    // Patch::Boundary(id) => {
                    // let boundary = match mesh.parents()[parent].clone() {
                    //     Parent::Boundary => boundary,
                    //     _ => panic!(
                    //         "Bad construction of parent_to_patch at index : {:?}",
                    //         parent
                    //     ),
                    // };
                    // let faces_loc = mesh
                    //     .he_from_parent(ParentIndex(parent))
                    //     .iter()
                    //     .map(|he| if mesh.he_to_parent()[he].1 == Some(0) {
                    //         he_to_face[he.0]
                    //     })
                    //     .collect();
                    // if id.0 != boundaries.len() {
                    //     panic!("Wrong construction of boundary");
                    // }
                    // boundaries.push(BoundaryPatch::new(boundary, faces_loc));
                }
                Patch::Cell(id) => {
                    match mesh.parents()[parent] {
                        Parent::Cell => (),
                        _ => panic!(
                            "Bad construction of parent_to_patch at index : {:?}",
                            parent
                        ),
                    };
                    let faces_loc = mesh
                        .he_from_parent(ParentIndex(parent))
                        .iter()
                        .map(|he| he_to_face[he.0])
                        .collect();
                    if id.0 != cells.len() {
                        panic!("Wrong construction of boundary");
                    }
                    cells.push(Cell::new(faces_loc, &vertices, &faces, CellIndex(cells.len())));
                }
            }
        }

        Self {
            cells,
            faces,
            boundaries: mesh.boundaries().clone(),
            vertices,
        }
    }

    pub unsafe fn manual_new(
        cells: Vec<Cell>,
        faces: Vec<Face>,
        boundaries: Vec<BoundaryPatch>,
        vertices: Vec<Point2<f64>>,
    ) -> Self {
        Self {
            cells,
            faces,
            boundaries,
            vertices,
        }
    }

    /// https://docs.vtk.org/en/latest/design_documents/VTKFileFormats.html#unstructuredgrid
    pub fn export(&self, path: String) -> io::Result<()> {
        let path = PathBuf::from(format!("{}/mesh.vtu", &path,));

        let mut file = File::create(&path)?;

        writeln!(
            file,
            "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">"
        )?;
        writeln!(file, "  <UnstructuredGrid>")?;
        writeln!(
            file,
            "    <Piece NumberOfPoints=\"{}\" NumberOfCells=\"{}\">",
            self.num_vertices(),
            self.num_cells()
        )?;
        writeln!(file, "      <Points>")?;
        writeln!(
            file,
            "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">"
        )?;
        write!(file, "          ")?;
        for vertex in self.vertices() {
            write!(file, "{} {} 0 ", vertex.x, vertex.y)?;
        }
        writeln!(file)?;
        writeln!(file, "        </DataArray>")?;
        writeln!(file, "      </Points>")?;

        writeln!(file, "      <Cells>")?;
        writeln!(
            file,
            "        <DataArray type=\"UInt64\" Name=\"connectivity\" format=\"ascii\">"
        )?;
        write!(file, "          ")?;
        for cell in self.cells() {
            for vertex_id in cell.vertices_id() {
                write!(file, "{} ", vertex_id.0)?;
            }
        }
        writeln!(file)?;
        writeln!(file, "        </DataArray>")?;
        writeln!(
            file,
            "        <DataArray type=\"UInt64\" Name=\"offsets\" format=\"ascii\">"
        )?;
        write!(file, "          ")?;
        let mut offset = 0;
        for cell in self.cells() {
            offset += cell.vertices_id().len();
            write!(file, "{} ", offset)?;
        }
        writeln!(file)?;
        writeln!(file, "        </DataArray>")?;
        writeln!(
            file,
            "        <DataArray type=\"UInt64\" Name=\"types\" format=\"ascii\">"
        )?;
        write!(file, "          ")?;
        for cell in self.cells() {
            if cell.vertices_id().len() == 3 {
                write!(file, "5 ")?;
            } else {
                unimplemented!();
            }
        }
        writeln!(file)?;
        writeln!(file, "        </DataArray>")?;
        writeln!(file, "      </Cells>")?;

        writeln!(file, "    </Piece>")?;
        writeln!(file, "  </UnstructuredGrid>")?;
        writeln!(file, "</VTKFile>")?;

        Ok(())
    }

    pub fn serialize_file(&self, path: &str) -> std::io::Result<()> {
        let mut file = File::create(path)?;
        bincode::serde::encode_into_std_write(self, &mut file, bincode::config::standard())
            .unwrap();
        Ok(())
    }

    pub fn deserialize_file(path: &str) -> std::io::Result<Computational2DMesh> {
        let mut file = File::open(path)?;
        Ok(bincode::serde::decode_from_std_read(&mut file, bincode::config::standard()).unwrap())
    }
}
