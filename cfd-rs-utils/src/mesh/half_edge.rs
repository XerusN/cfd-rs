use crate::errors::MeshError;
use indices::*;
use nalgebra::{Point2, Vector2};

use std::fmt::{self, Debug, Error};
use std::fs::File;
use std::io::{self, Write};

use super::computational_mesh::BoundaryPatch;

pub mod indices;

#[cfg(test)]
mod test;

/// Parent of a half_edge, either a cell or a boundary.
/// ```None``` is meant as an error or temporary state.
#[derive(Clone, Debug, Default, PartialEq)]
pub enum Parent {
    #[default]
    None,
    Cell,
    Boundary,
}

/// Array based Half-edge data-structure mesh representation
/// Supports meshes of up to a billion element.
/// Since the crate is built for cfd on a classic computer (not HPC) it is easily enough.
#[derive(Clone, Default, PartialEq)]
pub struct Base2DMesh {
    //```he``` is for Half-edge
    he_to_vertex: Vec<VertexIndex>,
    he_to_twin: Vec<HalfEdgeIndex>,
    he_to_next_he: Vec<HalfEdgeIndex>,
    he_to_prev_he: Vec<HalfEdgeIndex>,
    he_to_parent: Vec<(ParentIndex, Option<BoundaryPatchIndex>)>,

    vertices: Vec<Point2<f64>>,
    parents: Vec<Parent>,
    boundaries: Vec<BoundaryPatch>,

    parent_to_first_he: Vec<HalfEdgeIndex>,
}

impl Base2DMesh {
    /// Gets all the vertices of an half-edge.
    pub fn vertices_from_he(&self, he_id: HalfEdgeIndex) -> [VertexIndex; 2] {
        [
            self.he_to_vertex[he_id],
            self.he_to_vertex[self.he_to_twin[he_id]],
        ]
    }

    pub fn vertices(&self, v_id: VertexIndex) -> Point2<f64> {
        self.vertices[v_id]
    }

    /// Returns the number of vertices.
    pub fn vertices_len(&self) -> usize {
        self.vertices.len()
    }

    /// Returns the number of half-edges.
    pub fn he_len(&self) -> usize {
        self.he_to_vertex.len()
    }

    /// Returns the number of parents.
    pub fn parents_len(&self) -> usize {
        self.parents.len()
    }

    /// Gets all half-edges from a parent (Cell or boundary).
    pub fn he_from_parent(&self, parent_id: ParentIndex) -> Vec<HalfEdgeIndex> {
        let first_he = self.parent_to_first_he[parent_id];
        let mut result = vec![first_he];
        let mut current_he = self.he_to_next_he[first_he];

        while first_he != current_he {
            //println!("{:?}, {:?}", current_he, first_he);
            result.push(current_he);
            current_he = self.he_to_next_he[current_he];
        }

        result
    }

    /// Gets all vertices from a parent (Cell or boundary).
    pub fn vertices_from_parent(&self, parent_id: ParentIndex) -> Vec<VertexIndex> {
        self.he_from_parent(parent_id)
            .into_iter()
            .map(|he_id| self.he_to_vertex[he_id])
            .collect()
    }

    /// Gets the parent from an HalfEdge.
    pub fn he_to_parent(&self) -> &Vec<(ParentIndex, Option<BoundaryPatchIndex>)> {
        &self.he_to_parent
    }

    /// Gets the twin HalfEdge from an HalfEdge.
    pub fn he_to_twin(&self) -> &Vec<HalfEdgeIndex> {
        &self.he_to_twin
    }

    /// Gets the next HalfEdge from an HalfEdge.
    pub fn he_to_next_he(&self) -> &Vec<HalfEdgeIndex> {
        &self.he_to_next_he
    }

    /// Gets the previous HalfEdge from an HalfEdge.
    pub fn he_to_prev_he(&self) -> &Vec<HalfEdgeIndex> {
        &self.he_to_prev_he
    }

    /// Gets the parents adjacent to another.
    /// This may have strange behaviours when used on a boundary.
    pub fn neighbors_from_parent(
        &self,
        parent_id: ParentIndex,
    ) -> Vec<(ParentIndex, Option<BoundaryPatchIndex>)> {
        self.he_from_parent(parent_id)
            .into_iter()
            .map(|he_id| self.he_to_parent[self.he_to_twin[he_id]])
            .collect()
    }

    /// Gets the parent properties from its index.
    pub fn parents(&self) -> &Vec<Parent> {
        &self.parents
    }

    pub fn boundaries(&self) -> &Vec<BoundaryPatch> {
        &self.boundaries
    }

    /// Gets the half-edges connected to a vertex
    pub fn he_from_vertex(&self, vertex_id: VertexIndex) -> Vec<HalfEdgeIndex> {
        let mut result = Vec::new();
        for i in 0..self.he_len() {
            if self.he_to_vertex[HalfEdgeIndex(i)] == vertex_id {
                result.push(HalfEdgeIndex(i))
            }
        }
        result
    }

    /// Get the normal (unit vector) from an edge, oriented toward the exterior of the triangle.
    pub fn normal(&self, he_id: HalfEdgeIndex) -> Vector2<f64> {
        let vertices = self.vertices_from_he(he_id);
        let vertices = [self.vertices[vertices[0]], self.vertices[vertices[1]]];
        let segment = Vector2::new(vertices[1].x - vertices[0].x, vertices[1].y - vertices[0].y);
        Vector2::new(-segment.y, segment.x).normalize()
    }

    /// Returns the HalfEdge as a vector
    pub fn he_vector(&self, he_id: HalfEdgeIndex) -> Vector2<f64> {
        let vertices = self.vertices_from_he(he_id);
        let vertices = [self.vertices[vertices[0]], self.vertices[vertices[1]]];
        Vector2::new(vertices[1].x - vertices[0].x, vertices[1].y - vertices[0].y)
    }

    /// Check that the mesh topology is valid.
    /// Used to confirm the topology before switching to an immutable mesh and for test purpose.
    ///
    /// This function aims to be exhaustive but obviously it is not.
    /// For now it mostly checks coherence between arrays and that not indices are out of bound.
    ///
    /// If you find some cases not covered by this function, don't hesitate to submit a pull request or detail what you would like implemented.
    pub fn check_mesh(&self) -> Result<(), MeshError> {
        //Checks length coherence between HalfEdges arrays
        assert_eq!(self.he_to_next_he.len(), self.he_to_twin.len());
        assert_eq!(self.he_to_prev_he.len(), self.he_to_twin.len());
        assert_eq!(self.he_to_vertex.len(), self.he_to_twin.len());
        assert_eq!(self.he_to_parent.len(), self.he_to_twin.len());

        for vertex in &self.he_to_vertex {
            if *vertex >= VertexIndex(self.vertices.len()) {
                return Err(MeshError::VertexIndexOutOfBound {
                    got: *vertex,
                    len: self.vertices.len(),
                });
            }
        }

        for he in &self.he_to_twin {
            if *he >= HalfEdgeIndex(self.he_to_vertex.len()) {
                return Err(MeshError::HalfEdgeIndexOutOfBound {
                    got: *he,
                    len: self.he_to_vertex.len(),
                });
            }
        }

        for (i, he) in self.he_to_twin.iter().enumerate() {
            if self.he_to_twin[*he] != HalfEdgeIndex(i) {
                return Err(MeshError::TwinNotCorrect {
                    he: HalfEdgeIndex(i),
                    he_twin: *he,
                    he_twin_twin: self.he_to_twin[*he],
                });
            }
        }

        for he in &self.he_to_next_he {
            if *he >= HalfEdgeIndex(self.he_to_vertex.len()) {
                return Err(MeshError::HalfEdgeIndexOutOfBound {
                    got: *he,
                    len: self.he_to_vertex.len(),
                });
            }
        }

        for he in &self.he_to_prev_he {
            if *he >= HalfEdgeIndex(self.he_to_vertex.len()) {
                return Err(MeshError::HalfEdgeIndexOutOfBound {
                    got: *he,
                    len: self.he_to_vertex.len(),
                });
            }
        }

        // Simple check from prev and next
        for (i, next) in self.he_to_next_he.iter().enumerate() {
            let he = HalfEdgeIndex(i);
            if he != self.he_to_prev_he[*next] {
                return Err(MeshError::NextPrevNotCorrect {
                    he,
                    he_next: *next,
                    he_next_prev: self.he_to_prev_he[*next],
                });
            }
        }

        // Might be redundant with the previous check
        for i in 0..self.he_to_vertex.len() {
            let origin = HalfEdgeIndex(i);
            let mut j = 0;
            let mut current_he = self.he_to_next_he[origin];
            while current_he != origin {
                current_he = self.he_to_next_he[current_he];
                j += 1;
                // len/2 could work
                if j >= self.he_to_vertex.len() {
                    break;
                }
            }
        }

        for (parent, _) in &self.he_to_parent {
            if *parent >= ParentIndex(self.parents.len()) {
                return Err(MeshError::ParentIndexOutOfBound {
                    got: *parent,
                    len: self.parents.len(),
                });
            }
        }

        for i in 0..self.parents.len() {
            let parent = ParentIndex(i);
            for he in self.he_from_parent(parent) {
                if self.he_to_parent[he].0 != parent {
                    return Err(MeshError::ParentNotCorrect {
                        parent,
                        he,
                        he_parent: self.he_to_parent[he].0,
                    });
                }
            }
        }

        Ok(())
    }

    /// Written by chatGPT, proper export function will be made later
    pub fn export_vtk(&self, filename: &str) -> io::Result<()> {
        let mut file = File::create(filename)?;

        // Write VTK header
        writeln!(file, "# vtk DataFile Version 3.0")?;
        writeln!(file, "2D Mesh Example")?;
        writeln!(file, "ASCII")?;
        writeln!(file, "DATASET POLYDATA")?;

        // Write points
        writeln!(file, "POINTS {} float", self.vertices.len())?;
        for vertex in &self.vertices {
            writeln!(file, "{} {} 0.0", vertex.x, vertex.y)?;
        }

        // Write lines (edges)
        let num_edges = self.he_to_next_he.len();
        writeln!(file, "LINES {} {}", num_edges, num_edges * 3)?;
        for (he, &twin_he) in self.he_to_twin.iter().enumerate() {
            let start_vertex = self.he_to_vertex[he];
            let end_vertex = self.he_to_vertex[twin_he];
            writeln!(file, "2 {} {}", start_vertex, end_vertex)?;
        }

        // Write additional attributes (e.g., he_to_parent)
        writeln!(file, "CELL_DATA {}", num_edges)?;
        writeln!(file, "SCALARS he_to_parent int 1")?;
        writeln!(file, "LOOKUP_TABLE default")?;
        for &parent in &self.he_to_parent {
            writeln!(file, "{}", parent.0)?;
        }

        Ok(())
    }
}

impl Debug for Base2DMesh {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> Result<(), Error> {
        writeln!(f, "Base2DMesh:")?;
        writeln!(f, "------------------------------------")?;
        writeln!(f, "vertex-id | coords [x, y]")?;
        for (id, vertex) in self.vertices.iter().enumerate() {
            writeln!(f, "{id:9} | {vertex:?}")?;
        }
        writeln!(f, "------------------------------------")?;
        writeln!(
            f,
            "he-id | vertex-id | twin-he | parent-id | next-he | prev-he | boundary-id"
        )?;
        for (he, twin) in self.he_to_twin().iter().enumerate() {
            let boundary = match self.he_to_parent[he].1 {
                None => None,
                Some(value) => Some(value.0),
            };
            writeln!(
                f,
                "{he:5} | {:9} | {:7} | {:9} | {:7} | {:7} | {:?}",
                self.he_to_vertex[he].0,
                twin.0,
                self.he_to_parent[he].0 .0,
                self.he_to_next_he()[he].0,
                self.he_to_prev_he()[he].0,
                boundary
            )?;
        }
        writeln!(f, "------------------------------------")?;
        writeln!(f, "parent-id | parent")?;
        for (id, parent) in self.parents.iter().enumerate() {
            writeln!(f, "{id:9} | {parent:?}")?;
        }
        writeln!(f, "------------------------------------")?;
        writeln!(f, "boundary-id | boundaries")?;
        for (id, boundary) in self.boundaries().iter().enumerate() {
            writeln!(f, "{id:11} | {boundary:?}")?;
        }
        Ok(())
    }
}

/// Gives access to modifications from Base2DMesh
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Modifiable2DMesh(pub Base2DMesh);

/// Mesh with valid topology, can be safely used in computations
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Safe2DMesh(pub Base2DMesh);

impl Modifiable2DMesh {
    /// Creates a mesh with only boundaries defined.
    /// Expects the edges to be sorted (the next edge is starting with the same vertex as the last from the previous edge)
    /// Parents designates the parent at the right side from an edge
    ///
    /// Failing to comply with those invariants might result in unexpected behaviours.
    /// This function is done in a simple version for testing purpose, but will be changed in the future.
    /// If you have any suggestion, do not hesitate to reach out, send an issue or provide a fix.
    ///
    /// # Safety
    ///
    /// The function is marked as unsafe to warn about the very unstable API and the very specific input needed.
    pub unsafe fn new_from_boundary(
        vertices: Vec<Point2<f64>>,
        edge_to_vertices_and_parent: Vec<(
            VertexIndex,
            VertexIndex,
            (ParentIndex, Option<BoundaryPatchIndex>),
        )>,
        parents: Vec<Parent>,
        boundaries: Vec<BoundaryPatch>,
    ) -> Self {
        let mut parents = parents;

        let mut he_to_vertex =
            Vec::<VertexIndex>::with_capacity(edge_to_vertices_and_parent.len() * 2);
        let mut he_to_twin =
            Vec::<HalfEdgeIndex>::with_capacity(edge_to_vertices_and_parent.len() * 2);
        let mut he_to_next_he =
            Vec::<HalfEdgeIndex>::with_capacity(edge_to_vertices_and_parent.len() * 2);
        let mut he_to_prev_he =
            Vec::<HalfEdgeIndex>::with_capacity(edge_to_vertices_and_parent.len() * 2);
        let mut he_to_parent = Vec::<(ParentIndex, Option<BoundaryPatchIndex>)>::with_capacity(
            edge_to_vertices_and_parent.len() * 2,
        );

        let mut parent_to_first_he = Vec::<HalfEdgeIndex>::with_capacity(parents.len() + 1);

        // All arrays are needed to be built correctly

        let cell = ParentIndex(parents.len());
        parents.push(Parent::Cell);

        let mut prev_he_1 = HalfEdgeIndex(0);
        let mut prev_he_2 = HalfEdgeIndex(0);

        for (i, edge) in edge_to_vertices_and_parent.iter().enumerate() {
            let new_he_1 = HalfEdgeIndex(he_to_vertex.len());
            let new_he_2 = HalfEdgeIndex(he_to_vertex.len() + 1);
            he_to_vertex.push(edge.0);
            he_to_vertex.push(edge.1);
            he_to_twin.push(new_he_2);
            he_to_twin.push(new_he_1);

            he_to_parent.push((cell, None));
            he_to_parent.push(edge.2);

            he_to_next_he.push(HalfEdgeIndex(0));
            he_to_next_he.push(HalfEdgeIndex(0));
            he_to_prev_he.push(HalfEdgeIndex(0));
            he_to_prev_he.push(HalfEdgeIndex(0));

            if i == 0 {
                prev_he_1 = new_he_1;
                prev_he_2 = new_he_2;
                continue;
            }

            he_to_next_he[new_he_2] = prev_he_2;
            he_to_prev_he[new_he_1] = prev_he_1;

            he_to_next_he[prev_he_1] = new_he_1;
            he_to_prev_he[prev_he_2] = new_he_2;

            prev_he_1 = new_he_1;
            prev_he_2 = new_he_2;
        }

        he_to_next_he[HalfEdgeIndex(1)] = prev_he_2;
        he_to_prev_he[HalfEdgeIndex(0)] = prev_he_1;

        he_to_next_he[prev_he_1] = HalfEdgeIndex(0);
        he_to_prev_he[prev_he_2] = HalfEdgeIndex(1);

        for i in 0..parents.len() {
            let current_id = ParentIndex(i);
            for (j, he_parent) in he_to_parent.iter().enumerate() {
                if current_id == he_parent.0 {
                    parent_to_first_he.push(HalfEdgeIndex(j));
                    break;
                }
            }
        }

        Modifiable2DMesh(Base2DMesh {
            vertices,

            he_to_vertex,
            he_to_twin,
            he_to_next_he,
            he_to_prev_he,
            he_to_parent,

            parents,
            parent_to_first_he,
            boundaries,
        })
    }

    pub fn validate_topology(self) -> Result<Safe2DMesh, MeshError> {
        self.0.check_mesh()?;
        Ok(Safe2DMesh(self.0))
    }

    /// Gets a mutable reference to a vertex from its index.
    pub fn vertex_mut_from_index(&mut self, vertex_id: VertexIndex) -> &mut Point2<f64> {
        &mut self.0.vertices[vertex_id]
    }

    /// Gets a mutable reference to the parent properties from its index.
    pub fn parent_mut_from_index(&mut self, parent_id: ParentIndex) -> &mut Parent {
        &mut self.0.parents[parent_id]
    }

    /// Creates a new vertex on an half edge at a distance of ```distance_ratio``` (between 0. and 1.) the HalfEdge length
    pub fn split_edge(
        &mut self,
        he_id: HalfEdgeIndex,
        distance_ratio: f64,
    ) -> Result<(), MeshError> {
        if he_id >= HalfEdgeIndex(self.0.he_len()) {
            return Err(MeshError::HalfEdgeIndexOutOfBound {
                got: he_id,
                len: self.0.he_len(),
            });
        }

        if (distance_ratio >= 1.0) | (distance_ratio <= 0.0) {
            return Err(MeshError::WrongFloatValue {
                got: distance_ratio,
                expected: (0.0, 1.0),
            });
        }

        let new_vertex_id = VertexIndex(self.0.vertices.len());
        let new_vertex_pos: Point2<f64> = {
            let edge_vertices = self.0.vertices_from_he(he_id);
            let edge_vertices = (
                self.0.vertices[edge_vertices[0]],
                self.0.vertices[edge_vertices[1]],
            );
            edge_vertices.0.lerp(&edge_vertices.1, distance_ratio)
        };

        let he_ids = (he_id, self.0.he_to_twin[he_id]);

        self.0.vertices.push(new_vertex_pos);

        let new_he_ids = (
            HalfEdgeIndex(self.0.he_to_twin.len()),
            HalfEdgeIndex(self.0.he_to_twin.len() + 1),
        );

        self.0.he_to_vertex.push(new_vertex_id);
        self.0.he_to_vertex.push(new_vertex_id);

        self.0.he_to_twin.push(he_ids.1);
        self.0.he_to_twin.push(he_ids.0);
        self.0.he_to_twin[he_ids.0] = new_he_ids.1;
        self.0.he_to_twin[he_ids.1] = new_he_ids.0;

        let next = self.0.he_to_next_he[he_ids.0];
        self.0.he_to_next_he[he_ids.0] = new_he_ids.0;
        self.0.he_to_next_he.push(next);
        self.0.he_to_prev_he[next] = new_he_ids.0;
        let next = self.0.he_to_next_he[he_ids.1];
        self.0.he_to_next_he[he_ids.1] = new_he_ids.1;
        self.0.he_to_next_he.push(next);
        self.0.he_to_prev_he[next] = new_he_ids.1;

        self.0.he_to_prev_he.push(he_ids.0);
        self.0.he_to_prev_he.push(he_ids.1);

        self.0.he_to_parent.push(self.0.he_to_parent[he_ids.0]);
        self.0.he_to_parent.push(self.0.he_to_parent[he_ids.1]);

        Ok(())
    }

    /// Adds an edge between two vertices.
    /// The vertices must share a common parent.
    /// Returns the new parent.
    ///
    /// The parent must be a cell (not a boundary)
    ///
    /// # Safety
    ///
    /// This function is marked as unsafe to warn about the risk of creating edges crossing each other, leading to wrong parent links.
    /// The edge you are creating must not intersect with a cell boundary.
    ///
    /// This issue will be fixed later by introducing a line crossing check algorithm.
    pub unsafe fn trimming(
        &mut self,
        vertices: (VertexIndex, VertexIndex),
        parent: ParentIndex,
        boundary: Option<BoundaryPatchIndex>,
    ) -> Result<ParentIndex, MeshError> {
        if vertices.0 >= VertexIndex(self.0.vertices_len()) {
            return Err(MeshError::VertexIndexOutOfBound {
                got: vertices.0,
                len: self.0.vertices_len(),
            });
        }
        if vertices.1 >= VertexIndex(self.0.vertices_len()) {
            return Err(MeshError::VertexIndexOutOfBound {
                got: vertices.1,
                len: self.0.vertices_len(),
            });
        }

        for (i, twin) in self.0.he_to_twin.iter().enumerate() {
            let he = HalfEdgeIndex(i);
            if (self.0.he_to_vertex[he] == vertices.0) & (self.0.he_to_vertex[*twin] == vertices.1)
            {
                return Err(MeshError::AlreadyExists);
            }
        }

        let hes_to_vertex = self.0.he_from_vertex(vertices.0);

        let mut he_from_vertex_with_parent = None;
        for he in &hes_to_vertex {
            if self.0.he_to_parent[*he].0 == parent {
                he_from_vertex_with_parent = Some(*he);
                break;
            }
        }

        let he_from_vertex_with_parent = match he_from_vertex_with_parent {
            None => {
                return Err(MeshError::ParentDoesNotContainVertex {
                    vertex: vertices.0,
                    parent,
                })
            }
            Some(value) => value,
        };

        // let hes_to_vertex_2 = self.0.he_from_vertex(vertices.1);

        // let mut he_from_vertex_with_parent_2 = None;
        // for he in &hes_to_vertex_2 {
        //     if self.0.he_to_parent[*he].0 == parent {
        //         he_from_vertex_with_parent_2 = Some(*he);
        //         break;
        //     }
        // }

        // let he_from_vertex_with_parent_2 = match he_from_vertex_with_parent_2 {
        //     None => {
        //         return Err(MeshError::ParentDoesNotContainVertex {
        //             vertex: vertices.0,
        //             parent,
        //         })
        //     }
        //     Some(value) => value,
        // };

        // let mut boundary_index = None;

        // if let Parent::Boundary = self.0.parents[parent] {
        //     if let Some(b_id) = self.0.he_to_parent()[he_from_vertex_with_parent].1 {
        //         if let Some(b_id2) = self.0.he_to_parent()[he_from_vertex_with_parent_2].1 {
        //             if b_id == b_id2 {
        //                 boundary_index = Some(b_id)
        //             }
        //         }
        //         if let Some(b_id2) =
        //             self.0.he_to_parent()[self.0.he_to_next_he()[he_from_vertex_with_parent_2]].1
        //         {
        //             if b_id == b_id2 {
        //                 boundary_index = Some(b_id)
        //             }
        //         }
        //     }
        //     if let Some(b_id) =
        //         self.0.he_to_parent()[self.0.he_to_next_he()[he_from_vertex_with_parent]].1
        //     {
        //         if let Some(b_id2) = self.0.he_to_parent()[he_from_vertex_with_parent_2].1 {
        //             if b_id == b_id2 {
        //                 boundary_index = Some(b_id)
        //             }
        //         }
        //         if let Some(b_id2) =
        //             self.0.he_to_parent()[self.0.he_to_next_he()[he_from_vertex_with_parent_2]].1
        //         {
        //             if b_id == b_id2 {
        //                 boundary_index = Some(b_id)
        //             }
        //         }
        //     }
        //     if let None = boundary_index {
        //         boundary_index = self.0.he_to_parent()[he_from_vertex_with_parent].1;
        //     }
        // }

        let new_he = self.0.he_len();
        self.0.he_to_vertex.push(vertices.1);
        self.0.he_to_vertex.push(vertices.0);
        self.0.he_to_twin.push(HalfEdgeIndex(new_he + 1));
        self.0.he_to_twin.push(HalfEdgeIndex(new_he));

        let new_cell = self.0.parents_len();
        self.0.parents.push(Parent::Cell);
        self.0.he_to_parent.push((parent, boundary));
        self.0.he_to_parent.push((ParentIndex(new_cell), None));

        self.0.he_to_next_he.push(he_from_vertex_with_parent);
        self.0.he_to_next_he.push(HalfEdgeIndex(usize::MAX)); // Placeholder to be quite sure to have an error when checking the mesh if the value is not set correctly later in the function
        self.0.he_to_prev_he.push(HalfEdgeIndex(usize::MAX));
        self.0
            .he_to_prev_he
            .push(self.0.he_to_prev_he[he_from_vertex_with_parent]);
        self.0.he_to_next_he[self.0.he_to_prev_he[he_from_vertex_with_parent]] =
            HalfEdgeIndex(new_he + 1);
        self.0.he_to_prev_he[he_from_vertex_with_parent] = HalfEdgeIndex(new_he);

        let mut current_he = he_from_vertex_with_parent;

        let mut i = 0;
        while self.0.he_to_vertex[current_he] != vertices.1 {
            self.0.he_to_parent[current_he].0 = parent;
            current_he = self.0.he_to_next_he[current_he];
            i += 1;
            if i >= self.0.he_len() {
                return Err(MeshError::ParentDoesNotContainVertex {
                    vertex: vertices.1,
                    parent,
                });
            }
        }

        self.0.he_to_prev_he[HalfEdgeIndex(new_he)] = self.0.he_to_prev_he[current_he];
        self.0.he_to_next_he[HalfEdgeIndex(new_he + 1)] = current_he;
        self.0.he_to_next_he[self.0.he_to_prev_he[current_he]] = HalfEdgeIndex(new_he);
        self.0.he_to_prev_he[current_he] = HalfEdgeIndex(new_he + 1);

        let mut i = 0;
        while self.0.he_to_vertex[current_he] != vertices.0 {
            self.0.he_to_parent[current_he] = (ParentIndex(new_cell), None);
            current_he = self.0.he_to_next_he[current_he];
            i += 1;
            if i >= self.0.he_len() {
                // Should not be possible
                return Err(MeshError::ParentDoesNotContainVertex {
                    vertex: vertices.0,
                    parent,
                });
            }
        }

        self.0.parent_to_first_he[parent] = HalfEdgeIndex(new_he);

        self.0.parent_to_first_he.push(HalfEdgeIndex(new_he + 1));

        Ok(ParentIndex(new_cell))
    }

    /// Creates a triangle
    ///
    /// # Safety
    ///
    /// You must specify the right half_edge corresponding to the direction you want your new triangle.
    pub unsafe fn notching(
        &mut self,
        he: HalfEdgeIndex,
        pos: Point2<f64>,
    ) -> Result<ParentIndex, MeshError> {
        if he >= HalfEdgeIndex(self.0.he_len()) {
            return Err(MeshError::HalfEdgeIndexOutOfBound {
                got: he,
                len: self.0.he_len(),
            });
        }

        let parent = self.0.he_to_parent[self.0.he_to_twin[he]];
        let parent_2 = self.0.he_to_parent[self.0.he_to_twin[he]];

        let mut boundary = None;
        if let Parent::Boundary = self.0.parents[parent_2.0] {
            boundary = parent_2.1;
        }

        let vertices = self.0.vertices_from_he(he);
        let new_vertex = self.0.vertices_len();
        self.split_edge(he, 0.5)?;
        self.0.vertices[new_vertex] = pos;
        let new_parent;
        unsafe {
            new_parent = self.trimming((vertices[0], vertices[1]), parent.0, boundary)?;
        }

        Ok(new_parent)
    }

    /// Meant to implement mesh adjustement using delaunay condition.
    ///
    /// Parents needs to be adjacent triangles and the triangles must not have aligned edges
    ///
    /// Can't swap edges for a boundary (ParentIndex must not point to a boundary)
    pub fn swap_edge(&mut self, parents: (ParentIndex, ParentIndex)) -> Result<(), MeshError> {
        if self.0.vertices_from_parent(parents.0).len() != 3 {
            return Err(MeshError::ParentNotTriangle { parent: parents.0 });
        }
        if self.0.vertices_from_parent(parents.1).len() != 3 {
            return Err(MeshError::ParentNotTriangle { parent: parents.1 });
        }

        let mut common_he = None;
        'outer: for he_0 in self.0.he_from_parent(parents.0) {
            for he_1 in self.0.he_from_parent(parents.1) {
                if he_0 == self.0.he_to_twin()[he_1] {
                    common_he = Some(he_0);
                    break 'outer;
                }
            }
        }

        let common_he = match common_he {
            None => {
                return Err(MeshError::NoCommonEdge {
                    parent_0: parents.0,
                    parent_1: parents.1,
                })
            }
            Some(he) => he,
        };

        // Not implemented using an edge removal followed by a trimming due to all the modifications needed for an edge removal
        let hes_0 = [
            common_he,
            self.0.he_to_next_he()[common_he],
            self.0.he_to_next_he()[self.0.he_to_next_he()[common_he]],
        ];
        let hes_1 = [
            self.0.he_to_twin()[common_he],
            self.0.he_to_next_he()[self.0.he_to_twin()[common_he]],
            self.0.he_to_next_he()[self.0.he_to_next_he()[self.0.he_to_twin()[common_he]]],
        ];

        if let Parent::Boundary = self.0.parents[parents.0] {
            return Err(MeshError::BoundaryEdgeSwap {
                parent_0: parents.0,
                parent_1: parents.1,
            });
        }
        if let Parent::Boundary = self.0.parents[parents.1] {
            return Err(MeshError::BoundaryEdgeSwap {
                parent_0: parents.0,
                parent_1: parents.1,
            });
        }

        // Check that no edges are aligned to avoid triangle with no area creation
        let vec_0 = self.0.he_vector(hes_0[1]);
        let vec_1 = self.0.he_vector(hes_1[2]);
        if vec_0.x / vec_1.x == vec_0.y / vec_1.y {
            return Err(MeshError::AllignedEdges {
                parent_0: parents.0,
                parent_1: parents.1,
            });
        }

        // Check that no edges are aligned to avoid triangle with no area creation
        let vec_0 = self.0.he_vector(hes_0[2]);
        let vec_1 = self.0.he_vector(hes_1[1]);
        if vec_0.x / vec_1.x == vec_0.y / vec_1.y {
            return Err(MeshError::AllignedEdges {
                parent_0: parents.0,
                parent_1: parents.1,
            });
        }

        self.0.he_to_vertex[hes_0[0]] = self.0.he_to_vertex[hes_0[2]];
        self.0.he_to_vertex[hes_1[0]] = self.0.he_to_vertex[hes_1[2]];

        self.0.he_to_next_he[hes_0[0]] = hes_1[2];
        self.0.he_to_prev_he[hes_1[2]] = hes_0[0];

        self.0.he_to_next_he[hes_0[1]] = hes_0[0];
        self.0.he_to_prev_he[hes_0[0]] = hes_0[1];

        self.0.he_to_next_he[hes_0[2]] = hes_1[1];
        self.0.he_to_prev_he[hes_1[1]] = hes_0[2];

        self.0.he_to_next_he[hes_1[0]] = hes_0[2];
        self.0.he_to_prev_he[hes_0[2]] = hes_1[0];

        self.0.he_to_next_he[hes_1[1]] = hes_1[0];
        self.0.he_to_prev_he[hes_1[0]] = hes_1[1];

        self.0.he_to_next_he[hes_1[2]] = hes_0[1];
        self.0.he_to_prev_he[hes_0[1]] = hes_1[2];

        self.0.he_to_parent[hes_0[2]] = (parents.1, None);
        self.0.he_to_parent[hes_1[2]] = (parents.0, None);

        self.0.parent_to_first_he[parents.0] = hes_0[0];
        self.0.parent_to_first_he[parents.1] = hes_1[0];

        Ok(())
    }
}
