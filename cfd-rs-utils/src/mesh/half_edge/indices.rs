use std::fmt;
use std::ops::{Index, IndexMut};

use super::Parent;
use crate::mesh::computational_mesh::{BoundaryPatch, Cell};
use nalgebra::Point2;
use serde::{Deserialize, Serialize};

#[derive(Copy, Clone, Debug, Default, PartialEq, PartialOrd, Serialize, Deserialize, Ord, Eq)]
pub struct ParentIndex(pub usize);

impl Index<ParentIndex> for Vec<HalfEdgeIndex> {
    type Output = HalfEdgeIndex;

    fn index(&self, index: ParentIndex) -> &Self::Output {
        &self[index.0]
    }
}

impl IndexMut<ParentIndex> for Vec<HalfEdgeIndex> {
    fn index_mut(&mut self, index: ParentIndex) -> &mut Self::Output {
        &mut self[index.0]
    }
}

impl Index<ParentIndex> for Vec<Parent> {
    type Output = Parent;

    fn index(&self, index: ParentIndex) -> &Self::Output {
        &self[index.0]
    }
}

impl IndexMut<ParentIndex> for Vec<Parent> {
    fn index_mut(&mut self, index: ParentIndex) -> &mut Self::Output {
        &mut self[index.0]
    }
}

impl fmt::Display for ParentIndex {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}

#[derive(Copy, Clone, Debug, Default, PartialEq, PartialOrd, Serialize, Deserialize, Ord, Eq)]
pub struct HalfEdgeIndex(pub usize);

impl Index<HalfEdgeIndex> for Vec<HalfEdgeIndex> {
    type Output = HalfEdgeIndex;

    fn index(&self, index: HalfEdgeIndex) -> &Self::Output {
        &self[index.0]
    }
}

impl IndexMut<HalfEdgeIndex> for Vec<HalfEdgeIndex> {
    fn index_mut(&mut self, index: HalfEdgeIndex) -> &mut Self::Output {
        &mut self[index.0]
    }
}

impl Index<HalfEdgeIndex> for Vec<(ParentIndex, Option<BoundaryPatchIndex>)> {
    type Output = (ParentIndex, Option<BoundaryPatchIndex>);

    fn index(&self, index: HalfEdgeIndex) -> &Self::Output {
        &self[index.0]
    }
}

impl IndexMut<HalfEdgeIndex> for Vec<(ParentIndex, Option<BoundaryPatchIndex>)> {
    fn index_mut(&mut self, index: HalfEdgeIndex) -> &mut Self::Output {
        &mut self[index.0]
    }
}

impl Index<HalfEdgeIndex> for Vec<VertexIndex> {
    type Output = VertexIndex;

    fn index(&self, index: HalfEdgeIndex) -> &Self::Output {
        &self[index.0]
    }
}

impl IndexMut<HalfEdgeIndex> for Vec<VertexIndex> {
    fn index_mut(&mut self, index: HalfEdgeIndex) -> &mut Self::Output {
        &mut self[index.0]
    }
}

impl fmt::Display for HalfEdgeIndex {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}

#[derive(Copy, Clone, Debug, Default, PartialEq, PartialOrd, Serialize, Deserialize, Ord, Eq)]
pub struct VertexIndex(pub usize);

impl Index<VertexIndex> for Vec<Point2<f64>> {
    type Output = Point2<f64>;

    fn index(&self, index: VertexIndex) -> &Self::Output {
        &self[index.0]
    }
}

impl IndexMut<VertexIndex> for Vec<Point2<f64>> {
    fn index_mut(&mut self, index: VertexIndex) -> &mut Self::Output {
        &mut self[index.0]
    }
}

impl fmt::Display for VertexIndex {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}

#[derive(Copy, Clone, Debug, Default, PartialEq, PartialOrd, Serialize, Deserialize, Ord, Eq)]
pub struct FaceIndex(pub usize);

impl Index<FaceIndex> for Vec<Point2<f64>> {
    type Output = Point2<f64>;

    fn index(&self, index: FaceIndex) -> &Self::Output {
        &self[index.0]
    }
}

impl IndexMut<FaceIndex> for Vec<Point2<f64>> {
    fn index_mut(&mut self, index: FaceIndex) -> &mut Self::Output {
        &mut self[index.0]
    }
}

impl fmt::Display for FaceIndex {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}

#[derive(Copy, Clone, Debug, Default, PartialEq, PartialOrd, Serialize, Deserialize, Ord, Eq)]
pub struct CellIndex(pub usize);

impl Index<CellIndex> for Vec<Cell> {
    type Output = Cell;

    fn index(&self, index: CellIndex) -> &Self::Output {
        &self[index.0]
    }
}

impl IndexMut<CellIndex> for Vec<Cell> {
    fn index_mut(&mut self, index: CellIndex) -> &mut Self::Output {
        &mut self[index.0]
    }
}

impl fmt::Display for CellIndex {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}

#[derive(Copy, Clone, Debug, Default, PartialEq, PartialOrd, Serialize, Deserialize, Ord, Eq)]
pub struct BoundaryPatchIndex(pub usize);

impl Index<BoundaryPatchIndex> for Vec<BoundaryPatch> {
    type Output = BoundaryPatch;

    fn index(&self, index: BoundaryPatchIndex) -> &Self::Output {
        &self[index.0]
    }
}

impl IndexMut<BoundaryPatchIndex> for Vec<BoundaryPatch> {
    fn index_mut(&mut self, index: BoundaryPatchIndex) -> &mut Self::Output {
        &mut self[index.0]
    }
}

impl fmt::Display for BoundaryPatchIndex {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}
