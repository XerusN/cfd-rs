use thiserror::Error;

use crate::mesh::indices::*;

#[derive(Clone, Debug, Default, Error, PartialEq)]
pub enum MeshError {
    #[default]
    #[error(
        "An Unspecified error happened, you can blame the crate developer for the lack of details"
    )]
    Unspecified,
    #[error("Value not in the right range (expected in {expected:?}, got {got:?}")]
    WrongFloatValue { got: f64, expected: (f64, f64) },
    #[error("An HalfEdgeIndex is out of bound in array (got {got:?}, there are only {len:?} half-edges)")]
    HalfEdgeIndexOutOfBound { got: HalfEdgeIndex, len: usize },
    #[error(
        "An VertexIndex is out of bound in array (got {got:?}, there are only {len:?} vertices)"
    )]
    VertexIndexOutOfBound { got: VertexIndex, len: usize },
    #[error(
        "An ParentIndex is out of bound in array (got {got:?}, there are only {len:?} parents)"
    )]
    ParentIndexOutOfBound { got: ParentIndex, len: usize },
    #[error("Twins are badly set (origin : {he:?}, twin : {he_twin:?} and twin from twin : {he_twin_twin:?})")]
    TwinNotCorrect {
        he: HalfEdgeIndex,
        he_twin: HalfEdgeIndex,
        he_twin_twin: HalfEdgeIndex,
    },
    #[error("Links between parents and half-edges are badly set (origin : {parent:?}, half-edge : {he:?} and parent from half_edge : {he_parent:?})")]
    ParentNotCorrect {
        parent: ParentIndex,
        he: HalfEdgeIndex,
        he_parent: ParentIndex,
    },
    #[error("No cycle when starting from half-edge (origin : {he:?})")]
    WrongHalfEdgeLoop { he: HalfEdgeIndex },
    #[error("Wrong definition of next or previous half-edge (origin : {he:?}, next : {he_next:?} and prev from next : {he_next_prev:?})")]
    NextPrevNotCorrect {
        he: HalfEdgeIndex,
        he_next: HalfEdgeIndex,
        he_next_prev: HalfEdgeIndex,
    },
    #[error("The object already exists")]
    AlreadyExists,
    #[error("The vertices is not contained in parent (vertex : {vertex:?}, parent : {parent:?})")]
    ParentDoesNotContainVertex {
        vertex: VertexIndex,
        parent: ParentIndex,
    },
    #[error("The mesh passed as input is badly constructed")]
    WrongMeshInitialisation,
    #[error("No Element is creatable on this edge")]
    NoElementCreatable(HalfEdgeIndex),
    #[error("Maximum iteration reached ({max_it:?})")]
    MaxIterationReached { max_it: usize },
    #[error("Parent needs to be a triangle ({parent:?})")]
    ParentNotTriangle { parent: ParentIndex },
    #[error("No common edges ({parent_0:?} and {parent_1:?})")]
    NoCommonEdge {
        parent_0: ParentIndex,
        parent_1: ParentIndex,
    },
    #[error("Parents have aligned edges ({parent_0:?} and {parent_1:?})")]
    AllignedEdges {
        parent_0: ParentIndex,
        parent_1: ParentIndex,
    },
    #[error("Can't swap edge with a boundary ({parent_0:?} and {parent_1:?})")]
    BoundaryEdgeSwap {
        parent_0: ParentIndex,
        parent_1: ParentIndex,
    },
}
