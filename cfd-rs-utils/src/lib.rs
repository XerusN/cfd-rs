//! Lib used as support for cfd-rs and cfd-rs-mesh.
//! Mostly for type definition and common functions.
//!
//! For now only functionnalities to construct an advancing front algorithm are being created.
//! Mesh is represented using array based half-edge data structure

pub mod control;
pub mod errors;
pub mod geometry;
pub mod mesh;
