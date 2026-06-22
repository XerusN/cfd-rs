use cfd_rs_utils::mesh::assembled_mesh::MeshCore;

use crate::finite_volume::solvers::SolverCore;

pub mod grid;

pub trait OutputConfig: Sized {
    fn try_output<M: MeshCore>(&mut self, solver: &SolverCore<M>);
}