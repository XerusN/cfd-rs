use cfd_rs_utils::mesh::assembled_mesh::MeshCore;

use crate::{
    finite_volume::solvers::SolverCore,
    post_processing::{grid::GridPostprocConfig, plot::PlotPostprocConfig},
};

pub mod grid;
pub mod plot;

#[derive(Debug, Clone, PartialEq)]
pub enum OutputConfig {
    Grid(GridPostprocConfig),
    Plot(PlotPostprocConfig),
}
