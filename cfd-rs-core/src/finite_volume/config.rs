use cfd_rs_utils::control::OutputControl;
use hashbrown::HashMap;
use nalgebra::{Point2, Vector2};

use crate::finite_volume::gradients::{GradientInterp, GradientScheme};

use super::{
    boundary::FieldsBoundaryConditions,
    equation::{
        discretizations::{
            convection::ConvectionScheme, divergence::DivergenceScheme, laplacian::LaplacianScheme,
            time_schemes::TimeIntegration,
        },
        variables::Variable,
    },
};

// /// Switch to private fields
// #[derive(Clone, Debug, PartialEq)]
// pub struct CaseConfig {
//     pub geometry: GeometryConfig,
//     pub output: OutputConfig,
// }

#[derive(Clone, PartialEq, Debug)]
pub struct GeometryConfig {
    pub import_path: Option<String>,
    pub meshing: MeshingConfig,
}

#[derive(Clone, PartialEq, Debug)]
pub enum MeshingConfig {
    Cartesian {
        length: Vector2<f64>,
        n_elements: Vector2<usize>,
    },
    Line {
        length: f64,
        n_elements: usize,
    },
    AdvancingFront {
        element_size: f64,
    },
}

#[derive(Clone, PartialEq, Debug)]
pub struct Schemes {
    pub transient: TimeIntegration,
    pub laplacian: LaplacianScheme,
    pub convection: ConvectionScheme,
    pub divergence: DivergenceScheme,
}

#[derive(Clone, PartialEq, Debug, Default)]
pub struct SchemesConfig {
    pub transient: Option<TimeIntegration>,
    pub laplacian: Option<LaplacianScheme>,
    pub convection: Option<ConvectionScheme>,
    pub divergence: Option<DivergenceScheme>,
}

#[derive(Clone, PartialEq, Debug, Default)]
pub struct GradientConfig {
    pub scheme: Option<GradientScheme>,
    pub interp: Option<GradientInterp>,
}

#[derive(Clone, PartialEq, Debug)]
pub struct OutputConfig {
    pub control: OutputControl,
    pub directory: String,
}

#[derive(Clone)]
pub enum InitFunc<'a> {
    Scalar(&'a dyn Fn(&Point2<f64>) -> f64),
    Vector2(
        &'a dyn Fn(&Point2<f64>) -> f64,
        &'a dyn Fn(&Point2<f64>) -> f64,
    ),
}
