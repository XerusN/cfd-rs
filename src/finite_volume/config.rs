use cfd_rs_utils::control::OutputControl;
use hashbrown::HashMap;
use nalgebra::{Point2, Vector2};

use super::{
    boundary::FieldsBoundaryConditions,
    discretizations::{
        convection::ConvectionScheme, divergence::DivergenceScheme, laplacian::LaplacianScheme,
        time_schemes::TimeIntegration,
    },
    equation::Variable,
    gradients::GradientConfig,
};

/// Switch to private fields
#[derive(Clone, PartialEq, Debug)]
pub struct CaseConfig {
    pub schemes: Schemes,
    pub initial_fields: HashMap<Variable, InitFunc>,
    pub geometry: GeometryConfig,
    pub bc: FieldsBoundaryConditions,
    pub output: OutputConfig,
}

#[derive(Clone, PartialEq, Debug)]
pub struct GeometryConfig {
    pub import_path: Option<String>,
    pub element_size: f64,
}

#[derive(Clone, PartialEq, Debug)]
pub struct Schemes {
    pub transient: TimeIntegration,
    pub gradients: GradientConfig,
    pub laplacian: LaplacianScheme,
    pub convection: ConvectionScheme,
    pub divergence: DivergenceScheme,
}

#[derive(Clone, PartialEq, Debug)]
pub struct OutputConfig {
    pub control: OutputControl,
    pub directory: String,
}

#[derive(Clone, PartialEq, Debug)]
pub enum InitFunc {
    Scalar(fn(&Point2<f64>) -> f64),
    Vector2(fn(&Point2<f64>) -> f64, fn(&Point2<f64>) -> f64),
}
