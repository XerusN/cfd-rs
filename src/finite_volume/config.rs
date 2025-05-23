use cfd_rs_utils::control::OutputControl;

use super::{
    boundary::FieldsBoundaryConditions,
    discretizations::{
        convection::ConvectionScheme, divergence::DivergenceScheme, laplacian::LaplacianScheme,
        time_schemes::TimeIntegration,
    },
    gradients::GradientConfig,
};

/// Switch to private fields
#[derive(Clone, PartialEq, Debug)]
pub struct CaseConfig {
    pub schemes: Schemes,
    pub geometry: GeometryConfig,
    pub bc: FieldsBoundaryConditions,
    pub output: OutputConfig,
}

#[derive(Clone, PartialEq, Debug)]
pub struct GeometryConfig {}

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
