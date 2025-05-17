use super::{
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
