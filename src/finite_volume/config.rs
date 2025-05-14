use super::{
    discretizations::{
        convection::ConvectionScheme, divergence::DivergenceScheme, laplacian::LaplacianScheme,
        time_schemes::TimeIntegration,
    },
    gradients::GradientScheme,
    interpolations::GradientInterpConfig,
};

#[derive(Clone, PartialEq, Debug)]
pub enum CaseConfig {
    Simple(Schemes, GeometryConfig),
}

#[derive(Clone, PartialEq, Debug)]
pub struct GeometryConfig {}

#[derive(Clone, PartialEq, Debug)]
pub struct Schemes {
    pub transient: TimeIntegration,
    pub gradients: (GradientScheme, GradientInterpConfig),
    pub laplacian: LaplacianScheme,
    pub convection: ConvectionScheme,
    pub divergence: DivergenceScheme,
}
