use super::{discretizations::{time_schemes::TimeScheme, laplacian::LaplacianScheme, convection::ConvectionScheme, divergence::DivergenceScheme}, gradients::GradientScheme, interpolations::{DecompositionConfig, GradientInterpConfig, InterpolationConfig}};

#[derive(Clone, PartialEq, Debug)]
pub enum CaseConfig {
    Simple(Schemes, GeometryConfig),
}

#[derive(Clone, PartialEq, Debug)]
pub struct GeometryConfig {
    
}

#[derive(Clone, PartialEq, Debug)]
pub struct Schemes {
    pub transient: TimeScheme,
    pub gradients: (GradientScheme, GradientInterpConfig),
    pub laplacian: LaplacianScheme,
    pub convection: ConvectionScheme,
    pub divergence: DivergenceScheme,
}