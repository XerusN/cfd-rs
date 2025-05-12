use super::{discretisations::{SpaceDiscretizationConfig, TimeDiscretizationConfig}, gradients::GradientConfig, interpolations::{DecompositionConfig, GradientInterpConfig, InterpolationConfig}};

#[derive(Clone, PartialEq, Debug)]
pub enum CaseConfig {
    Simple(SimpleSchemesConfig, GeometryConfig),
}

#[derive(Clone, PartialEq, Debug)]
pub struct GeometryConfig {
    
}

#[derive(Clone, PartialEq, Debug)]
pub struct SimpleSchemesConfig {
    pub transient: TimeDiscretizationConfig,
    
    pub gradients: (GradientConfig, GradientInterpConfig),
    
    pub momentum_eq: MomentumEqConfig,
    
}

#[derive(Clone, PartialEq, Debug)]
pub struct MomentumEqConfig {
    pub convection: SpaceDiscretizationConfig,
    pub diffusion: (SpaceDiscretizationConfig, DecompositionConfig),
}

pub struct PressureCorrEq {
    pub speed: InterpolationConfig,
    pub 
}