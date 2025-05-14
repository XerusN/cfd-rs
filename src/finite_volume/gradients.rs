use nalgebra::Vector2;
use super::{case::Case, equation::*, interpolations::GradientInterpConfig};

#[derive(Clone, Debug, PartialEq)]
pub enum GradientScheme {
    GreenGaussCompact,
    GreenGaussExtended,
    LeastSquare,
}

#[derive(Clone, Debug, PartialEq)]
pub struct GradientConfig {
    pub schemes: GradientScheme,
    pub interp: GradientInterpConfig,
}