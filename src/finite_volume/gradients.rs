use crate::finite_volume::interpolations::GradientInterpConfig;

#[derive(Clone, Debug, PartialEq)]
pub enum GradientConfig {
    GreenGaussCompact(GradientInterpConfig),
    GreenGaussExtended(GradientInterpConfig),
    LeastSquare(GradientInterpConfig),
}