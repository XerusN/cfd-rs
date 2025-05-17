use hashbrown::HashMap;

use cfd_rs::finite_volume::{boundary::{BoundaryCondition, FieldsBoundaryConditions}, config::{GeometryConfig, Schemes}, discretizations::{convection::ConvectionScheme, divergence::DivergenceScheme, laplacian::LaplacianScheme, time_schemes::TimeIntegration}, equation::{Dimension, Variable}, gradients::{GradientConfig, GradientScheme}, interpolations::GradientInterpConfig};
use cfd_rs_utils::geometry;
use cfd_rs::finite_volume::solvers::poisson;

fn main() {
    
    let schemes = Schemes{
        transient: TimeIntegration::ForwardEuler,
        convection: ConvectionScheme::UpwindSecondOrder,
        laplacian: LaplacianScheme::Centered,
        divergence: DivergenceScheme::Centered,
        gradients: GradientConfig{
            scheme: GradientScheme::GreenGaussCompact,
            interp: GradientInterpConfig::AveragedCorrected,
        },
    };
    
    let geometry = GeometryConfig{};
    let bc = vec![BoundaryCondition::Neumann(300.), BoundaryCondition::Neumann(200.), BoundaryCondition::Neumann(0.), BoundaryCondition::Neumann(0.)];
    let mut bc_fields = HashMap::new();
    bc_fields.insert(Variable::new("P".to_string(), Dimension::Scalar), bc);
    let bc_fields = FieldsBoundaryConditions::new(bc_fields);
    
}