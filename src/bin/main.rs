use hashbrown::HashMap;

use cfd_rs::finite_volume::solvers::poisson;
use cfd_rs::finite_volume::{
    boundary::{BoundaryCondition, FieldsBoundaryConditions},
    case::Case,
    config::{CaseConfig, GeometryConfig, Schemes},
    discretizations::{
        convection::ConvectionScheme, divergence::DivergenceScheme, laplacian::LaplacianScheme,
        time_schemes::TimeIntegration,
    },
    equation::{Dimension, Variable},
    gradients::{GradientConfig, GradientScheme},
    interpolations::GradientInterpConfig,
    solvers::poisson::PoissonCase,
};
use cfd_rs_utils::geometry;

fn main() {
    let schemes = Schemes {
        transient: TimeIntegration::ForwardEuler,
        convection: ConvectionScheme::UpwindSecondOrder,
        laplacian: LaplacianScheme::Centered,
        divergence: DivergenceScheme::Centered,
        gradients: GradientConfig {
            scheme: GradientScheme::GreenGaussCompact,
            interp: GradientInterpConfig::AveragedCorrected,
        },
    };

    let geometry = GeometryConfig {};
    let bc = vec![
        BoundaryCondition::Neumann(300.),
        BoundaryCondition::Neumann(200.),
        BoundaryCondition::Neumann(0.),
        BoundaryCondition::Neumann(0.),
    ];
    let mut bc_fields = HashMap::new();
    bc_fields.insert(Variable::new("P".to_string(), Dimension::Scalar), bc);
    let bc_fields = FieldsBoundaryConditions::new(bc_fields);
    let config = CaseConfig {
        schemes,
        geometry,
        bc: bc_fields,
    };

    let mut case = PoissonCase::new(config);

    case.next_step();

    #[cfg(debug_assertions)]
    case.export("./target/exports");
}
