use cfd_rs::finite_volume::config::OutputConfig;
use cfd_rs::finite_volume::solvers::simple::SimpleCase;
use cfd_rs_utils::control::OutputControl;
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

fn main() {
    let schemes = Schemes {
        transient: TimeIntegration::ForwardEuler,
        convection: ConvectionScheme::UpwindSecondOrder,
        laplacian: LaplacianScheme::OrthogonalCorrection,
        divergence: DivergenceScheme::Centered,
        gradients: GradientConfig {
            scheme: GradientScheme::GreenGaussCompact,
            interp: GradientInterpConfig::AveragedCorrected,
        },
    };

    let geometry = GeometryConfig {
        import_path: Some("./target/exports/mesh.cfd".to_string()),
        element_size: 0.01,
    };
    let bc = vec![
        BoundaryCondition::Dirichlet(300.),
        BoundaryCondition::Dirichlet(200.),
        BoundaryCondition::Dirichlet(0.),
        BoundaryCondition::Dirichlet(0.),
    ];
    let output = OutputConfig {
        control: OutputControl::Iteration(1),
        directory: "./target/exports".to_string(),
    };
    let mut bc_fields = HashMap::new();
    bc_fields.insert(Variable::new("P".to_string(), Dimension::Scalar), bc);
    let bc_fields = FieldsBoundaryConditions::new(bc_fields);
    let config = CaseConfig {
        schemes,
        geometry,
        bc: bc_fields,
        output,
    };

    let mut case = SimpleCase::new(config);

    //case.mesh().serialize_file(&"./target/exports/mesh.cfd").unwrap();

    println!("{}", case.name());

    // #[cfg(debug_assertions)]
    // case.export().unwrap();

    case.next_step();

    // #[cfg(debug_assertions)]
    case.export().unwrap();
}
