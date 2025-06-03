use cfd_rs::finite_volume::boundary::BoundaryValue;
use cfd_rs::finite_volume::case::simple::SimpleCase;
use cfd_rs::finite_volume::config::{InitFunc, OutputConfig};
use cfd_rs_utils::control::OutputControl;
use hashbrown::HashMap;

use cfd_rs::finite_volume::case::{self, poisson};
use cfd_rs::finite_volume::{
    boundary::{BoundaryCondition, FieldsBoundaryConditions},
    case::poisson::PoissonCase,
    case::Case,
    config::{CaseConfig, GeometryConfig, Schemes},
    discretizations::{
        convection::ConvectionScheme, divergence::DivergenceScheme, laplacian::LaplacianScheme,
        time_schemes::TimeIntegration,
    },
    equation::{Dimension, Variable},
    gradients::{GradientConfig, GradientScheme},
    interpolations::GradientInterpConfig,
};
use nalgebra::{Point2, Vector2};

fn poisson() -> CaseConfig {
    let schemes = Schemes {
        transient: TimeIntegration::ForwardEuler,
        convection: ConvectionScheme::UpwindSecondOrder,
        laplacian: LaplacianScheme::OrthogonalCorrection,
        divergence: DivergenceScheme::Basic,
        gradients: GradientConfig {
            scheme: GradientScheme::GreenGaussCompact,
            interp: GradientInterpConfig::AveragedCorrected,
        },
    };

    let geometry = GeometryConfig {
        import_path: Some("./target/exports/mesh.cfd".to_string()),
        element_size: 0.01,
    };

    let output = OutputConfig {
        control: OutputControl::Iteration(1),
        directory: "./target/exports".to_string(),
    };

    let mut bc_fields = HashMap::new();
    let bc = vec![
        BoundaryCondition::Dirichlet(BoundaryValue::Scalar(0.)),
        BoundaryCondition::Dirichlet(BoundaryValue::Scalar(0.)),
        BoundaryCondition::Dirichlet(BoundaryValue::Scalar(200.)),
        BoundaryCondition::Dirichlet(BoundaryValue::Scalar(300.)),
    ];
    bc_fields.insert(Variable::new("T".to_string(), Dimension::Scalar), bc);
    let bc_fields = FieldsBoundaryConditions::new(bc_fields);

    let mut initial_fields = HashMap::new();
    initial_fields.insert(
        Variable::new("T".to_string(), Dimension::Scalar),
        InitFunc::Scalar(constant),
    );

    CaseConfig {
        schemes,
        geometry,
        bc: bc_fields,
        output,
        initial_fields,
    }
}

fn simple() -> CaseConfig {
    let schemes = Schemes {
        transient: TimeIntegration::ForwardEuler,
        convection: ConvectionScheme::UpwindSecondOrder,
        laplacian: LaplacianScheme::OrthogonalCorrection,
        divergence: DivergenceScheme::Basic,
        gradients: GradientConfig {
            scheme: GradientScheme::GreenGaussCompact,
            interp: GradientInterpConfig::AveragedCorrected,
        },
    };

    let geometry = GeometryConfig {
        import_path: Some("./target/exports/mesh.cfd".to_string()),
        element_size: 0.01,
    };

    let output = OutputConfig {
        control: OutputControl::Iteration(1),
        directory: "./target/exports".to_string(),
    };

    let mut bc_fields = HashMap::new();
    let bc = vec![
        BoundaryCondition::Dirichlet(BoundaryValue::Vector2(Vector2::new(0., 0.))),
        BoundaryCondition::Dirichlet(BoundaryValue::Vector2(Vector2::new(0., 0.))),
        BoundaryCondition::Dirichlet(BoundaryValue::Vector2(Vector2::new(1., 0.))),
        BoundaryCondition::Dirichlet(BoundaryValue::Vector2(Vector2::new(0., 0.))),
    ];
    bc_fields.insert(Variable::new("U".to_string(), Dimension::Vector2), bc);
    let bc = vec![
        BoundaryCondition::Neumann(BoundaryValue::Scalar(0.)),
        BoundaryCondition::Neumann(BoundaryValue::Scalar(0.)),
        BoundaryCondition::Neumann(BoundaryValue::Scalar(0.)),
        BoundaryCondition::Neumann(BoundaryValue::Scalar(0.)),
    ];
    bc_fields.insert(Variable::new("P".to_string(), Dimension::Scalar), bc);
    let bc_fields = FieldsBoundaryConditions::new(bc_fields);

    let mut initial_fields = HashMap::new();
    initial_fields.insert(
        Variable::new("Phi".to_string(), Dimension::Scalar),
        InitFunc::Scalar(constant),
    );

    CaseConfig {
        schemes,
        geometry,
        bc: bc_fields,
        output,
        initial_fields,
    }
}

fn convection_setup() -> CaseConfig {
    let schemes = Schemes {
        transient: TimeIntegration::ForwardEuler,
        convection: ConvectionScheme::UpwindSecondOrder,
        laplacian: LaplacianScheme::OrthogonalCorrection,
        divergence: DivergenceScheme::Basic,
        gradients: GradientConfig {
            scheme: GradientScheme::GreenGaussCompact,
            interp: GradientInterpConfig::AveragedCorrected,
        },
    };

    let geometry = GeometryConfig {
        import_path: None,
        element_size: 0.01,
    };

    let output = OutputConfig {
        control: OutputControl::Iteration(1),
        directory: "./target/exports".to_string(),
    };

    let mut bc_fields = HashMap::new();
    let bc = vec![
        BoundaryCondition::Dirichlet(BoundaryValue::Scalar(1.)),
        BoundaryCondition::Dirichlet(BoundaryValue::Scalar(0.)),
        BoundaryCondition::Neumann(BoundaryValue::Scalar(0.)),
    ];
    bc_fields.insert(Variable::new("Phi".to_string(), Dimension::Scalar), bc);
    let bc_fields = FieldsBoundaryConditions::new(bc_fields);

    let mut initial_fields = HashMap::new();
    initial_fields.insert(
        Variable::new("Phi".to_string(), Dimension::Scalar),
        InitFunc::Scalar(constant),
    );

    CaseConfig {
        schemes,
        geometry,
        bc: bc_fields,
        output,
        initial_fields,
    }
}

pub fn constant(_point: &Point2<f64>) -> f64 {
    0.
}

fn main() {
    let mut case = PoissonCase::new(poisson());
    //let mut case = simple();

    for _ in 0..2 {
        case.next_step();
        case.export().unwrap();
    }
}
