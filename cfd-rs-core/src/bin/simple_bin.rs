use core::f64;
use std::mem;
use std::ops::Deref;

use cfd_rs::finite_volume::base::Field;
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
    config::{CaseConfig, GeometryConfig, MeshingConfig, Schemes},
    discretizations::{
        convection::ConvectionScheme, divergence::DivergenceScheme, laplacian::LaplacianScheme,
        time_schemes::TimeIntegration,
    },
    equation::{Dimension, Variable},
    gradients::{GradientConfig, GradientScheme},
    interpolations::GradientInterpConfig,
};
use nalgebra::{Point2, Vector2};

fn simple(geometry: GeometryConfig) -> CaseConfig {
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

    let output = OutputConfig {
        control: OutputControl::Iteration(1),
        directory: "./target/exports".to_string(),
    };

    let mut bc_fields = HashMap::new();
    let bc;
    match geometry.meshing {
        MeshingConfig::AdvancingFront { .. } => {
            bc = vec![
                BoundaryCondition::Dirichlet(BoundaryValue::Vector2(Vector2::new(0., 0.))),
                BoundaryCondition::Dirichlet(BoundaryValue::Vector2(Vector2::new(0., 0.))),
                BoundaryCondition::Dirichlet(BoundaryValue::Vector2(Vector2::new(1., 0.))),
                BoundaryCondition::Dirichlet(BoundaryValue::Vector2(Vector2::new(0., 0.))),
            ];
        }
        _ => {
            bc = vec![
                BoundaryCondition::Dirichlet(BoundaryValue::Vector2(Vector2::new(0., 0.))),
                BoundaryCondition::Dirichlet(BoundaryValue::Vector2(Vector2::new(0., 0.))),
                BoundaryCondition::Dirichlet(BoundaryValue::Vector2(Vector2::new(0., 1.))),
                BoundaryCondition::Dirichlet(BoundaryValue::Vector2(Vector2::new(0., 0.))),
            ];
        }
    }
    bc_fields.insert(Variable::new("U".to_string(), Dimension::Vector2), bc);
    let bc = vec![
        BoundaryCondition::Neumann(BoundaryValue::Scalar(0.)),
        BoundaryCondition::Neumann(BoundaryValue::Scalar(0.)),
        BoundaryCondition::Neumann(BoundaryValue::Scalar(0.)),
        BoundaryCondition::Neumann(BoundaryValue::Scalar(0.)),
    ];
    bc_fields.insert(Variable::new("P".to_string(), Dimension::Scalar), bc);
    let bc = vec![
        BoundaryCondition::Neumann(BoundaryValue::Scalar(0.)),
        BoundaryCondition::Neumann(BoundaryValue::Scalar(0.)),
        BoundaryCondition::Neumann(BoundaryValue::Scalar(0.)),
        BoundaryCondition::Neumann(BoundaryValue::Scalar(0.)),
    ];
    bc_fields.insert(Variable::new("div(U)".to_string(), Dimension::Scalar), bc);
    let bc = vec![
        BoundaryCondition::Neumann(BoundaryValue::Vector2(Vector2::new(0., 0.))),
        BoundaryCondition::Neumann(BoundaryValue::Vector2(Vector2::new(0., 0.))),
        BoundaryCondition::Neumann(BoundaryValue::Vector2(Vector2::new(0., 0.))),
        BoundaryCondition::Neumann(BoundaryValue::Vector2(Vector2::new(0., 0.))),
    ];
    bc_fields.insert(Variable::new("grad(P)".to_string(), Dimension::Vector2), bc);
    let bc_fields = FieldsBoundaryConditions::new(bc_fields);

    let mut initial_fields = HashMap::new();
    initial_fields.insert(
        Variable::new("U".to_string(), Dimension::Vector2),
        InitFunc::Vector2(constant, constant),
    );
    initial_fields.insert(
        Variable::new("P".to_string(), Dimension::Scalar),
        InitFunc::Scalar(constant),
    );
    initial_fields.insert(
        Variable::new("div(U)".to_string(), Dimension::Scalar),
        InitFunc::Scalar(constant),
    );
    initial_fields.insert(
        Variable::new("grad(P)".to_string(), Dimension::Vector2),
        InitFunc::Vector2(constant, constant),
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

pub fn constant_1(_point: &Point2<f64>) -> f64 {
    1.
}

pub fn sinusoidal_1d(point: &Point2<f64>) -> f64 {
    if point.x < 0.25 {
        (point.x * f64::consts::PI * 4.).sin()
    } else {
        0.
    }
}

pub fn u_x(point: &Point2<f64>) -> f64 {
    let pi = f64::consts::PI;
    (point.x * 2. * pi).sin() * (2. * point.y * pi).sin()
}

pub fn u_y(point: &Point2<f64>) -> f64 {
    let pi = f64::consts::PI;
    (point.x * pi * 2.).cos() * (point.y * pi * 2.).cos()
}

pub fn solution(point: &Point2<f64>) -> f64 {
    let pi = f64::consts::PI;
    let s = Point2::new((pi * point.x).sin(), (pi * point.y).sin());
    s.x.powi(2) * s.y.powi(2)
}

pub fn divergence_source_x(point: &Point2<f64>) -> f64 {
    let pi = f64::consts::PI;
    let s2 = Point2::new((2. * pi * point.x).sin(), (2. * pi * point.y).sin());
    let s = Point2::new((pi * point.x).sin(), (pi * point.y).sin());
    1. * pi * s2.x * s.y.powi(2) + 0.3 * pi * s2.y * s.x.powi(2) + s.x.powi(2) * s.y.powi(2)
}

// pub fn divergence_source_y(point: &Point2<f64>) -> f64 {
//     let pi = f64::consts::PI;
//     let c = Point2::new((2.*pi*point.x).cos(), (2.*pi*point.y).cos());
//     let s = Point2::new((2.*pi*point.x).sin(), (2.*pi*point.y).sin());
//     pi*s.x*c.y + 0.5*pi*c.x*s.y + c.x*c.y
// }

fn main() {
    let geometry = GeometryConfig {
        import_path: Some("../meshes/mesh3.cfd".to_string()),
        meshing: MeshingConfig::AdvancingFront { element_size: 0.01 },
    };
    let geometry = GeometryConfig {
        import_path: None,
        meshing: MeshingConfig::Cartesian {
            length: Vector2::new(1., 1.),
            n_elements: Vector2::new(50, 50),
        },
    };

    let mut case = SimpleCase::new(simple(geometry));

    case.export_cell_centered().unwrap();

    for _ in 0..10 {
        case.next_step();
        // {
        //     let temp = case.field(&Variable::new("Phi".to_string(), Dimension::Scalar)).expect("");
        //     let field = if let Field::Scalar(values) = temp.deref() {
        //         values
        //     } else {
        //         panic!();
        //     };
        //     println!("{:?}", field.face_values());
        // }
        // if case.step() % 10 == 0 {
        //     case.export().unwrap();
        // }
        case.export().unwrap();
        println!("{:?}", case.time());
    }
}
