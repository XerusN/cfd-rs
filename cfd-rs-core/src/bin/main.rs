// use core::f64;
// use std::mem;
// use std::ops::Deref;

// use cfd_rs::finite_volume::base::Field;
// use cfd_rs::finite_volume::boundary::BoundaryValue;
// use cfd_rs::finite_volume::case::burger::BurgerCase;
// use cfd_rs::finite_volume::case::convection_diffusion_test::ConvectionDiffusionCase;
// use cfd_rs::finite_volume::case::convection_test::ConvectionCase;
// use cfd_rs::finite_volume::case::diffusion_test::DiffusionCase;
// use cfd_rs::finite_volume::case::divergence_test::DivergenceCase;
// use cfd_rs::finite_volume::case::divergence_test2::DivergenceCase2;
// use cfd_rs::finite_volume::case::simple::SimpleCase;
// use cfd_rs::finite_volume::config::{InitFunc, OutputConfig};
// use cfd_rs_utils::control::OutputControl;
// use hashbrown::HashMap;

// use cfd_rs::finite_volume::case::{self, poisson};
// use cfd_rs::finite_volume::{
//     boundary::{BoundaryCondition, FieldsBoundaryConditions},
//     case::poisson::PoissonCase,
//     case::Case,
//     config::{CaseConfig, GeometryConfig, Schemes},
//     discretizations::{
//         convection::ConvectionScheme, divergence::DivergenceScheme, laplacian::LaplacianScheme,
//         time_schemes::TimeIntegration,
//     },
//     equation::{Dimension, Variable},
//     gradients::{GradientConfig, GradientScheme},
//     interpolations::GradientInterpConfig,
// };
// use nalgebra::{Point2, Vector2};

// fn poisson() -> CaseConfig {
//     let schemes = Schemes {
//         transient: TimeIntegration::ForwardEuler,
//         convection: ConvectionScheme::UpwindSecondOrder,
//         laplacian: LaplacianScheme::OrthogonalCorrection,
//         divergence: DivergenceScheme::Basic,
//         gradients: GradientConfig {
//             scheme: GradientScheme::GreenGaussCompact,
//             interp: GradientInterpConfig::AveragedCorrected,
//         },
//     };

//     // let geometry = GeometryConfig {
//     //     import_path: Some("./target/exports/mesh.cfd".to_string()),
//     //     element_size: 0.01,
//     // };
//     let geometry = GeometryConfig {
//         import_path: None,
//         element_size: 0.01,
//     };

//     let output = OutputConfig {
//         control: OutputControl::Iteration(1),
//         directory: "./target/exports".to_string(),
//     };

//     let mut bc_fields = HashMap::new();
//     let bc = vec![
//         BoundaryCondition::Dirichlet(BoundaryValue::Scalar(0.)),
//         BoundaryCondition::Dirichlet(BoundaryValue::Scalar(0.)),
//         BoundaryCondition::Dirichlet(BoundaryValue::Scalar(200.)),
//         BoundaryCondition::Dirichlet(BoundaryValue::Scalar(300.)),
//     ];
//     bc_fields.insert(Variable::new("T".to_string(), Dimension::Scalar), bc);
//     let bc_fields = FieldsBoundaryConditions::new(bc_fields);

//     let mut initial_fields = HashMap::new();
//     initial_fields.insert(
//         Variable::new("T".to_string(), Dimension::Scalar),
//         InitFunc::Scalar(constant),
//     );

//     CaseConfig {
//         schemes,
//         geometry,
//         bc: bc_fields,
//         output,
//         initial_fields,
//     }
// }

// fn simple() -> CaseConfig {
//     let schemes = Schemes {
//         transient: TimeIntegration::ForwardEuler,
//         convection: ConvectionScheme::UpwindSecondOrder,
//         laplacian: LaplacianScheme::OrthogonalCorrection,
//         divergence: DivergenceScheme::Basic,
//         gradients: GradientConfig {
//             scheme: GradientScheme::GreenGaussCompact,
//             interp: GradientInterpConfig::AveragedCorrected,
//         },
//     };

//     let geometry = GeometryConfig {
//         import_path: Some("./target/exports/mesh.cfd".to_string()),
//         element_size: 0.01,
//     };
//     // let geometry = GeometryConfig {
//     //     import_path: None,
//     //     element_size: 0.02,
//     // };

//     let output = OutputConfig {
//         control: OutputControl::Iteration(1),
//         directory: "./target/exports".to_string(),
//     };

//     let mut bc_fields = HashMap::new();
//     let bc = vec![
//         BoundaryCondition::Dirichlet(BoundaryValue::Vector2(Vector2::new(0., 0.))),
//         BoundaryCondition::Dirichlet(BoundaryValue::Vector2(Vector2::new(0., 0.))),
//         BoundaryCondition::Dirichlet(BoundaryValue::Vector2(Vector2::new(1., 0.))),
//         BoundaryCondition::Dirichlet(BoundaryValue::Vector2(Vector2::new(0., 0.))),
//     ];
//     bc_fields.insert(Variable::new("U".to_string(), Dimension::Vector2), bc);
//     let bc = vec![
//         BoundaryCondition::Neumann(BoundaryValue::Scalar(0.)),
//         BoundaryCondition::Neumann(BoundaryValue::Scalar(0.)),
//         BoundaryCondition::Neumann(BoundaryValue::Scalar(0.)),
//         BoundaryCondition::Neumann(BoundaryValue::Scalar(0.)),
//     ];
//     bc_fields.insert(Variable::new("P".to_string(), Dimension::Scalar), bc);
//     let bc = vec![
//         BoundaryCondition::Neumann(BoundaryValue::Scalar(0.)),
//         BoundaryCondition::Neumann(BoundaryValue::Scalar(0.)),
//         BoundaryCondition::Neumann(BoundaryValue::Scalar(0.)),
//         BoundaryCondition::Neumann(BoundaryValue::Scalar(0.)),
//     ];
//     bc_fields.insert(Variable::new("div(U)".to_string(), Dimension::Scalar), bc);
//     let bc = vec![
//         BoundaryCondition::Neumann(BoundaryValue::Vector2(Vector2::new(0., 0.))),
//         BoundaryCondition::Neumann(BoundaryValue::Vector2(Vector2::new(0., 0.))),
//         BoundaryCondition::Neumann(BoundaryValue::Vector2(Vector2::new(0., 0.))),
//         BoundaryCondition::Neumann(BoundaryValue::Vector2(Vector2::new(0., 0.))),
//     ];
//     bc_fields.insert(Variable::new("grad(P)".to_string(), Dimension::Vector2), bc);
//     let bc_fields = FieldsBoundaryConditions::new(bc_fields);

//     let mut initial_fields = HashMap::new();
//     initial_fields.insert(
//         Variable::new("U".to_string(), Dimension::Vector2),
//         InitFunc::Vector2(constant, constant),
//     );
//     initial_fields.insert(
//         Variable::new("P".to_string(), Dimension::Scalar),
//         InitFunc::Scalar(constant),
//     );
//     initial_fields.insert(
//         Variable::new("div(U)".to_string(), Dimension::Scalar),
//         InitFunc::Scalar(constant),
//     );
//     initial_fields.insert(
//         Variable::new("grad(P)".to_string(), Dimension::Vector2),
//         InitFunc::Vector2(constant, constant),
//     );

//     CaseConfig {
//         schemes,
//         geometry,
//         bc: bc_fields,
//         output,
//         initial_fields,
//     }
// }

// fn pipe() -> CaseConfig {
//     let schemes = Schemes {
//         transient: TimeIntegration::ForwardEuler,
//         convection: ConvectionScheme::UpwindSecondOrder,
//         laplacian: LaplacianScheme::OrthogonalCorrection,
//         divergence: DivergenceScheme::Basic,
//         gradients: GradientConfig {
//             scheme: GradientScheme::GreenGaussCompact,
//             interp: GradientInterpConfig::AveragedCorrected,
//         },
//     };

//     // let geometry = GeometryConfig {
//     //     import_path: Some("./target/exports/mesh.cfd".to_string()),
//     //     element_size: 0.01,
//     // };
//     let geometry = GeometryConfig {
//         import_path: None,
//         element_size: 0.01,
//     };

//     let output = OutputConfig {
//         control: OutputControl::Iteration(1),
//         directory: "./target/exports".to_string(),
//     };

//     let mut bc_fields = HashMap::new();
//     let bc = vec![
//         BoundaryCondition::Dirichlet(BoundaryValue::Vector2(Vector2::new(1., 0.))),
//         BoundaryCondition::Dirichlet(BoundaryValue::Vector2(Vector2::new(0., 0.))),
//         BoundaryCondition::Neumann(BoundaryValue::Vector2(Vector2::new(0., 0.))),
//         BoundaryCondition::Dirichlet(BoundaryValue::Vector2(Vector2::new(0., 0.))),
//     ];
//     bc_fields.insert(Variable::new("U".to_string(), Dimension::Vector2), bc);
//     let bc = vec![
//         BoundaryCondition::Neumann(BoundaryValue::Scalar(0.)),
//         BoundaryCondition::Neumann(BoundaryValue::Scalar(0.)),
//         BoundaryCondition::Dirichlet(BoundaryValue::Scalar(0.)),
//         BoundaryCondition::Neumann(BoundaryValue::Scalar(0.)),
//     ];
//     bc_fields.insert(Variable::new("P".to_string(), Dimension::Scalar), bc);
//     let bc = vec![
//         BoundaryCondition::Neumann(BoundaryValue::Scalar(0.)),
//         BoundaryCondition::Neumann(BoundaryValue::Scalar(0.)),
//         BoundaryCondition::Neumann(BoundaryValue::Scalar(0.)),
//         BoundaryCondition::Neumann(BoundaryValue::Scalar(0.)),
//     ];
//     bc_fields.insert(Variable::new("div(U)".to_string(), Dimension::Scalar), bc);
//     let bc_fields = FieldsBoundaryConditions::new(bc_fields);

//     let mut initial_fields = HashMap::new();
//     initial_fields.insert(
//         Variable::new("U".to_string(), Dimension::Vector2),
//         InitFunc::Vector2(constant_1, constant),
//     );
//     initial_fields.insert(
//         Variable::new("P".to_string(), Dimension::Scalar),
//         InitFunc::Scalar(constant),
//     );
//     initial_fields.insert(
//         Variable::new("div(U)".to_string(), Dimension::Scalar),
//         InitFunc::Scalar(constant),
//     );

//     CaseConfig {
//         schemes,
//         geometry,
//         bc: bc_fields,
//         output,
//         initial_fields,
//     }
// }

// fn convection_setup() -> CaseConfig {
//     let schemes = Schemes {
//         transient: TimeIntegration::ForwardEuler,
//         convection: ConvectionScheme::UpwindSecondOrder,
//         laplacian: LaplacianScheme::OrthogonalCorrection,
//         divergence: DivergenceScheme::Basic,
//         gradients: GradientConfig {
//             scheme: GradientScheme::GreenGaussCompact,
//             interp: GradientInterpConfig::AveragedCorrected,
//         },
//     };

//     let geometry = GeometryConfig {
//         import_path: None,
//         element_size: 0.01,
//     };

//     let output = OutputConfig {
//         control: OutputControl::Iteration(1),
//         directory: "./target/exports".to_string(),
//     };

//     let mut bc_fields = HashMap::new();
//     let bc = vec![
//         BoundaryCondition::Dirichlet(BoundaryValue::Scalar(0.)),
//         BoundaryCondition::Dirichlet(BoundaryValue::Scalar(0.)),
//         BoundaryCondition::Neumann(BoundaryValue::Scalar(0.)),
//     ];
//     bc_fields.insert(Variable::new("Phi".to_string(), Dimension::Scalar), bc);

//     let bc = vec![
//         BoundaryCondition::Neumann(BoundaryValue::Vector2(Vector2::new(0., 0.))),
//         BoundaryCondition::Neumann(BoundaryValue::Vector2(Vector2::new(0., 0.))),
//         BoundaryCondition::Neumann(BoundaryValue::Vector2(Vector2::new(0., 0.))),
//     ];
//     bc_fields.insert(Variable::new("Speed".to_string(), Dimension::Vector2), bc);
//     let bc_fields = FieldsBoundaryConditions::new(bc_fields);

//     let mut initial_fields = HashMap::new();
//     initial_fields.insert(
//         Variable::new("Phi".to_string(), Dimension::Scalar),
//         InitFunc::Scalar(sinusoidal_1d),
//     );
//     initial_fields.insert(
//         Variable::new("Speed".to_string(), Dimension::Vector2),
//         InitFunc::Vector2(constant_1, constant),
//     );

//     CaseConfig {
//         schemes,
//         geometry,
//         bc: bc_fields,
//         output,
//         initial_fields,
//     }
// }

// fn convection_setup_2d() -> CaseConfig {
//     let schemes = Schemes {
//         transient: TimeIntegration::ForwardEuler,
//         convection: ConvectionScheme::UpwindSecondOrder,
//         laplacian: LaplacianScheme::OrthogonalCorrection,
//         divergence: DivergenceScheme::Basic,
//         gradients: GradientConfig {
//             scheme: GradientScheme::GreenGaussCompact,
//             interp: GradientInterpConfig::AveragedCorrected,
//         },
//     };

//     let geometry = GeometryConfig {
//         import_path: Some("./target/exports/mesh2.cfd".to_string()),
//         element_size: 0.01,
//     };

//     let output = OutputConfig {
//         control: OutputControl::Iteration(1),
//         directory: "./target/exports".to_string(),
//     };

//     let mut bc_fields = HashMap::new();
//     let bc = vec![
//         BoundaryCondition::Neumann(BoundaryValue::Scalar(0.)),
//         BoundaryCondition::Dirichlet(BoundaryValue::Scalar(0.)),
//         BoundaryCondition::Neumann(BoundaryValue::Scalar(0.)),
//         BoundaryCondition::Dirichlet(BoundaryValue::Scalar(0.)),
//     ];
//     bc_fields.insert(Variable::new("Phi".to_string(), Dimension::Scalar), bc);

//     let bc = vec![
//         BoundaryCondition::Neumann(BoundaryValue::Vector2(Vector2::new(0., 0.))),
//         BoundaryCondition::Neumann(BoundaryValue::Vector2(Vector2::new(0., 0.))),
//         BoundaryCondition::Neumann(BoundaryValue::Vector2(Vector2::new(0., 0.))),
//         BoundaryCondition::Neumann(BoundaryValue::Vector2(Vector2::new(0., 0.))),
//     ];
//     bc_fields.insert(Variable::new("Speed".to_string(), Dimension::Vector2), bc);
//     let bc_fields = FieldsBoundaryConditions::new(bc_fields);

//     let mut initial_fields = HashMap::new();
//     initial_fields.insert(
//         Variable::new("Phi".to_string(), Dimension::Scalar),
//         InitFunc::Scalar(sinusoidal_1d),
//     );
//     initial_fields.insert(
//         Variable::new("Speed".to_string(), Dimension::Vector2),
//         InitFunc::Vector2(constant_1, constant),
//     );

//     CaseConfig {
//         schemes,
//         geometry,
//         bc: bc_fields,
//         output,
//         initial_fields,
//     }
// }

// fn burger_setup_2d() -> CaseConfig {
//     let schemes = Schemes {
//         transient: TimeIntegration::ForwardEuler,
//         convection: ConvectionScheme::UpwindSecondOrder,
//         laplacian: LaplacianScheme::OrthogonalCorrection,
//         divergence: DivergenceScheme::Basic,
//         gradients: GradientConfig {
//             scheme: GradientScheme::GreenGaussCompact,
//             interp: GradientInterpConfig::AveragedCorrected,
//         },
//     };

//     let geometry = GeometryConfig {
//         import_path: Some("./target/exports/mesh2.cfd".to_string()),
//         element_size: 0.01,
//     };

//     let output = OutputConfig {
//         control: OutputControl::Iteration(1),
//         directory: "./target/exports".to_string(),
//     };

//     let mut bc_fields = HashMap::new();

//     let bc = vec![
//         BoundaryCondition::Neumann(BoundaryValue::Vector2(Vector2::new(0., 0.))),
//         BoundaryCondition::Neumann(BoundaryValue::Vector2(Vector2::new(0., 0.))),
//         BoundaryCondition::Neumann(BoundaryValue::Vector2(Vector2::new(0., 0.))),
//         BoundaryCondition::Dirichlet(BoundaryValue::Vector2(Vector2::new(0., 0.))),
//     ];
//     bc_fields.insert(Variable::new("Speed".to_string(), Dimension::Vector2), bc);
//     let bc_fields = FieldsBoundaryConditions::new(bc_fields);

//     let mut initial_fields = HashMap::new();
//     initial_fields.insert(
//         Variable::new("Speed".to_string(), Dimension::Vector2),
//         InitFunc::Vector2(sinusoidal_1d, constant),
//     );

//     CaseConfig {
//         schemes,
//         geometry,
//         bc: bc_fields,
//         output,
//         initial_fields,
//     }
// }

// fn diffusion_setup() -> CaseConfig {
//     let schemes = Schemes {
//         transient: TimeIntegration::ForwardEuler,
//         convection: ConvectionScheme::UpwindSecondOrder,
//         laplacian: LaplacianScheme::OrthogonalCorrection,
//         divergence: DivergenceScheme::Basic,
//         gradients: GradientConfig {
//             scheme: GradientScheme::GreenGaussCompact,
//             interp: GradientInterpConfig::AveragedCorrected,
//         },
//     };

//     let geometry = GeometryConfig {
//         import_path: None,
//         element_size: 0.01,
//     };

//     let output = OutputConfig {
//         control: OutputControl::Iteration(1),
//         directory: "./target/exports".to_string(),
//     };

//     let mut bc_fields = HashMap::new();
//     let bc = vec![
//         BoundaryCondition::Dirichlet(BoundaryValue::Scalar(0.)),
//         BoundaryCondition::Dirichlet(BoundaryValue::Scalar(0.)),
//         BoundaryCondition::Neumann(BoundaryValue::Scalar(0.)),
//     ];
//     bc_fields.insert(Variable::new("Phi".to_string(), Dimension::Scalar), bc);

//     let bc_fields = FieldsBoundaryConditions::new(bc_fields);

//     let mut initial_fields = HashMap::new();
//     initial_fields.insert(
//         Variable::new("Phi".to_string(), Dimension::Scalar),
//         InitFunc::Scalar(sinusoidal_1d),
//     );

//     CaseConfig {
//         schemes,
//         geometry,
//         bc: bc_fields,
//         output,
//         initial_fields,
//     }
// }

// fn diffusion_setup_2d() -> CaseConfig {
//     let schemes = Schemes {
//         transient: TimeIntegration::ForwardEuler,
//         convection: ConvectionScheme::UpwindSecondOrder,
//         laplacian: LaplacianScheme::OrthogonalCorrection,
//         divergence: DivergenceScheme::Basic,
//         gradients: GradientConfig {
//             scheme: GradientScheme::GreenGaussCompact,
//             interp: GradientInterpConfig::AveragedCorrected,
//         },
//     };

//     let geometry = GeometryConfig {
//         import_path: Some("./target/exports/mesh2.cfd".to_string()),
//         element_size: 0.01,
//     };

//     let output = OutputConfig {
//         control: OutputControl::Iteration(1),
//         directory: "./target/exports".to_string(),
//     };

//     let mut bc_fields = HashMap::new();
//     let bc = vec![
//         BoundaryCondition::Dirichlet(BoundaryValue::Scalar(1.)),
//         BoundaryCondition::Dirichlet(BoundaryValue::Scalar(2.)),
//         BoundaryCondition::Dirichlet(BoundaryValue::Scalar(0.)),
//         BoundaryCondition::Dirichlet(BoundaryValue::Scalar(0.)),
//     ];
//     bc_fields.insert(Variable::new("Phi".to_string(), Dimension::Scalar), bc);

//     let bc_fields = FieldsBoundaryConditions::new(bc_fields);

//     let mut initial_fields = HashMap::new();
//     initial_fields.insert(
//         Variable::new("Phi".to_string(), Dimension::Scalar),
//         InitFunc::Scalar(constant),
//     );

//     CaseConfig {
//         schemes,
//         geometry,
//         bc: bc_fields,
//         output,
//         initial_fields,
//     }
// }

// fn divergence_setup_2() -> CaseConfig {
//     let schemes = Schemes {
//         transient: TimeIntegration::ForwardEuler,
//         convection: ConvectionScheme::UpwindSecondOrder,
//         laplacian: LaplacianScheme::OrthogonalCorrection,
//         divergence: DivergenceScheme::Basic,
//         gradients: GradientConfig {
//             scheme: GradientScheme::GreenGaussCompact,
//             interp: GradientInterpConfig::AveragedCorrected,
//         },
//     };

//     let geometry = GeometryConfig {
//         import_path: None,
//         element_size: 0.01,
//     };

//     let output = OutputConfig {
//         control: OutputControl::Iteration(1),
//         directory: "./target/exports".to_string(),
//     };

//     let mut bc_fields = HashMap::new();
//     let bc = vec![
//         BoundaryCondition::Dirichlet(BoundaryValue::Vector2(Vector2::new(100., 100.))),
//         BoundaryCondition::Dirichlet(BoundaryValue::Vector2(Vector2::new(100., 100.))),
//         BoundaryCondition::Dirichlet(BoundaryValue::Vector2(Vector2::new(100., 100.))),
//         BoundaryCondition::Dirichlet(BoundaryValue::Vector2(Vector2::new(100., 100.))),
//     ];
//     bc_fields.insert(Variable::new("U".to_string(), Dimension::Vector2), bc);
//     let bc = vec![
//         BoundaryCondition::Dirichlet(BoundaryValue::Vector2(Vector2::new(0., 0.))),
//         BoundaryCondition::Dirichlet(BoundaryValue::Vector2(Vector2::new(0., 0.))),
//         BoundaryCondition::Dirichlet(BoundaryValue::Vector2(Vector2::new(0., 0.))),
//         BoundaryCondition::Dirichlet(BoundaryValue::Vector2(Vector2::new(0., 0.))),
//     ];
//     bc_fields.insert(Variable::new("f".to_string(), Dimension::Vector2), bc);

//     let bc_fields = FieldsBoundaryConditions::new(bc_fields);

//     let mut initial_fields = HashMap::new();
//     initial_fields.insert(
//         Variable::new("f".to_string(), Dimension::Vector2),
//         InitFunc::Vector2(solution, solution),
//     );
//     initial_fields.insert(
//         Variable::new("U".to_string(), Dimension::Vector2),
//         InitFunc::Vector2(constant_1, constant_1),
//     );

//     CaseConfig {
//         schemes,
//         geometry,
//         bc: bc_fields,
//         output,
//         initial_fields,
//     }
// }

// fn test_div(point: &Point2<f64>) -> f64 {
//     let (x, y) = (point.x, point.y);
//     x*(x-1.)*y*(y-1.)
// }

// fn divergence_setup_2d() -> CaseConfig {
//     let schemes = Schemes {
//         transient: TimeIntegration::ForwardEuler,
//         convection: ConvectionScheme::UpwindSecondOrder,
//         laplacian: LaplacianScheme::OrthogonalCorrection,
//         divergence: DivergenceScheme::Basic,
//         gradients: GradientConfig {
//             scheme: GradientScheme::GreenGaussCompact,
//             interp: GradientInterpConfig::AveragedCorrected,
//         },
//     };

//     let geometry = GeometryConfig {
//         import_path: None,
//         element_size: 0.01,
//     };

//     let output = OutputConfig {
//         control: OutputControl::Iteration(1),
//         directory: "./target/exports".to_string(),
//     };

//     let mut bc_fields = HashMap::new();
//     let bc = vec![
//         BoundaryCondition::Dirichlet(BoundaryValue::Vector2(Vector2::new(0., 0.))),
//         BoundaryCondition::Dirichlet(BoundaryValue::Vector2(Vector2::new(0., 0.))),
//         BoundaryCondition::Dirichlet(BoundaryValue::Vector2(Vector2::new(0., 0.))),
//         BoundaryCondition::Dirichlet(BoundaryValue::Vector2(Vector2::new(0., 0.))),
//     ];
//     bc_fields.insert(Variable::new("U".to_string(), Dimension::Vector2), bc);
//     let bc = vec![
//         BoundaryCondition::Dirichlet(BoundaryValue::Vector2(Vector2::new(0., 0.))),
//         BoundaryCondition::Dirichlet(BoundaryValue::Vector2(Vector2::new(0., 0.))),
//         BoundaryCondition::Dirichlet(BoundaryValue::Vector2(Vector2::new(0., 0.))),
//         BoundaryCondition::Dirichlet(BoundaryValue::Vector2(Vector2::new(0., 0.))),
//     ];
//     bc_fields.insert(Variable::new("f".to_string(), Dimension::Vector2), bc);

//     let bc_fields = FieldsBoundaryConditions::new(bc_fields);

//     let mut initial_fields = HashMap::new();
//     initial_fields.insert(
//         Variable::new("f".to_string(), Dimension::Vector2),
//         InitFunc::Vector2(solution, solution),
//     );
//     initial_fields.insert(
//         Variable::new("U".to_string(), Dimension::Vector2),
//         InitFunc::Vector2(constant_1, constant_1),
//     );

//     CaseConfig {
//         schemes,
//         geometry,
//         bc: bc_fields,
//         output,
//         initial_fields,
//     }
// }

// pub fn constant(_point: &Point2<f64>) -> f64 {
//     0.
// }

// pub fn constant_1(_point: &Point2<f64>) -> f64 {
//     1.
// }

// pub fn sinusoidal_1d(point: &Point2<f64>) -> f64 {
//     if point.x < 0.25 {
//         (point.x * f64::consts::PI * 4.).sin()
//     } else {
//         0.
//     }
// }

// pub fn u_x(point: &Point2<f64>) -> f64 {
//     let pi = f64::consts::PI;
//     (point.x * 2. * pi).sin() * (2. * point.y * pi).sin()
// }

// pub fn u_y(point: &Point2<f64>) -> f64 {
//     let pi = f64::consts::PI;
//     (point.x * pi * 2.).cos() * (point.y * pi * 2.).cos()
// }

// pub fn solution(point: &Point2<f64>) -> f64 {
//     let pi = f64::consts::PI;
//     let s = Point2::new((pi * point.x).sin(), (pi * point.y).sin());
//     s.x.powi(2) * s.y.powi(2)
// }

// pub fn divergence_source_x(point: &Point2<f64>) -> f64 {
//     let pi = f64::consts::PI;
//     let s2 = Point2::new((2. * pi * point.x).sin(), (2. * pi * point.y).sin());
//     let s = Point2::new((pi * point.x).sin(), (pi * point.y).sin());
//     1. * pi * s2.x * s.y.powi(2) + 0.3 * pi * s2.y * s.x.powi(2) + s.x.powi(2) * s.y.powi(2)
// }

// // pub fn divergence_source_y(point: &Point2<f64>) -> f64 {
// //     let pi = f64::consts::PI;
// //     let c = Point2::new((2.*pi*point.x).cos(), (2.*pi*point.y).cos());
// //     let s = Point2::new((2.*pi*point.x).sin(), (2.*pi*point.y).sin());
// //     pi*s.x*c.y + 0.5*pi*c.x*s.y + c.x*c.y
// // }

// fn main() {
//     //let mut case = ConvectionCase::new(convection_setup_2d());
//     // let mut case = PoissonCase::new(poisson());
//     // let mut case = SimpleCase::new(simple());
//     // let mut case = ConvectionDiffusionCase::new(convection_setup_2d());
//     //let mut case = BurgerCase::new(burger_setup_2d());
//     let mut case = SimpleCase::new(simple());
//     //let mut case = DiffusionCase::new(convection_setup_2d());
//     // let mut case = DivergenceCase::new(divergence_setup_2d());
//     // let mut case = SimpleCase::new(simple());

//     // //println!("{:?}", case.mesh().cells().iter().map(|cell| cell.volume()).collect::<Vec<f64>>());
//     // {
//     //     let temp = case.field(&Variable::new("Phi".to_string(), Dimension::Scalar)).expect("");
//     //     let field = if let Field::Scalar(values) = temp.deref() {
//     //         values
//     //     } else {
//     //         panic!();
//     //     };
//     //     println!("{:?}", field.face_values());
//     // }

//     case.export().unwrap();

//     for _ in 0..10 {
//         case.next_step();
//         // {
//         //     let temp = case.field(&Variable::new("Phi".to_string(), Dimension::Scalar)).expect("");
//         //     let field = if let Field::Scalar(values) = temp.deref() {
//         //         values
//         //     } else {
//         //         panic!();
//         //     };
//         //     println!("{:?}", field.face_values());
//         // }
//         // if case.step() % 10 == 0 {
//         //     case.export().unwrap();
//         // }
//         case.export().unwrap();
//         println!("{:?}", case.time());
//     }
// }
fn main() {}
