use std::path::{Path, PathBuf};

use super::super::super::examples_default_logger;
use cfd_rs::{
    finite_volume::{
        boundary::{BoundaryCondition, BoundaryValue},
        config::{
            GeometryConfig, GradientConfig, InitFunctionsTrait, MeshingConfig, SchemesConfig,
        },
        mesh::mesh,
        solvers::{
            poisson::{Config, MainEquations, MainVariables, PoissonCase, PoissonUserFunctions},
            EquationsEnumConfigTrait, VariablesEnumConfigTrait,
        },
    },
    post_processing::{
        grid::{GridPostprocConfig, GridPostprocType, VariablesExtracted},
        OutputConfig,
    },
};
use cfd_rs_utils::control::{clear_dump_folder, Frequency, OutputControl};
use log::{debug, info, trace, warn, LevelFilter};
use nalgebra::{Point2, Vector2};

fn file_directory() -> PathBuf {
    let mut directory = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    directory.pop();
    directory.push(file!());
    directory.pop();
    directory
}

#[allow(dead_code)]
#[derive(Debug)]
struct InitFunctions {
    var: MainVariables,
}

impl InitFunctionsTrait for InitFunctions {
    fn init_x(&self, point: &Point2<f64>) -> f64 {
        match &self.var {
            &MainVariables::T => point.x + point.y,
            _ => 0.,
        }
    }
}

#[allow(dead_code)]
struct VariablesEnumConfig;

impl VariablesEnumConfigTrait<MainVariables, InitFunctions> for VariablesEnumConfig {
    fn gradient_configs(&self, variable_enum: &MainVariables) -> GradientConfig {
        match variable_enum {
            _ => GradientConfig::default(),
        }
    }

    fn initial_fields(&self, variable_enum: &MainVariables) -> InitFunctions {
        InitFunctions {
            var: variable_enum.clone(),
        }
    }

    fn boundary_conditions(&self, variable_enum: &MainVariables) -> Vec<BoundaryCondition> {
        match variable_enum {
            MainVariables::T => {
                vec![
                    // bot | cart: left
                    BoundaryCondition::Dirichlet(BoundaryValue::Scalar(1.)),
                    // right | cart: bot
                    BoundaryCondition::Neumann(BoundaryValue::Scalar(0.)),
                    // top | cart: right
                    BoundaryCondition::Dirichlet(BoundaryValue::Scalar(2.)),
                    // left | cart: top
                    BoundaryCondition::Neumann(BoundaryValue::Scalar(0.)),
                ]
            }
            MainVariables::Laplacian => {
                vec![
                    BoundaryCondition::Neumann(BoundaryValue::Scalar(0.)),
                    BoundaryCondition::Neumann(BoundaryValue::Scalar(0.)),
                    BoundaryCondition::Neumann(BoundaryValue::Scalar(0.)),
                    BoundaryCondition::Neumann(BoundaryValue::Scalar(0.)),
                ]
            }
            MainVariables::Grad => {
                vec![
                    BoundaryCondition::Neumann(BoundaryValue::Vector2(Vector2::new(0., 0.))),
                    BoundaryCondition::Neumann(BoundaryValue::Vector2(Vector2::new(0., 0.))),
                    BoundaryCondition::Neumann(BoundaryValue::Vector2(Vector2::new(0., 0.))),
                    BoundaryCondition::Neumann(BoundaryValue::Vector2(Vector2::new(0., 0.))),
                ]
            }
        }
    }
}

#[allow(dead_code)]
struct EquationsEnumConfig;

impl EquationsEnumConfigTrait<MainEquations> for EquationsEnumConfig {
    fn schemes_configs(&self, equations_enum: &MainEquations) -> SchemesConfig {
        match equations_enum {
            _ => SchemesConfig::default(),
        }
    }
}

#[allow(dead_code)]
struct UserFunctions;

impl PoissonUserFunctions for UserFunctions {}

#[allow(dead_code)]
fn poisson() -> (Config, EquationsEnumConfig, VariablesEnumConfig) {
    let directory = file_directory();
    let dump = directory.join("dump");

    let geometry = GeometryConfig {
        import_path: Some("../meshes/mesh4.cfd".to_string()),
        meshing: MeshingConfig::AdvancingFront { element_size: 0.01 },
    };

    let postproc = vec![OutputConfig::Grid(GridPostprocConfig::new(
        "main".to_owned(),
        GridPostprocType::CellCentered,
        dump,
        VariablesExtracted::All,
        OutputControl::new(true, true, Frequency::Iteration(1)),
    ))];

    let schemes = EquationsEnumConfig {};

    let variables_config = VariablesEnumConfig {};

    let main_config = Config {
        run_name: "Poisson-test_1".to_owned(),
        geometry,
        output: postproc,
    };

    (main_config, schemes, variables_config)
}

#[test]
pub fn main() {
    let directory = file_directory();
    let dump = directory.join("dump");
    clear_dump_folder(&dump);

    examples_default_logger(&directory.join("logs").join("main.log"), LevelFilter::Debug).unwrap();

    let (config, schemes, variables_config) = poisson();
    let mesh = mesh(&config.geometry);
    let mut case = PoissonCase::new(config, mesh, schemes, variables_config, UserFunctions {});

    for _ in 0..10 {
        case.next_step();
        info!("Time: {:?}", case.core.time());
    }

    case.core.plot_2d()
}
