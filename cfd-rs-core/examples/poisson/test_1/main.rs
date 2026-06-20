use std::path::{Path, PathBuf};

use cfd_rs::finite_volume::{
    boundary::{BoundaryCondition, BoundaryValue},
    config::{GeometryConfig, GradientConfig, InitFunctionsTrait, MeshingConfig, OutputConfig, SchemesConfig},
    solvers::{EquationsEnumConfigTrait, VariablesEnumConfigTrait, poisson::{Config, MainEquations, MainVariables, PoissonUserFunctions, PoissonCase}},
    mesh::mesh,
};
use cfd_rs_utils::{control::OutputControl};
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
    var: MainVariables
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
        InitFunctions { var: variable_enum.clone() }
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
            },
            MainVariables::Laplacian => {
                vec![
                    BoundaryCondition::Neumann(BoundaryValue::Scalar(0.)),
                    BoundaryCondition::Neumann(BoundaryValue::Scalar(0.)),
                    BoundaryCondition::Neumann(BoundaryValue::Scalar(0.)),
                    BoundaryCondition::Neumann(BoundaryValue::Scalar(0.)),
                ]
            },
            MainVariables::Grad => {
                vec![
                    BoundaryCondition::Neumann(BoundaryValue::Vector2(Vector2::new(0., 0.))),
                    BoundaryCondition::Neumann(BoundaryValue::Vector2(Vector2::new(0., 0.))),
                    BoundaryCondition::Neumann(BoundaryValue::Vector2(Vector2::new(0., 0.))),
                    BoundaryCondition::Neumann(BoundaryValue::Vector2(Vector2::new(0., 0.))),
                ]
            },
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

impl PoissonUserFunctions for UserFunctions {
    
}

#[allow(dead_code)]
fn poisson() -> (Config, EquationsEnumConfig, VariablesEnumConfig) {

    let geometry = GeometryConfig {
        import_path: Some("../meshes/mesh4.cfd".to_string()),
        meshing: MeshingConfig::AdvancingFront { element_size: 0.01 },
    };

    let output = OutputConfig {
        control: OutputControl::Iteration(1),
        directory: "./exports".to_string(),
    };
    
    let schemes = EquationsEnumConfig{};

    let variables_config = VariablesEnumConfig{};

    let main_config = Config {
        run_name: "Poisson-test_1".to_owned(),
        geometry,
        output,
    };
    
    (main_config, schemes, variables_config)
}

#[test]
pub fn main() {
    let directory = file_directory();
    let dump = directory.join("dump");
    
    let (config, schemes, variables_config) = poisson();
    let mesh = mesh(&config.geometry);
    let mut case = PoissonCase::new(config, mesh, schemes, variables_config, UserFunctions{});
    
    case.core.export_cell_centered(&dump).unwrap();

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
        case.core.export_cell_centered(&dump).unwrap();
        println!("{:?}", case.core.time());
    }

    case.core.plot_2d()
}
