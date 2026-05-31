use cfd_rs::finite_volume::{
    boundary::{BoundaryCondition, BoundaryValue, FieldsBoundaryConditions},
    config::{CaseConfig, GeometryConfig, InitFunc, MeshingConfig, OutputConfig, Schemes},
    equation::{
        discretizations::{
            convection::ConvectionScheme, divergence::DivergenceScheme, laplacian::LaplacianScheme,
            time_schemes::TimeIntegration,
        },
        variables::{ControlVolumeType, Dimension, Variable},
    },
    gradients::{GradientConfig, GradientInterpConfig, GradientScheme},
    mesh::mesh,
    solvers::{diffusion::DiffusionCase, Case},
};
use cfd_rs_utils::control::OutputControl;
use hashbrown::HashMap;
use nalgebra::Point2;

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

    // let geometry = GeometryConfig {
    //     import_path: Some("./target/exports/mesh.cfd".to_string()),
    //     element_size: 0.01,
    // };
    let geometry = GeometryConfig {
        import_path: Some("../meshes/mesh_unstructured.cfd".to_string()),
        meshing: MeshingConfig::AdvancingFront { element_size: 0.01 },
    };

    let output = OutputConfig {
        control: OutputControl::Iteration(1),
        directory: "./exports".to_string(),
    };

    let t = Variable::new("T".to_string(), Dimension::Scalar, ControlVolumeType::Nodes);

    let mut bc_fields = HashMap::new();
    let bc = vec![
        // bot | cart: left
        BoundaryCondition::Dirichlet(BoundaryValue::Scalar(0.)),
        // right | cart: bot
        BoundaryCondition::Neumann(BoundaryValue::Scalar(0.)),
        // top | cart: right
        BoundaryCondition::Dirichlet(BoundaryValue::Scalar(0.)),
        // left | cart: top
        BoundaryCondition::Neumann(BoundaryValue::Scalar(0.)),
    ];
    // let bc = vec![
    //     // bot | cart: left
    //     BoundaryCondition::Dirichlet(BoundaryValue::Scalar(0.)),
    //     // right | cart: bot
    //     BoundaryCondition::Dirichlet(BoundaryValue::Scalar(0.)),
    //     // top | cart: right
    //     BoundaryCondition::Dirichlet(BoundaryValue::Scalar(0.)),
    //     // left | cart: top
    //     BoundaryCondition::Dirichlet(BoundaryValue::Scalar(0.)),
    // ];
    bc_fields.insert(t.clone(), bc);
    let bc_fields = FieldsBoundaryConditions::new(bc_fields);

    let mut initial_fields = HashMap::new();
    initial_fields.insert(t.clone(), InitFunc::Scalar(gaussian));

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

pub fn custom(point: &Point2<f64>) -> f64 {
    point.x + point.y
}

pub fn gaussian(point: &Point2<f64>) -> f64 {
    (-(point.x - 0.33).powi(2) / (2. * (0.1 as f64).powi(2))).exp()
}

fn main() {
    let config = poisson();
    let mesh = mesh(&config.geometry);
    let mut case = DiffusionCase::new(config, mesh);

    case.export_cell_centered().unwrap();

    for _ in 0..20 {
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
        case.export_cell_centered().unwrap();
        println!("{:?}", case.time());
    }
}
