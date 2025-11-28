use cfd_rs::finite_volume::{boundary::{BoundaryCondition, BoundaryValue, FieldsBoundaryConditions}, case::{Case, poisson::PoissonCase}, config::{CaseConfig, GeometryConfig, InitFunc, MeshingConfig, OutputConfig, Schemes}, equation::{discretizations::{convection::ConvectionScheme, divergence::DivergenceScheme, laplacian::LaplacianScheme, time_schemes::TimeIntegration}, variables::{ControlVolumeType, Dimension, Variable}}, gradients::{GradientConfig, GradientInterpConfig, GradientScheme}, mesh::mesh};
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
        import_path: Some("../meshes/mesh3.cfd".to_string()),
        meshing: MeshingConfig::AdvancingFront { element_size: 0.01 },
    };

    let output = OutputConfig {
        control: OutputControl::Iteration(1),
        directory: "./target/exports".to_string(),
    };
    
    let t = Variable::new("T".to_string(), Dimension::Scalar, ControlVolumeType::Cells);
    
    let mut bc_fields = HashMap::new();
    let bc = vec![
        BoundaryCondition::Dirichlet(BoundaryValue::Scalar(0.)),
        BoundaryCondition::Dirichlet(BoundaryValue::Scalar(0.)),
        BoundaryCondition::Dirichlet(BoundaryValue::Scalar(200.)),
        BoundaryCondition::Dirichlet(BoundaryValue::Scalar(300.)),
    ];
    bc_fields.insert(t.clone(), bc);
    let bc_fields = FieldsBoundaryConditions::new(bc_fields);

    let mut initial_fields = HashMap::new();
    initial_fields.insert(
        t.clone(),
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
    let config = poisson();
    let mesh = mesh(&config.geometry);
    let mut case = PoissonCase::new(config, mesh);

    case.export_cell_centered("".to_owned()).unwrap();

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
        case.export_cell_centered("".to_owned()).unwrap();
        println!("{:?}", case.time());
    }
}
