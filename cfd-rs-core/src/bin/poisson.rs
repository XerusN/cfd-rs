use std::ops::Deref;

use cfd_rs::finite_volume::{
    boundary::{BoundaryCondition, BoundaryValue, FieldsBoundaryConditions}, case::{Case, poisson::PoissonCase}, config::{CaseConfig, GeometryConfig, InitFunc, MeshingConfig, OutputConfig, Schemes}, equation::{
        discretizations::{
            convection::ConvectionScheme, divergence::DivergenceScheme, laplacian::LaplacianScheme,
            time_schemes::TimeIntegration,
        },
        variables::{ControlVolumeType, Dimension, Variable},
    }, fields::Field, gradients::{GradientConfig, GradientInterpConfig, GradientScheme}, mesh::mesh
};
use cfd_rs_utils::{control::OutputControl, mesh::assembled_mesh::MeshCore};
use hashbrown::HashMap;
use nalgebra::{Point2, Vector2};
use plotters::prelude::*;

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
    // let geometry = GeometryConfig {
    //     import_path: Some("../meshes/mesh_unstructured.cfd".to_string()),
    //     meshing: MeshingConfig::AdvancingFront { element_size: 0.01 },
    // };
    let geometry = GeometryConfig {
        import_path: Some("../meshes/mesh4.cfd".to_string()),
        meshing: MeshingConfig::AdvancingFront { element_size: 0.01 },
    };

    let output = OutputConfig {
        control: OutputControl::Iteration(1),
        directory: "./exports".to_string(),
    };

    let t = Variable::new("T".to_string(), Dimension::Scalar, ControlVolumeType::Nodes);
    let grad_t = Variable::new(
        "Grad T".to_string(),
        Dimension::Vector2,
        ControlVolumeType::Nodes,
    );
    let lap = Variable::new(
        "Laplacian".to_string(),
        Dimension::Scalar,
        ControlVolumeType::Nodes,
    );

    let mut bc_fields = HashMap::new();
    let bc = vec![
        // bot | cart: left
        BoundaryCondition::Dirichlet(BoundaryValue::Scalar(1.)),
        // right | cart: bot
        BoundaryCondition::Neumann(BoundaryValue::Scalar(0.)),
        // top | cart: right
        BoundaryCondition::Dirichlet(BoundaryValue::Scalar(2.)),
        // left | cart: top
        BoundaryCondition::Neumann(BoundaryValue::Scalar(0.)),
    ];
    bc_fields.insert(t.clone(), bc);
    let bc = vec![
        BoundaryCondition::Neumann(BoundaryValue::Vector2(Vector2::new(0., 0.))),
        BoundaryCondition::Neumann(BoundaryValue::Vector2(Vector2::new(0., 0.))),
        BoundaryCondition::Neumann(BoundaryValue::Vector2(Vector2::new(0., 0.))),
        BoundaryCondition::Neumann(BoundaryValue::Vector2(Vector2::new(0., 0.))),
    ];
    bc_fields.insert(grad_t.clone(), bc);
    let bc = vec![
        BoundaryCondition::Neumann(BoundaryValue::Scalar(0.)),
        BoundaryCondition::Neumann(BoundaryValue::Scalar(0.)),
        BoundaryCondition::Neumann(BoundaryValue::Scalar(0.)),
        BoundaryCondition::Neumann(BoundaryValue::Scalar(0.)),
    ];
    bc_fields.insert(lap.clone(), bc);
    let bc_fields = FieldsBoundaryConditions::new(bc_fields);

    let mut initial_fields = HashMap::new();
    initial_fields.insert(t.clone(), InitFunc::Scalar(custom));
    initial_fields.insert(lap.clone(), InitFunc::Scalar(constant));
    initial_fields.insert(grad_t.clone(), InitFunc::Vector2(constant, constant));

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

fn main() {
    let config = poisson();
    let mesh = mesh(&config.geometry);
    let mut case = PoissonCase::new(config, mesh);

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
        case.export_cell_centered().unwrap();
        println!("{:?}", case.time());
    }
    
    plot(&case)
}


fn plot<M: MeshCore>(case: &PoissonCase<M>) {
    for var in case.fields_list() {
        plot_var(case, var, &(1., 2.));
    }
}

fn plot_var<M: MeshCore>(case: &PoissonCase<M>, var: &Variable, limits: &(f64, f64)) {
    let path = "../figures/poisson_".to_string() + var.name() + ".jpeg";
    let root = BitMapBackend::new(&path, (1920, 1080)).into_drawing_area();
    
    root.fill(&WHITE).unwrap();
    
    let mut chart = ChartBuilder::on(&root)
        .caption("Test", ("sans-serif", 80))
        .margin(5)
        .x_label_area_size(40)
        .y_label_area_size(40)
        .build_cartesian_2d(0f64..1f64, 0f64..1f64).unwrap();
    
    chart
        .configure_mesh()
        .x_labels(5)
        .y_labels(5)
        .max_light_lines(4)
        .x_label_offset(35)
        .y_label_offset(25)
        .label_style(("sans-serif", 20))
        .draw().unwrap();
    
    let field = case.field(var).unwrap();
    let values = match field.deref() {
        Field::Scalar(scalar_field) => scalar_field.values(),
        Field::Vector2(vector_field) => vector_field.x.values(),
    };
    
    if let ControlVolumeType::Nodes = var.cvt() {
        chart.draw_series(
            values.iter().zip(0..).map(|(v, i)| {
                Polygon::new::<Vec<(f64, f64)>, RGBColor>(
                    case.mesh().nodes.cv_nodes()[i].iter().map(|p| (p.x, p.y)).collect::<Vec<_>>(),
                    ViridisRGB::get_color_normalized(*v, limits.0, limits.1).into()
                )
            })
        ).unwrap();
    } else {
        chart.draw_series(
            values.iter().zip(0..).map(|(v, i)| {
                Polygon::new::<Vec<(f64, f64)>, RGBColor>(
                    case.mesh().cells.neighboring_nodes()[i].iter().map(|node| (case.mesh().nodes.centers()[*node].x, case.mesh().nodes.centers()[*node].y)).collect::<Vec<_>>(),
                    ViridisRGB::get_color_normalized(*v, limits.0, limits.1).into()
                )
            })
        ).unwrap();
    }
    root.present().unwrap();
}