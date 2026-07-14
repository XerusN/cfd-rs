use std::{
    cell::{Ref, RefCell},
    fs::create_dir,
    io,
    ops::Deref,
    path::{Path, PathBuf},
};

use cfd_rs_utils::{
    control::{Frequency, OutputControl},
    mesh::assembled_mesh::{Mesh, MeshCore},
};
use log::debug;
use nalgebra::DVector;
use plotters::{
    backend::BitMapBackend,
    chart::ChartBuilder,
    coord::Shift,
    drawing::{DrawingArea, IntoDrawingArea},
    element::Polygon,
    style::{colors::colormaps::ViridisRGB, RGBColor, WHITE},
};
use serde::de::value;

use crate::finite_volume::{
    equation::variables::{ControlVolumeType, Dimension, Variable},
    fields::Field,
    solvers::SolverCore,
};

#[derive(Debug, Clone, PartialEq)]
pub enum VariablesExtracted {
    Some(Vec<Variable>),
    All,
}

#[derive(Debug, Clone, PartialEq)]
pub enum PlotPostprocType {
    Plane2D(Plane2DConfig),
    // add None, lines, planes, ...?
}

#[derive(Debug, Clone, PartialEq)]
pub enum Limits {
    None,
    Scalar(),
}

#[derive(Debug, Clone, PartialEq)]
pub struct PlotPostprocConfig {
    name: String,
    postproc_type: PlotPostprocType,
    directory: PathBuf,
    variables: VariablesExtracted,
    output_control: OutputControl,
    export_count: RefCell<Option<usize>>,
    last_export_step: RefCell<Option<usize>>,
    last_export_time: RefCell<Option<f64>>,
    var_limits: Option<()>,
}

#[derive(Clone, PartialEq, Debug, Default)]
pub struct Plane2DConfig {
    pub x_limits: (Option<f64>, Option<f64>),
    pub y_limits: (Option<f64>, Option<f64>),
}

impl PlotPostprocConfig {
    pub fn new(
        name: String,
        postproc_type: PlotPostprocType,
        directory: PathBuf,
        variables: VariablesExtracted,
        output_control: OutputControl,
        var_limits: Option<()>,
    ) -> Self {
        if let Some(_) = var_limits {
            todo!("Implement var limits in plots, find an appropriate structure to be able to declare limits for variables extracted and for vec magnitude/directions");
        }

        Self {
            name,
            postproc_type,
            directory,
            export_count: RefCell::new(None),
            variables,
            output_control,
            last_export_step: RefCell::new(None),
            last_export_time: RefCell::new(None),
            var_limits,
        }
    }

    pub fn name(&self) -> &str {
        &self.name
    }

    pub fn postproc_type(&self) -> &PlotPostprocType {
        &self.postproc_type
    }

    pub fn directory(&self) -> &Path {
        self.directory.as_path()
    }

    pub fn export_count<'a>(&'a self) -> Ref<'a, Option<usize>> {
        self.export_count.borrow()
    }

    pub fn variables(&self) -> &VariablesExtracted {
        &self.variables
    }

    pub fn output_control(&self) -> &OutputControl {
        &self.output_control
    }

    pub fn new_plots<M: MeshCore>(&self, core: &SolverCore<M>) -> io::Result<()> {
        let mut export_count = self.export_count.borrow_mut();
        let count;
        match export_count.as_mut() {
            None => {
                *export_count = Some(0);
                count = 0;
            }
            Some(value) => {
                *value += 1;
                count = *value;
            }
        }
        let mut last_export_time = self.last_export_time.borrow_mut();
        match last_export_time.as_mut() {
            None => {
                *last_export_time = Some(core.time());
            }
            Some(value) => {
                *value = core.time();
            }
        }
        let mut last_export_step = self.last_export_step.borrow_mut();
        match last_export_step.as_mut() {
            None => {
                *last_export_step = Some(core.step());
            }
            Some(value) => {
                *value = core.step();
            }
        }

        

        let output_vars = match self.variables() {
            VariablesExtracted::All => core.fields_list(),
            VariablesExtracted::Some(vars) => {
                let mut temp = vec![];
                for var in core.fields_list() {
                    if vars.contains(var) {
                        temp.push(var);
                    }
                }
                temp
            }
        };

        let directory = self
            .directory()
            .join(PathBuf::from(format!("{}", &self.name())));
        let result = create_dir(&directory);
        match result {
            Ok(_) => (),
            Err(error) => {
                match error.kind() {
                    io::ErrorKind::AlreadyExists => (),
                    _ => return Err(error),
                }
            },
        }
        
        
        debug!("Output of {:?} in {:?}", self.name(), directory);

        match self.postproc_type() {
            PlotPostprocType::Plane2D(config) => {
                for var in output_vars {
                    let path =
                        directory.join(PathBuf::from(format!("{}_{:06}.jpeg", var.name(), count)));
                    plot_var_2d(core, &path, var.name(), var, None, config)
                }
            }
        }

        Ok(())
    }

    pub fn try_output<M: MeshCore>(&self, core: &SolverCore<M>) -> io::Result<()> {
        let mut export = false;

        if (core.step() == 0) & self.output_control().initial_output() {
            export = true;
            debug!("Initial output")
        }

        if core.is_final_iter() & self.output_control().final_output() {
            export = true;
            debug!("Final output")
        }

        match self.output_control().frequency() {
            Frequency::None => (),
            Frequency::Iteration(it) => match *self.last_export_step.borrow() {
                None => export = true,
                Some(last_step) => {
                    if last_step == core.step() {
                        ()
                    } else if (core.step() - last_step) % it == 0 {
                        export = true;
                        debug!("Frequency output")
                    }
                }
            },
            Frequency::TimeStep(dt) => {
                let quotient = core.time().div_euclid(*dt);
                match *self.last_export_time.borrow() {
                    None => export = true,
                    Some(last_time) => {
                        if last_time == core.time() {
                            ()
                        } else if quotient * dt > last_time {
                            export = true;
                            debug!("Frequency output")
                        }
                    }
                }
            }
        }

        if export {
            self.new_plots(core)
        } else {
            Ok(())
        }
    }
}

fn plot_var_2d<M: MeshCore>(
    solver_core: &SolverCore<M>,
    path: &Path,
    title: &str,
    var: &Variable,
    limits: Option<()>,
    config: &Plane2DConfig,
) {
    let root = BitMapBackend::new(path, (1920, 1080)).into_drawing_area();

    root.fill(&WHITE).unwrap();
    let root = root.titled(title, ("sans-serif", 100)).unwrap();

    let roots = match var.dim() {
        Dimension::Scalar => vec![root],
        Dimension::Vector2 => {
            let (left, temp) = root.split_horizontally(640);
            let (middle, right) = temp.split_horizontally(640);
            vec![left, middle, right]
        }
    };

    let bounding_box = solver_core.mesh().bounding_box();
    // root.split_evenly(_);
    let x_limits = (
        match config.x_limits.0 {
            None => bounding_box.0.x,
            Some(lim) => lim,
        },
        match config.x_limits.1 {
            None => bounding_box.1.x,
            Some(lim) => lim,
        },
    );
    let y_limits = (
        match config.y_limits.0 {
            None => bounding_box.0.y,
            Some(lim) => lim,
        },
        match config.y_limits.1 {
            None => bounding_box.1.y,
            Some(lim) => lim,
        },
    );

    // let limits = (
    //     match limits.0 {
    //         None => values.min(),
    //         Some(lim) => lim,
    //     },
    //     match limits.1 {
    //         None => values.max(),
    //         Some(lim) => lim,
    //     }
    // );
    let field = solver_core.field(var).unwrap();
    match field.deref() {
        Field::Scalar(scalar_field, _) => {
            let values = scalar_field.values();
            let limits = (values.min(), values.max());

            plot_field(
                &roots[0],
                solver_core.mesh(),
                values,
                limits,
                var.cvt(),
                x_limits,
                y_limits,
                None,
            );
        }
        Field::Vector2(vector_field, _) => {
            for i in 0..3 {
                let values = match i {
                    0 => vector_field.x.values(),
                    1 => vector_field.y.values(),
                    2 => {
                        &vector_field.compute_mag()
                    }
                    _ => panic!("Should never happen"),
                };
                let title = match i {
                    0 => "X Component",
                    1 => "Y Component",
                    2 => "Magnitude (L2)",
                    _ => panic!("Should never happen"),
                };
                let limits = (values.min(), values.max());
                plot_field(
                    &roots[i],
                    solver_core.mesh(),
                    values,
                    limits,
                    var.cvt(),
                    x_limits,
                    y_limits,
                    Some(title),
                )
            }
        }
    };
}

fn plot_field<M: MeshCore>(
    plot_root: &DrawingArea<BitMapBackend, Shift>,
    mesh: &Mesh<M>,
    values: &DVector<f64>,
    mut limits: (f64, f64),
    cvt: &ControlVolumeType,
    x_limits: (f64, f64),
    y_limits: (f64, f64),
    title: Option<&str>,
) {
    let title = match title {
        None => "",
        Some(s) => s,
    };

    let mut chart = ChartBuilder::on(&plot_root)
        .caption(title, ("sans-serif", 80))
        .margin(5)
        .x_label_area_size(40)
        .y_label_area_size(40)
        .build_cartesian_2d(x_limits.0..x_limits.1, y_limits.0..y_limits.1)
        .unwrap();

    chart
        .configure_mesh()
        .x_labels(5)
        .y_labels(5)
        .max_light_lines(4)
        .x_label_offset(35)
        .y_label_offset(25)
        .label_style(("sans-serif", 20))
        .draw()
        .unwrap();
    
    if limits.0 == limits.1 {
        limits.1 += f64::EPSILON;
    }

    match cvt {
        ControlVolumeType::Nodes => {
            chart
                .draw_series(values.iter().zip(0..).map(|(v, i)| {
                    Polygon::new::<Vec<(f64, f64)>, RGBColor>(
                        mesh.nodes.cv_nodes()[i]
                            .iter()
                            .map(|p| (p.x, p.y))
                            .collect::<Vec<_>>(),
                        ViridisRGB::get_color_normalized(*v, limits.0, limits.1).into(),
                    )
                }))
                .unwrap();
        }
        ControlVolumeType::Cells => {
            chart
                .draw_series(values.iter().zip(0..).map(|(v, i)| {
                    Polygon::new::<Vec<(f64, f64)>, RGBColor>(
                        mesh.cells.neighboring_nodes()[i]
                            .iter()
                            .map(|node| {
                                (mesh.nodes.centers()[*node].x, mesh.nodes.centers()[*node].y)
                            })
                            .collect::<Vec<_>>(),
                        ViridisRGB::get_color_normalized(*v, limits.0, limits.1).into(),
                    )
                }))
                .unwrap();
        }
    }
    plot_root.present().unwrap();
}
