use std::{cell::{Ref, RefCell}, fs::File, io::{self, Write}, ops::Deref, path::{Path, PathBuf}};

use cfd_rs_utils::{control::{Frequency, OutputControl}, mesh::assembled_mesh::{Mesh, MeshCore}};
use log::warn;
use nalgebra::Vector2;

use crate::{finite_volume::{equation::variables::{ControlVolumeType, Variable}, fields::Field, solvers::SolverCore}, post_processing::OutputConfig};

#[derive(Debug, Clone, PartialEq)]
pub enum VariablesExtracted {
    None,
    Some(Vec<Variable>),
    All
}


#[derive(Debug, Clone, PartialEq)]
pub enum GridPostprocType {
    CellCentered,
    NodeCentered,
    // add None, lines, planes, ...?
}


#[derive(Debug, Clone, PartialEq)]
pub struct GridPostprocConfig {
    name: String,
    postproc_type: GridPostprocType,
    directory: PathBuf,
    variables: VariablesExtracted,
    output_control: OutputControl,
    export_count: RefCell<Option<usize>>,
    last_export_step: RefCell<Option<usize>>,
    last_export_time: RefCell<Option<f64>>,
}

impl GridPostprocConfig {
    
    pub fn new(name: String, postproc_type: GridPostprocType, directory: PathBuf, variables: VariablesExtracted, output_control: OutputControl) -> Self {
        GridPostprocConfig { name, postproc_type, directory, export_count: RefCell::new(None), variables, output_control, last_export_step: RefCell::new(None), last_export_time: RefCell::new(None) }
    }
    
    pub fn name(&self) -> &str {
        &self.name
    }
    
    pub fn postproc_type(&self) -> &GridPostprocType {
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
    
    pub fn new_file<M: MeshCore>(&self, core: &SolverCore<M>) -> io::Result<()> {
        let mut export_count = self.export_count.borrow_mut();
        let count;
        match *export_count {
            None => {
                *export_count = Some(0);
                count = 0;
            },
            Some(mut value) => {
                value += 1;
                count = value;
            }
        }
        let mut last_export_time = self.last_export_time.borrow_mut();
        match *last_export_time {
            None => {
                *last_export_time = Some(core.time());
            },
            Some(mut value) => {
                value = core.time();
            }
        }
        let mut last_export_step = self.last_export_step.borrow_mut();
        match *last_export_step {
            None => {
                *last_export_step = Some(core.step());
            },
            Some(mut value) => {
                value = core.step();
            }
        }
        
        let path = self.directory().join(PathBuf::from(format!(
            "{}_{:06}.vtu",
            &self.name(),
            count
        )));
        
        match self.postproc_type() {
            GridPostprocType::CellCentered => export_cell_centered(core, &path, self.variables())?,
            GridPostprocType::NodeCentered => export_node_centered(core, &path, self.variables())?,
        }
        
        Ok(())
        
    }

    pub fn try_output<M: MeshCore>(&self, core: &SolverCore<M>) -> io::Result<()> {
        
        let mut export = false;
        
        if (core.step() == 0) & self.output_control().initial_output() {
            export = true;
        }
        
        if core.is_final_iter() & self.output_control().final_output() {
            export = true;
        }
        
        match self.output_control().frequency() {
            Frequency::None => (),
            Frequency::Iteration(it) => {
                if core.step() % it == 0 {
                    export = true;
                }
            }
            Frequency::TimeStep(dt) => {
                let quotient = core.time().div_euclid(*dt);
                match *self.last_export_time.borrow() {
                    None => export = true,
                    Some(last_time) => {
                        if quotient*dt > last_time {
                            export = true;
                        }
                    }
                }
            }
        }
        
        if export {
            self.new_file(core)
        } else {
            Ok(())
        }
    }
}

/// https://docs.vtk.org/en/latest/vtk_file_formats/vtkxml_file_format.html#unstructuredgrid
fn export_cell_centered<M: MeshCore>(core: &SolverCore<M>, output_path: &Path, vars: &VariablesExtracted) -> io::Result<()> {
    
    println!("{:?}", &output_path);
    
    let mesh = core.mesh();
    
    let (vars_node, vars_cells) = match vars {
        VariablesExtracted::None => (vec!(), vec!()),
        VariablesExtracted::All => {
            let mut temp_node = vec!();
            let mut temp_cell = vec!();
            for var in core.fields_list() {
                match var.cvt() {
                    ControlVolumeType::Cells => temp_cell.push(var),
                    ControlVolumeType::Nodes => temp_node.push(var),
                }
            }
            (temp_cell, temp_node)
        },
        VariablesExtracted::Some(values) => {
            let mut temp_node = vec!();
            let mut temp_cell = vec!();
            for var in core.fields_list() {
                if values.contains(var) {
                    match var.cvt() {
                        ControlVolumeType::Cells => temp_cell.push(var),
                        ControlVolumeType::Nodes => temp_node.push(var),
                    }
                }
            }
            (temp_cell, temp_node)
        }
    };

    let mut file = File::create(output_path)?;

    writeln!(
        file,
        "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">"
    )?;

    export_mesh(&mut file, mesh, &ControlVolumeType::Cells)?;

    writeln!(file, "      <CellData>")?;

    export_scalar(&mut file, &mesh.cells.volumes(), "Cell volumes")?;

    writeln!(
        file,
        "        <DataArray type=\"Float64\" Name=\"{}\" format=\"ascii\">",
        "Cell id",
    )?;
    write!(file, "          ")?;
    for id in 0..mesh.cells.n {
        write!(file, "{} ", id)?;
    }
    writeln!(file)?;
    writeln!(file, "        </DataArray>")?;

    export_variables(&mut file, core, &ControlVolumeType::Cells, &vars_cells)?;

    writeln!(file, "      </CellData>")?;

    writeln!(file, "      <PointData>")?;

    export_scalar(&mut file, &mesh.nodes.volumes(), "Node volumes")?;

    export_variables(&mut file, core, &ControlVolumeType::Nodes, &vars_node)?;

    writeln!(file, "      </PointData>")?;

    writeln!(file, "    </Piece>")?;
    writeln!(file, "  </UnstructuredGrid>")?;
    writeln!(file, "</VTKFile>")?;

    Ok(())
}

/// https://docs.vtk.org/en/latest/vtk_file_formats/vtkxml_file_format.html#unstructuredgrid
/// https://vtk.org/doc/nightly/html/vtkCellType_8h_source.html
fn export_node_centered<M: MeshCore>(core: &SolverCore<M>, output_path: &Path, vars: &VariablesExtracted) -> io::Result<()> {

    let mesh = core.mesh();
    
    let vars_node = match vars {
        VariablesExtracted::None => vec!(),
        VariablesExtracted::All => {
            let mut temp_node = vec!();
            for var in core.fields_list() {
                match var.cvt() {
                    ControlVolumeType::Cells => (),
                    ControlVolumeType::Nodes => temp_node.push(var),
                }
            }
            temp_node
        },
        VariablesExtracted::Some(values) => {
            let mut temp_node = vec!();
            for var in core.fields_list() {
                if values.contains(var) {
                    match var.cvt() {
                        ControlVolumeType::Cells => panic!("Cannot export Cells variables in mode NodeCentered"),
                        ControlVolumeType::Nodes => temp_node.push(var),
                    }
                }
            }
            temp_node
        }
    };

    let mut file = File::create(output_path)?;

    writeln!(
        file,
        "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">"
    )?;
    
    writeln!(file, "      <CellData>")?;

    export_scalar(&mut file, &mesh.nodes.volumes(), "Node volumes")?;

    writeln!(
        file,
        "        <DataArray type=\"Float64\" Name=\"{}\" format=\"ascii\">",
        "Node id",
    )?;
    write!(file, "          ")?;
    for id in 0..mesh.nodes.n {
        write!(file, "{} ", id)?;
    }
    writeln!(file)?;
    writeln!(file, "        </DataArray>")?;
    
    export_variables(&mut file, core, &ControlVolumeType::Nodes, &vars_node)?;

    writeln!(file, "      </CellData>")?;

    writeln!(file, "    </Piece>")?;
    writeln!(file, "  </UnstructuredGrid>")?;
    writeln!(file, "</VTKFile>")?;

    Ok(())
}

fn export_mesh<M: MeshCore>(
    file: &mut File,
    mesh: &Mesh<M>,
    cvt: &ControlVolumeType,
) -> io::Result<()> {
    match *cvt {
        ControlVolumeType::Cells => {
            writeln!(file, "  <UnstructuredGrid>")?;
            writeln!(
                file,
                "    <Piece NumberOfPoints=\"{}\" NumberOfCells=\"{}\">",
                mesh.nodes.n, mesh.cells.n
            )?;
            writeln!(file, "      <Points>")?;
            writeln!(
                file,
                "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">"
            )?;
            write!(file, "          ")?;
            for node in mesh.nodes.centers() {
                write!(file, "{} {} 0 ", node.x, node.y)?;
            }
            writeln!(file)?;
            writeln!(file, "        </DataArray>")?;
            writeln!(file, "      </Points>")?;

            writeln!(file, "      <Cells>")?;
            writeln!(
                file,
                "        <DataArray type=\"UInt64\" Name=\"connectivity\" format=\"ascii\">"
            )?;
            write!(file, "          ")?;
            for nodes in mesh.cells.neighboring_nodes() {
                for i_node in nodes {
                    write!(file, "{} ", i_node)?;
                }
            }
            writeln!(file)?;
            writeln!(file, "        </DataArray>")?;
            writeln!(
                file,
                "        <DataArray type=\"UInt64\" Name=\"offsets\" format=\"ascii\">"
            )?;
            write!(file, "          ")?;
            let mut offset = 0;
            for faces in mesh.cells.neighboring_pairs() {
                offset += faces.len();
                write!(file, "{} ", offset)?;
            }
            writeln!(file)?;
            writeln!(file, "        </DataArray>")?;
            writeln!(
                file,
                "        <DataArray type=\"UInt64\" Name=\"types\" format=\"ascii\">"
            )?;
            write!(file, "          ")?;
            for faces in mesh.cells.neighboring_pairs() {
                if faces.len() == 3 {
                    // Triangle
                    write!(file, "5 ")?;
                } else if faces.len() >= 4 {
                    // Polygon
                    write!(file, "7 ")?;
                }
            }
            writeln!(file)?;
            writeln!(file, "        </DataArray>")?;
            writeln!(file, "      </Cells>")?;
        }
        ControlVolumeType::Nodes => {
            writeln!(file, "  <UnstructuredGrid>")?;
            let mut cv_points_number = 0;
            for cv_nodes in mesh.nodes.cv_nodes() {
                cv_points_number += cv_nodes.len();
            }
            writeln!(
                file,
                "    <Piece NumberOfPoints=\"{}\" NumberOfCells=\"{}\">",
                cv_points_number, mesh.nodes.n
            )?;
            writeln!(file, "      <Points>")?;
            writeln!(
                file,
                "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">"
            )?;
            write!(file, "          ")?;
            for cv_nodes in mesh.nodes.cv_nodes() {
                for node in cv_nodes {
                    write!(file, "{} {} 0 ", node.x, node.y)?;
                }
            }
            writeln!(file)?;
            writeln!(file, "        </DataArray>")?;
            writeln!(file, "      </Points>")?;

            writeln!(file, "      <Cells>")?;
            writeln!(
                file,
                "        <DataArray type=\"UInt64\" Name=\"connectivity\" format=\"ascii\">"
            )?;
            write!(file, "          ")?;
            let mut current_cv_node = 0;
            for cv_nodes in mesh.nodes.cv_nodes() {
                for _ in cv_nodes {
                    write!(file, "{} ", current_cv_node)?;
                    current_cv_node += 1;
                }
            }
            writeln!(file)?;
            writeln!(file, "        </DataArray>")?;
            writeln!(
                file,
                "        <DataArray type=\"UInt64\" Name=\"offsets\" format=\"ascii\">"
            )?;
            write!(file, "          ")?;
            let mut offset = 0;
            for cv_nodes in mesh.nodes.cv_nodes() {
                offset += cv_nodes.len();
                write!(file, "{} ", offset)?;
            }
            writeln!(file)?;
            writeln!(file, "        </DataArray>")?;
            writeln!(
                file,
                "        <DataArray type=\"UInt64\" Name=\"types\" format=\"ascii\">"
            )?;
            write!(file, "          ")?;
            for cv_nodes in mesh.nodes.cv_nodes() {
                if cv_nodes.len() == 3 {
                    write!(file, "5 ")?;
                } else if cv_nodes.len() >= 4 {
                    write!(file, "7 ")?;
                }
            }
            writeln!(file)?;
            writeln!(file, "        </DataArray>")?;
            writeln!(file, "      </Cells>")?;
        }
    }

    Ok(())
}

fn export_scalar(file: &mut File, data: &[f64], name: &str) -> io::Result<()> {
    writeln!(
        file,
        "        <DataArray type=\"Float64\" Name=\"{}\" format=\"ascii\">",
        name,
    )?;
    write!(file, "          ")?;
    for value in data {
        write!(file, "{} ", value)?;
    }
    writeln!(file)?;
    writeln!(file, "        </DataArray>")?;
    Ok(())
}

fn export_vector(file: &mut File, data: &[Vector2<f64>], name: &str) -> io::Result<()> {
    writeln!(
        file,
        "        <DataArray type=\"Float64\" Name=\"{}\" format=\"ascii\">",
        name,
    )?;
    write!(file, "          ")?;
    for value in data {
        write!(file, "{} ", value)?;
    }
    writeln!(file)?;
    writeln!(file, "        </DataArray>")?;
    Ok(())
}

fn export_variables<M: MeshCore>(
    file: &mut File,
    solver_core: &SolverCore<M>,
    cvt: &ControlVolumeType,
    vars: &[&Variable],
) -> io::Result<()> {
    for var in vars {
        if var.cvt() != cvt {
            panic!("Should not have wrong cvt variable here, checked before")
        }
        let temp = solver_core
            .field(*var)
            .expect("Incoherence between variable list and fields");
        match temp.deref() {
            Field::Scalar(scalar_field, _) => {
                writeln!(
                    file,
                    "        <DataArray type=\"Float64\" Name=\"{}\" format=\"ascii\">",
                    var.name(),
                )?;
                write!(file, "          ")?;
                let field = scalar_field.values();
                for value in field {
                    write!(file, "{} ", value)?;
                }
            }
            Field::Vector2(fields, _) => {
                writeln!(
                    file,
                    "        <DataArray NumberOfComponents=\"2\" type=\"Float64\" Name=\"{}\" format=\"ascii\">",
                    var.name(),
                )?;
                write!(file, "          ")?;
                let field_x = fields.x.values();
                let field_y = fields.y.values();
                for i in 0..field_x.len() {
                    write!(file, "{} {} ", field_x[i], field_y[i])?;
                }
            }
        };

        writeln!(file)?;
        writeln!(file, "        </DataArray>")?;
    }
    Ok(())
}