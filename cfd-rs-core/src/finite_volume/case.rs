use hashbrown::HashMap;
use nalgebra::Vector2;
use std::{
    cell::{Ref, RefCell, RefMut},
    fs::File,
    io::{self, Write},
    ops::Deref,
    path::PathBuf,
};

use cfd_rs_utils::mesh::assembled_mesh::{Mesh, MeshCore};

use crate::finite_volume::equation::Component;

use super::{
    config::{CaseConfig, Schemes},
    equation::{
        variables::{ControlVolumeType, Dimension, Variable},
        Equation, EquationSolver,
    },
    error::CfdError,
    fields::{Field, ScalarField},
};

// pub mod burger;
// pub mod convection_diffusion_test;
// pub mod convection_test;
// pub mod diffusion_test;
// pub mod divergence_test;
// pub mod divergence_test2;
pub mod poisson;
pub mod simple;

#[derive(Clone, Debug, PartialEq)]
pub struct GradRequirements {
    cell: bool,
    face: bool,
}

impl GradRequirements {
    pub fn new(cell: bool, face: bool) -> Self {
        Self { cell, face }
    }

    pub fn cell(&self) -> bool {
        self.cell
    }

    pub fn face(&self) -> bool {
        self.face
    }

    /// Returns the most restrictive requirement (true)
    pub fn update_requirements(&mut self, other: &Self) {
        self.cell = self.cell | other.cell;
        self.face = self.face | other.face;
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct VariableFields {
    pub map: HashMap<Variable, (RefCell<Field>, GradRequirements)>,
}

impl VariableFields {
    pub fn new<M: MeshCore>(
        equations: &CaseEquations,
        mesh: &Mesh<M>,
        config: &CaseConfig,
    ) -> Self {
        let mut fields = HashMap::new();
        for (var, grad_req) in equations.variables_requirements() {
            // print!("Field init for {}: ", var.name());
            match *var.dim() {
                Dimension::Scalar => {
                    let values = RefCell::new(Field::Scalar(ScalarField::new(
                        mesh,
                        &grad_req,
                        &var,
                        Component::X,
                        config,
                    )));
                    fields.insert(var, (values, grad_req));
                }
                Dimension::Vector2 => {
                    let values = RefCell::new(Field::Vector2(Vector2::new(
                        ScalarField::new(mesh, &grad_req, &var, Component::X, config),
                        ScalarField::new(mesh, &grad_req, &var, Component::Y, config),
                    )));
                    fields.insert(var, (values, grad_req));
                }
            }
        }
        VariableFields { map: fields }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct CaseEquations {
    pub map: HashMap<String, Equation>,
}

impl CaseEquations {
    pub fn new() -> Self {
        let map = HashMap::new();
        CaseEquations { map }
    }

    /// Do not initialize fields before adding all equations here.
    /// Will throw an error if an equation with the same name is already present.
    pub fn add_eq(&mut self, name: String, equation: Equation) -> Result<(), CfdError> {
        if let Err(_) = self.map.try_insert(name.clone(), equation) {
            Err(CfdError::EquationAlreadyAdded { name: name })
        } else {
            Ok(())
        }
    }

    pub fn variables_requirements(&self) -> HashMap<Variable, GradRequirements> {
        let mut variables_glob: HashMap<Variable, GradRequirements> = HashMap::new();

        for eq in self.map.values() {
            for (variable, grad) in eq.fields_required() {
                //println!("{:?}, {:?}", variable, grad);
                if variables_glob.contains_key(variable) {
                    variables_glob
                        .get_mut(variable)
                        .expect("?")
                        .update_requirements(grad);
                } else {
                    variables_glob.insert(variable.clone(), grad.clone());
                }
            }
        }

        variables_glob
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct SolversSet {
    node: EquationSolver,
    cell: EquationSolver,
}

impl SolversSet {
    pub fn new<M: MeshCore>(mesh: &Mesh<M>) -> Self {
        let node = EquationSolver::new(mesh, ControlVolumeType::Nodes);
        let cell = EquationSolver::new(mesh, ControlVolumeType::Cells);
        Self { node, cell }
    }
    
    pub fn get_node(&self) -> &EquationSolver {
        &self.node
    }
    
    pub fn get_node_mut(&mut self) -> &mut EquationSolver {
        &mut self.node
    }
    
    pub fn get_cell(&self) -> &EquationSolver {
        &self.cell
    }
    
    pub fn get_cell_mut(&mut self) -> &mut EquationSolver {
        &mut self.cell
    }
    
    pub fn get_from_cvt(&self, cvt: &ControlVolumeType) -> &EquationSolver {
        match cvt {
            ControlVolumeType::Cells => &self.cell,
            ControlVolumeType::Nodes => &self.node,
        }
    }
    
    pub fn get_from_cvt_mut(&mut self, cvt: &ControlVolumeType) -> &mut EquationSolver {
        match cvt {
            ControlVolumeType::Cells => &mut self.cell,
            ControlVolumeType::Nodes => &mut self.node,
        }
    }
}

pub trait Case<T: MeshCore>: Sized {
    fn name(&self) -> &str;

    fn step(&self) -> usize;

    fn time(&self) -> f64;

    fn time_step(&self) -> f64;

    fn next_step(&mut self);

    fn import_from_file(file_name: &str) -> io::Result<()>;

    /// https://docs.vtk.org/en/latest/vtk_file_formats/vtkxml_file_format.html#unstructuredgrid
    fn export_cell_centered(&self, path: String) -> io::Result<()> {
        let path = PathBuf::from(format!(
            "{}/{}_cells_nodes_{:06}.vtu",
            &self.config().output.directory,
            &self.name(),
            &self.step()
        ));
        let mesh = self.mesh();

        let mut file = File::create(path)?;

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

        export_variables(&mut file, self, &ControlVolumeType::Cells)?;

        writeln!(file, "      </CellData>")?;

        writeln!(file, "      <PointData>")?;

        export_scalar(&mut file, &mesh.nodes.volumes(), "Node volumes")?;

        export_variables(&mut file, self, &ControlVolumeType::Nodes)?;

        writeln!(file, "      </PointData>")?;

        writeln!(file, "    </Piece>")?;
        writeln!(file, "  </UnstructuredGrid>")?;
        writeln!(file, "</VTKFile>")?;

        Ok(())
    }

    /// https://docs.vtk.org/en/latest/vtk_file_formats/vtkxml_file_format.html#unstructuredgrid
    /// https://vtk.org/doc/nightly/html/vtkCellType_8h_source.html
    fn export_node_centered(&self, path: String) -> io::Result<()> {
        let path = PathBuf::from(format!(
            "{}/{}_nodes_{:06}.vtu",
            &self.config().output.directory,
            &self.name(),
            &self.step()
        ));
        let mesh = self.mesh();

        let mut file = File::create(path)?;

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

        export_variables(&mut file, self, &ControlVolumeType::Nodes)?;

        writeln!(file, "      </CellData>")?;

        writeln!(file, "    </Piece>")?;
        writeln!(file, "  </UnstructuredGrid>")?;
        writeln!(file, "</VTKFile>")?;

        Ok(())
    }

    fn config(&self) -> &CaseConfig;

    fn fields_list(&self) -> Vec<&Variable>;

    fn field(&self, var: &Variable) -> Option<Ref<Field>>;

    fn field_mut(&mut self, var: &Variable) -> Option<RefMut<Field>>;

    fn equations_list(&self) -> Vec<&String>;

    fn equation(&self, name: &str) -> Option<&Equation>;

    fn equation_mut(&mut self, name: &str) -> Option<&mut Equation>;

    fn solver(&self, cvt: &ControlVolumeType) -> &EquationSolver;

    fn solver_mut(&mut self, cvt: &ControlVolumeType) -> &mut EquationSolver;

    fn schemes(&self) -> &Schemes;

    fn mesh(&self) -> &Mesh<T>;

    fn equation_solver_borrow(
        &mut self,
    ) -> (
        &mut SolversSet,
        &mut VariableFields,
        &CaseEquations,
        &Mesh<T>,
        &CaseConfig,
    );

    fn new(config: CaseConfig, mesh: Mesh<T>) -> Self;
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

fn export_variables<M: MeshCore, T: Case<M>>(
    file: &mut File,
    case: &T,
    cvt: &ControlVolumeType,
) -> io::Result<()> {
    for var in case.fields_list() {
        if var.cvt() != cvt {
            continue;
        }
        let temp = case
            .field(var)
            .expect("Incoherence between variable list and fields");
        match temp.deref() {
            Field::Scalar(scalar_field) => {
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
            Field::Vector2(fields) => {
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
