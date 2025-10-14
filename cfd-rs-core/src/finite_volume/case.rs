use hashbrown::HashMap;
use nalgebra::Vector2;
use std::{
    cell::{Ref, RefCell, RefMut},
    fs::File,
    io::{self, Write},
    ops::Deref,
    path::PathBuf,
};

use cfd_rs_utils::mesh::computational_mesh::Computational2DMesh;

use crate::finite_volume::equation::Component;

use super::{
    base::{CellScalarField, Field},
    config::{self, CaseConfig, Schemes},
    equation::{Dimension, Equation, EquationSolver, Variable},
    error::CfdError,
};

pub mod burger;
pub mod convection_diffusion_test;
pub mod convection_test;
pub mod diffusion_test;
pub mod divergence_test;
pub mod divergence_test2;
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
    pub fn new(equations: &CaseEquations, mesh: &Computational2DMesh, config: &CaseConfig) -> Self {
        let mut fields = HashMap::new();
        for (var, grad_req) in equations.variables_requirements() {
            // print!("Field init for {}: ", var.name());
            match *var.dim() {
                Dimension::Scalar => {
                    let values = RefCell::new(Field::Scalar(CellScalarField::new(
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
                        CellScalarField::new(mesh, &grad_req, &var, Component::X, config),
                        CellScalarField::new(mesh, &grad_req, &var, Component::Y, config),
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

pub trait Case {
    fn name(&self) -> &str;

    fn step(&self) -> usize;

    fn time(&self) -> f64;

    fn time_step(&self) -> f64;

    fn next_step(&mut self);

    fn import_from_file(file_name: &str) -> io::Result<()>;

    /// https://docs.vtk.org/en/latest/vtk_file_formats/vtkxml_file_format.html#unstructuredgrid
    fn export(&self) -> io::Result<()> {
        let path = PathBuf::from(format!(
            "{}/{}_{:06}.vtu",
            &self.config().output.directory,
            &self.name(),
            &self.step()
        ));

        let mut file = File::create(&path)?;

        writeln!(
            file,
            "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">"
        )?;
        writeln!(file, "  <UnstructuredGrid>")?;
        writeln!(
            file,
            "    <Piece NumberOfPoints=\"{}\" NumberOfCells=\"{}\">",
            self.mesh().num_vertices(),
            self.mesh().num_cells()
        )?;
        writeln!(file, "      <Points>")?;
        writeln!(
            file,
            "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">"
        )?;
        write!(file, "          ")?;
        for vertex in self.mesh().vertices() {
            write!(file, "{} {} 0 ", vertex.x, vertex.y)?;
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
        for cell in self.mesh().cells() {
            for vertex_id in cell.vertices_id() {
                write!(file, "{} ", vertex_id.0)?;
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
        for cell in self.mesh().cells() {
            offset += cell.vertices_id().len();
            write!(file, "{} ", offset)?;
        }
        writeln!(file)?;
        writeln!(file, "        </DataArray>")?;
        writeln!(
            file,
            "        <DataArray type=\"UInt64\" Name=\"types\" format=\"ascii\">"
        )?;
        write!(file, "          ")?;
        for cell in self.mesh().cells() {
            if cell.vertices_id().len() == 3 {
                write!(file, "5 ")?;
            } else if cell.vertices_id().len() == 4 {
                write!(file, "7 ")?;
            }
        }
        writeln!(file)?;
        writeln!(file, "        </DataArray>")?;
        writeln!(file, "      </Cells>")?;

        writeln!(file, "      <CellData>")?;
        // Does not support vector fields yet
        for var in self.fields_list() {
            let temp = self
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

    fn solver(&self) -> &EquationSolver;

    fn solver_mut(&mut self) -> &mut EquationSolver;

    fn schemes(&self) -> &Schemes;

    fn mesh(&self) -> &Computational2DMesh;

    fn equation_solver_borrow(
        &mut self,
    ) -> (
        &mut EquationSolver,
        &mut VariableFields,
        &CaseEquations,
        &Computational2DMesh,
        &CaseConfig,
    );

    fn new(config: CaseConfig) -> Self;
}
