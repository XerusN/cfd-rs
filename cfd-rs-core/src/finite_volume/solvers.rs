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

use crate::finite_volume::{
    boundary::{BoundaryCondition},
    config::{GradientConfig, InitFunctionsTrait, SchemesConfig},
    equation::{Component, operations::Op},
    gradients::GradientMethods,
};

use super::{
    config::Schemes,
    equation::{
        variables::{ControlVolumeType, Dimension, Variable},
        Equation, EquationSolver,
    },
    error::CfdError,
    fields::{Field, ScalarField},
};

use std::fmt::Debug;
use std::hash::Hash;

use plotters::prelude::*;

// pub mod burger;
// pub mod convection_diffusion_test;
// pub mod convection_test;
// pub mod diffusion_test;
// pub mod divergence_test;
// pub mod divergence_test2;
// pub mod diffusion;
pub mod poisson;
// pub mod simple;

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
    pub fn new<M: MeshCore, I: InitFunctionsTrait>(
        equations: &EquationsSet,
        mesh: &Mesh<M>,
        variables_hashmaps: VariablesHashMaps<I>,
    ) -> Self {
        let (mut gradient_methods, mut initial_fields, mut boundary_conditions) =
            variables_hashmaps.deconstruct_mut();
        // println!("{:?}", &initial_fields);
        let mut fields = HashMap::new();
        for (var, grad_req) in equations.variables_requirements() {
            // print!("Field init for {}: ", var.name());
            match *var.dim() {
                Dimension::Scalar => {
                    let values = RefCell::new(Field::Scalar(
                        ScalarField::new(
                            mesh,
                            &grad_req,
                            &var,
                            Component::X,
                            &initial_fields.remove(&var).expect(&format!(
                                "An initial field function should be defined for var {var:?}"
                            )),
                            gradient_methods.remove(&var).expect(&format!(
                                "An initial field function should be defined for var {var:?}"
                            )),
                        ),
                        boundary_conditions.remove(&var).expect(&format!(
                            "A boundary condition should be defined for var {var:?}"
                        )),
                    ));
                    fields.insert(var, (values, grad_req));
                }
                Dimension::Vector2 => {
                    let values = RefCell::new(Field::Vector2(
                        Vector2::new(
                            ScalarField::new(
                                mesh,
                                &grad_req,
                                &var,
                                Component::X,
                                initial_fields.get(&var).expect(&format!(
                                    "An initial field function should be defined for var {var:?}"
                                )),
                                gradient_methods.get(&var).expect(&format!(
                                    "A gradient config should be defined for var {var:?}"
                                )).clone(),
                            ),
                            ScalarField::new(
                                mesh,
                                &grad_req,
                                &var,
                                Component::Y,
                                initial_fields.get(&var).expect(&format!(
                                    "An initial field function should be defined for var {var:?}"
                                )),
                                gradient_methods.remove(&var).expect(&format!(
                                    "A gradient config should be defined for var {var:?}"
                                )),
                            ),
                        ),
                        boundary_conditions.remove(&var).expect(&format!(
                            "A boundary condition should be defined for var {var:?}"
                        )),
                    ));
                    fields.insert(var, (values, grad_req));
                }
            }
        }
        Self { map: fields }
    }
}

pub fn find_var_in_fields<'a>(var: &'a Variable, fields: &'a VariableFields) -> &'a RefCell<Field> {
    &fields
        .map
        .get(var)
        .expect(&format!("Missing variable {var:?} in fields",))
        .0
}

#[derive(Clone, Debug, PartialEq)]
pub struct Parameters {
    pub map: HashMap<String, f64>,
}

impl Parameters {
    pub fn new(parameters: Vec<(String, f64)>) -> Self {
        let mut map = HashMap::new();
        for (name, value) in parameters {
            map.insert(name, value);
        }
        Self { map }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct EquationsSet {
    pub map: HashMap<String, Equation>,
}

impl EquationsSet {
    pub fn new() -> Self {
        let map = HashMap::new();
        Self { map }
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

    pub fn add_eq_from_enum<E: EquationsEnum, EC: EquationsEnumConfigTrait<E>>(
        &mut self,
        enum_variant: &E,
        equations_config: &EC,
    ) -> Result<(), CfdError> {
        let name = enum_variant.name();
        let equation = Equation::new(
            enum_variant.lhs(),
            enum_variant.rhs(),
            enum_variant.schemes(equations_config.schemes_configs(enum_variant)),
        )?;
        if let Err(_) = self.map.try_insert(name.to_owned(), equation) {
            Err(CfdError::EquationAlreadyAdded {
                name: name.to_owned(),
            })
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
pub struct EquationSolversSet {
    node: EquationSolver,
    cell: EquationSolver,
}

impl EquationSolversSet {
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

#[derive(Clone, Debug, PartialEq)]
pub struct SolverCore<M: MeshCore> {
    equation_solvers: EquationSolversSet,
    variable_fields: VariableFields,
    equations: EquationsSet,
    mesh: Mesh<M>,
    parameters_set: Parameters,

    name: String,
    step: usize,
    time: f64,
    time_step: f64,
}

impl<M: MeshCore> SolverCore<M> {
    pub fn name(&self) -> &str {
        &self.name
    }

    pub fn step(&self) -> usize {
        self.step
    }

    pub fn time(&self) -> f64 {
        self.time
    }

    pub fn time_step(&self) -> f64 {
        self.time_step
    }
    pub fn set_time_step(&mut self, time_step: f64) {
        self.time_step = time_step
    }

    pub fn increment(&mut self) {
        self.time += self.time_step;
        self.step += 1;
    }

    pub fn import_from_file(_file_name: &str) -> std::io::Result<()> {
        todo!()
    }

    pub fn fields_list(&self) -> Vec<&Variable> {
        self.variable_fields.map.keys().collect()
    }

    pub fn field(&self, var: &Variable) -> Option<Ref<Field>> {
        match self.variable_fields.map.get(var) {
            None => None,
            Some((field, _)) => Some(field.borrow()),
        }
    }

    pub fn field_mut(&mut self, var: &Variable) -> Option<RefMut<Field>> {
        match self.variable_fields.map.get_mut(var) {
            None => None,
            Some((field, _)) => Some(field.borrow_mut()),
        }
    }

    pub fn equations_list(&self) -> Vec<&String> {
        self.equations.map.keys().collect()
    }

    pub fn equation(&self, name: &str) -> Option<&Equation> {
        match self.equations.map.get(name) {
            None => None,
            Some(value) => Some(&value),
        }
    }

    pub fn equation_mut(&mut self, name: &str) -> Option<&mut Equation> {
        match self.equations.map.get_mut(name) {
            None => None,
            Some(value) => Some(value),
        }
    }

    pub fn solver(&self, cvt: &ControlVolumeType) -> &EquationSolver {
        &self.equation_solvers.get_from_cvt(cvt)
    }

    pub fn solver_mut(&mut self, cvt: &ControlVolumeType) -> &mut EquationSolver {
        self.equation_solvers.get_from_cvt_mut(cvt)
    }

    pub fn mesh(&self) -> &Mesh<M> {
        &self.mesh
    }

    pub fn equation_solver_borrow(
        &mut self,
    ) -> (
        &mut EquationSolversSet,
        &mut VariableFields,
        &EquationsSet,
        &Mesh<M>,
    ) {
        (
            &mut self.equation_solvers,
            &mut self.variable_fields,
            &self.equations,
            &self.mesh,
        )
    }

    // pub fn init_core(config: CaseConfig, mesh: Mesh<M>) -> Self {
    //     todo!()
    // }

    pub fn plot_2d(&self) {
        for var in self.fields_list() {
            plot_var(self, var, &(1., 2.));
        }
    }

    /// https://docs.vtk.org/en/latest/vtk_file_formats/vtkxml_file_format.html#unstructuredgrid
    pub fn export_cell_centered(&self, output_dir: &PathBuf) -> io::Result<()> {
        let path = output_dir.join(PathBuf::from(format!(
            "{}_cells_nodes_{:06}.vtu",
            &self.name(),
            &self.step()
        )));
        println!("{:?}", &path);
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
    pub fn export_node_centered(&self, output_dir: &str) -> io::Result<()> {
        let path = PathBuf::from(format!(
            "{}/{}_nodes_{:06}.vtu",
            output_dir,
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
}

fn plot_var<M: MeshCore>(solver_core: &SolverCore<M>, var: &Variable, limits: &(f64, f64)) {
    let path = "../figures/poisson_".to_string() + var.name() + ".jpeg";
    let root = BitMapBackend::new(&path, (1920, 1080)).into_drawing_area();

    root.fill(&WHITE).unwrap();

    let mut chart = ChartBuilder::on(&root)
        .caption("Test", ("sans-serif", 80))
        .margin(5)
        .x_label_area_size(40)
        .y_label_area_size(40)
        .build_cartesian_2d(0f64..1f64, 0f64..1f64)
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

    let field = solver_core.field(var).unwrap();
    let values = match field.deref() {
        Field::Scalar(scalar_field, _) => scalar_field.values(),
        Field::Vector2(vector_field, _) => vector_field.x.values(),
    };

    if let ControlVolumeType::Nodes = var.cvt() {
        chart
            .draw_series(values.iter().zip(0..).map(|(v, i)| {
                Polygon::new::<Vec<(f64, f64)>, RGBColor>(
                    solver_core.mesh().nodes.cv_nodes()[i]
                        .iter()
                        .map(|p| (p.x, p.y))
                        .collect::<Vec<_>>(),
                    ViridisRGB::get_color_normalized(*v, limits.0, limits.1).into(),
                )
            }))
            .unwrap();
    } else {
        chart
            .draw_series(values.iter().zip(0..).map(|(v, i)| {
                Polygon::new::<Vec<(f64, f64)>, RGBColor>(
                    solver_core.mesh().cells.neighboring_nodes()[i]
                        .iter()
                        .map(|node| {
                            (
                                solver_core.mesh().nodes.centers()[*node].x,
                                solver_core.mesh().nodes.centers()[*node].y,
                            )
                        })
                        .collect::<Vec<_>>(),
                    ViridisRGB::get_color_normalized(*v, limits.0, limits.1).into(),
                )
            }))
            .unwrap();
    }
    root.present().unwrap();
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
) -> io::Result<()> {
    for var in solver_core.fields_list() {
        if var.cvt() != cvt {
            continue;
        }
        let temp = solver_core
            .field(var)
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

// pub trait VariablesEnum: Eq + Hash + Clone + Debug {
//     fn name(&self) -> &'static str;
//     fn dim(&self) -> Dimension;
//     fn cvt(&self) -> ControlVolumeType;
//     fn var(&self) -> Variable;
// }

pub trait VariablesEnum: Clone + Debug + Eq + PartialEq + Hash {
    fn name(&self) -> &'static str;
    fn dim(&self) -> Dimension;
    fn cvt(&self) -> ControlVolumeType;
    fn var(&self) -> Variable;
    fn default_gradient_method(&self) -> GradientMethods;
    fn gradient_method(&self, gradient_config: GradientConfig) -> GradientMethods {
        let mut method = self.default_gradient_method();
        gradient_config
            .interp
            .is_some()
            .then(|| method.interp = gradient_config.interp.unwrap());
        gradient_config
            .scheme
            .is_some()
            .then(|| method.scheme = gradient_config.scheme.unwrap());
        method
    }
}

pub trait VariablesEnumConfigTrait<V: VariablesEnum, I: InitFunctionsTrait> {
    fn gradient_configs(&self, variable_enum: &V) -> GradientConfig;
    fn initial_fields(&self, variable_enum: &V) -> I;
    fn boundary_conditions(&self, variable_enum: &V) -> Vec<BoundaryCondition>;
}

pub trait EquationsEnumConfigTrait<E: EquationsEnum> {
    fn schemes_configs(&self, equations_enum: &E) -> SchemesConfig;
}

pub struct VariablesHashMaps<I: InitFunctionsTrait> {
    pub gradient_methods: HashMap<Variable, GradientMethods>,
    pub initial_fields: HashMap<Variable, I>,
    pub boundary_conditions: HashMap<Variable, Vec<BoundaryCondition>>,
}

impl<I: InitFunctionsTrait> VariablesHashMaps<I> {
    pub fn new() -> Self {
        Self {
            gradient_methods: HashMap::new(),
            initial_fields: HashMap::new(),
            boundary_conditions: HashMap::new(),
        }
    }

    pub fn add_var_from_enum<V: VariablesEnum, VC: VariablesEnumConfigTrait<V, I>>(
        &mut self,
        enum_variant: &V,
        var_enum_config: &VC,
    ) {
        self.gradient_methods.insert(
            enum_variant.var(),
            enum_variant.gradient_method(var_enum_config.gradient_configs(enum_variant)),
        );
        self.initial_fields.insert(
            enum_variant.var(),
            var_enum_config.initial_fields(enum_variant),
        );
        self.boundary_conditions.insert(
            enum_variant.var(),
            var_enum_config.boundary_conditions(enum_variant),
        );
    }

    pub fn deconstruct_mut(
        self,
    ) -> (
        HashMap<Variable, GradientMethods>,
        HashMap<Variable, I>,
        HashMap<Variable, Vec<BoundaryCondition>>,
    ) {
        (
            self.gradient_methods,
            self.initial_fields,
            self.boundary_conditions,
        )
    }
}

// Macro to generate getter functions for solver specific variables
#[macro_export]
macro_rules! generate_solver_variables {
    (
        $T:ty,
        $(
            $variant:ident => {
                name: $name:expr,
                dim: $dim:expr,
                cvt: $cvt:expr,
                default_gradient_config: $default_gradient_config:expr,
            }
        ),+ $(,)?
    ) => {
        impl VariablesEnum for $T {
            fn name(&self) -> &'static str {
                match self {
                    $(
                        Self::$variant => $name,
                    )+
                }
            }

            fn dim(&self) -> Dimension {
                match self {
                    $(
                        Self::$variant => $dim,
                    )+
                }
            }

            fn cvt(&self) -> ControlVolumeType {
                match self {
                    $(
                        Self::$variant => $cvt,
                    )+
                }
            }

            fn var(&self) -> Variable {
                match self {
                    $(
                        Self::$variant => Variable::new($name.to_owned(), $dim, $cvt),
                    )+
                }
            }

            fn default_gradient_method(&self) -> GradientMethods {
                match self {
                    $(
                        Self::$variant => $default_gradient_config,
                    )+
                }
            }
        }
    };
}

pub trait EquationsEnum: Clone + Debug + Eq + PartialEq + Hash {
    fn name(&self) -> &'static str;
    fn lhs(&self) -> Op;
    fn rhs(&self) -> Op;
    fn default_schemes(&self) -> Schemes;
    fn schemes(&self, schemes_config: SchemesConfig) -> Schemes {
        let mut schemes = self.default_schemes();
        schemes_config
            .transient
            .is_some()
            .then(|| schemes.transient = schemes_config.transient.unwrap());
        schemes_config
            .convection
            .is_some()
            .then(|| schemes.convection = schemes_config.convection.unwrap());
        schemes_config
            .laplacian
            .is_some()
            .then(|| schemes.laplacian = schemes_config.laplacian.unwrap());
        schemes_config
            .divergence
            .is_some()
            .then(|| schemes.divergence = schemes_config.divergence.unwrap());
        schemes
    }
}

// Macro to generate getter functions for solver specific variables
#[macro_export]
macro_rules! generate_solver_equations {
    (
        $T:ty,
        $(
            $variant:ident => {
                name: $name:expr,
                lhs: $lhs:expr,
                rhs: $rhs:expr,
                default_schemes: $default_schemes:expr,
            }
        ),+ $(,)?
    ) => {
        impl EquationsEnum for $T {
            fn name(&self) -> &'static str {
                match self {
                    $(
                        Self::$variant => $name,
                    )+
                }
            }

            fn lhs(&self) -> Op {
                match self {
                    $(
                        Self::$variant => $lhs,
                    )+
                }
            }

            fn rhs(&self) -> Op {
                match self {
                    $(
                        Self::$variant => $rhs,
                    )+
                }
            }

            fn default_schemes(&self) -> Schemes {
                match self {
                    $(
                        Self::$variant => $default_schemes,
                    )+
                }
            }


        }
    };
}
