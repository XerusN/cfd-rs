use std::{collections::HashMap, io};

use cfd_rs_utils::mesh::computational_mesh::Computational2DMesh;

use super::{
    base::Field,
    config::Schemes,
    equation::{System, Variable},
};

#[derive(Clone, Debug, PartialEq)]
pub struct GradRequirements {
    cell: bool,
    face: bool,
}

impl GradRequirements {
    pub fn new(cell: bool, face: bool) -> Self {
        Self{cell, face}
    }
    
    pub fn cell(&self) -> bool {
        self.cell
    }

    pub fn face(&self) -> bool {
        self.face
    }

    /// Returns the most restrictive requirement (true)
    pub fn update_requirements(&mut self, other: Self) {
        self.cell = self.cell | other.cell;
        self.face = self.face | other.face;
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct VariableFields {
    pub map: HashMap<Variable, (Field, GradRequirements)>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct CaseSystems {
    pub map: HashMap<String, System>,
}

pub trait Case {
    fn name(&self) -> &str;

    fn step(&self) -> usize;

    fn time(&self) -> f64;

    fn time_step(&self) -> f64;

    fn next_step(&mut self);

    fn import_from_file(file_name: &str) -> io::Result<()>;

    fn export(&self, directory: &str) -> io::Result<()>;

    fn fields_list(&self) -> Vec<&Variable>;

    fn field(&self, var: &Variable) -> Option<&Field>;

    fn field_mut(&mut self, var: &Variable) -> Option<&mut Field>;

    fn equations_list(&self) -> Vec<&String>;

    fn equation(&self, name: &str) -> Option<&System>;

    fn equation_mut(&mut self, name: &str) -> Option<&mut System>;

    fn schemes(&self) -> &Schemes;

    fn mesh(&self) -> &Computational2DMesh;
}
