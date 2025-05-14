use std::collections::HashMap;

use cfd_rs_utils::mesh::computational_mesh::*;
use nalgebra::Vector2;
use crate::finite_volume::{base::Field, case::{Case, CaseSystems, VariableFields}, config::CaseConfig};
use super::super::equation::{Variable, System};

#[derive(Clone, Debug, PartialEq)]
pub struct SimpleCase {
    name: String,
    step: usize,
    time: f64,
    time_step: f64,
    
    density: f64,
    kinematic_viscosity: f64,
    
    mesh: Computational2DMesh,
    
    fields: VariableFields,
    systems: CaseSystems,
}

impl Case for SimpleCase {
    fn name(&self) -> &str {
        &self.name
    }
    
    fn step(&self) -> usize {
        self.step
    }
    
    fn time(&self) -> f64 {
        self.time
    }
    
    fn time_step(&self) -> f64 {
        self.time_step
    }
    
    fn export(&self, directory: &str) -> std::io::Result<()> {
        todo!()
    }
    
    fn import_from_file(file_name: &str) -> std::io::Result<()> {
        todo!()
    }
    
    fn fields_list(&self) -> Vec<&Variable> {
        self.fields.map.keys().collect()
    }
    
    #[inline]
    fn field(&self, var: &Variable) -> Option<&Field> {
        self.fields.map.get(var)
    }
    
    #[inline]
    fn field_mut(&mut self, var: &Variable) -> Option<&mut Field> {
        self.fields.map.get_mut(var)
    }
    
    fn equations_list(&self) -> Vec<&String> {
        self.systems.map.keys().collect()
    }
    
    fn equation(&self, name: &str) -> Option<&System> {
        self.systems.map.get(name)
    }
    
    fn equation_mut(&mut self, name: &str) -> Option<&mut System> {
        self.systems.map.get_mut(name)
    }
    
    fn next_step(&mut self) {
        todo!()
    }
    
}

pub fn setup(config: CaseConfig) -> SimpleCase {
    
    
    
    
    todo!()
}