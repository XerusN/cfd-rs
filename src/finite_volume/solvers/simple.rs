use super::super::equation::{System, Variable};
use crate::finite_volume::{
    base::Field,
    case::{Case, CaseSystems, VariableFields},
    config::{CaseConfig, GeometryConfig, Schemes},
};
use cfd_rs_utils::mesh::computational_mesh::*;
use nalgebra::Vector2;

#[derive(Clone, Debug, PartialEq)]
pub struct SimpleCase {
    name: String,
    step: usize,
    time: f64,
    time_step: f64,

    config: CaseConfig,

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
        match self.fields.map.get(var) {
            None => None,
            Some((field, _)) => Some(field),
        }
    }

    #[inline]
    fn field_mut(&mut self, var: &Variable) -> Option<&mut Field> {
        match self.fields.map.get_mut(var) {
            None => None,
            Some((field, _)) => Some(field),
        }
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

    fn schemes(&self) -> &Schemes {
        &self.config.schemes
    }

    fn mesh(&self) -> &Computational2DMesh {
        &self.mesh
    }
    
    fn equation_solver_borrow(&mut self) -> (&mut CaseSystems, &mut VariableFields, &Computational2DMesh, &CaseConfig) {
        (&mut self.systems, &mut self.fields, &self.mesh, &self.config)
    }
    
    fn next_step(&mut self) {
        //
        //
        //
        //
        todo!()
    }
}

pub fn setup(config: CaseConfig) -> SimpleCase {
    todo!()
}
