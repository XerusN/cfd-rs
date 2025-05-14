use std::{collections::HashMap, io};

use super::{base::Field, equation::{Variable, System}};

#[derive(Clone, Debug, PartialEq)]
pub struct VariableFields{
    pub map: HashMap<Variable, Field>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct CaseSystems{
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
}