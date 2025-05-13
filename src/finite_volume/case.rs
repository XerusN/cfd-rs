use std::io;

use super::{base::Field, equation::Variable};

pub trait Case {
    fn name(&self) -> &str;
    
    fn step(&self) -> usize;
    
    fn time(&self) -> f64;
    
    fn time_step(&self) -> f64;
    
    fn next_step(&mut self);
    
    fn import_from_file(file_name: &str) -> io::Result<()>;
    
    fn export(&self, directory: &str) -> io::Result<()>;
    
    fn variable_field(&self, var: Variable) -> &Vec<Field>;
    
    fn variable_field_mut(&mut self, var: Variable) -> &mut Vec<Field>;
}