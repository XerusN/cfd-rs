use std::io;

pub trait Case {
    fn name(&self) -> &str;
    
    fn step(&self) -> usize;
    
    fn time(&self) -> f64;
    
    fn time_step(&self) -> f64;
    
    fn next_step(&mut self);
    
    fn import_from_file(file_name: &str) -> io::Result<()>;
    
    fn export(&self, directory: &str) -> io::Result<()>;
}