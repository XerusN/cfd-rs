use std::path::PathBuf;

use cfd_rs_utils::control::OutputControl;

use crate::finite_volume::equation::variables::Variable;

#[derive(Debug, Clone, PartialEq)]
pub enum VariablesExtracted {
    None,
    Some(Vec<Variable>),
    All
}


#[derive(Debug, Clone, PartialEq)]
pub struct GridPostprocConfig {
    directory: PathBuf,
    export_count: usize,
    variables: VariablesExtracted,
    output_control: OutputControl,
}

impl GridPostprocConfig {
    
    pub fn new(directory: PathBuf, variables: VariablesExtracted, output_control: OutputControl) -> Self {
        GridPostprocConfig { directory, export_count: 0, variables, output_control }
    }
    
    
}