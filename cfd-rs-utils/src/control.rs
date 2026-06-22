
#[derive(Clone, Debug, Default, PartialEq)]
pub enum Frequency {
    TimeStep(f64),
    Iteration(usize),
    #[default]
    None,
}

#[derive(Clone, Debug, PartialEq)]
pub struct OutputControl {
    initial_output: bool,
    final_output: bool,
    frequency: Frequency,
}

impl Default for OutputControl {
    fn default() -> Self {
        OutputControl { frequency: Frequency::default(), initial_output: true, final_output: false }
    }
}

impl OutputControl {
    pub fn new(initial_output: bool, final_output: bool, frequency: Frequency) -> Self {
        OutputControl { frequency, initial_output, final_output }
    }
    
    pub fn initial_output(&self) -> bool {
        self.initial_output
    }
    
    pub fn final_output(&self) -> bool {
        self.final_output
    }
    
    pub fn frequency(&self) -> &Frequency {
        &self.frequency
    }
}
