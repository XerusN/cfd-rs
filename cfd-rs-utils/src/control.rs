#[derive(Clone, Debug, Default, PartialEq)]
pub enum OutputControl {
    TimeStep(f64),
    Iteration(usize),
    #[default]
    Final,
    None,
}
