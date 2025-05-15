use crate::finite_volume::case::GradRequirements;

#[derive(Clone, Debug, PartialEq)]
pub enum DivergenceScheme {
    Centered,
}

impl DivergenceScheme {
    pub fn required_grads(&self) -> GradRequirements {
        match *self {
            Self::Centered => GradRequirements::new(false, false),
        }
    }
}
