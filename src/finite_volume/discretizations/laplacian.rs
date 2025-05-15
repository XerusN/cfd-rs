use crate::finite_volume::case::GradRequirements;

#[derive(Clone, Debug, PartialEq)]
pub enum LaplacianScheme {
    Centered,
}

impl LaplacianScheme {
    pub fn required_grads(&self) -> GradRequirements {
        match *self {
            Self::Centered => GradRequirements::new(false, false),
        }
    }
}
