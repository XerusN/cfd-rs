use thiserror::Error;

use super::equation::{Op, Variable};

#[derive(Clone, Debug, Default, Error, PartialEq)]
pub enum CfdError {
    #[default]
    #[error(
        "An Unspecified error happened, you can blame the crate developer for the lack of details"
    )]
    Unspecified,
    #[error(
        "Multiple Variables (e.g. {var1:?} and {var2:?}) are recognized as being the unknown for the following equation: {lhs:?} = {rhs:?}"
    )]
    EquationInconsistentUnknown {
        lhs: Op,
        rhs: Op,
        var1: Variable,
        var2: Variable,
    },
    #[error(
        "No Unknown was indentified (you need a least on implicit term or time derivative) for the following equation: {lhs:?} = {rhs:?}"
    )]
    EquationNoUnknown { lhs: Op, rhs: Op },
}
