use super::{
    base::{CellScalarField, Field},
    case::Case,
    equation::*,
    interpolations::GradientInterpConfig,
};
use cfd_rs_utils::mesh::computational_mesh::Computational2DMesh;
use nalgebra::Vector2;

#[derive(Clone, Debug, PartialEq)]
pub enum GradientScheme {
    GreenGaussCompact,
    GreenGaussExtended,
    LeastSquare,
}

#[derive(Clone, Debug, PartialEq)]
pub struct GradientConfig {
    pub scheme: GradientScheme,
    pub interp: GradientInterpConfig,
}

pub fn update_grads(field: &mut Field, mesh: &Computational2DMesh, config: &GradientConfig) {
    match field {
        Field::Scalar(field) => update_grad_scalar(field, mesh, config),
        Field::Vector2(field) => {
            update_grad_scalar(&mut field.x, mesh, config);
            update_grad_scalar(&mut field.y, mesh, config);
        }
    }
}

fn update_grad_scalar(
    field: &mut CellScalarField,
    mesh: &Computational2DMesh,
    config: &GradientConfig,
) {
    match config.scheme {
        GradientScheme::GreenGaussCompact => (),
        _ => panic!("GradientScheme not implemented for {:?}", config.scheme),
    }

    match config.interp {
        _ => panic!("GradientInterp not implemented for {:?}", config.interp),
    }
}
