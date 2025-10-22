use nalgebra::Vector2;
use std::ops::{Add, Div, Mul, Sub};

use crate::finite_volume::equation::FieldOperator;

#[derive(Clone, Debug, PartialEq)]
pub enum Op {
    Add(Box<(Op, Op)>),
    Sub(Box<(Op, Op)>),
    FieldOperator(FieldOperator),
    MulScalar(f64, Box<Op>),
    MulVector(Vector2<f64>, Box<Op>),
    DivScalar(f64, Box<Op>),
    Scalar(f64),
    Vector2(Vector2<f64>),
}

impl Op {
    pub fn collect_field_operators(&self, collector: &mut Vec<FieldOperator>) {
        match self {
            Op::Add(pair) | Op::Sub(pair) => {
                Self::collect_field_operators(&pair.0, collector);
                Self::collect_field_operators(&pair.1, collector);
            }
            Op::MulScalar(_, inner) | Op::DivScalar(_, inner) => {
                Self::collect_field_operators(inner, collector);
            }
            Op::MulVector(_, inner) => {
                Self::collect_field_operators(inner, collector);
            }
            Op::FieldOperator(f_op) => collector.push(f_op.clone()),
            Op::Scalar(_) | Op::Vector2(_) => {}
        }
    }
}

impl Add for Op {
    type Output = Op;

    fn add(self, rhs: Self) -> Self::Output {
        Op::Add(Box::new((self, rhs)))
    }
}

impl Sub for Op {
    type Output = Op;

    fn sub(self, rhs: Self) -> Self::Output {
        Op::Sub(Box::new((self, rhs)))
    }
}

impl Mul<f64> for Op {
    type Output = Op;

    fn mul(self, rhs: f64) -> Self::Output {
        Op::MulScalar(rhs, Box::new(self))
    }
}

impl Mul<Op> for f64 {
    type Output = Op;

    fn mul(self, rhs: Op) -> Self::Output {
        Op::MulScalar(self, Box::new(rhs))
    }
}

impl Mul<Vector2<f64>> for Op {
    type Output = Op;

    fn mul(self, rhs: Vector2<f64>) -> Self::Output {
        Op::MulVector(rhs, Box::new(self))
    }
}

impl Mul<Op> for Vector2<f64> {
    type Output = Op;

    fn mul(self, rhs: Op) -> Self::Output {
        Op::MulVector(self, Box::new(rhs))
    }
}

impl Div<f64> for Op {
    type Output = Op;

    fn div(self, rhs: f64) -> Self::Output {
        Op::DivScalar(rhs, Box::new(self))
    }
}
