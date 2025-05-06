

#[derive(Clone, Debug, PartialEq)]
pub enum SpaceDiscretizationConfig {
    CentralDifference,
    Upwind,
    SecondOrderUpwind,
    FROMM,
    Quick,
}

#[derive(Clone, Debug, PartialEq)]
pub enum TimeDiscretizationConfig {
    FirstOrderEuler,
}