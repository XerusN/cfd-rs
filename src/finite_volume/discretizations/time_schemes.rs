#[derive(Clone, Debug, PartialEq)]
pub enum TimeIntegration {
    ForwardEuler,
}

impl TimeIntegration {
    pub fn time_integration_type(&self) -> TimeIntegrationCategory {
        match *self {
            Self::ForwardEuler => TimeIntegrationCategory::Explicit,
            _ => panic!("TimeIntegrationCategory not implemented for {self:?}"),
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub enum TimeIntegrationCategory {
    Implicit,
    Explicit,
}
