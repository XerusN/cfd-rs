#[macro_export]
macro_rules! laplacian {
    ( $variable:expr, $integration:expr ) => {{
        let var: Variable = $variable;
        let int: IntegrationCategory = $integration;
        Op::FieldOperator(FieldOperator::DifferentialOperator(
            DifferentialOperator::Laplacian(var, int),
        ))
    }};
}

#[macro_export]
macro_rules! divergence {
    ( $variable:expr, $integration:expr ) => {{
        let var: Variable = $variable;
        let int: IntegrationCategory = $integration;
        Op::FieldOperator(FieldOperator::DifferentialOperator(
            DifferentialOperator::Divergence(var, int),
        ))
    }};
}

#[macro_export]
macro_rules! gradient {
    ( $variable:expr ) => {{
        let var: Variable = $variable;
        Op::FieldOperator(FieldOperator::Gradient(var))
    }};
}

#[macro_export]
macro_rules! convection {
    ( $variable:expr, $speed:expr, $integration:expr ) => {{
        let var: Variable = $variable;
        let int: IntegrationCategory = $integration;
        let speed: Variable = $speed;
        Op::FieldOperator(FieldOperator::DifferentialOperator(
            DifferentialOperator::Convection {
                var: var,
                speed: speed,
                integration: int,
            },
        ))
    }};
}

#[macro_export]
macro_rules! time_derivative {
    ( $variable:expr ) => {{
        let var: Variable = $variable;
        Op::FieldOperator(FieldOperator::DifferentialOperator(
            DifferentialOperator::TimeDerivative(var),
        ))
    }};
}
