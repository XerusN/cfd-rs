

#[macro_export]
macro_rules! laplacian {
    ( $variable:expr, $integration:expr ) => {
        {
            let var: &Variable = $variable;
            let int: IntegrationCategory = $integration;
            Op::FieldOperator(FieldOperator::DifferentialOperator(
            DifferentialOperator::Laplacian(var.clone(), int)))
        }
    };
}