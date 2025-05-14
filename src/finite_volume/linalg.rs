#[derive(Debug, Clone, PartialEq)]
pub enum LinaglConfig {
    Jacobi(ConvergenceConfig),
    ConjugateGadient(ConvergenceConfig),
}

/// The f64 is the tolerance
#[derive(Debug, Clone, PartialEq)]
pub enum ConvergenceConfig {
    EuclideanNorm(f64),
    MaxNorm(f64),
}
