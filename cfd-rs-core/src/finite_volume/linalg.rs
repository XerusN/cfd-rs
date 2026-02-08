use nalgebra::DVector;
use nalgebra_sparse::CsrMatrix;

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

pub fn easy_jacobi(matrix: &CsrMatrix<f64>, rhs: &DVector<f64>, x: &mut DVector<f64>) {
    let mut x_2 = x.clone();
    let mut r = x.clone();
    let mut norm = 1.;
    let mut old_norm = 1.;
    let mut i = 0;

    while norm > 1e-4 {
        // Update
        for (i, phi) in x.iter_mut().enumerate() {
            *phi = rhs[i];
            let row = matrix.get_row(i).unwrap();
            for j in row.col_indices() {
                if *j != i {
                    let a = row.get_entry(*j).unwrap();
                    let a_ij = match a {
                        nalgebra_sparse::SparseEntry::Zero => {
                            panic!("wtf");
                        }
                        nalgebra_sparse::SparseEntry::NonZero(value) => value,
                    };
                    *phi -= a_ij * x_2[*j];
                }
            }
            let a = row.get_entry(i).unwrap();
            let a_ii = match a {
                nalgebra_sparse::SparseEntry::Zero => {
                    panic!("ill conditioned system");
                }
                nalgebra_sparse::SparseEntry::NonZero(value) => value,
            };
            *phi /= a_ii;
        }

        // Residual
        for (i, phi) in r.iter_mut().enumerate() {
            *phi = rhs[i];
            let row = matrix.get_row(i).unwrap();
            for j in row.col_indices() {
                let a = row.get_entry(*j).unwrap();
                let a_ij = match a {
                    nalgebra_sparse::SparseEntry::Zero => {
                        panic!("wtf");
                    }
                    nalgebra_sparse::SparseEntry::NonZero(value) => value,
                };
                *phi -= a_ij * x[*j];
            }
        }

        x_2 = x.clone();

        norm = r.norm();
        if old_norm == norm {
            panic!("Did not converge")
        }
        old_norm = norm;

        i += 1;

        println!("Iter: {i} Norm: {norm}");
        // if i > 100000 {
        //     break;
        // }
    }
}
