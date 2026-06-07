/** \file
 *
 * \brief Eigenvalue functions
 *
 */

#ifndef EIGEN_FUNCS
#define EIGEN_FUNCS

/// Eigenvalue Decomposition Routine
void eigen( double** A, int n, double* d, double** V );

/// Computes all eigenvalues and eigenvalues of a real symmetric matrix
void jacobi( double** A, int n, double* d, double** V );

#endif
