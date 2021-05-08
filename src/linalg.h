#ifndef LINALG_H
#define LINALG_H
#include <vector>

namespace LinalAlg {

double **CreateMatrix(long nrow, long ncol);
void DeleteMatrix(double **m);

/**
 * @brief Tridiag: Sove the tridiagonal system `A x = r`, where `A =
 * Diag(a,b,c)`
 * @param a :  sub-diagonal, vector of n-elements, the firts element `a[0]` is
 * ignored
 * @param b : main diagonal, vector of n-elements
 * @param c : superdiagonal, vector of n-elements, the last element `c[n-1]` is
 * ignored
 * @param n : dimension of the system
 * @param r : right vector, n-elements vector
 * @param x : resulting vector, `x` and `u` can be the same vector (in-place)
 */
void Tridiag(const double *a, const double *b, const double *c, const int n,
             const double *r, double *x);

/**
 * @brief Tridiag2: Same as `Tridiag()`, but solves two systems in the same
 * procedure `A x1 = r1` and `A x2 = r2`.
 */
void Tridiag2(const double *a, const double *b, const double *c, const int n,
              const double *r1, const double *r2, double *x1, double *x2);

/**
 * @brief Tridiag3: Same as `Tridiag()` but solves tree systems in the same
 * procedure `A x1 = r1`, `A x2 = r2` and  `A x3 = r3`
 */
void Tridiag3(const double *a, const double *b, const double *c, const int n,
              const double *r1, const double *r2, const double *r3, double *x1,
              double *x2, double *x3);

/**
 * @brief Cyclic: Solve the cyclic system `A x = r`, where `A = Diag(a,b,c)`
 * except for `A[0,n-1] = a[0]` and `A[n-1,0] = c[n-1]`. This function modifies
 * the vector `b`!!!
 */
void CyclicPack(const double *a, double *b, const double *c, int n,
                const double *r, double *x);

void CyclicPackV1(const double *a, const double *b, const double *c, int n,
                  const double *r, double *x);

};  // namespace LinalAlg

#endif  // LINALG_H
