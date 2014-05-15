/* adi.c:
 *
 * Solves a 2D heat equation using the Alternating Direction Implicit
 * method (and the Peaceman-Rachford Scheme).
 *
 * Also include a serial and parallel (MPI) TDMA solver.
 *
 * Author: Adam M. Holmes
 */

#include "adi.h"


/* Perform a timestep using the ADI method on the given domain */
void
adi(mat_t *u, double h, double dt)
{
  // Spacial step to time step ratio
  double bmux = dt/(h*h);
  double bmuy = dt/(h*h);
  
  // Intermediate values used in Peaceman-Rachford difference equations
  mat_t *w = mat_init(u->xsize, u->ysize);
  
  // Y direction ---------------------------------------------------------------
  int m; for (m = 1; m < u->xsize-1; m++) {
    
    // a,b,c are tridiagonal bands, d is result vector
    double *a = calloc(u->xsize, sizeof(double));
    double *b = calloc(u->xsize, sizeof(double));
    double *c = calloc(u->xsize, sizeof(double));
    double *d = calloc(u->xsize, sizeof(double));
    
    // Set up Ax=b problem for every m line
    int el; for (el = 0; el < u->xsize; el++) {
      if (el == 0) {
        // Set A matrix band
        b[el] = 1;
        // Boundary condition at w(0,m)
        d[el] = 0;
      } else if (el == u->xsize-1) {
        // Set A matrix band
        b[el] = 1;
        // Boundary condition at w(L,m)
        d[el] = 0;
      } else {
        // Set A matrix bands
        a[el] = -bmux/2;
        b[el] = 1 + bmux;
        c[el] = -bmux/2;
        // Finite difference
        d[el] = (bmuy/2)*mat_get(el, m-1, u) + (1-bmuy)*mat_get(el, m, u) +
                (bmuy/2)*mat_get(el, m+1, u);
      }
    }
    
    // Solve this m
    TDMA_solve_serial(a, b, c, d, u->xsize);
    
    // Save result to w
    for (el = 1; el < u->xsize-1; el++) {
      mat_set(el, m, d[el], w);
    }
    
    free(a); free(b); free(c); free(d);
  }

  // X direction ---------------------------------------------------------------
  int el; for (el = 1; el < u->xsize-1; el++) {
    
    // a,b,c are tridiagonal bands, d is result vector
    double *a = calloc(u->xsize, sizeof(double));
    double *b = calloc(u->xsize, sizeof(double));
    double *c = calloc(u->xsize, sizeof(double));
    double *d = calloc(u->xsize, sizeof(double));
    
    // Set up Ax=b problem for every el line
    int m; for (m = 0; m < u->xsize; m++) {
      if (m == 0) {
        // Set A matrix band
        b[m] = 1;
        // Boundary condition at u(l,0)
        d[m] = 0;
      } else if (m == u->xsize-1) {
        // Set A matrix band
        b[m] = 1;
        // Boundary condition at u(l,M)
        d[m] = 0;
      } else {
        // Set A matrix bands
        a[m] = -bmuy/2;
        b[m] = 1 + bmuy;
        c[m] = -bmuy/2;
        // Finite difference
        d[m] = (bmux/2)*mat_get(el-1, m, w) + (1-bmux)*mat_get(el, m, w) +
               (bmux/2)*mat_get(el+1, m, w);
      }
    }
    
    // Solve this el
    TDMA_solve_serial(a, b, c, d, u->xsize);
    
    // Save result to w
    for (m = 1; m < u->xsize-1; m++) {
      mat_set(el, m, d[m], u);
    }
    
    free(a); free(b); free(c); free(d);
  }
  
  // Enforce boundary conditions
  int i; for (i = 0; i < u->xsize; i++) {
    mat_set(0, i, 0, u);
    mat_set(u->xsize-1, i, 0, u);
  }
  
  mat_free(w);
}


/* Solve a tridiagonal matrix in serial; put result in d */
void
TDMA_solve_serial(double *a, double *b, double *c, double *d, int n) 
{
  // http://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
  n--;
  c[0] /= b[0];
  d[0] /= b[0];
  int i; for (i = 1; i < n; i++) {
      c[i] /= b[i] - a[i]*c[i-1];
      d[i] = (d[i] - a[i]*d[i-1]) / (b[i] - a[i]*c[i-1]);
  }
  d[n] = (d[n] - a[n]*d[n-1]) / (b[n] - a[n]*c[n-1]);
  for (i = n; i-- > 0;) {
      d[i] -= c[i] * d[i+1];
  }
}



/* Solve a tridiagonal matrix in parallel */
double *
TDMA_solve_parallel(double *a, double *b, double *c, double *d, int len)
{
  double *x = calloc(len, sizeof(double));
  
  
  return x; 
}


