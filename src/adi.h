/* adi.h:
 *
 * Solves a 2D heat equation using the Alternating Direction Implicit
 * method (and the Peaceman-Rachford Scheme).
 *
 * Author: Adam M. Holmes
 */

#ifndef ADI_H
#define ADI_H

#include "mat.h"


/* Perform a timestep using the ADI method on the given domain */
void adi(mat_t *u, double h, double dt);

/* Solve a tridiagonal matrix in serial; put result in d */
void TDMA_solve_serial(double *a, double *b, double *c, double *d, int len);

/* Solve a tridiagonal matrix in parallel */
double * TDMA_solve_parallel(double *a, double *b, double *c, double *d, int len);


#endif