/* adi.h:
 *
 * Solves a 2D heat equation using the Alternating Direction Implicit
 * method (and the Peaceman-Rachford Scheme).
 *
 * Author: Adam M. Holmes
 */

#ifndef ADI_H
#define ADI_H

#include <mpi.h>
#include <string.h>
#include "mat.h"


typedef struct 
{
  int line_count;
  double **lines;
} mpi_data;


/* Perform a timestep using the ADI method on the given domain */
void adi(mat_t *u, double h, double dt);

/* Solve a line in the y direction */
double * y_direction(mat_t *u, double dt, double h, int m);

/* Solve a line in the x direction */
double * x_direction(mat_t *w, double dt, double h, int el);

/* Solve a tridiagonal matrix; put result in d */
void TDMA_solve(double *a, double *b, double *c, double *d, int len);


#endif
