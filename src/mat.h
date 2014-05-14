/* mat.h:
 *
 * Functions that represent a matrix. Also used as the domain
 * for thr ADI PDE solver.
 *
 * Author: Adam M. Holmes
 */

#ifndef MAT_H
#define MAT_H

#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>


typedef struct 
{
  double **data;
  int xsize;
  int ysize;
} mat_t;


/* Initialize and return an empty matrix of given size */
mat_t * mat_init(int xsize, int ysize);

/* Set a value in the matrix */
void mat_set(int x, int y, double val, mat_t *a);

/* Get a value in the matrix */
double mat_get(int x, int y, mat_t *a);

/* Print a matrix to the given output steam */
void mat_print(FILE *stream, mat_t *a);

/* Print a matrix to the given output steam in a gnuplot friendly way */
void mat_print_gnuplot(FILE *stream, mat_t *a, double h);

/* Free matrix from memory */
void mat_free(mat_t *a);


#endif
