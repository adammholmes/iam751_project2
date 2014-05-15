/* mat.c:
 *
 * Functions that represent a matrix. Also used as the domain
 * for thr ADI PDE solver.
 *
 * Author: Adam M. Holmes
 */

#include "mat.h"


/* Initialize and return an empty matrix of given size */
mat_t *
mat_init(int xsize, int ysize)
{
  mat_t *temp = calloc(1, sizeof(mat_t));
  assert(temp);
  temp->data = calloc(ysize, sizeof(double *));
  assert(temp->data);
  int i; for (i = 0; i < ysize; i++) {
    temp->data[i] = calloc(xsize, sizeof(double));
    assert(temp->data[i]);
  }
  temp->xsize = xsize;
  temp->ysize = ysize;
  return temp;
}


/* Set a value in the matrix */
void
mat_set(int x, int y, double val, mat_t *a)
{
  assert(x >= 0 && x < a->xsize);
  assert(y >= 0 && y < a->ysize);
  a->data[y][x] = val;
}


/* Get a value in the matrix */
double
mat_get(int x, int y, mat_t *a)
{
  assert(x >= 0); assert(x < a->xsize);
  assert(y >= 0); assert(y < a->ysize);
  return a->data[y][x];
}


/* Print a matrix to the given output stream */
void
mat_print(FILE *stream, mat_t *a)
{
  int y; for (y = 0; y < a->ysize; y++) {
    int x; for (x = 0; x < a->xsize; x++) {
      fprintf(stream, "%.8f ", a->data[y][x]);
    }
    fprintf(stream, "\n");
  }
}


/* Print a matrix to the given output stream in a gnuplot friendly way */
void
mat_print_gnuplot(FILE *stream, mat_t *a, double h)
{
  int x; for (x = 0; x < a->xsize; x++) {
    int y; for (y = 0; y < a->ysize; y++) {
      fprintf(stream, "%.8f ", x*h);
      fprintf(stream, "%.8f ", y*h);
      fprintf(stream, "%.8f ", a->data[y][x]);
      fprintf(stream, "\n");
    }
    fprintf(stream, "\n");
  }
}


/* Free matrix from memory */
void
mat_free(mat_t *a)
{
  int i; for (i = 0; i < a->xsize; i++) {
    free(a->data[i]);
  }
  free(a);
}
