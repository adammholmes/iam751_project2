/* adi.c:
 *
 * Solves a 2D heat equation using the Alternating Direction Implicit
 * method (and the Peaceman-Rachford Scheme).
 *
 * Author: Adam M. Holmes
 */

#include "adi.h"


/* Perform a timestep using the ADI method on the given domain */
void
adi(mat_t *u, double h, double dt)
{
  // MPI info
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  
  int el; int m;
  int lines = u->xsize;
  int l_per_proc = (lines / (size-1));
  
  assert(size > 1);
  assert(lines > size);
  
  // Intermediate values used in Peaceman-Rachford difference equations
  mat_t *w = mat_init(u->xsize, u->ysize);
  
  
  // Y direction ---------------------------------------------------------------
  
  if (rank == 0) {
    
    // Grab every line from every process
    if (size > 2) {
      int t; for (t = 1; t < size-1; t++) {
        int m; for (m = l_per_proc * (t-1); m < l_per_proc * t; m++) {
          if (t == 1 && m == 0) {m++;}
          double res[lines];
          MPI_Recv(&res, lines, MPI_DOUBLE, t, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
          // Save result to w
          for (el = 1; el < u->xsize-1; el++) {
            mat_set(el, m, res[el], w);
          }
        }
      }
    } else { // Basically a serial solution
      int t; for (t = 1; t < size-1; t++) {
        int m; for (m = l_per_proc * (t-1); m < size-1; m++) {
          if (t == 1 && m == 0) {m++;}
          double res[lines];
          MPI_Recv(&res, lines, MPI_DOUBLE, t, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
          // Save result to w
          for (el = 1; el < u->xsize-1; el++) {
            mat_set(el, m, res[el], w);
          }
        }
      }
    }
    
    // Grab the last group of lines
    if (size > 2) {
      int m; for (m = l_per_proc * (size-2); m < lines-1; m++) {
        double res[lines];
        MPI_Recv(&res, lines, MPI_DOUBLE, size-1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        // Save result to w
        for (el = 1; el < u->xsize-1; el++) {
          mat_set(el, m, res[el], w);
        }
      }
    }
    
  } else if (rank == size-1 && size > 2) {
    
    // Last grouping of lines
    int m; for (m = l_per_proc * (size-2); m < lines-1; m++) {
      double *x = y_direction(u, dt, h, m);
      MPI_Send(x, lines, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
      free(x);
    }
    
  } else if (rank == 1 && size == 2) {
    
    // Basically a serial solution
    int m; for (m = l_per_proc * (rank-1); m < lines-1; m++) {
      if (rank == 1 && m == 0) {m++;}
      double *x = y_direction(u, dt, h, m);
      MPI_Send(x, lines, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
      free(x);
    }
    
  } else {
    
    // Normal line groupings
    int m; for (m = l_per_proc * (rank-1); m < l_per_proc * rank; m++) {
      if (rank == 1 && m == 0) {m++;}
      double *x = y_direction(u, dt, h, m);
      MPI_Send(x, lines, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
      free(x);
    }
    
  }
  
  // Wait for Y passes to complete
  MPI_Barrier(MPI_COMM_WORLD);
  
  // Update every process with completed intermediate values
  int i; for (i = 0; i < w->xsize; i++) {
    MPI_Bcast(w->data[i], w->xsize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  }
  
  
  // X Direction ---------------------------------------------------------------
  
  if (rank == 0) {

    // Grab every line from every process
    if (size > 2) {
      int t; for (t = 1; t < size-1; t++) {
        int el; for (el = l_per_proc * (t-1); el < l_per_proc * t; el++) {
          if (t == 1 && el == 0) {el++;}
          double res[lines];
          MPI_Recv(&res, lines, MPI_DOUBLE, t, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
          // Save result to w
          for (m = 1; m < u->xsize-1; m++) {
            mat_set(el, m, res[m], u);
          }
        }
      }
    } else { // Basically a serial solution
      int t; for (t = 1; t < size-1; t++) {
        int el; for (el = l_per_proc * (t-1); el < lines-1; el++) {
          if (t == 1 && el == 0) {el++;}
          double res[lines];
          MPI_Recv(&res, lines, MPI_DOUBLE, t, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
          // Save result to w
          for (m = 1; m < u->xsize-1; m++) {
            mat_set(el, m, res[m], u);
          }
        }
      }
    }
    
    // Grab the last group of lines
    if (size > 2) {
      int el; for (el = l_per_proc * (size-2); el < lines-1; el++) {
        double res[lines];
        MPI_Recv(&res, lines, MPI_DOUBLE, size-1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        // Save result to w
        for (m = 1; m < u->xsize-1; m++) {
          mat_set(el, m, res[m], u);
        }
      }
    }
    
    // Enforce boundary conditions to complete step
    int i; for (i = 0; i < u->xsize; i++) {
      mat_set(0, i, 0, u);
      mat_set(u->xsize-1, i, 0, u);
    }
    
  } else if (rank == size-1 && size > 2) {
    
    // Last grouping of lines
    int el; for (el = l_per_proc * (size-2); el < lines-1; el++) {
      double *x = x_direction(w, dt, h, el);
      MPI_Send(x, lines, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
      free(x);
    }
    
  } else if (rank == 1 && size == 2) {
    
    // Basically a serial solution
    int el; for (el = l_per_proc * (rank-1); el < lines-1; el++) {
      if (rank == 1 && el == 0) {el++;}
      double *x = x_direction(w, dt, h, el);
      MPI_Send(x, lines, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
      free(x);
    }
    
  } else {
    
    // Normal line groupings
    int el; for (el = l_per_proc * (rank-1); el < l_per_proc * rank; el++) {
      if (rank == 1 && el == 0) {el++;}
      double *x = x_direction(w, dt, h, el);
      MPI_Send(x, lines, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
      free(x);
    }
    
  }
  
  // Wait for X passes to complete
  MPI_Barrier(MPI_COMM_WORLD);
  
  // Update every process with completed final values
  for (i = 0; i < u->xsize; i++) {
    MPI_Bcast(u->data[i], u->xsize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  }
  
  mat_free(w);
}


/* Solve a line in the y direction */
double *
y_direction(mat_t *u, double dt, double h, int m)
{
  // Spacial step to time step ratio
  double bmux = dt/(h*h);
  double bmuy = dt/(h*h);
  
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
  TDMA_solve(a, b, c, d, u->xsize);
  
  free(a); free(b); free(c);
  
  return d;
}


/* Solve a line in the x direction */
double *
x_direction(mat_t *w, double dt, double h, int el)
{
  // Spacial step to time step ratio
  double bmux = dt/(h*h);
  double bmuy = dt/(h*h);
  
  // a,b,c are tridiagonal bands, d is result vector
  double *a = calloc(w->xsize, sizeof(double));
  double *b = calloc(w->xsize, sizeof(double));
  double *c = calloc(w->xsize, sizeof(double));
  double *d = calloc(w->xsize, sizeof(double));
  
  // Set up Ax=b problem for every el line
  int m; for (m = 0; m < w->xsize; m++) {
    if (m == 0) {
      // Set A matrix band
      b[m] = 1;
      // Boundary condition at u(l,0)
      d[m] = 0;
    } else if (m == w->xsize-1) {
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
  TDMA_solve(a, b, c, d, w->xsize);
  
  free(a); free(b); free(c);
  
  return d;
}


/* Solve a tridiagonal matrix; put result in d */
void
TDMA_solve(double *a, double *b, double *c, double *d, int n) 
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
