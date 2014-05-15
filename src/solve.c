/* solve.c:
 *
 * u_t = u_xx + u_yy
 *
 * Solves a 2D heat equation using the Alternating Direction Implicit
 * method (and the Peaceman-Rachford Scheme).
 *
 * Calculated data is stored under '../data/'. This folder needs
 * to exist for the program to work. A data file is created for
 * every timestep.
 *
 * Each line of the ADI scheme is done in parallel using MPI.
 *
 * Author: Adam M. Holmes
 *
 * http://youtu.be/dpAJEny81cU
 */

#include <math.h>
#include "adi.h"

#define H 0.01
#define DT 0.0002
#define SIZE ((1 / H) * 2) + 1
#define PI 3.1415926535897932


/* Solution function */
double
func(double x, double y, double t, double h) {
  return exp(-t) * sin(PI * (x*h)) * sin(PI * (y*h));
}


/* Set initial conditions */
void
init(mat_t *a, double h) {
  int y; for(y = 0; y < a->ysize; y++) {
    int x; for(x = 0; x < a->xsize; x++) {
      double val = func(x, y, 0, h);
      mat_set(x, y, val, a);
    }
  }
}


int
main(int argc, char **argv)
{
  // Option exists to modify the default number of timesteps
  int timesteps = 5;
  if (argc > 1) {
    timesteps = atof(argv[1]);
  }
  
  // Create and initialize domain; create first data file
  mat_t *temp = mat_init(SIZE,SIZE);
  init(temp, H);
  FILE *f = fopen("data/0.dat", "w"); assert(f);
  mat_print_gnuplot(f, temp, H);
  fclose(f);
  
  MPI_Init(NULL, NULL);
  
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  
  // Step through time
  int t; for (t = 1; t < timesteps; t++) {
    // Perform ADI to advance one timestep
    adi(temp, H, DT);
    
    // First process should write results to a file
    if (rank == 0) {
      // Figure out name of this timestep's output file
      char filename[10 + (int)(log10(t)+1)];
      strcpy(filename, "data/");
      sprintf(&filename[5], "%d", t);
      strcat(filename, ".dat");
    
      // Create this timestep's data file; write to it; close file
      FILE *f = fopen(filename, "w"); assert(f);
      mat_print_gnuplot(f, temp, H);
      fclose(f);
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
  }
  
  MPI_Finalize();
  
  return 0;
}
