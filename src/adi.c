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
  
  double halfbmux = dt/(2*h*h);
  double halfbmuy = dt/(2*h*h);
  
  mat_t *w = mat_init(u->xsize, u->ysize);
  
  // Y direction ---------------------------------------------------------------
  int m; for(m = 1; m < u->xsize-1; m++) {
    
    double *p = calloc(u->xsize, sizeof(double)); assert(p);
    double *q = calloc(u->xsize, sizeof(double)); assert(q);
    q[1] = 0; // Boundary values at x(0) and y(m)
    
    int el; for(el = 1; el < u->xsize-1; el++) {
      double dd = mat_get(el, m, u) + halfbmuy*(mat_get(el, m-1, u) -
                  2*mat_get(el, m, u) + mat_get(el, m+1, u));
      double denom = 1 + halfbmux*(2 - p[el]);
      p[el+1] = halfbmux / denom;
      q[el+1] = (dd + halfbmux*q[el]) / denom;
    }
    
    for (el = u->xsize-2; el >= 0; el--) {
      double val = p[el+1] * mat_get(el+1, m, w) + q[el+1];
      mat_set(el, m, val, w);
    }
    
  }

  // X direction ---------------------------------------------------------------
  int el; for(el = 1; el < u->xsize-1; el++) {
    
    double *p = calloc(u->xsize, sizeof(double)); assert(p);
    double *q = calloc(u->xsize, sizeof(double)); assert(q);
    q[1] = 0; // Boundary values at x(el) and y(0)
  
    int m; for(m = 1; m < u->xsize-1; m++) {
      double dd = mat_get(el, m, w) + halfbmux*(mat_get(el-1, m, w) -
                  2*mat_get(el, m, w) + mat_get(el+1, m, w));
      double denom = 1 + halfbmuy*(2 - p[m]);
      p[m+1] = halfbmux / denom;
      q[m+1] = (dd + halfbmux*q[m]) / denom;
    }
    
    for (m = u->xsize-2; m >= 0; m--) {
      double val = p[m+1] * mat_get(el, m+1, u) + q[m+1];
      mat_set(el, m, val, u);
    }
    
  }
  
  mat_free(w);
}
