// area.c

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h> 

# define NPOINTS 2000
# define MAXITER 2000


struct complex{
  double real;
  double imag;
};

int main(){
  int i, j, iter, numoutside = 0;
  double area, error, ztemp, time1, time2, elapsed;
  struct complex z, c;

/*
 *   
 *
 *     Outer loops run over npoints, initialise z=c
 *
 *     Inner loop has the iteration z=z*z+c, and threshold test
 */
 time1 = omp_get_wtime();
#pragma omp parallel default (none) \
                         private(i, j, iter, c, z, ztemp) \
                         shared(numoutside)
{
#pragma omp single
{
  for (i=0; i<NPOINTS; i++) {
    #pragma omp task
    {
      for (j=0; j<NPOINTS; j++) {
        c.real = -2.0+2.5*(double)(i)/(double)(NPOINTS)+1.0e-7;
        c.imag = 1.125*(double)(j)/(double)(NPOINTS)+1.0e-7;
        z=c;
        for (iter=0; iter<MAXITER; iter++){
          ztemp=(z.real*z.real)-(z.imag*z.imag)+c.real;
          z.imag=z.real*z.imag*2+c.imag; 
          z.real=ztemp; 
          if ((z.real*z.real+z.imag*z.imag)>4.0e0) {
#pragma omp atomic            
            numoutside++; 
            break;
          }
        }
} /* task */
      }
    }
}
}
time2 = omp_get_wtime();
elapsed = time2 - time1;
printf("Elsapsed time: %f\n", elapsed);

/*
 *  Calculate area and error and output the results
 */

      area=2.0*2.5*1.125*(double)(NPOINTS*NPOINTS-numoutside)/(double)(NPOINTS*NPOINTS);
      error=area/(double)NPOINTS;

      printf("Area of Mandlebrot set = %12.8f +/- %12.8f\n",area,error);

  }
