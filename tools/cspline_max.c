
/***************************************************
 * cspline_max.c
 *
 * Estimate the maximum of a curve using c-spline
 * interpolation.
 *
 * By Yuan-Chung Cheng <yccheng@berkeley.edu>
 *
 *
 ***************************************************/

#define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_min.h>

/* convergence criterion */
#define EPSABS (0.01)

/* Maxinum number of points, size of cache. */
#define MAX_NPOINTS  10000

/* global variables; the cspline variables are made global
   so that the minimization kernel can use the cspline function. */
gsl_interp_accel *cspline_max_acc;
gsl_spline *cspline_max_spline;

/* minimization kernel; a cspline representation 
   of the given data curve. */
double func_y (double x, void * params)
{
  return -1.0*gsl_spline_eval(cspline_max_spline,x,cspline_max_acc);
}

int main()
{
  char buffer[256];

  double x[MAX_NPOINTS],y[MAX_NPOINTS];
  double max;

  double xmin;

  int i,N;
  int kmax;

  // read data straight from stdin
  N=0;
  while(1) {
    fscanf(stdin,"%s",buffer);
    if(feof(stdin)) break;
    x[N]=strtod(buffer,NULL);
    buffer[0]='\0'; /* clear buffer in case the next read fails */
    fscanf(stdin,"%s",buffer);
    y[N]=strtod(buffer,NULL);
    N++;
  }

  if(N<=3) {
    printf("Not enough number of data points!!\n");
    exit(EXIT_FAILURE);
  }

#ifdef DEBUGDEBUG
  printf("Read %d points:\n",N);
  for(i=0;i<N;i++) {
    printf("x=%f, y=%f\n",x[i],y[i]);
  }
#endif

  /* initialize cspline function */
  cspline_max_acc=gsl_interp_accel_alloc();
  cspline_max_spline=gsl_spline_alloc(gsl_interp_cspline, N);
  gsl_spline_init(cspline_max_spline, x, y, N);

#ifdef DEBUGDEBUG
  {
    double xx,dxx;
    xx=x[0];
    dxx=(x[N-1]-x[0])/100.0;
    printf("C-Spline approximation:\n");
    for(i=0;i<100;i++) {
      printf("xx=%f, yy=%f\n",xx,gsl_spline_eval(cspline_max_spline,xx,cspline_max_acc));
      xx=xx+dxx;
    }
  }
#endif

  /* minimization procedure, we first find the 
     first guess, the maximal point in the set */
  kmax=0;
  max=y[0];
  for(i=0;i<N;i++) {
    if(y[i]>max) {
      kmax=i;
      max=y[i];
    }
  }
  xmin=x[kmax]; // initial guess of minimum/maximum

  /* now ready to initialize and invoke the GSL minimization process */
  {
    int status;
    int iter=0,max_iter=100;
    const gsl_min_fminimizer_type *T;
    gsl_min_fminimizer *s;
    gsl_function F;

    double xlower,xupper;

    F.function=&func_y;
    F.params=NULL;

    T=gsl_min_fminimizer_brent;
    s=gsl_min_fminimizer_alloc (T); 
    gsl_min_fminimizer_set (s, &F, xmin, x[0], x[N-1]);
    
    do {
      iter++;
      status=gsl_min_fminimizer_iterate(s);
      
      xmin=gsl_min_fminimizer_x_minimum(s);
      xlower=gsl_min_fminimizer_x_lower(s);
      xupper=gsl_min_fminimizer_x_upper(s);

      status=gsl_min_test_interval(xlower,xupper,EPSABS,0.0);

    } while (status == GSL_CONTINUE && iter < max_iter);
  }

  /* print (x,y) pair at the maximum */
  printf("%f %f\n",xmin,gsl_spline_eval(cspline_max_spline,xmin,cspline_max_acc));

  /* done; clean up */
  gsl_spline_free(cspline_max_spline);
  gsl_interp_accel_free(cspline_max_acc);

  return 1;

}

/*
 * $Log$
 * Revision 1.2  2007/03/07 22:44:19  platin
 *   - not show debug details.
 *
 * Revision 1.1  2006/07/12 17:32:42  platin
 *
 *   - cspline_max: use cspline interpolation to estimate maximum points.
 *   - normalize: normalize data set so that maximum of y is 1.0.
 *
 *
 */
