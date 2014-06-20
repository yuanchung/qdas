
/***************************************************
 * gau_max.c
 *
 * Estimate the maximum of a curve using least-square
 * fitting to a Gaussian function.
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
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multimin.h>

/* convergence criterion */
#define EPSABS (0.001)

/* Maxinum number of points, size of cache. */
#define MAX_NPOINTS  10000

/* global variables; the data points are made global
   so that the minimization kernel can use them. */
gsl_vector *data_x;
gsl_vector *data_y;

/* we fit to a gaussian f(x)=A0*exp(-(x-x0)*(x-x0)/2/sigma/sigma) */
inline double gau_fx(double x,double a0, double x0, double sigma) {
  return a0*exp(-(x-x0)*(x-x0)/2.0/sigma/sigma);
}

/* target function; we minimize \sum_{i=1}^n |Gau(x_i,A0,x0,sigma)-y_i|^2
   to obtain best fit A0,x0, and sigma */
double target_func (const gsl_vector *v, void *empty_params)
{
  int n,i;
  double A0,x0,sigma;
  double sum,xi,yi,fx;

  A0=gsl_vector_get(v,0);
  x0=gsl_vector_get(v,1);
  sigma=gsl_vector_get(v,2);

  n=data_x->size;
  sum=0.0;
  for(i=0;i<n;i++) {
    xi=gsl_vector_get(data_x,i);
    yi=gsl_vector_get(data_y,i);
    fx=gau_fx(xi,A0,x0,sigma);
    sum=sum+(fx-yi)*(fx-yi);
  }
  return sum;
}

int main()
{
  char buffer[256];

  double x[MAX_NPOINTS],y[MAX_NPOINTS];
  double max;

  double a0min,x0min,sigmamin;

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

  /* initialize data vector */
  data_x=gsl_vector_alloc(N);
  data_y=gsl_vector_alloc(N);
  for(i=0;i<N;i++) {
    gsl_vector_set(data_x,i,x[i]);
    gsl_vector_set(data_y,i,y[i]);
  }

  /* minimization procedure, we first find the 
     first guess for x0, the maximal point in the set */
  kmax=0;
  max=y[0];
  for(i=0;i<N;i++) {
    if(y[i]>max) {
      kmax=i;
      max=y[i];
    }
  }

  /* now ready to initialize and invoke the GSL minimization process */
  {
    size_t np = 3; // number of variables to fit     

    const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex;
    gsl_multimin_fminimizer *s = NULL;
    gsl_vector *ss, *xx;
    gsl_multimin_function minex_func;
     
    size_t iter = 0;
    int status;
    double size;
     
    /* Initial vertex size vector */
    ss = gsl_vector_alloc (np);
    /* Set all step sizes to 1 */
    gsl_vector_set_all (ss, 1.0);
     
    /* Starting point */
    xx = gsl_vector_alloc (np);
    gsl_vector_set (xx, 0, y[kmax]); // initial guess for A0
    gsl_vector_set (xx, 1, x[kmax]); // initial guess for x0
    gsl_vector_set (xx, 2, 30.0); // initial guess for sigma
     
    /* Initialize method and iterate */
    minex_func.f = &target_func;
    minex_func.n = np;
    minex_func.params = NULL;
     
    s = gsl_multimin_fminimizer_alloc (T, np);
    gsl_multimin_fminimizer_set (s, &minex_func, xx, ss);
    do {
      iter++;
      status = gsl_multimin_fminimizer_iterate(s);
     
      if (status)
	break;
      
      size = gsl_multimin_fminimizer_size (s);
      status = gsl_multimin_test_size (size, EPSABS);
     
#ifdef DEBUGDEBUG
      if (status == GSL_SUCCESS) {
	printf ("converged to minimum at\n");
      }
      printf ("%5d ", iter);
      for (i = 0; i < np; i++) {
	printf ("%10.3e ", gsl_vector_get (s->x, i));
      }
      printf ("f() = %7.3f size = %.3f\n", s->fval, size);
#endif
    } while (status == GSL_CONTINUE && iter < 1000);
    
    a0min=gsl_vector_get (s->x,0);
    x0min=gsl_vector_get (s->x,1);
    sigmamin=gsl_vector_get (s->x,2);

    gsl_vector_free(xx);
    gsl_vector_free(ss);
    gsl_multimin_fminimizer_free (s);
     
  }

  /* print (x,y) pair at the maximum */
  printf("%f %f\n",x0min,gau_fx(x0min,a0min,x0min,sigmamin));

  /* done; clean up */
  gsl_vector_free(data_x);
  gsl_vector_free(data_y);

  return 1;

}

/*
 * $Log$
 * Revision 1.2  2007/03/26 19:10:34  platin
 *   - remove unnessary debug message.
 *
 * Revision 1.1  2006/12/05 07:33:18  platin
 * *** empty log message ***
 *
 *
 */
