/***************************************************
 * qng2d.c
 * 
 * By Yuan-Chung Cheng <yccheng@mit.edu>
 *
 * functions to do numerical integration by using
 * quadrature methods. These functions do not really 
 * do the intergations, but call numerical integration
 * functions in GNU Scientific Library to process the 
 * integration. Actually these are wrappers...
 *
 ***************************************************/

#include <stdio.h>
#include <math.h>
#include <errno.h>
#include <gsl/gsl_integration.h>
#include "qng2d.h"

/*
 * Functional definition
 */


/* 
 *  global variables in this file
 */


/***************************************************
    START OF FUNCTIONS
***************************************************/

/*
 * 2D version
 */

/*
  
  int qng2d(DoubleFunc1d func, 
            double xlower, double xupper,
	    double ylower, double yupper,
	    double abseps,double releps,
	    double *result, double *abserr);

  This function calculate the following integral by GSL qng method:

                 yupper xupper
	             / /
     *result =       | | func(x,y) dxdy
                     / /
                 ylower xlower

  Return:
    
 */

/* construce the integrands */

DoubleFunc2d Func_QNG2D; /* the to-be-integrate 2d function */

double F_Qng2d_FixY(double x, void * params) 
{
  double y;

  y=*(double *)params;

#ifdef DEBUG_QNG2D
  printf("Func_QNG2D(%18.12f,%18.12f) = %18.12f\n",x,y,Func_QNG2D(x,y));
#endif

  return Func_QNG2D(x,y);
}

/* 2D integrand,, integrate over X on a fixed Y value */
/* 
 * params should give the parameters to integrate over X:
 * params[0]: xlower;
 * params[1]: xupper;
 * params[2]: abseps;
 * params[3]: releps;
 *
 */
double F_Qng2d(double y, void * params)
{
  double *array;
  double result,abserr;
  size_t neval;

  double valy;
  size_t returnval;
  gsl_function F;
  
  valy=y;
  array=(double *)params;
  F.function=F_Qng2d_FixY;
  F.params=&valy;

  gsl_integration_qng (&F, array[0], array[1], array[2], array[3], &result, &abserr, &neval);

#ifdef DEBUG_QNG2D
  fprintf(stderr,"F_Qng2d: returnval =  %d\n", returnval);
  fprintf(stderr,"F_Qng2d: intervals =  %d\n", neval);
  fprintf(stderr,"F_Qng2d: result    =  % .18f\n", result);
  fprintf(stderr,"F_Qng2d: abserr    =  % .18f\n", abserr);
#endif

  return result;
}


size_t qng2d(DoubleFunc2d func, double xlower, double xupper,
	     double ylower, double yupper,
	     double abseps,double releps, double *result, double *abserr)
{
  gsl_function F;
  
  int returnval;
  size_t neval;
  double params[4];
  
  /* set the 2d integrand globally */
  Func_QNG2D = func;

  /* init the function */
  F.function=F_Qng2d;
  params[0]=xlower;
  params[1]=xupper;
  params[2]=abseps;
  params[3]=releps;
  F.params=params;

  returnval=gsl_integration_qng (&F, ylower, yupper, abseps, releps, result, abserr, &neval);
  
#ifdef DEBUG_QNG2D
  fprintf(stderr,"qng2d: returnval =  %d\n", returnval);
  fprintf(stderr,"qng2d: intervals =  %d\n", neval);
  fprintf(stderr,"qng2d: result    =  % .18f\n", *result);
  fprintf(stderr,"qng2d: abserr    =  % .18f\n", *abserr);
#endif
  
  return returnval;
}

/*
 * TEST SECTION
 *
 */

#ifdef MAIN
double ff(double x, double y) {
  return cos(5.0*x)*cos(3.0*y)*exp(-1.0*x*x-1.0*y*y);
}

int
main ()
{
  double result, error;
  double expected = 0.000160292230954833645;

  qng2d(ff, 0, M_PI, 0, M_PI, 0.0, 1e-3, &result, &error);
  
  printf("result          = % .18f\n", result);
  printf("exact result    = % .18f\n", expected);
  printf("estimated error = % .18f\n", error);
  printf("actual error    = % .18f\n",  result - expected);
}

#endif

/*
 * $Log$
 * Revision 1.1  2007/01/03 07:30:00  platin
 *   - add a new package in which pulse duration integral is carried out
 *     to include pulse overlap effect using impulsive response
 *     function formalism.
 *
 *
 */

