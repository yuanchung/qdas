/***************************************************
 * 1cpeps.c
 * 
 * By Yuan-Chung Cheng <yccheng@berkeley.edu>
 *
 * Compute one-color PEPS as a function of 
 * population time T. Using the OHMART bath module.
 *
 ***************************************************/

#define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_integration.h>

#include "bathtype_mt99art.h"
#include "bath_ohmart.h"

#define CM2FS (5309.1)
#define NWSPACE (100000)
#define EPSABS (1e-7)
#define EPSREL (1e-7)

/* global variables */
const char *program_name;

/* Global variable given in the command line */
double beta, gamma0, wc;

/* Global variable for the convenience of minimization */
/* population time */
double T;
/* bath function */
mt99art_bath_func bf;

/* Display usage information and exit.  */
static void usage ()
{
  printf ("\
Usage: %s beta gamma0 wc \n\
\n",
          program_name);
}

/* kernel for the integrated 1c photon-echo signal
   the expression is given in p. 41 on notebook 1
*/
double f_Sint(double tt,void *params)
{
  double tau,TT;
  double f0,q0;
  double *p = (double *)params;
  double t;

  tau=p[0];
  TT=p[1];
  t=tt/CM2FS;

  /* integrand is exp(-2*f0)*[cos(q0)^2],
     where f0=Paa(tau)-Paa(tau+T+t)+Paa(tau+T)+Paa(T+t)-Paa(T)+Paa(t)
           q0=Qaa(t)+Qaa(T)-Qaa(T+t) */

  f0=bath_mt99art_gt_r(&bf,tau)-bath_mt99art_gt_r(&bf,tau+TT+t)+bath_mt99art_gt_r(&bf,tau+TT)
    +bath_mt99art_gt_r(&bf,TT+t)-bath_mt99art_gt_r(&bf,TT)+bath_mt99art_gt_r(&bf,t);
  q0=bath_mt99art_gt_i(&bf,t)+bath_mt99art_gt_i(&bf,TT)-bath_mt99art_gt_i(&bf,TT+t);

  //  printf("tau=%ffs, T=%ffs, t=%ffs, f0=%f, q0=%f\n",
  //	 tau*CM2FS,TT*CM2FS,t*CM2FS,f0,q0);

  return (exp(-2.0*f0)*(cos(q0)*cos(q0)));
}

/* the integrated photon-echo signal at tau>0 and T>0

             /infinity
  S(tau,T) = | dt |R2(tau,T,t)+R3(tau,T,t)|^2
            /0
*/
double S_int(double ttau, void * null_params)
{
  gsl_integration_workspace *w;
  double result, error;
  double params[2];
  gsl_function F;
  double ra=0.0,rb=2000.0; /* integration range for t, in fs */
  double tau,TT;

  tau=ttau/CM2FS;
  TT=T/CM2FS;

  w=gsl_integration_workspace_alloc(NWSPACE);

  /* fill in params, tau and T */
  params[0]=tau;
  params[1]=TT;
  F.function = &f_Sint;
  F.params = params;

  gsl_integration_qags(&F, ra, rb, EPSABS/100.0, EPSREL/100.0, NWSPACE, w, &result, &error);
  gsl_integration_workspace_free(w);

  //  printf("S1\n");

  return -1.0*result; // return negative value
}

/* the approximated target S function as a function of tau, parameter is T.
   see p. 42 of Notebook 1 for the expression for S */
double S_approx(double ttau, void * params)
{
  double TT;
  double A,B,tau,tau2,tau4;
  double dtmp;
  double ret;

  tau=ttau/CM2FS;
  tau2=tau*tau;
  tau4=tau2*tau2;

  TT=T/CM2FS;

  A = GSL_REAL(bath_mt99art_ct(&bf,0.0)) +
    GSL_REAL(bath_mt99art_ct(&bf,TT)) - GSL_REAL(bath_mt99art_ct(&bf,TT+tau)) +
    bath_mt99art_ht_i(&bf,TT)*bath_mt99art_ht_i(&bf,TT);

  B = bath_mt99art_ht_r(&bf,TT+tau) - bath_mt99art_ht_r(&bf,TT);
  
  //  printf("T=%f, tau=%f, A=%f, B=%f\n",T,ttau,A,B);

  dtmp=-2.0*bath_mt99art_gt_r(&bf,tau)+B*B/A;

  // we return the value in negative sign because we will do "minimization"
  // to find the maximum
  ret=-1.0*(exp(dtmp)/sqrt(A)*(gsl_sf_erf(B/sqrt(A))+1));

  //  printf("S: tau=%f, S=%f\n",tau,ret);
  return ret;
}

/* minimize S and return the peak-shift value */
double peps(double (* S)(double,void *))
{                                                                          
  int status;                                                              
  int iter = 0, max_iter = 1000;                                            
  const gsl_min_fminimizer_type *MT;                                        
  gsl_min_fminimizer *s;                                                   
  double a = 0.0, b = 100.0; // locate minimum in this range
  gsl_function F; 

  double Sa, Sb, Stau0, tau0;

  // we first find a initial tau so that S(tau0) < S(a) && S(tau0) < S(b)
  Sa=S(a,NULL);
  Sb=S(b,NULL);
  tau0=30.0;
  Stau0=S(tau0,NULL);
  while (Stau0 > Sa || Stau0 > Sb) {
    if (Stau0 > Sa) {
      tau0 = (tau0 + a)/2.0;
    } else {
      tau0 = (tau0 + b)/2.0;
    }
    Stau0=S(tau0,NULL);
  }

  F.function = S;                                                       
  F.params = NULL;  
  MT = gsl_min_fminimizer_brent;
  s = gsl_min_fminimizer_alloc (MT);
  //  printf("HERE1\n");
  gsl_min_fminimizer_set (s, &F, tau0, a, b);
  //  printf("HERE2\n");
  
  //  printf ("using %s method\n",gsl_min_fminimizer_name (s));
  do {                                                                      
    iter++;                                                              
    //    printf("iter=%d, %f, %f, %f\n",iter,tau0,a,b);
    status = gsl_min_fminimizer_iterate (s);                             
    
    tau0 = gsl_min_fminimizer_x_minimum (s);                                
    a = gsl_min_fminimizer_x_lower (s);                                  
    b = gsl_min_fminimizer_x_upper (s);                                  
    
    status = gsl_min_test_interval (a, b, 0.01, 0.0);
    if (status == GSL_SUCCESS) {
      gsl_min_fminimizer_free(s);
      return tau0;
    }
  } while (status == GSL_CONTINUE && iter < max_iter);                       
  
  gsl_min_fminimizer_free(s);
  return 0.0;
}

/* Main program */
int main(int argc, char *argv[])
{
  int i;

  double tau;
  double shift,C0,CT,fT;
  gsl_matrix *tmp;

  tmp=gsl_matrix_alloc(2,2); // actually redundant...

  /* Set program name for messages.  */
  program_name = argv[0]; 

  if(argc != 4) {
    usage();
    exit(EXIT_FAILURE);
  }
 
  beta=atof(argv[1]);
  gamma0=atof(argv[2]);
  wc=atof(argv[3]);

  printf("\n");
  printf("beta   = %f (%f K)\n",beta,1.4387/beta);
  printf("gamma0 = %f\n",gamma0);
  printf("wc     = %f\n",wc);
  printf("\n");
 
  bath_ohmart_ct_init(&bf,beta,gamma0,wc,tmp);

  T=0.0;
  for(tau=-10.0;tau<100.0;tau=tau+1.0) {
    printf("tau=%f, S=%f\n",tau,S_approx(tau,NULL));
  }

  T=0.0;
  for(tau=-10.0;tau<100.0;tau=tau+1.0) {
    printf("tau=%f, Sint=%20.12f\n",tau,S_int(tau,NULL));
  }

  C0=GSL_REAL(bath_mt99art_ct(&bf,0.0));
  printf("C(0)=%f\n",C0);
  printf("\n");
  for(T=0.0;T<1000;T=T+5.0) {
    CT=GSL_REAL(bath_mt99art_ct(&bf,T/CM2FS));
    fT=bath_mt99art_ht_i(&bf,T/CM2FS)*bath_mt99art_ht_i(&bf,T/CM2FS);
    shift=CT/sqrt(M_PI)/C0/sqrt(C0+fT)*CM2FS; // directly from the correlation function
    printf("T=%f PEPS=%20.12f aPEPS=%20.12f aaPEPS=%20.12f\n",T,peps(S_int),peps(S_approx),shift);
  }

  return 0;
}
