/***************************************************
 * 2c3ppes-hj.c
 * 
 * By Yuan-Chung Cheng <yccheng@berkeley.edu>
 *
 * Compute two-color three-pulse photon-echo signal 
 * at Hohjai's negative T' time-ordering condition 
 * (i.e. 750, 800, 750) as a function of coherence 
 * time tau at a given population time T. See p. 8 on 
 * Notebook 2.
 * Using the JS03 bath module.
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
#include <gsl/gsl_integration.h>

#include "bath_js03.h"

#define CM2FS (5309.1)
#define NWSPACE (100000)
#define EPSABS (1e-5)
#define EPSREL (1e-5)

/* global variables */
const char *program_name;

/* Global variable given in the command line */
double beta, gamma1, wc1, gamma2, wc2;

/* Global variables that defines the system */
double theta; /* mixing angle between two chromophores */
/* two bath functions */
js03_bath_func bf1,bf2;
js03_lineshape_func gf1,gf2;
/* population time */
double T;

/* Display usage information and exit.  */
static void usage ()
{
  printf ("\
Usage: %s time T theta gamma1 wc1 gamma2 wc2 \n\
\n\
       time is the population time slice to be calculated (in fs).\n\
       T is the temperature (in K).\n\
       theta is the mixing angle between the two chromophores.\n\
       gamma1 and wc1 define bath for the lower-energy chromophore.\n\
       gamma2 and wc2 define bath for the higher-energy one.\n\
       Pulse sequence is scrambled downhill: w2, w1, w2.\n\
\n",
          program_name);
}

/* this is the kernel; the integrated photon-echo signal
             /infinity
  S(T,tau) = | dt |R2(t,T,tau)-R1*(t,T,tau)|^2
            /0
  the expression using lineshape functions are defined in p. 8 on Notebook 2.
*/
double f_S(double t,void *params)
{
  double tau;
  double Qab_t;
  double f0;
  double s4,c4,s2c2;
  double *p = (double *)params;

  tau=p[0];

  /* some useful constant factors */
  s4=sin(theta)*sin(theta)*sin(theta)*sin(theta);
  c4=cos(theta)*cos(theta)*cos(theta)*cos(theta);
  s2c2=sin(theta)*sin(theta)*cos(theta)*cos(theta);

  /* integrand is exp(-2*f0)*[1-cos(2*Qab(t))],
     where f0=Pbb(T+tau)-Pab(t+T+tau)+Pab(tau)+Pab(t)-Pab(T)+Paa(T+t)
     a -> alpha -> the lower energy eigenstate 
     b -> beta  -> the higher energy eigenstate
     Pulse sequence at this impulsive limit is (Eb,Ea,Eb) */
  f0=0;
  // Pbb(T+tau)
  f0=f0+s4*bath_js03_gt_r_cached(&gf1,T+tau)+c4*bath_js03_gt_r_cached(&gf2,T+tau);
  // -Pab(t+T+tau)
  f0=f0-s2c2*(bath_js03_gt_r_cached(&gf1,t+T+tau)+bath_js03_gt_r_cached(&gf2,t+T+tau));
  // Pab(tau)
  f0=f0+s2c2*(bath_js03_gt_r_cached(&gf1,tau)+bath_js03_gt_r_cached(&gf2,tau));
  // Pab(t)
  f0=f0+s2c2*(bath_js03_gt_r_cached(&gf1,t)+bath_js03_gt_r_cached(&gf2,t));
  // -Pab(T)
  f0=f0-s2c2*(bath_js03_gt_r_cached(&gf1,T)+bath_js03_gt_r_cached(&gf2,T));
  // Paa(T+t)
  f0=f0+c4*bath_js03_gt_r_cached(&gf1,T+t)+s4*bath_js03_gt_r_cached(&gf2,T+t);
  // img part
  Qab_t=s2c2*(bath_js03_gt_i_cached(&gf1,t)+bath_js03_gt_i_cached(&gf2,t));

  //  exp(-2*f0)*[1-cos(2*Qab(t))]
  //  printf("T=%.12f, tau=%.12f, t=%.12f, Integrand=%.12f (f0=%.12f,Qab_t=%.12f)\n",
  //	 T*CM2FS,tau*CM2FS,t*CM2FS,exp(-2.0*f0)*(1-cos(2.0*Qab_t)),f0,Qab_t);
  return (exp(-2.0*f0)*(1-cos(2.0*Qab_t)));
}

double S(double tau)
{
  gsl_integration_workspace *w;
  double result, error;
  double params[1];
  gsl_function F;

  double ra=0.0,rb=90.0; /* integration range for t */

  w=gsl_integration_workspace_alloc(NWSPACE);

  /* fill in params */
  params[0]=tau;
  F.function = &f_S;
  F.params = params;

  gsl_integration_qags(&F, ra, rb, EPSABS/100.0, EPSREL/100.0, NWSPACE, w, &result, &error);
  gsl_integration_workspace_free(w);

  return result;
}

/* Main program */
int main(int argc, char *argv[])
{
  int i;

  double tau;
  double shift,C0,CT,fT;

  /* Set program name for messages.  */
  program_name = argv[0]; 

  if(argc != 8) {
    usage();
    exit(EXIT_FAILURE);
  }
 
  T=atof(argv[1])/CM2FS;
  beta=1.4387/atof(argv[2]);
  theta=atof(argv[3]);
  gamma1=atof(argv[4]);
  wc1=atof(argv[5]);
  gamma2=atof(argv[6]);
  wc2=atof(argv[7]);

  printf("\n");
  printf("# beta   = %f (%f K)\n",beta,1.4387/beta);
  printf("# theta  = %f\n",theta);
  printf("# gamma1 = %f\n",gamma1);
  printf("# wc1    = %f\n",wc1);
  printf("# gamma2 = %f\n",gamma2);
  printf("# wc2    = %f\n",wc2);
  printf("\n");
 
  // initialize bath functions
  bath_js03_ct_init(&bf1,beta,gamma1,wc1);
  bath_js03_gt_init(&gf1,&bf1);
  bath_js03_ct_init(&bf2,beta,gamma2,wc2);
  bath_js03_gt_init(&gf2,&bf2);

  // tau from 0fs to 300fs with step size 3fs
  for(tau=0.0;tau<=300.0;tau=tau+3.0) {
    printf("tau=%.2f, T=%.2f, S=%.18f\n",tau,T*CM2FS,S(tau/CM2FS));
  }

}
