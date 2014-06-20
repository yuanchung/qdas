/***************************************************
 * bath_htgau.c
 *
 * Part of the qdas package;
 * High-Temperature GAUssian bath module.
 * This module provides necessary bath correlation
 * functions for HT Gaussian bath in the excitonic basis.
 *
 * By Yuan-Chung Cheng <yuanchung@ntu.edu.tw>
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
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf.h>

#include "qdas.h"
#include "aux.h"
#include "bath_htgau.h"


// convert time from cm unit to fs
#define CM2FS (5309.1)

/* Global variables that defines the Gaussian energy-gap 
   correlation function for each excitation */
double BathHTGauLambda[256];
double BathHTGauTau[256];

/* F(x) is a n approximation for exp(x^2)*[1-erf(x)] */
inline double bath_htgau_Fx(double x)
{
  return (2.0/(2.32*x+sqrt(1.52*x*x+4.0)));
}
/* the temperature factor is embedded in the wsquare factor */
double bath_htgau_wsquare(double lambda0,double tau0,double beta)
{
  double u,ret;
  u=beta/tau0;
  ret=2.0*lambda0/beta*(5.0*u/sqrt(M_PI)-2.0*u*u*bath_htgau_Fx(u)-4.0*u*u*bath_htgau_Fx(2.0*u)+bath_htgau_Fx(2.5*u));
  return ret;
}


/* lineshape functions */
inline double bath_htgau_gn_r(double t, double lambda, double tau, double beta)
{
  double dtmp1, dtmp2,wsquare;

  dtmp1=tau*tau;
  dtmp2=-0.5*dtmp1 + 0.5*dtmp1*exp(-1.0*t*t/tau/tau) + 
    0.5*sqrt(M_PI)*tau*t*gsl_sf_erf(t/tau);
  wsquare=bath_htgau_wsquare(lambda,tau,beta);
  return (wsquare*dtmp2);
}

inline double bath_htgau_gn_i(double t, double lambda, double tau)
{
  return (-0.5*sqrt(M_PI)*lambda*tau*gsl_sf_erf(t/tau));
}

/* vibrational contributions */
inline double bath_htgau_gv_r(double t, vibmodes *vibs,double beta)
{
  int j;
  double sum;
  double wj,Sj;

  sum=0.0;
  for(j=0;j<vibs->n;j++) {
    Sj=vibs->S[j];
    wj=vibs->omega[j];
    sum=sum+Sj*(coth(beta*wj/2.0)*(1.0-cos(wj*t)));
  }
  return sum;
}

inline double bath_htgau_gv_i(double t, vibmodes *vibs)
{
  int j;
  double sum;
  double wj,Sj;

  sum=0.0;
  for(j=0;j<vibs->n;j++) {
    Sj=vibs->S[j];
    wj=vibs->omega[j];
    sum=sum+Sj*sin(wj*t);//+Sj*wj*t; // last term takes care of vib. S shift
  }
  return sum;
}

/* interface functions */
int bath_htgau_init_params(const size_t nsize, const double beta, 
			   const size_t bath_nparams, const double *bath_params)
{
  int i,idx;

  printf("Use JSF's Gaussian Ansatz (high-temperature approx.):\n");
  printf("\n");
  printf("                 /t                    /t   /t1\n");
  printf("g(t) = -I*lambda*| M(t1)*dt1 + <dw^2>*|dt1 | M(t2)*dt2\n");
  printf("                 /0                   /0   /0\n");
  printf("\n");
  printf("M(t) = exp(-t^2/tau^2)\n");
  printf("\n");
  printf("For each site:\n");
  printf("%4s %12s %12s %12s\n","Site","lambda (cm^-1)","tau (fs)","<dw^2>");
  for(i=0;i<nsize;i++) {
    idx=i*2;
    BathHTGauLambda[i]=bath_params[idx];
    BathHTGauTau[i]=bath_params[idx+1];
    printf("%4d %12.4f %12.4f %12.4f\n",i+1,BathHTGauLambda[i],BathHTGauTau[i]*CM2FS,
	   bath_htgau_wsquare(BathHTGauLambda[i],BathHTGauTau[i],beta));
  }
  printf("\n");

  return 0;
}

int bath_htgau_free_params()
{
  return 0;
}

/*
 * $Log$
 * Revision 1.2  2007/01/30 08:33:15  platin
 *   - use more reasonable coupling constant factor in the bath_htgau module,
 *     note that this is still a high-T approximation.
 *   - more informative messages in 2c3ppes-htgau.c
 *
 * Revision 1.1.1.1  2006/05/24 00:42:18  platin
 *
 *   - initial import of the qdas package.
 *   - qdas stands for quantum dynamics and spectroscopy.
 *   - basic code inherited from the "lineshape" package.
 *
 *
 *
 */
