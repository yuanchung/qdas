/***************************************************
 * 2c3ppes.c
 * 
 * By Yuan-Chung Cheng <yccheng@berkeley.edu>
 *
 * Compute full two-color three-pulse photon-echo signal 
 * including different time-ordering conditions
 * (i.e. 750, 800, 750) as a function of coherence 
 * time tau at a given population time T. See p. Eq. 1-3
 * on Notebook 2.
 *
 * Using the JS03 bath module.
 *
 * Notations:
 * T: population time in the "naturally ordered" exp.
 * tau: coherence time in the naturally ordered exp.
 * t1: time delay between the first and the second pulse, t1>0
 * t2: time delay between the second and the third pulse, t2>0
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
#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp.h>

#include "aux.h"
#include "bath_js03.h"

#define CM2FS (5309.1)
#define NWSPACE (100000)
#define EPSABS (1e-6)
#define EPSREL (1e-6)
#define CSCALE (1e3)
#define LSF_CACHE_MAXPOINTS (30000)
#define LSF_CACHE_MIN (0.0)
#define LSF_CACHE_MAX (5.0)
#define LSF_CACHE_STEPFACTOR (0.0001)

/* typedef for lineshape functions, including all contributions */
typedef struct {
  js03_lineshape_func gn; /* bath lineshape function, the homogeneous contribution */
  double S,w,gg; /* coherent vibrational mode with Huang-Rhys factor S and frequency w and decay constant gg */
  double Delta; /* std. dev. of Gaussian static disorder */
} lineshape_func;

/* global variables */
const char *program_name;

/* Global variable given in the command line */
double theta; /* mixing angle between two chromophores */
double beta, gamma1, wc1, gamma2, wc2,gammaj,wcj;

/* Global variables that defines the system */
/* two bath functions */
js03_bath_func bf1,bf2,bfj;
lineshape_func gf1,gf2,gfj;
/* transition dipole moments |a>->|B>, |b>->|H> */

/* BL and HL from GD paper */
const double mu_a[3]={0.78550,  0.44743,  0.42755};
const double mu_b[3]={0.18975*0.95,  -0.96874*0.95,  -0.15979*0.95};
/* BM and HM from GD paper */
//const double mu_a[3]={-0.9544328,   0.0099420,   0.2982603};
//const double mu_b[3]={0.28130,  -0.32148,  -0.90417};

/* kinetic rates; note that these are only effective in TG calculations */
double Kbh; /* |H> -> |B> */
double Kpb; /* |B> -> |P> */

/* Display usage information and exit.  */
static void usage ()
{
  printf ("\
Usage: %s temp theta gamma1 wc1 Delta1 gamma2 wc2 Delta2 [gammaj wcj Deltaj S1 w1 gg1 S2 w2 gg2 Sj wj ggj Rbh Rpb]\n\
\n\
       temp is the temperature (in K).\n\
       theta is the mixing angle between the two chromophores.\n\
       gamma1 and wc1 define bath for the lower-energy chromophore.\n\
       gamma2 and wc2 define bath for the higher-energy one.\n\
       Delta1 and Delta2 are standard dev. of Gaussian static disorder\n\
       for chromophore 1 and 2, respectively.\n\
       Optional argvs:\n\
       gammaj, wcj, and Deltaj are for off-diagonal couplings.\n\
       S[12] and w[12] are vibrational Huang-Rhys factor and mode\n\
       frequency for chromophore 1 and 2, respectively.\n\
       gg[12] are the decay rate (in cm^-1) for the vib. modes.\n\
       Sj, wj, and ggj are for off-diagonal couplings.\n\
       Rbh is |H> -> |B> population transfer rate in fs^-1.\n\
       Rpb is |B> -> |P> population transfer rate in fs^-1.\n\
\n",
          program_name);
}

/* compute transition dipole moment of mixed states,
   dm = c1*mu_a + c2*mu_b, where mu_a and mu_b are
   defined as global variables.

*/
void compute_mixed_dipole(double c1,double c2,double *dm)
{
  int i;

  for(i=0;i<3;i++) {
    dm[i]=c1*mu_a[i]+c2*mu_b[i];
  }
}

inline double vec_dot_product(double *a, double *b)
{
  int i;
  double sum;

  sum=0.0;
  for(i=0;i<3;i++) {
    sum=sum+a[i]*b[i];
  }
  return sum;
}

/* orientational ensemble average of dipole factors,
   <MaMbMcMd> = 1/15 * { (Ma*Mb)(Mc*Md) + (Ma*Mc)(Mb*Md) + (Ma*Md)(Mb*Mc) } */
double compute_average_dipole(double *Ma, double *Mb, double *Mc, double *Md)
{
  double ab,cd,ac,bd,ad,bc;

  ab=vec_dot_product(Ma,Mb);
  cd=vec_dot_product(Mc,Md);
  ac=vec_dot_product(Ma,Mc);
  bd=vec_dot_product(Mb,Md);
  ad=vec_dot_product(Ma,Md);
  bc=vec_dot_product(Mb,Mc);
  // we do not apply the 1/15 factor
  return ((ab*cd+ac*bd+ad*bc));
}

/* compute ensemble angle averaged dipole factors needed for 
   peak shift calculations;
   exciton states: |alpha> = |e1>, |beta>=|e2>, |f>=|e1e2>
*/
double Me1Me1Me2Me2;
double Me2Me2Me2fMe2f;
double Me1Me2Me1fMe2f;
/* transition dipole factors needed for PEPS calculations:
   dfA = d_alpha^4*d_beta^4+d_alpha^2*d_beta^2*d_alpha->mu^2*d_beta->mu^2
   dfB = d_alpha^3*d_beta^3*d_alpha->mu*d_beta->mu
*/
double dfA_1,dfB_1; // for region 1
double dfA_2,dfB_2; // for region 2
double dfA_3,dfB_3; // for region 3
void compute_dipole_factors()
{
  double Me1[3],Me2[3],Me1f[3],Me2f[3];

  compute_mixed_dipole(cos(theta),sin(theta),Me1);
  compute_mixed_dipole(-1.0*sin(theta),cos(theta),Me2);
  compute_mixed_dipole(sin(theta),cos(theta),Me1f);
  compute_mixed_dipole(cos(theta),-1.0*sin(theta),Me2f);

  Me1Me1Me2Me2 = compute_average_dipole(Me1,Me1,Me2,Me2);
  Me2Me2Me2fMe2f = compute_average_dipole(Me2,Me2,Me2f,Me2f);
  Me1Me2Me1fMe2f = compute_average_dipole(Me1,Me2,Me1f,Me2f);

  /* scale for numerical fittness */
  Me1Me1Me2Me2=Me1Me1Me2Me2*CSCALE;
  Me2Me2Me2fMe2f=Me2Me2Me2fMe2f*CSCALE;
  Me1Me2Me1fMe2f=Me1Me2Me1fMe2f*CSCALE;

  dfA_1=Me1Me1Me2Me2*Me1Me1Me2Me2+Me2Me2Me2fMe2f*Me2Me2Me2fMe2f;
  dfB_1=Me1Me1Me2Me2*Me2Me2Me2fMe2f;

  dfA_2=Me1Me1Me2Me2*Me1Me1Me2Me2+Me1Me2Me1fMe2f*Me1Me2Me1fMe2f;
  dfB_2=Me1Me1Me2Me2*Me1Me2Me1fMe2f;

  dfA_3=dfA_2;
  dfB_3=dfB_2;
}

/* lineshape functions including inhomogeneous and vib. contributions */
double gt_r(lineshape_func *gt,double t)
{
  double ret;
  double S,w,Delta,g,tmp;

  S=gt->S;
  w=gt->w;
  g=gt->gg;
  Delta=gt->Delta;

  //  printf("S=%20.12f, w=%20.12f, Delta=%20.12f\n",S,w,Delta);
  ret=bath_js03_gt_r_cached(&(gt->gn),t);
  //  ret=ret+S*(coth(beta*w/2.0)*(1.0-cos(w*t)));
  // ad-hoc dissipative oscillator
  tmp=exp(-1.0*g*t);
  ret=ret+S*w*w*coth(beta*w/2.0)*pow(g*g+w*w,-2.0)*(w*w+g*w*w*t-w*w*tmp*cos(w*t)-2.0*g*tmp*w*sin(w*t)+g*g*g*t+g*g*tmp*cos(w*t)-g*g);
  ret=ret+0.5 * Delta * Delta * t * t;
  //  printf("gt_r=%20.12f, t=%20.12f\n",ret,t);
  return ret;
}

double gt_i(lineshape_func *gt,double t)
{
  double ret;
  double S,w,g,tmp;

  S=gt->S;
  w=gt->w;
  g=gt->gg;

  ret=bath_js03_gt_i_cached(&(gt->gn),t);
  //ret=ret+S*(sin(w*t)-w*t);
  // ad-hoc dissipative oscillator
  tmp=exp(-1.0*g*t);
  ret=ret+-1.0*S*w*w*pow(g*g+w*w,-2.0)*(-2.0*g*w+w*g*g*t+w*w*w*t+2*w*tmp*g*cos(w*t)-w*w*tmp*sin(w*t)+g*g*tmp*sin(w*t));

  //  printf("gt_i=%20.12f, t=%20.12f\n",ret,t);
  return ret;
}

/* cached lineshape functions for exciton states */
#define Paa(tt) gsl_spline_eval(Paa_Spline,tt,Paa_Spline_Acc)
gsl_spline *Paa_Spline;
gsl_interp_accel *Paa_Spline_Acc;
#define Qaa(tt) gsl_spline_eval(Qaa_Spline,tt,Qaa_Spline_Acc)
gsl_spline *Qaa_Spline;
gsl_interp_accel *Qaa_Spline_Acc;
#define Pbb(tt) gsl_spline_eval(Pbb_Spline,tt,Pbb_Spline_Acc)
gsl_spline *Pbb_Spline;
gsl_interp_accel *Pbb_Spline_Acc;
#define Qbb(tt) gsl_spline_eval(Qbb_Spline,tt,Qbb_Spline_Acc)
gsl_spline *Qbb_Spline;
gsl_interp_accel *Qbb_Spline_Acc;
#define Pab(tt) gsl_spline_eval(Pab_Spline,tt,Pab_Spline_Acc)
gsl_spline *Pab_Spline;
gsl_interp_accel *Pab_Spline_Acc;
#define Qab(tt) gsl_spline_eval(Qab_Spline,tt,Qab_Spline_Acc)
gsl_spline *Qab_Spline;
gsl_interp_accel *Qab_Spline_Acc;
/* initialize [PQ]{aa,bb,ab} spline functions;
   call this to prepare all [PQ]{aa,bb,ab}(t) functions */
void exciton_PQ_init()
{
  double s4,c4,s2c2;

  double xx[LSF_CACHE_MAXPOINTS];
  double yy[LSF_CACHE_MAXPOINTS];
  double tt;
  size_t i,npoints;

  printf("Initializing cached lineshape functions for exciton states:\n");

  /* some useful constant factors */
  s4=sin(theta)*sin(theta)*sin(theta)*sin(theta);
  c4=cos(theta)*cos(theta)*cos(theta)*cos(theta);
  s2c2=sin(theta)*sin(theta)*cos(theta)*cos(theta);

  /* set up grid points to cache */
  i=0;tt=LSF_CACHE_MIN;
  while(tt<=LSF_CACHE_MAX && i<LSF_CACHE_MAXPOINTS) {
      xx[i]=tt;
      i++;
      if(tt<=0.2) {
	/* 1/LSF_CACHE_STEPFACTOR points between 0.0-0.2; this is about 1000 fs */
        tt = tt + LSF_CACHE_STEPFACTOR*0.2;
      } else {
	/* exponential step size after tt>0.2 */
        tt = tt*1.05;
      }
  }
  npoints=i;

  /* Paa(t) */
  printf("  Initializing Paa(t)...\n");
  for(i=0;i<npoints;i++) {
    yy[i]=c4*gt_r(&gf1,xx[i])+s4*gt_r(&gf2,xx[i])+4.0*s2c2*gt_r(&gfj,xx[i]);
  }
  Paa_Spline_Acc=gsl_interp_accel_alloc();
  Paa_Spline=gsl_spline_alloc (gsl_interp_cspline,npoints);
  gsl_spline_init(Paa_Spline,xx,yy,npoints);
  /* Qaa(t) */
  printf("  Initializing Qaa(t)...\n");
  for(i=0;i<npoints;i++) {
    yy[i]=c4*gt_i(&gf1,xx[i])+s4*gt_i(&gf2,xx[i])+4.0*s2c2*gt_i(&gfj,xx[i]);
  }
  Qaa_Spline_Acc=gsl_interp_accel_alloc();
  Qaa_Spline=gsl_spline_alloc (gsl_interp_cspline,npoints);
  gsl_spline_init(Qaa_Spline,xx,yy,npoints);
  /* Pbb(t) */
  printf("  Initializing Pbb(t)...\n");
  for(i=0;i<npoints;i++) {
    yy[i]=s4*gt_r(&gf1,xx[i])+c4*gt_r(&gf2,xx[i])+4.0*s2c2*gt_r(&gfj,xx[i]);
  }
  Pbb_Spline_Acc=gsl_interp_accel_alloc();
  Pbb_Spline=gsl_spline_alloc (gsl_interp_cspline,npoints);
  gsl_spline_init(Pbb_Spline,xx,yy,npoints);
  /* Qbb(t) */
  printf("  Initializing Qbb(t)...\n");
  for(i=0;i<npoints;i++) {
    yy[i]=s4*gt_i(&gf1,xx[i])+c4*gt_i(&gf2,xx[i])+4.0*s2c2*gt_i(&gfj,xx[i]);
  }
  Qbb_Spline_Acc=gsl_interp_accel_alloc();
  Qbb_Spline=gsl_spline_alloc (gsl_interp_cspline,npoints);
  gsl_spline_init(Qbb_Spline,xx,yy,npoints);
  /* Pab(t) */
  printf("  Initializing Pab(t)...\n");
  for(i=0;i<npoints;i++) {
    yy[i]=s2c2*(gt_r(&gf1,xx[i])+gt_r(&gf2,xx[i]))-4.0*s2c2*gt_r(&gfj,xx[i]);
  }
  Pab_Spline_Acc=gsl_interp_accel_alloc();
  Pab_Spline=gsl_spline_alloc (gsl_interp_cspline,npoints);
  gsl_spline_init(Pab_Spline,xx,yy,npoints);
  /* Qab(t) */
  printf("  Initializing Qab(t)...\n");
  for(i=0;i<npoints;i++) {
    yy[i]=s2c2*(gt_i(&gf1,xx[i])+gt_i(&gf2,xx[i]))-4.0*s2c2*gt_r(&gfj,xx[i]);
  }
  Qab_Spline_Acc=gsl_interp_accel_alloc();
  Qab_Spline=gsl_spline_alloc (gsl_interp_cspline,npoints);
  gsl_spline_init(Qab_Spline,xx,yy,npoints);

  printf("\n");

}

/* kernel for the integrated photon-echo signal at tau>0 and T>0
   t1=tau, t2=T
  the expression using lineshape functions are defined in Eq.3 on Notebook 2.
*/
double f_S1(double t,void *params)
{
  double t1,t2;
  double Qab_t2pt;
  double Qab_t2;
  double f0;
  double *p = (double *)params;

  t1=p[0];
  t2=p[1];

  /* integrand is exp(-2*f0)*[dfA-2*dfB*cos(2*(Qab(t2+t)-Qab(t2)))],
     where f0=Pbb(t1)-Pab(t+t1+t2)+Pab(t1+t2)+Pab(t2+t)-Pab(t2)+Paa(t)
     a -> alpha -> the lower energy eigenstate 
     b -> beta  -> the higher energy eigenstate
     Pulse sequence at this impulsive limit is (Eb,Eb,Ea) */
  f0=Pbb(t1)-Pab(t+t1+t2)+Pab(t1+t2)+Pab(t2+t)-Pab(t2)+Paa(t);

  // img part
  Qab_t2=Qab(t2);
  Qab_t2pt=Qab(t+t2);

  //  printf("f0=%20.12f\n",f0);
  return (exp(-2.0*f0)*(dfA_1-2.0*dfB_1*cos(2.0*(Qab_t2pt-Qab_t2))));
}

/* the integrated photon-echo signal at tau>0 and T>0

             /infinity
  S(t1,t2) = | dt |R3(t1,t2,t)-R1*(t1,t2,t)|^2
            /0
*/
double S1(double t1, double t2)
{
  gsl_integration_workspace *w;
  double result, error;
  double params[2];
  gsl_function F;

  double ra=0.0,rb=2.0; /* integration range for t */

  w=gsl_integration_workspace_alloc(NWSPACE);

  /* fill in params */
  params[0]=t1;
  params[1]=t2;
  F.function = &f_S1;
  F.params = params;

  gsl_integration_qags(&F, ra, rb, EPSABS/100.0, EPSREL/100.0, NWSPACE, w, &result, &error);
  gsl_integration_workspace_free(w);

  //  printf("S1\n");

  return result;
}

/* kernel for the integrated photon-echo signal at T<0, and tau > |T|
   t1=tau-|T|, t2=|T|
  the expression using lineshape functions are defined in Eq.1 on Notebook 2.
*/
double f_S2(double t,void *params)
{
  double t1,t2;
  double Qab_t;
  double f0;
  double *p = (double *)params;

  t1=p[0];
  t2=p[1];

  /* integrand is exp(-2*f0)*[dfA-2*dfB*cos(2*Qab(t))],
     where f0=Pbb(t1+t2)-Pab(t+t1+t2)+Pab(t1)+Pab(t)-Pab(t2)+Paa(t2+t)
     a -> alpha -> the lower energy eigenstate 
     b -> beta  -> the higher energy eigenstate
     Pulse sequence at this impulsive limit is (Eb,Ea,Eb) */
  f0=Pbb(t1+t2)-Pab(t+t1+t2)+Pab(t1)+Pab(t)-Pab(t2)+Paa(t2+t);

  // img part
  Qab_t=Qab(t);

  return (exp(-2.0*f0)*(dfA_2-2.0*dfB_2*cos(2.0*Qab_t)));
}

/* the integrated photon-echo signal at T<0, and tau > |T|
   t1=tau-|T|, t2=|T|
             /infinity
  S(t1,t2) = | dt |R2(t1,t2,t)-R1*(t1,t2,t)|^2
            /0
*/
double S2(double t1, double t2)
{
  gsl_integration_workspace *w;
  double result, error;
  double params[2];
  gsl_function F;

  double ra=0.0,rb=2.0; /* integration range for t */

  w=gsl_integration_workspace_alloc(NWSPACE);

  /* fill in params */
  params[0]=t1;
  params[1]=t2;
  F.function = &f_S2;
  F.params = params;

  gsl_integration_qags(&F, ra, rb, EPSABS/100.0, EPSREL/100.0, NWSPACE, w, &result, &error);
  gsl_integration_workspace_free(w);

  //  printf("S2\n");

  return result;
}

/* kernel for the integrated photon-echo signal at tau>0, T<0, and tau < |T|
   t1=|T|-tau, t2=tau
  the expression using lineshape functions are defined in Eq.2 on Notebook 2.
*/
double f_S3(double t,void *params)
{
  double t1,t2;
  double Qab_t;
  double f0;
  double *p = (double *)params;

  t1=p[0];
  t2=p[1];

  /* integrand is exp(-2*f0)*[dfA-2*dfB*cos(2*Qab(t))],
     where f0=Pbb(t2)+Pab(t1)+Pab(t)-Pab(t2+t)-Pab(t1+t2)+Paa(t1+t2+t)
     a -> alpha -> the lower energy eigenstate 
     b -> beta  -> the higher energy eigenstate
     Pulse sequence at this impulsive limit is (Ea,Eb,Eb) */
  f0=Pbb(t2)+Pab(t1)+Pab(t)-Pab(t2+t)-Pab(t1+t2)+Paa(t1+t2+t);

  // img part
  Qab_t=Qab(t);

  return (exp(-2.0*f0)*(dfA_3-2.0*dfB_3*cos(2.0*Qab_t)));
}

/* the integrated photon-echo signal at tau>0, T<0, and tau < |T|
   t1=|T|-tau, t2=tau
             /infinity
  S(t1,t2) = | dt |R2(t1,t2,t)-R1*(t1,t2,t)|^2
            /0
*/
double S3(double t1, double t2)
{
  gsl_integration_workspace *w;
  double result, error;
  double params[2];
  gsl_function F;

  double ra=0.0,rb=2.0; /* integration range for t */

  w=gsl_integration_workspace_alloc(NWSPACE);

  /* fill in params */
  params[0]=t1;
  params[1]=t2;
  F.function = &f_S3;
  F.params = params;

  gsl_integration_qags(&F, ra, rb, EPSABS/100.0, EPSREL/100.0, NWSPACE, w, &result, &error);
  gsl_integration_workspace_free(w);

  //  printf("S3\n");

  return result;
}

/* Sint determines time-ordering and compute 
   integrated signal using different pathways */
double Sint(double tau, double T)
{
  double Sval;

  if(tau>=0 && T>=0) {
    // region 1, natural ordering (see p. 10 on notebook 2)
    Sval=S1(tau,T);
  } else if (tau>=0 && T<=0 && tau>=fabs(T)) {
    // region 2
    Sval=S2(tau-fabs(T),fabs(T));
  } else if (tau>=0 && T<=0 && tau<fabs(T)) {
    // region 3
    Sval=S3(fabs(T)-tau,tau);
  } else {
    Sval=0.0;
  }

  return Sval;
}

/* kernel for the integrated TG signal; tau=0, t2=T
*/
double f_TG(double t,void *params)
{
  double t2;
  double Qab_t2pt;
  double Qab_t2;
  double f0;
  double GabT,GbbT;
  double *p = (double *)params;
  double dF1,dF2;

  t2=p[0];

  /* conditional probablity, p. 19 on notebook2 */
  if(Kbh>0.0 || Kpb>0.0) {
    GabT=Kbh/(Kbh-Kpb)*(exp(-1.0*Kpb*t2)-exp(-1.0*Kbh*t2));
    GbbT=exp(-1.0*Kbh*t2);
  } else {
    GabT=0.0;
    GbbT=1.0;
  }

  /* integrand is exp(-2*f0)*[dF1-2*dF2*cos(2*(Qab(t2+t)-Qab(t2)))],
     where f0=Paa(t)
     a -> alpha -> the lower energy eigenstate 
     b -> beta  -> the higher energy eigenstate
     Pulse sequence at this impulsive limit is (Eb,Eb,Ea) */
  f0=Paa(t);

  dF1=pow(Me1Me1Me2Me2*(1.0+GabT),2.0)+pow(Me2Me2Me2fMe2f*GbbT,2.0);
  dF2=Me1Me1Me2Me2*(1.0+GabT)*Me2Me2Me2fMe2f*GbbT;
  // modified to make TG-> at t->infinity
  //  dF1=pow(Me1Me1Me2Me2*(exp(-5.0*t2)+GabT),2.0)+pow(Me2Me2Me2fMe2f*GbbT,2.0);
  //  dF2=Me1Me1Me2Me2*(exp(-5.0*t2)+GabT)*Me2Me2Me2fMe2f*GbbT;

  // img part
  Qab_t2=Qab(t2);
  Qab_t2pt=Qab(t2+t);

  //  printf("f0=%20.12f\n",f0);
  return (exp(-2.0*f0)*(dF1-2.0*dF2*cos(2.0*(Qab_t2pt-Qab_t2))));
}

/* the integrated photon-echo signal at tau>0 and T>0

                /infinity
  TG(t1=0,t2) = | dt |R3(t1=0,t2,t)-Gbb(T)*R1*(t1=0,t2,t)+Gab(T)*R2~(t1=0,t2,t)|^2
                /0

		Note that t2 must be in cm unit.
*/
double TG(double t2)
{
  gsl_integration_workspace *w;
  double result, error;
  double params[1];
  gsl_function F;

  double ra=0.0,rb=2.0; /* integration range for t */

  w=gsl_integration_workspace_alloc(NWSPACE);

  /* fill in params */
  params[0]=t2;
  F.function = &f_TG;
  F.params = params;

  gsl_integration_qags(&F, ra, rb, EPSABS/100.0, EPSREL/100.0, NWSPACE, w, &result, &error);
  gsl_integration_workspace_free(w);

  //  printf("S1\n");

  return result;
}

// compute Sint(tau,T) at fixed T for 
// tau from -300fs to 300fs with step size 3fs;
// show result and normalized data
void Sint_fixedT(double T)
{
  double tau,Smax;
  double Stau[4096];
  int i;

  // tau from -300fs to 300fs with step size 3fs
  Smax=0.0;
  for(i=0,tau=-300.0;tau<=300.0;tau=tau+3.0,i++) {
    Stau[i]=Sint(tau/CM2FS,T/CM2FS);
    if(Stau[i]>Smax) Smax=Stau[i];
    printf("tau=%.2f, T=%.2f, Sint=%.18f\n",tau,T,Stau[i]);
  }
  printf("tau=, T=, Sint=\n");

  // also output result nornamized to maximum=1.0
  for(i=0,tau=-300.0;tau<=300.0;tau=tau+3.0,i++) {
    printf("tau=%.2f, T=%.2f, Sint_norm=%.18f\n",tau,T,Stau[i]/Smax);
  }
  printf("tau=, T=, Sint_norm=\n");
}

/* Main program */
int main(int argc, char *argv[])
{
  double T;

  /* Set program name for messages.  */
  program_name = argv[0]; 

  if(argc != 9 && argc != 12 && argc != 15 && argc != 18 && argc != 21 && argc != 22 && argc != 23) {
    usage();
    exit(EXIT_FAILURE);
  }
  gf1.S=0.0;gf1.w=1000.0;gf1.gg=0.0; 
  gf2.S=0.0;gf2.w=1000.0;gf2.gg=0.0; 
  gfj.S=0.0;gfj.w=1000.0;gfj.gg=0.0; 
  beta=1.4387/atof(argv[1]);
  theta=atof(argv[2]);
  gamma1=atof(argv[3]);
  wc1=atof(argv[4]);
  gf1.Delta=atof(argv[5]);
  gamma2=atof(argv[6]);
  wc2=atof(argv[7]);
  gf2.Delta=atof(argv[8]);

  /* no population transfer by default */
  Kbh=0.0;
  Kpb=0.0;

  if(argc >= 12) {
    gammaj=atof(argv[9]);
    wcj=atof(argv[10]);
    gfj.Delta=atof(argv[11]);
  }

  if(argc >= 15) {
    gf1.S=atof(argv[12]);
    gf1.w=atof(argv[13]);
    gf1.gg=atof(argv[14]);
  }

  if(argc >= 18) {
    gf2.S=atof(argv[15]);
    gf2.w=atof(argv[16]);
    gf2.gg=atof(argv[17]);
  }

  if(argc >= 21) {
    gfj.S=atof(argv[18]);
    gfj.w=atof(argv[19]);
    gfj.gg=atof(argv[20]);
  }

  if(argc >= 22) {
    Kbh=CM2FS*atof(argv[21]);
  }
  if(argc >= 23) {
    Kpb=CM2FS*atof(argv[22]);
  }

  printf("\n");
  printf("# beta   = %f (%f K)\n",beta,1.4387/beta);
  printf("# theta  = %f\n",theta);
  printf("# gamma1 = %f\n",gamma1);
  printf("# wc1    = %f\n",wc1);
  printf("# delta1 = %f\n",gf1.Delta);
  printf("# gamma2 = %f\n",gamma2);
  printf("# wc2    = %f\n",wc2);
  printf("# delta2 = %f\n",gf2.Delta);
  printf("# gammaj = %f\n",gammaj);
  printf("# wcj    = %f\n",wcj);
  printf("# deltaj = %f\n",gfj.Delta);
  printf("# S1     = %f\n",gf1.S);
  printf("# w1     = %f\n",gf1.w);
  printf("# gg1    = %f\n",gf1.gg);
  printf("# S2     = %f\n",gf2.S);
  printf("# w2     = %f\n",gf2.w);
  printf("# gg2    = %f\n",gf2.gg);
  printf("# Sj     = %f\n",gfj.S);
  printf("# wj     = %f\n",gfj.w);
  printf("# ggj    = %f\n",gfj.gg);
  printf("# Rbh    = %f\n",Kbh);
  printf("# Rpb    = %f\n",Kpb);
  printf("\n");
 
  // compute transition dipole factors for various pathways
  compute_dipole_factors();
  printf("# Me1Me1Me2Me2   = %f\n",Me1Me1Me2Me2/CSCALE);
  printf("# Me2Me2Me2fMe2f = %f\n",Me2Me2Me2fMe2f/CSCALE);
  printf("# Me1Me2Me1fMe2f = %f\n",Me1Me2Me1fMe2f/CSCALE);
  printf("\n");

  // initialize bath functions
  bath_js03_ct_init(&bf1,beta,gamma1,wc1);
  bath_js03_gt_init(&gf1.gn,&bf1);
  bath_js03_ct_init(&bf2,beta,gamma2,wc2);
  bath_js03_gt_init(&gf2.gn,&bf2);
  bath_js03_ct_init(&bfj,beta,gammaj,wcj);
  bath_js03_gt_init(&gfj.gn,&bfj);

  // initialize cached exciton lineshape functions
  exciton_PQ_init();

  // We turnoff the error handler, which means that the number generated
  // might have errors in them...
  gsl_set_error_handler_off();

  for(T=300.0;T>0.0;T=T-25.0) {
    Sint_fixedT(T);
  }

  for(T=0.0;T>-150.0;T=T-5.0) {
    Sint_fixedT(T);
  }
  for(T=-150.0;T>=-300.0;T=T-10.0) {
    Sint_fixedT(T);
  }

  // TG signal
  for(T=0.0;T<1000.0;T=T+5.0) {
    printf("T=%.2f, TG=%.18f\n",T,TG(T/CM2FS));
  }
  for(T=1000.0;T<5000.0;T=T+50.0) {
    printf("T=%.2f, TG=%.18f\n",T,TG(T/CM2FS));
  }

  gsl_set_error_handler (NULL);

  return 1;
}
