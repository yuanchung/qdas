/***************************************************
 * 2c3ppes-pi.c
 * 
 * By Yuan-Chung Cheng <yccheng@berkeley.edu>
 *
 * Compute full two-color three-pulse photon-echo signal 
 * including different time-ordering conditions
 * (i.e. 750, 800, 750) as a function of coherence 
 * time tau at a given population time T. See p. Eq. 1-3
 * on Notebook 2.
 *
 * Including finite pulse duration integration for the first
 * two pulses.
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
#include "qng2d.h"
#include "bath_js03.h"

#define CM2FS (5309.1)
#define NWSPACE (100000)
#define EPSABS (0.0)
#define EPSREL (1e-3)
#define CSCALE (1e6)
#define LSF_CACHE_MAXPOINTS (30000)
#define LSF_CACHE_MIN (0.0)
#define LSF_CACHE_MAX (5.0)
#define LSF_CACHE_STEPFACTOR (0.0001)

/* pulse duration (intensity FWHM) is 40 fs */
#define PULSE_FWHM (40.0/CM2FS)

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
double beta, gamma1, wc1, gamma2, wc2;

/* Global variables that defines the system */
/* two bath functions */
js03_bath_func bf1,bf2;
lineshape_func gf1,gf2;

/* transition dipole moments |a>->|B>, |b>->|H> */
/* BL and HL from GD paper */
const double mu_a[3]={0.78550,  0.44743,  0.42755};
const double mu_b[3]={0.18975*0.6,  -0.96874*0.6,  -0.15979*0.6};
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
Usage: %s temp theta gamma1 wc1 Delta1 gamma2 wc2 Delta2 [S1 w1 gg1 S2 w2 gg2 Rbh Rpb]\n\
\n\
       temp is the temperature (in K).\n\
       theta is the mixing angle between the two chromophores.\n\
       gamma1 and wc1 define bath for the lower-energy chromophore.\n\
       gamma2 and wc2 define bath for the higher-energy one.\n\
       Delta1 and Delta2 are standard dev. of Gaussian static disorder\n\
       for chromophore 1 and 2, respectively.\n\
       Optional argvs:\n\
       Sj and wj are vibrational Huang-Rhys factor and mode\n\
       frequency for chromophore 1 and 2, respectively.\n\
       ggj are the decay rate (in cm^-1) for the vib. modes.\n\
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
}

/* lineshape functions for localized states,
   including inhomogeneous and vib. contributions */
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
    yy[i]=c4*gt_r(&gf1,xx[i])+s4*gt_r(&gf2,xx[i]);
  }
  Paa_Spline_Acc=gsl_interp_accel_alloc();
  Paa_Spline=gsl_spline_alloc (gsl_interp_cspline,npoints);
  gsl_spline_init(Paa_Spline,xx,yy,npoints);
  /* Qaa(t) */
  printf("  Initializing Qaa(t)...\n");
  for(i=0;i<npoints;i++) {
    yy[i]=c4*gt_i(&gf1,xx[i])+s4*gt_i(&gf2,xx[i]);
  }
  Qaa_Spline_Acc=gsl_interp_accel_alloc();
  Qaa_Spline=gsl_spline_alloc (gsl_interp_cspline,npoints);
  gsl_spline_init(Qaa_Spline,xx,yy,npoints);
  /* Pbb(t) */
  printf("  Initializing Pbb(t)...\n");
  for(i=0;i<npoints;i++) {
    yy[i]=s4*gt_r(&gf1,xx[i])+c4*gt_r(&gf2,xx[i]);
  }
  Pbb_Spline_Acc=gsl_interp_accel_alloc();
  Pbb_Spline=gsl_spline_alloc (gsl_interp_cspline,npoints);
  gsl_spline_init(Pbb_Spline,xx,yy,npoints);
  /* Qbb(t) */
  printf("  Initializing Qbb(t)...\n");
  for(i=0;i<npoints;i++) {
    yy[i]=s4*gt_i(&gf1,xx[i])+c4*gt_i(&gf2,xx[i]);
  }
  Qbb_Spline_Acc=gsl_interp_accel_alloc();
  Qbb_Spline=gsl_spline_alloc (gsl_interp_cspline,npoints);
  gsl_spline_init(Qbb_Spline,xx,yy,npoints);
  /* Pab(t) */
  printf("  Initializing Pab(t)...\n");
  for(i=0;i<npoints;i++) {
    yy[i]=s2c2*(gt_r(&gf1,xx[i])+gt_r(&gf2,xx[i]));
  }
  Pab_Spline_Acc=gsl_interp_accel_alloc();
  Pab_Spline=gsl_spline_alloc (gsl_interp_cspline,npoints);
  gsl_spline_init(Pab_Spline,xx,yy,npoints);
  /* Qab(t) */
  printf("  Initializing Qab(t)...\n");
  for(i=0;i<npoints;i++) {
    yy[i]=s2c2*(gt_i(&gf1,xx[i])+gt_i(&gf2,xx[i]));
  }
  Qab_Spline_Acc=gsl_interp_accel_alloc();
  Qab_Spline=gsl_spline_alloc (gsl_interp_cspline,npoints);
  gsl_spline_init(Qab_Spline,xx,yy,npoints);

  printf("\n");

}

/* Impulsive response function at region 1 -- t1>0 and t2>0
   RR1 = R3 - R1*
   the third pulse is centered at t=0
   Pulse sequence at this impulsive limit is (Eb,Eb,Ea)
   the expression used are defined in Eq.5 on p.28-p.30, Notebook 2.
*/
double RR1_r(double t1,double t2,double t)
{
  double f0,a1,a2;
  double tmp1,tmp2;

  f0=Paa(t)-Pab(t+t1+t2)+Pab(t1+t2)+Pab(t+t2)-Pab(t2)+Pbb(t1);

  // exp(I*..) factor
  // a1 = -Qaa(t)-Qab(t+t1+t2)+Qab(t1+t2)+Qab(t+t2)-Qab(t2)+Qbb(t1)
  // a2 = -Qaa(t)-Qab(t+t1+t2)+Qab(t1+t2)-Qab(t+t2)+Qab(t2)+Qbb(t1)
  tmp1 = -1.0*Qaa(t)-Qab(t+t1+t2)+Qab(t1+t2)+Qbb(t1);
  tmp2 = Qab(t+t2)-Qab(t2);
  a1=tmp1+tmp2;
  a2=tmp1-tmp2;

  //  printf("f0=%20.12f\n",f0);
  return (exp(-1.0*f0)*(Me1Me1Me2Me2*cos(a1)-Me2Me2Me2fMe2f*cos(a2)));
}

double RR1_i(double t1,double t2,double t)
{
  double f0,a1,a2;
  double tmp1,tmp2;

  f0=Paa(t)-Pab(t+t1+t2)+Pab(t1+t2)+Pab(t2+t)-Pab(t2)+Pbb(t1);

  // exp(I*..) factor
  // a1 = -Qaa(t)-Qab(t+t1+t2)+Qab(t1+t2)+Qab(t+t2)-Qab(t2)+Qbb(t1)
  // a2 = -Qaa(t)-Qab(t+t1+t2)+Qab(t1+t2)-Qab(t+t2)+Qab(t2)+Qbb(t1)
  tmp1 = -1.0*Qaa(t)-Qab(t+t1+t2)+Qab(t1+t2)+Qbb(t1);
  tmp2 = Qab(t+t2)-Qab(t2);
  a1=tmp1+tmp2;
  a2=tmp1-tmp2;

  //  printf("f0=%20.12f\n",f0);
  return (exp(-1.0*f0)*(Me1Me1Me2Me2*sin(a1)-Me2Me2Me2fMe2f*sin(a2)));
}

/* Impulsive response function at region 2 -- t1>0, t2<0, and t1 > |t2|
   RR2 = R2 - R1*
   the third pulse is centered at t=0
   Pulse sequence at this impulsive limit is (Eb,Ea,Eb)
   the expression used are defined in Eq.5 on p.28-p.30, Notebook 2.
*/
double RR2_r(double t1,double t2,double t)
{
  double f0,a1,a2;
  double tmp1,tmp2;

  f0=Paa(t+t2)-Pab(t+t1+t2)+Pab(t1)+Pab(t)-Pab(t2)+Pbb(t1+t2);

  // exp(I*..) factor
  // a1 = -Qaa(t+t2)-Qab(t+t1+t2)+Qab(t1)+Qab(t)+Qab(t2)+Qbb(t1+t2)
  // a2 = -Qaa(t+t2)-Qab(t+t1+t2)+Qab(t1)-Qab(t)+Qab(t2)+Qbb(t1+t2)
  tmp1 = -1.0*Qaa(t+t2)-Qab(t+t1+t2)+Qab(t1)+Qab(t2)+Qbb(t1+t2);
  tmp2 = Qab(t);
  a1=tmp1+tmp2;
  a2=tmp1-tmp2;

  //  printf("f0=%20.12f\n",f0);
  return (exp(-1.0*f0)*(Me1Me1Me2Me2*cos(a1)-Me1Me2Me1fMe2f*cos(a2)));
}

double RR2_i(double t1,double t2,double t)
{
  double f0,a1,a2;
  double tmp1,tmp2;

  f0=Paa(t+t2)-Pab(t+t1+t2)+Pab(t1)+Pab(t)-Pab(t2)+Pbb(t1+t2);

  // exp(I*..) factor
  // a1 = -Qaa(t+t2)-Qab(t+t1+t2)+Qab(t1)+Qab(t)+Qab(t2)+Qbb(t1+t2)
  // a2 = -Qaa(t+t2)-Qab(t+t1+t2)+Qab(t1)-Qab(t)+Qab(t2)+Qbb(t1+t2)
  tmp1 = -1.0*Qaa(t+t2)-Qab(t+t1+t2)+Qab(t1)+Qab(t2)+Qbb(t1+t2);
  tmp2 = Qab(t);
  a1=tmp1+tmp2;
  a2=tmp1-tmp2;

  //  printf("f0=%20.12f\n",f0);
  return (exp(-1.0*f0)*(Me1Me1Me2Me2*sin(a1)-Me1Me2Me1fMe2f*sin(a2)));
}


/* Impulsive response function at region 3 -- t1>0, t2<0, and t1 < |t2|
   RR3 = R1 - R2*
   the third pulse is centered at t=0
   Pulse sequence at this impulsive limit is (Ea,Eb,Eb)
   the expression used are defined in Eq.5 on p.28-p.30, Notebook 2.
*/
double RR3_r(double t1,double t2,double t)
{
  double f0,a1,a2;
  double tmp1,tmp2;

  f0=Paa(t+t1+t2)-Pab(t+t2)+Pab(t1)+Pab(t)-Pab(t1+t2)+Pbb(t2);

  // exp(I*..) factor
  // a1 = -Qaa(t+t1+t2)-Qab(t+t2)-Qab(t1)+Qab(t)+Qab(t1+t2)+Qbb(t2)
  // a2 = -Qaa(t+t1+t2)-Qab(t+t2)-Qab(t1)-Qab(t)+Qab(t1+t2)+Qbb(t2)
  tmp1 = -1.0*Qaa(t+t1+t2)-Qab(t+t2)-Qab(t1)+Qab(t1+t2)+Qbb(t2);
  tmp2 = Qab(t);
  a1=tmp1+tmp2;
  a2=tmp1-tmp2;

  //  printf("f0=%20.12f\n",f0);
  return (exp(-1.0*f0)*(Me1Me1Me2Me2*cos(a1)-Me1Me2Me1fMe2f*cos(a2)));
}

double RR3_i(double t1,double t2,double t)
{
  double f0,a1,a2;
  double tmp1,tmp2;

  f0=Paa(t+t1+t2)-Pab(t+t2)+Pab(t1)+Pab(t)-Pab(t1+t2)+Pbb(t2);

  // exp(I*..) factor
  // a1 = -Qaa(t+t1+t2)-Qab(t+t2)-Qab(t1)+Qab(t)+Qab(t1+t2)+Qbb(t2)
  // a2 = -Qaa(t+t1+t2)-Qab(t+t2)-Qab(t1)-Qab(t)+Qab(t1+t2)+Qbb(t2)
  tmp1 = -1.0*Qaa(t+t1+t2)-Qab(t+t2)-Qab(t1)+Qab(t1+t2)+Qbb(t2);
  tmp2 = Qab(t);
  a1=tmp1+tmp2;
  a2=tmp1-tmp2;

  //  printf("f0=%20.12f\n",f0);
  return (exp(-1.0*f0)*(Me1Me1Me2Me2*sin(a1)-Me1Me2Me1fMe2f*sin(a2)));
}

/* Impulsive response function at region 6 -- t1<0, t2>0, and t1 < |t2|
   RR6 = R4 - R2*
   the third pulse is centered at t=0
   Pulse sequence at this impulsive limit is (Eb,Eb,Ea)
   the expression used are defined in Eq.5 on p.28-p.30, Notebook 2.
*/
double RR6_r(double t1,double t2,double t)
{
  double f0,a1,a2;
  double tmp1,tmp2;

  f0=Paa(t)-Pab(t+t2)+Pab(t+t1+t2)+Pab(t2)-Pab(t1+t2)+Pbb(t1);

  // exp(I*..) factor
  // a1 = -Qaa(t)+Qab(t+t2)-Qab(t+t1+t2)-Qab(t2)+Qab(t1+t2)-Qbb(t1)
  // a2 = -Qaa(t)-Qab(t+t2)-Qab(t+t1+t2)+Qab(t2)+Qab(t1+t2)-Qbb(t1)
  tmp1 = -1.0*Qaa(t)-Qab(t+t1+t2)+Qab(t1+t2)-Qbb(t1);
  tmp2 = Qab(t+t2)-Qab(t2);
  a1=tmp1+tmp2;
  a2=tmp1-tmp2;

  //  printf("f0=%20.12f\n",f0);
  return (exp(-1.0*f0)*(Me1Me1Me2Me2*cos(a1)-Me1Me2Me1fMe2f*cos(a2)));
}

double RR6_i(double t1,double t2,double t)
{
  double f0,a1,a2;
  double tmp1,tmp2;

  f0=Paa(t)-Pab(t+t2)+Pab(t+t1+t2)+Pab(t2)-Pab(t1+t2)+Pbb(t1);

  // exp(I*..) factor
  // a1 = -Qaa(t)+Qab(t+t2)-Qab(t+t1+t2)-Qab(t2)+Qab(t1+t2)-Qbb(t1)
  // a2 = -Qaa(t)-Qab(t+t2)-Qab(t+t1+t2)+Qab(t2)+Qab(t1+t2)-Qbb(t1)
  tmp1 = -1.0*Qaa(t)-Qab(t+t1+t2)+Qab(t1+t2)-Qbb(t1);
  tmp2 = Qab(t+t2)-Qab(t2);
  a1=tmp1+tmp2;
  a2=tmp1-tmp2;

  //  printf("f0=%20.12f\n",f0);
  return (exp(-1.0*f0)*(Me1Me1Me2Me2*sin(a1)-Me1Me2Me1fMe2f*sin(a2)));
}

/* The following is the real 2D functions to be integrated */

/* Global parameters tau, T and t --> center of pulses;
   these parameters must be assigned before calling f2d_*
   Note that the center of the third pulse is set to t=0;
   these parameters denote the measuring time t (t>0);
   center of pulse -k1 is param_center1 (param_center1<0),
   center of pulse k2 is param_center2 (param_center2<0);
   order of param_center1 and param_center2 can be exchanged
*/
double param_center1;
double param_center2;
double param_t;

/* for region 1 and region  6; Eb-Eb-Ea */
double f2d_P16_r(double tau1, double tau2)
{
  double E1E2;

  /* Gaussian pulse field factor */
  E1E2=exp(-1.385*(tau1-param_center1)*(tau1-param_center1)/PULSE_FWHM/PULSE_FWHM)
    * exp(-1.385*(tau2-param_center2)*(tau2-param_center2)/PULSE_FWHM/PULSE_FWHM);

  if(tau1<=tau2)
    return E1E2*RR1_r(tau2-tau1,-1.0*tau2,param_t);
  else
    return E1E2*RR6_r(tau1-tau2,-1.0*tau1,param_t);
}
double f2d_P16_i(double tau1, double tau2)
{
  double E1E2;
  /* Gaussian pulse field factor */
  E1E2=exp(-1.385*(tau1-param_center1)*(tau1-param_center1)/PULSE_FWHM/PULSE_FWHM)
    * exp(-1.385*(tau2-param_center2)*(tau2-param_center2)/PULSE_FWHM/PULSE_FWHM);

  if(tau1<=tau2)
    return E1E2*RR1_i(tau2-tau1,-1.0*tau2,param_t);
  else
    return E1E2*RR6_i(tau1-tau2,-1.0*tau1,param_t);
}
/* this function sets up param_* and call gng2d to 
   obtain integrated polarization P(tau,T,t);
   |P(tau,T,t)|^2 is returned.
*/
double P16_square(double t, void *params)
{
  double lower1,lower2,upper1,upper2;
  double P_r,P_i,error;
  double tau,T;
  double *p = (double *)params;

  tau=p[0];
  T=p[1];

  /* we assume T>0 and tau>-T without any check */
  param_center1=-1.0*T-tau;
  param_center2=-1.0*T;
  param_t=t;

  /* integration area; only use (+-) 2*FWHM pulse */
  lower1=param_center1-2.0*PULSE_FWHM;
  lower2=param_center2-2.0*PULSE_FWHM;
  /* upper bound is more tricky; we make sure only tau[12] < 0 is integrated */
  upper1=param_center1+2.0*PULSE_FWHM;
  if(upper1>0.0) upper1=0.0;
  upper2=param_center2+2.0*PULSE_FWHM;
  if(upper2>0.0) upper2=0.0;

  /* call qng2d to perform 2D numerical integration */
  //  printf("HERE1\n");
  qng2d(&f2d_P16_r, lower1, upper1, lower2, upper2, EPSABS, EPSREL/10.0, &P_r, &error);
  qng2d(&f2d_P16_i, lower1, upper1, lower2, upper2, EPSABS, EPSREL/10.0, &P_i, &error);

  /* return |P(t)|^2 */
  return (P_r*P_r+P_i*P_i);
}

/* for region 2 and region  3; Eb-Ea-Eb
   Note that in this case
   center of pulse -k1 is param_center1 (param_center1<0),
   center of pulse k3 is param_center2 (param_center2<0);
   the k3 pulse is Ea. these sequence is defined in P23_square()
 */
double f2d_P23_r(double tau1, double tau2)
{
  double E1E2;
  /* Gaussian pulse field factor */
  E1E2=exp(-1.385*(tau1-param_center1)*(tau1-param_center1)/PULSE_FWHM/PULSE_FWHM)
    * exp(-1.385*(tau2-param_center2)*(tau2-param_center2)/PULSE_FWHM/PULSE_FWHM);

  if(tau1<=tau2)
    return E1E2*RR2_r(tau2-tau1,-1.0*tau2,param_t);
  else
    return E1E2*RR3_r(tau1-tau2,-1.0*tau1,param_t);
}
double f2d_P23_i(double tau1, double tau2)
{
  double E1E2;
  /* Gaussian pulse field factor */
  E1E2=exp(-1.385*(tau1-param_center1)*(tau1-param_center1)/PULSE_FWHM/PULSE_FWHM)
    * exp(-1.385*(tau2-param_center2)*(tau2-param_center2)/PULSE_FWHM/PULSE_FWHM);

  if(tau1<=tau2)
    return E1E2*RR2_i(tau2-tau1,-1.0*tau2,param_t);
  else
    return E1E2*RR3_i(tau1-tau2,-1.0*tau1,param_t);
}
/* this function sets up param_* and call gng2d to 
   obtain integrated polarization P(tau,T,t);
   |P(tau,T,t)|^2 is returned.
*/
double P23_square(double t, void *params)
{
  double lower1,lower2,upper1,upper2;
  double P_r,P_i,error;
  double tau,T;
  double *p = (double *)params;

  tau=p[0];
  T=p[1];

  /* we assume T<0 and tau>0 without any check */
  param_center1=-1.0*tau;
  param_center2=T;
  param_t=t;

  /* integration area; only use (+-) 2*FWHM pulse */
  lower1=param_center1-2.0*PULSE_FWHM;
  lower2=param_center2-2.0*PULSE_FWHM;
  /* upper bound is more tricky; we make sure only tau[12] < 0 is integrated */
  upper1=param_center1+2.0*PULSE_FWHM;
  if(upper1>0.0) upper1=0.0;
  upper2=param_center2+2.0*PULSE_FWHM;
  if(upper2>0.0) upper2=0.0;

  /* call qng2d to perform 2D numerical integration */
  qng2d(f2d_P23_r, lower1, upper1, lower2, upper2, EPSABS, EPSREL/10.0, &P_r, &error);
  qng2d(f2d_P23_i, lower1, upper1, lower2, upper2, EPSABS, EPSREL/10.0, &P_i, &error);

  /* return |P(t)|^2 */
  return (P_r*P_r+P_i*P_i);
}


/* the integrated photon-echo signal at (tau,T)

             /infinity
  S(tau,T) = | dt |P(tau,T,t)|^2
            /0
*/
double Sint(double tau, double T)
{
  //  gsl_integration_workspace *w;
  double result, error;
  size_t neval;
  double params[2];
  gsl_function F;

  double ra=0.0,rb=2.0; /* integration range for t; this is in cm unit so it is > 10000 fs */

  /* fill in params */
  params[0]=tau;
  params[1]=T;
  F.params = params;
  /* here we determine the region of (tau,T), so that we can call correct P(t) */
  if(T>=0.0 && tau>=(-1.0*T)) {
    // region 1 or 6
    //    printf("CASE 1: tau=%f, T=%f\n",tau,T);
    F.function = &P16_square;
  } else if(T<0.0 && tau >=0.0) {
    // region 2 or 3
    //    printf("CASE 2: tau=%f, T=%f\n",tau,T);
    F.function = &P23_square;
  } else {
    // we ignore anything in region 4 and 5
    //    printf("CASE 3: tau=%f, T=%f\n",tau,T);
    return 0.0;
  }

  //  w=gsl_integration_workspace_alloc(NWSPACE);
  //  gsl_integration_qags(&F, ra, rb, EPSABS, EPSREL, NWSPACE, w, &result, &error);
  //  gsl_integration_workspace_free(w);
  gsl_integration_qng (&F, ra, rb, EPSABS, EPSREL, &result, &error, &neval);

  return result;
}

/* kernel for the integrated TG signal; tau=0, t2=T
*/
double f_TG(double t,void *params)
{
  double t2;
  double Qab_t2pt;
  double Qab_t2;
  double f0;
  double s4,c4,s2c2;
  double GabT,GbbT;
  double *p = (double *)params;
  double dF1,dF2;

  t2=p[0];

  /* some useful constant factors */
  s4=sin(theta)*sin(theta)*sin(theta)*sin(theta);
  c4=cos(theta)*cos(theta)*cos(theta)*cos(theta);
  s2c2=sin(theta)*sin(theta)*cos(theta)*cos(theta);

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
  f0=c4*gt_r(&gf1,t)+s4*gt_r(&gf2,t);

  dF1=pow(Me1Me1Me2Me2*(1.0+GabT),2.0)+pow(Me2Me2Me2fMe2f*GbbT,2.0);
  dF2=Me1Me1Me2Me2*(1.0+GabT)*Me2Me2Me2fMe2f*GbbT;
  // modified to make TG-> at t->infinity
  //  dF1=pow(Me1Me1Me2Me2*(exp(-5.0*t2)+GabT),2.0)+pow(Me2Me2Me2fMe2f*GbbT,2.0);
  //  dF2=Me1Me1Me2Me2*(exp(-5.0*t2)+GabT)*Me2Me2Me2fMe2f*GbbT;

  // img part
  Qab_t2=s2c2*(gt_i(&gf1,t2)+gt_i(&gf2,t2));
  Qab_t2pt=s2c2*(gt_i(&gf1,t2+t)+gt_i(&gf2,t2+t));

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

  double ra=0.0,rb=90.0; /* integration range for t */

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
// tau from -300fs to 300fs with step size 5fs;
// show result and normalized data
void Sint_fixedT(double T)
{
  double tau,Smax;
  double Stau[4096];
  int i;

  // tau from -300fs to 300fs with step size 3fs
  Smax=0.0;
  for(i=0,tau=-300.0;tau<=300.0;tau=tau+5.0,i++) {
    Stau[i]=Sint(tau/CM2FS,T/CM2FS);
    if(Stau[i]>Smax) Smax=Stau[i];
    printf("tau=%.2f, T=%.2f, Sint=%.18f\n",tau,T,Stau[i]);
  }
  printf("tau=, T=, Sint=\n");

  // also output result nornamized to maximum=1.0
  for(i=0,tau=-300.0;tau<=300.0;tau=tau+5.0,i++) {
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

  if(argc != 9 && argc != 12 && argc != 15 && argc != 16 && argc != 17) {
    usage();
    exit(EXIT_FAILURE);
  }
  gf1.S=0.0;gf1.w=1000.0;gf1.gg=0.0; 
  gf2.S=0.0;gf2.w=1000.0;gf2.gg=0.0; 
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
    gf1.S=atof(argv[9]);
    gf1.w=atof(argv[10]);
    gf1.gg=atof(argv[11]);
  }

  if(argc >= 15) {
    gf2.S=atof(argv[12]);
    gf2.w=atof(argv[13]);
    gf2.gg=atof(argv[14]);
  }

  if(argc >= 16) {
    Kbh=CM2FS*atof(argv[15]);
  }
  if(argc >= 17) {
    Kpb=CM2FS*atof(argv[16]);
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
  printf("# S1     = %f\n",gf1.S);
  printf("# w1     = %f\n",gf1.w);
  printf("# gg1    = %f\n",gf1.gg);
  printf("# S2     = %f\n",gf2.S);
  printf("# w2     = %f\n",gf2.w);
  printf("# gg2    = %f\n",gf2.gg);
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

  // initialize cached exciton lineshape functions
  exciton_PQ_init();

  // We turnoff the error handler, which means that the number generated
  // might have errors in them...
  gsl_set_error_handler_off();

  for(T=0.0;T>=-300.0;T=T-10.0) {
    Sint_fixedT(T);
  }
  
  gsl_set_error_handler (NULL);

  // TG signal
  //  for(T=0.0;T<1000.0;T=T+5.0) {
  //    printf("T=%.2f, TG=%.18f\n",T,TG(T/CM2FS));
  //  }
  //  for(T=1000.0;T<5000.0;T=T+50.0) {
  //    printf("T=%.2f, TG=%.18f\n",T,TG(T/CM2FS));
  //  }

  return 1;
}

/*
 * $Log$
 * Revision 1.2  2007/01/12 07:19:59  platin
 *   - major changes in 2c3ppes* files.
 *   - add a 2c3ppes program that uses HTGAU bath.
 *
 *
 */
