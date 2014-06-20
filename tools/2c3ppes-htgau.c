/***************************************************
 * 2c3ppes-htgau.c
 * 
 * By Yuan-Chung Cheng <yccheng@berkeley.edu>
 *
 * Compute two-color three-pulse photon-echo signal 
 * including different time-ordering conditions
 * (i.e. 750, 800, 750) as a function of time t1 and t2;
 * note that the ordering here is 750-800-750
 *
 * Using the HTGau bath module.
 *
 * Notations:
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
#include <gsl/gsl_min.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "aux.h"
#include "bath_htgau.h"
#include "spectrum.h"

#define CM2FS (5309.1)
#define NWSPACE (100000)
#define EPSABS (1e-6)
#define EPSREL (1e-6)
#define CSCALE (1e3)
#define LSF_CACHE_MAXPOINTS (30000)
#define LSF_CACHE_MIN (0.0)
#define LSF_CACHE_MAX (5.0)
#define LSF_CACHE_STEPFACTOR (0.0001)

/* global variables */
const char *program_name;

/* Global variable given in the command line */
double theta; /* mixing angle between two chromophores */
/* define the bath */
double beta, tau0, lambda1, lambda2, corr_coeff;
/* static disorder */
double Delta1,Delta2;
/* vib. modes */
double S1,w1,gg1,phi1,S2,w2,gg2,phi2,Sj,wj,ggj,phij;

/* Global variables that defines the system */
/* transition dipole moments |a>->|B>, |b>->|H> */

/* BL and HL from GD paper */
const double mu_a[3]={0.78550,  0.44743,  0.42755};
//const double mu_b[3]={0.18975*0.95,  -0.96874*0.95,  -0.15979*0.95};
const double mu_b[3]={0.18975*0.8,  -0.96874*0.8,  -0.15979*0.8};

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
Usage: %s temp theta tau0 lambda1 lambda2 c Delta1 Delta2 [S1 w1 gg1 phi1 S2 w2 gg2 phi2 Sj wj ggj phij Rbh Rpb]\n\
\n\
       temp is the temperature (in K).\n\
       theta is the mixing angle between the two chromophores.\n\
       tau0 is the decay time of the Gaussian correlation function (in fs):\n\
            M(t) = lambda*exp(-t^2/tau0^2)\n\
       lambda1 define bath for the lower-energy chromophore (in cm^-1).\n\
       lambda2 define bath for the higher-energy one (in cm^-1).\n\
       c is the cross-correlation coefficient for the two bath (c between -1 .. 1).\n\
       Delta[12] define static disorder (not correlated; in cm^-1).\n\
       Optional argvs:\n\
       S[12] and w[12] are vibrational Huang-Rhys factor and mode\n\
       frequency in cm^-1 for chromophore 1 and 2, respectively.\n\
       gg[12] are the decay rate (in cm^-1) for the vib. modes.\n\
       phi[12] are the initial phase (in rad) for the vib. modes.\n\
       Sj, wj, ggj, phij are for off-diagonal vib. couplings.\n\
       Note: these vib. modes are not correlated (via c).\n\
       Rbh is |H> -> |B> population transfer rate in fs^-1.\n\
       Rpb is |B> -> |P> population transfer rate in fs^-1.\n\
\n",
          program_name);
}

/* functions used to find the maximum point of (xx[],yy[]) */
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

void find_max(double *xmax, double *ymax, 
	      const double *xx, const double *yy, const size_t npoints, const double xinit)
{

  int status;
  int iter=0,max_iter=100;
  const gsl_min_fminimizer_type *T;
  gsl_min_fminimizer *s;
  gsl_function F;

  double xlower,xupper;

  /* initialize cspline function */
  cspline_max_acc=gsl_interp_accel_alloc();
  cspline_max_spline=gsl_spline_alloc(gsl_interp_cspline, npoints);
  gsl_spline_init(cspline_max_spline, xx, yy, npoints);

  /* now ready to initialize and invoke the GSL minimization process */
  F.function=&func_y;
  F.params=NULL;

  T=gsl_min_fminimizer_brent;
  s=gsl_min_fminimizer_alloc (T); 
  gsl_min_fminimizer_set (s, &F, xinit, xx[0], xx[npoints-1]);
  
  do {
    iter++;
    status=gsl_min_fminimizer_iterate(s);
    
    *xmax=gsl_min_fminimizer_x_minimum(s);
    xlower=gsl_min_fminimizer_x_lower(s);
    xupper=gsl_min_fminimizer_x_upper(s);
    
    status=gsl_min_test_interval(xlower,xupper,0.0,0.01);
    
  } while (status == GSL_CONTINUE && iter < max_iter);

  // eval the maximum value
  *ymax=gsl_spline_eval(cspline_max_spline,*xmax,cspline_max_acc);
  
  /* done; clean up */
  gsl_spline_free(cspline_max_spline);
  gsl_interp_accel_free(cspline_max_acc);
  gsl_min_fminimizer_free(s);
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

/* lineshape function, vib. contributions */
double gv_r(double S, double w, double g,double phi,double t)
{
  double ret,tmp;

  // ret=ret+S*(coth(beta*w/2.0)*(1.0-cos(w*t)));
  // plus an ad-hoc dissipative oscillator
  tmp=exp(-1.0*g*t);
  ret=-1.0*S*w*w*coth(beta*w/2.0)*pow(g*g+w*w,-2.0)*(sin(phi)*t*w*w*w-g*cos(phi)*t*w*w-w*w*cos(phi)+w*w*tmp*cos(w*t+phi)-2.0*g*w*sin(phi)+sin(phi)*t*w*g*g+2.0*g*tmp*w*sin(w*t+phi)-g*g*g*cos(phi)*t-g*g*tmp*cos(w*t+phi)+g*g*cos(phi));
  return ret;
}

double gv_i(double S, double w, double g,double phi,double t)
{
  double ret,tmp;

  // ret=ret+S*(sin(w*t)-w*t);
  // ad-hoc dissipative oscillator
  tmp=exp(-1.0*g*t);
  ret=-1.0*S*w*w*pow(g*g+w*w,-2.0)*(-2.0*g*w*cos(phi)+w*w*sin(phi)-g*g*sin(phi)+w*cos(phi)*g*g*t+w*w*w*cos(phi)*t+g*g*g*sin(phi)*t+g*sin(phi)*t*w*w+2*w*tmp*g*cos(w*t+phi)-w*w*tmp*sin(w*t+phi)+g*g*tmp*sin(w*t+phi));

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

  double g_11,g_22,g_12,g_j;

  double xx[LSF_CACHE_MAXPOINTS];
  double aa[LSF_CACHE_MAXPOINTS];
  double bb[LSF_CACHE_MAXPOINTS];
  double ab[LSF_CACHE_MAXPOINTS];
  double t;
  size_t i,npoints;

  printf("Initializing cached lineshape functions for exciton states:\n");

  /* some useful constant factors */
  s4=sin(theta)*sin(theta)*sin(theta)*sin(theta);
  c4=cos(theta)*cos(theta)*cos(theta)*cos(theta);
  s2c2=sin(theta)*sin(theta)*cos(theta)*cos(theta);

  /* set up grid points to cache */
  i=0;t=LSF_CACHE_MIN;
  while(t<=LSF_CACHE_MAX && i<LSF_CACHE_MAXPOINTS) {
      xx[i]=t;
      i++;
      if(t<=0.2) {
	/* 1/LSF_CACHE_STEPFACTOR points between 0.0-0.2; this is about 1000 fs */
        t = t + LSF_CACHE_STEPFACTOR*0.2;
      } else {
	/* exponential step size after t>0.2 */
        t = t*1.05;
      }
  }
  npoints=i;

  /* All real parts: Paa(t), Pbb(t), Pab(t) */
  printf("  Initializing Paa(t)/Pbb(t)/Pab(t)...\n");
  for(i=0;i<npoints;i++) {
    t=xx[i];
    g_11=bath_htgau_gn_r(t,lambda1,tau0,beta)+gv_r(S1,w1,gg1,phi1,t)+0.5*Delta1*Delta1*t*t;
    g_22=bath_htgau_gn_r(t,lambda2,tau0,beta)+gv_r(S2,w2,gg2,phi2,t)+0.5*Delta2*Delta2*t*t;
    /* only vib. contribution to g_j */
    g_j=gv_r(Sj,wj,ggj,phij,t);
    /* cross-correlation term, note that coherent vib. modes 
       are not correlated */
    g_12=bath_htgau_gn_r(t,corr_coeff*sqrt(lambda1*lambda2),tau0,beta)+0.5*corr_coeff*Delta1*Delta2*t*t;
    aa[i]=c4*g_11+s4*g_22+4.0*s2c2*g_j+2.0*s2c2*g_12;
    bb[i]=s4*g_11+c4*g_22+4.0*s2c2*g_j+2.0*s2c2*g_12;
    ab[i]=s2c2*(g_11+g_22)-4.0*s2c2*g_j+(s4+c4)*g_12;
  }
  Paa_Spline_Acc=gsl_interp_accel_alloc();
  Paa_Spline=gsl_spline_alloc (gsl_interp_cspline,npoints);
  gsl_spline_init(Paa_Spline,xx,aa,npoints);
  Pbb_Spline_Acc=gsl_interp_accel_alloc();
  Pbb_Spline=gsl_spline_alloc (gsl_interp_cspline,npoints);
  gsl_spline_init(Pbb_Spline,xx,bb,npoints);
  Pab_Spline_Acc=gsl_interp_accel_alloc();
  Pab_Spline=gsl_spline_alloc (gsl_interp_cspline,npoints);
  gsl_spline_init(Pab_Spline,xx,ab,npoints);

  /* Qaa(t)/Qbb(t)/Qab(t) */
  printf("  Initializing Qaa(t)/Qbb(t)/Qab(t)...\n");
  for(i=0;i<npoints;i++) {
    t=xx[i];
    g_11=bath_htgau_gn_i(t,lambda1,tau0)+gv_i(S1,w1,gg1,phi1,t);
    g_22=bath_htgau_gn_i(t,lambda2,tau0)+gv_i(S2,w2,gg2,phi2,t);
    g_j=gv_i(Sj,wj,ggj,phij,t);
    g_12=bath_htgau_gn_i(t,corr_coeff*sqrt(lambda1*lambda2),tau0);
    aa[i]=c4*g_11+s4*g_22+4.0*s2c2*g_j+2.0*s2c2*g_12;
    bb[i]=s4*g_11+c4*g_22+4.0*s2c2*g_j+2.0*s2c2*g_12;
    ab[i]=s2c2*(g_11+g_22)-4.0*s2c2*g_j+(s4+c4)*g_12;
  }
  Qaa_Spline_Acc=gsl_interp_accel_alloc();
  Qaa_Spline=gsl_spline_alloc (gsl_interp_cspline,npoints);
  gsl_spline_init(Qaa_Spline,xx,aa,npoints);
  Qbb_Spline_Acc=gsl_interp_accel_alloc();
  Qbb_Spline=gsl_spline_alloc (gsl_interp_cspline,npoints);
  gsl_spline_init(Qbb_Spline,xx,bb,npoints);
  Qab_Spline_Acc=gsl_interp_accel_alloc();
  Qab_Spline=gsl_spline_alloc (gsl_interp_cspline,npoints);
  gsl_spline_init(Qab_Spline,xx,ab,npoints);

  printf("\n");

}

/* kernel for the integrated photon-echo signal at region 1; tau>0 and T>0
   t1=tau, t2=T
  the expression using lineshape functions are defined in Eq.3 on Notebook 2.
*/
double f_S_reg1(double t,void *params)
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
double S_reg1(double t1, double t2)
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
  F.function = &f_S_reg1;
  F.params = params;

  gsl_integration_qags(&F, ra, rb, EPSABS/100.0, EPSREL/100.0, NWSPACE, w, &result, &error);
  gsl_integration_workspace_free(w);

  //  printf("S_reg1\n");

  return result;
}

/* kernel for the integrated photon-echo signal at T<0, and tau > |T|
   t1=tau-|T|, t2=|T|
  the expression using lineshape functions are defined in Eq.1 on Notebook 2.
*/
double f_S_reg2(double t,void *params)
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
double S_reg2(double t1, double t2)
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
  F.function = &f_S_reg2;
  F.params = params;

  gsl_integration_qags(&F, ra, rb, EPSABS/100.0, EPSREL/100.0, NWSPACE, w, &result, &error);
  gsl_integration_workspace_free(w);

  return result;
}

/* kernel for the integrated photon-echo signal at tau>0, T<0, and tau < |T|
   t1=|T|-tau, t2=tau
  the expression using lineshape functions are defined in Eq.2 on Notebook 2.
*/
double f_S_reg3(double t,void *params)
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
double S_reg3(double t1, double t2)
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
  F.function = &f_S_reg3;
  F.params = params;

  gsl_integration_qags(&F, ra, rb, EPSABS/100.0, EPSREL/100.0, NWSPACE, w, &result, &error);
  gsl_integration_workspace_free(w);

  return result;
}

/* Sint determines time-ordering and compute 
   integrated signal using different pathways;
   Note that we use t1/t2 ordering instead of tau/T */
double Sint(double t1, double t2)
{
  double Sval;

  if(t1>=0 && t2>=0) {
    // region 2, HJ's coherence PE ordering (see p. 10 on notebook 2)
    Sval=S_reg2(t1,t2);
  } else if (t1<=0 && t2>=0 && t2>=fabs(t1)) {
    // region 3
    Sval=S_reg3(-1.0*t1,t2+t1);
  } else if (t1>=0 && t2<=0 && t1>=fabs(t2)) {
    // region 1
    Sval=S_reg1(t1+t2,fabs(t2));
  } else {
    Sval=0.0;
  }

  return Sval;
}

// compute Sint(t1,t2) at fixed t2 for 
// t1 from -300fs to 300fs with step size 3fs;
// show result and normalized data
void Sint_fixedT(double t2)
{
  double t1,Smax;
  double St1[4096];
  double tt1[4096];
  double maxt1;
  int i;

  // range to show and step size, in fs
  double t1min=-50.0;
  double t1max=200.0;
  double t1step=5.0;
  size_t npoints;

  // t1 from -500fs to 500fs with step size 5fs
  Smax=0.0;maxt1=-300.0;
  for(i=0,t1=t1min;t1<=t1max;t1=t1+t1step,i++) {
    St1[i]=Sint(t1/CM2FS,t2/CM2FS);
    tt1[i]=t1;
    if(St1[i]>Smax) {
      Smax=St1[i];
      maxt1=t1;
    }
    printf("ta=%.2f, tb=%.2f, Sint=%.18f\n",t1,t2,St1[i]);
  }
  printf("ta=, tb=, Sint=\n");
  npoints=i;

  // also output result nornamized to maximum=1.0
  for(i=0,t1=t1min;t1<=t1max;t1=t1+t1step,i++) {
    printf("ta=%.2f, tb=%.2f, Sint_norm=%.18f\n",t1,t2,St1[i]/Smax);
  }
  printf("ta=, tb=, Sint_norm=\n");

  // find peak shift
  find_max(&maxt1, &Smax, tt1,St1,npoints,maxt1); 

  // the max value at this t2; peak-shift and the projection to t2 axis here
  printf("ta=%.2f, tb=%.2f, Smax=%.18f\n",maxt1,t2,Smax);

}

/* linear absorption lineshape */

/* kernel exp(-g(t)) for the two excitations */
void ExpG_e1(double t, double *ret_r,double *ret_i)
{
  double rp,theta;

  rp=exp(-1.0*Paa(t));
  theta=-1.0*Qaa(t);

  *ret_r=rp*cos(theta);
  *ret_i=rp*sin(theta);
}
void ExpG_e2(double t, double *ret_r,double *ret_i)
{
  double rp,theta;

  rp=exp(-1.0*Pbb(t));
  theta=-1.0*Qbb(t);

  *ret_r=rp*cos(theta);
  *ret_i=rp*sin(theta);
}

/* generate linear abs. spectrum */
void generate_lineshape(spectrum *spec, void (* ExpG)(double, double *,double *),
			const double sigma, const double mu2)
{
  size_t i,niter;
  double delta;
  gsl_vector_complex *Fw = gsl_vector_complex_alloc(FFT_N);
  double *ww = malloc(FFT_N*sizeof(double));
  spectrum *stmp,*spec0;

  gsl_rng *r=gsl_rng_alloc(gsl_rng_ranlux);
  gsl_rng_set (r,(unsigned)time(NULL));

  stmp=spectrum_alloc(1.0,(double)FFT_N,1.0);
  spec0=spectrum_alloc(1.0,(double)FFT_N,1.0);

  // Now we are ready to generate linear absorption spectrum for excitation 1
  Ft_hFourier_forward(ExpG, ww, Fw);
  for(i=0;i<FFT_N;i++) {
    spec0->w[i]=ww[i]; // note that we do not apply eigenvalue shift
    spec0->intensity[i]=GSL_REAL(gsl_vector_complex_get(Fw,i)); // nor do we apply dipole moment factor here
    //    printf("ww=%f %f\n",spec0->w[i],spec0->intensity[i]);
    // sync all spectra coord.
    spec->w[i]=ww[i];
    stmp->w[i]=ww[i];
  }
  // static disorder:
  // assume the static disorder is small so that the eigen wavefunction is not changed,
  // we just average 500 shifted spectrum
  niter = (sigma>1.0) ? 500 : 1;
  spectrum_copy(spec,spec0);
  for(i=0;i<niter;i++) {
    delta=gsl_ran_gaussian(r,sigma);
    spectrum_copy(stmp,spec0);
    spectrum_shift_w(stmp,delta);
    spectrum_addto(spec,stmp);
  }
  // divide by niter+1 to restore the amplitude, also apply the dipole moment factor, and done!
  for(i=0;i<FFT_N;i++) {
    spec->intensity[i]=mu2*spec->intensity[i]/(double)(niter+1);
  }

  /* clean up */
  free(ww);
  gsl_vector_complex_free(Fw);
  gsl_rng_free(r);
  spectrum_free(stmp);
  spectrum_free(spec0);
}

/* show abs. spectra  for both excitations;
   note that center excitation energy is set at 0 for both peaks */
void show_spectra()
{
  double Me1[3],Me2[3];
  double mu2_e1,mu2_e2;
  double Imax,wmax;
  spectrum *spec1;
  spectrum *spec2;
  spectrum *spec_all;
  size_t i;

  spec1=spectrum_alloc(1.0,(double)FFT_N,1.0);
  spec2=spectrum_alloc(1.0,(double)FFT_N,1.0);
  spec_all=spectrum_alloc(1.0,(double)FFT_N,1.0);

  // square of dipole moments; for weighting
  compute_mixed_dipole(cos(theta),sin(theta),Me1);
  compute_mixed_dipole(-1.0*sin(theta),cos(theta),Me2);
  mu2_e1=vec_dot_product(Me1,Me1);
  mu2_e2=vec_dot_product(Me2,Me2);

  generate_lineshape(spec1, ExpG_e1,Delta1,mu2_e1);
  generate_lineshape(spec2, ExpG_e2,Delta1,mu2_e2);

  // normalize to the maximun of B peak
  Imax=-10.0;wmax=0.0;
  for(i=0;i<FFT_N;i++) {
    if(spec1->intensity[i]>Imax) {
      Imax=spec1->intensity[i];
      wmax=spec1->w[i];
    }
  }
  for(i=0;i<FFT_N;i++) {
    spec1->intensity[i]=spec1->intensity[i]/Imax;
    spec2->intensity[i]=spec2->intensity[i]/Imax;
  }

  // print the result
  printf("\n");
  printf("\n");
  printf("Absorption lineshape (w in cm^-1):\n");
  for(i=0;i<FFT_N;i++) {
    printf("wb=%8.2f, Ib=%12.6f, wh=%8.2f, Ih=%12.6f\n",
	   spec1->w[i],spec1->intensity[i],
	   spec2->w[i],spec2->intensity[i]);
  }
  printf("\n");

  spectrum_free(spec1);
  spectrum_free(spec2);
  spectrum_free(spec_all);
}

/* Main program */
int main(int argc, char *argv[])
{
  double T;

  double Me1[3],Me2[3],Me1f[3],Me2f[3];

  /* Set program name for messages.  */
  program_name = argv[0]; 

  if(argc != 9 && argc != 13 && argc != 17 && argc != 21 && argc != 22 && argc != 23) {
    usage();
    exit(EXIT_FAILURE);
  }
  S1=0.0;w1=1000.0;gg1=0.0;phi1=0.0;
  S2=0.0;w2=1000.0;gg2=0.0;phi2=0.0;
  Sj=0.0;wj=1000.0;ggj=0.0;phij=0.0;
  beta=1.4387/atof(argv[1]);
  theta=atof(argv[2]);
  tau0=atof(argv[3])/CM2FS;
  lambda1=atof(argv[4]);
  lambda2=atof(argv[5]);
  corr_coeff=atof(argv[6]);
  Delta1=atof(argv[7]);
  Delta2=atof(argv[8]);

  /* no population transfer by default */
  Kbh=0.0;
  Kpb=0.0;

  if(argc >= 13) {
    S1=atof(argv[9]);
    w1=atof(argv[10]);
    gg1=atof(argv[11]);
    phi1=atof(argv[12]);
  }

  if(argc >= 17) {
    S2=atof(argv[13]);
    w2=atof(argv[14]);
    gg2=atof(argv[15]);
    phi2=atof(argv[16]);
  }

  if(argc >= 21) {
    Sj=atof(argv[17]);
    wj=atof(argv[18]);
    ggj=atof(argv[19]);
    phij=atof(argv[20]);
  }

  if(argc >= 22) {
    Kbh=CM2FS*atof(argv[21]);
  }
  if(argc >= 23) {
    Kpb=CM2FS*atof(argv[22]);
  }

  printf("\n");
  printf("        2C3PPES-HTGAU\n");
  printf("        Impulsive-limit three-pulse\n");
  printf("        photon-echo signals for a dimer system.\n");
  printf("\n");
  printf("        Part of the QDAS package.\n");
  printf("\n");
  printf("        Copyright(C) 2006.\n");
  printf("        Yuan-Chung Cheng <yccheng@berkeley.edu>.\n");
  printf("\n");
  printf("\n");
  printf("2C-3PPES using JSF's Gaussian Ansatz (high-temperature approx.):\n");
  printf("\n");
  printf("                 /t                    /t   /t1\n");
  printf("g(t) = -I*lambda*| M(t1)*dt1 + <dw^2>*|dt1 | M(t2)*dt2\n");
  printf("                 /0                   /0   /0\n");
  printf("\n");
  printf("M(t) = exp(-t^2/tau^2)\n");
  printf("\n");
  printf("For each site:\n");
  printf("%6s %12s %12s %12s\n","Site","lambda (cm^-1)","tau (fs)","<dw^2>");
  printf("%6d %12.4f %12.4f %12.4f\n",1,lambda1,tau0*CM2FS,
	 bath_htgau_wsquare(lambda1,tau0,beta));
  printf("%6d %12.4f %12.4f %12.4f\n",2,lambda2,tau0*CM2FS,
	 bath_htgau_wsquare(lambda2,tau0,beta));
  printf("%6s %12.4f %12.4f %12.4f\n","cross",corr_coeff*sqrt(lambda1*lambda2),tau0*CM2FS,
	 bath_htgau_wsquare(corr_coeff*sqrt(lambda1*lambda2),tau0,beta));
  printf("\n");
  printf("\n");
  printf("Parameter list:\n");
  printf("# beta    = %f (%f K)\n",beta,1.4387/beta);
  printf("# theta   = %f\n",theta);
  printf("# tau0    = %f\n",tau0*CM2FS);
  printf("# lambda1 = %f\n",lambda1);
  printf("# lambda2 = %f\n",lambda2);
  printf("# c       = %f\n",corr_coeff);
  printf("# Delta1  = %f\n",Delta1);
  printf("# Delta2  = %f\n",Delta2);
  printf("# S1      = %f\n",S1);
  printf("# w1      = %f\n",w1);
  printf("# gg1     = %f\n",gg1);
  printf("# phi1    = %f\n",phi1);
  printf("# S2      = %f\n",S2);
  printf("# w2      = %f\n",w2);
  printf("# gg2     = %f\n",gg2);
  printf("# phi2    = %f\n",phi2);
  printf("# Sj      = %f\n",Sj);
  printf("# wj      = %f\n",wj);
  printf("# ggj     = %f\n",ggj);
  printf("# phij    = %f\n",phij);
  printf("# Rbh     = %f\n",Kbh);
  printf("# Rpb     = %f\n",Kpb);
  printf("\n");
 
  // compute transition dipole factors for various pathways

  compute_mixed_dipole(cos(theta),sin(theta),Me1);
  compute_mixed_dipole(-1.0*sin(theta),cos(theta),Me2);
  compute_mixed_dipole(sin(theta),cos(theta),Me1f);
  compute_mixed_dipole(cos(theta),-1.0*sin(theta),Me2f);
  compute_dipole_factors();
  printf("\n");
  printf("Transition dipole moments for exciton states:\n");
  printf("# Me1Me1   = %f\n",vec_dot_product(Me1,Me1));
  printf("# Me2Me2   = %f\n",vec_dot_product(Me2,Me2));
  printf("# Me1fMe1f = %f\n",vec_dot_product(Me1f,Me1f));
  printf("# Me2fMe2f = %f\n",vec_dot_product(Me2f,Me2f));
  printf("# Me1Me1Me2Me2   = %f\n",Me1Me1Me2Me2/CSCALE);
  printf("# Me2Me2Me2fMe2f = %f\n",Me2Me2Me2fMe2f/CSCALE);
  printf("# Me1Me2Me1fMe2f = %f\n",Me1Me2Me1fMe2f/CSCALE);
  printf("\n");
  printf("\n");

  // initialize cached exciton lineshape functions
  exciton_PQ_init();

  // show linear absorption spectra for both excitations
  show_spectra();

  // We turnoff the error handler, which means that the number generated
  // might have errors in them...
  gsl_set_error_handler_off();

  //  for(T=300.0;T>0.0;T=T-25.0) {
  //    Sint_fixedT(T);
  //  }

  printf("Three-pulse photon-echo signals:\n");
  for(T=70.0;T<=400.0;T=T+5.0) {
    Sint_fixedT(T);
  }

  gsl_set_error_handler (NULL);

  return 1;
}

/*
 * $Log$
 * Revision 1.5  2007/02/21 06:13:05  platin
 *   - minor changes.
 *
 * Revision 1.4  2007/02/02 08:07:22  platin
 *
 *   - make static disorder correlated.
 *
 *
 */
