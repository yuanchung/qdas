/***************************************************
 * bath_js03.c
 *
 * Bath functions that implement the JS03 BChl
 * bath correlation functions for general use.
 *
 * By Yuan-Chung Cheng <yuanchung@ntu.edu.tw>
 *
 *
 ***************************************************/

#define _GNU_SOURCE
//#define HAVE_INLINE

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
#include <gsl/gsl_integration.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_fft_complex.h>

#include "qdas.h"
#include "aux.h"
#include "bath_js03.h"

/* Use the following inline gsl_ functions to optimize performance */
inline gsl_complex fgsl_matrix_complex_get(const gsl_matrix_complex * m, 
					   const size_t i, const size_t j)
{
  return *(gsl_complex *)(m->data + 2*(i * m->tda + j)) ;
} 

inline void fgsl_matrix_complex_set(gsl_matrix_complex * m, 
				    const size_t i, const size_t j, const gsl_complex x)
{
  *(gsl_complex *)(m->data + 2*(i * m->tda + j)) = x ;
}

inline double fgsl_matrix_get(const gsl_matrix * m, const size_t i, const size_t j)
{
  return m->data[i * m->tda + j] ;
} 

inline void fgsl_matrix_set(gsl_matrix * m, const size_t i, const size_t j, const double x)
{
  m->data[i * m->tda + j] = x ;
}

// Defined Constants

#define NWSPACE (100000)
#define EPSABS (1e-5)
#define EPSREL (1e-5)

/* define grid points in the Ct table; note that the range should >> 1/wc */
#define TablePointsMax 30000
#define TableMax  (100.0)
// the following constant controls the step size in short times.
// The initial step size will be TableStepFactor*(1/wc).
#define TableStepFactor (0.001)

// Global variables

// system-bath interacting terms and bath correlation 
// functions in the localized basis; we handle at most
// 256 terms here.
// |n><m|*<B(tau)B(0)>, note the maximal number of terms
// for a N state system is N*(N+1)/2
int BATH_JS03OpKet[256];
int BATH_JS03OpBra[256];
js03_bath_func BATH_JS03BathFunc[256];
// Q is a static matrix with dimension (n^4)*num_bath_ops, which
// holds the coefficient that converts the correlations in the
// site-localized basis to the exciton matrix.
gsl_matrix *BATH_JS03OpQ;
gsl_vector *BATH_JS03Ct;

// normalized spectral function taken from J&S
inline double bath_js03_Jw(double omega, double wc)
{
  return ((0.5*omega+0.58*omega*omega/wc)*exp(-1.0*omega/wc));
}

/* 
   bath correlation functions; note that either Ct_r nor Ct_i here is
   scaled by the coupling strength gamma. This enables us to save
   unscaled values in the cache files.
*/
/* real part of C(t) */
inline double bath_js03_f_Ct_r(double w, void *params)
{
  double aw;
  double tau,beta,wc;
  double *p = (double *) params;;
  double cothx;

  if (w == 0.0) return 0.0;

  tau = p[0];
  beta = p[1];
  wc = p[2];

  cothx=((1.0+exp(-1.0*beta*w))/(1.0-exp(-1.0*beta*w)));
  aw=(bath_js03_Jw(w,wc)*cothx*cos(w*tau));
  //    printf("HERE:%f %f %f %f %f\n",w,wc,aw,bath_js03_Jw(w,wc),coth(0.5*beta*w));
  return aw;
}

double bath_js03_Ct_r(double tau, double beta, double wc)
{
  gsl_integration_workspace *w;
  double result, error;
  double params[3];
  gsl_function F;

  double ra=0.0,rb=10000.0; /* integration range */

  w=gsl_integration_workspace_alloc(NWSPACE);

  /* fill in params */
  params[0]=tau;
  params[1]=beta;
  params[2]=wc;

  F.function = &bath_js03_f_Ct_r;
  F.params = params;

  gsl_integration_qags(&F, ra, rb, EPSABS/100.0, EPSREL/100.0, NWSPACE, w, &result, &error);
  gsl_integration_workspace_free(w);

  return result;
}

/* imaginary part */
double bath_js03_Ct_i(double tau, double wc)
{
  double d1,d2,d3;

  /* We got analytical form over w for this one,
     for Jang and Silbey's spectral function,
     J(w)=(0.5*w+0.58*w*w/wc)*exp(-w/wc),
     Ci(t) = 0.16*wc^3*t*(-28+wc^2*t^2)/(1+wc^2*t^2)^3 */
  d1=wc*wc*wc*tau;
  d2=wc*wc*tau*tau;
  d3=1+d2;
  d3=d3*d3*d3;
  return (0.16*d1*(d2-28.0)/d3);
}

/* the following functions are used to cache and speed up the evaluation
   of the bath correlation function */
void bath_js03_ct_init(js03_bath_func *f, double beta, double gamma, double wc)
{
  double t[TablePointsMax],ct_r[TablePointsMax],ct_i[TablePointsMax];
  char fname[256];
  FILE *fh;
  size_t i;
  size_t npoints;

  // assign gamma, beta, wc
  f->gamma=gamma;
  f->beta=beta;
  f->wc=wc;

  // the cache file name
  snprintf(fname,255,"js03_ct_table_%09.6f_%06.2f.cache",beta,wc);

#ifdef BATH_DEBUG
  printf("Initializing bath correlation function ...\n");
  printf("  Beta = %f cm (%5.2f K)\n",beta,1.4387/beta);
  printf("  J&S Spectral Density:\n");
  printf("    Wc = %f\n",wc);
  printf("    Coupling Strength Gamma = %f\n",gamma);
  printf("    Cache file: %s\n",fname);
  printf("\n");
#endif

  fh=fopen(fname,"r");
  if(fh) {
    // found the cache file, read from it
#ifdef BATH_DEBUG
    printf("Found the cache file, use data points in it.\n");
    fflush(stdout);
#endif
    i=0;
    while(!feof(fh)) {
      fscanf(fh,"%lf %lf %lf\n",t+i,ct_r+i,ct_i+i);
      //      printf("%8.2f\t%8.4f\t%8.4f\n",w[i],ct_r[i],ct_i[i]);
      i++;
    }
    // update the number of points
    npoints = i;
    fclose(fh);
  } else {
    double tau;
    // cannot find the cache file, need to re-generate the data
    printf("Cache file not found, re-generate data points.\n");
    printf("  t\t  Re[C(t)]\t  Im[C(t)]\n");
    fflush(stdout);

    // We turnoff the error handler, which means that the number generated
    // might have errors in them...
    gsl_set_error_handler_off();

    tau=0.0;i=0;
    while(tau<=TableMax && i<TablePointsMax) {
      t[i]=tau;
      ct_r[i]=bath_js03_Ct_r(tau,beta,wc);
      ct_i[i]=bath_js03_Ct_i(tau,wc);
      printf("%f %f %f\n",tau,ct_r[i],ct_i[i]);

      i++;
      if(tau<(10.0/wc)) {
	tau=tau + TableStepFactor/wc;
      } else {
	tau=tau*1.1;
      }
      fflush(stdout);
    }
    gsl_set_error_handler (NULL);

    npoints=i;

    // keep a copy in the cache file
    printf("Writing cache file %s...\n",fname);
    fh=fopen(fname,"w");
    for(i=0;i<npoints;i++) {
      fprintf(fh,"%.12f %.12f %.12f\n",t[i],ct_r[i],ct_i[i]);
    }
    fclose(fh);
  }

  // prepare the GSL spline functions, and done
  f->BATH_JS03_Ct_r_Spline_Acc=gsl_interp_accel_alloc();
  f->BATH_JS03_Ct_i_Spline_Acc=gsl_interp_accel_alloc();

  f->BATH_JS03_Ct_r_Spline=gsl_spline_alloc (gsl_interp_cspline, npoints);
  gsl_spline_init(f->BATH_JS03_Ct_r_Spline,t,ct_r,npoints);
  f->BATH_JS03_Ct_i_Spline=gsl_spline_alloc (gsl_interp_cspline, npoints);
  gsl_spline_init(f->BATH_JS03_Ct_i_Spline,t,ct_i,npoints);

}

/* these two functions do the real work */
double bath_js03_ct_r_cached(js03_bath_func *f, double tt)
{
  double t;
  t=fabs(tt);
  return (f->gamma*gsl_spline_eval(f->BATH_JS03_Ct_r_Spline,t,f->BATH_JS03_Ct_r_Spline_Acc));
}

double bath_js03_ct_i_cached(js03_bath_func *f, double t)
{
  if(t>=0.0) {
    return (f->gamma*gsl_spline_eval(f->BATH_JS03_Ct_i_Spline,t,f->BATH_JS03_Ct_i_Spline_Acc));
  } else {
    return -1.0*(f->gamma*gsl_spline_eval(f->BATH_JS03_Ct_i_Spline,-1.0*t,f->BATH_JS03_Ct_i_Spline_Acc));
  }
}

/* bath lineshape function,
          /t    /t1
   g(t) = | dt1 | dt2 C(t1-t2)
          /0    /0

          /t     
        = | dtau C(tau)*(t-tau)
          /0     

   Note that either gt_r nor gt_i is scaled.
*/
// use this temporary container to hold required parameters for integrand
struct tmp_bft_struct 
{
  js03_bath_func *bf;
  double t;
};
typedef struct tmp_bft_struct tmp_bft;

double bath_js03_f_gt_r(double tau, void *params)
{
  js03_bath_func *bf;
  double Ctau,t;
  tmp_bft *bft=params;

  bf=bft->bf;
  t=bft->t;

  Ctau=gsl_spline_eval(bf->BATH_JS03_Ct_r_Spline,fabs(tau),bf->BATH_JS03_Ct_r_Spline_Acc);
  return (Ctau*(t-tau));
}

double bath_js03_gt_r(double t, js03_bath_func *bf)
{
  gsl_integration_workspace *w;
  double result, error;
  gsl_function F;
  tmp_bft bft;

  w=gsl_integration_workspace_alloc(NWSPACE);

  /* fill in params */
  bft.bf=bf;
  bft.t=t;
  F.function = &bath_js03_f_gt_r;
  F.params = &bft;

  /* integrate up to time t */
  gsl_integration_qags (&F, 0.0, t, EPSABS, EPSREL, NWSPACE, w, &result, &error);
  gsl_integration_workspace_free(w);

  return result;
}

double bath_js03_f_gt_i(double tau, void *params)
{
  js03_bath_func *bf;
  double Ctau,t;
  tmp_bft *bft=params;

  bf=bft->bf;
  t=bft->t;

  if(tau>=0) {
    Ctau=gsl_spline_eval(bf->BATH_JS03_Ct_i_Spline,tau,bf->BATH_JS03_Ct_i_Spline_Acc);
  } else {
    Ctau=-1.0*gsl_spline_eval(bf->BATH_JS03_Ct_i_Spline,-1.0*tau,bf->BATH_JS03_Ct_i_Spline_Acc);
  }
  return (Ctau*(t-tau));
}

double bath_js03_gt_i(double t, js03_bath_func *bf)
{
  gsl_integration_workspace *w;
  double result, error;
  gsl_function F;
  tmp_bft bft;

  w=gsl_integration_workspace_alloc(NWSPACE);

  /* fill in params */
  bft.bf=bf;
  bft.t=t;
  F.function = &bath_js03_f_gt_i;
  F.params = &bft;

  /* integrate up to time t */
  gsl_integration_qags (&F, 0.0, t, EPSABS, EPSREL, NWSPACE, w, &result, &error);
  gsl_integration_workspace_free(w);

  return result;
}

void bath_js03_gt_init(js03_lineshape_func *f, js03_bath_func *bf)
{
  double t[TablePointsMax],gt_r[TablePointsMax],gt_i[TablePointsMax];
  char fname[256];
  FILE *fh;
  size_t i;
  size_t npoints;
  double wc,beta;

  wc=bf->wc;
  beta=bf->beta;
  f->bf=bf;
  f->gamma=bf->gamma;

  // the cache file name
  snprintf(fname,255,"js03_gt_table_%09.6f_%06.2f.cache",beta,wc);

  fh=fopen(fname,"r");
  if(fh) {
    // found the cache file, read from it
#ifdef BATH_DEBUG
    printf("Found the cache file, use data points in it.\n");
    fflush(stdout);
#endif
    i=0;
    while(!feof(fh)) {
      fscanf(fh,"%lf %lf %lf\n",t+i,gt_r+i,gt_i+i);
      //      printf("%8.2f\t%8.4f\t%8.4f\n",w[i],gt_r[i],gt_i[i]);
      i++;
    }
    // update the number of points
    npoints = i;
    fclose(fh);
  } else {
    double tau;
    // cannot find the cache file, need to re-generate the data
    printf("Cache file for g(t) not found, re-generate data points.\n");
    printf("  t\t  Re[g(t)]\t  Im[g(t)]\n");
    fflush(stdout);
    // We turnoff the error handler, which means that the number generated
    // might have errors in them...
    gsl_set_error_handler_off();

    i=0; tau=0.0;
    while(tau<=TableMax && i<TablePointsMax) {
      t[i]=tau;
      gt_r[i]=bath_js03_gt_r(tau,bf);
      gt_i[i]=bath_js03_gt_i(tau,bf);
      printf("%f %f %f\n",tau,gt_r[i],gt_i[i]);

      i++;
      if(tau<(10.0/wc)) {
	tau=tau + TableStepFactor/wc;
      } else {
	tau=tau*1.1;
      }
      fflush(stdout);
    }
    gsl_set_error_handler (NULL);
    npoints=i;
    // keep a copy in the cache file
    printf("Writing cache file %s...\n",fname);
    fh=fopen(fname,"w");
    for(i=0;i<npoints;i++) {
      fprintf(fh,"%.12f %.12f %.12f\n",t[i],gt_r[i],gt_i[i]);
    }
    fclose(fh);
  }

  // prepare the GSL spline functions, and done
  f->BATH_JS03_gt_r_Spline_Acc=gsl_interp_accel_alloc();
  f->BATH_JS03_gt_i_Spline_Acc=gsl_interp_accel_alloc();

  f->BATH_JS03_gt_r_Spline=gsl_spline_alloc (gsl_interp_cspline, npoints);
  gsl_spline_init(f->BATH_JS03_gt_r_Spline,t,gt_r,npoints);
  f->BATH_JS03_gt_i_Spline=gsl_spline_alloc (gsl_interp_cspline, npoints);
  gsl_spline_init(f->BATH_JS03_gt_i_Spline,t,gt_i,npoints);

}

void js03_lineshape_func_free(js03_lineshape_func *f)
{
  gsl_interp_accel_free(f->BATH_JS03_gt_r_Spline_Acc);
  gsl_interp_accel_free(f->BATH_JS03_gt_i_Spline_Acc);
  gsl_spline_free(f->BATH_JS03_gt_r_Spline);
  gsl_spline_free(f->BATH_JS03_gt_i_Spline);
}

/* these two functions do the real work for g(t) */
double bath_js03_gt_r_cached(js03_lineshape_func *f, double tt)
{
  double t;
  t=fabs(tt);
  return (f->gamma*gsl_spline_eval(f->BATH_JS03_gt_r_Spline,t,f->BATH_JS03_gt_r_Spline_Acc));
}

double bath_js03_gt_i_cached(js03_lineshape_func *f, double t)
{
  if(t>=0.0) {
    return (f->gamma*gsl_spline_eval(f->BATH_JS03_gt_i_Spline,t,f->BATH_JS03_gt_i_Spline_Acc));
  } else {
    return (f->gamma*gsl_spline_eval(f->BATH_JS03_gt_i_Spline,-1.0*t,f->BATH_JS03_gt_i_Spline_Acc));
  }
}

/* interface functions */
int bath_js03_init_params(const size_t nsize, const double beta, 
			  const size_t bath_nparams, const double *bath_params)
{
  int i,idx;
  int nterms;

  // initialize system-op arrays; 
  // zero can be used to indicate the end of the coupling terms
  for(i=0;i<256;i++) {
    BATH_JS03OpBra[i]=0;
    BATH_JS03OpKet[i]=0;
  }

  printf("Use J&S's BChl spectral function:\n");
  printf("J(w)=gamma*(0.5*w+0.58*w^2/wc)*exp(-w/wc)\n");
  printf("\n");
  printf("System-Op and J(w) for each site:\n");
  printf("\n");
  printf("  %16s %12s %12s\n","System-Op","gamma","wc (cm^-1)");

  // input in the mod_params array should be in a table form
  // n m gamma wc
  // parse them to the global variables
  nterms=bath_nparams/4; // number of system-bath coupling terms
  for(i=0;i<nterms;i++) {
    double gamma,wc;
    int na,nb;
    idx=i*4;
    // system OP represented by ket state |n> and bra state <m|, i.e. |n><m|.
    // We make sure n<=m
    na=bath_params[idx];
    nb=bath_params[idx+1];
    BATH_JS03OpKet[i]=((na < nb) ? na : nb);
    BATH_JS03OpBra[i]=((na < nb) ? nb : na);
    gamma=bath_params[idx+2];
    wc=bath_params[idx+3];
    if(BATH_JS03OpBra[i]<=0 || BATH_JS03OpKet[i]<=0) {
      printf("Error: can't parse system-bath interaction terms.\n");
      exit(EXIT_FAILURE);
    }
    bath_js03_ct_init(BATH_JS03BathFunc+i,beta,gamma,wc);

    // print out the parameters...
    if(BATH_JS03OpBra[i] == BATH_JS03OpKet[i]) {
      // diagonal e-ph couplings
      printf("  |%d><%d|            %12.4f %12.4f\n",
	     BATH_JS03OpKet[i],BATH_JS03OpBra[i],gamma,wc);
    } else {
      // off-diagonal terms, note that we will automatically 
      // symmetrize the system part in all evaluations later.
      printf("  |%d><%d|+|%d><%d|   %12.4f %12.4f\n",
	     BATH_JS03OpKet[i],BATH_JS03OpBra[i],BATH_JS03OpBra[i],BATH_JS03OpKet[i],
	     gamma,wc);
    }
  }
  printf("\n");

  return 0;
}

int bath_js03_free_params()
{
  int i;
  
  i=0;
  while(BATH_JS03OpBra[i]>0) {
    gsl_interp_accel_free(BATH_JS03BathFunc[i].BATH_JS03_Ct_r_Spline_Acc);
    gsl_interp_accel_free(BATH_JS03BathFunc[i].BATH_JS03_Ct_i_Spline_Acc);
    gsl_spline_free(BATH_JS03BathFunc[i].BATH_JS03_Ct_r_Spline);
    gsl_spline_free(BATH_JS03BathFunc[i].BATH_JS03_Ct_i_Spline);
    i++;
  }

  return 0;
}

/* The following functions handles the correlation 
   tensor in the exciton matrix */
void bath_js03_init_OpQ(int ndim, gsl_matrix *U)
{
  int nops,i;
  size_t n,m,a,b,c,d;
  size_t idxa,idxb,idxc;
  double dd;

  /* we first count the number of bath operators */
  nops=0;
  while(BATH_JS03OpBra[nops]>0) {
    nops++;
  }

  /* allocate the matrix, totally n^2*n^2 rows... */
  BATH_JS03OpQ=gsl_matrix_alloc(ndim*ndim*ndim*ndim,nops);
  gsl_matrix_set_zero(BATH_JS03OpQ);
  BATH_JS03Ct=gsl_vector_alloc(nops); // as a workspace

  /* construct matrix elements, no tricks here */
  for(i=0;i<nops;i++) {
    n=BATH_JS03OpKet[i]-1;
    m=BATH_JS03OpBra[i]-1; // note: m>=n
    for(a=0,idxa=0;a<ndim;a++,idxa+=ndim*ndim*ndim) {
      for(b=0,idxb=0;b<ndim;b++,idxb+=ndim*ndim) {
	for(c=0,idxc=0;c<ndim;c++,idxc+=ndim) {
	  for(d=0;d<ndim;d++) {
	    //	    printf("(n,m,a,b,c,d) = (%d,%d,%d,%d,%d,%d)\n",n,m,a,b,c,d);
	    if(n==m) {
	      // diagonal term
	      dd=fgsl_matrix_get(U,n,a)*fgsl_matrix_get(U,n,b)*
		fgsl_matrix_get(U,n,c)*fgsl_matrix_get(U,n,d);
	      fgsl_matrix_set(BATH_JS03OpQ,idxa+idxb+idxc+d,i,dd);
	    } else {
	      // off-diagonal term
	      dd=fgsl_matrix_get(U,n,a)*fgsl_matrix_get(U,m,b)*
		fgsl_matrix_get(U,n,c)*fgsl_matrix_get(U,m,d);
	      dd+=fgsl_matrix_get(U,n,a)*fgsl_matrix_get(U,m,b)*
		fgsl_matrix_get(U,m,c)*fgsl_matrix_get(U,n,d);
	      dd+=fgsl_matrix_get(U,m,a)*fgsl_matrix_get(U,n,b)*
		fgsl_matrix_get(U,m,c)*fgsl_matrix_get(U,n,d);
	      dd+=fgsl_matrix_get(U,m,a)*fgsl_matrix_get(U,n,b)*
		fgsl_matrix_get(U,n,c)*fgsl_matrix_get(U,m,d);
	      fgsl_matrix_set(BATH_JS03OpQ,idxa+idxb+idxc+d,i,dd);
	    }
	  }
	}
      }
    }
  }

  printf("OpQ =\n");
  gsl_matrix_lprint(BATH_JS03OpQ);
  printf("\n");
  printf("\n");

}

void bath_js03_free_OpQ()
{
  gsl_matrix_free(BATH_JS03OpQ);
  gsl_vector_free(BATH_JS03Ct);
}

/* evaluate the exciton-basis correlation tensor (in an array form) at time t,
   input: t: time
          qr: a n^4 long arrays, q_r, for the real part
              of Q(t)=<q_ab(t)q_cd(0)>, note that the 
              element (a,b; c,d) is stored at array index a*n^3+b*n^2+c*n+d. 
    details see Notebook 1, p. 68. */
void bath_js03_Qt_r(gsl_vector *qr, double t)
{
  int i;
  
  /* C(t) */
  i=0;
  while(BATH_JS03OpBra[i]>0) {
    gsl_vector_set(BATH_JS03Ct,i,bath_js03_ct_r_cached(BATH_JS03BathFunc+i,t));
    i++;
  }
  /* Q(t) = OpQ * C(t) */
  gsl_blas_dgemv (CblasNoTrans, 1.0, BATH_JS03OpQ, BATH_JS03Ct, 0.0, qr);
  //  printf("OpQ =\n");
  //  matrix_print(BATH_JS03OpQ);
  //  printf("\n");
  //  printf("Ct =\n");
  //  vector_print(BATH_JS03Ct);
  //  printf("\n");
  //  printf("qr =\n");
  ///  vector_print(qr);
    //  printf("\n");
}
/* imaginary part */
void bath_js03_Qt_i(gsl_vector *qi, double t)
{
  int i;
  
  /* C(t) */
  i=0;
  while(BATH_JS03OpBra[i]>0) {
    gsl_vector_set(BATH_JS03Ct,i,bath_js03_ct_i_cached(BATH_JS03BathFunc+i,t));
    i++;
  }
  /* Q(t) = OpQ * C(t) */
  gsl_blas_dgemv (CblasNoTrans, 1.0, BATH_JS03OpQ, BATH_JS03Ct, 0.0, qi);
}

/* real part and imagnary part of the vibrational contributions */
double bath_js03_gv_r(double t, vibmodes *vibs,double beta)
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

double bath_js03_gv_i(double t, vibmodes *vibs)
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


#ifdef MAIN
int main()
{
  // unit in cm^-1
  double beta=0.018685; // 77K
  double wc=100.0;
  //  double beta=0.2;
  //  double wc=50.0;

  double t,dt;

  return 1;

}
#endif

/*
 * $Log$
 * Revision 1.11  2007/07/18 23:11:14  platin
 *   - remove unused code.
 *
 * Revision 1.10  2007/06/21 22:11:03  platin
 *
 *   - code clean up.
 *
 * Revision 1.9  2007/01/12 07:20:32  platin
 *
 *   - increase percision.
 *
 * Revision 1.8  2006/10/31 23:24:47  platin
 *
 *   - increase the size of cached table.
 *
 * Revision 1.7  2006/10/31 22:17:12  platin
 *
 *   - fix negative time part in cached g(t).
 *   - save unscaled g(t) in cache file too.
 *
 * Revision 1.6  2006/10/31 19:35:32  platin
 *
 *   - add support for lineshape function g(t) in the bath_js03 module.
 *
 * Revision 1.5  2006/10/27 23:02:59  platin
 *
 *   - commit modifications up to 10/27/2006.
 *
 * Revision 1.4  2006/06/28 17:14:11  platin
 *
 *   - add the dm3pes module that computes the three-pulse photon-echo
 *     signals using Domcke's density-matrix based method.
 *   - aux functions used to implement dm3pes.
 *
 * Revision 1.3  2006/06/08 17:49:11  platin
 *
 *   - add density matrix handeling functions (from qbdyn).
 *
 * Revision 1.2  2006/05/26 19:19:30  platin
 *
 *   - add dynamics module.
 *   - revise the params.c, use a more reasonable scheme to handle
 *     parameters for different modules.
 *
 * Revision 1.1.1.1  2006/05/24 00:42:19  platin
 *
 *   - initial import of the qdas package.
 *   - qdas stands for quantum dynamics and spectroscopy.
 *   - basic code inherited from the "lineshape" package.
 *
 *
 *
 */
