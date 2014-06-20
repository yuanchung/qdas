/***************************************************
 * bathtype_mt99art.c
 *
 * A bathtype module that provides shared variables 
 * and functions for bath modules that implement Meier
 * and Tannor's exponential-parameterized artificial
 * bath correlation functions.
 *
 * By Yuan-Chung Cheng <yuanchung@ntu.edu.tw>
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
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_blas.h>

#include "qdas.h"
#include "aux.h"
#include "bath.h"
#include "bathtype_mt99art.h"

// Defined Constants

// Global variables

// system-bath interacting terms and bath correlation 
// functions in the localized basis; we handle at most
// 256 terms here.
// |n><m|*<B(tau)B(0)>, note the maximal number of terms
// for a N state system is N*(N+1)/2
size_t BATH_MT99ART_NOps=0; /* number of bath operator terms */
mt99art_bath_func BATH_MT99ARTBathFunc[256];
// interacting system operators in the exciton basis
gsl_matrix_complex *BATH_MT99ARTOpSbar[256];

/* 
   bath correlation functions, return a complex number:
   The definition adapted here is different from Meier and Tannor.
   here:
   C_r(t) = \int_0^\infy dw J(w)*coth(beta*w/2)*cos(w*t),
   C_i(t) = \int_0^\infy dw J(w)*sin(w*t).
   So, C(t) here = PI*(a(t)+i*b(t)). with a(t) and b(t) defined in M&T.
*/

gsl_complex bath_mt99art_ct(mt99art_bath_func *f, double tt)
{
  double t;
  gsl_complex at,bt,ztmp;
  int i;
  
  t=fabs(tt);

  at=gsl_complex_rect(0.0,0.0);
  for(i=0;i<f->Nr;i++) {
    ztmp=gsl_complex_mul_real(f->Gr[i],t);
    ztmp=gsl_complex_exp(ztmp);
    ztmp=gsl_complex_mul(f->Ar[i],ztmp);
    at=gsl_complex_add(at,ztmp);
  }
  bt=gsl_complex_rect(0.0,0.0);
  for(i=0;i<f->Ni;i++) {
    ztmp=gsl_complex_mul_real(f->Gi[i],t);
    ztmp=gsl_complex_exp(ztmp);
    ztmp=gsl_complex_mul(f->Ai[i],ztmp);
    bt=gsl_complex_add(bt,ztmp);
  }

  if(tt<0.0) {
    // extra PI and minus sign factor because different definition of C(t) in M&T
    return gsl_complex_rect(GSL_REAL(at)*M_PI,GSL_REAL(bt)*M_PI);
  } else {
    return gsl_complex_rect(GSL_REAL(at)*M_PI,-1.0*GSL_REAL(bt)*M_PI);
  }

}

/* bath lineshape function,
   g(t)= \int_0^t dt1 \int_0^t1 dt2 C(t1-t2).
   Note that definition of C(t) is different from M&T.
   bath_mt99art_gt_r(): real part
   bath_mt99art_gt_i(): imaginary part
*/
double bath_mt99art_gt_r(mt99art_bath_func *f, double tt)
{
  gsl_complex Gr,Grt;
  gsl_complex gt_r,ztmp;
  double t;
  int i;
  
  t=fabs(tt);

  gt_r=gsl_complex_rect(0.0,0.0);
  
  for(i=0;i<f->Nr;i++) {
    Gr=f->Gr[i];
    Grt=gsl_complex_mul_real(Gr,t);
    ztmp=gsl_complex_exp(Grt);
    ztmp=gsl_complex_sub(ztmp,Grt);
    ztmp=gsl_complex_sub(ztmp,gsl_complex_rect(1.0,0.0));
    ztmp=gsl_complex_div(ztmp,gsl_complex_mul(Gr,Gr));
    ztmp=gsl_complex_mul(f->Ar[i],ztmp);
    gt_r=gsl_complex_add(gt_r,ztmp);
  }

  //  printf("G(t)_r = %f %f\n",GSL_REAL(gt_r),GSL_IMAG(gt_r));

  // extra PI factor because different definition of C(t) in M&T
  return GSL_REAL(gt_r)*M_PI;
}

double bath_mt99art_gt_i(mt99art_bath_func *f, double tt)
{
  gsl_complex Gi,Git;
  gsl_complex gt_i,ztmp;
  double t;

  int i;

  t=fabs(tt);
  
  gt_i=gsl_complex_rect(0.0,0.0);
  
  for(i=0;i<f->Ni;i++) {
    Gi=f->Gi[i];
    Git=gsl_complex_mul_real(Gi,t);
    ztmp=gsl_complex_exp(Git);
    ztmp=gsl_complex_sub(ztmp,Git);
    ztmp=gsl_complex_sub(ztmp,gsl_complex_rect(1.0,0.0));
    ztmp=gsl_complex_div(ztmp,gsl_complex_mul(Gi,Gi));
    ztmp=gsl_complex_mul(f->Ai[i],ztmp);
    gt_i=gsl_complex_add(gt_i,ztmp);
  }

  //  printf("G(t)_i = %f %f\n",GSL_REAL(gt_i),GSL_IMAG(gt_i));

  if(tt<0.0) {
    // extra PI factor because different definition of C(t) in M&T
    return GSL_REAL(gt_i)*M_PI;
  } else {
    return -1.0*GSL_REAL(gt_i)*M_PI;
  }
}

/* integrated correlation.
   h(t)= \int_0^t dt1 C(t1).
   Note that definition of C(t) is different from M&T.
   bath_mt99art_ht_r(): real part
   bath_mt99art_ht_i(): imaginary part
*/
double bath_mt99art_ht_r(mt99art_bath_func *f, double tt)
{
  gsl_complex Gr,Grt;
  gsl_complex ht_r,ztmp;
  double t;
  int i;
  
  t=fabs(tt);

  ht_r=gsl_complex_rect(0.0,0.0);
  
  for(i=0;i<f->Nr;i++) {
    Gr=f->Gr[i];
    Grt=gsl_complex_mul_real(Gr,t);
    ztmp=gsl_complex_exp(Grt);
    ztmp=gsl_complex_sub(ztmp,gsl_complex_rect(1.0,0.0));
    ztmp=gsl_complex_div(ztmp,Gr);
    ztmp=gsl_complex_mul(f->Ar[i],ztmp);
    ht_r=gsl_complex_add(ht_r,ztmp);
  }

  //  printf("H(t)_r = %f %f\n",GSL_REAL(ht_r),GSL_IMAG(ht_r));

  // extra PI factor because different definition of C(t) in M&T
  return GSL_REAL(ht_r)*M_PI;
}

double bath_mt99art_ht_i(mt99art_bath_func *f, double tt)
{
  gsl_complex Gi,Git;
  gsl_complex ht_i,ztmp;
  double t;
  int i;
  
  t=fabs(tt);

  ht_i=gsl_complex_rect(0.0,0.0);
  
  for(i=0;i<f->Ni;i++) {
    Gi=f->Gi[i];
    Git=gsl_complex_mul_real(Gi,t);
    ztmp=gsl_complex_exp(Git);
    ztmp=gsl_complex_sub(ztmp,gsl_complex_rect(1.0,0.0));
    ztmp=gsl_complex_div(ztmp,Gi);
    ztmp=gsl_complex_mul(f->Ai[i],ztmp);
    ht_i=gsl_complex_add(ht_i,ztmp);
  }

  //  printf("H(t)_i = %f %f\n",GSL_REAL(ht_i),GSL_IMAG(ht_i));

  if(tt<0) {
    // extra PI factor because different definition of C(t) in M&T
    return GSL_REAL(ht_i)*M_PI;
  } else {
    return -1.0*GSL_REAL(ht_i)*M_PI;
  }

}

/* The following function initializes the interacting system operators
   in the exciton basis */
void bath_mt99art_init_OpSbar(const qdas_keys *keys, const gsl_matrix *U)
{
  int i;
  gsl_matrix_complex *UU=gsl_matrix_complex_alloc(keys->nsize,keys->nsize);
  gsl_matrix_complex *mtmp=gsl_matrix_complex_alloc(keys->nsize,keys->nsize);
  gsl_matrix *Stmp=gsl_matrix_alloc(keys->nsize,keys->nsize);;
  double dtmp;
  size_t ndim;
  int n,m;

  for(i=0;i<keys->bathgroup_nterms;i++) {
    int bath_mod=keys->bathgroup[i].bath_mod;
    if(strcmp(bath_get_keyword(bath_mod),"JS03ART") != 0 &&
       strcmp(bath_get_keyword(bath_mod),"MBOART") != 0 &&
       strcmp(bath_get_keyword(bath_mod),"MT99ART") != 0 &&
       strcmp(bath_get_keyword(bath_mod),"MT99OHM") != 0 &&
       strcmp(bath_get_keyword(bath_mod),"OHMART") != 0 ) {
      printf("Bath-type \"%s\" not supported!\n",bath_get_keyword(bath_mod));
      exit(EXIT_FAILURE);
    }
  }

  ndim=keys->nsize;

  /* convert U to complex form */
  gsl_matrix_complex_copy_real(UU,U);

  /* we allocate/initialize the Sbar matrices for all system-bath interaction terms */
  for(i=0;i<BATH_MT99ART_NOps;i++) {
    /* get the real bath Op matrix in the site basis */
    gsl_matrix_memcpy(Stmp,BATH_MT99ARTBathFunc[i].S);
    /* fill in two-exciton matrix elements in Stmp, now in localized basis */
    gsl_matrix_fill_2es_from_1es(Stmp, keys->ntes, keys->tes_list);
    /* convert into a complex matrix */
    BATH_MT99ARTOpSbar[i]=gsl_matrix_complex_alloc(ndim,ndim);
    gsl_matrix_complex_copy_real(BATH_MT99ARTOpSbar[i],Stmp);
    /* and do the transform to the eigenbasis; Sbar=U^dagger*S*U */
    gsl_blas_zgemm(CblasConjTrans, CblasNoTrans, gsl_complex_rect(1.0,0.0), 
		   UU, BATH_MT99ARTOpSbar[i], gsl_complex_rect(0.0,0.0), mtmp);
    gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, gsl_complex_rect(1.0,0.0), 
		   mtmp, UU, gsl_complex_rect(0.0,0.0), BATH_MT99ARTOpSbar[i]);

    if(keys->printlv >= 3) {
      printf("\n");
      printf("S[%d]=\n",i);
      for(n=0;n<ndim;n++) {
	for(m=0;m<ndim;m++) {
	  dtmp=GSL_REAL(gsl_matrix_complex_get(BATH_MT99ARTOpSbar[i],n,m));
	  printf("%12.6f + %12.6fI ",
		 GSL_REAL(gsl_matrix_complex_get(BATH_MT99ARTOpSbar[i],n,m)),
		 GSL_IMAG(gsl_matrix_complex_get(BATH_MT99ARTOpSbar[i],n,m)));
	}
	printf("\n");
      }    
      printf("\n");
    }
  } // for(i=0;i<BATH_MT99ART_NOps;i++) 
  
  gsl_matrix_complex_free(mtmp);
  gsl_matrix_complex_free(UU);
  gsl_matrix_free(Stmp);

}

void bath_mt99art_free_OpSbar()
{
  int i;
  for(i=0;i<BATH_MT99ART_NOps;i++) {
    gsl_matrix_complex_free(BATH_MT99ARTOpSbar[i]);
  }
}

int bath_mt99art_free_params()
{
  int i;
  for(i=0;i<BATH_MT99ART_NOps;i++) {
    gsl_matrix_free(BATH_MT99ARTBathFunc[i].S);
  }
  BATH_MT99ART_NOps=0;

  return 0;
}

/*
 * $Log$
 * Revision 1.4  2007/03/15 08:20:50  platin
 *   - add direct support for MT99 artificial bath.
 *
 * Revision 1.3  2007/03/15 02:36:26  platin
 *
 *   - add multi-mode brownian oscillator model.
 *
 * Revision 1.2  2006/12/11 19:59:18  platin
 *
 *   - a bug in evaluating Sn in TES basis is fixed.
 *
 * Revision 1.1  2006/11/02 19:29:31  platin
 *
 *   - implement the bathtype interface.
 *
 *
 */
