 /***************************************************
 * bath_mboart.c
 *
 * Bath functions that implement Multimode Brownian 
 * Oscillator model, this can be used in Meier and Tannor's 
 * time-nonlocal EOM.
 *
 * By Yuan-Chung Cheng <yuanchung@ntu.edu.tw>
 *
 * Ref: mbo_artificial.mw
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
#include "bath_mboart.h"
#include "bathtype_mt99art.h"

// Defined Constants

// Global variables

/* Over-damped brownian oscillator model,
   spectral density.
   J(w)=2*lambda*w*gamma0/(w^2+gamma0^2)
*/

double bath_mboart_J(double omega,double lambda, double gamma0)
{
  return (2.0*lambda*omega*gamma0/(omega*omega+gamma0*gamma0));
}

/* initialize f using assigned temperature and parameters */
/* C(t) = a(t) - i*b(t)

               /\inf
   a(t) = 1/Pi | dw  J(w)*cos(w*t)*coth(beta*w/2)
               /0

        = lambda*gamma0*cot(beta*gamma0/2)*exp(-gamma0*t)
          + 4*gamma0*lambda/beta*sum(nu_n*exp(-nu_n*t)/(nu_n^2-gamma0^2),n=1..inf)
          where nu_n = 2*pi*n/beta

               /\inf
   b(t) = 1/Pi | dw  J(w)*sin(w*t)
               /0

        = lambda*gamma0*exp(-gamma0*t)

   also note we parameterize
   a(t) = \sum_k alpha_k*exp(gamma_k*t)
   b(t) = \sum_k alpha_k*exp(gamma_k*t)
*/
void bath_mboart_ct_init(mt99art_bath_func *f, 
			 const double beta, const double lambda, const double gamma0,
			 const gsl_matrix *S)
{
  double dtmp;
  gsl_complex vren;
  double max;
  size_t nr,ni;
  int i;

  // assign beta, S ... etc.
  f->beta=beta;
  f->S=gsl_matrix_alloc(S->size1,S->size1);
  gsl_matrix_memcpy(f->S,S);

#ifdef BATH_DEBUG
  printf("\n");
  printf("Initializing bath correlation function ...\n");
  printf("  Beta = %f cm (%5.2f K)\n",beta,1.4387/beta);
  printf("  Overdamped browian oscillator spectral sensity:\n");
  printf("    lambda0 = %f\n",lambda);
  printf("    gamma0  = %f\n",gamma0);
  printf("\n");
#endif

  /* the real part, a(t) */
  // first term from the Browanian oscillator parameterization
  nr=0;
  f->Ar[nr]=gsl_complex_rect(lambda*gamma0/tan(beta*gamma0/2),0.0);
  f->Gr[nr]=gsl_complex_rect(-1.0*gamma0,0.0);
  nr++;
  // then Matsubara terms, at most 200 terms
  max=-10.0;
  for(i=0;i<200;i++) {
    double nu_i,c_i;
    nu_i=2.0*M_PI/beta*(double)(i+1);
    c_i=4.0*gamma0*lambda/beta*nu_i/(nu_i*nu_i-gamma0*gamma0);
    if(fabs(c_i) > max) {
      max=fabs(c_i);
    }
    if (fabs(c_i) < (0.001*max) || fabs(c_i) < 1e-9) {
      // the following terms are too small
      break;
    }
    // accept this term
    // minus sign in Ar comes from an extra i in J(iw) = i*Ji(w)
    f->Ar[nr]=gsl_complex_rect(c_i,0.0);
    f->Gr[nr]=gsl_complex_rect(-1.0*nu_i,0.0);
    nr++;
  }
  /* the negative imaginary part, b(t) */
  ni=0;
  f->Ai[ni]=gsl_complex_rect(-1.0*lambda*gamma0,0.0); // minus sign for the definition in M&T
  f->Gi[ni]=gsl_complex_rect(-1.0*gamma0,0.0);
  ni++;

  f->Nr=nr;
  f->Ni=ni;

  /* for consistency, we compute the renormalization 
     energy from the imaginary parameterization. The results 
     should be close to lambda0 for browanian oscillator */
  /* mu = 2 * \sum_k Ai[k]/Gi[k] (remember the extra minus 
     sign in the definition of Ai[k] here.) */
  vren=gsl_complex_div(f->Ai[0],f->Gi[0]);
  vren=gsl_complex_mul(vren,gsl_complex_rect(2.0,0.0));
  f->Vren=GSL_REAL(vren);

#ifdef BATH_DEBUG
  printf("  M&T's parameterization:\n");
  printf("\n");
  printf("  (Nr,Ni) = (%d,%d)\n",nr,ni);
  printf("\n");
  printf("  Vren = %f\n",f->Vren);
  printf("\n");
  printf("  a(t) = \\sum_k alpha_k*exp(gamma_k*t)\n");
  printf("\n");
  printf("    %30s %30s\n","alpha_k","gamma_k");
  for(i=0;i<nr;i++) {
    printf("    %14f%+14f   %14f%+14f\n",
	   GSL_REAL(f->Ar[i]),GSL_IMAG(f->Ar[i]),GSL_REAL(f->Gr[i]),GSL_IMAG(f->Gr[i]));
  }
  printf("\n");
  printf("\n");
  printf("  b(t) = \\sum_k alpha_k*exp(gamma_k*t)\n");
  printf("\n");
  printf("    %30s %30s\n","alpha_k","gamma_k");
  for(i=0;i<ni;i++) {
    printf("    %14f%+14f   %14f%+14f\n",
	   GSL_REAL(f->Ai[i]),GSL_IMAG(f->Ai[i]),GSL_REAL(f->Gi[i]),GSL_IMAG(f->Gi[i]));
  }
  printf("\n");
  printf("\n");
#endif

  // done!
}

/* interface functions */
int bath_mboart_init_params(const size_t nsize, const double beta, 
			    const size_t bath_nparams, const double *bath_params)
{
  int i,idx;
  int nterms;
  size_t ndim;

  double lambda0[256],gamma0[256];

  int opbra[256];
  int opket[256];
  gsl_matrix *S;

  ndim=nsize;
  S=gsl_matrix_alloc(ndim,ndim);

  // input in the mod_params array should be in a table form
  // n m lambda0 gamma0
  // lambda0 and gamma0 in cm^-1
  // parse them to the global variables

  nterms=bath_nparams/4;   // number of simple system-bath coupling terms;
  for(i=0;i<nterms;i++) {
    int na,nb;
    idx=i*4;
    // system OP represented by ket state |n> and bra state <m|, i.e. |n><m|.
    // We make sure n<=m
    na=(int)bath_params[idx];
    nb=(int)bath_params[idx+1];
    opket[i]=((na < nb) ? na : nb);
    opbra[i]=((na < nb) ? nb : na);
    lambda0[i]=bath_params[idx+2]; 
    gamma0[i]=bath_params[idx+3];
    if(opbra[i]<=0 || opket[i]<=0) {
      printf("Error: can't parse system-bath interaction terms.\n");
      exit(EXIT_FAILURE);
    }
    /* construct the S operator */
    gsl_matrix_set_zero(S);
    gsl_matrix_set(S,opket[i]-1,opbra[i]-1,1.0);
    gsl_matrix_set(S,opbra[i]-1,opket[i]-1,1.0);
    bath_mboart_ct_init(BATH_MT99ARTBathFunc+BATH_MT99ART_NOps+i,beta,lambda0[i],gamma0[i],S);
  }

  /* show what we have just done */
  printf("Use over-damped Brownian oscillator model:\n");
  printf("J(w)=2*lambda0*gamma0*w/(w^2+gamma0^2)\n");
  printf("\n");
  printf("System-Op and J(w) for each site:\n");
  printf("\n");
  printf("  %16s %16s %16s %12s\n","System-Op","lambda0 (cm^-1)","gamma0 (cm^-1)","(Nr,Ni)");
  for(i=0;i<nterms;i++) {
    // print out the parameters...
    if(opbra[i] == opket[i]) {
      // diagonal e-ph couplings
      printf("  |%d><%d|            %16.4f %16.4f (%4d,%4d)\n",
	     opket[i],opbra[i],lambda0[i],gamma0[i],
	     BATH_MT99ARTBathFunc[BATH_MT99ART_NOps+i].Nr,BATH_MT99ARTBathFunc[BATH_MT99ART_NOps+i].Ni);
    } else {
      // off-diagonal terms, note that we will automatically 
      // symmetrize the system part in all evaluations later.
      printf("  |%d><%d|+|%d><%d|   %16.4f %16.4f (%4d,%4d)\n",
	     opket[i],opbra[i],
	     opbra[i],opket[i],
	     lambda0[i],gamma0[i],
	     BATH_MT99ARTBathFunc[BATH_MT99ART_NOps+i].Nr,BATH_MT99ARTBathFunc[BATH_MT99ART_NOps+i].Ni);
    }
  }
  printf("\n");

  /* add to the number of MT99 bath terms */
  BATH_MT99ART_NOps+=nterms;

  gsl_matrix_free(S);
  return 0;
}

int bath_mboart_free_params()
{
  int i;
  for(i=0;i<BATH_MT99ART_NOps;i++) {
    gsl_matrix_free(BATH_MT99ARTBathFunc[i].S);
  }
  return 0;
}

#ifdef MAIN
int main()
{
  // unit in cm^-1
  double beta=0.018685; // 77K
  //  double beta=0.0048281; // 298K
  double lambda0=80.0;
  double gamma0=100.0;
  //  double beta=0.2;
  //  double wc=50.0;

  double w,t;
  gsl_complex ct;
  mt99art_bath_func bf;
  gsl_matrix *S=gsl_matrix_alloc(2,2);
  int i=0;

  bath_mboart_ct_init(&bf,beta,lambda0,gamma0,S);
  t=0.0;
  for(i=0;i<2000;i++) {
    ct=bath_mt99art_ct(&bf,t);
    printf("t=%f, Ct=%f + %fI\n",t,GSL_REAL(ct),GSL_IMAG(ct));
    printf("t=%f, Gt=%f + %fI\n",t,
	   bath_mt99art_gt_r(&bf,t),bath_mt99art_gt_i(&bf,t));
    printf("t=%f, Ht=%f + %fI\n",t,
	   bath_mt99art_ht_r(&bf,t),bath_mt99art_ht_i(&bf,t));

    t=t+0.0001;
  }

  gsl_matrix_free(S);
  return 1;

}
#endif

/*
 * $Log$
 * Revision 1.2  2007/03/15 08:20:50  platin
 *   - add direct support for MT99 artificial bath.
 *
 * Revision 1.1  2007/03/15 02:36:26  platin
 *
 *   - add multi-mode brownian oscillator model.
 *
 *
 */
