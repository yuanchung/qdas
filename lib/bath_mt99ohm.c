 
/***************************************************
 * bath_mt99ohm.c
 *
 * Bath functions that implement Meier and Tannor's
 * parameterization of Ohmic spectral density with 
 * exponential cutoff.
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
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_blas.h>

#include "qdas.h"
#include "aux.h"
#include "bath_mt99ohm.h"
#include "bathtype_mt99art.h"

// Defined Constants

/* Meier and Tannor's parameterization of the Ohmic
   spectral density.
   J(w)=gamma*w*exp(-w/wc)
       = \sum_{k=1..3} pk*w/{[(w+Omega_k)^2+Gamma_k^2]*[(w-Omega_k)^2+Gamma_k^2]}
   where

   k    pk/(gamma*wc^4)  Gamma_k/wc  Omega_k/wc
   1    12.0677        2.2593      0.2378
   2   -19.9762        5.4377      0.0888 
   3     0.1834        0.8099      0.0482
*/

const double mt99ohm_ohmic_params[3][3] = 
  {{ 12.0677, 2.2593, 0.2378},
   {-19.9762, 5.4377, 0.0888},
   {  0.1834, 0.8099, 0.0482}};

// Global variables

/* Meier and Tannor's parameterization of the Ohmic
   spectral density.
   J(w)=gamma*w*exp(-w/wc)
       = \sum_{k=1..3} pk*w/{[(w+Omega_k)^2+Gamma_k^2]*[(w-Omega_k)^2+Gamma_k^2]}
*/
double bath_mt99ohm_J(double omega, double gamma, double wc)
{
  double pk,Gk,Ok;
  int i;
  double sum;

  sum=0.0;
  for(i=0;i<3;i++) {
    pk=mt99ohm_ohmic_params[i][0]*gamma*wc*wc*wc*wc;
    Gk=mt99ohm_ohmic_params[i][1]*wc;
    Ok=mt99ohm_ohmic_params[i][2]*wc;
    sum=sum+pk*omega/(((omega+Ok)*(omega+Ok)+Gk*Gk)*((omega-Ok)*(omega-Ok)+Gk*Gk));
  }

  return sum;
}

/* J(iw) = i*Ji(w) */
double bath_mt99ohm_Ji(double omega, double gamma, double wc)
{
  double pk,Gk,Ok,dtmp;
  int i;
  double sum;

  sum=0.0;
  for(i=0;i<3;i++) {
    pk=mt99ohm_ohmic_params[i][0]*gamma*wc*wc*wc*wc;
    Gk=mt99ohm_ohmic_params[i][1]*wc;
    Ok=mt99ohm_ohmic_params[i][2]*wc;
    dtmp=Gk*Gk+Ok*Ok-omega*omega;
    //    if(fabs(dtmp) < wc*wc) {
    //    if(fabs(Gk-omega) < wc) {
      // FIXME: abort when near-resonance Matsubara freq causes divergent term!
    //      printf("ERROR: near resonance Matsubara terms (nu_k=%f)!!\n",omega);
    //      exit(EXIT_FAILURE);
    //    }
    sum=sum+pk*omega/(dtmp*dtmp+omega*omega*Ok*Ok);
  }

  return sum;
}

/* initialize f using assigned temperature and parameters */
void bath_mt99ohm_ct_init(mt99art_bath_func *f, 
			  const double beta, const double gamma, const double wc,
			  const gsl_matrix *S)
{
  double pk,Gk,Ok;
  gsl_complex ztmp;
  double max;
  size_t nr,ni;
  int i;

  // assign gamma, beta, wc, S ... etc.
  f->beta=beta;
  f->S=gsl_matrix_alloc(S->size1,S->size1);
  gsl_matrix_memcpy(f->S,S);
  f->Vren=2.0*gamma*wc/M_PI;

#ifdef BATH_DEBUG
  printf("\n");
  printf("Initializing bath correlation function ...\n");
  printf("  Beta = %f cm (%5.2f K)\n",beta,1.4387/beta);
  printf("  Ohmic spectral sensity:\n");
  printf("    Wc = %f\n",wc);
  printf("    Coupling Strength Gamma = %f\n",gamma);
  printf("\n");
#endif

  /* the real part, a(t) */
  // first 6 terms from the Ohmic parameterization
  nr=0;
  for(i=0;i<3;i++) {
    pk=mt99ohm_ohmic_params[i][0]*gamma*wc*wc*wc*wc;
    Gk=mt99ohm_ohmic_params[i][1]*wc;
    Ok=mt99ohm_ohmic_params[i][2]*wc;
    ztmp=gsl_complex_coth(gsl_complex_rect(beta*Ok/2.0,beta*Gk/2.0));
    ztmp=gsl_complex_mul_real(ztmp,pk/Ok/Gk/8.0);
    if(gsl_complex_abs(ztmp) > 1e-6) {
      f->Ar[nr]=ztmp;
      f->Gr[nr]=gsl_complex_rect(-1.0*Gk,Ok);
      nr++;
    }
    ztmp=gsl_complex_coth(gsl_complex_rect(beta*Ok/2.0,-1.0*beta*Gk/2.0));
    ztmp=gsl_complex_mul_real(ztmp,pk/Ok/Gk/8.0);
    if(gsl_complex_abs(ztmp) > 1e-6) {
      f->Ar[nr]=ztmp;
      f->Gr[nr]=gsl_complex_rect(-1.0*Gk,-1.0*Ok);
      nr++;
    }
  }
  // then Matsubara terms, at most 200 terms
  max=-10.0;
  for(i=0;i<200;i++) {
    double nu_k,Jnu_k;
    nu_k=2.0*M_PI/beta*(double)(i+1);
    Jnu_k=bath_mt99ohm_Ji(nu_k,gamma,wc);
    if(fabs(Jnu_k) > max) {
      max=fabs(Jnu_k);
    }
    if (fabs(Jnu_k) < (0.001*max) || fabs(Jnu_k/beta) < 1e-9) {
      // the following terms are too small
      break;
    }
    // accept this term
    // minus sign in Ar comes from an extra i in J(iw) = i*Ji(w)
    f->Ar[nr]=gsl_complex_rect(-2.0*Jnu_k/beta,0.0);
    f->Gr[nr]=gsl_complex_rect(-1.0*nu_k,0.0);
    nr++;
  }
  /* the negative imaginary part, b(t) */
  // totally 6 terms from the Ohmic parameterization
  ni=0;
  for(i=0;i<3;i++) {
    pk=mt99ohm_ohmic_params[i][0]*gamma*wc*wc*wc*wc;
    Gk=mt99ohm_ohmic_params[i][1]*wc;
    Ok=mt99ohm_ohmic_params[i][2]*wc;
    if(fabs(pk/Ok/Gk/8.0) > 1e-6) {
      f->Ai[ni]=gsl_complex_rect(0.0,pk/Ok/Gk/8.0);
      f->Gi[ni]=gsl_complex_rect(-1.0*Gk,Ok);
      ni++;
      f->Ai[ni]=gsl_complex_rect(0.0,-1.0*pk/Ok/Gk/8.0);
      f->Gi[ni]=gsl_complex_rect(-1.0*Gk,-1.0*Ok);
      ni++;
    }
  }

  f->Nr=nr;
  f->Ni=ni;

#ifdef BATH_DEBUG
  printf("  M&T's parameterization:\n");
  printf("\n");
  printf("  (Nr,Ni)=(%lu,%lu)\n",nr,ni);
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
int bath_mt99ohm_init_params(const size_t nsize, const double beta, 
			     const size_t bath_nparams, const double *bath_params)
{
  int i,idx;
  int nterms;
  size_t ndim;

  double gamma[256],wc[256];

  int opbra[256];
  int opket[256];
  gsl_matrix *S;

  ndim=nsize;
  S=gsl_matrix_alloc(ndim,ndim);

  // input in the mod_params array should be in a table form
  // n m gamma wc
  // parse them to the global variables
  nterms=bath_nparams/4; // number of simple system-bath coupling terms
  for(i=0;i<nterms;i++) {
    int na,nb;
    idx=i*4;
    // system OP represented by ket state |n> and bra state <m|, i.e. |n><m|.
    // We make sure n<=m
    na=(int)bath_params[idx];
    nb=(int)bath_params[idx+1];
    opket[i]=((na < nb) ? na : nb);
    opbra[i]=((na < nb) ? nb : na);
    gamma[i]=bath_params[idx+2];
    wc[i]=bath_params[idx+3];
    if(opbra[i]<=0 || opket[i]<=0) {
      printf("Error: can't parse system-bath interaction terms.\n");
      exit(EXIT_FAILURE);
    }
    /* construct the S operator */
    gsl_matrix_set_zero(S);
    gsl_matrix_set(S,opket[i]-1,opbra[i]-1,1.0);
    gsl_matrix_set(S,opbra[i]-1,opket[i]-1,1.0);
    bath_mt99ohm_ct_init(BATH_MT99ARTBathFunc+BATH_MT99ART_NOps+i,beta,gamma[i],wc[i],S);
  }

  /* show what we have just done */
  printf("Use M&T's parameterization of Ohmic spectral function:\n");
  printf("J(w)=gamma*w*exp(-w/wc)\n");
  printf("\n");
  printf("System-Op and J(w) for each site:\n");
  printf("\n");
  printf("  %16s %12s %12s %12s\n","System-Op","gamma","wc (cm^-1)","(Nr,Ni)");
  for(i=0;i<nterms;i++) {
    // print out the parameters...
    if(opbra[i] == opket[i]) {
      // diagonal e-ph couplings
      printf("  |%d><%d|            %12.4f %12.4f (%4lu,%4lu)\n",
	     opket[i],opbra[i],gamma[i],wc[i],
	     BATH_MT99ARTBathFunc[BATH_MT99ART_NOps+i].Nr,BATH_MT99ARTBathFunc[BATH_MT99ART_NOps+i].Ni);
    } else {
      // off-diagonal terms, note that we will automatically 
      // symmetrize the system part in all evaluations later.
      printf("  |%d><%d|+|%d><%d|   %12.4f %12.4f (%4lu,%4lu)\n",
	     opket[i],opbra[i],
	     opbra[i],opket[i],
	     gamma[i],wc[i],
	     BATH_MT99ARTBathFunc[BATH_MT99ART_NOps+i].Nr,BATH_MT99ARTBathFunc[BATH_MT99ART_NOps+i].Ni);
    }
  }
  printf("\n");

  /* add to the number of MT99 bath terms */
  BATH_MT99ART_NOps+=nterms;

  gsl_matrix_free(S);
  return 0;
}


#ifdef MAIN
int main()
{
  // unit in cm^-1
  double beta=0.018685; // 77K
  //  double beta=0.0048281; // 298K
  double gamma=1.0;
  double wc=100.0;
  //  double beta=0.2;
  //  double wc=50.0;

  double Gk,Ok,pk;

  double w,t;
  gsl_complex ct;
  mt99art_bath_func bf;
  int i=0;

  printf("pk      Gk      Ok\n");
  for(i=0;i<3;i++) {
    pk=mt99ohm_ohmic_params[i][0]*gamma*wc*wc*wc*wc;
    Gk=mt99ohm_ohmic_params[i][1]*wc;
    Ok=mt99ohm_ohmic_params[i][2]*wc;
    printf("%f %f %f\n",pk,Gk,Ok);
  }
  
  w=0.0;
  for(i=0;i<3000;i++) {
    printf("w=%f, Jw=%f\n",w,bath_mt99ohm_J(w,gamma,wc));
    w=w+1.0;
  }
  
  bath_mt99ohm_ct_init(&bf,beta,gamma,wc);
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

  for(wc=100.0;wc<250;wc=wc+1.0) {
    bath_mt99ohm_ct_init(&bf,beta,gamma,wc);
    ct=bath_mt99art_ct(&bf,0.0);
    printf("wc=%f, CtZero=%f + %fI\n",wc,GSL_REAL(ct),GSL_IMAG(ct));    
  }


  return 1;

}
#endif

/*
 * $Log$
 * Revision 1.15  2007/03/12 23:08:27  platin
 *   - fixed MT99 bath modules to allow more general bath functions,
 *     now baths not defined by gamma and wc can be used.
 *
 * Revision 1.14  2006/11/02 19:29:31  platin
 *
 *   - implement the bathtype interface.
 *
 * Revision 1.13  2006/11/01 23:41:54  platin
 *
 *   - modify the mt99ohm bath model to fit into the current
 *     tnl-kernel structure.
 *
 * Revision 1.12  2006/08/23 22:30:26  platin
 *
 *   - add support for two-exciton states in codes regarding H and
 *     transition dipoles. The "Assign" format for transition dipole
 *     input can also be used to include effects of excited state absorption.
 *
 *   - basic support for TESLIST keyword.
 *
 * Revision 1.11  2006/08/10 16:20:18  platin
 *
 *   - minor changes.
 *
 * Revision 1.10  2006/08/08 19:10:29  platin
 *
 *   - correct renormalization energy term.
 *
 * Revision 1.9  2006/08/03 02:15:27  platin
 *
 *   - fix the pulse period setup; this should handle the FWHM correctly.
 *   - in the js03art module, save some time when gamma is effectively zero.
 *
 * Revision 1.8  2006/08/01 19:10:51  platin
 *
 *   - include artificial JS03 bath model that uses four terms for
 *     M&T's parameterization of bath correlation functions.
 *
 * Revision 1.7  2006/07/20 20:18:46  platin
 *
 *   - use more Matsubara terms...
 *
 * Revision 1.6  2006/07/20 17:01:30  platin
 *
 *   - minor change on output formats.
 *
 * Revision 1.5  2006/07/15 07:51:55  platin
 *
 *   - more messages.
 *
 * Revision 1.4  2006/07/15 07:47:39  platin
 *
 *   - add support for g(t) and h(t).
 *   - add a PI factor to C(t) return value; i.e. use a definition of C(t)
 *     that is different from M&T.
 *
 * Revision 1.3  2006/07/13 23:37:29  platin
 *
 *   - add keyword TPRINT that adjusts time step of output.
 *   - minor message changes.
 *
 * Revision 1.2  2006/07/11 18:14:41  platin
 *
 *   - avoid bath terms that are too small.
 *
 * Revision 1.1  2006/07/11 17:22:42  platin
 *
 *   - MT99OHM bath model.
 *
 *
 */
