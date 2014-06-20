 
/***************************************************
 * bath_js03art.c
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
#include "bath_js03art.h"
#include "bathtype_mt99art.h"

// Defined Constants

#define NTERMS_MAX (256)

/* Parameterization of the JS03 BChl correlation functions that's
   suitable for Meier and Tannor's TNL calculations
   J(w)=gamma*(0.5*w+0.58*w^2/wc)*exp(-w/wc)
   C_r(t) = 1/Pi * \int_0^\infy dw J(w)*coth(beta*w/2)*cos(w*t),
   C_i(t) = 1/Pi * \int_0^\infy dw J(w)*sin(w*t).
          = a(t) - ib(t)
   where

   a(t) = \sum_k alpha_k*exp(gamma_k*t)
   b(t) = \sum_k alpha_k*exp(gamma_k*t)

*/

#define NTERMS_Wc100T77K_R (5)
const double js03art_At_wc100T77K_params[NTERMS_Wc100T77K_R][3] =  // @ 77K and wc=100 cm^-1
  { /* alpha_k,   Re(gamma_k), Im(gamma_k) */
    {  7180.60,  -301.422,  153.682},
    {  7180.60,  -301.422, -153.682},
    {-3855.965,  -545.696,  384.234},
    {-3855.965,  -545.696, -384.234},
    {  -79.676,  -65.9197,    0.0}};

#define NTERMS_Wc100T77K_I (5)
const double js03art_Bt_wc100T77K_params[NTERMS_Wc100T77K_I][3] =  // @ 77K and wc=100 cm^-1
  { /* alpha_k,   Re(gamma_k), Im(gamma_k) */
    {  -6015.01,   -159.622,  36.8461},
    {  -6015.01,   -159.622, -36.8461},
    {  3840.285,   -337.215,  337.13},
    {  3840.285,   -337.215, -337.13},
    {  4349.45,    -168.166,    0.0}};

// @77K wc=120, updated 08092006
#define NTERMS_Wc120T77K_R (5)
const double js03art_At_wc120T77K_params[NTERMS_Wc120T77K_R][3] =  // @ 77K and wc=120 cm^-1
  { /* alpha_k,   Re(gamma_k), Im(gamma_k) */
    {  10072.45,  -348.156,  187.649},
    {  10072.45,  -348.156, -187.649},
    {  -5063.8 ,  -658.184,  474.8  },
    {  -5063.8 ,  -658.184, -474.8  },
    {  -1060.14,  -140.636,    0.0 }};
//    {  -0.02717793604, -0.001,0.0}};

#define NTERMS_Wc120T77K_I (5)
const double js03art_Bt_wc120T77K_params[NTERMS_Wc120T77K_I][3] =  // @ 77K and wc=120 cm^-1
  { /* alpha_k,   Re(gamma_k), Im(gamma_k) */
    { -8148.815,  -204.008,  42.9763},
    { -8148.815,  -204.008, -42.9763},
    {   5709.15,  -404.94 ,  397.681 },
    {   5709.15,  -404.94 , -397.681 },
    {   4879.33,  -234.575,    0.0 }};
    //    {   41.33817042, -1.0, 0.0}};

// Global variables

/* initialize f using assigned temperature and parameters */
void bath_js03art_ct_init(mt99art_bath_func *f, 
			  const double beta, const double gamma, const double wc,
			  const gsl_matrix *S)
{
  int i,j;
  int nterms_r,nterms_i;

  double bath_At_params[NTERMS_MAX][3],bath_Bt_params[NTERMS_MAX][3];

  gsl_complex vren;

  /* FIXME: implement auto-fitting functions!! 
     at this point, we only support hand-fitted C(t) */
  // only accept 77+-5K and wc=100+-2 or 120+-2
  if(fabs(1.4387/beta-77.0)>5.0) {
    printf("Sorry, parameterizations at T=%.2fK not supported yet!!\n",1.4387/beta);
    exit(EXIT_FAILURE);
  }
  if(fabs(wc-100.0)<2.0) {
    nterms_r=NTERMS_Wc100T77K_R;
    for(i=0;i<nterms_r;i++) {
      for(j=0;j<nterms_r;j++) {
	bath_At_params[i][j]=js03art_At_wc100T77K_params[i][j];
      }
    }
    nterms_i=NTERMS_Wc100T77K_I;
    for(i=0;i<nterms_i;i++) {
      for(j=0;j<nterms_i;j++) {
	bath_Bt_params[i][j]=js03art_Bt_wc100T77K_params[i][j];
      }
    }
  } else if(fabs(wc-120.0)<2.0) {
    nterms_r=NTERMS_Wc120T77K_R;
    for(i=0;i<nterms_r;i++) {
      for(j=0;j<nterms_r;j++) {
	bath_At_params[i][j]=js03art_At_wc120T77K_params[i][j];
      }
    }
    nterms_i=NTERMS_Wc120T77K_I;
    for(i=0;i<nterms_i;i++) {
      for(j=0;j<nterms_i;j++) {
	bath_Bt_params[i][j]=js03art_Bt_wc120T77K_params[i][j];
      }
    }
  } else {
    printf("Sorry, parameterizations at wc=%.2fcm-1 not supported yet!!\n",wc);
    exit(EXIT_FAILURE);
  }

  // assign gamma, beta, wc, S
  f->beta=beta;
  f->S=gsl_matrix_alloc(S->size1,S->size1);
  gsl_matrix_memcpy(f->S,S);

#ifdef BATH_DEBUG
  printf("\n");
  printf("Initializing bath correlation function ...\n");
  printf("  Beta = %f cm (%5.2f K)\n",beta,1.4387/beta);
  printf("  JS03 BChl spectral sensity:\n");
  printf("    Wc = %f\n",wc);
  printf("    Coupling Strength Gamma = %f\n",gamma);
  printf("\n");
#endif

  /* the real part, a(t) */
  // totally nterms terms from the parameterization
  for(i=0;i<nterms_r;i++) {
    f->Ar[i]=gsl_complex_rect(gamma*bath_At_params[i][0],0.0);
    f->Gr[i]=gsl_complex_rect(bath_At_params[i][1],bath_At_params[i][2]);
  }
  /* the negative imaginary part, b(t) */
  for(i=0;i<nterms_i;i++) {
    f->Ai[i]=gsl_complex_rect(gamma*bath_Bt_params[i][0],0.0);
    f->Gi[i]=gsl_complex_rect(bath_Bt_params[i][1],bath_Bt_params[i][2]);
  }

  if(gamma>1e-12) {
    f->Nr=nterms_r;
    f->Ni=nterms_i;
  } else {
    f->Nr=0;
    f->Ni=0;
  }

  /* for consistency, we compute the renormalization 
     energy from the imaginary parameterization. The results 
     should be close to mu=2.16*gamma*wc/M_PI for JS03 spectral density */
  /* mu = 2 * \sum_k Ai[k]/Gi[k] (remember the extra minus 
     sign in the definition of Ai[k] here.) */
  vren=gsl_complex_rect(0.0,0.0);
  for(i=0;i<nterms_i;i++) {
    vren=gsl_complex_add(vren,gsl_complex_div(f->Ai[i],f->Gi[i]));
  }
  vren=gsl_complex_mul(vren,gsl_complex_rect(2.0,0.0));
  f->Vren=GSL_REAL(vren);

#ifdef BATH_DEBUG
  printf("  M&T's parameterization:\n");
  printf("\n");
  printf("Vren = %f\n",f->Vren);
  printf("\n");
  printf("  a(t) = \\sum_k alpha_k*exp(gamma_k*t)\n");
  printf("\n");
  printf("    %30s %30s\n","alpha_k","gamma_k");
  for(i=0;i<f->Nr;i++) {
    printf("    %14f%+14f   %14f%+14f\n",
	   GSL_REAL(f->Ar[i]),GSL_IMAG(f->Ar[i]),GSL_REAL(f->Gr[i]),GSL_IMAG(f->Gr[i]));
  }
  printf("\n");
  printf("\n");
  printf("  b(t) = \\sum_k alpha_k*exp(gamma_k*t)\n");
  printf("\n");
  printf("    %30s %30s\n","alpha_k","gamma_k");
  for(i=0;i<f->Ni;i++) {
    printf("    %14f%+14f   %14f%+14f\n",
	   GSL_REAL(f->Ai[i]),GSL_IMAG(f->Ai[i]),GSL_REAL(f->Gi[i]),GSL_IMAG(f->Gi[i]));
  }
  printf("\n");
  printf("\n");
#endif

  // done!
}

/* interface functions */
int bath_js03art_init_params(const size_t nsize, const double beta, 
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

  nterms=bath_nparams/4; // number of simple system-bath coupling terms;
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
    bath_js03art_ct_init(BATH_MT99ARTBathFunc+BATH_MT99ART_NOps+i,beta,gamma[i],wc[i],S);
  }

  /* show what we have just done */
  printf("Use artificial C(t) parameterization of JS03 BChl J(w):\n");
  printf("J(w)=gamma*(0.5*w+0.58*w^2/wc)*exp(-w/wc)\n");
  printf("\n");
  printf("System-Op and J(w) for each site:\n");
  printf("\n");
  printf("  %16s %12s %12s %12s\n","System-Op","gamma","wc (cm^-1)","(Nr,Ni)");
  for(i=0;i<nterms;i++) {
    // print out the parameters...
    if(opbra[i] == opket[i]) {
      // diagonal e-ph couplings
      printf("  |%d><%d|            %12.4f %12.4f (%4d,%4d)\n",
	     opket[i],opbra[i],gamma[i],wc[i],
	     BATH_MT99ARTBathFunc[BATH_MT99ART_NOps+i].Nr,BATH_MT99ARTBathFunc[BATH_MT99ART_NOps+i].Ni);
    } else {
      // off-diagonal terms, note that we will automatically 
      // symmetrize the system part in all evaluations later.
      printf("  |%d><%d|+|%d><%d|   %12.4f %12.4f (%4d,%4d)\n",
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
  //double beta=0.0048281; // 298K
  double gamma=1.0;
  double wc=100.0;

  double w,t;
  gsl_complex ct;
  mt99art_bath_func bf;
  gsl_matrix *S;
S=gsl_matrix_alloc(2,2);
  int i=0;

  bath_js03art_ct_init(&bf,beta,gamma,wc,S);
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


  return 1;

}
#endif

/*
 * $Log$
 * Revision 1.13  2007/03/12 23:08:27  platin
 *   - fixed MT99 bath modules to allow more general bath functions,
 *     now baths not defined by gamma and wc can be used.
 *
 * Revision 1.12  2006/11/08 06:23:50  platin
 *
 *   - minor changes and adding a artificial ohmic bath module that
 *     uses more efficient artificial fit to Ohmic correlation
 *     functions.
 *
 * Revision 1.11  2006/11/02 19:29:31  platin
 *
 *   - implement the bathtype interface.
 *
 * Revision 1.10  2006/11/01 23:41:54  platin
 *
 *   - modify the mt99ohm bath model to fit into the current
 *     tnl-kernel structure.
 *
 * Revision 1.9  2006/10/27 23:02:59  platin
 *
 *   - commit modifications up to 10/27/2006.
 *
 * Revision 1.8  2006/09/01 17:53:24  platin
 *
 *   - commit improved parameters for js03art bath @ wc=100 and 77K.
 *
 * Revision 1.7  2006/08/26 05:01:17  platin
 *
 *   - include primitive two-exciton state support for dynamics and
 *     optical transitions.
 *
 * Revision 1.6  2006/08/15 23:05:06  platin
 *
 *   - remove ad-hoc fix in tnl-* dynamics; now both programs correctly
 *     handle renormalization term (Hren) and eliminate erroneous dynamics,
 *     such as nonpositive dynamics and non-zero long time coherence terms...
 *
 * Revision 1.5  2006/08/11 06:01:03  platin
 *
 *   - add a smal term in the parameterization for imaginary part
 *     of the TCF to preserve positivity.
 *
 * Revision 1.4  2006/08/10 16:22:51  platin
 *
 *   - reparameterized the 77Kwc120 C(t), a new exponential term is added
 *     to provide a better fit. This params set satisfies detail balance up
 *     to about DeltaE=800 cm^-1
 *
 * Revision 1.3  2006/08/08 19:10:29  platin
 *
 *   - correct renormalization energy term.
 *
 * Revision 1.2  2006/08/03 02:15:27  platin
 *
 *   - fix the pulse period setup; this should handle the FWHM correctly.
 *   - in the js03art module, save some time when gamma is effectively zero.
 *
 * Revision 1.1  2006/08/01 19:10:51  platin
 *
 *   - include artificial JS03 bath model that uses four terms for
 *     M&T's parameterization of bath correlation functions.
 *
 *
 */
