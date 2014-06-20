 /***************************************************
 * bath_mt99art.c
 *
 * Bath functions that directly implement Meier and Tannor's
 * artificicial bath ansatz. This model actually
 * interpolate between overdamped oscillator to coherent
 * vibration bath.
 *
 * By Yuan-Chung Cheng <yuanchung@ntu.edu.tw>
 *
 * Ref: mt99_artificial.mw
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
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_blas.h>

#include "qdas.h"
#include "aux.h"
#include "bath_mt99art.h"
#include "bathtype_mt99art.h"

// Defined Constants

#define MATSUBARA_REL_CUTOFF (1e-3)
#define MATSUBARA_ABS_CUTOFF (1e-9)

// Global variables

/* Over-damped brownian oscillator model,
   spectral density.
   J(w)=4*lambda0*Gamma0^3*w/((w+Omega0)^2+Gamma0^2)/((w-Omega0)^2+Gamma0^2)
   where lambda0 is related to the reorganization energy, 
         Omega0 is the characteristic frequency, and
         Gamma0 is the decay rate of bath correlation.

   the reorganization energy is
   
                         2
                   Gamma0
   lambda0 * -----------------------
              /      2         2\ 
              \Omega0  + Gamma0 / 

*/

double bath_mt99art_Jw(double w,double lambda0,double omega0,double gamma0)
{
  double p,Op,Om;
  p=4.0*lambda0*gamma0*gamma0*gamma0;
  Op=w+omega0;
  Om=w-omega0;
  return (p*w/(Op*Op+gamma0*gamma0)/(Om*Om+gamma0*gamma0));
}

// complex version of Jw; used in the evaluation of Matsubara terms
gsl_complex bath_mt99art_Jw_complex(gsl_complex w,double lambda0,double omega0,double gamma0)
{
  gsl_complex ret;
  gsl_complex p,Op,Om;
  gsl_complex ztmp;

  p=gsl_complex_rect(4.0*lambda0*gamma0*gamma0*gamma0,0.0);
  Op=gsl_complex_add_real(w,omega0);
  Om=gsl_complex_sub_real(w,omega0);

  ret=gsl_complex_mul(p,w);
  ztmp=gsl_complex_add_real(gsl_complex_mul(Op,Op),gamma0*gamma0);
  ret=gsl_complex_div(ret,ztmp);
  ztmp=gsl_complex_add_real(gsl_complex_mul(Om,Om),gamma0*gamma0);
  ret=gsl_complex_div(ret,ztmp);

  return (ret);
}

/* initialize f using assigned temperature and parameters */
/* C(t) = a(t) - i*b(t)

               /\inf
   a(t) = 1/Pi | dw  J(w)*cos(w*t)*coth(beta*w/2)
               /0

        =   lambda0*gamma0^2/2/omega0*coth(beta*Om/2)*exp(-I*Om*t)
          + lambda0*gamma0^2/2/omega0*coth(beta*Op/2)*exp(I*Op*t)
          + sum(2*I/beta*Jw(I*nu_n)*exp(-nu_n*t),n=1..NN);
          
          where nu_n = 2*pi*n/beta

               /\inf
   b(t) = 1/Pi | dw  J(w)*sin(w*t)
               /0

        = I*lambda0*gamma0^2/2/omega0*(exp(-I*Om*t)-exp(I*Op*t))

	where Op=omega0+I*gamma0, and Om=omega0-I*gamma0.

   also note we parameterize
   a(t) = \sum_k alpha_k*exp(gamma_k*t)
   b(t) = \sum_k alpha_k*exp(gamma_k*t)
*/
void bath_mt99art_ct_init(mt99art_bath_func *f, 
			  const double beta, 
			  const double lambda0, const double omega0,const double gamma0,
			  const gsl_matrix *S)
{
  double dtmp;
  gsl_complex Op,Om,vren;
  double max;
  size_t nr,ni;
  int i;

  // assign beta, S ... etc.
  f->beta=beta;
  f->S=gsl_matrix_alloc(S->size1,S->size1);
  gsl_matrix_memcpy(f->S,S);

  Op=gsl_complex_rect(omega0,gamma0);
  Om=gsl_complex_rect(omega0,-1.0*gamma0);

#ifdef BATH_DEBUG
  printf("\n");
  printf("Initializing bath correlation function ...\n");
  printf("  Beta = %f cm (%5.2f K)\n",beta,1.4387/beta);
  printf("  MT99 artificial bath model:\n");
  printf("    lambda0 = %f\n",lambda0);
  printf("    omega0  = %f\n",omega0);
  printf("    gamma0  = %f\n",gamma0);
  printf("\n");
#endif

  /* the real part, a(t) */
  // first two terms from the MT99 parameterization
  nr=0;
  // the Om term
  f->Ar[nr]=gsl_complex_mul_real(gsl_complex_coth(gsl_complex_rect(beta*omega0/2.0,-1.0*beta*gamma0/2.0)),
				 lambda0*gamma0*gamma0/2.0/omega0);
  f->Gr[nr]=gsl_complex_rect(-1.0*gamma0,-1.0*omega0);
  nr++;
  // the Op term
  f->Ar[nr]=gsl_complex_mul_real(gsl_complex_coth(gsl_complex_rect(beta*omega0/2.0,beta*gamma0/2.0)),
				 lambda0*gamma0*gamma0/2.0/omega0);
  f->Gr[nr]=gsl_complex_rect(-1.0*gamma0,omega0);
  nr++;
  // then Matsubara terms, at most 200 terms
  max=fabs(GSL_REAL(f->Ar[0]));
  for(i=0;i<200;i++) {
    gsl_complex nu_i,c_i,Jw_i;
    nu_i=gsl_complex_rect(0.0,2.0*M_PI/beta*(double)(i+1));
    Jw_i=bath_mt99art_Jw_complex(nu_i,lambda0,omega0,gamma0);
    c_i=gsl_complex_mul(gsl_complex_rect(0.0,2.0/beta),Jw_i);
    if(gsl_complex_abs(c_i) > max) {
      max=gsl_complex_abs(c_i);
    }
    if (gsl_complex_abs(c_i) < (MATSUBARA_REL_CUTOFF*max) || gsl_complex_abs(c_i) < MATSUBARA_ABS_CUTOFF) {
      // the following terms are too small
      break;
    }
    // accept this term
    f->Ar[nr]=gsl_complex_rect(GSL_REAL(c_i),0.0);
    f->Gr[nr]=gsl_complex_rect(-2.0*M_PI/beta*(double)(i+1),0.0);
    nr++;
  }
  /* the negative imaginary part, b(t); 
   notice the extra -1 factor */
  ni=0;
  f->Ai[ni]=gsl_complex_rect(0.0,-1.0*lambda0*gamma0*gamma0/2.0/omega0); // Om term
  f->Gi[ni]=gsl_complex_rect(-1.0*gamma0,-1.0*omega0);
  ni++;
  f->Ai[ni]=gsl_complex_rect(0.0,lambda0*gamma0*gamma0/2.0/omega0); // Op term
  f->Gi[ni]=gsl_complex_rect(-1.0*gamma0,omega0);
  ni++;

  f->Nr=nr;
  f->Ni=ni;

  /* for consistency, we compute the renormalization 
     energy from the imaginary parameterization. The results 
     should be close to 2 * lambda0 * gamma0^2/(omega0^2+gamma0^2) for 
     MT99 bath.
     mu = 2 * \sum_k Ai[k]/Gi[k] (remember the extra minus 
     sign in the definition of Ai[k] here.) */
  vren=gsl_complex_rect(0.0,0.0);
  for(i=0;i<ni;i++) {
    vren=gsl_complex_add(vren,gsl_complex_div(f->Ai[i],f->Gi[i]));
  }
  vren=gsl_complex_mul(vren,gsl_complex_rect(2.0,0.0));
  //  f->Vren=GSL_REAL(vren);
  f->Vren=GSL_REAL(vren);
  //f->Vren=0.0;

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
int bath_mt99art_init_params(const size_t nsize, const double beta, 
			     const size_t bath_nparams, const double *bath_params)
{
  int i,idx;
  int nterms;
  size_t ndim;

  double lambda0[256],omega0[256],gamma0[256];
  int opbra[256];
  int opket[256];
  gsl_matrix *S;

  ndim=nsize;
  S=gsl_matrix_alloc(ndim,ndim);

  // input in the mod_params array should be in a table form
  // n m lambda0 omega0 gamma0
  // lambda0, omega0, and gamma0 all in cm^-1
  // parse them to the global variables


  nterms=bath_nparams/5;   // number of simple system-bath coupling terms;
  for(i=0;i<nterms;i++) {
    int na,nb;
    idx=i*5;
    // system OP represented by ket state |n> and bra state <m|, i.e. |n><m|.
    // We make sure n<=m
    na=(int)bath_params[idx];
    nb=(int)bath_params[idx+1];
    opket[i]=((na < nb) ? na : nb);
    opbra[i]=((na < nb) ? nb : na);
    lambda0[i]=bath_params[idx+2]; 
    omega0[i]=bath_params[idx+3]; 
    gamma0[i]=bath_params[idx+4];
    if(opbra[i]<=0 || opket[i]<=0) {
      printf("Error: can't parse system-bath interaction terms.\n");
      exit(EXIT_FAILURE);
    }
    /* construct the S operator and initialize the bath C(t) */
    gsl_matrix_set_zero(S);
    gsl_matrix_set(S,opket[i]-1,opbra[i]-1,1.0);
    gsl_matrix_set(S,opbra[i]-1,opket[i]-1,1.0);
    bath_mt99art_ct_init(BATH_MT99ARTBathFunc+i+BATH_MT99ART_NOps,beta,lambda0[i],omega0[i],gamma0[i],S);
  }

  /* show what we have just done */
  printf("Use Meier and Tannor's artificial bath model:\n");
  printf("J(w)=4*lambda0*gamma0^3*w/((w+omega0)^2+gamma0^2)/((w-omega0)^2+gamma0^2)\n");
  printf("\n");
  printf("System-Op and J(w) for each site (in cm^-1):\n");
  printf("\n");
  printf("  %16s %14s %14s %14s %12s\n","System-Op","lambda0","omega0","gamma0","(Nr,Ni)");
  for(i=0;i<nterms;i++) {
    // print out the parameters...
    if(opbra[i] == opket[i]) {
      // diagonal e-ph couplings
      printf("  |%d><%d|            %14.4f %14.4f %14.4f (%4d,%4d)\n",
	     opket[i],opbra[i],lambda0[i],omega0[i],gamma0[i],
	     BATH_MT99ARTBathFunc[BATH_MT99ART_NOps+i].Nr,BATH_MT99ARTBathFunc[BATH_MT99ART_NOps+i].Ni);
    } else {
      // off-diagonal terms, note that we will automatically 
      // symmetrize the system part in all evaluations later.
      printf("  |%d><%d|+|%d><%d|   %14.4f %14.4f %14.4f (%4d,%4d)\n",
	     opket[i],opbra[i],
	     opbra[i],opket[i],
	     lambda0[i],omega0[i],gamma0[i],
	     BATH_MT99ARTBathFunc[BATH_MT99ART_NOps+i].Nr,BATH_MT99ARTBathFunc[BATH_MT99ART_NOps+i].Ni);
    }
  }
  printf("\n");

  /* now treat the non-local bath terms */
  /* FIXME: no support yet!
  for(i=0;i<keys->nlbath_nterms;i++) {
    double nlgamma,nlwc;
    gamma=keys->nlbath[i].gamma;
    wc=keys->nlbath[i].wc;
    bath_mt99art_ct_init(BATH_MT99ARTBathFunc+nterms,
			 keys->beta,gamma,wc,keys->nlbath[i].S);
    nterms++;
  }
  */

  /* save the number of terms */
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
  // three MT99 modes
  double lambda1=500.0;
  double omega1=250.0;
  double gamma1=20.0;
  double lambda2=150.0;
  double omega2=10.0;
  double gamma2=20.0;
  double lambda3=0.0;
  double omega3=3.0;
  double gamma3=10.0;
  //  double beta=0.2;
  //  double wc=50.0;

  double w,t;
  gsl_complex ct;
  mt99art_bath_func bf1,bf2,bf3;
  gsl_matrix *S=gsl_matrix_alloc(2,2);
  int i=0;

  printf("lambda1 = %f\n",lambda1);
  printf("omega1  = %f\n",omega1);
  printf("gamma1  = %f\n",gamma1);
  printf("lambda2 = %f\n",lambda2);
  printf("omega2  = %f\n",omega2);
  printf("gamma2  = %f\n",gamma2);
  printf("lambda3 = %f\n",lambda3);
  printf("omega3  = %f\n",omega3);
  printf("gamma3  = %f\n",gamma3);

  for(w=0;w<1000.0;w=w+1.0) {
    printf("w=%16.6f, Jw=%16.6f\n",w,
	   bath_mt99art_Jw(w,lambda1,omega1,gamma1)
	   +bath_mt99art_Jw(w,lambda2,omega2,gamma2)
	   +bath_mt99art_Jw(w,lambda3,omega3,gamma3));
  }

  bath_mt99art_ct_init(&bf1,beta,lambda1,omega1,gamma1,S);
  bath_mt99art_ct_init(&bf2,beta,lambda2,omega2,gamma2,S);
  bath_mt99art_ct_init(&bf3,beta,lambda3,omega3,gamma3,S);
  t=0.0;
  for(i=0;i<2000;i++) {
    ct=gsl_complex_add(bath_mt99art_ct(&bf1,t),bath_mt99art_ct(&bf2,t));
    ct=gsl_complex_add(ct,bath_mt99art_ct(&bf3,t));
    printf("t=%f, Ct=%f + %fI\n",t,GSL_REAL(ct)/M_PI,GSL_IMAG(ct)/M_PI);
    printf("t=%f, Ht=%f + %fI\n",t,
	   bath_mt99art_ht_r(&bf1,t)/M_PI+bath_mt99art_ht_r(&bf2,t)/M_PI+bath_mt99art_ht_r(&bf3,t)/M_PI,
	   bath_mt99art_ht_i(&bf1,t)/M_PI+bath_mt99art_ht_i(&bf2,t)/M_PI+bath_mt99art_ht_i(&bf3,t)/M_PI);
    printf("t=%f, Gt=%f + %fI\n",t,
	   bath_mt99art_gt_r(&bf1,t)/M_PI+bath_mt99art_gt_r(&bf2,t)/M_PI+bath_mt99art_gt_r(&bf3,t)/M_PI,
	   bath_mt99art_gt_i(&bf1,t)/M_PI+bath_mt99art_gt_i(&bf2,t)/M_PI+bath_mt99art_gt_i(&bf3,t)/M_PI);
    t=t+0.0001;
  }

  gsl_matrix_free(S);
  return 1;

}
#endif

/*
 * $Log$
 * Revision 1.4  2007/08/09 23:50:18  platin
 *   - implemented site localized extra population dynamics in TNL.
 *   - implemented the DIAGONALDYNAMICS method for population-only
 *     incoherent dynamics.
 *
 * Revision 1.3  2007-07-27 22:50:57  platin
 *
 *   - minor correction; test two model spectral density in bath_mt99art.c.
 *
 * Revision 1.2  2007/03/26 19:11:22  platin
 *
 *   - minor changes.
 *
 * Revision 1.1  2007/03/15 08:20:50  platin
 *
 *   - add direct support for MT99 artificial bath.
 *
 *
 */
