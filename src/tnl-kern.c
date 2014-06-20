/***************************************************
 * tnl-kern.c
 * 
 * By Yuan-Chung Cheng <yuanchung@ntu.edu.tw>
 *
 * Simulating three-pulse photon-echo signal using
 * Domcke's density-matrix based method.
 * Domcke's method in JCP 123, 164112 (2005).
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
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_sf_bessel.h>

#include "tnl-kern.h"
#include "params.h"
#include "aux.h"
#include "dmop.h"
#include "bath.h"
#include "bathtype_mt99art.h"

// modules

/* Global constants */

/* global variables */

/* Macros */

/* Display usage information and exit.  */

/* Equations of motion */

/* computate steady-state bath aux density matrix, real part:
   Gr*sr(t) = i[Sn,sigma(t)]
   use mtmp as working space; field-interactions are ignored
   Note that this is steady-state in the interaction picture...
*/
void compute_bath_aux_dm_r_ss(gsl_matrix_complex *sr, 
			      const gsl_matrix_complex *sigma, const gsl_matrix_complex *Sn,
			      const gsl_complex Gr)
{
  /* i[Sn,sigma(t)] commutator term
     = i*Sn*sigma(t) - i*sigma(t)*Sn */
  gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,gsl_complex_rect(0.0,1.0), Sn, sigma,
		 gsl_complex_rect(0.0,0.0), sr); // sr = i*Sn*sigma(t)
  gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,gsl_complex_rect(0.0,-1.0), sigma, Sn,
		 gsl_complex_rect(1.0,0.0), sr); // sr = sr - i*sigma(t)*Sn
  gsl_matrix_complex_scale(sr,gsl_complex_inverse(Gr)); // sr = sr(t)/Gr
  
  // DONE!
}

/* computate steady-state bath aux density matrix, imag part:
   Gi*si(t) = -[Sn,sigma(t)]
   use mtmp as working space; field-interactions are ignored
*/
void compute_bath_aux_dm_i_ss(gsl_matrix_complex *si, 
			      const gsl_matrix_complex *sigma, const gsl_matrix_complex *Ss,
			      const gsl_complex Gi)
{
  /* [Sn - phi,sigma(t)]_+ anticommutator term
     = -[Ss*sigma(t) + sigma(t)*Ss] */
  gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,gsl_complex_rect(-1.0,0.0), Ss, sigma,
		 gsl_complex_rect(0.0,0.0), si); // si = - Ss*sigma(t)
  gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,gsl_complex_rect(-1.0,0.0), sigma, Ss,
		 gsl_complex_rect(1.0,0.0), si); // si = si - sigma(t)*Ss
  
  gsl_matrix_complex_scale(si,gsl_complex_inverse(Gi)); // si = si(t)/Gi
  
  // DONE!
}

/* compute equilibrium bath aux density matrix using SCF equation, real part:
   Gr*sr(t) = i[Sn,sigma(t)] + i[W(t),sr(t)]
 */
void compute_bath_aux_dm_r_eq(gsl_matrix_complex *sr, 
			      const gsl_matrix_complex *sigma, const gsl_matrix_complex *Sn,
			      const gsl_matrix_complex *Ht, gsl_complex Gr)
{
  int i;
  size_t ndim;
  gsl_matrix_complex *Sn_sigmat,*mtmp;

  ndim=sr->size1;

  Sn_sigmat=gsl_matrix_complex_alloc(ndim,ndim);
  mtmp=gsl_matrix_complex_alloc(ndim,ndim);

  /* i[Sn,sigma(t)] commutator term
     = i*Sn*sigma(t) - i*sigma(t)*Sn */
  gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,gsl_complex_rect(0.0,1.0), Sn, sigma,
		 gsl_complex_rect(0.0,0.0), Sn_sigmat); // Sn_sigmat = i*Sn*sigma(t)
  gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,gsl_complex_rect(0.0,-1.0), sigma, Sn,
		 gsl_complex_rect(1.0,0.0), Sn_sigmat); // Sn_sigmat = Sn_sigmat - i*sigma(t)*Sn

  /* start from sr=0 */
  gsl_matrix_complex_set_zero(sr);

  // FIXME: should check convergence */
  for(i=0;i<100;i++) {

    /* i[W(t),sr(t)] = i*W(t)*sr(t) - i*sr(t)*W(t) */
    gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,gsl_complex_rect(0.0,1.0), Ht, sr,
		   gsl_complex_rect(0.0,0.0), mtmp); // mtmp = i*W(t)*sr(t)
    gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,gsl_complex_rect(0.0,-1.0), sr, Ht,
		   gsl_complex_rect(1.0,0.0), mtmp); // mtmp = mtmp - i*sr(t)*W(t)
    //    printf("TEST %d\n",i);
    gsl_matrix_complex_memcpy(sr,Sn_sigmat);
    gsl_matrix_complex_add(sr,mtmp);
    gsl_matrix_complex_scale(sr,gsl_complex_inverse(Gr));
  }

  gsl_matrix_complex_free(Sn_sigmat);
  gsl_matrix_complex_free(mtmp);

#ifdef DEBUGDEBUG
  printf("\n");
  printf("Init sigma = \n");
  gsl_matrix_complex_print(sigma);
  printf("Init Sn = \n");
  gsl_matrix_complex_print(Sn);
  printf("Init Ht = \n");
  gsl_matrix_complex_print(Ht);
  printf("Init Gr = %20.16f, %20.16f\n",GSL_REAL(Gr),GSL_IMAG(Gr));
  printf("Init Gr^-1 = %20.16f, %20.16f\n",GSL_REAL(gsl_complex_inverse(Gr)),GSL_IMAG(gsl_complex_inverse(Gr)));
  printf("Init dm_r = \n");
  gsl_matrix_complex_print(sr);
#endif

  // DONE!
}

/* compute equilibrium bath aux density matrix using SCF equation, imaginary part:
   Gi*si(t) = -[Sn,sigma(t)]_+ + i[W(t),si(t)]
 */
void compute_bath_aux_dm_i_eq(gsl_matrix_complex *si, 
			      const gsl_matrix_complex *sigma, const gsl_matrix_complex *Sn,
			      const gsl_matrix_complex *Ht, gsl_complex Gi)
{
  int i;
  size_t ndim;
  gsl_matrix_complex *Sn_sigmat,*mtmp;

  ndim=si->size1;

  Sn_sigmat=gsl_matrix_complex_alloc(ndim,ndim);
  mtmp=gsl_matrix_complex_alloc(ndim,ndim);

  /* -[Sn,sigma(t)]_+ commutator term
     = -Sn*sigma(t) - sigma(t)*Sn */
  gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,gsl_complex_rect(-1.0,0.0), Sn, sigma,
		 gsl_complex_rect(0.0,0.0), Sn_sigmat); // Sn_sigmat = -Sn*sigma(t)
  gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,gsl_complex_rect(-1.0,0.0), sigma, Sn,
		 gsl_complex_rect(1.0,0.0), Sn_sigmat); // Sn_sigmat = Sn_sigmat - sigma(t)*Sn

  /* start from si=0 */
  gsl_matrix_complex_set_zero(si);

  // FIXME: should check convergence */
  for(i=0;i<100;i++) {

    /* i[W(t),si(t)] = i*W(t)*si(t) - i*si(t)*W(t) */
    gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,gsl_complex_rect(0.0,1.0), Ht, si,
		   gsl_complex_rect(0.0,0.0), mtmp); // mtmp = i*W(t)*si(t)
    gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,gsl_complex_rect(0.0,-1.0), si, Ht,
		   gsl_complex_rect(1.0,0.0), mtmp); // mtmp = mtmp - i*si(t)*W(t)

    gsl_matrix_complex_memcpy(si,Sn_sigmat);
    gsl_matrix_complex_add(si,mtmp);
    gsl_matrix_complex_scale(si,gsl_complex_inverse(Gi));
  }

  gsl_matrix_complex_free(Sn_sigmat);
  gsl_matrix_complex_free(mtmp);

  // DONE!
}

/* compute derivative of a bath aux density matrix, real part:
   d sr(t) /dt = -i[Sn,sigma(t)] - i[W(t),sr(t)] + Gr*sr(t)
   use mtmp as working space
 */
void compute_bath_aux_dm_r_dt(gsl_matrix_complex *dsrdt, const gsl_matrix_complex *sr, 
			      const gsl_matrix_complex *sigma, const gsl_matrix_complex *Sn,
			      const gsl_matrix_complex *Ht, gsl_complex Gr)
{
  gsl_matrix_complex_memcpy(dsrdt,sr);
  gsl_matrix_complex_scale(dsrdt,Gr); // dsrdt = Gr*sr(t)

  /* I don't know why the zsymm and zhemm version is not faster
  gsl_blas_zsymm(CblasLeft, CblasUpper, gsl_complex_rect(0.0,-1.0), Sn, sigma, 
		 gsl_complex_rect(1.0,0.0), dsrdt); // dsrdt = dsrdt - i*Sn*sigma(t)
  gsl_blas_zsymm(CblasRight, CblasUpper, gsl_complex_rect(0.0,1.0), Sn, sigma, 
		 gsl_complex_rect(1.0,0.0), dsrdt); // dsrdt = dsrdt + i*sigma(t)*Sn
  
  gsl_blas_zhemm(CblasLeft, CblasUpper, gsl_complex_rect(0.0,-1.0), Ht, sr, 
		 gsl_complex_rect(1.0,0.0), dsrdt); // dsrdt = dsrdt - i*W(t)*sr(t)
  gsl_blas_zhemm(CblasRight, CblasUpper, gsl_complex_rect(0.0,1.0), Ht, sr, 
		 gsl_complex_rect(1.0,0.0), dsrdt); // dsrdt = dsrdt + i*sr(t)*W(t)
  */

  /* -i[Sn,sigma(t)] commutator term
     = -i*Sn*sigma(t) + i*sigma(t)*Sn ; note that Sn is a symmetric matrix */
  gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,gsl_complex_rect(0.0,-1.0), Sn, sigma, 
		 gsl_complex_rect(1.0,0.0), dsrdt); // dsrdt = dsrdt - i*Sn*sigma(t)
  gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,gsl_complex_rect(0.0,1.0), sigma, Sn,
		 gsl_complex_rect(1.0,0.0), dsrdt); // dsrdt = dsrdt + i*sigma(t)*Sn
  
  /* -i[W(t),sr(t)] = -i*W(t)*sr(t) + i*sr(t)*W(t) */
  gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,gsl_complex_rect(0.0,-1.0), Ht, sr, 
		 gsl_complex_rect(1.0,0.0), dsrdt); // dsrdt = dsrdt - i*W(t)*sr(t)
  gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,gsl_complex_rect(0.0,1.0), sr, Ht, 
		 gsl_complex_rect(1.0,0.0), dsrdt); // dsrdt = dsrdt + i*sr(t)*W(t)

  // DONE!
}

/* compute derivative of a bath aux density matrix, imagnary part:
   d si(t) /dt = [Sn,sigma(t)]_+ - i[W(t),si(t)] + Gi*si(t)
 */
void compute_bath_aux_dm_i_dt(gsl_matrix_complex *dsidt, const gsl_matrix_complex *si, 
			      const gsl_matrix_complex *sigma, const gsl_matrix_complex *Ss,
			      const gsl_matrix_complex *Ht, gsl_complex Gi)
{
  gsl_matrix_complex_memcpy(dsidt,si);
  gsl_matrix_complex_scale(dsidt,Gi); // dsidt = Gi*si(t)

  /* I don't know why the zsymm and zhemm version is not faster
  gsl_blas_zsymm(CblasLeft, CblasUpper, gsl_complex_rect(1.0,0.0), Ss, sigma, 
		 gsl_complex_rect(1.0,0.0), dsidt); // dsidt = dsidt + Ss*sigma(t)
  gsl_blas_zsymm(CblasRight, CblasUpper, gsl_complex_rect(1.0,0.0), Ss, sigma, 
		 gsl_complex_rect(1.0,0.0), dsidt); // dsidt = dsidt + sigma(t)*Ss
  gsl_blas_zhemm(CblasLeft, CblasUpper, gsl_complex_rect(0.0,-1.0), Ht, si, 
		 gsl_complex_rect(1.0,0.0), dsidt); // dsidt = dsidt - i*W(t)*si(t)
  gsl_blas_zhemm(CblasRight, CblasUpper, gsl_complex_rect(0.0,1.0), Ht, si, 
		 gsl_complex_rect(1.0,0.0), dsidt); // dsidt = dsidt + i*si(t)*W(t)
  */
  /* [Sn - phi,sigma(t)]_+ anticommutator term
     = Ss*sigma(t) + sigma(t)*Ss */
  gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,gsl_complex_rect(1.0,0.0), Ss, sigma, 
		 gsl_complex_rect(1.0,0.0), dsidt); // dsidt = dsidt + Ss*sigma(t)
  gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,gsl_complex_rect(1.0,0.0), sigma, Ss,
		 gsl_complex_rect(1.0,0.0), dsidt); // dsidt = dsidt + sigma(t)*Ss

  /* -i[W(t),si(t)] = -i*W(t)*si(t) + i*si(t)*W(t) */
  gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,gsl_complex_rect(0.0,-1.0), Ht, si, 
		 gsl_complex_rect(1.0,0.0), dsidt); // dsidt = dsidt - i*W(t)*si(t)
  gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,gsl_complex_rect(0.0,1.0), si, Ht,
		 gsl_complex_rect(1.0,0.0), dsidt); // dsidt = dsidt + i*si(t)*W(t)

  // DONE!
}

/* Compute the derivative of density matrix using given tfs and ops
   Equation of motion is in p. 76 of Notebook 1;  */
void compute_dsigma_dt(gsl_matrix_complex *dsigmadt, tnl_aux_funcs *tfs,
		       const gsl_matrix_complex *Weff_t,const qdas_keys *keys)
{
  size_t ndim,n,i;

  gsl_matrix_complex *mtmp,*sumOp;

  gsl_complex ztmp,dpdt;
  size_t a,src,dest;

  ndim=tfs->sigma->size1;

  mtmp=gsl_matrix_complex_alloc(ndim,ndim);
  sumOp=gsl_matrix_complex_alloc(ndim,ndim);

  gsl_matrix_complex_set_zero(dsigmadt);

  /* bath terms */
  for(n=0;n<tfs->nops;n++) {
    gsl_matrix_complex_set_zero(sumOp);
    for(i=0;i<tfs->bfs[n].Nr;i++) {
      gsl_matrix_complex_memcpy(mtmp,tfs->bfs[n].dm_r[i]);
      gsl_matrix_complex_scale(mtmp,tfs->bfs[n].Ar[i]);
      gsl_matrix_complex_add(sumOp,mtmp);
    }
    for(i=0;i<tfs->bfs[n].Ni;i++) {
      gsl_matrix_complex_memcpy(mtmp,tfs->bfs[n].dm_i[i]);
      gsl_matrix_complex_scale(mtmp,tfs->bfs[n].Ai[i]);
      gsl_matrix_complex_add(sumOp,mtmp);
    }
    gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,gsl_complex_rect(0.0,-1.0), tfs->bfs[n].Sbar, sumOp,
		   gsl_complex_rect(0.0,0.0), mtmp); // mtmp = -i*Sn*sumOp
    gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,gsl_complex_rect(0.0,1.0), sumOp, tfs->bfs[n].Sbar,
		   gsl_complex_rect(1.0,0.0), mtmp); // mtmp = mtmp + i*sumOp*Sn
    gsl_matrix_complex_add(dsigmadt,mtmp);
  }

  if(keys->diagonaldynamics == QDAS_SET_DIAGONAL_DYNAMICS) {
    /* ad-hoc version of the population-only dynamics is requested,
       so we neglect all changes in the coherence density matrix terms. 
       Note that this is not really a secular approximation (which does 
       not exist for TNL dynamics), and only makes scense when the 
       initial density matrix does not have coherence terms. */
    gsl_matrix_complex_memcpy(mtmp,dsigmadt);
    gsl_matrix_complex_set_zero(dsigmadt);
    for(i=0;i<ndim;i++) {
      gsl_matrix_complex_set(dsigmadt,i,i,
			     gsl_matrix_complex_get(mtmp,i,i));
    }
  }

  /* coherent commutator term dsigmadt=-i[Weff_t,sigma] */
  gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,gsl_complex_rect(0.0,-1.0), Weff_t, tfs->sigma,
		 gsl_complex_rect(1.0,0.0), dsigmadt); // dsigmadt = dsigmadt - I*Weff_t*sigma
  gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,gsl_complex_rect(0.0,1.0), tfs->sigma, Weff_t,
		 gsl_complex_rect(1.0,0.0), dsigmadt); //dsigmadt = dsigmadt + I*sigma*Weff_t

  /* extra population transfer terms; using Lindblad theorem */
  if(keys->nptrans>0) {
    double k;
    for(n=0;n<keys->nptrans;n++) {
      /* add k*rho(src,src)*dt to rho(dest,dest), and 
	 sub the same amount from rho(src,src) */
      src=keys->extra_ptrans[n].src;
      dest=keys->extra_ptrans[n].dest;
      k=keys->extra_ptrans[n].k;
      dpdt=gsl_complex_mul_real(gsl_matrix_complex_get(tfs->sigma,src,src),k);
      ztmp=gsl_matrix_complex_get(dsigmadt,src,src);
      gsl_matrix_complex_set(dsigmadt,src,src,gsl_complex_sub(ztmp,dpdt));
      ztmp=gsl_matrix_complex_get(dsigmadt,dest,dest);
      gsl_matrix_complex_set(dsigmadt,dest,dest,gsl_complex_add(ztmp,dpdt));
      // off-diagonal terms
      for(a=0;a<ndim;a++) {
	if(a != src && a != dest) {
	  dpdt=gsl_complex_mul_real(gsl_matrix_complex_get(tfs->sigma,src,a),0.5*k);
	  ztmp=gsl_matrix_complex_get(dsigmadt,src,a);
	  gsl_matrix_complex_set(dsigmadt,src,a,gsl_complex_sub(ztmp,dpdt));
	  dpdt=gsl_complex_mul_real(gsl_matrix_complex_get(tfs->sigma,a,src),0.5*k);
	  ztmp=gsl_matrix_complex_get(dsigmadt,a,src);
	  gsl_matrix_complex_set(dsigmadt,a,src,gsl_complex_sub(ztmp,dpdt));
	}
      }
    }
  }

  /* extra site-representation population transfer terms; using Lindblad theorem */
  if(keys->nsiteptrans>0) {
    double k;
    gsl_matrix_complex *drdt_n;
    drdt_n=gsl_matrix_complex_alloc(ndim,ndim);
    /* drho/dt = \sum_{k}\gamma_{k}\left(A_{k}\rho_{s}A_{k}^{\dagger}-
       \frac{1}{2}A_{k}^{\dagger}A_{k}\rho_{s}-
       \frac{1}{2}\rho_{s}A_{k}^{\dagger}A_{k}\right). */
    for(n=0;n<keys->nsiteptrans;n++) {
      k=keys->extra_siteptrans[n].k;
      /* first term A*rho*A^\dagger */
      gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,gsl_complex_rect(1.0,0.0), 
		     keys->extra_siteptrans[n].Op1, tfs->sigma,
		     gsl_complex_rect(0.0,0.0),mtmp); // mtmp = A*sigma      
      gsl_blas_zgemm(CblasNoTrans,CblasTrans,gsl_complex_rect(1.0,0.0), 
		     mtmp,keys->extra_siteptrans[n].Op1,
		     gsl_complex_rect(0.0,0.0),drdt_n); // drdt_n = mtmp*A^\dagger
      gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,gsl_complex_rect(1.0,0.0), 
		     keys->extra_siteptrans[n].Op2, tfs->sigma,
		     gsl_complex_rect(0.0,0.0),mtmp); // mtmp = A^\dagger*A*sigma      
      gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,gsl_complex_rect(1.0,0.0), 
		     tfs->sigma,keys->extra_siteptrans[n].Op2,
		     gsl_complex_rect(1.0,0.0),mtmp); // mtmp = mtmp + sigma*A^\dagger*A
      gsl_matrix_complex_scale(mtmp,gsl_complex_rect(0.5,0.0));
      gsl_matrix_complex_sub(drdt_n,mtmp); // drdt_n = drdt_n - mtmp
      gsl_matrix_complex_scale(drdt_n,gsl_complex_rect(k,0.0)); // rate constant
      gsl_matrix_complex_add(dsigmadt,drdt_n);
    }
    gsl_matrix_complex_free(drdt_n);
  }

  /* extra nonlinear pairwise transfer terms; also using Lindblad theorem */
  if(keys->nnlptrans>0) {
    double k;
    for(n=0;n<keys->nnlptrans;n++) {
      /* sub k*rho(src,src)*rho(src,src)*dt from rho(src,src), and 
	 add 1/2 that amount to rho(dest,dest) */
      src=keys->extra_nlptrans[n].src;
      dest=keys->extra_nlptrans[n].dest;
      k=keys->extra_nlptrans[n].k;
      // assigned nonlinear population transfer; exciton-exciton anil...
      // scale k by population of the source state; d/dt p(src) = -k*p(src)*p(src)
      k=k*GSL_REAL(gsl_matrix_complex_get(tfs->sigma,src,src));
      dpdt=gsl_complex_mul_real(gsl_matrix_complex_get(tfs->sigma,src,src),k);
      ztmp=gsl_matrix_complex_get(dsigmadt,src,src);
      gsl_matrix_complex_set(dsigmadt,src,src,gsl_complex_sub(ztmp,dpdt));
      // d/dt p(dest) = 1/2*k*p(src)*p(src)
      dpdt=gsl_complex_mul_real(gsl_matrix_complex_get(tfs->sigma,src,src),0.5*k);
      ztmp=gsl_matrix_complex_get(dsigmadt,dest,dest);
      gsl_matrix_complex_set(dsigmadt,dest,dest,gsl_complex_add(ztmp,dpdt));
      for(a=0;a<ndim;a++) {
	if(a != src && a != dest) {
	  dpdt=gsl_complex_mul_real(gsl_matrix_complex_get(tfs->sigma,src,a),0.5*k);
	  ztmp=gsl_matrix_complex_get(dsigmadt,src,a);
	  gsl_matrix_complex_set(dsigmadt,src,a,gsl_complex_sub(ztmp,dpdt));
	  dpdt=gsl_complex_mul_real(gsl_matrix_complex_get(tfs->sigma,a,src),0.5*k);
	  ztmp=gsl_matrix_complex_get(dsigmadt,a,src);
	  gsl_matrix_complex_set(dsigmadt,a,src,gsl_complex_sub(ztmp,dpdt));
	}
      }
    }
  }

  // done!! clean up
  gsl_matrix_complex_free(mtmp);
  gsl_matrix_complex_free(sumOp);
}

/* initialize a bath aux function from mt99art parameters */
void bath_aux_funcs_init(bath_aux_funcs *f, const gsl_matrix_complex *sigma, 
			 gsl_matrix_complex *Sbar, const gsl_matrix_complex *H0, 
			 const gsl_vector *lambda,
			 mt99art_bath_func *bf)
{
  int i;
  size_t ndim;

  ndim=Sbar->size1;

  f->Sbar=Sbar;
  f->Nr=bf->Nr;
  f->Ni=bf->Ni;
  f->Vren=bf->Vren;

  /* fill in each bath funcs */
  for(i=0;i<f->Nr;i++) {
    f->Ar[i]=bf->Ar[i];
    f->Gr[i]=bf->Gr[i]; 
    f->dm_r[i]=gsl_matrix_complex_alloc(ndim,ndim);
    // initialize bath aux function using the equilibrium results
    //    compute_bath_aux_dm_r_eq(f->dm_r[i], sigma, Sbar, H0, f->Gr[i]);
    //    compute_bath_aux_dm_r_ss(f->dm_r[i], sigma, Sbar, f->Gr[i]);
    // initialize as zero matrix
    gsl_matrix_complex_set_zero(f->dm_r[i]);
    // get the initial derivative; initial field is zero
    f->ddm_rdt[i]=gsl_matrix_complex_alloc(ndim,ndim);
    compute_bath_aux_dm_r_dt(f->ddm_rdt[i], f->dm_r[i], sigma, Sbar, H0, f->Gr[i]);
  }
  for(i=0;i<f->Ni;i++) {
    f->Ai[i]=bf->Ai[i];
    f->Gi[i]=bf->Gi[i];
    f->dm_i[i]=gsl_matrix_complex_alloc(ndim,ndim);
    // initialize bath aux function using the equilibrium results
    //    compute_bath_aux_dm_i_eq(f->dm_i[i], sigma, Sbar, H0, f->Gi[i]);
    //    compute_bath_aux_dm_i_ss(f->dm_i[i], sigma, Sbar, f->Gi[i]);
    // initialize as zero matrix
    gsl_matrix_complex_set_zero(f->dm_i[i]);
    // get the initial derivative; initial field is zero
    f->ddm_idt[i]=gsl_matrix_complex_alloc(ndim,ndim);
    compute_bath_aux_dm_i_dt(f->ddm_idt[i], f->dm_i[i], sigma, f->Sbar, H0, f->Gi[i]);
  }

}

void bath_aux_funcs_free(bath_aux_funcs *f)
{
  int i;
 
 for(i=0;i<f->Nr;i++) {
   gsl_matrix_complex_free(f->dm_r[i]);
   gsl_matrix_complex_free(f->ddm_rdt[i]);
 }
 for(i=0;i<f->Ni;i++) {
   gsl_matrix_complex_free(f->dm_i[i]);
   gsl_matrix_complex_free(f->ddm_idt[i]);
 }
}

/* initialize a TNL density matrix and bath aux functions,
   use the bath information from the mt99art module. */
void tnl_aux_funcs_init(tnl_aux_funcs *f, size_t ndim,
			const gsl_matrix_complex *initdm, 
			const gsl_matrix_complex *Heff0, 
			const gsl_matrix_complex *H0, 
			const gsl_vector *lambda, const gsl_matrix *U,
			const qdas_keys *keys)
{
  size_t i;
  gsl_matrix_complex *mtmp=gsl_matrix_complex_alloc(ndim,ndim);
  gsl_matrix_complex *Sbar=gsl_matrix_complex_alloc(ndim,ndim);
  gsl_matrix_complex *UU=gsl_matrix_complex_alloc(ndim,ndim);

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

  gsl_matrix_complex_copy_real(UU,U);

  f->sigma=gsl_matrix_complex_alloc(ndim,ndim);
  f->dsigmadt=gsl_matrix_complex_alloc(ndim,ndim);

  gsl_matrix_complex_memcpy(f->sigma,initdm);
  
  /* handle dynamical bath parts */
  for(i=0;i<BATH_MT99ART_NOps;i++) {
    bath_aux_funcs_init(f->bfs+i,f->sigma,BATH_MT99ARTOpSbar[i],
			H0,lambda,BATH_MT99ARTBathFunc+i);
  }
  f->nops=BATH_MT99ART_NOps;

  // finally fill in the initial derivative too...
  compute_dsigma_dt(f->dsigmadt, f, Heff0,keys);

  gsl_matrix_complex_free(mtmp);
  gsl_matrix_complex_free(Sbar);
  gsl_matrix_complex_free(UU);
}

void tnl_aux_funcs_free(tnl_aux_funcs *f)
{
  int i;
  gsl_matrix_complex_free(f->sigma);
  gsl_matrix_complex_free(f->dsigmadt);
  for(i=0;i<f->nops;i++) {
    bath_aux_funcs_free(f->bfs+i);
  }
}

/* set f1=f2; note that only time dependent parts are copied */
void copy_tnl_aux_funcs(tnl_aux_funcs *f1, const tnl_aux_funcs *f2)
{
  int n,i;

  gsl_matrix_complex_memcpy(f1->sigma,f2->sigma);
  gsl_matrix_complex_memcpy(f1->dsigmadt,f2->dsigmadt);

  for(n=0;n<f1->nops;n++) {
    for(i=0;i<f1->bfs[n].Nr;i++) {
      gsl_matrix_complex_memcpy(f1->bfs[n].dm_r[i],f2->bfs[n].dm_r[i]);
      gsl_matrix_complex_memcpy(f1->bfs[n].ddm_rdt[i],f2->bfs[n].ddm_rdt[i]);
    }
    for(i=0;i<f1->bfs[n].Ni;i++) {
      gsl_matrix_complex_memcpy(f1->bfs[n].dm_i[i],f2->bfs[n].dm_i[i]);
      gsl_matrix_complex_memcpy(f1->bfs[n].ddm_idt[i],f2->bfs[n].ddm_idt[i]);
    }
  }

}

/* 
   setup the eigenbasis and dipole matrix in the eigen basis,
   this function should be used before anything else in this
   tnl-kernel module.

   return variables:
   U: coefficient matrix for eigenstates, U^-1*H*U diagonalized H.
   lambda: eigenvalues
   mX: the upper-half of the dipole operator; mX=\sum_{a}V_{ga}*|g><a|
   Note that mX is actually a vector; we return mX_x, mX_y, mX_z, mX_abs when
   the memory space is allocated (pointer not NULL)...

   This function also adjusts the indices in poptrans terms in keys->ptrans,
   which can make sure that level crossover in disordered system (MC or GH runs) 
   does not affect the assignment of poptrans terms.
 */
void setup_eigenbasis(gsl_matrix *U, gsl_vector *lambda, 
		      gsl_matrix *mX_x, gsl_matrix *mX_y, gsl_matrix *mX_z, gsl_matrix *mX_abs,
		      qdas_keys *keys, const gsl_matrix *H)
{
  gsl_eigen_symmv_workspace *w;
  gsl_matrix *mtmp=gsl_matrix_alloc(keys->nsize,keys->nsize);
  gsl_matrix *matX=gsl_matrix_alloc(keys->nsize,keys->nsize);
  gsl_matrix *matY=gsl_matrix_alloc(keys->nsize,keys->nsize);
  gsl_matrix *matZ=gsl_matrix_alloc(keys->nsize,keys->nsize);
  gsl_vector *vtmp=gsl_vector_alloc(keys->nsize);

  int i;
  size_t ndim;
  size_t a,b;
  size_t map[BUFFER_SIZE];

  printf("\n");
  printf("Setting up and transforming operators to the eigenbasis\n");

  ndim=H->size1;

  /* diagonalize H to obtain eigenstates and U */
  gsl_matrix_memcpy(mtmp,H);
  w=gsl_eigen_symmv_alloc(keys->nsize);
  gsl_eigen_symmv(mtmp, lambda, U, w);
  gsl_eigen_symmv_sort(lambda,U,GSL_EIGEN_SORT_VAL_ASC);
  gsl_eigen_symmv_free(w);

  /* adjust indices in eigenbasis so that the order is consistent with the 
     zero-th order eigenstates.
     basically we compare eigen states of two Hamiltonians and adjust the 
     indices for the current instance of H. The reason of doing this is that 
     the assignment of poptrans states in the keyword file and the initial 
     density density matrix in the input file assume the state
     labeling of the zeroth order (averaged) Hamiltonian given in the .key file.
     When static disorder is present, the state index can change. Therefore, we
     adjust the indices of poptrans terms to
     make sure that level crossover in disordered system (MC or GH runs) 
     does not affect the assignment of poptrans terms.
  */
  compute_eigen_map(map, keys->He, H);
  gsl_vector_memcpy(vtmp,lambda);
  gsl_matrix_memcpy(mtmp,U);
  for(a=0;a<ndim;a++) {
    //    printf("%d --> %d\n",map[a],a);
    if(map[a] != a) {
      printf("Eigenstate level crossover detected: reassign %lu-th to %lu-th.\n",map[a],a);
    }
    gsl_vector_set(lambda,a,gsl_vector_get(vtmp,map[a]));
    for(b=0;b<ndim;b++) {
      gsl_matrix_set(U,b,a,gsl_matrix_get(mtmp,b,map[a]));
    }
  }

  /* setup the upper-half of dipole operator in the site-localized basis... */
  gsl_matrix_set_zero(matX);
  gsl_matrix_set_zero(matY);
  gsl_matrix_set_zero(matZ);
  for(i=0;i<keys->ndipoles;i++) {
    gsl_matrix_set(matX,keys->mu[i].n,keys->mu[i].m,keys->mu[i].vec[0]); // note n<m
    gsl_matrix_set(matY,keys->mu[i].n,keys->mu[i].m,keys->mu[i].vec[1]);
    gsl_matrix_set(matZ,keys->mu[i].n,keys->mu[i].m,keys->mu[i].vec[2]);
  }
  /* U^\daggger*X*U to convert mat[XYZ] to exciton basis */
  gsl_blas_dgemm (CblasTrans, CblasNoTrans, 1.0, U, matX, 0.0, mtmp);
  gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, mtmp, U, 0.0, matX);
  gsl_blas_dgemm (CblasTrans, CblasNoTrans, 1.0, U, matY, 0.0, mtmp);
  gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, mtmp, U, 0.0, matY);
  gsl_blas_dgemm (CblasTrans, CblasNoTrans, 1.0, U, matZ, 0.0, mtmp);
  gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, mtmp, U, 0.0, matZ);

  if(mX_x) {
    gsl_matrix_memcpy(mX_x,matX);
  }
  if(mX_y) {
    gsl_matrix_memcpy(mX_y,matY);
  }
  if(mX_z) {
    gsl_matrix_memcpy(mX_z,matZ);
  }
  if(mX_abs) {
    gsl_matrix_set_zero(mX_abs);
    for(i=0;i<ndim;i++) {
      double mx,my,mz;
      int j;
      mx=0.0;my=0.0;mz=0.0;
      for(j=0;j<ndim;j++) {
	mx=gsl_matrix_get(matX,i,j);
	my=gsl_matrix_get(matY,i,j);
	mz=gsl_matrix_get(matZ,i,j);
	gsl_matrix_set(mX_abs,i,j,sqrt(mx*mx+my*my+mz*mz));
      }
    }
  }

  if(keys->printlv >= 3) {  
    printf("\n");
    printf("lambda = \n");
    gsl_vector_print(lambda);
    printf("\n");
    printf("U = \n");
    gsl_matrix_print(U);
    printf("\n");
    printf("mX_abs = \n");
    gsl_matrix_print(mX_abs);
    printf("\n");
  }

  // fill in the Lindblad operator in the exciton basis for
  // all site-representation population transfer terms
  for(i=0;i<keys->nsiteptrans;i++) {
    gsl_matrix_set_zero(matX);
    //    gsl_matrix_set(matX,keys->extra_siteptrans[i].src,keys->extra_siteptrans[i].dest,1.0);
    gsl_matrix_set(matX,keys->extra_siteptrans[i].dest,keys->extra_siteptrans[i].src,1.0); /* FIXME: which setup for lindblad op is right? kind of confused... */
    /* U^\daggger*A*U to convert A to exciton basis */
    gsl_blas_dgemm (CblasTrans, CblasNoTrans, 1.0, U, matX, 0.0, mtmp);
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, mtmp, U, 0.0, matX);
    /* save in a complex matrix */
    gsl_matrix_complex_copy_real(keys->extra_siteptrans[i].Op1,matX);
    /* also A^\dagger*A */
    gsl_matrix_memcpy(mtmp,matX);
    gsl_blas_dgemm (CblasTrans, CblasNoTrans, 1.0, mtmp, matX, 0.0, matY);
    gsl_matrix_complex_copy_real(keys->extra_siteptrans[i].Op2,matY);
  }

  printf("Eigenbasis set and now we are ready to propagate the dynamics.\n");
  printf("\n");

  gsl_matrix_free(mtmp);
  gsl_matrix_free(matX);
  gsl_matrix_free(matY);
  gsl_matrix_free(matZ);

}

/* update using Crank-Nicholson scheme:
   f(t+dt) = f(t) + dt/2 * (df(t)/dt + df(t+dt)/dt) */
void cn_update(gsl_matrix_complex *ftpdt, const gsl_matrix_complex *ft, 
	       const gsl_matrix_complex *dftdt, const gsl_matrix_complex *dftpdtdt,
	       double delta_t, gsl_matrix_complex *mtmp)
{
  gsl_matrix_complex_memcpy(ftpdt,ft);
  gsl_matrix_complex_memcpy(mtmp,dftdt);
  gsl_matrix_complex_add(mtmp,dftpdtdt);
  gsl_matrix_complex_scale(mtmp,gsl_complex_rect(delta_t/2.0,0.0));
  gsl_matrix_complex_add(ftpdt,mtmp);
}


/* update the tnl aux functions using the iterative 
Crank-Nicholson scheme, the previous tfs is given to 
save computation time */
void advance_tnl_aux_funcs(tnl_aux_funcs *tfs, const tnl_aux_funcs *tfs_old,
			   const gsl_matrix_complex *Weff_t, const gsl_matrix_complex *W_t, 
			   const gsl_vector *lambda, const qdas_keys *keys, double delta_t)
{
  size_t ndim,n,i,niter;
  gsl_matrix_complex *mtmp,*sigma_new;
  double conv;

  ndim=tfs->sigma->size1;

  mtmp=gsl_matrix_complex_alloc(ndim,ndim);
  sigma_new=gsl_matrix_complex_alloc(ndim,ndim);

  /* construct first guess from previous time step */
  copy_tnl_aux_funcs(tfs,tfs_old);
  cn_update(tfs->sigma,tfs_old->sigma,tfs->dsigmadt,tfs_old->dsigmadt,delta_t,mtmp);
  
  /* start iteration */
  niter=1;
  do {
    /* update aux functions and time derivatives */
    for(n=0;n<tfs->nops;n++) {
      for(i=0;i<tfs->bfs[n].Nr;i++) {
	cn_update(tfs->bfs[n].dm_r[i],tfs_old->bfs[n].dm_r[i],
		  tfs->bfs[n].ddm_rdt[i],tfs_old->bfs[n].ddm_rdt[i],delta_t,mtmp);
	compute_bath_aux_dm_r_dt(tfs->bfs[n].ddm_rdt[i],tfs->bfs[n].dm_r[i],
				 tfs->sigma, tfs->bfs[n].Sbar, W_t, tfs->bfs[n].Gr[i]);
      }
      for(i=0;i<tfs->bfs[n].Ni;i++) {
	cn_update(tfs->bfs[n].dm_i[i],tfs_old->bfs[n].dm_i[i],
		  tfs->bfs[n].ddm_idt[i],tfs_old->bfs[n].ddm_idt[i],delta_t,mtmp);
	compute_bath_aux_dm_i_dt(tfs->bfs[n].ddm_idt[i],tfs->bfs[n].dm_i[i],
				 tfs->sigma, tfs->bfs[n].Sbar, W_t, tfs->bfs[n].Gi[i]);
      }
    }
    compute_dsigma_dt(tfs->dsigmadt, tfs, Weff_t,keys);

    /* recalculate the density matrix and check convergence */
    cn_update(sigma_new,tfs_old->sigma,tfs->dsigmadt,tfs_old->dsigmadt,delta_t,mtmp);
    gsl_matrix_complex_memcpy(mtmp,sigma_new);
    gsl_matrix_complex_sub(mtmp,tfs->sigma);
    conv=gsl_matrix_complex_absmax(mtmp);
    if(conv < keys->cnconv) {
      //      printf("Convergence achieved after %ld iterations (maxabs=%20.16f)\n",niter,conv);
      gsl_matrix_complex_memcpy(tfs->sigma,sigma_new);
      break;
    }

    /* not converged, update data and do next run */
    gsl_matrix_complex_memcpy(tfs->sigma,sigma_new);
    //    printf("conv=%20.16f, Next run!\n",conv);
    niter++;
  } while(1);

  // DONE
  gsl_matrix_complex_free(sigma_new);
  gsl_matrix_complex_free(mtmp);
}


/*
 * $Log$
 * Revision 1.19  2007/08/09 23:50:46  platin
 *   - implemented site localized extra population dynamics in TNL.
 *   - implemented the DIAGONALDYNAMICS method for population-only
 *     incoherent dynamics.
 *
 * Revision 1.18  2007/06/23 00:54:49  platin
 *
 *   - bump to 0.8p3.
 *   - add support for index adjustment so that poptrans terms can be
 *     consistent when level energy crossover happens at disordered jobs.
 *
 * Revision 1.17  2007/06/01 18:00:02  platin
 *
 *   - bump to 0.8p1
 *   - improved the Gaussian-Hermite module; now users can select between
 *     1-20 point rules can be used.
 *   - added keyword to select between MC and GH samplings.
 *   - static disorder average loop in both tnl-dm3pes and tnl-dynamics
 *     is fixed.
 *
 * Revision 1.16  2007/03/26 19:15:25  platin
 *
 *   - bump to v0.8.
 *   - fix the pulse duration factor; the fwhm value is the field fwhm.
 *
 * Revision 1.15  2007/03/15 08:21:06  platin
 *
 *   - support MT99ART bath.
 *
 * Revision 1.14  2007/03/15 02:37:05  platin
 *
 *   - add support for the multimode browanian oscillator bath.
 *
 * Revision 1.13  2007/03/06 23:51:49  platin
 *
 *   - clean up BS disorder_aux_xxx code.
 *
 * Revision 1.12  2007/03/06 23:32:01  platin
 *
 *   - remove useless -markov codes.
 *   - use SS initial condition for all bath aux functions.
 *
 * Revision 1.11  2007/02/21 06:14:59  platin
 *
 *   - add comment on the markov type,, it is not working.
 *
 * Revision 1.10  2007/02/20 22:02:51  platin
 *
 *   - implemented Markov approximation in tnl-kern.c
 *   - tnl-dm3pes-2d now used the Markov model.
 *
 * Revision 1.9  2006/12/22 01:50:02  platin
 *
 *   - fixed stupid errors.
 *
 * Revision 1.8  2006/12/22 00:30:14  platin
 *
 *   - support nonlinear pairwise rxn via Linblad formalism.
 *
 * Revision 1.7  2006/12/20 07:24:49  platin
 *
 *   - bug fix.
 *
 * Revision 1.6  2006/12/20 07:21:23  platin
 *
 *   - nonlinear poptrans stuff.
 *
 * Revision 1.5  2006/12/20 07:12:38  platin
 *
 *   - add a indicative mesage for poptrans.
 *
 * Revision 1.4  2006/12/20 07:07:22  platin
 *
 *   - remove redudant code.
 *   - add options of using nonlinear poptrans.
 *
 * Revision 1.3  2006/11/29 23:17:15  platin
 *
 *   - add a decay term in the pseudo- static disorder bath functions,
 *     this might serve as an approximation to include static disorder...
 *
 * Revision 1.2  2006/11/02 19:29:06  platin
 *
 *   - change the propagation kernel to use the new bathtype interface.
 *     the programs can handle different MT99ART type baths now.
 *
 * Revision 1.1  2006/10/27 22:24:54  platin
 *
 *   - commit re-organized tnl-dm3pes and tnl-dynamics code.
 *
 *
 */
