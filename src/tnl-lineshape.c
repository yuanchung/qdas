
/***************************************************
 * tnl-lineshape.c
 * 
 * By Yuan-Chung Cheng <yuanchung@ntu.edu.tw>
 *
 * Simulating linear absorption spectrum of a molecular
 * system using density-matrix based method.
 * See Notebok 2, p. 49.
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

#include "params.h"
#include "aux.h"
#include "dmop.h"
#include "bath.h"
#include "bathtype_mt99art.h"
#include "tnl-kern.h"
#include "gauss-hermite.h"
#include "samplelist.h"

/* Global constants */

/* global variables */
const char *program_name;

/* Macros */

/* Display usage information and exit.  */
static void usage ()
{
  printf ("\
Usage: %s [OPTIONS] KEYFILE \n\
Use the parameters described in the KEYFILE to\n\
simulate absorption spectrum of the system.\n\
\n\
OPTIONS:\n\
\n\
    -nnn   No options yet.\n",
	  program_name);
}


/* Use density-matrix based method to propagate the
   polarization in the eigenbasis; the dissipation part is 
   implemented using the QME method in Notebook 2, page 49-50 */
void propagate(const double tstop, const double tstep, 
	       gsl_vector_complex *Pt,
	       qdas_keys *keys, const gsl_matrix *H)
{
  gsl_matrix *U=gsl_matrix_alloc(keys->nsize,keys->nsize);
  gsl_matrix *mX=gsl_matrix_alloc(keys->nsize,keys->nsize);
  gsl_vector *lambda=gsl_vector_alloc(keys->nsize);

  tnl_aux_funcs *sigma_t,*sigma_tpdt; /* sigma(t) and sigma(t+dt) */

  gsl_matrix_complex *mtmp=gsl_matrix_complex_alloc(keys->nsize,keys->nsize);
  gsl_matrix_complex *initdm=gsl_matrix_complex_alloc(keys->nsize,keys->nsize);
  gsl_matrix_complex *zU=gsl_matrix_complex_alloc(keys->nsize,keys->nsize);
  gsl_matrix_complex *Weff_t=gsl_matrix_complex_alloc(keys->nsize,keys->nsize);
  gsl_matrix_complex *Hren=gsl_matrix_complex_alloc(keys->nsize,keys->nsize);

  gsl_matrix_complex *Heigen=gsl_matrix_complex_alloc(keys->nsize,keys->nsize);

  size_t i,j,ndim;

  double t;
  double dtmp;
  size_t a,t_idx;

  size_t pcounter,poutput;

  // fill global vars
  ndim=keys->nsize;

  /* obtain U/lambda, and mX_abs, mX will be used to 
     construct the initial sigma too */
  setup_eigenbasis(U,lambda,NULL,NULL,NULL,mX,keys,H);
  /* setup_eigenbasis() returns the upper-diagonal part of mX,
     we want full one here. */
  for(i=0;i<ndim;i++) {
    for(j=0;j<i;j++) {
      gsl_matrix_set(mX,i,j,gsl_matrix_get(mX,j,i));
    }
  }
  printf("mX = \n");
  gsl_matrix_print(mX);
  printf("\n");

  /* initial bath in the exciton basis */
  bath_mt99art_init_OpSbar(keys,U);

  /* Hamiltonian in the eigen basis */
  gsl_matrix_complex_set_zero(Heigen);
  for(a=0;a<ndim;a++) {
    dtmp=gsl_vector_get(lambda,a);
    gsl_matrix_complex_set(Heigen,a,a,gsl_complex_rect(dtmp,0.0));
  }

  /* initial density matrix is |g><g|*mX */
  gsl_matrix_complex_set_zero(mtmp);
  gsl_matrix_complex_set(mtmp,0,0,gsl_complex_rect(1.0,0.0)); // |g><g| in site basis
  // U'* |g><g| *U to exciton basis, use Weff_t as tmp vars
  gsl_matrix_complex_copy_real(zU,U);
  gsl_blas_zgemm(CblasConjTrans,CblasNoTrans,gsl_complex_rect(1.0,0.0), zU, mtmp,
		 gsl_complex_rect(0.0,0.0), Weff_t); // Weff_t = U'*sigma
  gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,gsl_complex_rect(1.0,0.0), Weff_t, zU,
		 gsl_complex_rect(0.0,0.0), mtmp); // mtmp = Weff_t*U
  // use Weff_t as tmp vars
  gsl_matrix_complex_copy_real(Weff_t,mX);
  gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,gsl_complex_rect(1.0,0.0), mtmp, Weff_t,
		 gsl_complex_rect(0.0,0.0), initdm); // initdm = mtmp*mX
  

  /* the renormalization part,
     Hren = \sum_n mu_n/2 * Sn*Sn; 
     note we initialize sigma(t) in order to obtain bath 
     functions, but then free and re-initialize sigma(t) 
     to include the effect of Hren */
  sigma_t=(tnl_aux_funcs *)malloc(sizeof(tnl_aux_funcs));
  tnl_aux_funcs_init(sigma_t,ndim,initdm,Heigen,Heigen,lambda,U,keys);
  gsl_matrix_complex_set_zero(Hren);
  for(a=0;a<sigma_t->nops;a++) {
    gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,gsl_complex_rect(1.0,0.0), 
  		   sigma_t->bfs[a].Sbar, sigma_t->bfs[a].Sbar, 
  		   gsl_complex_rect(0.0,0.0), mtmp); 
    gsl_matrix_complex_scale(mtmp,gsl_complex_rect(sigma_t->bfs[a].Vren/2.0,0.0));
    gsl_matrix_complex_add(Hren,mtmp);
  }
  tnl_aux_funcs_free(sigma_t);
  free(sigma_t);
  printf("Hren = \n");
  gsl_matrix_complex_print(Hren);

  /* aux density matrices */
  gsl_matrix_complex_memcpy(mtmp,Hren); 
  gsl_matrix_complex_add(mtmp,Heigen);
  sigma_t=(tnl_aux_funcs *)malloc(sizeof(tnl_aux_funcs));
  tnl_aux_funcs_init(sigma_t,ndim,initdm,mtmp,Heigen,lambda,U,keys);
  sigma_tpdt=(tnl_aux_funcs *)malloc(sizeof(tnl_aux_funcs));
  tnl_aux_funcs_init(sigma_tpdt,ndim,initdm,mtmp,Heigen,lambda,U,keys);

  // start prograte in time
  pcounter=0; 
  poutput=(size_t)rint(keys->tprint/tstep);

  t_idx=0;
  for(t=0.0;t<=tstop;t+=tstep) {

    tnl_aux_funcs *taf_ptr;
    
    if((pcounter % poutput) == 0) {
      gsl_complex Pval;
      /* tprint check points, save the P(t);
	 P(t) = Tr{ mX*sigma(t)} */
      // use Weff_t as tmp vars
      gsl_matrix_complex_copy_real(Weff_t,mX);
      gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,gsl_complex_rect(1.0,0.0), Weff_t, sigma_t->sigma,
		     gsl_complex_rect(0.0,0.0), mtmp);
      Pval=gsl_matrix_complex_trace(mtmp);
      gsl_vector_complex_set(Pt,t_idx,Pval);
      t_idx++;
      /*
      printf("t=%12.4f, Pt=%18.12f %18.12f\n",t*TIME_CM2FS,GSL_REAL(Pval),GSL_IMAG(Pval));      
      printf("t=%12.4f, sigma=",t*TIME_CM2FS);      
      dmop_lvprint(sigma_t->sigma);
      printf("\n");
      */
    } //if((pcounter % poutput) == 0)
    pcounter++;

    /* prepare for the t+dt step */

    /* Weff_t is the normalized effective potential,
       Weff(t) = H0 + \sum mu/2*Sn*Sn */
    gsl_matrix_complex_memcpy(Weff_t,Heigen);
    gsl_matrix_complex_add(Weff_t,Hren);
    advance_tnl_aux_funcs(sigma_tpdt,sigma_t,Weff_t, Heigen,lambda,keys,tstep);

    /* done, ready for next run (swap pointers) */
    taf_ptr=sigma_t;
    sigma_t=sigma_tpdt;
    sigma_tpdt=taf_ptr;
  } // for(t...

  /* clean up */
  bath_mt99art_free_OpSbar();

  gsl_matrix_free(U);
  gsl_matrix_free(mX);
  gsl_vector_free(lambda);
  gsl_matrix_complex_free(mtmp);
  gsl_matrix_complex_free(initdm);
  gsl_matrix_complex_free(zU);
  gsl_matrix_complex_free(Weff_t);
  gsl_matrix_complex_free(Heigen);
  gsl_matrix_complex_free(Hren);

  tnl_aux_funcs_free(sigma_t);
  tnl_aux_funcs_free(sigma_tpdt);
  free(sigma_t);
  free(sigma_tpdt);
}

/* Main program */
int main(int argc, char *argv[])
{

  /* variables for the options */
  char *key_file_name=NULL;
  
  int i;
  
  qdas_keys keys;
  
  /* Set program name for messages.  */
  program_name = argv[0]; 
  
  /* HEADER MESSAGE */
  printf("\n");
  printf("        TNL-LINESHAPE\n");
//  printf("        SVN revision: %s.\n","$Revision: 259 $");
  printf("        Linear absorption spectrum using\n");
  printf("        TNL Non-markovian quantum master equation method.\n");
  printf("\n");
  printf("        Part of the QDAS package, version %s.\n",QDAS_VERSION);
  printf("\n");
  printf("        Copyright(C) 2007.\n");
  printf("        Yuan-Chung Cheng <yuanchung@ntu.edu.tw>.\n");
  printf("\n");
  printf("\n");
  printf("Job started: ");
  print_time_stamp();
  printf("\n");
  
  /* parse the arguments... */
  for(i=1;i<argc;i++){
    if(argv[i][0] != '-') {
      if(key_file_name == NULL) {
	key_file_name=(char*)strdup(argv[i]);
      } else {
	usage();
	exit(EXIT_FAILURE);
      }
    } else {
      /* Options go here */
      usage();
      exit(EXIT_FAILURE);
    }
  }

  if(key_file_name == NULL) {
    usage();
    exit(EXIT_FAILURE);
  }

  /* initialize bath modules before we read the keyword file; 
     the bath parameters will be read and prepared in params.c */
  bath_mod_init();

  /* use the key file to initialize all parameters */
  params_init(key_file_name,&keys);

  if(keys.npulses > 0) {
    printf("Pulse assignment(s) not allowed in abs. spectrum simulation.\n");
    exit(EXIT_FAILURE);
  }
  if(!keys.ndipoles) {
    printf("Required keyword \"DIPOLE\" not found!\n");
    exit(EXIT_FAILURE);
  }
  printf("\n");
  printf("\n");
  printf("------------------------\n");
  printf("| Dynamical Simulation |\n");
  printf("------------------------\n");
  printf("\n");
  printf("Propagate time-domain linear response using Crank-Nicholson scheme.\n");
  printf("Iteration convegence threshold is %g.\n",keys.cnconv);
  printf("\n");
  printf("Simulate trajectory using time increment %f fs.\n",keys.tstep*TIME_CM2FS);
  printf("\n");
  printf("Will sample time domain polarization from t=0.0 to t=%f fs,\n",
	 (keys.tstop)*TIME_CM2FS);
  printf("with time step size delta_t=%f fs. This setup will be\n",keys.tprint*TIME_CM2FS);
  printf("used in FFT to obtain the spectrum\n");
  printf("\n");

  /* Parameters read and initialized, ready to do real work */
  {
    gsl_vector_complex *Pt=NULL;
    gsl_vector_complex *Pt_sum=NULL;
    size_t nn,narray;
    qdas_sample_item *joblist=NULL;
    size_t iter;

    /* use multiple trajectories to account for static disorder is required,
       so we need to store the sum of P(t) trajectories. */
    /* we use FFT so the size of P(t) array should be 2^n */
    
    nn=(size_t)ceil((fabs(keys.tstop))/keys.tprint)+1;
    narray=1024; // start with 1024 poins, and intensionally make room for padding...
    while(narray<=(nn+1)) 
      narray=narray*2;
    Pt=gsl_vector_complex_alloc(narray);
    Pt_sum=gsl_vector_complex_alloc(narray);

    gsl_vector_complex_set_zero(Pt);
    gsl_vector_complex_set_zero(Pt_sum);

    joblist=prepare_job_list(&keys);

    iter=0;
    while(joblist) {

      printf("Sample #%lu (weight=%20.16f): ",iter+1,joblist->weight);
      print_time_stamp();
      printf("\n");
#ifdef DEBUG
      printf("\n");
      printf("Hamiltonian = \n");
      gsl_matrix_print(joblist->H);
      printf("\n");
#endif
      /* now indeed propagate the dynamics and obtain the polarization; note that we ignore joblist->P */
      propagate(keys.tstop,keys.tstep,Pt,&keys,joblist->H);
      gsl_vector_complex_add(Pt_sum,Pt);
      /* fflush stdout to make the message more prompt */
      fflush(stdout);
      iter++;
      joblist=joblist->next;
    }

    /* print out the averaged results */
    {
      double t;
      int t_idx;
      gsl_complex value;

      // scale the sum to get the average
      gsl_vector_complex_scale(Pt_sum,gsl_complex_rect(1.0/((double)iter),0.0));
      /* show averaged P(t) */
      // This corresponds to the layout in propagate()
      printf("\n");
      printf("\n");
      printf("-------------------------\n");
      printf("| Linear Polarizability |\n");
      printf("-------------------------\n");
      printf("\n");
      printf("Time evolution of the averaged linear polarizability:\n");
      printf("\n");
      for(t=0.0,t_idx=0;t<=(keys.tstop+keys.tstep/10.0);t+=keys.tprint,t_idx++) { 
	value=gsl_vector_complex_get(Pt_sum,t_idx);
	printf("t=%12.4f, Pavg=%20.16f  %20.16f, |Pavg|=%20.16f\n",
	       t*TIME_CM2FS,GSL_REAL(value),GSL_IMAG(value),gsl_complex_abs(value));
      }
      printf("\n");
    }

    /* FFT to obtain the spectrum */
    {
      gsl_vector *omega=gsl_vector_alloc(Pt_sum->size);
      gsl_vector_complex *Fomega=gsl_vector_complex_alloc(Pt_sum->size);
      printf("\n");
      printf("\n");
      printf("-----------------\n");
      printf("| Abs. Spectrum |\n");
      printf("-----------------\n");
      printf("\n");
      printf("Construct spectrum using FFT with %lu points.\n",omega->size);
      hFourier_forward(Pt_sum , keys.tprint, omega,Fomega);
      printf("\n");
      printf("Linear absorption spectrum:\n");
      // print backwards because we actually need backward FFT
      for(i=omega->size;i>0;i--) {
	if(-1.0*gsl_vector_get(omega,i-1) >= keys.spec_start &&
	   -1.0*gsl_vector_get(omega,i-1) <= keys.spec_end) {
	  printf("w=%12.4f, I=%20.16f\n",
		 -1.0*gsl_vector_get(omega,i-1),
		 GSL_REAL(gsl_vector_complex_get(Fomega,i-1)));
	}
      }
      printf("\n");
      gsl_vector_free(omega);
      gsl_vector_complex_free(Fomega);
    }

    gsl_vector_complex_free(Pt);
    gsl_vector_complex_free(Pt_sum);
  }

  
  /* every thing done; clean up */
  params_close();

  printf("\n");
  printf("Job finished: ");
  print_time_stamp();
  printf("\n");

  return 0;
}

/*
 * $Log$
 *
 */
