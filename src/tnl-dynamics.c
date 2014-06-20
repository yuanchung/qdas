
/***************************************************
 * tnl-dynamics.c
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

#include "tnl-dynamics.h"
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
Usage: %s [OPTIONS] KEYFILE [INPUTFILE] \n\
Use the parameters described in the KEYFILE to\n\
propagate the initial density matrix defined in \n\
the INPUTFILE. If no INPUTFILE is given, standard \n\
input is read. \n\
\n\
OPTIONS:\n\
\n\
    -site  The input density matrix is in the site basis.\n",
	  program_name);
}


/* build V(t) = \sum l_i * Vi_bar(t) + c.c.
   where Vi_bar(t) = E_i(t-tau_i) * \sum_a exp{I*w_i*t} * V_{ga} * |g><a|.
   See p. 73 of Notebook 1.
   note that mX_abs is the upper-half only, and the output V is the full matrix;
   also note that field polarizations are ignored here... */
void build_V_matrix(gsl_matrix_complex *V, 
		    const double t, const gsl_matrix *mX_abs,
		    const qdas_keys *keys)
{
  int i,a,b;
  gsl_complex ztmp1,ztmp2;
  double E0,tau0,tau_fwhm,w0;
  double Vab,field;
  size_t ndim;

  ndim=V->size1;

  gsl_matrix_complex_set_zero(V);

  for(i=0;i<keys->npulses;i++) {
    E0=keys->pulse_seq[i].E0;
    tau0=keys->pulse_seq[i].tau0;
    tau_fwhm=keys->pulse_seq[i].fwhm_tau;
    w0=keys->pulse_seq[i].w0;
    if(fabs(t-tau0) > 2.0*tau_fwhm) continue; // only use (+-) 2*FWHM pulse
    // FWHM of the field
    field=E0*exp(-2.0*1.385*(t-tau0)*(t-tau0)/tau_fwhm/tau_fwhm);
    /* consider full transitions, not only |g>->|e> transitions */
    for(a=0;a<ndim;a++) {
      for(b=a+1;b<ndim;b++) { // Vag=0 when a=b, we assume zero perment dipole.
	Vab=gsl_matrix_get(mX_abs,a,b);
	ztmp1=gsl_complex_polar(1.0,w0*(t-tau0));
	ztmp2=gsl_complex_mul_real(ztmp1,field*Vab);
	/* Vn_bar term */
	ztmp1=gsl_matrix_complex_get(V,a,b);
	gsl_matrix_complex_set(V,a,b,gsl_complex_add(ztmp1,ztmp2));
      } // b
    } // a 
  }
  /* now set the complex conjugate lower-half terms */
  for(a=0;a<ndim;a++) {
    for(b=0;b<a;b++) { // Vag=0 when a=b, we assume zero perment dipole.
      gsl_matrix_complex_set(V,a,b,gsl_complex_conjugate(gsl_matrix_complex_get(V,b,a)));
    }
  }

  // done
}

/* Use TNL-QME to propagate the density-matrix 
   in the eigenbasis; the dissipation part is 
   implemented using the QME method in Notebook 1, page 65-69 */
void propagate(const double tstop, const double tstep, 
	       const gsl_matrix_complex *dm0,
	       gsl_matrix_complex **dm_traj,
	       qdas_keys *keys, const gsl_matrix *H, int jcount, int siteinitdm)
{
  gsl_matrix *U=gsl_matrix_alloc(keys->nsize,keys->nsize);
  gsl_matrix *mX_abs=gsl_matrix_alloc(keys->nsize,keys->nsize);
  gsl_vector *lambda=gsl_vector_alloc(keys->nsize);

  tnl_aux_funcs *sigma_t,*sigma_tpdt; /* sigma(t) and sigma(t+dt) */

  gsl_matrix_complex *initdm=gsl_matrix_complex_alloc(keys->nsize,keys->nsize);
  gsl_matrix_complex *V=gsl_matrix_complex_alloc(keys->nsize,keys->nsize);
  gsl_matrix_complex *W_t=gsl_matrix_complex_alloc(keys->nsize,keys->nsize);
  gsl_matrix_complex *Weff_t=gsl_matrix_complex_alloc(keys->nsize,keys->nsize);
  gsl_matrix_complex *Hren=gsl_matrix_complex_alloc(keys->nsize,keys->nsize);

  gsl_matrix_complex *Heigen=gsl_matrix_complex_alloc(keys->nsize,keys->nsize);

  size_t ndim;

  double t;
  double dtmp;
  double tf;
  size_t a,t_idx;

  size_t pcounter,poutput;

  // fill global vars
  ndim=keys->nsize;

  // propagate tstop fs longer than the center of the third pulse
  tf=tstop;

  /* obtain U/lambda, and mX_abs */
  setup_eigenbasis(U,lambda,NULL,NULL,NULL,mX_abs,keys,H);
 
  /* transform initial density matrix if needed */
  if(siteinitdm == 1) {
    // U' * dm0 * U to eigen basis,  use V,W_t,Weff_t as tmp vars
    gsl_matrix_complex_copy_real(V,U);
    gsl_blas_zgemm(CblasConjTrans,CblasNoTrans,gsl_complex_rect(1.0,0.0), V, dm0,
		   gsl_complex_rect(0.0,0.0), W_t); // W_t = U*sigma(t)
    gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,gsl_complex_rect(1.0,0.0), W_t, V,
		   gsl_complex_rect(0.0,0.0), Weff_t); // Weff_t = W_t*U'
    gsl_matrix_complex_memcpy(initdm,Weff_t);
  } else {
    gsl_matrix_complex_memcpy(initdm,dm0);
  }

  /* initial bath in the exciton basis */
  bath_mt99art_init_OpSbar(keys,U);

  /* Hamiltonian in the eigen basis */
  gsl_matrix_complex_set_zero(Heigen);
  for(a=0;a<ndim;a++) {
    dtmp=gsl_vector_get(lambda,a);
    gsl_matrix_complex_set(Heigen,a,a,gsl_complex_rect(dtmp,0.0));
  }  

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
  		   gsl_complex_rect(0.0,0.0), V); // use V as tmp storage
    gsl_matrix_complex_scale(V,gsl_complex_rect(sigma_t->bfs[a].Vren/2.0,0.0));
    gsl_matrix_complex_add(Hren,V);
  }
  tnl_aux_funcs_free(sigma_t);
  free(sigma_t);
  printf("Hren = \n");
  gsl_matrix_complex_print(Hren);

  /* aux density matrices */
  gsl_matrix_complex_memcpy(V,Hren); // again use V as tmp variable to store Weff(0)
  gsl_matrix_complex_add(V,Heigen);
  sigma_t=(tnl_aux_funcs *)malloc(sizeof(tnl_aux_funcs));
  tnl_aux_funcs_init(sigma_t,ndim,initdm,V,Heigen,lambda,U,keys);
  sigma_tpdt=(tnl_aux_funcs *)malloc(sizeof(tnl_aux_funcs));
  tnl_aux_funcs_init(sigma_tpdt,ndim,initdm,V,Heigen,lambda,U,keys);

  // start prograte in time
  printf("\n");
  printf("Propagate the dynamics from t=0.00 to t=%.2f.\n",tf*TIME_CM2FS);
  printf("\n");

  /* propagate from time zero, note that the 
     first pulse should be in the positive time.
     FIXME: remove these restrictions!! */
  pcounter=0; 
  poutput=(size_t)rint(keys->tprint/tstep);

  t_idx=0;
  for(t=0.0;t<=tf;t+=tstep) {

    tnl_aux_funcs *taf_ptr;
    
    if(jcount>1 && (pcounter % poutput) == 0) {
      /* if this is a static-disorder multiple-iteration job, 
	 we store the density matrix in the dm_traj array 
	 and continue without printing anything unless in debugdebug mode */
      gsl_matrix_complex_add(dm_traj[t_idx],sigma_t->sigma);
      t_idx++;
    }

    if(
#ifndef DEBUGDEBUG
       jcount<=1 &&  // don't output detail trajectory info when static average is on, unless in DEBUG mode
#endif
       (pcounter % poutput) == 0) {
      /* print density matrix rho(t) */
      printf("t=%12.4f, rho=",t*TIME_CM2FS);      
      dmop_lvprint(sigma_t->sigma);
      printf("\n");
      // print Psite,, use V,W_t,Weff_t as tmp vars
      printf("t=%12.4f, rho_site=",t*TIME_CM2FS);
      gsl_matrix_complex_copy_real(V,U);
      // U*sigma*U' to site basis
      gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,gsl_complex_rect(1.0,0.0), V, sigma_t->sigma,
		     gsl_complex_rect(0.0,0.0), W_t); // W_t = U*sigma(t)
      gsl_blas_zgemm(CblasNoTrans,CblasConjTrans,gsl_complex_rect(1.0,0.0), W_t, V,
		     gsl_complex_rect(0.0,0.0), Weff_t); // Weff_t = W_t*U'
      dmop_lvprint(Weff_t);
      printf("\n");

#ifdef DEBUGDEBUGDEBUG
      printf("t=%12.4f, Psite=",t*TIME_CM2FS);
      for(a=0;a<ndim;a++) {
	gsl_complex ztmp;
	ztmp=gsl_matrix_complex_get(Weff_t,a,a);
	printf(" %14.10f",GSL_REAL(ztmp));
      }
      printf("\n");
      
      printf("Fsigma = \n");gsl_matrix_complex_print(sigma_t->sigma);
      printf("Fdsigmadt = \n");gsl_matrix_complex_print(sigma_t->dsigmadt);
      int n,i;
      for(n=0;n<sigma_t->nops;n++) {
	for(i=0;i<sigma_t->bfs[n].Nr;i++) {
	  printf("dm_r[%d] = \n",i);gsl_matrix_complex_print(sigma_t->bfs[n].dm_r[i]);
	}
	for(i=0;i<sigma_t->bfs[n].Ni;i++) {
	  printf("dm_i[%d] = \n",i);gsl_matrix_complex_print(sigma_t->bfs[n].dm_i[i]);
	}
      }
#endif
    } //if((pcounter % poutput) == 0)
    pcounter++;

    /* prepare for the t+dt step */
    
    /* update aux density matrices */

    /* first construct modulation parts Hm and save in W_t */
    if(keys->nsinmoduls>0) {
      gsl_complex ztmp;
      // use V, Weff_t as tmp vars, and save result in W_t
      gsl_matrix_complex_set_zero(W_t);
      // construct \\sum |n><m|*Jm*sin(2*pi*(t+dt)/Tm+Phi) in site basis
      for(a=0;a<keys->nsinmoduls;a++) {
	ztmp=gsl_matrix_complex_get(W_t,keys->smodul[a].n,keys->smodul[a].m);
	ztmp=gsl_complex_add(ztmp,
			     gsl_complex_rect(keys->smodul[a].Jm*sin(2.0*M_PI*(t+tstep)/keys->smodul[a].Tm + keys->smodul[a].Phi),0.0));
	gsl_matrix_complex_set(W_t,keys->smodul[a].n,keys->smodul[a].m,ztmp);
	if(keys->smodul[a].n != keys->smodul[a].m) {
	  gsl_matrix_complex_set(W_t,keys->smodul[a].m,keys->smodul[a].n,gsl_complex_conjugate(ztmp));
	}
      }
      // now transform into eigenbasis using U'*Hm*U
      gsl_matrix_complex_copy_real(V,U);
      gsl_blas_zgemm(CblasConjTrans,CblasNoTrans,gsl_complex_rect(1.0,0.0), V, W_t,
		     gsl_complex_rect(0.0,0.0), Weff_t); // Weff_t = U'*Hm
      gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,gsl_complex_rect(1.0,0.0), Weff_t, V,
		     gsl_complex_rect(0.0,0.0), W_t); // W_t = Weff_t*U
      //      printf("W(t) = \n");gsl_matrix_complex_print(W_t);
    } else {
      gsl_matrix_complex_set_zero(W_t);
    }

    /* here V=V(t+dt) is full media-field op */
    build_V_matrix(V,t+tstep,mX_abs,keys);
    //    printf("V = \n");gsl_matrix_complex_print(V);

    /* W_t is the time dependent Hamiltonian, W_t =  \sum |><|*Jm*sin(t/Tm) + H - V */
    gsl_matrix_complex_add(W_t,Heigen);
    gsl_matrix_complex_sub(W_t,V);
    /* Weff_t is the normalized effective potential,
       Weff(t) = W(t) + \sum mu/2*Sn*Sn */
    gsl_matrix_complex_memcpy(Weff_t,W_t);
    gsl_matrix_complex_add(Weff_t,Hren);

    advance_tnl_aux_funcs(sigma_tpdt,sigma_t,Weff_t, W_t,lambda,keys,tstep);

    /* done, ready for next run (swap pointers) */
    taf_ptr=sigma_t;
    sigma_t=sigma_tpdt;
    sigma_tpdt=taf_ptr;
  } // for(t...

  /* clean up */
  bath_mt99art_free_OpSbar();

  gsl_matrix_free(U);
  gsl_matrix_free(mX_abs);
  gsl_vector_free(lambda);
  gsl_matrix_complex_free(initdm);
  gsl_matrix_complex_free(V);
  gsl_matrix_complex_free(W_t);
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
  char *input_file_name=NULL;
  
  int siteinitdm = 0;

  size_t dmdim;
  gsl_matrix_complex *initdm;

  int i;
  
  qdas_keys keys;
  
  /* Set program name for messages.  */
  program_name = argv[0]; 
  
  /* HEADER MESSAGE */
  printf("\n");
  printf("        TNL-DYNAMICS\n");
//  printf("        SVN revision: %s.\n","$Revision: 269 $");
  printf("        Time-Nonlocal non-Markovia quantum dynamics\n");
  printf("        using quantum master equation method.\n");
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
      } else if(input_file_name == NULL) {
	input_file_name=(char*)strdup(argv[i]);
      } else {
	usage();
	exit(EXIT_FAILURE);
      }
    } else {
      if(!strncmp(argv[i],"-site",5)) {
	/* initial density matrix is in site basis */
	siteinitdm = 1;
      } else {
	/* new options go here */
	usage();
	exit(EXIT_FAILURE);
      }
    }
  }

  if(key_file_name == NULL || input_file_name == NULL) {
    usage();
    exit(EXIT_FAILURE);
  }

  /* read the input density matrix */
  dmdim=dmop_read(input_file_name,&initdm);

  printf("Input file: %s\n",input_file_name);
  printf("\n");
  printf("Dimension of the system: %ld\n",dmdim);
  printf("\n");
  printf("Initial density matrix at t=0:\n");
  if(siteinitdm == 1)
    printf("  *** Initial density matrix in the site basis ***\n");
  else
    printf("\n");
  dmop_dmprint(initdm);
  printf("\n");

  /* use the key file to initialize all parameters */
  params_init(key_file_name,&keys);

  /* Check the dynamics module and 
     print out parameters specifically for this module */
  if(keys.nsize != dmdim) {
    printf("NSIZE in %s is inconsistent with the dimension of density matrix!\n",key_file_name);
    exit(EXIT_FAILURE);
  }

  printf("\n");
  printf("\n");
  printf("------------------------\n");
  printf("| Dynamical Simulation |\n");
  printf("------------------------\n");
  printf("\n");
  printf("Propagate density matrix using Crank-Nicholson iterative scheme.\n");
  printf("Iteration convegence threshold is %g.\n",keys.cnconv);
  printf("\n");
  printf("Simulate trajectory using time increment %f fs.\n",keys.tstep*TIME_CM2FS);
  printf("\n");
  printf("Will output simulated trajectory from t=0.0 to t=%f fs,\n",
	 (keys.tstop)*TIME_CM2FS);
  printf("with time step size delta_t=%f fs.\n",keys.tprint*TIME_CM2FS);
  printf("\n");

  /* Parameters read and initialized, ready to do real work */
  {
    gsl_matrix_complex **dm_sum=NULL;
    gsl_matrix_complex **dm_traj=NULL;
    size_t narray;
    qdas_sample_item *joblist=NULL;
    size_t iter;
    size_t JJCOUNT;

    /* use multiple trajectories to account for static disorder is required,
       so we need to store the sum of dm trajectories. */
    narray=(size_t)ceil((fabs(keys.tstop))/keys.tprint)+1;
    dm_sum=(gsl_matrix_complex **)malloc(sizeof(gsl_matrix_complex *)*narray);
    dm_traj=(gsl_matrix_complex **)malloc(sizeof(gsl_matrix_complex *)*narray);
    for(i=0;i<narray;i++) {
      dm_sum[i]=gsl_matrix_complex_alloc(dmdim,dmdim);
      gsl_matrix_complex_set_zero(dm_sum[i]);
      dm_traj[i]=gsl_matrix_complex_alloc(dmdim,dmdim);
      gsl_matrix_complex_set_zero(dm_traj[i]);
    }

    joblist=prepare_job_list(&keys);
    JJCOUNT=job_list_count(joblist);

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
      /* now indeed propagate the dynamics; note that we ignore joblist->P */
      propagate(keys.tstop,keys.tstep,initdm,dm_traj,&keys, joblist->H,JJCOUNT,siteinitdm);
#ifdef DEBUGDEBUG
      printf("t=, rho=\n"); // this is for convenient gnuplot plot
      printf("t=, Psite=\n"); // this is for convenient gnuplot plot
#endif
      for(i=0;i<narray;i++) {
	gsl_matrix_complex_scale(dm_traj[i],gsl_complex_rect(joblist->weight,0.0));
	gsl_matrix_complex_add(dm_sum[i],dm_traj[i]);
      }
      fflush(stdout);
      iter++;
      joblist=joblist->next;
    }

    if(iter>1) {
      double t;
      size_t pcounter,poutput,t_idx;
      gsl_complex ztmp,ss;

      /* normalize dm_sum and show averaged dm */
      ztmp=gsl_complex_rect(0.0,0.0);
      for(i=0;i<dmdim;i++) {
	ztmp=gsl_complex_add(ztmp,gsl_matrix_complex_get(dm_sum[0],i,i));
      }
      printf("ztmp=%18.12f + %18.12f\n\n\n\n",GSL_REAL(ztmp),GSL_IMAG(ztmp));
      ss=gsl_complex_inverse(ztmp); // use ss to normalize density matrices
      printf("ss=%18.12f + %18.12f\n\n\n\n",GSL_REAL(ss),GSL_IMAG(ss));
      pcounter=0; 
      poutput=(size_t)rint(keys.tprint/keys.tstep);
      t_idx=0;
      // output the averaged dm
      for(t=0.0;t<=keys.tstop;t+=keys.tstep) { // This corresponds to the layout in propagate()
	if((pcounter % poutput) == 0) {
	  gsl_matrix_complex_scale(dm_sum[t_idx],ss);
	  printf("t=%12.4f, rho_avg=",t*TIME_CM2FS);
	  dmop_lvprint(dm_sum[t_idx]);
	  printf("\n");
	  t_idx++;
	}
	pcounter++;
      }
    } // if(iter>1)

    for(i=0;i<narray;i++) {
      gsl_matrix_complex_free(dm_sum[i]); // release memory space too
      gsl_matrix_complex_free(dm_traj[i]);
    }
    free(dm_sum);
    free(dm_traj);
  }

  /* every thing done; clean up */
  gsl_matrix_complex_free(initdm);
  params_close();

  printf("\n");
  printf("Job finished: ");
  print_time_stamp();
  printf("\n");

  return 0;
}

/*
 * $Log$
 * Revision 1.15  2007/07/10 06:16:48  platin
 *   - use a second generation ranlux random number generator.
 *
 * Revision 1.14  2007/06/23 00:54:49  platin
 *
 *   - bump to 0.8p3.
 *   - add support for index adjustment so that poptrans terms can be
 *     consistent when level energy crossover happens at disordered jobs.
 *
 * Revision 1.13  2007/06/20 19:04:25  platin
 *
 *   - time stamp the start of job propagation and the finish of job.
 *   - fflush stdout to make the message more prompt.
 *
 * Revision 1.12  2007/06/18 19:28:16  platin
 *
 *   - add cvs revision tags.
 *
 * Revision 1.11  2007/06/14 18:08:46  platin
 *
 *   - support for RSEED.
 *   - bump to 0.8p2
 *
 * Revision 1.10  2007/06/01 18:00:02  platin
 *
 *   - bump to 0.8p1
 *   - improved the Gaussian-Hermite module; now users can select between
 *     1-20 point rules can be used.
 *   - added keyword to select between MC and GH samplings.
 *   - static disorder average loop in both tnl-dm3pes and tnl-dynamics
 *     is fixed.
 *
 * Revision 1.9  2007/03/26 19:15:25  platin
 *
 *   - bump to v0.8.
 *   - fix the pulse duration factor; the fwhm value is the field fwhm.
 *
 * Revision 1.8  2007/03/15 02:37:05  platin
 *
 *   - add support for the multimode browanian oscillator bath.
 *
 * Revision 1.7  2007/03/06 23:32:01  platin
 *
 *   - remove useless -markov codes.
 *   - use SS initial condition for all bath aux functions.
 *
 * Revision 1.6  2007/02/28 04:16:06  platin
 *
 *   - implement correct phase-locking
 *
 * Revision 1.5  2007/02/20 23:31:13  platin
 *
 *   - hide unnecessary debug messages.
 *
 * Revision 1.4  2007/02/20 22:02:51  platin
 *
 *   - implemented Markov approximation in tnl-kern.c
 *   - tnl-dm3pes-2d now used the Markov model.
 *
 * Revision 1.3  2006/11/08 06:54:11  platin
 *
 *   - support for new ohmart module.
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
 * Revision 1.19  2006/08/26 05:01:17  platin
 *
 *   - include primitive two-exciton state support for dynamics and
 *     optical transitions.
 *
 * Revision 1.18  2006/08/23 22:30:26  platin
 *
 *   - add support for two-exciton states in codes regarding H and
 *     transition dipoles. The "Assign" format for transition dipole
 *     input can also be used to include effects of excited state absorption.
 *
 *   - basic support for TESLIST keyword.
 *
 * Revision 1.17  2006/08/15 23:05:06  platin
 *
 *   - remove ad-hoc fix in tnl-* dynamics; now both programs correctly
 *     handle renormalization term (Hren) and eliminate erroneous dynamics,
 *     such as nonpositive dynamics and non-zero long time coherence terms...
 *
 * Revision 1.16  2006/08/11 07:06:31  platin
 *
 *   - progress commit.
 *
 * Revision 1.15  2006/08/11 06:00:06  platin
 *
 *   - bump to 0.5p2
 *   - now propagates in S-picture, remove interacting picture code.
 *   - add a "Markovian-like" term in bath functions to preserve positivity.
 *
 * Revision 1.14  2006/08/10 16:19:21  platin
 *
 *   - remove useless pulse time code.
 *
 * Revision 1.13  2006/08/03 02:15:27  platin
 *
 *   - fix the pulse period setup; this should handle the FWHM correctly.
 *   - in the js03art module, save some time when gamma is effectively zero.
 *
 * Revision 1.12  2006/08/01 22:25:28  platin
 *
 *   - makeshift modification, use js03art bath instead.
 *
 * Revision 1.11  2006/08/01 20:22:11  platin
 *
 *   - fix bug.
 *
 * Revision 1.10  2006/07/31 08:34:32  platin
 *
 *   - align laser pulse phase as 0 at center of pulse.
 *
 * Revision 1.9  2006/07/31 03:27:33  platin
 *
 *   - pre-compute Sn(t) to avoid repeating evaluations; this gives
 *     a slight performance increase.
 *
 * Revision 1.8  2006/07/31 02:19:52  platin
 *
 *   - progressive commit.
 *
 * Revision 1.7  2006/07/28 23:31:08  platin
 *
 *   - major change. Use short-time Crank-Nicholson propagator.
 *
 * Revision 1.6  2006/07/27 17:54:18  platin
 *
 *   - minor fix.
 *
 * Revision 1.5  2006/07/27 17:48:03  platin
 *
 *   - add support of given initial condition.
 *   - remove restriction on the number of pulses.
 *
 * Revision 1.4  2006/07/20 20:19:47  platin
 *
 *   - minor changes.
 *   - use -DASYMMETRIC_FIELD to use non-Gaussian pulse shape in tnl-dm3pes.
 *
 * Revision 1.3  2006/07/20 18:41:46  platin
 *
 *   - sign problem should be fixed.
 *   - bump version number for tnl-dm3pes to 0.4.
 *
 * Revision 1.2  2006/07/20 17:03:05  platin
 *
 *   - apply an artificial minus sign to correct the dynamics;
 *     this must be fixed later.
 *   - update tnl-dynamics.c.
 *
 * Revision 1.1  2006/07/11 16:46:32  platin
 *
 *   - import the tnl-dynamics bit; still requires modification to take
 *     initial states etc.
 *
 *
 */
