/***************************************************
 * tnl-dm3pes.c
 * 
 * By Yuan-Chung Cheng <yuanchung@ntu.edu.tw>
 *
 * Simulating three-pulse photon-echo signal using
 * Domcke's density-matrix based method.
 * Domcke's method in JCP 123, 164112 (2005).
 * Update and use the more efficient version in   
 * Egorova et al., JCP 126, 074314 (2007).
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

#include "tnl-dm3pes.h"
#include "params.h"
#include "aux.h"
#include "dmop.h"
#include "bath.h"
#include "bathtype_mt99art.h"
#include "tnl-kern.h"
#include "gauss-hermite.h"
#include "linpolar.h"
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
calculate three-pulse photon-echo signals.\n\
\n\
OPTIONS:\n\
\n\
    -nocache  Do not cache polarizability data.\n\
    -resume   Resume calculation if a .qcache file is found.\n",
	  program_name);
}


/* This function applys the detection window in keys and
   save the screened dipole matrix in mX_new. mX_new
   can be used to generated screened polarization output. */
void apply_detection_window(gsl_matrix *mX_new,
			    const gsl_matrix *mX,
			    const gsl_vector *lambda,
			    const qdas_keys *keys)
{
  double center,halfwidth,Eab,delta;
  size_t a,b,ndim;

  center=(keys->detect.upper+keys->detect.lower)/2.0;
  halfwidth=(keys->detect.upper-keys->detect.lower)/2.0;
  ndim=mX->size1;

  for(a=0;a<ndim;a++) {
    for(b=0;b<ndim;b++) {
      /* mX(a,b) denotes the transition dipole from b->a */
      Eab=gsl_vector_get(lambda,b)-gsl_vector_get(lambda,a);
      delta=fabs(fabs(Eab)-fabs(center));
      if(delta<halfwidth) {
	gsl_matrix_set(mX_new,a,b,
		       gsl_matrix_get(mX,a,b));
      } else {
      /* we apply a squared detection window by eliminating all
	 transitions outside the window. */
	gsl_matrix_set(mX_new,a,b,0.0);
      }
    }
  }

#ifdef DEBUG
  printf("Detection mX = \n");
  gsl_matrix_print(mX_new);
  printf("\n");
#endif

}

/* build V(t) = l1 * V1_bar(t) + l2 * V2^\dagger_bar(t) + l3 * V3^\dagger_bar(t)
   where Vi_bar(t) = E_i(t-tau_i) * \sum_a exp{I*(w_i-Ea+Eg)*t} * V_{ga} * |g><a|.
   See p. 73 of Notebook 1. */
void build_V_matrix(gsl_matrix_complex *V, int l1, int l2, int l3,
		    const double t, 
		    gsl_matrix *mX_p1, gsl_matrix *mX_p2, gsl_matrix *mX_p3,
		    const qdas_keys *keys)
{
  int i,a,b;
  gsl_complex ztmp1,ztmp2;
  double E0,tau0,tau_fwhm,w0;
  double Vab,field;
  size_t ndim;
  int ll[3];

  gsl_matrix *mX[3];
  mX[0]=mX_p1;
  mX[1]=mX_p2;
  mX[2]=mX_p3;

  ll[0]=l1;ll[1]=l2;ll[2]=l3; // use to decide which term to pickup

  ndim=V->size1;

  gsl_matrix_complex_set_zero(V);

  for(i=0;i<3;i++) { // three pulses
    if(ll[i] == 1) {
      E0=keys->pulse_seq[i].E0;
      tau0=keys->pulse_seq[i].tau0;
      tau_fwhm=keys->pulse_seq[i].fwhm_tau;
      w0=keys->pulse_seq[i].w0;
      /* Gaussian field is default */
      if(fabs(t-tau0) > 2.0*tau_fwhm) continue; // only use (+-) 2*FWHM pulse
      // FWHM of the field
      field=E0*exp(-2.0*1.385*(t-tau0)*(t-tau0)/tau_fwhm/tau_fwhm);
      /* consider full transitions, not only |g>->|e> transitions */
      for(a=0;a<ndim;a++) {
	for(b=0;b<ndim;b++) { // Vag=0 when a=b, we assume zero perment dipole.
	  Vab=gsl_matrix_get(mX[i],a,b);
	  // t-tau0 phase term is right, should also work for 3PPEPS
	  ztmp1=gsl_complex_polar(1.0,w0*(t-tau0)); 
	  ztmp2=gsl_complex_mul_real(ztmp1,field*Vab);
	  /* V1_bar term */
	  if(i==0) {
	    ztmp1=gsl_matrix_complex_get(V,a,b);
	    gsl_matrix_complex_set(V,a,b,gsl_complex_add(ztmp1,ztmp2));
	  } else {
	    /* V2^\dagger_bar and V3^\dagger_bar term */
	    ztmp1=gsl_matrix_complex_get(V,b,a);
	    gsl_matrix_complex_set(V,b,a,gsl_complex_add(ztmp1,gsl_complex_conjugate(ztmp2)));
	  }
	}
      } // a
    } // if ll[i] == 1
  }
  // done
}

// compute polarization for aux density matrices, see p. 74 */
// P3(t) = Tr{ X *(rho1-rho2-rho3)}
// we take the ks= k3 + k2 - k1 signal, and ignore the -ks term (which is P(ks,t)^*)
// also note that we use mtmp to calculate the third order density matrix
gsl_complex compute_polarization(tnl_aux_funcs *tfs[],
				 gsl_matrix *mX_det, gsl_matrix_complex *mtmp)
{
  double Vab;
  gsl_complex ztmp,sigma_ba;
  gsl_complex sum;
  int a,b;
  size_t ndim;

  ndim=mX_det->size1;

  // mtmp = rho1-rho2-rho3
  gsl_matrix_complex_memcpy(mtmp,tfs[0]->sigma);
  gsl_matrix_complex_sub(mtmp,tfs[1]->sigma);
  gsl_matrix_complex_sub(mtmp,tfs[2]->sigma);

  /* Tr{X*rho} */
  sum=gsl_complex_rect(0.0,0.0);
  for(a=0;a<ndim;a++) {
    for(b=0;b<ndim;b++) {
      Vab=gsl_matrix_get(mX_det,a,b);
      sigma_ba=gsl_matrix_complex_get(mtmp,b,a);
      ztmp=gsl_complex_mul_real(sigma_ba,Vab);
      sum=gsl_complex_add(sum,ztmp);
    }
  }

  return sum;
}

/* constuct ground state in the exciton basis. this is not 
   redundent because the excited state can have negative 
   (compared to groundstate) transition frequency after RWA.
   So, the lowest eigenstate is not necessarily the groundstate.
   Rather, the first state in the site basis is the ground state, by convention. */
void construct_groundstate(gsl_matrix_complex *dm, gsl_matrix *U)
{
  gsl_matrix *dm0;
  gsl_matrix *mtmp;
  size_t ndim;

  ndim=dm->size1;
  dm0=gsl_matrix_alloc(ndim,ndim);
  mtmp=gsl_matrix_alloc(ndim,ndim);
  gsl_matrix_set_zero(dm0);
  gsl_matrix_set(dm0,0,0,1.0); // ground state in the site representation

  /* transform to exciton basis using U^\dagger*dm*U */
  gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,U,dm0,0.0,mtmp);
  gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,mtmp,U,0.0,dm0);

  /* copy to the complex matrix */
  gsl_matrix_complex_copy_real(dm,dm0);

  gsl_matrix_free(dm0);
  gsl_matrix_free(mtmp);

}

/* Use density-matrix based method to propagate the 3rd-order 
   polarization in the eigenbasis; the dissipation part is 
   implemented using the QME method in Notebook 1, page 65-69;
   if a complex vector Pt_traj is given, the trajectory of P(t) will be saved */
void propagate(gsl_vector_complex *Pt_traj,
	       const double tstop, const double tstep, 
	       qdas_keys *keys, const gsl_matrix *H, const linpolar_4vecs *Pol, int Jcount)
{
  gsl_matrix *U=gsl_matrix_alloc(keys->nsize,keys->nsize);
  gsl_matrix *mX_x=gsl_matrix_alloc(keys->nsize,keys->nsize);
  gsl_matrix *mX_y=gsl_matrix_alloc(keys->nsize,keys->nsize);
  gsl_matrix *mX_z=gsl_matrix_alloc(keys->nsize,keys->nsize);
  gsl_matrix *mX_abs=gsl_matrix_alloc(keys->nsize,keys->nsize);
  gsl_matrix *mX_p1=gsl_matrix_alloc(keys->nsize,keys->nsize);
  gsl_matrix *mX_p2=gsl_matrix_alloc(keys->nsize,keys->nsize);
  gsl_matrix *mX_p3=gsl_matrix_alloc(keys->nsize,keys->nsize);
  gsl_matrix *mX_p4=gsl_matrix_alloc(keys->nsize,keys->nsize);
  gsl_matrix *mX_det=gsl_matrix_alloc(keys->nsize,keys->nsize);
  gsl_vector *lambda=gsl_vector_alloc(keys->nsize);

  tnl_aux_funcs *sigma_t[3]; // totally three aux density matrices
  tnl_aux_funcs *sigma_tpdt[3];

  gsl_matrix_complex *V=gsl_matrix_complex_alloc(keys->nsize,keys->nsize);
  gsl_matrix_complex *W_t=gsl_matrix_complex_alloc(keys->nsize,keys->nsize);
  gsl_matrix_complex *Weff_t=gsl_matrix_complex_alloc(keys->nsize,keys->nsize);
  gsl_matrix_complex *Heigen=gsl_matrix_complex_alloc(keys->nsize,keys->nsize);
  gsl_matrix_complex *Hren=gsl_matrix_complex_alloc(keys->nsize,keys->nsize);

  gsl_matrix_complex *groundstate=gsl_matrix_complex_alloc(keys->nsize,keys->nsize);

  gsl_complex Pt;

  size_t ndim;

  double t;
  double t1,t2,t3;
  double dtmp;
  double Sint; // integrated signal
  double tf;
  size_t i,a,t_idx;

  size_t pcounter,poutput;

  // handy variables
  ndim=keys->nsize;

  // pulse time
  t1=keys->pulse_seq[0].tau0;
  t2=keys->pulse_seq[1].tau0;
  t3=keys->pulse_seq[2].tau0;

  // propagate tstop fs longer than the center of the third pulse
  tf=t3+tstop+tstep/10.0;

  /* obtain U/lambda, and mX */
  setup_eigenbasis(U,lambda,mX_x,mX_y,mX_z,mX_abs,keys,H);
 
  if(Pol) {
    /* a set of polarization vectors is given; we need to consider
       polarization of the field and each pulse should have
       a different interaction matrices */
    /* generate mX matrix corresponding to each pulses */
    compute_pulse_interactions(mX_p1,mX_p2,mX_p3,mX_p4,
			       mX_x,mX_y,mX_z, Pol);
    /* apply detection window and save the screened mX_abs in mX_det */
    apply_detection_window(mX_det,mX_p4,lambda,keys);
  } else {
    /* no polarization stuff; use the abs interactions */
    gsl_matrix_memcpy(mX_p1,mX_abs);
    gsl_matrix_memcpy(mX_p2,mX_abs);
    gsl_matrix_memcpy(mX_p3,mX_abs);
    /* apply detection window and save the screened mX_abs in mX_det */
    apply_detection_window(mX_det,mX_abs,lambda,keys);
  }
  
  /* no need for these primitive mXs now, put them away... */
  gsl_matrix_free(mX_x);
  gsl_matrix_free(mX_y);
  gsl_matrix_free(mX_z);
  gsl_matrix_free(mX_p4);
  gsl_matrix_free(mX_abs);
 
  /* initialize bath in the exciton basis */
  bath_mt99art_init_OpSbar(keys,U);

  /* construct the ground state density matrix in the exciton basis */
  construct_groundstate(groundstate,U);

#ifdef DEBUGDEBUG
  printf("\n");
  printf("Groundstate = \n");
  gsl_matrix_complex_print(groundstate);
  printf("\n");
#endif

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
  sigma_t[0]=(tnl_aux_funcs *)malloc(sizeof(tnl_aux_funcs));
  tnl_aux_funcs_init(sigma_t[0],ndim,groundstate,Heigen,Heigen,lambda,U,keys);
  gsl_matrix_complex_set_zero(Hren);
  for(a=0;a<sigma_t[0]->nops;a++) {
    gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,gsl_complex_rect(1.0,0.0), 
                  sigma_t[0]->bfs[a].Sbar, sigma_t[0]->bfs[a].Sbar, 
                  gsl_complex_rect(0.0,0.0), V); // use V as tmp storage
    gsl_matrix_complex_scale(V,gsl_complex_rect(sigma_t[0]->bfs[a].Vren/2.0,0.0));
    gsl_matrix_complex_add(Hren,V);
  }
  tnl_aux_funcs_free(sigma_t[0]);
  free(sigma_t[0]);

  if(keys->printlv >= 3) {
    printf("Hren = \n");
    gsl_matrix_complex_print(Hren);
  }

  /* aux density matrices, initially everything is in the ground state */
  gsl_matrix_complex_memcpy(V,Hren); // again use V as tmp variable to store Weff(0)
  gsl_matrix_complex_add(V,Heigen);
  for(i=0;i<3;i++) {
    sigma_t[i]=(tnl_aux_funcs *)malloc(sizeof(tnl_aux_funcs));
    tnl_aux_funcs_init(sigma_t[i],ndim,groundstate,V,Heigen,lambda,U,keys);
    sigma_tpdt[i]=(tnl_aux_funcs *)malloc(sizeof(tnl_aux_funcs));
    tnl_aux_funcs_init(sigma_tpdt[i],ndim,groundstate,V,Heigen,lambda,U,keys);
  }

  // prograte in time
  printf("\n");
  printf("Propagate the 3PE polarization from t=%.2f to t=%.2f.\n",-1.0*t3*TIME_CM2FS,(tf-t3)*TIME_CM2FS);
  printf("\n");
  printf("Note that the value of P(t) shown is scaled by 1e6.\n");
  printf("\n");

  /* propagate from time zero, note that the 
     first pulse should be in the positive time.
     FIXME: remove these restrictions!! */
  Sint=0.0;
  pcounter=0; 
  poutput=(size_t)rint(keys->tprint/tstep);

  //  printf("poutput=%ld\n",poutput);
  t_idx=0;
  for(t=0.0;t<=tf;t+=tstep) {

    tnl_aux_funcs *taf_ptr;

    /* compute polarization */
    Pt=compute_polarization(sigma_t,mX_det,V); // use V as tmp variable
    // V contains the third-order density matrix when returned

    //    if(t>=t3) {
    // compute approximated integrated signal
    dtmp=gsl_complex_abs(Pt)*(1e6); // scale by 1M
    Sint=Sint+dtmp*dtmp*tstep;
    //    }

    // make sure we output the t=t3 point
    if(fabs(t-t3)<1e-6) pcounter=0;

    if(
#ifndef DEBUGDEBUG
       // FIXME: should use a more intuitive method to control output details
       Jcount<=1 &&  // don't output detail trajectory info when average is on, unless in DEBUG mode
#endif
       (pcounter % poutput) == 0) { // print P(t) every poutput runs
      printf("t=%12.4f, P=%20.16f  %20.16f, |P|=%20.16f\n",
	     (t-t3)*TIME_CM2FS,GSL_REAL(Pt)*(1e6),GSL_IMAG(Pt)*(1e6),gsl_complex_abs(Pt)*(1e6));
#ifdef DEBUG
      fflush(stdout);
#endif
#ifdef DEBUGDEBUG
      {
	gsl_complex ztmp;
	int nn;
	// V contains the third-order density matrix
	printf("t=%12.4f, rho3=",(t-t3)*TIME_CM2FS);
	dmop_lvprint(V);
	printf("\n");
	// use W_t and V as tmp space
	gsl_matrix_complex_copy_real(W_t,mX_det);
	for(nn=0;nn<3;nn++) {
	  printf("t=%12.4f, s%1d=",(t-t3)*TIME_CM2FS,nn+1);
	  dmop_lvprint(sigma_t[nn]->sigma);
	  printf("\n");
	  gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,gsl_complex_rect(1.0,0.0), 
			 W_t, sigma_t[nn]->sigma, 
			 gsl_complex_rect(0.0,0.0), V); // use V as tmp storage
	  ztmp=gsl_matrix_complex_trace(V);
	  printf("t=%12.4f, <X*s%1d>= %20.12f + %20.12fI\n",
		 (t-t3)*TIME_CM2FS,nn+1,GSL_REAL(ztmp),GSL_IMAG(ztmp));
	}
	printf("\n");
      }
#endif
    }
    pcounter++;

    // save the sum of P(t) when static disorder sampling is on
    if(Pt_traj) {
      gsl_vector_complex_set(Pt_traj,t_idx,Pt);
      t_idx++;
    }

    /* prepare for the t+dt step */

    /* update aux density matrices */
    build_V_matrix(V,1,1,1,t+tstep,mX_p1,mX_p2,mX_p3,keys);
    /* W_t is the time dependent Hamiltonian, W_t = H-V */
    gsl_matrix_complex_memcpy(W_t,Heigen);
    gsl_matrix_complex_sub(W_t,V);
    /* Weff_t is the normalized effective potential,
       Weff(t) = H - V + \sum mu/2*Sn*Sn */
    gsl_matrix_complex_memcpy(Weff_t,W_t);gsl_matrix_complex_add(Weff_t,Hren);
    advance_tnl_aux_funcs(sigma_tpdt[0],sigma_t[0],Weff_t,W_t,lambda,keys,tstep);

    build_V_matrix(V,1,1,0,t+tstep,mX_p1,mX_p2,mX_p3,keys);
    gsl_matrix_complex_memcpy(W_t,Heigen);
    gsl_matrix_complex_sub(W_t,V);
    gsl_matrix_complex_memcpy(Weff_t,W_t);gsl_matrix_complex_add(Weff_t,Hren);
    advance_tnl_aux_funcs(sigma_tpdt[1],sigma_t[1],Weff_t,W_t,lambda,keys,tstep);

    build_V_matrix(V,1,0,1,t+tstep,mX_p1,mX_p2,mX_p3,keys);
    gsl_matrix_complex_memcpy(W_t,Heigen);
    gsl_matrix_complex_sub(W_t,V);
    gsl_matrix_complex_memcpy(Weff_t,W_t);gsl_matrix_complex_add(Weff_t,Hren);
    advance_tnl_aux_funcs(sigma_tpdt[2],sigma_t[2],Weff_t,W_t,lambda,keys,tstep);

    /* done, ready for next run */
    for(i=0;i<3;i++) {
    taf_ptr=sigma_t[i];
    sigma_t[i]=sigma_tpdt[i];
    sigma_tpdt[i]=taf_ptr;
    }
  }

  printf("\n");
  printf("Int. Signal (approx.): tau=%9.2f, T=%9.2f, S*1e12=%26.16f\n",
	 (t2-t1)*TIME_CM2FS,(t3-t2)*TIME_CM2FS,Sint*TIME_CM2FS);
  printf("\n");

  /* clean up */
  bath_mt99art_free_OpSbar();

  gsl_matrix_free(U);
  gsl_matrix_free(mX_p1);
  gsl_matrix_free(mX_p2);
  gsl_matrix_free(mX_p3);
  gsl_matrix_free(mX_det);
  gsl_vector_free(lambda);
  gsl_matrix_complex_free(V);
  gsl_matrix_complex_free(W_t);
  gsl_matrix_complex_free(Weff_t);
  gsl_matrix_complex_free(Heigen);
  gsl_matrix_complex_free(Hren);
  gsl_matrix_complex_free(groundstate);

  for(i=0;i<3;i++) {
    tnl_aux_funcs_free(sigma_t[i]);
    free(sigma_t[i]);
    tnl_aux_funcs_free(sigma_tpdt[i]);
    free(sigma_tpdt[i]);
  }
  // DONE!
}


/* the following functions handle cache files */
/* read the cache file; if not successful, return 0; *iter=0 too, if not succeful */
int read_ptavg_cache(char *fname, size_t *iter, gsl_vector_complex *Pt) 
{
  FILE *input_file;
  
  input_file=fopen(fname, "r");
  if(input_file == NULL) {
    *iter=0;
    return 0;
  }
  fread(iter, sizeof(size_t), 1, input_file);
  gsl_vector_complex_fread(input_file,Pt);
  fclose(input_file);

  return *iter;
}

void save_ptavg_cache(char *fname, size_t iter, gsl_vector_complex *Pt) 
{
  FILE *output_file;
  
  output_file=fopen(fname, "w");
  if(output_file == NULL) {
    fprintf(stderr, "error while opening \"%s\" for writing: %s \n",
            fname,strerror(errno));
    exit(errno);
  }

  fwrite(&iter, sizeof(size_t), 1, output_file);
  gsl_vector_complex_fwrite(output_file,Pt);
  fclose(output_file);
}


/* Main program */
int main(int argc, char *argv[])
{
  int i;

  char *key_file_name=NULL;
  char *cache_file_name=NULL;

  /* variables for the options */
  int resume=0;
  int cache=1;

  /* main keyword structure */  
  qdas_keys keys;
  
  /* Set program name for messages.  */
  program_name = argv[0]; 
  
  /* HEADER MESSAGE */
  printf("\n");
  printf("%s\n",program_name);
  printf("\n");
  printf("        TNL-DM3PES\n");
//  printf("        SVN revision: %s.\n","$Revision: 259 $");
  printf("        Density-matrix based method for three-pulse\n");
  printf("        photon-echo signals. The dissipative dynamics\n");
  printf("        is simulated by a TNL non-Markovian quantum\n");
  printf("        master equation method.\n");
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
	cache_file_name=(char *)malloc((strlen(key_file_name)+8)*sizeof(char));
	strcpy(cache_file_name,key_file_name);
	strcpy(cache_file_name+strlen(key_file_name),".qcache");
      } else {
	usage();
	exit(EXIT_FAILURE);
      }
    } else {
      if(!strncmp(argv[i],"-resume",7)) {
	/* resume previously interrupted calculation when
	   a cache file is found */
	resume = 1;
      } else if(!strncmp(argv[i],"-nocache",6)) {
	/* do not save cache file. */
	cache = 0;
      } else {
	/* new options go here */
	usage();
	exit(EXIT_FAILURE);
      }
    }
  }

  /* Check if we have a key file , 
     if yes, use the key file to initialize all parameters */
  if(key_file_name != NULL){
    params_init(key_file_name,&keys);
  } else {
    usage();
    exit(EXIT_FAILURE);
  }

  /* Check the dynamics module and 
     print out parameters specifically for this module */
  if(keys.npulses != 3) {
    printf("Number of pulses has to be 3 in a 3PES experiment.\n");
    exit(EXIT_FAILURE);
  }
  if(!keys.ndipoles) {
    printf("Required keyword \"DIPOLE\" not found!\n");
    exit(EXIT_FAILURE);
  }

  printf("-------------\n");
  printf("| Detection |\n");
  printf("-------------\n");
  printf("\n");
  printf("The reported signals include only contributions from\n");
  printf("electronic transitions in the energy window:\n");
  printf("  From %f to %f.\n",
	 keys.detect.lower,keys.detect.upper);
  printf("\n");

  printf("\n");
  printf("\n");
  printf("------------------------\n");
  printf("| Dynamical Simulation |\n");
  printf("------------------------\n");
  printf("\n");
  printf("Propagate density matrix using Crank-Nicholson iterative scheme.\n");
  printf("Iteration convegence threshold is %g.\n",keys.cnconv);
  printf("\n");
  printf("Simulate trajectory using time increment of %f fs.\n",keys.tstep*TIME_CM2FS);
  printf("\n");
  printf("Will output P_3P(t) from t=%.2f to t=%.2f fs,\n",
	 -1.0*(keys.pulse_seq[2].tau0)*TIME_CM2FS,(keys.tstop)*TIME_CM2FS);
  printf("with time step size delta_t=%f fs.\n",keys.tprint*TIME_CM2FS);
  printf("\n");

  /* Parameters read and initialized, ready to do real work */
  {
    gsl_vector_complex *Pt_avg=NULL;
    gsl_vector_complex *Pt_traj=NULL;
    size_t narray;
    qdas_sample_item *joblist=NULL;
    size_t iter;
    size_t JJCOUNT;

    /* use multiple trajectories to account for static disorder may be required,
       so we need to store the sum of Pt trajectories. */
    narray=(size_t)ceil(fabs(keys.tstop+keys.pulse_seq[2].tau0)/keys.tstep)+1;
    Pt_avg=gsl_vector_complex_alloc(narray);
    gsl_vector_complex_set_zero(Pt_avg);
    Pt_traj=gsl_vector_complex_alloc(narray);
    gsl_vector_complex_set_zero(Pt_traj);

    joblist=prepare_job_list(&keys);
    JJCOUNT=job_list_count(joblist);

    iter=0;
    /* attemp to resume from a previous calculation if possible */
    if(resume && read_ptavg_cache(cache_file_name,&iter,Pt_avg) > 0) {
      printf("\n");
      printf("Found cache file %s,\n",cache_file_name);
      printf("will resume the calculation from iteration %lu.\n",iter+1);
      printf("\n");
      printf("\n");
      for(i=0;i<iter;i++)
	joblist=joblist->next; // skip these runs...
    }
    while(joblist) {

      printf("Sample #%lu (weight=%20.16f): ",iter+1,joblist->weight);
      print_time_stamp();
      printf("\n");
#ifdef DEBUG
      printf("\n");
      printf("Hamiltonian = \n");
      gsl_matrix_print(joblist->H);
      printf("\n");
      if(keys.polar) {
	printf("Polarization = \n");
	gsl_matrix_print(joblist->P);
	printf("\n");
      }
#endif
      /* now indeed propagate the dynamics */
      propagate(Pt_traj,keys.tstop,keys.tstep,&keys, joblist->H, joblist->P,JJCOUNT);
#ifdef DEBUGDEBUG
      printf("t=, P=, |P|=\n"); // this is for convenient gnuplot plot
#endif
      gsl_vector_complex_scale(Pt_traj,gsl_complex_rect(joblist->weight,0.0));
      gsl_vector_complex_add(Pt_avg,Pt_traj);
      fflush(stdout);
      if(cache == 1) {
	// this run finished, save the cache file and proceed to the next one
	save_ptavg_cache(cache_file_name,iter+1,Pt_avg);
      }
      iter++;
      joblist=joblist->next;
    }
    
    if(iter>1) {
      /* print out the averaged results */
      double t;
      double Sint,dtmp;
      size_t pcounter,poutput,t_idx;
      gsl_complex Pt;

      // scale the sum to get the average
      gsl_vector_complex_scale(Pt_avg,gsl_complex_rect(1.0/((double)iter),0.0));
      /* show averaged P(t) */
      pcounter=0; 
      poutput=(size_t)rint(keys.tprint/keys.tstep);
      t_idx=0;
      Sint=0.0;
      // This corresponds to the layout in propagate()
      for(t=0.0;t<=(keys.tstop+keys.pulse_seq[2].tau0+keys.tstep/10.0);t+=keys.tstep) { 
	Pt=gsl_vector_complex_get(Pt_avg,t_idx);
	if(fabs(t-keys.pulse_seq[2].tau0)<1e-6) pcounter=0;
	if((pcounter % poutput) == 0) {
	  printf("t=%12.4f, Pavg=%20.16f  %20.16f, |Pavg|=%20.16f\n",
		 (t-keys.pulse_seq[2].tau0)*TIME_CM2FS,GSL_REAL(Pt)*(1e6),GSL_IMAG(Pt)*(1e6),gsl_complex_abs(Pt)*(1e6));
	}
	pcounter++;

	// compute approximated integrated signal
	dtmp=gsl_complex_abs(Pt)*(1e6); // scale by 1M
	Sint=Sint+dtmp*dtmp*keys.tstep;
	t_idx++;
      }
      printf("\n");
      printf("Averaged integrated signal (approx.): tau=%9.2f, T=%9.2f, S*1e12=%32.16f\n",
	     (keys.pulse_seq[1].tau0-keys.pulse_seq[0].tau0)*TIME_CM2FS,
	     (keys.pulse_seq[2].tau0-keys.pulse_seq[1].tau0)*TIME_CM2FS,Sint*TIME_CM2FS);
      printf("\n");

    } // if(iter>1)

    /* clean up */
    gsl_vector_complex_free(Pt_avg);
    gsl_vector_complex_free(Pt_traj);
    
  } // end of ready to do real work

  /* done; clean up */
  params_close();

  printf("\n");
  printf("Job finished: ");
  print_time_stamp();
  printf("\n");

  return 0;
}

/*
 * $Log$
 * Revision 1.30  2007/07/15 04:51:18  platin
 *   - add info on RNG.
 *
 * Revision 1.29  2007-07-10 06:16:48  platin
 *
 *   - use a second generation ranlux random number generator.
 *
 * Revision 1.28  2007/06/23 00:54:49  platin
 *
 *   - bump to 0.8p3.
 *   - add support for index adjustment so that poptrans terms can be
 *     consistent when level energy crossover happens at disordered jobs.
 *
 * Revision 1.27  2007/06/20 19:04:25  platin
 *
 *   - time stamp the start of job propagation and the finish of job.
 *   - fflush stdout to make the message more prompt.
 *
 * Revision 1.26  2007/06/18 19:28:16  platin
 *
 *   - add cvs revision tags.
 *
 * Revision 1.25  2007/06/14 18:08:46  platin
 *
 *   - support for RSEED.
 *   - bump to 0.8p2
 *
 * Revision 1.24  2007/06/04 23:51:21  platin
 *
 *   - adjust groundstate convention. Now tnl-dm3pes can correctly handle
 *     negative excited state energy.
 *
 * Revision 1.23  2007/06/01 18:00:02  platin
 *
 *   - bump to 0.8p1
 *   - improved the Gaussian-Hermite module; now users can select between
 *     1-20 point rules can be used.
 *   - added keyword to select between MC and GH samplings.
 *   - static disorder average loop in both tnl-dm3pes and tnl-dynamics
 *     is fixed.
 *
 * Revision 1.22  2007/05/31 01:29:46  platin
 *
 *   - implemented the more efficient verison in Egorova et al.
 *   - bump to 0.8p1.
 *
 * Revision 1.21  2007/05/31 01:11:34  platin
 *
 *   - more reasonable debug info.
 *
 * Revision 1.20  2007/05/07 23:51:04  platin
 *
 *   - minor changes.
 *
 * Revision 1.19  2007/03/26 19:15:25  platin
 *
 *   - bump to v0.8.
 *   - fix the pulse duration factor; the fwhm value is the field fwhm.
 *
 * Revision 1.18  2007/03/15 02:37:05  platin
 *
 *   - add support for the multimode browanian oscillator bath.
 *
 * Revision 1.17  2007/03/11 17:57:12  platin
 *
 *   - support two-exciton states in Gauss-Hermite...
 *
 * Revision 1.16  2007/03/11 08:20:37  platin
 *
 *   - monor debug info change.
 *
 * Revision 1.15  2007/03/11 00:18:27  platin
 *
 *   - ignore small weighted terms in GH integration.
 *
 * Revision 1.14  2007/03/09 08:59:46  platin
 *
 *   - include Gauss-Hermite Quadrature code for static disorder of two DOF.
 *   - bump version.
 *
 * Revision 1.13  2007/03/06 23:30:51  platin
 *
 *   - merge changes in -2d to main file, the (t-tau0) phase factor
 *     should be correct for all.
 *
 * Revision 1.12  2007/02/28 04:16:06  platin
 *
 *   - implement correct phase-locking
 *
 * Revision 1.11  2007/02/20 22:02:51  platin
 *
 *   - implemented Markov approximation in tnl-kern.c
 *   - tnl-dm3pes-2d now used the Markov model.
 *
 * Revision 1.10  2007/02/20 06:43:03  platin
 *
 *   - simple message that tracks program name.
 *
 * Revision 1.9  2007/02/19 20:25:25  platin
 *
 *  - tryy to fix a round-off error,
 *
 * Revision 1.8  2007/02/19 19:51:30  platin
 *
 *   - bugfix.
 *
 * Revision 1.7  2007/02/19 19:49:20  platin
 *
 *   - fixed t output in averaged signal.
 *
 * Revision 1.6  2007/01/19 20:28:15  platin
 *
 *   - change def. of t in dm3pes according to normal convention.
 *
 * Revision 1.4  2006/12/20 07:07:22  platin
 *
 *   - remove redudant code.
 *   - add options of using nonlinear poptrans.
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
 * Revision 1.21  2006/08/26 05:01:17  platin
 *
 *   - include primitive two-exciton state support for dynamics and
 *     optical transitions.
 *
 * Revision 1.20  2006/08/23 22:30:26  platin
 *
 *   - add support for two-exciton states in codes regarding H and
 *     transition dipoles. The "Assign" format for transition dipole
 *     input can also be used to include effects of excited state absorption.
 *
 *   - basic support for TESLIST keyword.
 *
 * Revision 1.19  2006/08/15 23:05:06  platin
 *
 *   - remove ad-hoc fix in tnl-* dynamics; now both programs correctly
 *     handle renormalization term (Hren) and eliminate erroneous dynamics,
 *     such as nonpositive dynamics and non-zero long time coherence terms...
 *
 * Revision 1.18  2006/08/11 07:05:36  platin
 *
 *   - bump to 0.5p1
 *   - not using interaction picture anymore; add ad-hoc correction to
 *     the positivity problem.
 *
 * Revision 1.17  2006/08/03 02:15:27  platin
 *
 *   - fix the pulse period setup; this should handle the FWHM correctly.
 *   - in the js03art module, save some time when gamma is effectively zero.
 *
 * Revision 1.16  2006/08/01 22:27:39  platin
 *
 *   - makeshift change, use js03art bath module instead; should fix
 *     this and use general bath modules...
 *
 * Revision 1.15  2006/08/01 20:22:11  platin
 *
 *   - fix bug.
 *
 * Revision 1.14  2006/07/31 08:36:43  platin
 *
 *   - align laser pulse phase as 0 at center of pulse.
 *   - collect total signal, instead of t>t3, to reflect expr. condition.
 *   - use Crank-Nicholson split-operator short-time propagation scheme.
 *
 * Revision 1.13  2006/07/31 02:19:52  platin
 *
 *   - progressive commit.
 *
 * Revision 1.12  2006/07/20 22:34:43  platin
 *
 *   - adjust coefficient for asymmetric pulse form.
 *
 * Revision 1.11  2006/07/20 20:19:47  platin
 *
 *   - minor changes.
 *   - use -DASYMMETRIC_FIELD to use non-Gaussian pulse shape in tnl-dm3pes.
 *
 * Revision 1.10  2006/07/20 18:41:46  platin
 *
 *   - sign problem should be fixed.
 *   - bump version number for tnl-dm3pes to 0.4.
 *
 * Revision 1.9  2006/07/20 17:03:05  platin
 *
 *   - apply an artificial minus sign to correct the dynamics;
 *     this must be fixed later.
 *   - update tnl-dynamics.c.
 *
 * Revision 1.8  2006/07/15 00:02:10  platin
 *
 *   - use |mu| of excited states, instead of mu_x. Note the orientation
 *     factor is still ignored.
 *
 * Revision 1.7  2006/07/13 23:37:29  platin
 *
 *   - add keyword TPRINT that adjusts time step of output.
 *   - minor message changes.
 *
 * Revision 1.6  2006/07/13 23:07:46  platin
 *
 *   - revert to the "good" short-time propagator in interaction picture.
 *
 * Revision 1.4  2006/07/13 19:49:04  platin
 *
 *   - instead output the real part of Pt, output the full complex Pt
 *     associated with exp(i*ks*r), where ks= k3 + k2 - k1.
 *
 * Revision 1.3  2006/07/11 16:46:32  platin
 *
 *   - import the tnl-dynamics bit; still requires modification to take
 *     initial states etc.
 *
 * Revision 1.2  2006/07/11 00:01:08  platin
 *
 *   - corrected sign problem in tnl-dm3pes.
 *   - implemented density matrix dynamics in tnl-dm3pes; plan to
 *     move it to tnl-dynamics module
 *
 * Revision 1.1  2006/07/10 20:25:41  platin
 *
 *   - a possible sign problem in dm3pes.c,
 *   - fist import of the tnl-dm2pes module, which calculates 3P PE signal
 *     using non-Markovian TNL method from Meier and Tannor.
 *
 *
 */
