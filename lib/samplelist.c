/***************************************************
 * samplelist.c
 * 
 * By Yuan-Chung Cheng <yuanchung@ntu.edu.tw>
 *
 * Part of the qdas package;
 * Operations related to the handeling of 
 * static disorder samplings...
 *
 ***************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "params.h"
#include "aux.h"
#include "gauss-hermite.h"
#include "linpolar.h"
#include "samplelist.h"

/* this function prepared a list of Hamiltonians and polarizations
   for sampling runs; the number of items in the list will be returned in nitems */
qdas_sample_item * prepare_job_list(const qdas_keys *keys)
{
    int i;
    size_t nsd;
    size_t iter;
    qdas_sample_item *joblist;
    qdas_sample_item **this;

    gsl_rng *r_sd=gsl_rng_alloc(gsl_rng_ranlxd2);
    gsl_rng *r_pol=gsl_rng_alloc(gsl_rng_ranlxd2);    
    if(keys->polm_mc_rseed)
      gsl_rng_set (r_pol,keys->polm_mc_rseed);
    else
      gsl_rng_set (r_pol,(unsigned)time(NULL));
    if(keys->sdm_mc_rseed)
      gsl_rng_set (r_sd,keys->sdm_mc_rseed);
    else
      gsl_rng_set (r_sd,(unsigned)time(NULL));

    // how should we propagate and perhaps average trajectories?
    nsd=count_static_disorder_terms(keys->disorder);    
    if(nsd==0 && keys->polar == NULL) {
      // no sampling is needed...
      joblist=(qdas_sample_item *)malloc(sizeof(qdas_sample_item));
      joblist->weight=1.0;
      joblist->H=gsl_matrix_alloc(keys->nsize,keys->nsize);
      gsl_matrix_memcpy(joblist->H,keys->He);
      joblist->P=NULL;
      joblist->next=NULL;
    } else if (nsd == 0 ) {
      /* no static disorder but we need to sample pulse polarizations for an isotropic sample */
      printf("\n");
      printf("Use Monte-Carlo sampling to account for field polarization factors.\n");
      printf("Will use totally %lu iterations.\n",keys->polm_mc_niter);
      printf("\n");
      this=&joblist;
      for(i=0;i<keys->polm_mc_niter;i++) {
	*this=(qdas_sample_item *)malloc(sizeof(qdas_sample_item));
	(*this)->weight=1.0;
	(*this)->H=gsl_matrix_alloc(keys->nsize,keys->nsize);
	gsl_matrix_memcpy((*this)->H,keys->He);
	(*this)->P=linpolar_4vecs_alloc();
	linpolar_4vecs_set((*this)->P,keys->polar->angle);
	/* random rotation on the field polarization vectors */
	linpolar_4vecs_ranrot((*this)->P,keys->align,r_pol);
	(*this)->next=NULL;
	this=(qdas_sample_item **)&((*this)->next);
      }

      // now we handle static disorder in the following; note that if static disorder is on
      // the polm_mc_niter will be disregarded and the number of iterations follows the
      // sd sampling...

    } else if (keys->sdmethod == QDAS_SDMETHOD_GH) {
      // use Gauss-Hermite integration

      gh_sdisorder_series *ghs;
      ghs=gh_sdisorder_series_alloc(keys->disorder,keys->sdm_gh_order);
      gh_sdisorder_series_set(ghs,keys);

      printf("\n");
      printf("Use Gauss-Hermite %lu-point rule to account for static disorder.\n",keys->sdm_gh_order);
      if(keys->polar)
	printf("Also use Monte-Carlo sampling for field polarization factors concurrently.\n");
      printf("Will use at most %lu iterations.\n",ghs->npoints);
      printf("\n");
#ifdef DEBUG
      gh_sdisorder_series_print(ghs);      
#endif
      this=&joblist;
      for(iter=0;iter<ghs->npoints;iter++) {
	double weight;
	weight=ghs->weights[iter];
	/* FIXME: move this to the gauss-hermite module! */
	if(weight>keys->sdm_gh_wmin) { // ignore small terms
	  *this=(qdas_sample_item *)malloc(sizeof(qdas_sample_item));
	  (*this)->weight=ghs->weights[iter];
	  (*this)->H=gsl_matrix_alloc(keys->nsize,keys->nsize);
	  gsl_matrix_memcpy((*this)->H,keys->He);
	  gsl_matrix_add((*this)->H,ghs->deltaH[iter]);
	  if(keys->polar) {
	    (*this)->P=linpolar_4vecs_alloc();
	    linpolar_4vecs_set((*this)->P,keys->polar->angle);
	    /* random rotation on the field polarization vectors */
	    linpolar_4vecs_ranrot((*this)->P,keys->align,r_pol);
	  } else {
	    (*this)->P=NULL;
	  }
	  (*this)->next=NULL;
	  this=(qdas_sample_item **)&((*this)->next);
	}
      }
      gh_sdisorder_series_free(ghs);
    } else if(keys->sdmethod == QDAS_SDMETHOD_MC) {
      // use Monte-Carlo sampling

      printf("\n");
      printf("Use Monte-Carlo sampling to account for static disorder.\n");
      if(keys->polar)
	printf("Also use Monte-Carlo sampling for field polarization factors concurrently.\n");
      printf("Will use totally %lu iterations.\n",keys->sdm_mc_niter);
      printf("\n");

      this=&joblist;
      for(iter=0;iter<keys->sdm_mc_niter;iter++) {
	  *this=(qdas_sample_item *)malloc(sizeof(qdas_sample_item));
	  (*this)->weight=1.0;
	  (*this)->H=gsl_matrix_alloc(keys->nsize,keys->nsize);
	  construct_Hamiltonian((*this)->H,keys,r_sd);
	  if(keys->polar) {
	    (*this)->P=linpolar_4vecs_alloc();
	    linpolar_4vecs_set((*this)->P,keys->polar->angle);
	    /* random rotation on the field polarization vectors */
	    linpolar_4vecs_ranrot((*this)->P,keys->align,r_pol);
	  } else {
	    (*this)->P=NULL;
	  }
	  (*this)->next=NULL;
	  this=(qdas_sample_item **)&((*this)->next);
      }
    } // done preparation of samplings

    // DONE
    gsl_rng_free(r_sd);
    gsl_rng_free(r_pol);

    return joblist;
}

int job_list_count(qdas_sample_item *list)
{
  int count=0;
  while(list) {
    count++;
    list=list->next;
  }
  return count;
}


/*
 * $Log$
 */
