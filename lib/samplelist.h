/***************************************************
 * samplelist.h
 * 
 * By Yuan-Chung Cheng <yccheng@mit.edu>
 *
 * Header file for samplelist.c.
 *
 ***************************************************/

#ifndef _SAMPLELIST_H
#define _SAMPLELIST_H 1

#include "qdas.h"

/* types */

/* linked-list structure for sampling jobs */
typedef struct {
  gsl_matrix *H; /* Hamilton */
  gsl_matrix *P; /* polarization 4 vectors */
  double weight; /* a weight for this sample; used in non-uniform sampling methods */
  void *next;
} qdas_sample_item;

/* Macros */

/* exported functions */
qdas_sample_item * prepare_job_list(const qdas_keys *keys);
int job_list_count(qdas_sample_item *list);


#endif /* samplelist.h */

/*
 * $Log$
 */
