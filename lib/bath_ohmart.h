/***************************************************
 * bath_ohmart.h
 * 
 * Header file for bath_ohmart.c, used to provide
 * correlation functions for the TNL-* modules.
 *
 * By Yuan-Chung Cheng <yuanchung@ntu.edu.tw>
 *
 *
 ***************************************************/

#ifndef BATH_OHMART_H
#define BATH_OHMART_H 1

#include "qdas.h"
#include "spectrum.h"
#include <gsl/gsl_complex.h>
#include "bathtype_mt99art.h"

/* exported functions */

/* functions provided by this bath module */
int bath_ohmart_init_params(const size_t nsize, const double beta, 
			    const size_t bath_nparams, const double *bath_params);
void bath_ohmart_ct_init(mt99art_bath_func *f, 
			 const double beta, const double gamma, const double wc,
			 const gsl_matrix *S);

#endif /* bath_ohmart.h */

/*
 * $Log$
 * Revision 1.2  2006/11/29 23:15:23  platin
 *   - minor change.
 *
 * Revision 1.1  2006/11/08 06:23:50  platin
 *
 *   - minor changes and adding a artificial ohmic bath module that
 *     uses more efficient artificial fit to Ohmic correlation
 *     functions.
 *
 *
 */
