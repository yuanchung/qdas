/***************************************************
 * bath_js03art.h
 * 
 * Header file for bath_js03art.c, used to provide
 * correlation functions for the TNL-* modules.
 *
 * By Yuan-Chung Cheng <yuanchung@ntu.edu.tw>
 *
 *
 ***************************************************/

#ifndef BATH_JS03ART_H
#define BATH_JS03ART_H 1

#include "qdas.h"
#include "spectrum.h"
#include <gsl/gsl_complex.h>
#include "bathtype_mt99art.h"

/* exported functions */

/* functions provided by this bath module */
int bath_js03art_init_params(const size_t nsize, const double beta, 
			     const size_t bath_nparams, const double *bath_params);

#endif /* bath_js03art.h */

/*
 * $Log$
 * Revision 1.4  2006/11/02 19:29:31  platin
 *   - implement the bathtype interface.
 *
 * Revision 1.3  2006/10/27 23:02:59  platin
 *
 *   - commit modifications up to 10/27/2006.
 *
 * Revision 1.2  2006/08/26 05:01:17  platin
 *
 *   - include primitive two-exciton state support for dynamics and
 *     optical transitions.
 *
 * Revision 1.1  2006/08/01 19:10:51  platin
 *
 *   - include artificial JS03 bath model that uses four terms for
 *     M&T's parameterization of bath correlation functions.
 *
 *
 */
