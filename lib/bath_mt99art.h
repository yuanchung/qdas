/***************************************************
 * bath_mt99art.h
 * 
 * Header file for bath_mt99art.c, used to provide
 * correlation functions and time-domaine lineshape
 * for the YCC06 module.
 *
 * By Yuan-Chung Cheng <yuanchung@ntu.edu.tw>
 *
 *
 ***************************************************/

#ifndef BATH_MT99ART_H
#define BATH_MT99ART_H 1

#include "qdas.h"
#include "spectrum.h"
#include <gsl/gsl_complex.h>
#include "bathtype_mt99art.h"

/* exported functions */

/* functions provided by this bath module */
int bath_mt99art_init_params(const size_t nsize, const double beta, 
			     const size_t bath_nparams, const double *bath_params);

#endif /* bath_mt99art.h */

/*
 * $Log$
 * Revision 1.1  2007/03/15 08:20:50  platin
 *   - add direct support for MT99 artificial bath.
 *
 *
 */
