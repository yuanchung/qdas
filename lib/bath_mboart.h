/***************************************************
 * bath_mboart.h
 * 
 * Header file for bath_mboart.c, used to provide
 * correlation functions and time-domaine lineshape
 * for the YCC06 module.
 *
 * By Yuan-Chung Cheng <yuanchung@ntu.edu.tw>
 *
 *
 ***************************************************/

#ifndef BATH_MBOART_H
#define BATH_MBOART_H 1

#include "qdas.h"
#include "spectrum.h"
#include <gsl/gsl_complex.h>
#include "bathtype_mt99art.h"

/* exported functions */

/* functions provided by this bath module */
int bath_mboart_init_params(const size_t nsize, const double beta, 
			    const size_t bath_nparams, const double *bath_params);

#endif /* bath_mboart.h */

/*
 * $Log$
 * Revision 1.1  2007/03/15 02:36:26  platin
 *   - add multi-mode brownian oscillator model.
 *
 *
 */
