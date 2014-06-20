/***************************************************
 * bath_mt99ohm.h
 * 
 * Header file for bath_mt99ohm.c, used to provide
 * correlation functions and time-domaine lineshape
 * for the YCC06 module.
 *
 * By Yuan-Chung Cheng <yuanchung@ntu.edu.tw>
 *
 *
 ***************************************************/

#ifndef BATH_MT99OHM_H
#define BATH_MT99OHM_H 1

#include "qdas.h"
#include "spectrum.h"
#include <gsl/gsl_complex.h>
#include "bathtype_mt99art.h"

/* exported functions */

/* functions provided by this bath module */
int bath_mt99ohm_init_params(const size_t nsize, const double beta, 
			     const size_t bath_nparams, const double *bath_params);

#endif /* bath_mt99ohm.h */

/*
 * $Log$
 * Revision 1.4  2006/11/02 19:29:31  platin
 *   - implement the bathtype interface.
 *
 * Revision 1.3  2006/11/01 23:41:54  platin
 *
 *   - modify the mt99ohm bath model to fit into the current
 *     tnl-kernel structure.
 *
 * Revision 1.2  2006/07/15 07:47:39  platin
 *
 *   - add support for g(t) and h(t).
 *   - add a PI factor to C(t) return value; i.e. use a definition of C(t)
 *     that is different from M&T.
 *
 * Revision 1.1  2006/07/11 17:22:42  platin
 *
 *   - MT99OHM bath model.
 *
 *
 */
