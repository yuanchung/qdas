/***************************************************
 * bath_htgau.h
 * 
 * Header file for bath_htgau.c.
 *
 * By Yuan-Chung Cheng <yuanchung@ntu.edu.tw>
 *
 *
 ***************************************************/

#ifndef _BATH_HTGAU_H
#define _BATH_HTGAU_H 1

#include "qdas.h"

/* Shared global variables */
extern double BathHTGauLambda[256];
extern double BathHTGauTau[256];
extern double BathHTGauFRatio[256];

/* exported functions */

int bath_htgau_init_params(const size_t nsize, const double beta, 
			   const size_t bath_nparams, const double *bath_params);
int bath_htgau_free_params();
double bath_htgau_wsquare(double lambda0,double tau0,double beta);
double bath_htgau_gn_r(double t, double lambda, double tau, double beta);
double bath_htgau_gn_i(double t, double lambda, double tau);
double bath_htgau_gv_r(double t, vibmodes *vibs,double beta);
double bath_htgau_gv_i(double t, vibmodes *vibs);

#endif /* bath_htgau.h */

/*
 * $Log$
 * Revision 1.2  2007/01/30 08:33:15  platin
 *   - use more reasonable coupling constant factor in the bath_htgau module,
 *     note that this is still a high-T approximation.
 *   - more informative messages in 2c3ppes-htgau.c
 *
 * Revision 1.1.1.1  2006/05/24 00:42:18  platin
 *
 *   - initial import of the qdas package.
 *   - qdas stands for quantum dynamics and spectroscopy.
 *   - basic code inherited from the "lineshape" package.
 *
 *
 */
