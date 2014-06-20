/***************************************************
 * bath_js03.h
 * 
 * Header file for bath_js03.c, used to provide
 * correlation functions and time-domaine lineshape
 * for the YCC06 module.
 *
 * By Yuan-Chung Cheng <yuanchung@ntu.edu.tw>
 *
 *
 ***************************************************/

#ifndef BATH_JS03_H
#define BATH_JS03_H 1

#include "qdas.h"
#include "spectrum.h"
#include <gsl/gsl_complex.h>

/* exported functions */

/* data structure to hold a "bath function" */
struct js03_bath_func_struct 
{
  double gamma; /* gamma is a scaling factor */
  double beta; /* inverse temperature */
  double wc; /* cut-off frequency */
  gsl_spline *BATH_JS03_Ct_r_Spline;
  gsl_interp_accel *BATH_JS03_Ct_r_Spline_Acc;
  gsl_spline *BATH_JS03_Ct_i_Spline;
  gsl_interp_accel *BATH_JS03_Ct_i_Spline_Acc;
};

typedef struct js03_bath_func_struct js03_bath_func;

/* data structure to hold a "bath lineshape function" */
struct js03_lineshape_func_struct 
{
  double gamma; /* gamma is a scaling factor */
  js03_bath_func *bf;
  gsl_spline *BATH_JS03_gt_r_Spline;
  gsl_interp_accel *BATH_JS03_gt_r_Spline_Acc;
  gsl_spline *BATH_JS03_gt_i_Spline;
  gsl_interp_accel *BATH_JS03_gt_i_Spline_Acc;
};

typedef struct js03_lineshape_func_struct js03_lineshape_func;

/* Global variables in the JS03 bath module */
extern int BATH_JS03OpBra[];
extern int BATH_JS03OpKet[];
extern js03_bath_func BATH_JS03BathFunc[];

/* functions provided by this bath module */
int bath_js03_init_params(const size_t nsize, const double beta, 
			  const size_t bath_nparams, const double *bath_params);
int bath_js03_free_params();
double bath_js03_ct_r_cached(js03_bath_func *f, double t);
double bath_js03_ct_i_cached(js03_bath_func *f, double t);
double bath_js03_gt_r_cached(js03_lineshape_func *f, double t);
double bath_js03_gt_i_cached(js03_lineshape_func *f, double t);
double bath_js03_gv_r(double t, vibmodes *vibs,double beta);
double bath_js03_gv_i(double t, vibmodes *vibs);
void bath_js03_init_OpQ(int ndim, gsl_matrix *U);
void bath_js03_free_OpQ();
void bath_js03_Qt_r(gsl_vector *qr, double t);
void bath_js03_Qt_i(gsl_vector *qi, double t);
void bath_js03_ct_init(js03_bath_func *f, double beta, double gamma, double wc);
void bath_js03_gt_init(js03_lineshape_func *f, js03_bath_func *bf);

#endif /* bath_js03.h */

/*
 * $Log$
 * Revision 1.5  2006/11/28 01:24:15  platin
 *   - add parameter for ohmart bath at wc=100 and T=298K.
 *
 * Revision 1.4  2006/10/31 22:17:12  platin
 *
 *   - fix negative time part in cached g(t).
 *   - save unscaled g(t) in cache file too.
 *
 * Revision 1.3  2006/10/31 19:35:32  platin
 *
 *   - add support for lineshape function g(t) in the bath_js03 module.
 *
 * Revision 1.2  2006/05/26 19:19:30  platin
 *
 *   - add dynamics module.
 *   - revise the params.c, use a more reasonable scheme to handle
 *     parameters for different modules.
 *
 * Revision 1.1.1.1  2006/05/24 00:42:19  platin
 *
 *   - initial import of the qdas package.
 *   - qdas stands for quantum dynamics and spectroscopy.
 *   - basic code inherited from the "lineshape" package.
 *
 *
 *
 */
