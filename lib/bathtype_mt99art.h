/***************************************************
 * bathtype_mt99art.h
 * 
 * Header file for bathtype_mt99art.c, used to provide
 * correlation functions for the TNL-* modules.
 *
 * By Yuan-Chung Cheng <yuanchung@ntu.edu.tw>
 *
 *
 ***************************************************/

#ifndef BATHTYPE_MT99ART_H
#define BATHTYPE_MT99ART_H 1

#include "qdas.h"
#include "spectrum.h"
#include <gsl/gsl_complex.h>

/* exported functions */

/* Meier and Tannor's parameterization of bath correlation functions
   C(t) = a(t) - ib(t)
   where a(t) and b(t) are represented by superposition of 
   exponentials:
   a(t) = \sum_{k=1..Nr} Ar[k]*exp(Gr[k]*t)
   b(t) = \sum_{k=1..Ni} Ai[k]*exp(Gi[k]*t)
*/
struct mt99art_bath_func_struct 
{

  double beta; /* inverse temperature */

  size_t Nr; /* number of terms in a(t) */
  size_t Ni; /* number of terms in b(t) */
  double Vren; /* renormalization term Vren = 2/pi * \int_0^infinity J(w)/w * dw */
  
  gsl_complex Ar[256]; /* coefficients for a(t), at most 256 terms */
  gsl_complex Gr[256];
  gsl_complex Ai[256]; /* coefficients for b(t), at most 256 terms */
  gsl_complex Gi[256];

  gsl_matrix *S; /* the system Op of the system-bath interacting term */
};

typedef struct mt99art_bath_func_struct mt99art_bath_func;

/* Global variables in the MT99ART bath module */
extern size_t BATH_MT99ART_NOps; /* number of bath operator terms */
extern mt99art_bath_func BATH_MT99ARTBathFunc[];
extern gsl_matrix_complex *BATH_MT99ARTOpSbar[];

/* functions provided by this bath module */
gsl_complex bath_mt99art_ct(mt99art_bath_func *f, double t);
double bath_mt99art_gt_r(mt99art_bath_func *f, double t);
double bath_mt99art_gt_i(mt99art_bath_func *f, double t);
double bath_mt99art_ht_r(mt99art_bath_func *f, double t);
double bath_mt99art_ht_i(mt99art_bath_func *f, double t);
void bath_mt99art_init_OpSbar(const qdas_keys *keys, const gsl_matrix *U);
void bath_mt99art_free_OpSbar();
int bath_mt99art_free_params();

/* and now include the bath modules in this type */
#include "bath_mt99ohm.h"
#include "bath_js03art.h"
#include "bath_ohmart.h"

#endif /* bathtype_mt99art.h */

/*
 * $Log$
 * Revision 1.3  2007/03/12 23:08:27  platin
 *   - fixed MT99 bath modules to allow more general bath functions,
 *     now baths not defined by gamma and wc can be used.
 *
 * Revision 1.2  2006/11/08 06:54:35  platin
 *
 *   - add ohmart in mt99art bathtype.
 *
 * Revision 1.1  2006/11/02 19:29:31  platin
 *
 *   - implement the bathtype interface.
 *
 *
 */
