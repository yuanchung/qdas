/***************************************************
 * gauss-hermite.h
 * 
 * By Yuan-Chung Cheng <yuanchung@ntu.edu.tw>
 *
 * Header file for gauss-hermite.c.
 *
 ***************************************************/

#ifndef _GAUSS_HERMITE_H
#define _GAUSS_HERMITE_H 1

#include "qdas.h"
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

/* Global variables */

/* Macros */

/* Typedef */
/* data structure to hold a Gauss-Hermite static disorder series */
struct gh_sdisorder_series_struct 
{
  size_t order; /* n point Gauss-Hermite rule */
  double *xtab; /* abscissa points and weights for this GH rule */
  double *wtab;

  size_t npoints; /* number of points */
  gsl_matrix **deltaH; /* deltaH for each instance of static disorder */
  double *weights; /* weights for each deltaH */
};

typedef struct gh_sdisorder_series_struct gh_sdisorder_series;

/* exported functions */
gh_sdisorder_series *gh_sdisorder_series_alloc(gsl_matrix *disorder, size_t nrule);
void gh_sdisorder_series_free(gh_sdisorder_series *ghs);
void gh_sdisorder_series_set(gh_sdisorder_series *ghs, const qdas_keys *keys);
void gh_sdisorder_series_print(gh_sdisorder_series *ghs);

#endif /* gauss-hermite.h */

/*
 * $Log$
 * Revision 1.4  2007/06/01 17:58:55  platin
 *   - improved the Gaussian-Hermite module; now users can select between
 *     1-20 point rules can be used.
 *   - added keyword to select between MC and GH samplings.
 *
 * Revision 1.3  2007/03/11 17:55:28  platin
 *
 *   - support two-exciton states in Gauss-Hermite deltaH.
 *
 * Revision 1.2  2007/03/09 08:58:56  platin
 *
 *   - more supporting gsl_* functions.
 *   - remove obsoleted code in params.c.
 *   - minor changes in gauss-hermite*
 *
 * Revision 1.1  2007/03/09 06:01:52  platin
 *
 *   - import new Gauss-Hermite Quadrature module that can be used
 *     for static disorder average.
 *
 *
 */
