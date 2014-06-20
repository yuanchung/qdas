/***************************************************
 * spectrum2d.h
 * 
 * Header file for spectrum2d.c, provides definition
 * and interfaces of 2D spectrum2d functions.
 *
 * By Yuan-Chung Cheng <yuanchung@ntu.edu.tw>
 *
 ***************************************************/

#ifndef _SPECTRUM2D_H
#define _SPECTRUM2D_H 1

#include <gsl/gsl_spline.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_matrix.h>

/* data structure to hold a spectrum2d */

struct spectrum2d_struct 
{
  /* number of points in the spectrum2d */
  size_t size1;
  size_t size2;
  double *w1;
  double *w2;
  gsl_matrix_complex *data;
};

typedef struct spectrum2d_struct spectrum2d;

/* exported functions */
spectrum2d *spectrum2d_alloc(double w1_start, double w1_end, double w1_step,
			     double w2_start, double w2_end, double w2_step);
void spectrum2d_free(spectrum2d *spec);
spectrum2d *spectrum2d_dup(const spectrum2d *spec);
spectrum2d *spectrum2d_set_zero(spectrum2d *spec);
void spectrum2d_copy(spectrum2d *dest, const spectrum2d *source);
void spectrum2d_normalize(spectrum2d *spec);
void spectrum2d_scale(spectrum2d *spec, gsl_complex s);
void spectrum2d_add(spectrum2d *spec1, const spectrum2d *spec2);
void spectrum2d_print(spectrum2d *spec);

#endif /* spectrum2d.h */

/*
 * $Log$
 * Revision 1.2  2007/06/18 20:41:53  platin
 *   - minor changes; add the normalize function.
 *
 * Revision 1.1  2007/06/13 17:22:03  platin
 *
 *   - basic support for 2D spectrum data type.
 *
 *
 */
