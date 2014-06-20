/***************************************************
 * spectrum.h
 * 
 * Header file for spectrum.c, provides definition
 * and interfaces of spectrum functions.
 *
 * By Yuan-Chung Cheng <yuanchung@ntu.edu.tw>
 *
 *
 ***************************************************/

#ifndef _SPECTRUM_H
#define _SPECTRUM_H 1

#include <gsl/gsl_spline.h>

/* data structure to hold a spectrum */

struct spectrum_struct 
{
  size_t npoints; /* number of points in the spectrum */
  double *w;
  double *intensity;
  gsl_spline *spline;
  gsl_interp_accel *spline_acc;
};

typedef struct spectrum_struct spectrum;

/* exported functions */
spectrum *spectrum_alloc(double fstart, double fend, double fstep);
void spectrum_free(spectrum *spec);
spectrum *spectrum_dup(spectrum *spec);
spectrum *spectrum_set_zero(spectrum *spec);
spectrum *spectrum_normalize(spectrum *spec);
void spectrum_copy(spectrum *dest, spectrum *source);
void spectrum_addto(spectrum *dest, spectrum *source);
void spectrum_shift_w(spectrum *spec, double delta);
void spectrum_spline_init(spectrum *spec);
double spectrum_spline_eval(spectrum *spec, double w);
void spectrum_spline_reset(spectrum *spec);

#endif /* spectrum.h */

/*
 * $Log$
 * Revision 1.1  2006/05/24 00:42:19  platin
 * Initial revision
 *
 *
 */
