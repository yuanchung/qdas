/***************************************************
 * spectrum.c
 *
 * Provide functions to handle spectra.
 *
 * By Yuan-Chung Cheng <yuanchung@ntu.edu.tw>
 *
 *
 ***************************************************/

#define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_spline.h>

#include "spectrum.h"

// Defined Constants

// Global variables

// Functions

/* allocate space for a spectrum in range fstart..fend with spacing fstep */
spectrum *spectrum_alloc(double fstart, double fend, double fstep)
{
  int size;
  int i;
  double w;

  spectrum *spec;

  size=(int)(ceil(fabs(fend-fstart)/fabs(fstep))+1.0);

  spec=(spectrum *)malloc(sizeof(spectrum));
  spec->npoints=size;
  spec->w=(double *)malloc(size*sizeof(double));
  spec->intensity=(double *)malloc(size*sizeof(double));

  w=fstart;
  for(i=0;i<size;i++) {
    spec->w[i]=w;
    spec->intensity[i]=0.0;
    w=w+fstep;
  }

  spec->spline=NULL;
  spec->spline_acc=NULL;

  return spec;
}

/* duplicate a spectrum and return a pointer to the new spectrum */
spectrum *spectrum_dup(spectrum *spec)
{
  int size;
  int i;

  spectrum *ret;

  size=spec->npoints;

  ret=(spectrum *)malloc(sizeof(spectrum));
  ret->npoints=size;
  ret->w=(double *)malloc(size*sizeof(double));
  ret->intensity=(double *)malloc(size*sizeof(double));

  for(i=0;i<size;i++) {
    ret->w[i]=spec->w[i];
    ret->intensity[i]=spec->intensity[i];
  }

  ret->spline=NULL;
  ret->spline_acc=NULL;
  
  return ret;
}

void spectrum_free(spectrum *spec)
{
  spec->npoints=0;
  free(spec->w);
  free(spec->intensity);
  if(spec->spline) {
    gsl_spline_free(spec->spline);
    gsl_interp_accel_free(spec->spline_acc);
  }
  free(spec);
}

/* set the intensity of all frequencies to zero */
spectrum *spectrum_set_zero(spectrum *spec)
{
  int i;

  for(i=0;i<spec->npoints;i++) {
    spec->intensity[i]=0.0;
  }
  if(spec->spline) {
    gsl_spline_init(spec->spline,spec->w,spec->intensity,spec->npoints);
    gsl_interp_accel_free(spec->spline_acc);
    spec->spline_acc=gsl_interp_accel_alloc();
  }
  return spec;
}

/* normalize the maximum of the spectrum to 1.0 */
spectrum *spectrum_normalize(spectrum *spec)
{
  int i;
  double max;

  max=-1000.0;
  for(i=0;i<spec->npoints;i++) {
    if(spec->intensity[i] > max)
      max=spec->intensity[i];
  }
  for(i=0;i<spec->npoints;i++) {
    spec->intensity[i]=spec->intensity[i]/max;
  }

  if(spec->spline) {
    gsl_spline_init(spec->spline,spec->w,spec->intensity,spec->npoints);
    gsl_interp_accel_free(spec->spline_acc);
    spec->spline_acc=gsl_interp_accel_alloc();
  }

  return spec;
}

/* plain copy, the two spectra must have the same 
   number of points */
void spectrum_copy(spectrum *dest, spectrum *source)
{
  int i;

  if(dest->npoints != source->npoints) {
    printf("Spectra copy failure, number of points must be the same!\n");
    exit(EXIT_FAILURE);
  }
  for(i=0;i<dest->npoints;i++) {
    dest->w[i]=source->w[i];
    dest->intensity[i]=source->intensity[i];
  }

  if(dest->spline) {
    gsl_spline_init(dest->spline,dest->w,dest->intensity,dest->npoints);
    gsl_interp_accel_free(dest->spline_acc);
    dest->spline_acc=gsl_interp_accel_alloc();
  }
}

/* add the source spectrum to the dest */
void spectrum_addto(spectrum *dest, spectrum *source)
{
  int i;

  if(source->spline) {
    spectrum_spline_reset(source);
  } else {
    spectrum_spline_init(source);
  }

  for(i=0;i<dest->npoints;i++) {
    dest->intensity[i]=dest->intensity[i]+
      spectrum_spline_eval(source,dest->w[i]);
  }

  if(dest->spline) {
    gsl_spline_init(dest->spline,dest->w,dest->intensity,dest->npoints);
    gsl_interp_accel_free(dest->spline_acc);
    dest->spline_acc=gsl_interp_accel_alloc();
  }
}

/* w' = w + delta */
void spectrum_shift_w(spectrum *spec, double delta)
{
  int i;

  for(i=0;i<spec->npoints;i++) {
    spec->w[i]=spec->w[i]+delta;
  }

  if(spec->spline) {
    gsl_spline_init(spec->spline,spec->w,spec->intensity,spec->npoints);
    gsl_interp_accel_free(spec->spline_acc);
    spec->spline_acc=gsl_interp_accel_alloc();
  }

}

/* activate the spline function in the spectrum, so 
   that the intensity of any w value can be evaluated */
void spectrum_spline_init(spectrum *spec)
{
  spec->spline=gsl_spline_alloc(gsl_interp_cspline, spec->npoints);
  spec->spline_acc=gsl_interp_accel_alloc();
  gsl_spline_init(spec->spline,spec->w,spec->intensity,spec->npoints);
}

double spectrum_spline_eval(spectrum *spec, double w)
{
  return  gsl_spline_eval(spec->spline, w, spec->spline_acc);

}

void spectrum_spline_reset(spectrum *spec)
{
  gsl_spline_free(spec->spline);
  gsl_interp_accel_free(spec->spline_acc);
  spec->spline=gsl_spline_alloc(gsl_interp_cspline, spec->npoints);
  spec->spline_acc=gsl_interp_accel_alloc();
  gsl_spline_init(spec->spline,spec->w,spec->intensity,spec->npoints);
}


/*
 * $Log$
 * Revision 1.2  2007/01/30 21:38:29  platin
 *   - also copy w; it makes more sense this way.
 *
 * Revision 1.1.1.1  2006/05/24 00:42:19  platin
 *
 *   - initial import of the qdas package.
 *   - qdas stands for quantum dynamics and spectroscopy.
 *   - basic code inherited from the "lineshape" package.
 *
 *
 * Revision 1.4  2006/04/05 17:07:55  platin
 *
 *   - modular structure.
 *
 * Revision 1.3  2006/03/24 18:51:40  platin
 *
 *   - add spectrum_shift_w();
 *
 * Revision 1.2  2006/03/17 02:48:54  platin
 *
 *   - use Small's formula to weight vib. sidebands.
 *   - use Gaussian spectral function.
 *
 * Revision 1.1  2006/03/13 03:14:47  platin
 *
 *   - add spectrum module that can be used to handle spectra.
 *   - basically finished the lineshape part. need to test results.
 *
 *
 */
