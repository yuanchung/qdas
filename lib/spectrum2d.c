/***************************************************
 * spectrum2d.c
 *
 * Provide functions to handle 2D spectra.
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
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_matrix.h>

#include "spectrum2d.h"

// Defined Constants
#define EPSABS (1e-6)

// Global variables

// Functions

/* allocate space for a 2D spectrum in range 
   w1_start..w1_end with spacing w1_step in the first dimension,
   and w2_start..w2_end with spacing w2_step in the second dimension. */
spectrum2d *spectrum2d_alloc(double w1_start, double w1_end, double w1_step,
			     double w2_start, double w2_end, double w2_step)
{
  int size1,size2;
  int i;
  double w1,w2;

  spectrum2d *spec;

  size1=(int)(rint(fabs(w1_end-w1_start)/fabs(w1_step))+1.0);
  size2=(int)(rint(fabs(w2_end-w2_start)/fabs(w2_step))+1.0);

  spec=(spectrum2d *)malloc(sizeof(spectrum2d));
  spec->size1=size1;
  spec->size2=size2;
  spec->w1=(double *)malloc(size1*sizeof(double));
  spec->w2=(double *)malloc(size2*sizeof(double));
  spec->data=gsl_matrix_complex_alloc(size1,size2);

  w1=w1_start;
  for(i=0;i<size1;i++) {
    spec->w1[i]=w1;
    w1=w1+w1_step;
  }
  w2=w2_start;
  for(i=0;i<size2;i++) {
    spec->w2[i]=w2;
    w2=w2+w2_step;
  }

  gsl_matrix_complex_set_zero(spec->data);

  return spec;
}

/* duplicate a spectrum2d and return a pointer to the new spectrum2d */
spectrum2d *spectrum2d_dup(const spectrum2d *spec)
{
  int size1,size2;

  spectrum2d *ret;

  size1=spec->size1;
  size2=spec->size2;

  ret=spectrum2d_alloc(spec->w1[0],spec->w1[spec->size1-1],spec->w1[1]-spec->w1[0],
		       spec->w2[0],spec->w2[spec->size2-1],spec->w2[1]-spec->w2[0]);

  gsl_matrix_complex_memcpy(ret->data,spec->data);
  
  return ret;
}

void spectrum2d_free(spectrum2d *spec)
{
  free(spec->w1);
  free(spec->w2);
  gsl_matrix_complex_free(spec->data);
  free(spec);
}

/* set the intensity of all frequencies to zero */
spectrum2d *spectrum2d_set_zero(spectrum2d *spec)
{
  gsl_matrix_complex_set_zero(spec->data);
  return spec;
}

/* plain copy, the two spectra must have the same 
   number of points */
void spectrum2d_copy(spectrum2d *dest, const spectrum2d *source)
{
  int i;

  if((dest->size1 != source->size1) || (dest->size2 != source->size2)) {
    printf("2D spectra copy failure, number of points must be the same!\n");
    exit(EXIT_FAILURE);
  }
  for(i=0;i<dest->size1;i++) {
    dest->w1[i]=source->w1[i];
  }
  for(i=0;i<dest->size2;i++) {
    dest->w2[i]=source->w2[i];
  }

  gsl_matrix_complex_memcpy(dest->data,source->data);

}

/* normalize the spectrum so that the maximun of abs(S) is 1.0 */
void spectrum2d_normalize(spectrum2d *spec)
{
  double sabs,max;
  size_t n1,n2;

  max=0.0;
  for(n1=0;n1<spec->size1;n1++) {
    for(n2=0;n2<spec->size2;n2++) {
      sabs=gsl_complex_abs(gsl_matrix_complex_get(spec->data,n1,n2));
      if(sabs>max) max=sabs;
    }
  }
  gsl_matrix_complex_scale(spec->data,gsl_complex_rect(1.0/max,0.0));

}

void spectrum2d_scale(spectrum2d *spec, gsl_complex s)
{
  gsl_matrix_complex_scale(spec->data,s);
}

void spectrum2d_add(spectrum2d *spec1, const spectrum2d *spec2)
{
  int i;
  if(spec1->size1 != spec2->size1 ||spec1->size1 != spec2->size1) {
    fprintf(stderr, "Cannot combine two spectra with different dimensions.\n");
    exit(EXIT_FAILURE);
  }
  for(i=0;i<spec1->size1;i++) {
    if(fabs(spec1->w1[i] - spec2->w1[i]) > EPSABS) {
      fprintf(stderr, "Cannot combine two spectra with different w1 axises.\n");
      exit(EXIT_FAILURE);
    }
  }
  for(i=0;i<spec1->size2;i++) {
    if(fabs(spec1->w2[i] - spec2->w2[i]) > EPSABS) {
      fprintf(stderr, "Cannot combine two spectra with different w2 axises.\n");
      exit(EXIT_FAILURE);
    }
  }
  gsl_matrix_complex_add(spec1->data,spec2->data);
}

/* print a 2D spectrum in gnuplot-able format */
void spectrum2d_print(spectrum2d *spec)
{
  int n1,n2;

  for(n1=0;n1<spec->size1;n1++) {
    for(n2=0;n2<spec->size2;n2++) {
      printf("wtau=%f wt=%f Sr=%f Si=%f\n",spec->w1[n1],spec->w2[n2],
	     GSL_REAL(gsl_matrix_complex_get(spec->data,n1,n2)),
	     GSL_IMAG(gsl_matrix_complex_get(spec->data,n1,n2)));
    }
    printf("Sr=\n");
  }

}

/*
 * $Log$
 * Revision 1.3  2007/06/18 20:43:23  platin
 *   - add const qualitifiers.
 *
 * Revision 1.2  2007/06/18 20:41:53  platin
 *
 *   - minor changes; add the normalize function.
 *
 * Revision 1.1  2007/06/13 17:22:03  platin
 *
 *   - basic support for 2D spectrum data type.
 *
 *
 */
