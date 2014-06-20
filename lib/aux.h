/***************************************************
 * aux.h
 * 
 * Header file for aux.c, useful aux functions.
 *
 * By Yuan-Chung Cheng <yuanchung@ntu.edu.tw>
 *
 *
 ***************************************************/

#ifndef _AUX_H
#define _AUX_H 1

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>

#include "qdas.h"

// FIXME: put these in params.c !! 
// N=8192 and delta_t=0.0001 translates into about 8 cm^-1 step in the spectrum
#define FFT_2POWER (13)
#define FFT_N (1 << FFT_2POWER)
#define FFT_DELTA_T (0.0001)

#define FFT_REAL(z,i) ((z)[2*(i)])
#define FFT_IMAG(z,i) ((z)[2*(i)+1])

/* exported functions */
double coth(double x);
int power(int x,int n);
double factorial(int n);
int print_time_stamp();
size_t count_static_disorder_terms(gsl_matrix *m);
size_t gsl_matrix_count_nonzero(gsl_matrix *m);
void gsl_matrix_print(gsl_matrix *m);
void gsl_matrix_vprint(gsl_matrix *m);
void gsl_matrix_lprint(gsl_matrix *m);
void gsl_matrix_complex_print(gsl_matrix_complex *m);
void gsl_vector_print(gsl_vector *v);
void gsl_vector_complex_print(gsl_vector_complex *v);
double dot_product(double *v1, double *v2);
double gsl_matrix_trace(gsl_matrix *m);
gsl_complex gsl_matrix_complex_trace(gsl_matrix_complex *m);
double gsl_matrix_complex_absmax(gsl_matrix_complex *m);
double gsl_matrix_complex_realmax(gsl_matrix_complex *m);
double gsl_matrix_complex_imagmax(gsl_matrix_complex *m);

int gsl_matrix_complex_copy_real(gsl_matrix_complex *a, const gsl_matrix *b);
int gsl_matrix_copy_complex_real(gsl_matrix *a,const gsl_matrix_complex *b);
int gsl_matrix_copy_complex_imag(gsl_matrix *a,const gsl_matrix_complex *b);
int gsl_matrix_copy_complex_abs(gsl_matrix *a,const gsl_matrix_complex *b);

gsl_matrix_complex *gsl_matrix_complex_invert(gsl_matrix_complex *invM, gsl_matrix_complex *M);

void hFourier_forward(gsl_vector_complex *Ft , double delta_t,
		      gsl_vector *w, gsl_vector_complex *Fw);
void Ft_hFourier_forward(void (* Func)(double, double *,double *),
			 double *w, gsl_vector_complex *Fw);

void gsl_matrix_fill_2es_from_1es(gsl_matrix *M, size_t ntes, const te_state *tes_list);
void construct_Hamiltonian(gsl_matrix *H, const qdas_keys *keys, gsl_rng *r);
void show_eigen(const gsl_matrix *H, gsl_matrix *mu);
void sort_double_array_decend(size_t *ilist, const double *array, size_t N);
void compute_eigen_map(size_t *map, const gsl_matrix *H0, const gsl_matrix *H1);

#endif /* aux.h */

/*
 * $Log$
 * Revision 1.14  2007/07/18 23:11:14  platin
 *   - remove unused code.
 *
 * Revision 1.13  2007/06/23 00:54:23  platin
 *
 *   - add support for index adjustment so that poptrans terms can be
 *     consistent when level energy crossover happens at disordered jobs.
 *
 * Revision 1.12  2007/06/21 00:35:32  platin
 *
 *   - support functions for finding maximum value of the real (or imag) part
 *     of a gsl_matrix_complex.
 *   - support wc=40 and T=77K ohmic bath.
 *
 * Revision 1.11  2007/06/20 18:14:03  platin
 *
 *   - added time stamp function.
 *
 * Revision 1.10  2007/06/01 21:36:25  platin
 *
 *   - implement sorted weights in gauss-hermite module.
 *
 * Revision 1.9  2007/05/10 01:45:28  platin
 *
 *   - add routines for converting complex matrix into real one.
 *
 * Revision 1.8  2007/03/09 08:58:56  platin
 *
 *   - more supporting gsl_* functions.
 *   - remove obsoleted code in params.c.
 *   - minor changes in gauss-hermite*
 *
 * Revision 1.7  2007/03/09 06:01:52  platin
 *
 *   - import new Gauss-Hermite Quadrature module that can be used
 *     for static disorder average.
 *
 * Revision 1.6  2006/10/27 23:02:59  platin
 *
 *   - commit modifications up to 10/27/2006.
 *
 * Revision 1.5  2006/07/28 23:32:02  platin
 *
 *   - add a aux function that returns the absmax of a complex matrix.
 *
 * Revision 1.4  2006/06/28 17:14:11  platin
 *
 *   - add the dm3pes module that computes the three-pulse photon-echo
 *     signals using Domcke's density-matrix based method.
 *   - aux functions used to implement dm3pes.
 *
 * Revision 1.3  2006/05/26 23:13:54  platin
 *
 *   - bug fix and minor changes in makefiles.
 *
 * Revision 1.2  2006/05/26 19:19:30  platin
 *
 *   - add dynamics module.
 *   - revise the params.c, use a more reasonable scheme to handle
 *     parameters for different modules.
 *
 * Revision 1.1.1.1  2006/05/24 00:42:18  platin
 *
 *   - initial import of the qdas package.
 *   - qdas stands for quantum dynamics and spectroscopy.
 *   - basic code inherited from the "lineshape" package.
 *
 */
