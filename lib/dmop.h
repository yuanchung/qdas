/***************************************************
 * dmop.h
 * 
 * By Yuan-Chung Cheng <yccheng@mit.edu>
 *
 * Header file for dmop.c.
 *
 ***************************************************/

#ifndef _DMOP_H
#define _DMOP_H 1

#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

/* Global variables */

/* Macros */

/* exported functions */

/* dmop.c */
gsl_matrix_complex *dmop_alloc(size_t n);
void dmop_free(gsl_matrix_complex *dm);
size_t dmop_read(char *fname, gsl_matrix_complex **m);
size_t dmop_save(char *fname, gsl_matrix_complex *m);
void dmop_product(gsl_matrix_complex *m1, gsl_matrix_complex *m2, gsl_matrix_complex *res);
void dmop_copy(gsl_matrix_complex *m1, gsl_matrix_complex *m2);
double dmop_overlap(gsl_matrix_complex *m1, gsl_matrix_complex *m2);
double dmop_purity(gsl_matrix_complex *m);
gsl_matrix_complex *dmop_vec2mat(gsl_matrix_complex *rho, gsl_vector_complex *psi);
void dmop_dmprint(gsl_matrix_complex *rho);
void dmop_lfdmprint(gsl_matrix_complex *m);
void dmop_lvprint(gsl_matrix_complex *m);
gsl_vector_complex *dmop_vecnormalize(gsl_vector_complex *psi);
gsl_matrix_complex *dmop_normalize(gsl_matrix_complex *rho);
gsl_matrix_complex *dmop_dmscale(gsl_matrix_complex *dm, double s);

#endif /* dmop.h */

/*
 * $Log$
 * Revision 1.1  2006/06/08 17:49:11  platin
 *   - add density matrix handeling functions (from qbdyn).
 *
 *
 */
