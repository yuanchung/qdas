/***************************************************
 * linpolar.h
 * 
 * By Yuan-Chung Cheng <yccheng@mit.edu>
 *
 * Header file for linpolar.c.
 *
 ***************************************************/

#ifndef _LINPOLAR_H
#define _LINPOLAR_H 1

#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

/* types */
/* a linear polarization vector */
typedef gsl_vector linpolar_vec;
typedef gsl_matrix linpolar_4vecs; /* four linpolar_vec s; a 3x4 matrix with 4 column vectors */
typedef gsl_matrix linpolar_rotmat; /* 3x3 rotation matrix */

/* Macros */

/* exported functions */
linpolar_vec * linpolar_vec_alloc();
void linpolar_vec_free(linpolar_vec *vec);
linpolar_4vecs * linpolar_4vecs_alloc();
void linpolar_4vecs_free(linpolar_4vecs *mat);
linpolar_rotmat * linpolar_rotmat_alloc();
void linpolar_rotmat_free(linpolar_rotmat *mat);
void linpolar_vec_set(linpolar_vec *vec, double angle);
void linpolar_4vecs_set(linpolar_4vecs *mat, double angles[]);
void linpolar_vec_print(linpolar_vec *vec);
void linpolar_4vecs_print(linpolar_4vecs *vecs);
void linpolar_4vecs_ranrot(linpolar_4vecs *mat, double align, gsl_rng *r);
void compute_pulse_interactions(gsl_matrix *mIp1, gsl_matrix *mIp2,gsl_matrix *mIp3,gsl_matrix *mIp4,
				 gsl_matrix *mX,  gsl_matrix *mY,  gsl_matrix *mZ,
				const linpolar_4vecs *polar);

/* linpolar.c */

#endif /* linpolar.h */

/*
 * $Log$
 */
