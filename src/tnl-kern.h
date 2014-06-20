/***************************************************
 * tnl-kern.h
 * 
 * By Yuan-Chung Cheng <yuanchung@ntu.edu.tw>
 *
 * Header file for tnl-kern.c, some global constants...
 *
 ***************************************************/

#ifndef _TNL_KERN_H
#define _TNL_KERN_H 1

#include "gsl/gsl_matrix.h"
#include "qdas.h"

/* Global constants */
#define TNL_KERNEL_VERSION "0.7"

/* Global variables */

/* type def */

/* stuctures used for bath and aux density matrices */
struct bath_aux_funcs_struct 
{
  gsl_matrix_complex *Sbar; /* interacting operator in the exciton basis */
  size_t Nr; /* number of terms in a(t) */
  size_t Ni; /* number of terms in b(t) */
  double Vren; /* renormalization term Vren = 2/pi * \int_0^infinity J(w)/w * dw */
  double Vstatic; /* static term = Tr{S*exp(-beta*H)}/Zs */
  
  gsl_complex Ar[256]; /* coefficients for a(t), at most 256 terms */
  gsl_complex Gr[256];
  gsl_complex Ai[256]; /* coefficients for b(t), at most 256 terms */
  gsl_complex Gi[256];

  /* aux density matrices */
  gsl_matrix_complex *dm_r[256];
  gsl_matrix_complex *dm_i[256];

  /* derivatives */
  gsl_matrix_complex *ddm_rdt[256];
  gsl_matrix_complex *ddm_idt[256];

};

typedef struct bath_aux_funcs_struct bath_aux_funcs;

struct tnl_aux_funcs_struct 
{

  gsl_matrix_complex *sigma; /* system reduced density matrix */
  gsl_matrix_complex *dsigmadt; /* and its derivatives */

  size_t nops; /* number of system-bath interacting operators; at most 256 */
  bath_aux_funcs bfs[256]; /* bath aux functions, including bath aux density matrices */

};

typedef struct tnl_aux_funcs_struct tnl_aux_funcs;

/* Function prototypes */
void setup_eigenbasis(gsl_matrix *U, gsl_vector *lambda, 
		      gsl_matrix *mX_x, gsl_matrix *mX_y, gsl_matrix *mX_z, gsl_matrix *mX_abs,
		      qdas_keys *keys, const gsl_matrix *H);
void tnl_aux_funcs_init(tnl_aux_funcs *f, size_t ndim,
			const gsl_matrix_complex *initdm, 
			const gsl_matrix_complex *Weff0, 
			const gsl_matrix_complex *W0, 
			const gsl_vector *lambda, const gsl_matrix *U,
			const qdas_keys *keys);
void tnl_aux_funcs_free(tnl_aux_funcs *f);
void advance_tnl_aux_funcs(tnl_aux_funcs *tfs, const tnl_aux_funcs *tfs_old,
			   const gsl_matrix_complex *Weff_t, const gsl_matrix_complex *W_t, 
			   const gsl_vector *lambda, const qdas_keys *keys, double delta_t);
void advance_tnl_aux_funcs_markov(tnl_aux_funcs *tfs, const tnl_aux_funcs *tfs_old,
				  const gsl_matrix_complex *Weff_t, const gsl_matrix_complex *W_t, 
				  const gsl_vector *lambda, const qdas_keys *keys, double delta_t);
#endif /* tnl-kern.h */

/*
 * $Log$
 * Revision 1.4  2007/06/23 00:54:49  platin
 *   - bump to 0.8p3.
 *   - add support for index adjustment so that poptrans terms can be
 *     consistent when level energy crossover happens at disordered jobs.
 *
 * Revision 1.3  2007/03/06 23:51:49  platin
 *
 *   - clean up BS disorder_aux_xxx code.
 *
 * Revision 1.2  2007/02/20 22:02:51  platin
 *
 *   - implemented Markov approximation in tnl-kern.c
 *   - tnl-dm3pes-2d now used the Markov model.
 *
 * Revision 1.1  2006/10/27 22:24:54  platin
 *
 *   - commit re-organized tnl-dm3pes and tnl-dynamics code.
 *
 *
 */
