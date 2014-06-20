/***************************************************
 * dmop.c
 * 
 * By Yuan-Chung Cheng <yccheng@mit.edu>
 *
 * Operations related to the handeling of 
 * density matrices used by qbdyn package.
 *
 ***************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_cblas.h>

#include "dmop.h"

/*
 * Use GSL gsl_matrix_complex to implement density matrices
 */

/* Allocate the space for a n*n dimensional density matrix */
gsl_matrix_complex *dmop_alloc(size_t n) {
  return gsl_matrix_complex_alloc(n,n);
}

/* free the space allocated by dmop_allocate */
void dmop_free(gsl_matrix_complex *dm) {
  gsl_matrix_complex_free(dm);
}

/* read density matrix from file "fname" , allocate the required 
   memory space, put the result into m, and return
   the dimension of the dm. */
/* Format of the .inp file:
   n
   n*n matrix --> real part of DM
   n*n matrix --> img part of the DM.
*/
size_t dmop_read(char *fname, gsl_matrix_complex **m)
{
  FILE *input_file;
  int i,j;
  char buffer[256];
  size_t dim;
  gsl_complex ztmp;

  if(fname == NULL) {
    /* resd from stdio */
    input_file=stdin;
  } else {
    input_file=fopen(fname, "r");
    if(input_file == NULL){
      fprintf(stderr, "error while opening \"%s\" for reading: %s \n",
	      fname,strerror(errno));
      exit(errno);
    }
  }

  /* the first element in the file should be the 
     dimension of the density matrix */
  buffer[0]='\0';
  fscanf(input_file,"%s",buffer);
  dim=strtoul(buffer,NULL,0);

  if(dim == 0) {
    fprintf(stderr,"error in dmop.c: error reading input dimension.\n");
    exit(EXIT_FAILURE);
  }

  /* allocate space for the density matrix */
  *m=dmop_alloc(dim);
  
  /* just plainly read in the real part; better
     representations should be considered... */
  for(i=0;i<dim;i++) {
    for(j=0;j<dim;j++) {
      if (feof(input_file)) {
	fprintf(stderr,"error in dmop.c: error reading input density matrix.\n");
	exit(EXIT_FAILURE);
      }
      fscanf(input_file,"%s",buffer);
      gsl_matrix_complex_set(*m,i,j,gsl_complex_rect(strtod(buffer,NULL),0.0));
    }
  }

  buffer[0]='\0'; /* clear buffer in case there is leftovers */
  /* read the img part until EOF or finished ... */
  for(i=0;i<dim;i++) {
    for(j=0;j<dim;j++) {
      if (! feof(input_file)) {
	fscanf(input_file,"%s",buffer);
	ztmp=gsl_matrix_complex_get(*m,i,j);
	gsl_matrix_complex_set(*m,i,j,gsl_complex_rect(GSL_REAL(ztmp),strtod(buffer,NULL)));
      }
    }
  }
  
  if(fname != NULL) {
    fflush(input_file);
    fclose(input_file);
  }

  return dim;
}

/* save the dm m into a density matrix file */
size_t dmop_save(char *fname, gsl_matrix_complex *m)
{
  FILE *output_file;
  int i,j;

  if(fname == NULL) {
    /* print to stdout */
    output_file=stdout;
  } else {
    output_file=fopen(fname, "w");
    if(output_file == NULL){
      fprintf(stderr, "error while opening \"%s\" for writing: %s \n",
	      fname,strerror(errno));
      exit(errno);
    }
  }

  fprintf(output_file, "%d\n", m->size1);
  
  /* real part of the density matrix; better
     representations should be considered... */
  for(i=0;i<m->size1;i++) {
    for(j=0;j<m->size1;j++) {
      fprintf(output_file,"%20.17f ",GSL_REAL(gsl_matrix_complex_get(m,i,j)));
    }
    fprintf(output_file,"\n");
  }
  /* then img part */
  fprintf(output_file,"\n");
  for(i=0;i<m->size1;i++) {
    for(j=0;j<m->size1;j++) {
      fprintf(output_file,"%20.17f ",GSL_IMAG(gsl_matrix_complex_get(m,i,j)));
    }
    fprintf(output_file,"\n");
  }
  if(fname != NULL)
    fclose(output_file);
  
  return m->size1;
}

/* copy the content of m2 into m1, i.e. m1=m2 */
void dmop_copy(gsl_matrix_complex *m1, gsl_matrix_complex *m2)
{
  gsl_matrix_complex_memcpy(m1,m2);
}

/* calculate the overlap between two density matrices.
   F = Tr(r1*r2) */
double dmop_overlap(gsl_matrix_complex *m1, gsl_matrix_complex *m2)
{
  size_t i,j,dim;
  gsl_complex sum;
  gsl_complex ztmp;

  dim=m1->size1;

  sum=gsl_complex_rect(0.0,0.0);

  for(i=0;i<dim;i++) {
    for(j=0;j<dim;j++) {
      ztmp=gsl_complex_mul(gsl_matrix_complex_get(m1,i,j),gsl_matrix_complex_get(m2,j,i));
      sum=gsl_complex_add(sum,ztmp);
    }
  }

  return GSL_REAL(sum);
}

/* calculate the purity of a density matrix,
   P = Tr(r*r) */
double dmop_purity(gsl_matrix_complex *m)
{
  size_t i,j,dim;
  gsl_complex sum;
  gsl_complex ztmp;

  dim=m->size1;

  sum=gsl_complex_rect(0.0,0.0);

  for(i=0;i<dim;i++) {
    for(j=0;j<dim;j++) {
      ztmp=gsl_complex_mul(gsl_matrix_complex_get(m,i,j),gsl_matrix_complex_get(m,j,i));
      sum=gsl_complex_add(sum,ztmp);
    }
  }

  return GSL_REAL(sum);
}

/* give a pure state vector described by psi, and compute
   the corresponding density matrix: rho = psi' * psi */
gsl_matrix_complex *dmop_vec2mat(gsl_matrix_complex *rho, gsl_vector_complex *psi) 
{
  int i,j;
  gsl_complex rho_ij;

  for(i=0;i<psi->size;i++) {
    for(j=0;j<psi->size;j++) {
      rho_ij=gsl_complex_mul(gsl_complex_conjugate(gsl_vector_complex_get(psi,j)),
			     gsl_vector_complex_get(psi,i));
      gsl_matrix_complex_set(rho,i,j,rho_ij);
    }
  }

  return rho;
}

void dmop_dmprint(gsl_matrix_complex *rho)
{
  int i,j;
  int dim;
  gsl_complex rho_ij;

  dim=rho->size1;

  printf("%u\n",dim);
  
  for(i=0;i<dim;i++) {
    for(j=0;j<dim;j++) {
      rho_ij=gsl_matrix_complex_get(rho,i,j);
      printf("%12.9f + %12.9fI ",GSL_REAL(rho_ij),GSL_IMAG(rho_ij));
    }
    printf("\n");
  }
}

/* print the density matrix in rho,long format */
void dmop_lfdmprint(gsl_matrix_complex *rho)
{
  size_t i,j;
  size_t ndim;
  gsl_complex rho_ij;

  ndim=rho->size1;

  for(i=0;i<ndim;i++) {
    for(j=0;j<ndim;j++) {
      rho_ij=gsl_matrix_complex_get(rho,i,j);
      printf("%18.15f + %18.15fI ",GSL_REAL(rho_ij),GSL_IMAG(rho_ij));
    }
    printf("\n");
  }
}

/* print the density matrix in m, as a Liouville vector */
void dmop_lvprint(gsl_matrix_complex *rho)
{
  size_t i,j;
  size_t ndim;
  gsl_complex rho_ij;

  ndim=rho->size1;

  for(i=0;i<ndim;i++) {
    for(j=0;j<ndim;j++) {
      rho_ij=gsl_matrix_complex_get(rho,i,j);
      printf("%18.15f + %18.15fI ",GSL_REAL(rho_ij),GSL_IMAG(rho_ij));
    }
  }
}

/* Normalize a vector */
gsl_vector_complex *dmop_vecnormalize(gsl_vector_complex *psi)
{
  int i;
  size_t dim;

  gsl_complex ctmp;
  double sum;

  dim=psi->size;

  sum=0.0;
  for(i=0;i<dim;i++) {
    sum=sum+gsl_complex_abs2(gsl_vector_complex_get(psi,i));
  }
  sum=sqrt(sum);
  for(i=0;i<dim;i++) {
    ctmp=gsl_complex_div_real(gsl_vector_complex_get(psi,i),sum);
    gsl_vector_complex_set(psi,i,ctmp);
  }

  return psi;
}

/* Normalize a density matrix */
gsl_matrix_complex *dmop_normalize(gsl_matrix_complex *rho)
{
  int i,j;
  size_t dim;

  gsl_complex ctmp;
  double sum;

  dim=rho->size1;

  sum=0.0;
  for(i=0;i<dim;i++) {
    sum=sum+GSL_REAL(gsl_matrix_complex_get(rho,i,i));
  }
  for(i=0;i<dim;i++) {
    for(j=0;j<dim;j++) {
      ctmp=gsl_complex_div_real(gsl_matrix_complex_get(rho,i,j),sum);
      gsl_matrix_complex_set(rho,i,j,ctmp);
    }
  }

  return rho;
}

gsl_matrix_complex *dmop_dmscale(gsl_matrix_complex *dm, double s) 
{
  int dim,i;
  gsl_vector_complex_view row;
  
  dim=dm->size1;

  for(i=0;i<dim;i++) {
    row=gsl_matrix_complex_row(dm,i);
    gsl_blas_zdscal(s,&row.vector);
  }

  return dm;
}

/*
 * $Log$
 * Revision 1.2  2006/07/20 17:01:30  platin
 *   - minor change on output formats.
 *
 * Revision 1.1  2006/06/08 17:49:11  platin
 *
 *   - add density matrix handeling functions (from qbdyn).
 *
 *
 */
