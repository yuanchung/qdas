/***************************************************
 * aux.c
 *
 * Useful aux functions.
 *
 * By Yuan-Chung Cheng <yuanchung@ntu.edu.tw>
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
#include <gsl/gsl_errno.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_fft_complex.h>

#include "qdas.h"
#include "aux.h"

/* more numerical functions */
inline double coth(double x)
{
  return ((1.0+exp(-2.0*x))/(1.0-exp(-2.0*x)));
}

/* x^n */
int power(int x,int n) 
{
  int i;
  int ret;

  ret=1;
  for(i=0;i<n;i++) {
    ret=ret*x;
  }
  return ret;
}

/* n! */
double factorial(int n) 
{
  int i;
  double ret;
  ret=1.0;
  for(i=1;i<=n;i++) {
    ret=ret*((double)i);
  }
  return ret;
}

/* print the local time */
int print_time_stamp()
{
  time_t timer;

  timer=time(NULL);
  printf("%s",asctime(localtime(&timer)));
  return 0;
}

/* find the number of non-zero static disorder terms in a static disorder gsl_matrix */
size_t count_static_disorder_terms(gsl_matrix *m)
{
  size_t i,j,ndim;
  size_t count;

  ndim=m->size1;

  count=0;
  for(i=0;i<ndim;i++) {
    if(gsl_matrix_get(m,i,i) > DBL_EPSILON)
      count++;
  }
  /* only count upper half off-diagonal terms */
  for(i=0;i<ndim;i++) {
    for(j=i+1;j<ndim;j++) {
      if(gsl_matrix_get(m,i,j) > DBL_EPSILON)
	count++;
    }
  }
  return count;
}

/* find the number of non-zero elements in a gsl_matrix */
size_t gsl_matrix_count_nonzero(gsl_matrix *m)
{
  size_t i,j,ndim1,ndim2;
  size_t count;

  ndim1=m->size1;
  ndim2=m->size2;

  count=0;
  for(i=0;i<ndim1;i++) {
    for(j=0;j<ndim2;j++) {
      if(gsl_matrix_get(m,i,j) > DBL_EPSILON)
	count++;
    }
  }

  return count;
}

/* More matrix handeling functions */
void gsl_matrix_print(gsl_matrix *m)
{
  size_t i,j,ndim1,ndim2;

  ndim1=m->size1;
  ndim2=m->size2;

  for(i=0;i<ndim1;i++) {
    for(j=0;j<ndim2;j++) {
      printf("% 12.6f ",gsl_matrix_get(m,i,j));
    }
    printf("\n");
  }
}

/* print a gsl_matrix in a vector form */
void gsl_matrix_vprint(gsl_matrix *m)
{
  size_t i,j,ndim1,ndim2;

  ndim1=m->size1;
  ndim2=m->size2;

  for(i=0;i<ndim1;i++) {
    for(j=0;j<ndim2;j++) {
      printf("% 10.4f ",gsl_matrix_get(m,i,j));
    }
  }
}

void gsl_matrix_lprint(gsl_matrix *m)
{
  size_t i,j,ndim1,ndim2;

  ndim1=m->size1;
  ndim2=m->size2;

  for(i=0;i<ndim1;i++) {
    for(j=0;j<ndim2;j++) {
      printf("% 18.12f ",gsl_matrix_get(m,i,j));
    }
    printf("\n");
  }
}

void gsl_matrix_complex_print(gsl_matrix_complex *m)
{
  size_t i,j,ndim1,ndim2;
  double real,imag;

  ndim1=m->size1;
  ndim2=m->size2;

  for(i=0;i<ndim1;i++) {
    for(j=0;j<ndim2;j++) {
      real=GSL_REAL(gsl_matrix_complex_get(m,i,j));
      imag=GSL_IMAG(gsl_matrix_complex_get(m,i,j));
      // The following format is suitable for loading in octave
      printf("(% 5.2f,% 5.2f) ",real,imag);
    }
    printf("\n");
  }
}

void gsl_vector_complex_print(gsl_vector_complex *v)
{
  size_t i,ndim1;
  double real,imag;

  ndim1=v->size;

  for(i=0;i<ndim1;i++) {
    real=GSL_REAL(gsl_vector_complex_get(v,i));
    imag=GSL_IMAG(gsl_vector_complex_get(v,i));
    // The following format is suitable for loading in octave
    printf("(% 9.6f,% 9.6f)\n",real,imag);
  }
}

void gsl_vector_print(gsl_vector *v)
{
  size_t i,ndim1;

  ndim1=v->size;

  for(i=0;i<ndim1;i++) {
    printf("% 9.6f\n",gsl_vector_get(v,i));
  }
}

double dot_product(double *v1, double *v2)
{
  return (v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2]);
}

double gsl_matrix_trace(gsl_matrix *m)
{
  size_t i,ndim;
  double sum;

  ndim=m->size1;

  sum=0.0;
  for(i=0;i<ndim;i++) {
    sum=sum+gsl_matrix_get(m,i,i);
  }

  return sum;
}

gsl_complex gsl_matrix_complex_trace(gsl_matrix_complex *m)
{
  size_t i,ndim;
  gsl_complex sum;

  ndim=m->size1;

  sum=gsl_complex_rect(0.0,0.0);
  for(i=0;i<ndim;i++) {
    sum=gsl_complex_add(sum,gsl_matrix_complex_get(m,i,i));
  }

  return sum;
}

/* maxmal real value in a complex matrix */
double gsl_matrix_complex_realmax(gsl_matrix_complex *m)
{
  size_t i,j,ndim;
  double max,v;

  ndim=m->size1;

  max=GSL_REAL(gsl_matrix_complex_get(m,0,0));
  for(i=0;i<ndim;i++) {
    for(j=0;j<ndim;j++) { 
      v=GSL_REAL(gsl_matrix_complex_get(m,i,j));
      if(v>max) max=v;
    }
  }

  return max;
}

double gsl_matrix_complex_imagmax(gsl_matrix_complex *m)
{
  size_t i,j,ndim;
  double max,v;

  ndim=m->size1;

  max=GSL_IMAG(gsl_matrix_complex_get(m,0,0));
  for(i=0;i<ndim;i++) {
    for(j=0;j<ndim;j++) { 
      v=GSL_IMAG(gsl_matrix_complex_get(m,i,j));
      if(v>max) max=v;
    }
  }

  return max;
}

/* maxmal abs() value in a complex matrix */
double gsl_matrix_complex_absmax(gsl_matrix_complex *m)
{
  size_t i,j,ndim;
  double max,v;

  ndim=m->size1;

  max=gsl_complex_abs(gsl_matrix_complex_get(m,0,0));
  for(i=0;i<ndim;i++) {
    for(j=0;j<ndim;j++) { 
      v=gsl_complex_abs(gsl_matrix_complex_get(m,i,j));
      if(v>max) max=v;
    }
  }

  return max;
}

/* copy a real matrix into a complex matrix a=b */
int gsl_matrix_complex_copy_real(gsl_matrix_complex *a,const gsl_matrix *b) 
{
  int i,j;
  for(i=0;i<a->size1;i++) {
    for(j=0;j<a->size2;j++) {
      gsl_matrix_complex_set(a,i,j,gsl_complex_rect(gsl_matrix_get(b,i,j),0.0));
    }
  }
  return 0;
}

/* copy the real part of a complex matrix into a real matrix, a=Real(b) */
int gsl_matrix_copy_complex_real(gsl_matrix *a,const gsl_matrix_complex *b) 
{
  int i,j;
  for(i=0;i<a->size1;i++) {
    for(j=0;j<a->size2;j++) {
      gsl_matrix_set(a,i,j,GSL_REAL(gsl_matrix_complex_get(b,i,j)));
    }
  }
  return 0;
}

/* copy the imaginary part of a complex matrix into a real matrix, a=Imag(b) */
int gsl_matrix_copy_complex_imag(gsl_matrix *a,const gsl_matrix_complex *b) 
{
  int i,j;
  for(i=0;i<a->size1;i++) {
    for(j=0;j<a->size2;j++) {
      gsl_matrix_set(a,i,j,GSL_IMAG(gsl_matrix_complex_get(b,i,j)));
    }
  }
  return 0;
}

/* copy the absolute value of a complex matrix into a real matrix, a=abs(b) */
int gsl_matrix_copy_complex_abs(gsl_matrix *a,const gsl_matrix_complex *b) 
{
  int i,j;
  for(i=0;i<a->size1;i++) {
    for(j=0;j<a->size2;j++) {
      gsl_matrix_set(a,i,j,gsl_complex_abs(gsl_matrix_complex_get(b,i,j)));
    }
  }
  return 0;
}

/* invert M and save the result in invM, also return a pointer to invM */
gsl_matrix_complex *gsl_matrix_complex_invert(gsl_matrix_complex *invM, gsl_matrix_complex *M)
{
  int s;
  gsl_permutation *perm=gsl_permutation_alloc(M->size1);

  gsl_permutation_init(perm);

  //  printf("M Matrix:\n");
  //  matrix_complex_print(M);

  /* invert M */
  gsl_linalg_complex_LU_decomp(M,perm,&s);
  gsl_linalg_complex_LU_invert(M,perm,invM);

  //  printf("Inverse M Matrix:\n");
  //  matrix_complex_print(invM);

  gsl_permutation_free(perm);

  return invM;
}

struct index_value_pair_struct 
{
  size_t index;
  double value;
};

typedef struct index_value_pair_struct index_value_pair;

int index_value_pair_compfunc(const void *x, const void *y)
{
  index_value_pair *xx,*yy;
  int t;

  xx = (index_value_pair *)x;
  yy = (index_value_pair *)y;

  if (xx->value < yy->value) t = 1;
  else
    if (xx->value == yy->value) t = 0;
    else  t = -1;
  return t;
}

/* use qsort to sort an array of N double numbers and return a list 
   of integers that contains the indices of the elements in
   decending order.
   
   Input:
         ilist: pointer to an array that will be used to return the list.
         array: 
         N:     number of elements in array

   Return:
          in ilist

 */
void sort_double_array_decend(size_t *ilist, const double *array, size_t N)
{
  index_value_pair *data;
  int i;

  data=(index_value_pair *)malloc(N*sizeof(index_value_pair));
  for(i=0;i<N;i++) {
    data[i].index=i;
    data[i].value=array[i];
  }
  qsort(data, N, sizeof(index_value_pair), index_value_pair_compfunc);  
  for(i=0;i<N;i++) {
    ilist[i]=data[i].index;
  }

  free(data);
}

/* The following function computes the half-Fourier transformation
   of Ft

          /inf
   F(w) = |   exp(Iwt)*F(t)*dw
          /0

   where F(t) is given from t=0 to t=delta_t*(N-1) in a complex 
   vector on evenly distributed discrete grid points, and delta_t 
   is the time period between points.

   This function uses GSL's FFT routine to compute approximately
   F(w) and return the result in a gsl_vector_complex array */
void hFourier_forward(gsl_vector_complex *Ft , double delta_t,
		      gsl_vector *w, gsl_vector_complex *Fw)
{
  size_t N;
  double *data;
  double sqN;
  int i;

  N = Ft->size; // number of FFT points

  if(N != Fw->size || N != w->size) {
    printf("hFourier_forward: array size mismatch, abort!\n");
    exit(EXIT_FAILURE);
  }

  sqN = sqrt((double)N);
  data=malloc(2*N*sizeof(double));

  // prepare data array for GSL FFT
  for (i = 0; i < N; i++){
    FFT_REAL(data,i) = GSL_REAL(gsl_vector_complex_get(Ft,i));
    FFT_IMAG(data,i) = GSL_IMAG(gsl_vector_complex_get(Ft,i));
  }                                                                           

  gsl_fft_complex_radix2_forward (data, 1, N);                                
                                                                                
  /* Now convert the index to the frequence domain variable */
  for(i=N/2+1;i<N;i++) {
    /* from back of the array is j=-N/2...-1, N/2-1 points */
    int idx,j;
    double Fw_r,Fw_i;
    idx=i-N/2-1;
    j=i-N; // notice j<0, so we put it in the end of the array; N-idx
    gsl_vector_set(w,N-idx-1,-2.0*M_PI/delta_t/(double)N*(double)j);
    Fw_r=FFT_REAL(data,i)/sqN;
    Fw_i=FFT_IMAG(data,i)/sqN;
    gsl_vector_complex_set(Fw,N-idx-1,gsl_complex_rect(Fw_r,Fw_i));
  }
  for(i=0;i<=N/2;i++) {
    /* first N/2+1 points is j=0..N/2 */
    int idx,j;
    double Fw_r,Fw_i;
    idx=N/2-1+i;
    j=i;
    gsl_vector_set(w,N-idx-1,-2.0*M_PI/delta_t/(double)N*(double)j);
    Fw_r=FFT_REAL(data,i)/sqN;
    Fw_i=FFT_IMAG(data,i)/sqN;
    gsl_vector_complex_set(Fw,N-idx-1,gsl_complex_rect(Fw_r,Fw_i));
  }

  free(data);
}

/* The following function computes the half-Fourier transformation
   of a functipn F(t) = Fr(t) + i*Fi(t)

          /inf
   F(w) = |   exp(Iwt)*F(t)*dw
          /0

   Note: Fr and Fi should defined using "gsl_function" format.
         
   This function uses GSL's FFT routine to compute approximately
   F(w) and return the result in a gsl_vector_complex array */
void Ft_hFourier_forward(void (* Func)(double, double *,double *),
			 double *w, gsl_vector_complex *Fw)
{
  size_t N;
  double *data;
  double t_i;
  double sqN;
  int i;

  N = FFT_N; // number of FFT points

  if(N != Fw->size) {
    printf("hFourier_forward: array size and power2 mismatch, abort!\n");
    exit(EXIT_FAILURE);
  }

  sqN = sqrt((double)N);
  data=malloc(2*N*sizeof(double));

  t_i=0.0;
  for (i = 0; i < N; i++){
    double Fr,Fi;
    Func(t_i,&Fr,&Fi);
    FFT_REAL(data,i) = Fr;
    FFT_IMAG(data,i) = Fi;
    t_i=t_i+FFT_DELTA_T;
  }                                                                           

  gsl_fft_complex_radix2_forward (data, 1, N);                                
                                                                                
  /* Now convert the index to the frequence domain variable */
  for(i=N/2+1;i<N;i++) {
    /* from back of the array is j=-N/2...-1, N/2-1 points */
    int idx,j;
    double Fw_r,Fw_i;
    idx=i-N/2-1;
    j=i-N; // notice j<0, so we put it in the end of the array; N-idx
    w[N-idx-1]=-2.0*M_PI/FFT_DELTA_T/(double)N*(double)j;
    Fw_r=FFT_REAL(data,i)/sqN;
    Fw_i=FFT_IMAG(data,i)/sqN;
    gsl_vector_complex_set(Fw,N-idx-1,gsl_complex_rect(Fw_r,Fw_i));
  }
  for(i=0;i<=N/2;i++) {
    /* first N/2+1 points is j=0..N/2 */
    int idx,j;
    double Fw_r,Fw_i;
    idx=N/2-1+i;
    j=i;
    w[N-idx-1]=-2.0*M_PI/FFT_DELTA_T/(double)N*(double)j;
    Fw_r=FFT_REAL(data,i)/sqN;
    Fw_i=FFT_IMAG(data,i)/sqN;
    gsl_vector_complex_set(Fw,N-idx-1,gsl_complex_rect(Fw_r,Fw_i));
  }

  free(data);
}

/* fill in the two-exciton part of a matrix from its one-exciton elements;
   only 2es block terms will be changed. See Notebook 1, p. 88
   all off-block diagonal terms (1es and 2es cross terms) will not ne changed.
 */
void gsl_matrix_fill_2es_from_1es(gsl_matrix *M, size_t ntes, const te_state *tes_list)
{
  int s1,s2;
  int l1,m1,n1,l2,m2,n2;
  double dtmp;

  for(s1=0;s1<ntes;s1++) {
    l1=tes_list[s1].label;
    m1=tes_list[s1].site1;
    n1=tes_list[s1].site2;
    for(s2=0;s2<ntes;s2++) {
      l2=tes_list[s2].label;
      m2=tes_list[s2].site1;
      n2=tes_list[s2].site2;
      // see p. 88 of notebook 1
      dtmp=0.0;
      if(n1 == n2) {
	dtmp=dtmp+gsl_matrix_get(M,m1,m2);
      }
      if(m1 == n2) {
	dtmp=dtmp+gsl_matrix_get(M,n1,m2);
      }
      if(n1 == m2) {
	dtmp=dtmp+gsl_matrix_get(M,m1,n2);
      }
      if(m1 == m2) {
	dtmp=dtmp+gsl_matrix_get(M,n1,n2);
      }
      gsl_matrix_set(M,l1,l2,dtmp);
    }
  }
  // DONE
}

/* construct an instance of Monte-Carlo sampled Hamiltonian, the result will be saved in H */
void construct_Hamiltonian(gsl_matrix *H, const qdas_keys *keys, gsl_rng *r)
{
  int i,j;
  double sigma_ij,dtmp;
  gsl_matrix *deltaH;

  deltaH=gsl_matrix_alloc(keys->nsize,keys->nsize);
  gsl_matrix_set_zero(deltaH);

  /* generate Gaussian disorder for each matrix elements */
  for(i=0;i<keys->nsize;i++) {
    for(j=i;j<keys->nsize;j++) {
      sigma_ij=gsl_matrix_get(keys->disorder,i,j);
      if(sigma_ij<1e-6) {
	dtmp=0.0; /* no static disorder */
      } else {
	dtmp=gsl_ran_gaussian(r,sigma_ij);
      }
      gsl_matrix_set(deltaH,i,j,dtmp);
      gsl_matrix_set(deltaH,j,i,dtmp);
    }
  }
  /* now update the two-exciton states too */
  gsl_matrix_fill_2es_from_1es(deltaH,keys->ntes,keys->tes_list);
  /* H = He + deltaH */
  gsl_matrix_memcpy(H,keys->He);
  gsl_matrix_add(H,deltaH);
  // DONE
}

/* this function diagonalizes two matrices (H0 and H1) and returns the
   overlap matrix of the eigen vectors of the two matrices in m.
   If we denote eigenvectors of H0 as a0,a1,a2,...,aN and
   eigenvectors of H1 as b0,b1,b2,...,bN (all sorted), then
   m(i,j) = <a_i|b_j>
*/
void compute_eigen_overlap(gsl_matrix *m, const gsl_matrix *H0, const gsl_matrix *H1)
{
  gsl_eigen_symmv_workspace *w;
  gsl_matrix *mtmp1=gsl_matrix_alloc(H0->size1,H0->size1);
  gsl_matrix *U0=gsl_matrix_alloc(H0->size1,H0->size1);
  gsl_matrix *U1=gsl_matrix_alloc(H0->size1,H0->size1);
  gsl_vector *lambda=gsl_vector_alloc(H0->size1);

  /* diagonalize H0 and H1 to obtain eigenstates and U */
  w=gsl_eigen_symmv_alloc(H0->size1);
  gsl_matrix_memcpy(mtmp1,H0);
  gsl_eigen_symmv(mtmp1, lambda, U0, w);
  gsl_eigen_symmv_sort(lambda,U0,GSL_EIGEN_SORT_VAL_ASC);
  gsl_matrix_memcpy(mtmp1,H1);
  gsl_eigen_symmv(mtmp1, lambda, U1, w);
  gsl_eigen_symmv_sort(lambda,U1,GSL_EIGEN_SORT_VAL_ASC);
  gsl_eigen_symmv_free(w);

  /* m=U0^\dagger * U1 */
  gsl_blas_dgemm (CblasTrans, CblasNoTrans, 1.0, U0, U1, 0.0, m);

  // DONE
  gsl_matrix_free(mtmp1);
  gsl_matrix_free(U0);
  gsl_matrix_free(U1);
  gsl_vector_free(lambda);
}

/* this function diagonalizes two matrices (H0 and H1) and returns the
   overlap matrix of the population vectors of the two matrices in m.
   If we denote eigenvectors of H0 as a0,a1,a2,...,aN and
   eigenvectors of H1 as b0,b1,b2,...,bN (all sorted), then
   m(i,j) = sum(a_i(n)^2 * b_j(n)^2,n)
*/
void compute_eigenpop_overlap(gsl_matrix *m, const gsl_matrix *H0, const gsl_matrix *H1)
{
  gsl_eigen_symmv_workspace *w;
  gsl_matrix *mtmp1=gsl_matrix_alloc(H0->size1,H0->size1);
  gsl_matrix *U0=gsl_matrix_alloc(H0->size1,H0->size1);
  gsl_matrix *U1=gsl_matrix_alloc(H0->size1,H0->size1);
  gsl_vector *lambda=gsl_vector_alloc(H0->size1);

  /* diagonalize H0 and H1 to obtain eigenstates and U */
  w=gsl_eigen_symmv_alloc(H0->size1);
  gsl_matrix_memcpy(mtmp1,H0);
  gsl_eigen_symmv(mtmp1, lambda, U0, w);
  gsl_eigen_symmv_sort(lambda,U0,GSL_EIGEN_SORT_VAL_ASC);
  gsl_matrix_memcpy(mtmp1,H1);
  gsl_eigen_symmv(mtmp1, lambda, U1, w);
  gsl_eigen_symmv_sort(lambda,U1,GSL_EIGEN_SORT_VAL_ASC);
  gsl_eigen_symmv_free(w);

  // population overlap; we first construct "population" matrices from U; P = U .* U
  gsl_matrix_mul_elements(U0,U0);
  gsl_matrix_mul_elements(U1,U1);

  /* m=U0^\dagger * U1 */
  gsl_blas_dgemm (CblasTrans, CblasNoTrans, 1.0, U0, U1, 0.0, m);

  //  printf("Pop-overlap = \n");
  //  gsl_matrix_print(m);

  // DONE
  gsl_matrix_free(mtmp1);
  gsl_matrix_free(U0);
  gsl_matrix_free(U1);
  gsl_vector_free(lambda);
}

/* this function generates a map that relates the eigen states of H0
   to those of H1; basically finding the most overlapped pairs of the 
   eigenstates of the two Hamiltonians. the returned map is a vector
   so that the i-th eigenstate in H0 corresponds to the map[i]-th eigenstate
   in H1. */
void compute_eigen_map(size_t *map, const gsl_matrix *H0, const gsl_matrix *H1)
{
  size_t n1,n2,dim;
  gsl_matrix *overlap=gsl_matrix_alloc(H0->size1,H0->size1);
  
  size_t midx;
  double max,val;

  dim=H0->size1;
  // compute the population engen vector overlap map;
  // compared to the amplitude overlap,
  // this is a better measure of the "similarity" of two quantum states.
  compute_eigenpop_overlap(overlap,H0,H1);

  // find maximally overlapped pairs
  for(n1=0;n1<dim;n1++) {
    // find maximum and index in this row
    max=fabs(gsl_matrix_get(overlap,n1,0));
    midx=0;
    for(n2=1;n2<dim;n2++) {
      val=fabs(gsl_matrix_get(overlap,n1,n2));
      if(val > max) {
	midx=n2;
	max=val;
      }
    }
    map[n1]=midx;
    // midx is taken; set the overlap of this column to zero to avoid taking it twice
    for(n2=0;n2<dim;n2++)
      gsl_matrix_set(overlap,n2,midx,0.0); 
  }

  gsl_matrix_free(overlap);
}

/* Show eigenvalue and eigenvectors for a given Hamiltonian */
void show_eigen(const gsl_matrix *H, gsl_matrix *dipole)
{
  int i;
  gsl_eigen_symmv_workspace *w;
  gsl_matrix *mtmp1=gsl_matrix_alloc(H->size1,H->size1);
  gsl_matrix *U=gsl_matrix_alloc(H->size1,H->size1);
  gsl_matrix *mu;
  gsl_vector *lambda=gsl_vector_alloc(H->size1);

  if(dipole) {
    mu=dipole;
  } else {
    mu=gsl_matrix_alloc(H->size1,3);
    gsl_matrix_set_zero(mu);
  }

  /* diagonalize H to obtain eigenstates and U */
  gsl_matrix_memcpy(mtmp1,H);
  w=gsl_eigen_symmv_alloc(H->size1);
  gsl_eigen_symmv(mtmp1, lambda, U, w);
  gsl_eigen_symmv_sort(lambda,U,GSL_EIGEN_SORT_VAL_ASC);
  gsl_eigen_symmv_free(w);

  printf("\n");
  printf("------------------\n");
  printf("| Exciton States |\n");
  printf("------------------\n");
  printf("\n");
  printf("%12s %12s\n","Energy","Dipole Strength (|mu|^2)");
  printf("\n");
  for(i=0;i<H->size1;i++) {
    double mx,my,mz;
    int j;

    mx=my=mz=0.0;
    for(j=0;j<H->size1;j++) {
      mx=mx+gsl_matrix_get(U,j,i)*gsl_matrix_get(mu,j,0);
      my=my+gsl_matrix_get(U,j,i)*gsl_matrix_get(mu,j,1);
      mz=mz+gsl_matrix_get(U,j,i)*gsl_matrix_get(mu,j,2);
    }
    printf("%12.4f %12.4f\n",gsl_vector_get(lambda,i),mx*mx+my*my+mz*mz);
  }
  printf("\n");
  printf("U = \n");
  gsl_matrix_lprint(U);
  printf("\n");

  gsl_matrix_free(mtmp1);
  gsl_matrix_free(U);
  gsl_vector_free(lambda);
  if(!dipole) {
    gsl_matrix_free(mu);
  }

}

/*
 * $Log$
 * Revision 1.13  2007/06/23 00:54:23  platin
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
 *
 *
 */
