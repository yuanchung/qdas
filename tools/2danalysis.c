/***************************************************
 * 2danalysis.c
 * 
 * By Yuan-Chung Cheng <yccheng@berkeley.edu>
 *
 * Analysis tool for 2d spectrum generated using 2dspec.
 *
 * Format of the input (read from stdin):
 * 
 * Four columns: omega1, omega2, realS, imagS
 *
 * Convention: omega1 is the coherence frequency, 
 *         and omega2 is the rephasing frequency.
 *
 * ordering of omega1_ and omega2_ is as follows:
 *
 * omega1_1 omega2_1 Real[S(omega1_1,omega2_1)] Imag[S(omega1_1,omega2_1)]
 * omega1_1 omega2_2 Real[S(omega1_1,omega2_2)] Imag[S(omega1_1,omega2_2)]
 * .....
 * omega1_1 omega2_n Real[S(omega1_1,omega2_n)] Imag[S(omega1_1,omega2_n)]
 *
 * omega1_2 omega2_1 Real[S(omega1_2,omega2_1)] Imag[S(omega1_2,omega2_1)]
 * omega1_2 omega2_2 Real[S(omega1_2,omega2_2)] Imag[S(omega1_2,omega2_2)]
 * .....
 * omega1_2 omega2_n Real[S(omega1_2,omega2_n)] Imag[S(omega1_2,omega2_n)]
 *
 * omega1_3 omega2_1 Real[S(omega1_3,omega2_1)] Imag[S(omega1_3,omega2_1)]
 * omega1_3 omega2_2 Real[S(omega1_3,omega2_2)] Imag[S(omega1_3,omega2_2)]
 * .....
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
#include <gsl/gsl_errno.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_interp.h>
#include <fftw3.h>

#include "qdas.h"
#include "aux.h"
#include "spectrum2d.h"

#define CM2FS (5309.1)
#define EPSABS (1e-6)

/* global variables */
const char *program_name;

/* Display usage information and exit.  */
static void usage ()
{
  printf ("Usage: %s command [command-arguments] [FILE...]\n\
\n\
       where FILE contains the 2D spectral data to be processed,\n\
       possible command and command-arguments are listed below:\n\
\n\
       read w1 w2:\n\
           show the value of spectrum at (w1,w2). Note that data interpolation\n\
           is not implemented yet; this function returns value at the nearest\n\
           point.\n\
\n\
       rmax:\n\
           show the maximum value of the real part spectrum.\n\
\n\
       imax:\n\
           show the maximum value of the imaginary part spectrum.\n\
\n\
       absmax:\n\
           show the maximum value of the absolute value spectrum.\n\
\n\
       localrmax w1_lower w1_upper w2_lower w2_upper:\n\
           show the maximum of the real part spectrum in the rectangular\n\
           region defined by [w1_lower:w1_upper] [w2_lower:w2_upper].\n\
\n\
       localrmin w1_lower w1_upper w2_lower w2_upper:\n\
           show the minimum of the real part spectrum in the rectangular\n\
           region defined by [w1_lower:w1_upper] [w2_lower:w2_upper].\n\
\n\
       dslice:\n\
           show the diagonal slice of the spectrum.\n\
\n\
       adslice w:\n\
           show the anti-diagonal slice of the spectrum, through the (w,w) point.\n\
\n\
       vslice w1:\n\
           show a vertical slice of the spectrum at omega_tau=w1.\n\
\n\
       hslice w2:\n\
           show a horizontal slice of the spectrum at omega_t=w2.\n\
\n\
       integrate w1_lower w1_upper w2_lower w2_upper:\n\
           calculate the volumn under the rectangular region\n\
           defined by [w1_lower:w1_upper] [w2_lower:w2_upper].\n\
\n\
       combine weight1 file1 weight2 file2:\n\
           add the spectra in file1 and file2 using the weights. The\n\
           two spectra have to have the same dimensionality and axises.\n\
\n\
       Note that all numbers given will be round to the closest grid point.\n\
       We also assume that the spectral file contains spectral points on a\n\
       equally-spaced 2D square of grid points.\n\
\n",
          program_name);
}

/* find the index of the element n in v so that v[n] <= val < v[n+1] */
size_t gsl_vector_find_lower(gsl_vector *v,double val)
{
  size_t n;
  for(n=0;n<v->size;n++) {
    if(gsl_vector_get(v,n) > val) break;
  }
  return (n-1);
}

/* find the index of the element n in v so that v[n-1] < val <= v[n] */
size_t gsl_vector_find_upper(gsl_vector *v,double val)
{
  size_t n;
  for(n=0;n<v->size;n++) {
    if(gsl_vector_get(v,n) >= val) break;
  }
  return (n-1);
}

/* find the index of the element n in v so that v[n] <= val < v[n+1] */
size_t double_array_find_lower(double *v,size_t size,double val)
{
  size_t n;
  for(n=0;n<size;n++) {
    if(v[n] > val) break;
  }
  return (n-1);
}

/* find the index of the element n in v so that v[n-1] < val <= v[n] */
size_t double_array_find_upper(double *v,size_t size,double val)
{
  size_t n;
  for(n=0;n<size;n++) {
    if(v[n] >= val) break;
  }
  return (n-1);
}

/* find the index of the element n in v so that v[n] is nearest to val */
size_t double_array_find_nearest(double *v,size_t size,double val)
{
  size_t n,nearest;
  double delta,min_delta;
  nearest=0;
  min_delta=fabs(val-v[0]);
  for(n=0;n<size;n++) {
    delta=fabs(v[n]-val);
    if(min_delta > delta) {
      nearest=n;
      min_delta=delta;
    }
  }

  return nearest;
}

/* read the spectrum from the input file and save the 
   result in a spectrum2d structure; note that we assume the spectrum
   is given in a square 2d grid. */
spectrum2d* read_spectrum(char *fname)
{
  FILE *input_file;
  size_t ntotal,ndim;
  size_t n1,n2;
  spectrum2d *spec;

  /* first a dry run to compute the dimension of the input spectrum */
  input_file=fopen(fname, "r");
  if(input_file == NULL){
    fprintf(stderr, "error while opening \"%s\" for reading: %s \n",
            fname,strerror(errno));
    exit(errno);
  }
  ntotal=0;
  while (! feof(input_file)) {
    double w1,w2,Sr,Si;
    fscanf(input_file,"%lf %lf %lf %lf",&w1,&w2,&Sr,&Si);
    ntotal++;
  }
  ntotal--;
  fclose(input_file);
  ndim=(size_t)rint(sqrt((double)ntotal));
  if(ntotal != ndim*ndim) {
    fprintf(stderr, "Input format error, the matrix does not seem to be square.\n");
    exit(EXIT_FAILURE);
  }

  /* now we know the dimension of the input matrix; we can allocate a 2D spectrum accordingly */
  spec=spectrum2d_alloc(1.0,(double)ndim,1.0,1.0,(double)ndim,1.0);

  /* now ready to really read the file */
  input_file=fopen(fname, "r");
  if(input_file == NULL){
    fprintf(stderr, "error while opening \"%s\" for reading: %s \n",
            fname,strerror(errno));
    exit(errno);
  }
  printf("\n");
  printf("Will read %lux%lu complex elements from %s.\n",ndim,ndim,fname);
  printf("\n");
  for(n1=0;n1<ndim;n1++) {
    for(n2=0;n2<ndim;n2++) {
      double w1,w2,Sr,Si;
      if (! feof(input_file)) {
	fscanf(input_file,"%lf %lf %lf %lf",&w1,&w2,&Sr,&Si);
	// we assume the data is equally spaced, ordered, and aligned in either w1 or w2 dim.
	spec->w1[n1]=w1;
	spec->w2[n2]=w2;
	gsl_matrix_complex_set(spec->data,n1,n2,gsl_complex_rect(Sr,Si));
      } else {
	fprintf(stderr, "EOF while reading \"%s\", might be mismatch matrix dimensions.\n",
		fname);
	exit(EXIT_FAILURE);
      }
    }
  }
  fclose(input_file);

#ifdef DEBUG2DANALYSIS
  printf("Read:\n");
  for(n1=0;n1<ndim;n1++) {
    for(n2=0;n2<ndim;n2++) {
      printf("w1=%f w2=%f Sr=%f Si=%f\n",spec->w1[n1],spec->w2[n2],
	     GSL_REAL(gsl_matrix_complex_get(spec->data,n1,n2)),
	     GSL_IMAG(gsl_matrix_complex_get(spec->data,n1,n2)));
    }
    printf("Sr=\n");
  }
#endif

  // DONE!!
  return spec;
}

/* output the maximum real part value of the spectrum. */
void process_rmax(char *fname,spectrum2d *input)
{
  double max;
  max=gsl_matrix_complex_realmax(input->data);
  printf("Maximum real part value = %f\n",max);
}

/* output the maximum imaginary part value of the spectrum. */
void process_imax(char *fname,spectrum2d *input)
{
  double max;
  max=gsl_matrix_complex_imagmax(input->data);
  printf("Maximum imaginary part value = %f\n",max);
}

/* output the maximum abs value of the spectrum. */
void process_absmax(char *fname,spectrum2d *input)
{
  double max;
  max=gsl_matrix_complex_absmax(input->data);
  printf("Maximum absolute value = %f\n",max);
}

/* output the diagonal slice of the spectrum; we assume the spectrum is given
   on a square grid with equally-spaced points */
void process_dslice(char *fname,spectrum2d *input)
{
  size_t size;
  size_t n;
  gsl_complex S;

  size=input->size1;
  for(n=0;n<size;n++) {
    S=gsl_matrix_complex_get(input->data,n,n);
    printf("Diagonal slice (%s): w_tau=%20.12f, w_t=%20.12f, Sr=%20.12f, Si=%20.12f\n",
	   fname,input->w1[n],input->w2[n],GSL_REAL(S),GSL_IMAG(S));
  }
  printf("Diagonal slice (): w_tau=, w_t=, Sr=, Si=\n");
}

/* take a anti-diagonal slice of the spectrum, about through
   the omega point */
void process_adslice(char *fname,spectrum2d *input,double omega)
{
  size_t size;
  size_t n0,n1,n2;
  gsl_complex S;
  double w0,w1,w2,dw;

  size=input->size1;

  // this returns the nearest w1 point on the spectrum
  n0=double_array_find_nearest(input->w1,size,omega); 
  w0=input->w1[n0];

  dw=fabs(input->w1[1]-input->w1[0])/10.0; // this is the crition for "no freq. difference"

  size=input->size2;
  for(n1=0;n1<size;n1++) {
    for(n2=0;n2<size;n2++) {
      w1=input->w1[n1];
      w2=input->w2[n2];
      if(fabs(w1+w2-2.0*w0) < dw) {
	/* anti-diagonal slice at w1 + w2 = 2*w0 */
	S=gsl_matrix_complex_get(input->data,n1,n2);
	printf("Antidiagonal slice (%s): w_tau=%20.12f, w_t=%20.12f, Sr=%20.12f, Si=%20.12f\n",
	       fname,input->w1[n1],input->w2[n2],GSL_REAL(S),GSL_IMAG(S));
      }
    }
  }
  printf("Antidiagonal slice (): w_tau=, w_t=, Sr=, Si=\n");
}

/* output a horizontal slice of the spectrum at fixed w2~~omega */
void process_hslice(char *fname,spectrum2d *input,double omega)
{
  size_t size;
  size_t n,nslice;
  gsl_complex S;

  size=input->size2;

  // this returns the nearest w2 point on the spectrum
  nslice=double_array_find_nearest(input->w2,size,omega); 

  size=input->size1;
  for(n=0;n<size;n++) {
    S=gsl_matrix_complex_get(input->data,n,nslice);
    printf("Horizontal slice (%s): w_tau=%20.12f, w_t=%20.12f, Sr=%20.12f, Si=%20.12f\n",
	   fname,input->w1[n],input->w2[nslice],GSL_REAL(S),GSL_IMAG(S));
  }
  printf("Horizontal slice (): w_tau=, w_t=, Sr=, Si=\n");
}

/* output a vertical slice of the spectrum at fixed w2~~omega */
void process_vslice(char *fname,spectrum2d *input,double omega)
{
  size_t size;
  size_t n,nslice;
  gsl_complex S;

  size=input->size1;

  // this returns the nearest w1 point on the spectrum
  nslice=double_array_find_nearest(input->w1,size,omega); 

  size=input->size2;
  for(n=0;n<size;n++) {
    S=gsl_matrix_complex_get(input->data,nslice,n);
    printf("Vertical slice (%s): w_tau=%20.12f, w_t=%20.12f, Sr=%20.12f, Si=%20.12f\n",
	   fname,input->w1[nslice],input->w2[n],GSL_REAL(S),GSL_IMAG(S));
  }
  printf("Vertical slice (): w_tau=, w_t=, Sr=, Si=\n");
}

/* fill a vector with Simpson 1D weight */
void fill_simpson_weight(gsl_vector *vec)
{
  size_t dim,i;
  dim=vec->size;
  if(!(dim % 2)) {
    printf("Error: even number of intervals required by Simpson rule!\n");
    exit(EXIT_FAILURE);
  }
  /* Simpson weight = (1,4,2,4,2,4,2,...,2,4,1) */
  gsl_vector_set(vec,0,1.0);
  gsl_vector_set(vec,dim-1,1.0);
  for(i=1;i<dim-1;i=i+2) {
    gsl_vector_set(vec,i,4.0);    
  }
  for(i=2;i<dim-1;i=i+2) {
    gsl_vector_set(vec,i,2.0);    
  }
}

/* sum over all elements ofa matrix using the 2D Simpson weight.
   note that this multiplied by h*k/9 equals to the volume cover 
   by the 2d matrix.

   Ref: http://mathews.ecs.fullerton.edu/n2003/SimpsonsRule2DMod.html
*/
double simpson2d_sum(gsl_matrix *m)
{
  double sum=0.0;
  size_t i,j;  
  size_t dim1,dim2;
  gsl_vector *v1,*v2;

  dim1=m->size1;
  dim2=m->size2;

  v1=gsl_vector_alloc(dim1);
  v2=gsl_vector_alloc(dim2);

  /* we first construct arrays of 1d Simpson weights */
  fill_simpson_weight(v1);
  fill_simpson_weight(v2);

#ifdef DEBUG2DANALYSIS
  printf(" v1 =\n");
  gsl_vector_print(v1);
  printf("\n");
  printf(" v2 =\n");
  gsl_vector_print(v2);
  printf("\n");
#endif

  /* sum elements with the weights */
  sum=0.0;
  for(i=0;i<dim1;i++) {
    for(j=0;j<dim2;j++) {
      sum=sum+gsl_matrix_get(m,i,j)*gsl_vector_get(v1,i)*gsl_vector_get(v2,j);
    }
  }

  // DONE
  gsl_vector_free(v1);
  gsl_vector_free(v2);

  return sum;
}

void process_integrate(char *fname,spectrum2d *spec,double w1lower, double w1upper,
		       double w2lower,double w2upper)
{
  int n1lower,n1upper,n2lower,n2upper;

  if(spec->size1<1 || spec->size2<1 || w1upper<=w1lower || w2upper<=w2lower) {
    usage();
    exit(EXIT_FAILURE);
  }

  printf("\n");
  printf("The reference integration region:\n");
  printf("  w1 = %f ... %f\n",w1lower,w1upper);
  printf("  w2 = %f ... %f\n",w2lower,w2upper);
  printf("\n");

  // now find the indices for the integration region.
  n1lower=double_array_find_lower(spec->w1,spec->size1,w1lower+fabs(w1lower*1e-6));
  n1upper=double_array_find_upper(spec->w1,spec->size1,w1upper-fabs(w1upper*1e-6));
  if((n1upper-n1lower) % 2) n1upper=n1upper+1; /* Simpson rule requires even number of intervals */
  n2lower=double_array_find_lower(spec->w2,spec->size2,w2lower+fabs(w2lower*1e-6));
  n2upper=double_array_find_upper(spec->w2,spec->size2,w2upper-fabs(w2upper*1e-6));
  if((n2upper-n2lower) % 2) n2upper=n2upper+1;
  
  printf("The actual integration region:\n");
  printf("  w1 = %f ... %f\n",spec->w1[n1lower],spec->w1[n1upper]);
  printf("  w2 = %f ... %f\n",spec->w2[n2lower],spec->w2[n2upper]);
  printf("\n");

  // perform 2D sum using the Simpson 2D rule
  {
    double dw1,dw2;
    double sum_real,sum_imag,sum_abs;

    gsl_matrix *mtmp;
    gsl_matrix_complex *data;

    int n1,n2;

    dw1=spec->w1[1]-spec->w1[0];
    dw2=spec->w2[1]-spec->w2[0];

    /* copy the integration region into a smaller matrix */
    data=gsl_matrix_complex_alloc(n1upper-n1lower+1,n2upper-n2lower+1);
    mtmp=gsl_matrix_alloc(n1upper-n1lower+1,n2upper-n2lower+1);
    for(n1=0;n1<n1upper-n1lower+1;n1++) {
      for(n2=0;n2<n2upper-n2lower+1;n2++) {
	gsl_matrix_complex_set(data,n1,n2,gsl_matrix_complex_get(spec->data,n1lower+n1,n2lower+n2));
      }
    }

    printf("\n");
    printf("Integration using 2D Simpson rule...\n");
    gsl_matrix_copy_complex_real(mtmp,data);
    sum_real=simpson2d_sum(mtmp);
    gsl_matrix_copy_complex_imag(mtmp,data);
    sum_imag=simpson2d_sum(mtmp);
    gsl_matrix_copy_complex_abs(mtmp,data);
    sum_abs=simpson2d_sum(mtmp);
    printf("Integrated volume (%s): Ir=%20.12f, Ii=%20.12f, Iabs=%20.12f \n",
	   fname,sum_real*dw1*dw2/9.0,sum_imag*dw1*dw2/9.0,sum_abs*dw1*dw2/9.0);
    printf("\n");
    gsl_matrix_free(mtmp);
    gsl_matrix_complex_free(data);
  }
}

void process_combine(char *fname1, spectrum2d *spec1, double weight1,
		     char *fname2, spectrum2d *spec2, double weight2)
{
  /* scale the spectrum by weights */
  spectrum2d_scale(spec1,gsl_complex_rect(weight1,0.0));
  spectrum2d_scale(spec2,gsl_complex_rect(weight2,0.0));

  /* then add them up in spec1 */
  spectrum2d_add(spec1,spec2);

  printf("2D Map: %f*%s + %f*%s = \n",weight1,fname1,weight2,fname2);
  spectrum2d_print(spec1);
}

/* find maximal real amplitude within a rectangular regine */
void process_local_rmax(char *fname,spectrum2d *spec,double w1lower, double w1upper,
		       double w2lower,double w2upper)
{
  int n1lower,n1upper,n2lower,n2upper;

  if(spec->size1<1 || spec->size2<1 || w1upper<=w1lower || w2upper<=w2lower) {
    usage();
    exit(EXIT_FAILURE);
  }

  printf("\n");
  printf("The reference spectral region:\n");
  printf("  w1 = %f ... %f\n",w1lower,w1upper);
  printf("  w2 = %f ... %f\n",w2lower,w2upper);
  printf("\n");

  // now find the indices for the integration region.
  n1lower=double_array_find_lower(spec->w1,spec->size1,w1lower+fabs(w1lower*1e-6));
  n1upper=double_array_find_upper(spec->w1,spec->size1,w1upper-fabs(w1upper*1e-6));
  n2lower=double_array_find_lower(spec->w2,spec->size2,w2lower+fabs(w2lower*1e-6));
  n2upper=double_array_find_upper(spec->w2,spec->size2,w2upper-fabs(w2upper*1e-6));
  
  printf("The actual spectral region:\n");
  printf("  w1 = %f ... %f\n",spec->w1[n1lower],spec->w1[n1upper]);
  printf("  w2 = %f ... %f\n",spec->w2[n2lower],spec->w2[n2upper]);
  printf("\n");

  // Search for the maximal value
  {
    size_t i,j,n1,n2;
    double max,val;
    max=GSL_REAL(gsl_matrix_complex_get(spec->data,n1lower,n2lower));
    n1=n1lower;
    n2=n2lower;
    for(i=n1lower;i<=n1upper;i++) {
      for(j=n2lower;j<=n2upper;j++) { 
	val=GSL_REAL(gsl_matrix_complex_get(spec->data,i,j));
	if(val>max) {
	  max=val;
	  n1=i;
	  n2=j;
	}
      }
    }
    printf("\n");
    printf("Maximum real part value in the region = %.0f\n",max);
    printf("Maximum real part value found at: w1 = %f, w2 = %f\n",spec->w1[n1],spec->w2[n2]);
  }

}

/* find minimal real amplitude within a rectangular regine */
void process_local_rmin(char *fname,spectrum2d *spec,double w1lower, double w1upper,
		       double w2lower,double w2upper)
{
  int n1lower,n1upper,n2lower,n2upper;

  if(spec->size1<1 || spec->size2<1 || w1upper<=w1lower || w2upper<=w2lower) {
    usage();
    exit(EXIT_FAILURE);
  }

  printf("\n");
  printf("The reference spectral region:\n");
  printf("  w1 = %f ... %f\n",w1lower,w1upper);
  printf("  w2 = %f ... %f\n",w2lower,w2upper);
  printf("\n");

  // now find the indices for the integration region.
  n1lower=double_array_find_lower(spec->w1,spec->size1,w1lower+fabs(w1lower*1e-6));
  n1upper=double_array_find_upper(spec->w1,spec->size1,w1upper-fabs(w1upper*1e-6));
  n2lower=double_array_find_lower(spec->w2,spec->size2,w2lower+fabs(w2lower*1e-6));
  n2upper=double_array_find_upper(spec->w2,spec->size2,w2upper-fabs(w2upper*1e-6));
  
  printf("The actual spectral region:\n");
  printf("  w1 = %f ... %f\n",spec->w1[n1lower],spec->w1[n1upper]);
  printf("  w2 = %f ... %f\n",spec->w2[n2lower],spec->w2[n2upper]);
  printf("\n");

  // Search for the minimal value
  {
    size_t i,j,n1,n2;
    double min,val;
    min=GSL_REAL(gsl_matrix_complex_get(spec->data,n1lower,n2lower));
    n1=n1lower;
    n2=n2lower;
    for(i=n1lower;i<=n1upper;i++) {
      for(j=n2lower;j<=n2upper;j++) { 
	val=GSL_REAL(gsl_matrix_complex_get(spec->data,i,j));
	if(val<min) {
	  min=val;
	  n1=i;
	  n2=j;
	}
      }
    }
    printf("\n");
    printf("Minimum real part value in the region = %.0f\n",min);
    printf("Minimum real part value found at: w1 = %f, w2 = %f\n",spec->w1[n1],spec->w2[n2]);
  }

}

void process_read(char *fname,spectrum2d *spec,double w1, double w2)
{
  int n1,n2;

  double rval,ival,absval;
  if(spec->size1<1 || spec->size2<1) {
    usage();
    exit(EXIT_FAILURE);
  }

  printf("\n");
  printf("Read the reference spectral point: (%f, %f)",w1,w2);
  printf("\n");

  // now find the indices for the nearest point.
  n1=double_array_find_nearest(spec->w1,spec->size1,w1);
  n2=double_array_find_nearest(spec->w2,spec->size2,w2);
  
  rval=GSL_REAL(gsl_matrix_complex_get(spec->data,n1,n2));
  ival=GSL_IMAG(gsl_matrix_complex_get(spec->data,n1,n2));
  absval=gsl_complex_abs(gsl_matrix_complex_get(spec->data,n1,n2));
  printf("Spectral value: Sr=%f, Si=%f, Sabs=%f\n",rval,ival,absval);
  printf("Spectral value read at: w1 = %f, w2 = %f\n",spec->w1[n1],spec->w2[n2]);

}

/* Main program */
int main(int argc, char *argv[])
{
  char *fname;
  char *command;

  spectrum2d *input;

  /* Set program name for messages.  */
  program_name = argv[0]; 

  /* HEADER MESSAGE */
  printf("\n");
  printf("%s\n",program_name);
  printf("\n");
  printf("        2DANALYSIS\n");
  printf("        CVS revision: %s.\n","$Revision: 258 $");
  printf("        Two-dimensional spectra ANALYSIS program.\n");
  printf("\n");
  printf("        Part of the QDAS package, version %s.\n",QDAS_VERSION);
  printf("\n");
  printf("        Copyright(C) 2007.\n");
  printf("        Yuan-Chung Cheng <yccheng@berkeley.edu>.\n");
  printf("\n");
  
  if(argc < 3) {
    usage();
    exit(EXIT_FAILURE);
  }

  /* last argument is the file name, the second one is the command. */
  fname=argv[argc-1];
  command=argv[1];

  /* we read the input file first */
  input=read_spectrum(fname);

  /* now process command and options */
  if(!strncmp(command,"rmax",4) && argc == 3) {
    process_rmax(fname,input);
  } else if(!strncmp(command,"imax",4) && argc == 3) {
    process_imax(fname,input);
  } else if(!strncmp(command,"absmax",6) && argc == 3) {
    process_absmax(fname,input);
  } else if(!strncmp(command,"dslice",6) && argc == 3) {
    process_dslice(fname,input);
  } else if(!strncmp(command,"adslice",7) && argc == 4) {
    process_adslice(fname,input,atof(argv[2]));
  } else if(!strncmp(command,"vslice",6) && argc == 4) {
    process_vslice(fname,input,atof(argv[2]));
  } else if(!strncmp(command,"hslice",6) && argc == 4) {
    process_hslice(fname,input,atof(argv[2]));
  } else if(!strncmp(command,"integrate",3) && argc == 7) {
    process_integrate(fname,input,atof(argv[2]),atof(argv[3]),atof(argv[4]),atof(argv[5]));
  } else if(!strncmp(command,"localrmax",9) && argc == 7) {
    process_local_rmax(fname,input,atof(argv[2]),atof(argv[3]),atof(argv[4]),atof(argv[5]));
  } else if(!strncmp(command,"localrmin",9) && argc == 7) {
    process_local_rmin(fname,input,atof(argv[2]),atof(argv[3]),atof(argv[4]),atof(argv[5]));
  } else if(!strncmp(command,"read",4) && argc == 5) {
    process_read(fname,input,atof(argv[2]),atof(argv[3]));
  } else if(!strncmp(command,"combine",7) && argc == 6) {
    char *fname2;
    spectrum2d *input2;
    fname2=argv[3];
    input2=read_spectrum(fname2);
    process_combine(fname,input,atof(argv[4]),fname2,input2,atof(argv[2]));
    spectrum2d_free(input2);
  } else {
    fprintf(stderr, "\n");
    fprintf(stderr, "Unknown command and/or arguments!!\n");
    fprintf(stderr, "\n");
    usage();
    exit(EXIT_FAILURE);
  }
 
  spectrum2d_free(input);

  return 1;
}
