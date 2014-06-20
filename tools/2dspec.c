/***************************************************
 * 2dspec.c
 * 
 * By Yuan-Chung Cheng <yccheng@berkeley.edu>
 *
 * Take P(tau,t) map generated from *dm3pes output 
 * and compute 2D electronic spectrum. Note: at fixed T!!
 *
 * Format of the input (read from stdin):
 * 
 * Four columns: tau, t, realP, imagP
 *
 * ordering of tau and t is as follows:
 *
 * tau1 t1 Real[P(tau1,t1)] Imag[P(tau1,t1)]
 * tau1 t2 Real[P(tau1,t2)] Imag[P(tau1,t2)]
 * .....
 * tau1 tn Real[P(tau1,tn)] Imag[P(tau1,tn)]
 *
 * tau2 t1 Real[P(tau2,t1)] Imag[P(tau2,t1)]
 * tau2 t2 Real[P(tau2,t2)] Imag[P(tau2,t2)]
 * .....
 * tau2 tn Real[P(tau2,tn)] Imag[P(tau2,tn)]
 *
 * tau3 t1 Real[P(tau3,t1)] Imag[P(tau3,t1)]
 * tau3 t2 Real[P(tau3,t2)] Imag[P(tau3,t2)]
 * .....
 *
 *
 * For example, if the dm3pes output for T=0
 * is saved in prefix_*_0000.out, where *=tau,
 * then the following shell line will do the 
 * trick (works for gnu grep/tr/cut):

   grep -A 10000 -E "^t=[- ]+0\.0000," prefix_*_0000.out | \
        sed 's/--/P=/g' | grep P= | cut -d'_' -f2- | tr '_,' ' ' | \
        tr -d 'a-zA-Z=:|,_' | tr -s ' ' ' ' | \
        cut -d' ' -f1,3,4,5 > prefix_0000.tt

 *
 *
 * In the output:
 *
 *         omega1 is the coherence frequency, 
 *         and omega2 is the rephasing frequency.
 *
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

/* global variables and constants... */
const char *program_name;
#define TOOL_2DSPEC_VERSION "0.7"

//#define DEBUG2DSPEC 1

/* values that can be controlled by options */
/* 
   Sizes of FFT matrix which we will fill by doing 
   zero padding after data points.
   If 1024 points are used in each dimension, and the 
   step size in tau/t is 5fs, this gives about 6.5 cm^-1 resolution.
   default is 1024 points...
*/
size_t FFT_Size1=1024;
size_t FFT_Size2=1024;

/* uncomment the following to use exponential padding;
   this is still not working in this implementation. */
//#define USE_EXP_PADDING 1
// number of padding exp decay steps
//#define PAD_CONST (5.0)

/* window of the output spectrum */
double w1_min,w1_max;
double w2_min,w2_max;

/* Display usage information and exit.  */
static void usage ()
{
  printf ("Usage: %s [OPTIONS] N_tau N_t FILE\n\
\n\
       N_tau and N_t defines the size of the matrix to be processed.\n\
       FILE is the data file contain time domain traces.\n\
\n\
       Available options are:\n\
\n\
       -w1range w1min w1max\n\
       -w2range w2min w2max\n\
            min..max define the frequency range for the output window.\n\
\n\
       -fftsize n1 n2\n\
            size of the matrix for FFT. This also defines the resolution\n\
            of the output spectrum.\n\
\n",
          program_name);
}

/* 2d FFT using fftw on a complex GSL matrix; the input time domain 2d spectrum t_spec
   will be copied to the frequency domain f_spec and then FFT; 
   the frequency axes will be returned in f_spec too */
void proc_2dfft(spectrum2d *f_spec,const spectrum2d *t_spec)
{
  size_t dim1,dim2,n1,n2;
  size_t n1half,n2half;
  fftw_complex *in, *out;
  fftw_plan p;

  gsl_matrix_complex *m;
  double *freq1;
  double *freq2;

  double dt,dtau;
  double t0,tau0;

  // copy data to f_spec; t_spec will be left untouched
  spectrum2d_copy(f_spec,t_spec);
  m=f_spec->data;
  freq1=f_spec->w1;
  freq2=f_spec->w2;

  dim1=m->size1;
  dim2=m->size2;
  n1half=dim1/2;
  n2half=dim2/2;

  // FFTW data
  in = (fftw_complex*) fftw_malloc(dim1*dim2*sizeof(fftw_complex));
  out = (fftw_complex*) fftw_malloc(dim1*dim2*sizeof(fftw_complex));
  // FFTW plan; notice that we use forward FFT for both tau and t here
  p = fftw_plan_dft_2d(dim1, dim2,in,out,FFTW_FORWARD, FFTW_ESTIMATE);

  // fill in data and go
  for(n1=0;n1<dim1;n1++) {
    for(n2=0;n2<dim2;n2++) {
      in[n2+dim2*n1][0]=GSL_REAL(gsl_matrix_complex_get(m,n1,n2));
      in[n2+dim2*n1][1]=GSL_IMAG(gsl_matrix_complex_get(m,n1,n2));
    }
  }

  fftw_execute(p);

  // store the result back to m; in the order that 
  // the zero frequency is in the center
  // start from the second half of the array, they are negative freqs
  // we also flip the second dimension intensionally so that is is in order
  // when multiplied by -1
  for(n1=n1half+1;n1<dim1;n1++) {
    freq1[n1-(n1half+1)]=-1.0*(double)(dim1-n1)/(double)(dim1);
    for(n2=n2half+1;n2<dim2;n2++) {
      freq2[dim2-1-(n2-(n2half+1))]=-1.0*(double)(dim2-n2)/(double)(dim2);
      gsl_matrix_complex_set(m,n1-(n1half+1),dim2-1-(n2-(n2half+1)),
			     gsl_complex_rect(out[n2+dim2*n1][0],out[n2+dim2*n1][1]));
    }
    for(n2=0;n2<=n2half;n2++) {
      freq2[dim2-1-(n2+dim2-(n2half+1))]=(double)(n2)/(double)(dim2);
      gsl_matrix_complex_set(m,n1-(n1half+1),dim2-1-(n2+dim2-(n2half+1)),
			     gsl_complex_rect(out[n2+dim2*n1][0],out[n2+dim2*n1][1]));
    }
  }
  for(n1=0;n1<=n1half;n1++) {
    freq1[n1+dim1-(n1half+1)]=(double)(n1)/(double)(dim1);
    for(n2=n2half+1;n2<dim2;n2++) {
      freq2[dim2-1-(n2-(n2half+1))]=-1.0*(double)(dim2-n2)/(double)(dim2);
      gsl_matrix_complex_set(m,n1+dim1-(n1half+1),dim2-1-(n2-(n2half+1)),
			     gsl_complex_rect(out[n2+dim2*n1][0],out[n2+dim2*n1][1]));
    }
    for(n2=0;n2<=n2half;n2++) {
      freq2[dim2-1-(n2+dim2-(n2half+1))]=(double)(n2)/(double)(dim2);
      gsl_matrix_complex_set(m,n1+dim1-(n1half+1),dim2-1-(n2+dim2-(n2half+1)),
			     gsl_complex_rect(out[n2+dim2*n1][0],out[n2+dim2*n1][1]));
    }
  }

  fftw_destroy_plan(p);
  fftw_free(in); 
  fftw_free(out);

  // FFT done, now we adjust the spectrum according to the 2D electronic spectroscopy convention.
  // convert to angular frequency
  dtau=fabs(t_spec->w1[1]-t_spec->w1[0])/TIME_CM2FS;
  for(n1=0;n1<FFT_Size1;n1++) {
    freq1[n1]=2.0*M_PI*freq1[n1]/dtau;
  }
  dt=fabs(t_spec->w2[1]-t_spec->w2[0])/TIME_CM2FS;
  for(n2=0;n2<FFT_Size2;n2++) {
    // multiply w_t by -1 to get backward transform axis.
    freq2[n2]=-2.0*M_PI*freq2[n2]/dt;
  }

  /* now apply the phase factor to adjust to time shifts */
  tau0=t_spec->w1[0]/TIME_CM2FS;
  t0=t_spec->w2[0]/TIME_CM2FS;
  printf("\n");
  printf("Applying pahse factor exp(-w_tau*tau0)*exp(w_t*t0) for zero-time adjustments:\n");
  printf("  tau0 = %f\n",tau0*TIME_CM2FS);
  printf("  t0   = %f\n",t0*TIME_CM2FS);
  for(n1=0;n1<FFT_Size1;n1++) {
    for(n2=0;n2<FFT_Size2;n2++) {
      gsl_complex f0,ff;
      double wtau,wt;
      wtau=freq1[n1];
      wt=freq2[n2];
      //      f0=gsl_complex_polar(1.0,-2.0*M_PI*wtau*tau0-2.0*M_PI*wt*t0);
      f0=gsl_complex_polar(1.0,-1.0*wtau*tau0+wt*t0);
      ff=gsl_matrix_complex_get(m,n1,n2);
      gsl_matrix_complex_set(m,n1,n2,gsl_complex_mul(f0,ff));
    }
  }
  printf("\n");

  // DONE!!
}

spectrum2d* read_input(char *fname,int dim1,int dim2)
{
  int n1,n2;
  FILE *input_file;

  spectrum2d *spec;

  /*allocate memory space for the spectrum in the time domain; and read input */
  spec=spectrum2d_alloc(1.0,(double)FFT_Size1,1.0,1.0,(double)FFT_Size2,1.0);
  gsl_matrix_complex_set_zero(spec->data); // zero padding to the size defined by FFT_Size1 x FFT_Size2
  input_file=fopen(fname, "r");
  if(input_file == NULL){
    fprintf(stderr, "error while opening \"%s\" for reading: %s \n",
            fname,strerror(errno));
    exit(errno);
  }
  printf("\n");
  printf("Reading %dx%d complex elements from %s.\n",dim1,dim2,fname);
  printf("\n");
  for(n1=0;n1<dim1;n1++) {
    for(n2=0;n2<dim2;n2++) {
      double tau,t,Pr,Pi;
      if (! feof(input_file)) {
	fscanf(input_file,"%lf %lf %lf %lf",&tau,&t,&Pr,&Pi);
	spec->w1[n1]=tau;
	spec->w2[n2]=t;
	// we store Es(tau,T,t) = i*P(tau,T,t)
	gsl_matrix_complex_set(spec->data,n1,n2,gsl_complex_rect(-1.0*Pi,Pr));
	//	gsl_matrix_complex_set(spec->data,n1,n2,gsl_complex_rect(-1.0*Pi,0.0));
	//gsl_matrix_complex_set(spec->data,n1,n2,gsl_complex_rect(0.0,Pr));
#ifdef DEBUG2DSPEC
	printf("Read:");
	printf("tau=%f t=%f Er=%f Ei=%f\n",spec->w1[n1],spec->w2[n2],
	       GSL_REAL(gsl_matrix_complex_get(spec->data,n1,n2)),
	       GSL_IMAG(gsl_matrix_complex_get(spec->data,n1,n2)));
#endif
      } else {
	fprintf(stderr, "EOF while reading \"%s\", might be mismatch matrix dimensions.\n",
		fname);
	exit(EXIT_FAILURE);
      }
      //      printf("%lf %lf %lf %lf\n",tau,t,Pr,Pi);
    }
  }
  fclose(input_file);

  /* test data integrity */
  if(spec->w1[1]<=spec->w1[0]) {
    fprintf(stderr, "Incorrect ordering of tau points; data must be given in ascending tau points.\n");
    exit(EXIT_FAILURE);
  }
  if(spec->w2[1]<=spec->w2[0]) {
    fprintf(stderr, "Incorrect ordering of t points; data must be given in ascending t points.\n");
    exit(EXIT_FAILURE);
  }

#ifndef USE_EXP_PADDING
  printf("\n");
  printf("Zero-padding data points to form %lux%lu matrix for FFT.\n",FFT_Size1,FFT_Size2);
  printf("\n");
#else
  // pad extra points using a exponential tail
  printf("\n");
  printf("Exponentially pad data points to form %dx%d matrix for FFT.\n",FFT_Size1,FFT_Size2);
  printf("\n");
  for(n1=0;n1<dim1;n1++) { // pad dim2 direction
    gsl_complex f0,ff;
    f0=gsl_matrix_complex_get(spec->data,n1,dim2-1); // f0 is the last elements
    for(n2=dim2;n2<FFT_Size2;n2++) {
      ff=gsl_complex_mul_real(f0,exp(-1.0*(double)(n2-dim2+1)/PAD_CONST));
      gsl_matrix_complex_set(spec->data,n1,n2,ff);
    }
  }
  for(n2=0;n1<FFT_Size2;n1++) { // pad dim1 direction
    gsl_complex f0,ff;
    f0=gsl_matrix_complex_get(spec->data,dim1-1,n2); // f0 is the last elements
    for(n1=dim1;n1<FFT_Size1;n1++) {
      ff=gsl_complex_mul_real(f0,exp(-1.0*(double)(n1-dim1+1)/PAD_CONST));
      gsl_matrix_complex_set(spec->data,n1,n2,ff);
    }
  }
#endif

#ifdef DEBUG2DSPEC
  printf("After padding:\n");
  for(n1=0;n1<FFT_Size1;n1++) {
    for(n2=0;n2<FFT_Size2;n2++) {
      printf("tau=%6d t=%6d PrP=%f PiP=%f\n",n1,n2,
	     GSL_REAL(gsl_matrix_complex_get(spec->data,n1,n2)),
	     GSL_IMAG(gsl_matrix_complex_get(spec->data,n1,n2)));
    }
    printf("PrP=\n");
  }
#endif

  return spec;
}


/* Main program */
int main(int argc, char *argv[])
{
  int dim1,dim2,n1,n2;
  char *fname;

  int npt;

  spectrum2d *input;
  spectrum2d *output;

  /* Set program name for messages.  */
  program_name = argv[0]; 

  /* HEADER MESSAGE */
  printf("\n");
  printf("%s\n",program_name);
  printf("\n");
  printf("        2DSPEC\n");
  printf("        Version %s.\n",TOOL_2DSPEC_VERSION);
  printf("        CVS revision: %s.\n","$Revision: 258 $");
  printf("        Two-dimensional Fast Fourier Transform program.\n");
  printf("\n");
  printf("        Part of the QDAS package, version %s.\n",QDAS_VERSION);
  printf("\n");
  printf("        Copyright(C) 2007.\n");
  printf("        Yuan-Chung Cheng <yccheng@berkeley.edu>.\n");
  printf("\n");

  /* default values for the options */
  w1_min=-1500.0;  
  w1_max=1500.0;
  w2_min=-1500.0;  
  w2_max=1500.0;

  /* handle options */
  npt=1; // this points to the "free" position in the argv[] array
  while(1) {
    if(!strncmp(argv[npt],"-w1range",6)) {
      // -w1range for the window of output
      if(argc < npt+3) {
	usage();
	exit(EXIT_FAILURE);
      } else {
	w1_min=atof(argv[npt+1])<=atof(argv[npt+2]) ? atof(argv[npt+1]) : atof(argv[npt+2]);
	w1_max=atof(argv[npt+1])<=atof(argv[npt+2]) ? atof(argv[npt+2]) : atof(argv[npt+1]);
	npt=npt+3;
      }
    } else if(!strncmp(argv[npt],"-w2range",6)) {
      // -w1range for the window of output
      if(argc < npt+3) {
	usage();
	exit(EXIT_FAILURE);
      } else {
	w2_min=atof(argv[npt+1])<=atof(argv[npt+2]) ? atof(argv[npt+1]) : atof(argv[npt+2]);
	w2_max=atof(argv[npt+1])<=atof(argv[npt+2]) ? atof(argv[npt+2]) : atof(argv[npt+1]);
	npt=npt+3;
      }      
    } else if(!strncmp(argv[npt],"-fftsize",4)) {
      // -w1range for the window of output
      if(argc < npt+3) {
	usage();
	exit(EXIT_FAILURE);
      } else {
	FFT_Size1=atof(argv[npt+1]);
	FFT_Size2=atof(argv[npt+2]);
	npt=npt+3;
      }      
    } else {
      // no more options found, break out of the loop
      break;
    }
  }

  /* infomation of the options */
  printf("\n");
  printf("FFT using %lux%lu matrix points.\n",FFT_Size1,FFT_Size2);
  printf("\n");
  printf("Will output the 2D spectrum in the following window:.\n");
  printf("  w1 = %f ... %f\n",w1_min,w1_max);
  printf("  w2 = %f ... %f\n",w2_min,w2_max);
  printf("\n");

  /* now handel the dimension and file arguments */
  if(argc != npt+3) {
    usage();
    exit(EXIT_FAILURE);
  }
  dim1=atoi(argv[npt]); // dim1 is for tau
  dim2=atoi(argv[npt+1]); // dim2 is for t
  fname=argv[npt+2];

  if(dim1<1 || dim2<1) {
    usage();
    exit(EXIT_FAILURE);
  }

  input=read_input(fname,dim1,dim2);
  
  // allocate memory space for frequency domain 2d spectrum then perform 2d FFT
  printf("\n");
  printf("Performing 2D FFT.\n");
  printf("\n");
  output=spectrum2d_alloc(1.0,(double)FFT_Size1,1.0,1.0,(double)FFT_Size2,1.0);
  proc_2dfft(output,input);

  printf("2D Map:\n");
  for(n1=0;n1<FFT_Size1;n1++) {
    double w1;
    w1=output->w1[n1];
    // we only output stuff in the defined window (default is +-1500 cm^-1)
    if(w1>=w1_min && w1<=w1_max) {
      for(n2=0;n2<FFT_Size2;n2++) {
	double w2;
	w2=output->w2[n2];
	if(w2>=w2_min && w2<=w2_max) {
	  printf("wtau=%f wt=%f Sr=%f Si=%f\n",w1,w2,
		 GSL_REAL(gsl_matrix_complex_get(output->data,n1,n2)),
		 GSL_IMAG(gsl_matrix_complex_get(output->data,n1,n2)));
	}
      }
      printf("Sr=\n");
    }
  }

  // normalize the abs of the spectral points to 1
  spectrum2d_normalize(output);
  printf("Normalized 2D Map:\n");
  for(n1=0;n1<FFT_Size1;n1++) {
    double w1;
    w1=output->w1[n1];
    // we only output stuff in the defined window (default is +-1500 cm^-1)
    if(w1>=w1_min && w1<=w1_max) {
      for(n2=0;n2<FFT_Size2;n2++) {
	double w2;
	w2=output->w2[n2];
	if(w2>=w2_min && w2<=w2_max) {
	  printf("wtau=%f wt=%f Sr_norm=%f Si_norm=%f\n",w1,w2,
		 GSL_REAL(gsl_matrix_complex_get(output->data,n1,n2)),
		 GSL_IMAG(gsl_matrix_complex_get(output->data,n1,n2)));
	}
      }
      printf("Sr_norm=\n");
    }
  }

  // DONE; clean up
  spectrum2d_free(input);
  spectrum2d_free(output);

  return 1;
}

