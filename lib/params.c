/***************************************************
 * params.c
 * 
 * By Yuan-Chung Cheng <yuanchung@ntu.edu.tw>
 *
 * Read key file and initialize the parameters
 * for the qdas package. Also provides interface to access
 * the parameters.
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
#include "gsl/gsl_matrix.h"
#include "parser.h"

#include "qdas.h"
#include "params.h"
#include "aux.h"
#include "bath.h"

/* Variables for the parameters */
qdas_keys *Keys;

/* print out the parameters read */
void print_input()
{

  int i,j;
  int n;

  n=Keys->nsize;

  printf("\n");
  printf("NSIZE = %d\n",n);
  printf("\n");
  printf("\n");
  printf("-----------------------\n");
  printf("| Hamiltonian (cm^-1) |\n");
  printf("-----------------------\n");
  printf("\n");
  gsl_matrix_print(Keys->He);
  printf("\n");
  printf("\n");
  if(Keys->nctes_elems > 0) {
    printf("----------------------------\n");
    printf("| Two-exciton Construction |\n");
    printf("----------------------------\n");
    printf("\n");
    printf("The two-exciton block is constructed from %lu one-exciton states.\n",Keys->nctes_elems);
    printf("List of one-exciton states in the construction:\n  ");
    for(i=0;i<Keys->nctes_elems;i++) {
      printf("%lu ",Keys->ctes_elems_list[i]+1);
    }
    printf("\n");
    printf("\n");
  }
  printf("----------------------\n");
  printf("| Two-exciton States |\n");
  printf("----------------------\n");
  printf("\n");
  printf("%lu two-exciton states are assigned.\n",Keys->ntes);
  printf("\n");
  if(Keys->nctes_elems == 0) {
    printf("These assignments will be applied to adjuest system-bath interactions.\n");
    printf("Note that neither Hamiltonian nor dipole moments are adjuested by \n");
    printf("these settings; to simulate absorption of these states, extra DIPOLE\n");
    printf("assignments are needed. \n");
    printf("\n");
  }
  if(Keys->ntes>0) {
    for(i=0;i<Keys->ntes;i++) {
      printf("|%2lu> = |%2lu> x |%2lu>\n",
	     Keys->tes_list[i].label+1,Keys->tes_list[i].site1+1,Keys->tes_list[i].site2+1);
    }
  }
  printf("\n");
  printf("\n");
  printf("------------------------------\n");
  printf("| External Sinus Modulations |\n");
  printf("------------------------------\n");
  printf("\n");
  printf("%lu sinus modulations will be added to the Hamiltonian.\n",Keys->nsinmoduls);
  printf("\n");
  printf("Format: |n><m|*Jm*sin(2*pi*t/Tm + Phi) and c.c. (if applicable).\n");
  printf("\n");
  printf(" %4s %4s %12s %12s %12s\n","n","m","Jm","Tm (fs)","Phi");
  for(i=0;i<Keys->nsinmoduls;i++) {
    printf(" %4lu %4lu %12.4f %12.4f %12.4f\n",
	   Keys->smodul[i].n+1,Keys->smodul[i].m+1,
	   Keys->smodul[i].Jm,Keys->smodul[i].Tm*TIME_CM2FS,Keys->smodul[i].Phi);
  }
  printf("\n");
  printf("\n");
  printf("-------------------\n");
  printf("| Static Disorder |\n");
  printf("-------------------\n");
  printf("\n");
  switch(Keys->sdmethod) {
  case QDAS_SDMETHOD_MC:
    printf("Use Monte-Carlo method for ensemble average over static disorder.\n");
    printf("NITER = %lu samples will be used in the average.\n",Keys->sdm_mc_niter);
    if(Keys->sdm_mc_rseed) {
      printf("Will seed the random number generator with %lu.\n",Keys->sdm_mc_rseed);      
    } else {
      printf("Will seed the random number generator with the system time.\n");      
    }
    break;
  case QDAS_SDMETHOD_GH:
    printf("Use Gauss-Hermite %lu point rule for ensemble average over static disorder.\n",
	   Keys->sdm_gh_order);
    printf("All weights smaller than %20.16f will be ignored.\n",Keys->sdm_gh_wmin);
    break;
  default:
    printf("Unknow method for ensemble average over static disorder!\n");
    exit(EXIT_FAILURE);
  }
  printf("\n");
  printf("Std. dev. of Gaussian disorder for each He elements:\n");
  printf("\n");
  gsl_matrix_print(Keys->disorder);
  printf("\n");
  printf("\n");
  printf("---------------------------\n");
  printf("| Bath Spectral Functions |\n");
  printf("---------------------------\n");
  printf("\n");
  printf("Beta: %f cm (%5.2f K)\n",Keys->beta,1.4387/Keys->beta);
  printf("\n");
  bath_init_params(Keys);

  printf("-------------------------\n");
  printf("| Vibrational Couplings |\n");
  printf("-------------------------\n");
  printf("\n");
  printf("%12s %9s %9s\n","Excitation","Omega","S");
  for(i=0;i<n;i++) {
    for(j=0;j<Keys->vibrations[i].n;j++) {
      printf("%12d %9.2f %9.2f\n",i+1,Keys->vibrations[i].omega[j],Keys->vibrations[i].S[j]);
    }
  }
  printf("\n");
  printf("\n");

  if(Keys->nlbath_nterms>0) {
    printf("------------------------\n");
    printf("| Non-local Bath Terms |\n");
    printf("------------------------\n");
    printf("\n");
    for(i=0;i<Keys->nlbath_nterms;i++) {
      printf("  gamma=%f, wc=%f\n",Keys->nlbath[i].gamma,Keys->nlbath[i].wc);
      printf("  S =\n");
      gsl_matrix_print(Keys->nlbath[i].S);
      printf("\n");
    }
    printf("\n");
    printf("\n");
  }

/*   if(Keys->lineshape_mod) { */
/*     printf("--------------------\n"); */
/*     printf("| Lineshape Theory |\n"); */
/*     printf("--------------------\n"); */
/*     printf("\n"); */
/*     printf("Calculate linear absorption spectrum using:\n"); */
/*     printf("\n"); */
/*     if( !strncasecmp(Keys->lineshape_mod,"JSF00",10) ) { */
/*       printf("  High-temperature Gaussian lineshape model of JSF.\n"); */
/*     } else if ( !strncasecmp(Keys->lineshape_mod,"YCC06",10) ) { */
/*       printf("  Y.C. Cheng's time-domain ILT and JS's BChl J(w).\n"); */
/*     } else if ( !strncasecmp(Keys->lineshape_mod,"YCC06S",10) ) { */
/*       printf("  Y.C. Cheng's 2nd-order time-domain ILT and JS's BChl J(w).\n"); */
/*     } else if ( !strncasecmp(Keys->lineshape_mod,"YCC06M",10) ) { */
/*       printf("  Y.C. Cheng's Makovian time-domain ILT and JS's BChl J(w).\n"); */
/*     } else { */
/*       printf("Lineshape module \"%s\" does not exist!!\n",Keys->lineshape_mod); */
/*       exit(EXIT_FAILURE); */
/*     } */
/*     printf("\n"); */
/*     printf("\n"); */
/*   } */

  printf("------------------------\n");
  printf("| Abs. Spectrum Output |\n");
  printf("------------------------\n");
  printf("\n");
  printf("When applicable, will output simulated absorption spectrum in range:\n");
  printf("Start: %9.2f cm^-1\n",Keys->spec_start);
  printf("End  : %9.2f cm^-1\n",Keys->spec_end);
  printf("Step : %9.2f cm^-1\n",Keys->spec_step);
  printf("\n");
  printf("\n");

  if(Keys->npulses>0) {
    printf("----------------\n");
    printf("| Laser Pulses |\n");
    printf("----------------\n");
    printf("\n");
    printf("Time dependent external fields represented by Gaussian pulses:\n");
    printf("\n");
    printf("E(t)=E0*exp(-4*ln(2)*(t-tau0)^2/delta^2), pulse frequency is w0\n");
    printf("\n");
    printf("Note that the labeling does not correspond to time ordering,\n");
    printf("instead, the labeling represents the k-vector numbering.\n");
    printf("\n");
    printf("%4s %14s %14s %14s %14s\n","i","E0 (cm^-1)","tau0 (fs)","delta (FWHM, fs)","w0 (cm^-1)");
    for(i=0;i<Keys->npulses;i++) {
      printf("%4d %14.6f %14.2f %14.2f %14.2f\n",
	     i+1,Keys->pulse_seq[i].E0,
	     Keys->pulse_seq[i].tau0*TIME_CM2FS,
	     Keys->pulse_seq[i].fwhm_tau*TIME_CM2FS,Keys->pulse_seq[i].w0);
    }
    printf("\n");
    printf("\n");
  }

  printf("----------------------\n");
  printf("| Laser Polarization |\n");
  printf("----------------------\n");
  printf("\n");
  if(Keys->polar) {
    printf("Pulse linear polarization angles:\n");
    printf("  pulse 1   -->  %f\n",Keys->polar->angle[0]);
    printf("  pulse 2   -->  %f\n",Keys->polar->angle[1]);
    printf("  pulse 3   -->  %f\n",Keys->polar->angle[2]);
    printf("  detection -->  %f\n",Keys->polar->angle[3]);
    printf("\n");
    printf("Use Monte-Carlo method for isotropic averaging.\n");
    printf("%lu samples will be used if no static disorder is assigned.\n",Keys->polm_mc_niter);
    if(Keys->polm_mc_rseed) {
      printf("Will seed the random number generator with %lu.\n",Keys->polm_mc_rseed);
    } else {
      printf("Will seed the random number generator with the system time.\n");      
    }
  } else {
    printf("No polarization angle is defined; the laser polarization will be ignored.\n");
  }
  printf("\n");
  printf("\n");

  if(Keys->align != 1.0) {
    printf("-----------------------------------\n");
    printf("| Anisotropic Sample |\n");
    printf("-----------------------------------\n");
    printf("\n");
    printf("Anisotropic sample with alignment factor = %f\n",Keys->align);
    printf("\n");
    printf("\n");
  }

  if(Keys->nptrans>0) {
    printf("-----------------------------------\n");
    printf("| Extra Population Transfer Terms |\n");
    printf("-----------------------------------\n");
    printf("\n");
    printf("Extra population transfer terms handled using\n");
    printf("Lindblad theorem (labeling in exciton basis):\n");
    printf("\n");
    for(i=0;i<Keys->nptrans;i++) {
      printf("  |%lu>  -->  |%lu>  with time constant of %f fs\n",
	     Keys->extra_ptrans[i].src+1,Keys->extra_ptrans[i].dest+1,
	     TIME_CM2FS/Keys->extra_ptrans[i].k);
    }
    printf("\n");
    printf("\n");
  }

  if(Keys->nsiteptrans>0) {
    printf("-------------------------------------------------------\n");
    printf("| Extra Site-representation Population Transfer Terms |\n");
    printf("-------------------------------------------------------\n");
    printf("\n");
    printf("Extra population transfer terms handled using\n");
    printf("Lindblad theorem (labeling in site basis):\n");
    printf("\n");
    for(i=0;i<Keys->nsiteptrans;i++) {
      printf("  |%lu>  -->  |%lu>  with time constant of %f fs\n",
	     Keys->extra_siteptrans[i].src+1,Keys->extra_siteptrans[i].dest+1,
	     TIME_CM2FS/Keys->extra_siteptrans[i].k);
    }
    printf("\n");
    printf("\n");
  }

  if(Keys->nnlptrans>0) {
    printf("----------------------------------------------\n");
    printf("| Extra Non-linear Population Transfer Terms |\n");
    printf("----------------------------------------------\n");
    printf("\n");
    printf("Extra non-linear pairwise transfer terms handled using\n");
    printf("Lindblad theorem (labeling in exciton basis):\n");
    printf("\n");

    for(i=0;i<Keys->nnlptrans;i++) {
      printf("  |%lu> + |%lu>  -->  |%lu>  with time constant of %f fs\n",
	     Keys->extra_nlptrans[i].src+1,Keys->extra_nlptrans[i].src+1,Keys->extra_nlptrans[i].dest+1,
	     TIME_CM2FS/Keys->extra_nlptrans[i].k);
    }
    printf("\n");
    printf("\n");
  }

  if(Keys->diagonaldynamics == QDAS_SET_DIAGONAL_DYNAMICS) {
    printf("----------------------------\n");
    printf("| Population-only dynamics |\n");
    printf("----------------------------\n");
    printf("\n");
    printf("Propagating the dynamics of the population terms only;\n");
    printf("all derivatives on the off-diagonal density matrix will be\n");
    printf("ignored!!\n");
    printf("Note that this is not a secular approximation and is only\n");
    printf("sensible when the initial density matrix does not contain\n");
    printf("non-zero off-diagonal elements!!\n");
    printf("\n");
  }

  fflush(stdout);

}


/* functions for parser */

int f_nsize_()
{
  size_t n;
  int i;

  Keys->nsize=parser_next_int();

  if(Keys->nsize < 1) {
    printf("Error: NSIZE must be >= 1.\n");
    exit(EXIT_FAILURE);
  }

  /* Now we can fill in default values for some optional parameters */
  n=Keys->nsize;
  Keys->disorder = gsl_matrix_alloc(n,n);
  gsl_matrix_set_zero(Keys->disorder);

  Keys->vibrations = malloc(n*sizeof(vibmodes));
  for(i=0;i<n;i++) {
    Keys->vibrations[i].n=0;
  }

  return 0;
}

/* aux function used to read in a matrix in MATRIX or ASSIGN 
   format; this is used to read in Hamiltonian and static disorders */
void read_matrix(gsl_matrix *m, char *description)
{
  char *str;
  int i,j;
  int n;

  n=Keys->nsize;
  if(n < 1) {
    printf("Error: NSIZE must be assigned before %s.\n",description);
    exit(EXIT_FAILURE);
  }
  
  if(!(str=parser_next_token())) {
    fprintf(stderr,"Error while reading %s!\n",description);
    exit(EXIT_FAILURE);
  }
  
  if(! strncasecmp(str,"MAT",3)) {
    /* matrix type input, read in all n*n elements */
    for(i=0;i<n;i++) {
      for(j=0;j<n;j++) {
	gsl_matrix_set(m,i,j,parser_next_double());
      }
    }
    str=parser_next_token();
    if(strncasecmp(str,"END",3)) {
      /* not end in a block ending "END" */
      fprintf(stderr,"Error while reading %s!\n",description);
      exit(EXIT_FAILURE);
    }
  } else if (! strncasecmp(str,"ASSI",4)) {
    /* assign type input, read in pair assigned elements:
       site1 site2 matrix_elements */
    
    gsl_matrix_set_zero(m);
    
    while(1) {
      str=parser_next_token();
      if(! strncasecmp(str,"END",3)) {
	break;
      } else {
	int n1,n2;
	double tmp;
	n1=atoi(str);
	n2=parser_next_int();
	tmp=parser_next_double();
	//	printf("(%lu,%lu) -> %f\n",n1,n2,tmp);
	if((n1<1 || n2<1) || (n1>n || n2>n)) {
	  fprintf(stderr,"Error reading %s, site number (%d,%d) not in 1..%d!\n",
		  description,n1,n2,n);
	  exit(EXIT_FAILURE);
	}
	/* Note: index in keyword file starts from 1, 
	   and the matrix has to be symmetric */
	gsl_matrix_set(m,n1-1,n2-1,tmp);
	gsl_matrix_set(m,n2-1,n1-1,tmp);
      }
    }
  } else {
    fprintf(stderr,"Error: unknown format type while reading %s!\n",description);
    exit(EXIT_FAILURE);
  }
}

/* read the Hamiltonian of the system */
int f_hamiltonian_()
{
  if(Keys->nsize < 1) {
    printf("Error: NSIZE must be assigned before HAMILTONIAN.\n");
    exit(EXIT_FAILURE);
  }  
  Keys->He = gsl_matrix_alloc(Keys->nsize,Keys->nsize);
  read_matrix(Keys->He,"HAMILTONIAN");
  return 0;
}

/* read the degree of disorder of the system */
int f_disorder_()
{
  size_t i,j;
  double mij,mji;

  if(Keys->nsize < 1) {
    printf("Error: NSIZE must be assigned before DISORDER.\n");
    exit(EXIT_FAILURE);
  }
  read_matrix(Keys->disorder,"DISORDER");
  for(i=0;i<Keys->nsize;i++) {
    for(j=0;j<Keys->nsize;j++) {
      mij=gsl_matrix_get(Keys->disorder,i,j);
      mji=gsl_matrix_get(Keys->disorder,j,i);
      if(mij != mji) {
	printf("Error: assignment of DISORDER is not symmetric.\n");
	exit(EXIT_FAILURE);
      }
    }
  }
  return 0;
}

int f_beta_()
{
  Keys->beta=parser_next_double();
  if(Keys->beta <= 0.0) {
    printf("Error: BETA must be positive.\n");
    exit(EXIT_FAILURE);
  }
  return 0;
}

int f_bath_()
{
  char *str;
  size_t n;
  int idx;
  bathgroup_func *bath;

  n=Keys->nsize;
  if(n < 1) {
    printf("Error: NSIZE must be assigned before BATH.\n");
    exit(EXIT_FAILURE);
  }

  // we handle multiple BATH keywords now...
  bath = Keys->bathgroup + Keys->bathgroup_nterms;

  /* bath's name following the "bath" keyword */
  str=parser_next_token();
  bath->bath_mod = bath_lookup_id(str);

  if(bath->bath_mod < 0) {
    printf("Error: cannot recognize bath type \"%s\".\n",str);
    exit(EXIT_FAILURE);
  }
  /* we then read everything up to "END" and pile them up */
  idx=0;
  while(1) {
    str=parser_next_token();
    if(str == NULL) {
      printf("Error: parse error in BATH.\n");
      exit(EXIT_FAILURE);
    }
    if(!strncasecmp(str,"END",3)) {
      break;
    }
    bath->bath_params[idx] = strtod(str,NULL);
    idx++;
  }
  bath->bath_nparams = idx;
  //printf("BATH::%lu\n",idx);
  // step up bathgroup counter
  Keys->bathgroup_nterms++;

  return 0;
}

// FIXME: nlbath is not working now
int f_nlbath_()
{
  char *str;
  size_t n;

  n=Keys->nsize;
  if(n < 1) {
    printf("Error: NSIZE must be assigned before BATH.\n");
    exit(EXIT_FAILURE);
  }
  
  /* bath's name following the "nlbath" keyword */
  str=parser_next_token();
  //  Keys->bath_mod = bath_lookup_id(str);

  //  if(Keys->bath_mod < 0) {
  //    printf("Error: cannot recognize bath type \"%s\".\n",str);
  //    exit(EXIT_FAILURE);
  //  }
  /* followed by gamma and wc */
  Keys->nlbath[Keys->nlbath_nterms].gamma=parser_next_double();
  Keys->nlbath[Keys->nlbath_nterms].wc=parser_next_double();

  /* we then read everything up to "END" into the S 
     matrix using either MAT or ASSI method */
  Keys->nlbath[Keys->nlbath_nterms].S=gsl_matrix_alloc(Keys->nsize,Keys->nsize);
  read_matrix(Keys->nlbath[Keys->nlbath_nterms].S,"NLBATH");
  Keys->nlbath_nterms++;

  return 0;
}

/* vibrational modes */
int f_vibrations_()
{
  char *str;
  int n,nvibs;

  n=Keys->nsize;
  if(n < 1) {
    printf("Error: NSIZE must be assigned before VIBRATIONS.\n");
    exit(EXIT_FAILURE);
  }

  while(1) {
    str=parser_next_token();
    if(! strncasecmp(str,"END",3)) {
      break;
    } else {
      int idx;
      double omega,S;
      vibmodes *vib;

      /* vibrations assigned as:
	 excitation_number  omega  S */
      idx=atoi(str);
      omega=parser_next_double();
      S=parser_next_double();
      if(idx<1 || idx>n) {
	fprintf(stderr,"Error reading VIBRATIONS, excitation index %d not in 1..%d!\n",
		idx,n);
	exit(EXIT_FAILURE);
      }
      if(omega<0.0) {
	fprintf(stderr,"Error in mode frequency when reading VIBRATIONS!\n");
	exit(EXIT_FAILURE);
      }
      if(S<0.0) {
	fprintf(stderr,"Error in Huang-Rhys factor when reading VIBRATIONS!\n");
	exit(EXIT_FAILURE);
      }
      /* Note: index in keyword file starts from 1 */
      vib=Keys->vibrations+(idx-1);
      vib->n=vib->n+1;
      nvibs=vib->n;
      vib->omega[nvibs-1]=omega;
      vib->S[nvibs-1]=S;
    }
  }

  return 0;
}

int f_dipole_()
{
  char *str;
  size_t nterms;
  int i,j;
  double fscale;

  if(Keys->nsize < 1) {
    printf("Error: NSIZE must be set before assigning DIPOLE.\n");
    exit(EXIT_FAILURE);
  }

  str=parser_next_token();

  if (! strncasecmp(str,"ASSI",4)) {
    /* the Assigned type, format is
       n m fscale x y z */
    nterms=0;
    str=parser_next_token();
    while(strncasecmp(str,"END",3)) {
      int n,m;
      n=atoi(str);
      m=parser_next_int();
      Keys->mu[nterms].n = ((n < m) ? n : m) - 1; // always keep n < m, also minus 1 for C convention
      Keys->mu[nterms].m = ((n < m) ? m : n) - 1;
      if(Keys->mu[nterms].n < 0 || Keys->mu[nterms].m > Keys->nsize) {
	printf("Error: assigningd state number in DIPOLE out of range.\n");
	exit(EXIT_FAILURE);	
      }
      fscale=parser_next_double();
      for(j=0;j<3;j++) {
	Keys->mu[nterms].vec[j]=fscale*parser_next_double();
      }
      nterms++;
      str=parser_next_token();
    } // while
    /* end reading, fill the number of terms */
    Keys->ndipoles=nterms;
  } else {
    /* default is matrix type input, read in all n*3 elements */
    if(! strncasecmp(str,"MAT",3)) str=parser_next_token();
    if(Keys->nsize > DIPOLES_MAX) {
      printf("Error: only support at most %d DIPOLE terms.\n",DIPOLES_MAX);
      exit(EXIT_FAILURE);
    }

    Keys->ndipoles=Keys->nsize;

    /* read in nsize amplitude and coordinates */
    for(i=0;i<Keys->ndipoles;i++) {
      fscale=atof(str); // should be the scaling factor for |\mu|
      Keys->mu[i].n=0; // assume the leading state is the groundstate
      Keys->mu[i].m=i;
      for(j=0;j<3;j++) {
	Keys->mu[i].vec[j]=fscale*parser_next_double();
      }
      str=parser_next_token(); // this one should be END or the scaling factor
    }
    if(strncasecmp(str,"END",3)) {
      /* not end in a block ending "END" */
      fprintf(stderr,"Error while reading DIPOLE!\n");
      exit(EXIT_FAILURE);
    }
  } // end treating default format

  return 0;
}

int f_sinmodul_()
{
  char *str;
  size_t nterms;

  if(Keys->nsize < 1) {
    printf("Error: NSIZE must be set before assigning SINMODUL.\n");
    exit(EXIT_FAILURE);
  }

  str=parser_next_token();

  /* sin() type modulations, format is
     n m Jm Tm */
  nterms=0;
  while(strncasecmp(str,"END",3)) {
    int n,m;
    n=atoi(str);
    m=parser_next_int();
    Keys->smodul[nterms].n = ((n < m) ? n : m) - 1; // always keep n < m, also minus 1 for C convention
    Keys->smodul[nterms].m = ((n < m) ? m : n) - 1;
    if(Keys->smodul[nterms].n < 0 || Keys->smodul[nterms].m >= Keys->nsize) {
      printf("Error: assigningd state number in SINMODUL out of range.\n");
      exit(EXIT_FAILURE);
    }
    Keys->smodul[nterms].Jm=parser_next_double();
    Keys->smodul[nterms].Tm=(parser_next_double()/TIME_CM2FS);
    Keys->smodul[nterms].Phi=parser_next_double();
    nterms++;
    str=parser_next_token();
  } // while
    /* end reading, fill the number of terms */
  Keys->nsinmoduls=nterms;

  return 0;
}

int f_poptrans_()
{
  char *str;
  size_t nterms;

  if(Keys->nsize < 1) {
    printf("Error: NSIZE must be set before assigning POPTRANS.\n");
    exit(EXIT_FAILURE);
  }

  str=parser_next_token();

  /* extra population transfer terms; these terms are 
     handled using Lindblad theorem. Format is
     src dest 1/rate */
  nterms=0;
  while(strncasecmp(str,"END",3)) {
    int a,b;
    a=atoi(str);
    b=parser_next_int();
    Keys->extra_ptrans[nterms].src = a - 1; // extra population transfer from a to b
    Keys->extra_ptrans[nterms].dest = b - 1;
    if(Keys->extra_ptrans[nterms].src < 0 || Keys->extra_ptrans[nterms].src >= Keys->nsize) {
      printf("Error: assigningd state number in POPTRANS out of range.\n");
      exit(EXIT_FAILURE);
    }
    if(Keys->extra_ptrans[nterms].dest < 0 || Keys->extra_ptrans[nterms].dest >= Keys->nsize) {
      printf("Error: assigningd state number in POPTRANS out of range.\n");
      exit(EXIT_FAILURE);
    }
    if(Keys->extra_ptrans[nterms].src == Keys->extra_ptrans[nterms].dest) {
      printf("Error: source and destination states in POPTRANS must be different.\n");
      exit(EXIT_FAILURE);
    }
    /* input is the time constant for rate in fs; we inverse it to rate in energy unit */
    Keys->extra_ptrans[nterms].k=1.0/(parser_next_double()/TIME_CM2FS);
    nterms++;
    str=parser_next_token();
  } // while
    /* end reading, fill the number of terms */
  Keys->nptrans=nterms;

  return 0;
}

int f_sitepoptrans_()
{
  char *str;
  size_t nterms;

  if(Keys->nsize < 1) {
    printf("Error: NSIZE must be set before assigning SITEPOPTRANS.\n");
    exit(EXIT_FAILURE);
  }

  str=parser_next_token();

  /* extra site-representation population transfer terms; these terms are 
     handled using Lindblad theorem. Format is
     src dest 1/rate */
  nterms=0;
  while(strncasecmp(str,"END",3)) {
    int a,b;
    a=atoi(str);
    b=parser_next_int();
    Keys->extra_siteptrans[nterms].src = a - 1; // extra population transfer from a to b
    Keys->extra_siteptrans[nterms].dest = b - 1;
    if(Keys->extra_ptrans[nterms].src < 0 || Keys->extra_siteptrans[nterms].src >= Keys->nsize) {
      printf("Error: assigningd state number in SITEPOPTRANS out of range.\n");
      exit(EXIT_FAILURE);
    }
    if(Keys->extra_siteptrans[nterms].dest < 0 || Keys->extra_siteptrans[nterms].dest >= Keys->nsize) {
      printf("Error: assigningd state number in SITEPOPTRANS out of range.\n");
      exit(EXIT_FAILURE);
    }
    if(Keys->extra_siteptrans[nterms].src == Keys->extra_siteptrans[nterms].dest) {
      printf("Error: source and destination states in SITEPOPTRANS must be different.\n");
      exit(EXIT_FAILURE);
    }
    /* input is the time constant for rate in fs; we inverse it to rate in energy unit */
    Keys->extra_siteptrans[nterms].k=1.0/(parser_next_double()/TIME_CM2FS);
    /* and allocate space for the operator matrix */
    Keys->extra_siteptrans[nterms].Op1=gsl_matrix_complex_alloc(Keys->nsize,Keys->nsize);
    Keys->extra_siteptrans[nterms].Op2=gsl_matrix_complex_alloc(Keys->nsize,Keys->nsize);
    nterms++;
    str=parser_next_token();
  } // while
    /* end reading, fill the number of terms */
  Keys->nsiteptrans=nterms;

  return 0;
}

int f_nlpoptrans_()
{
  char *str;
  size_t nterms;

  if(Keys->nsize < 1) {
    printf("Error: NSIZE must be set before assigning NLPOPTRANS.\n");
    exit(EXIT_FAILURE);
  }

  str=parser_next_token();

  /* extra nonlinear population transfer terms; these terms are 
     handled using Lindblad theorem. Format is
     src dest 1/rate (in fs) */
  nterms=0;
  while(strncasecmp(str,"END",3)) {
    int a,b;
    a=atoi(str);
    b=parser_next_int();
    Keys->extra_nlptrans[nterms].src = a - 1; // extra population transfer from a to b
    Keys->extra_nlptrans[nterms].dest = b - 1;
    if(Keys->extra_nlptrans[nterms].src < 0 || Keys->extra_nlptrans[nterms].src >= Keys->nsize) {
      printf("Error: assigningd state number in NLPOPTRANS out of range.\n");
      exit(EXIT_FAILURE);
    }
    if(Keys->extra_nlptrans[nterms].dest < 0 || Keys->extra_nlptrans[nterms].dest >= Keys->nsize) {
      printf("Error: assigningd state number in NLPOPTRANS out of range.\n");
      exit(EXIT_FAILURE);
    }
    if(Keys->extra_nlptrans[nterms].src == Keys->extra_nlptrans[nterms].dest) {
      printf("Error: source and destination states in NLPOPTRANS must be different.\n");
      exit(EXIT_FAILURE);
    }
    /* input is the time constant for rate in fs; we inverse it to rate in energy unit */
    Keys->extra_nlptrans[nterms].k=1.0/(parser_next_double()/TIME_CM2FS);
    nterms++;
    str=parser_next_token();
  } // while
    /* end reading, fill the number of terms */
  Keys->nnlptrans=nterms;

  return 0;
}

int f_teslist_()
{
  char *str;
  size_t nterms;
  int i;

  if(Keys->nsize < 1) {
    printf("Error: NSIZE must be set before assigning TESLIST.\n");
    exit(EXIT_FAILURE);
  }

  nterms=0;
  str=parser_next_token();
  while(strncasecmp(str,"END",3)) {
    // format of the list is 
    // lable n m
    // read: state lable in the Hamiltonian is a two-exciton state of state n and m.
    int n,m;
    Keys->tes_list[nterms].label = atoi(str)-1;
    n=parser_next_int();
    m=parser_next_int();
    Keys->tes_list[nterms].site1 = ((n < m) ? n : m) - 1; // always keep n < m
    Keys->tes_list[nterms].site2 = ((n < m) ? m : n) - 1;
    if(Keys->tes_list[nterms].label < 0 || Keys->tes_list[nterms].label >= Keys->nsize) {
      printf("Error: assigned state number in TESLIST out of range.\n");
      exit(EXIT_FAILURE);
    }
    if(Keys->tes_list[nterms].site1 < 0 || Keys->tes_list[nterms].site2 > Keys->nsize) {
      printf("Error: assigned site number in TESLIST out of range.\n");
      exit(EXIT_FAILURE);
    }
    if(Keys->tes_list[nterms].site1 == Keys->tes_list[nterms].site2) {
      printf("Error: incorrect two-exciton state assignment in TESLIST.\n");
      exit(EXIT_FAILURE);
    }
    // check repeated label
    for(i=0;i<nterms;i++) {
      if(Keys->tes_list[nterms].label == Keys->tes_list[i].label) {
	printf("Error: double assignment of state %lu in TESLIST.\n",Keys->tes_list[nterms].label + 1);
	exit(EXIT_FAILURE);
      }
    }
    // check repeated states
    for(i=0;i<nterms;i++) {
      if(Keys->tes_list[nterms].site1 == Keys->tes_list[i].site1 &&
	 Keys->tes_list[nterms].site2 == Keys->tes_list[i].site2 ) {
	printf("Error: double assignment of sites (%lu,%lu) in TESLIST.\n",
	       Keys->tes_list[nterms].site1 + 1,
	       Keys->tes_list[nterms].site2 + 1);
	exit(EXIT_FAILURE);
      }
    }
    nterms++;
    if(nterms>TE_STATE_MAX) {
      printf("Error: number of assignments in TESLIST exceeds the program maximum.\n");
      exit(EXIT_FAILURE);
    }
    str=parser_next_token();
  } // while
  /* end reading, fill the number of terms */
  
  Keys->ntes=nterms;

  return 0;
}

int f_cteslist_()
{

  char *str;
  size_t nterms;
  int i;

  if(Keys->nsize < 1) {
    printf("Error: NSIZE must be set before assigning TESLIST.\n");
    exit(EXIT_FAILURE);
  }

  nterms=0;
  str=parser_next_token();
  while(strncasecmp(str,"END",3)) {
    // this should just be a list of numbers indicating the
    // one-exciton states...
    Keys->ctes_elems_list[nterms] = atoi(str)-1;
    if(Keys->ctes_elems_list[nterms] < 0 || Keys->ctes_elems_list[nterms] >= Keys->nsize) {
      printf("Error: assigned site number in CTESLIST out of range.\n");
      exit(EXIT_FAILURE);
    }
    // check repeated state
    for(i=0;i<nterms;i++) {
      if(Keys->ctes_elems_list[nterms] == Keys->ctes_elems_list[i]) {
	printf("Error: double assignment of state %lu in CTESLIST.\n",Keys->ctes_elems_list[nterms] + 1);
	exit(EXIT_FAILURE);
      }
    }
    nterms++;
    if( ((nterms-1)*nterms/2) > TE_STATE_MAX) {
      printf("Error: number of assignments in CTESLIST exceeds the program capacity.\n");
      exit(EXIT_FAILURE);
    }
    str=parser_next_token();
  } // while

  if( nterms < 2) {
    printf("Error: at least two one-exciton states are required in CTESLIST.\n");
    exit(EXIT_FAILURE);
  }

  /* end reading, fill the number of terms */  
  Keys->nctes_elems=nterms;

  return 0;
}

int f_sdmethod_()
{
  char *str;

  /* name of the static-disorder module should follow the "sdmethod" keyword */
  str=parser_next_token();
  
  if(! strncasecmp(str,"MC",2)) {
    /* Monte-Carlo method, the following is a integer for the number of interations */
    Keys->sdmethod=QDAS_SDMETHOD_MC;
    Keys->sdm_mc_niter=parser_next_int();
  } else if(! strncasecmp(str,"GH",2)) {
    /* Gauss-Hermite method, the following is a integer for the 
       number of points and a double for the minimum weight accepted */
    Keys->sdmethod=QDAS_SDMETHOD_GH;
    Keys->sdm_gh_order=parser_next_int();
    Keys->sdm_gh_wmin=parser_next_double();
  } else {
    fprintf(stderr,"Error in \"SDMETHOD\"; unknown method name %s!\n",str);
    exit(EXIT_FAILURE);
  }

  return 0;
}

int f_sdrseed_()
{
  Keys->sdm_mc_rseed=(size_t)parser_next_int();
  return 0;
}

int f_polmethod_()
{
  char *str;

  /* name of the polarization sampling module should follow the "polmethod" keyword */
  str=parser_next_token();
  
  if(! strncasecmp(str,"MC",2)) {
    /* Monte-Carlo method, the following is a integer for the number of interations */
    Keys->polmethod=QDAS_POLMETHOD_MC;
    Keys->polm_mc_niter=parser_next_int();
  } else {
    fprintf(stderr,"Error in \"POLMETHOD\"; unknown method name %s!\n",str);
    exit(EXIT_FAILURE);
  }

  return 0;
}

int f_polrseed_()
{
  Keys->polm_mc_rseed=(size_t)parser_next_int();
  return 0;
}

/* int f_lineshape_() */
/* { */
/*   char *str; */

/*   /\* name of the lineshape module should follow the "lineshape" keyword *\/ */
/*   str=parser_next_token(); */
/*   Keys->lineshape_mod = strdup(str); */

/*   return 0; */
/* } */

int f_spectrum_()
{
  // spectral range: start, end, step
  Keys->spec_start=parser_next_double();
  Keys->spec_end=parser_next_double();
  Keys->spec_step=parser_next_double();

  if(Keys->spec_start >= Keys->spec_end ||
     Keys->spec_step <= 0.0) {
    fprintf(stderr,"Error in parsing \"SPECTRUM\".\n");
    exit(EXIT_FAILURE);
  }

  return 0;
}

int f_pulseseq_()
{
  char *str;
  size_t n;

  n=0;
  while(1) {
    str=parser_next_token();
    if(strncasecmp(str,"END",3)) {
      Keys->pulse_seq[n].E0=atof(str);
      Keys->pulse_seq[n].tau0=(parser_next_double()/TIME_CM2FS);
      Keys->pulse_seq[n].fwhm_tau=(parser_next_double()/TIME_CM2FS);
      Keys->pulse_seq[n].w0=parser_next_double();
      n++;
    } else {
      break;
    }
  }
  Keys->npulses=n;

  return 0;
}


int f_detect_()
{

  char *str;

  str=parser_next_token();
  while(strncasecmp(str,"END",3)) {
    /* parse keyword in this sub-block */
    if(! strncasecmp(str,"WINDOW",6)) {
      /* detection window; lower_limit followed by upper_limit */
      Keys->detect.lower=parser_next_double();
      Keys->detect.upper=parser_next_double();
      if(Keys->detect.lower > Keys->detect.upper) {
	printf("%f %f\n",Keys->detect.lower,Keys->detect.upper);
	printf("Error: lower limit of the detection window exceeds the upper limit.\n");
	exit(EXIT_FAILURE);
      }
      str=parser_next_token();
    } else {
      fprintf(stderr,"Error in \"DETECT\"; unknown keyword %s!\n",str);
      exit(EXIT_FAILURE);
    }
  }
  return 0;
}

int f_linpolar_()
{
  Keys->polar=malloc(sizeof(polar_3ppe));
  // the following 4 numbers should indicate the angles in radius
  Keys->polar->angle[0]=parser_next_double();
  Keys->polar->angle[1]=parser_next_double();
  Keys->polar->angle[2]=parser_next_double();
  Keys->polar->angle[3]=parser_next_double();

  return 0;
}

int f_align_()
{
  // alignment factor
  Keys->align=parser_next_double();
  return 0;
}

int f_tstop_()
{
  // input is in fs unit, we convert it to cm unit and store in Keys.
  Keys->tstop=(parser_next_double()/TIME_CM2FS);
  return 0;
}

int f_tstep_()
{
  // input is in fs unit, we convert it to cm unit and store in Keys.
  Keys->tstep=(parser_next_double()/TIME_CM2FS);
  return 0;
}

int f_tprint_()
{
  // input is in fs unit, we convert it to cm unit and store in Keys.
  Keys->tprint=(parser_next_double()/TIME_CM2FS);
  return 0;
}

int f_printlv_()
{
  // print detail level
  Keys->printlv=parser_next_int();
  return 0;
}

int f_cnconv_()
{
  // convegence threshold for the maximal density matrix element
  Keys->cnconv=fabs(parser_next_double());
  return 0;
}

int f_diagonaldynamics_()
{
  // turn on the populational-only dynamics
  Keys->diagonaldynamics=QDAS_SET_DIAGONAL_DYNAMICS;
  return 0;
}

/* get dipole of a 1es state */
void get_dipole(double *res, size_t state, size_t ndipoles, dipole *mu)
{
  int i;
  res[0]=0.0;res[1]=0.0;res[2]=0.0;

  for(i=0;i<ndipoles;i++) {
    // mu[] has the transition dipole from n -> m
    if(state == mu[i].m) {
      res[0]=mu[i].vec[0];
      res[1]=mu[i].vec[1];
      res[2]=mu[i].vec[2];
    }
  }
}

/* this function taks a QDAS keywork structure as an input and modify it
   to include two-exciton states; the following elements in keys will be updated:
   nsize, He, dipole, disorder, vibrations, tes_list.
   nlbath and sinmodul are broken at this point.
*/
void add_two_exciton_states(qdas_keys *keys)
{
  size_t ndim0,ndim;
  size_t nn,count;
  size_t i,j,n,m;

  gsl_matrix *old;
  vibmodes *old_vibs;

  nn=keys->nctes_elems;

  ndim0=keys->nsize; // the original dimension of the system
  ndim=ndim0 + (nn*(nn-1)/2); // the new dimension of the system

  /* we first construct the tes_list */
  count=ndim0;
  keys->ntes=0;
  for(i=0;i<nn;i++) {
    n=keys->ctes_elems_list[i];
    for(j=i+1;j<nn;j++) {
      m=keys->ctes_elems_list[j];
      keys->tes_list[keys->ntes].label=count;
      keys->tes_list[keys->ntes].site1= n < m ? n : m; // we make sure site1 < site2
      keys->tes_list[keys->ntes].site2= n > m ? n : m; // we make sure site1 < site2
      keys->ntes++;
      count++;
    }
  }

  /* update the Hamiltonian */
  old=keys->He;
  keys->He=gsl_matrix_alloc(ndim,ndim);
  gsl_matrix_set_zero(keys->He);
  // keep 1ex states
  for(i=0;i<ndim0;i++) {
    for(j=0;j<ndim0;j++) {
      gsl_matrix_set(keys->He,i,j,
		     gsl_matrix_get(old,i,j));
    }
  }
  gsl_matrix_free(old);
  // expand to 2es states
  gsl_matrix_fill_2es_from_1es(keys->He,keys->ntes, keys->tes_list);

  /* update the static disorder matrix */
  old=keys->disorder;
  keys->disorder=gsl_matrix_alloc(ndim,ndim);
  gsl_matrix_set_zero(keys->disorder);
  for(i=0;i<ndim0;i++) {
    for(j=0;j<ndim0;j++) {
      gsl_matrix_set(keys->disorder,i,j,
		     gsl_matrix_get(old,i,j));
    }
  }
  gsl_matrix_free(old);

  /* update the vibration assignments */
  old_vibs=keys->vibrations;
  keys->vibrations=malloc(ndim*sizeof(vibmodes));
  for(i=0;i<ndim0;i++) {
    keys->vibrations[i].n=old_vibs[i].n;
    for(j=0;j<keys->vibrations[i].n;j++) {
      keys->vibrations[i].omega[j] = old_vibs[i].omega[j];
      keys->vibrations[i].S[j] = old_vibs[i].S[j];
    }
  }
  for(i=ndim0;i<ndim;i++) {
    // FIXME: we ignore vibrational couplings to the two-exciton states
    keys->vibrations[i].n=0;
  }
  free(old_vibs);

  /* update the dipoles; i.e. adding 1es -> 2es transitions; 
     this is only possible in the site basis */
  for(i=0;i<keys->ntes;i++) {
    // two extra dipole terms for each tes
    // |k> -> |kl>, the transition dipole equals to that of |g> -> |l>
    keys->mu[keys->ndipoles].n = keys->tes_list[i].site1;
    keys->mu[keys->ndipoles].m = keys->tes_list[i].label;
    get_dipole(keys->mu[keys->ndipoles].vec,keys->tes_list[i].site2,keys->ndipoles,keys->mu);
    keys->ndipoles++;
    // |l> -> |kl>, the transition dipole equals to that of |g> -> |k>
    keys->mu[keys->ndipoles].n = keys->tes_list[i].site2;
    keys->mu[keys->ndipoles].m = keys->tes_list[i].label;
    get_dipole(keys->mu[keys->ndipoles].vec,keys->tes_list[i].site1,keys->ndipoles,keys->mu);
    keys->ndipoles++;
  }

  // finally update the nsize
  keys->nsize=ndim;

  // DONE
}

/* initialize the parameters using the file "filename" */
void params_init(char *fname, qdas_keys *key)
{
  FILE *input_file;
  parser_klist *list;
  int i,j;

  /* assign the global Keys */
  Keys = key;

  /* initialize bath modules before we read the keyword file; 
     the bath parameters will be read and prepared in params.c */
  bath_mod_init();

  /* initialize optional parameters */
  Keys->nsize=0;
  Keys->He = NULL;
  Keys->bathgroup_nterms = 0;
  for(i=0;i<BATHGROUP_MAX;i++) {
    Keys->bathgroup[i].bath_mod=-1;
    Keys->bathgroup[i].bath_nparams=0;
    for(j=0;j<BUFFER_SIZE;j++) {
      Keys->bathgroup[i].bath_params[j] = -100.0; // negative values mean un-initialized
    }
  }
  Keys->nlbath_nterms=0;

  Keys->sdmethod=QDAS_SDMETHOD_MC;
  Keys->sdm_mc_niter=1;
  Keys->sdm_mc_rseed=0;
  Keys->sdm_gh_order=11; // default is 11 point G-H Rule
  Keys->sdm_gh_wmin=1e-9; // minimum weight

  //  Keys->lineshape_mod = NULL;
  Keys->vibrations = NULL;
  Keys->ndipoles = 0;

  Keys->ntes = 0;
  Keys->nctes_elems=0;

  Keys->spec_start=-1500.0;
  Keys->spec_end=1500.0;
  Keys->spec_step=5.0;

  Keys->nsinmoduls=0;

  Keys->nptrans=0;
  Keys->nsiteptrans=0;
  Keys->nnlptrans=0;

  Keys->npulses=0;

  /* huge default value for the detection window, also 
     ignore output field polarization by default */
  Keys->detect.upper=1e6;
  Keys->detect.lower=-1e6;

  Keys->polar=NULL;
  Keys->polmethod=QDAS_POLMETHOD_MC;
  Keys->polm_mc_niter=1;
  Keys->polm_mc_rseed=0;

  // Alignment factor for anisotropic samples (NG)
  Keys->align=1.0;

  // Dynamical variables and convergence threshold for CN iteration
  Keys->tprint=1.0;
  Keys->tstep=0.1;
  Keys->tstop=100.0;
  Keys->cnconv=1e-12; // double machine precision is ~2E-16

  Keys->diagonaldynamics=0; // turn off population-only dynamics

  // about output detail level; default is 1, output limited stuff
  Keys->printlv=1;

  /* the parameter file */
  printf("\n");
  printf("Parameter file: %s\n",fname);
  printf("\n");
  input_file=fopen(fname, "r");

  if(input_file == NULL){
    fprintf(stderr, "error while opening \"%s\" for reading: %s \n",
	    fname,strerror(errno));
    exit(errno);
  }
  
  /* required keywords */
  list=parser_klist_alloc(300);
  parser_keyword_add(list,"NSIZE_",PARSER_REQUIRED,f_nsize_);
  parser_keyword_add(list,"HAMILTONIAN_",PARSER_REQUIRED,f_hamiltonian_);
  parser_keyword_add(list,"BETA_",PARSER_REQUIRED,f_beta_);
  parser_keyword_add(list,"BATH_",PARSER_MULTIPLE,f_bath_);

  // FIXME: nlbath is not working now
  //  parser_keyword_add(list,"NLBATH_",PARSER_MULTIPLE,f_nlbath_);
  parser_keyword_add(list,"DISORDER_",PARSER_OPTIONAL,f_disorder_);
  parser_keyword_add(list,"SDMETHOD_",PARSER_OPTIONAL,f_sdmethod_);
  parser_keyword_add(list,"SDRSEED_",PARSER_OPTIONAL,f_sdrseed_);

  /* parameters that defines the system */
  //  parser_keyword_add(list,"LINESHAPE_",PARSER_OPTIONAL,f_lineshape_);
  parser_keyword_add(list,"VIBRATIONS_",PARSER_OPTIONAL,f_vibrations_);
  parser_keyword_add(list,"DIPOLE_",PARSER_OPTIONAL,f_dipole_);
  parser_keyword_add(list,"TESLIST_",PARSER_OPTIONAL,f_teslist_);
  parser_keyword_add(list,"CTESLIST_",PARSER_OPTIONAL,f_cteslist_);

  parser_keyword_add(list,"SPECTRUM_",PARSER_OPTIONAL,f_spectrum_);

  /* parameters for the dynamics module */
  parser_keyword_add(list,"POPTRANS_",PARSER_OPTIONAL,f_poptrans_);
  parser_keyword_add(list,"SITEPOPTRANS_",PARSER_OPTIONAL,f_sitepoptrans_);
  parser_keyword_add(list,"NLPOPTRANS_",PARSER_OPTIONAL,f_nlpoptrans_);
  parser_keyword_add(list,"SINMODUL_",PARSER_OPTIONAL,f_sinmodul_);
  parser_keyword_add(list,"PULSESEQ_",PARSER_OPTIONAL,f_pulseseq_);

  parser_keyword_add(list,"LINPOLAR_",PARSER_OPTIONAL,f_linpolar_);
  parser_keyword_add(list,"POLMETHOD_",PARSER_OPTIONAL,f_polmethod_);
  parser_keyword_add(list,"POLRSEED_",PARSER_OPTIONAL,f_polrseed_);

  parser_keyword_add(list,"ALIGN_",PARSER_OPTIONAL,f_align_);

  parser_keyword_add(list,"DETECT_",PARSER_OPTIONAL,f_detect_);
  parser_keyword_add(list,"TSTOP_",PARSER_OPTIONAL,f_tstop_);
  parser_keyword_add(list,"TSTEP_",PARSER_OPTIONAL,f_tstep_);
  parser_keyword_add(list,"TPRINT_",PARSER_OPTIONAL,f_tprint_);
  parser_keyword_add(list,"PRINTLV_",PARSER_OPTIONAL,f_printlv_);
  parser_keyword_add(list,"CNCONV_",PARSER_OPTIONAL,f_cnconv_);
  parser_keyword_add(list,"DIAGONALDYNAMICS_",PARSER_OPTIONAL,f_diagonaldynamics_);

  if(!parser_parse_input(input_file,list)) {
    printf("Input file parsing error!\n");
    exit(EXIT_FAILURE);
  }

  /* done reading and parsing */
  fclose(input_file);
  parser_klist_free(list);

  /* we now do some post-processing */
  /* We need at least one BATH keywork */
  if(Keys->bathgroup_nterms<1) {
      printf("BATH assignment not found!\n");
      exit(EXIT_FAILURE);
      printf("\n");
  }

  /* handle "CTESLIST" */
  if(Keys->nctes_elems > 0) {
    if(Keys->ntes > 0) {
      printf("\n");
      printf("TESLIST and CTESLIST are mutually exclusive keywords!\n");
      printf("Please correct your keyword file.\n");
      exit(EXIT_FAILURE);
      printf("\n");
    }
    add_two_exciton_states(Keys);
  }

  print_input(fname);

}   

void params_close()
{
  int i;
  gsl_matrix_free(Keys->He);
  gsl_matrix_free(Keys->disorder);
  for(i=0;i<Keys->nsiteptrans;i++) {
    gsl_matrix_complex_free(Keys->extra_siteptrans[i].Op1);
    gsl_matrix_complex_free(Keys->extra_siteptrans[i].Op2);
  }
  if(Keys->polar)
    free(Keys->polar);

  // also close the bath module
  for(i=0;i<Keys->bathgroup_nterms;i++)
    bath_free_params(Keys->bathgroup[i].bath_mod);
}

/*
 * $Log$
 * Revision 1.19  2007/08/09 23:50:18  platin
 *   - implemented site localized extra population dynamics in TNL.
 *   - implemented the DIAGONALDYNAMICS method for population-only
 *     incoherent dynamics.
 *
 * Revision 1.18  2007/06/23 00:54:23  platin
 *
 *   - add support for index adjustment so that poptrans terms can be
 *     consistent when level energy crossover happens at disordered jobs.
 *
 * Revision 1.17  2007/06/14 18:07:01  platin
 *
 *   - add keyword "RSEED" that allows using a consistent random sequence
 *     in the MC sampling. This should help with the consistency of the
 *     results across jobs.
 *
 * Revision 1.16  2007/06/01 17:58:55  platin
 *
 *   - improved the Gaussian-Hermite module; now users can select between
 *     1-20 point rules can be used.
 *   - added keyword to select between MC and GH samplings.
 *
 * Revision 1.15  2007/03/09 08:58:56  platin
 *
 *   - more supporting gsl_* functions.
 *   - remove obsoleted code in params.c.
 *   - minor changes in gauss-hermite*
 *
 * Revision 1.14  2007/02/21 07:50:07  platin
 *
 *   - typo fixed.
 *
 * Revision 1.13  2007/02/21 06:12:26  platin
 *
 *   - remove leftover debug line.
 *
 * Revision 1.12  2007/01/18 18:50:18  platin
 *
 *   - bug fix.
 *
 * Revision 1.11  2006/12/22 00:30:03  platin
 *
 *     - support nonlinear pairwise rxn via Linblad formalism.
 *
 * Revision 1.10  2006/11/02 21:55:44  platin
 *
 *   - remove unnessary .h dependency.
 *
 * Revision 1.9  2006/10/27 23:02:59  platin
 *
 *   - commit modifications up to 10/27/2006.
 *
 * Revision 1.8  2006/08/26 05:01:17  platin
 *
 *   - include primitive two-exciton state support for dynamics and
 *     optical transitions.
 *
 * Revision 1.7  2006/08/23 22:30:26  platin
 *
 *   - add support for two-exciton states in codes regarding H and
 *     transition dipoles. The "Assign" format for transition dipole
 *     input can also be used to include effects of excited state absorption.
 *
 *   - basic support for TESLIST keyword.
 *
 * Revision 1.6  2006/07/13 23:37:29  platin
 *
 *   - add keyword TPRINT that adjusts time step of output.
 *   - minor message changes.
 *
 * Revision 1.5  2006/06/29 00:20:13  platin
 *
 *   - implemented the 3P photon-echo signal module.
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
 * Revision 1.1.1.1  2006/05/24 00:42:19  platin
 *
 *   - initial import of the qdas package.
 *   - qdas stands for quantum dynamics and spectroscopy.
 *   - basic code inherited from the "lineshape" package.
 *
 *
 *
 */
