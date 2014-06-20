/***************************************************
 * bath_ohmart.c
 *
 * Bath functions that implement Meier and Tannor's
 * parameterization of Ohmic spectral density with 
 * exponential cutoff.
 *
 * By Yuan-Chung Cheng <yuanchung@ntu.edu.tw>
 *
 *
 ***************************************************/

#define _GNU_SOURCE
//#define HAVE_INLINE

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
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_blas.h>

#include "qdas.h"
#include "aux.h"
#include "bath_ohmart.h"
#include "bathtype_mt99art.h"

// Defined Constants

#define NTERMS_MAX (256)

/* Parameterization of the Ohmic correlation functions that's
   suitable for Meier and Tannor's TNL calculations
   J(w)=gamma*w*exp(-w/wc)
   C_r(t) = 1/Pi * \int_0^\infy dw J(w)*coth(beta*w/2)*cos(w*t),
   C_i(t) = 1/Pi * \int_0^\infy dw J(w)*sin(w*t).
          = a(t) - ib(t)
   where

   a(t) = \sum_k alpha_k*exp(gamma_k*t)
   b(t) = \sum_k alpha_k*exp(gamma_k*t)

*/

#define NTERMS_Wc40T298K_R (4)
const double ohmart_At_wc40T298K_params[NTERMS_Wc40T298K_R][3] =  // @ 298K and wc=40 cm^-1
  { /* alpha_k,   Re(gamma_k), Im(gamma_k) */
    { 6523.82,  -48.7526,  0.0},
    {  900.086,  -12.4042,  0.0},
    {-1055.5,  -151.931,  101.771},
    {-1055.5,  -151.931, -101.771}};

#define NTERMS_Wc40T298K_I (5)
const double ohmart_Bt_wc40T298K_params[NTERMS_Wc40T298K_I][3] =  // @ 298K and wc=40 cm^-1
  { /* alpha_k,   Re(gamma_k), Im(gamma_k) */
    { -92.455,  -33.8384,   31.1647},
    { -92.455,  -33.8384,  -31.1647},
    {265.5385,  -110.199,  108.749},
    {265.5385,  -110.199, -108.749},
    {-343.251,  -28.4544,    0.0}};

#define NTERMS_Wc40T77K_R (4)
const double ohmart_At_wc40T77K_params[NTERMS_Wc40T77K_R][3] =  // @ 77K and wc=40 cm^-1
  { /* alpha_k,   Re(gamma_k), Im(gamma_k) */
    {  23633.1698208212/2.0, -92.9416444359049,  0.147416660335236},
    {  23633.1698208212/2.0, -92.9416444359049, -0.147416660335236},
    {  -22640.2546449242, -97.6879549363807, 0.0},
    {  485.244352702437, -17.6817762327653, 0.0}};

#define NTERMS_Wc40T77K_I (5)
const double ohmart_Bt_wc40T77K_params[NTERMS_Wc40T77K_I][3] =  // @ 77K and wc=40 cm^-1
  { /* alpha_k,   Re(gamma_k), Im(gamma_k) */
    { -5757.92204560217/2.0, -126.206533479494,  0.161564408927534},
    { -5757.92204560217/2.0, -126.206533479494, -0.161564408927534},
    {  6043.31365909399/2.0, -127.822966309534, 41.5093105580789},
    {  6043.31365909399/2.0, -127.822966309534, -41.5093105580789},
    { -285.391613491813, -29.7028573890488, 0.0}};

#define NTERMS_Wc100T77K_R (5)
const double ohmart_At_wc100T77K_params[NTERMS_Wc100T77K_R][3] =  // @ 77K and wc=100 cm^-1
  { /* alpha_k,   Re(gamma_k), Im(gamma_k) */
    {  12443.6863218718,  -366.071568410724,  0.936604846614617},
    {  12443.6863218718,  -366.071568410724, -0.936604846614617},
    { -10670.5415377028,  -431.449141105434,  139.187403982732},
    { -10670.5415377028,  -431.449141105434, -139.187403982732},
    {  1278.30209120449,  -49.2518957854801, 0.0}};

#define NTERMS_Wc100T298K_R (4)
const double ohmart_At_wc100T298K_params[NTERMS_Wc100T298K_R][3] =  // @ 298K and wc=100 cm^-1
  { /* alpha_k,   Re(gamma_k), Im(gamma_k) */
    {  16250.1,  -121.014,  0.0},
    {  2076.48,  -29.3882,  0.0},
    { -2313.26,  -419.346,  300.032},
    { -2313.26,  -419.346, -300.032}};

#define NTERMS_Wc120T77K_R (5)
const double ohmart_At_wc120T77K_params[NTERMS_Wc120T77K_R][3] =  // @ 77K and wc=120 cm^-1
  { /* alpha_k,   Re(gamma_k), Im(gamma_k) */
    {  12876.8863535821,   -419.816921778719,   0.25279760068966},
    {  12876.8863535821,   -419.816921778719,  -0.25279760068966},
    { -10236.9429716212,   -531.356937443002,   208.184518003959},
    { -10236.9429716212,   -531.356937443002,  -208.184518003959},
    {  1086.15878524594,   -49.6282482683452,  0.0}};

#define NTERMS_Wc150T77K_R (5)
const double ohmart_At_wc150T77K_params[NTERMS_Wc150T77K_R][3] =  // @ 77K and wc=150 cm^-1
  { /* alpha_k,   Re(gamma_k), Im(gamma_k) */
    {  13657.5813615899,  -467.291739780798,  72.9956905702684},
    {  13657.5813615899,  -467.291739780798, -72.9956905702684},
    {  -9450.56281771260,  -672.566891876374, 337.342206189664},
    {  -9450.56281771260,  -672.566891876374,-337.342206189664},
    {  704.341141979774,  -42.6533694369846, 0.0}};

#define NTERMS_Wc150T100K_R (4)
const double ohmart_At_wc150T100K_params[NTERMS_Wc150T100K_R][3] =  // @ 100K and wc=150 cm^-1
  { /* alpha_k,   Re(gamma_k), Im(gamma_k) */
    {  27267.1565787438,  -460.389383779187,   0.0},
    { -18900.6495899502/2.0,  -673.398319375812, 309.235251478611},
    { -18900.6495899502/2.0,  -673.398319375812,-309.235251478611},
    {  1748.22727249094,  -60.2622979556219,  0.0}};

#define NTERMS_Wc150T150K_R (4)
const double ohmart_At_wc150T150K_params[NTERMS_Wc150T150K_R][3] =  // @ 150K and wc=150 cm^-1
  { /* alpha_k,   Re(gamma_k), Im(gamma_k) */
    {  27632.4729345049,  -357.953337785357,  0.0},
    { -18496.6050517295/2.0,  -592.676072606874,  2.05837951935921},
    { -18496.6050517295/2.0,  -592.676072606874, -2.05837951935921},
    {  3384.79561318171,  -66.017262704827,  0.0}};

#define NTERMS_Wc150T200K_R (4)
const double ohmart_At_wc150T200K_params[NTERMS_Wc150T200K_R][3] =  // @ 200K and wc=150 cm^-1
  { /* alpha_k,   Re(gamma_k), Im(gamma_k) */
    {  28353.8808192245,  -310.832458437857,  0.0},
    {  -17675.4083737494/2.0,  -552.568992300968,   2.22008328593291},
    {  -17675.4083737494/2.0,  -552.568992300968,  -2.22008328593291},
    {  4699.58266357656,  -66.7839857926646,  0.0}};

#define NTERMS_Wc150T250K_R (4)
const double ohmart_At_wc150T250K_params[NTERMS_Wc150T250K_R][3] =  // @ 250K and wc=150 cm^-1
  { /* alpha_k,   Re(gamma_k), Im(gamma_k) */
    {  29378.1932289507,  -272.516260043532,  0.0},
    {  -16485.543608721/2.0,  -545.697289265353,   2.43479915244978},
    {  -16485.543608721/2.0,  -545.697289265353,  -2.43479915244978},
    {  5442.50169576429,  -64.2142176766774,  0.0}};

#define NTERMS_Wc150T300K_R (3)
const double ohmart_At_wc150T300K_params[NTERMS_Wc150T300K_R][3] =  // @ 300K and wc=150 cm^-1
  { /* alpha_k,   Re(gamma_k), Im(gamma_k) */
    {  28426.2654667669,  -225.50461995441,  0.0},
    { -12206.3313935179,  -628.959878567182,  0.0},
    {  5084.89788316839,  -56.8411674422759,  0.0}};

#define NTERMS_Wc200T77K_R (5)
const double ohmart_At_wc200T77K_params[NTERMS_Wc200T77K_R][3] =  // @ 77K and wc=200 cm^-1
  { /* alpha_k,   Re(gamma_k), Im(gamma_k) */
    {  12090.0484320200,  -403.263337489127,  212.470144060889},
    {  12090.0484320200,  -403.263337489127, -212.470144060889},
    { -4855.26833686285,  -969.634834028302,  724.384579770432},
    { -4855.26833686285,  -969.634834028302, -724.384579770432},
    {  427.124303849066,  -37.1258539341996, 0.0}};

/* Now fitted imaginary part for each cut-off frequencies;
   note that the imaginary part of the correlation is 
   temperature independent */

#define NTERMS_Wc100_I (5)
const double ohmart_Bt_wc100_params[NTERMS_Wc100_I][3] =  // @ wc=100 cm^-1
  { /* alpha_k,   Re(gamma_k), Im(gamma_k) */
    {-930.955,  -92.4259,  63.7218},
    {-930.955,  -92.4259, -63.7218},
    { 1687.22,  -272.607,  268.786},
    { 1687.22,  -272.607, -268.786},
    {-1512.53,  -62.6094,    0.0}};

#ifdef NEVER_DEFINED
// the set above seems to be better at the tail part...
const double ohmart_Bt_wc100_params[NTERMS_Wc100_I][3] =  // @ wc=100 cm^-1
  { /* alpha_k,   Re(gamma_k), Im(gamma_k) */
    { -2448.82491591573,  -288.366283909217,  1.02490600213789},
    { -2448.82491591573,  -288.366283909217, -1.02490600213789},
    {  3840.16549949778,  -297.187818100805,  202.596872614433},
    {  3840.16549949778,  -297.187818100805, -202.596872614433},
    { -2782.68116716410,  -87.893592837696, 0.0}};
#endif

#define NTERMS_Wc120_I (5)
const double ohmart_Bt_wc120_params[NTERMS_Wc120_I][3] =  // @ wc=120 cm^-1
  { /* alpha_k,   Re(gamma_k), Im(gamma_k) */
    { -209.131824263686,   -71.3069098857846,   122.914894389386},
    { -209.131824263686,   -71.3069098857846,  -122.914894389386},
    {  2433.28488928941,   -328.896365011599,   327.059175954643},
    {  2433.28488928941,   -328.896365011599,  -327.059175954643},
    { -4448.30613005145,   -102.350485884965,  0.0}};

#define NTERMS_Wc150_I (5)
const double ohmart_Bt_wc150_params[NTERMS_Wc150_I][3] =  // @ wc=150 cm^-1
  { /* alpha_k,   Re(gamma_k), Im(gamma_k) */
    { -2210.62193900496,  -137.487359661568,  93.2130852433601},
    { -2210.62193900496,  -137.487359661568, -93.2130852433601},
    {  3806.62635423718,  -412.051731007188, 403.458011149567},
    {  3806.62635423718,  -412.051731007188,-403.458011149567},
    { -3192.00883046444,  -90.8248736586544, 0.0}};

#define NTERMS_Wc200_I (5)
const double ohmart_Bt_wc200_params[NTERMS_Wc200_I][3] =  // @ wc=200 cm^-1
  { /* alpha_k,   Re(gamma_k), Im(gamma_k) */
    { -3658.79319122295,  -181.347994217028,  129.146120449445},
    { -3658.79319122295,  -181.347994217028, -129.146120449445},
    {  6692.09724882265,  -549.066387239708,  540.036800836489},
    {  6692.09724882265,  -549.066387239708, -540.036800836489},
    { -6066.60811519935,  -124.418778267288, 0.0}};

/* initialize f using assigned temperature and parameters */
void bath_ohmart_ct_init(mt99art_bath_func *f, 
			  const double beta, const double gamma, const double wc,
			  const gsl_matrix *S)
{
  int i,j;
  int nterms_r=0,nterms_i=0;

  double bath_At_params[NTERMS_MAX][3],bath_Bt_params[NTERMS_MAX][3];

  gsl_complex vren;

  /* FIXME: implement auto-fitting functions!! 
     at this point, we only support hand-fitted C(t) */
  if(fabs(1.4387/beta-298.0)<2.0 && fabs(wc-40.0)<2.0) {
    // accept 298+-2K and wc=40
    nterms_r=NTERMS_Wc40T298K_R;
    for(i=0;i<nterms_r;i++) {
      for(j=0;j<nterms_r;j++) {
	bath_At_params[i][j]=ohmart_At_wc40T298K_params[i][j];
      }
    }
    nterms_i=NTERMS_Wc40T298K_I;
    for(i=0;i<nterms_i;i++) {
      for(j=0;j<nterms_i;j++) {
	bath_Bt_params[i][j]=ohmart_Bt_wc40T298K_params[i][j];
      }
    }
  } else if(fabs(1.4387/beta-77.0)<2.0 && fabs(wc-40.0)<2.0) {
    // accept 77+-2K and wc=40
    nterms_r=NTERMS_Wc40T77K_R;
    for(i=0;i<nterms_r;i++) {
      for(j=0;j<nterms_r;j++) {
	bath_At_params[i][j]=ohmart_At_wc40T77K_params[i][j];
      }
    }
    nterms_i=NTERMS_Wc40T77K_I;
    for(i=0;i<nterms_i;i++) {
      for(j=0;j<nterms_i;j++) {
	bath_Bt_params[i][j]=ohmart_Bt_wc40T77K_params[i][j];
      }
    }
  } else if((
	     fabs(1.4387/beta-77.0)<2.0 || // 77K
	     fabs(1.4387/beta-298.0)<2.0 // 298K
	     ) && fabs(wc-100.0)<2.0) {
    // Wc = 100
    if(fabs(1.4387/beta-77.0)<2.0) {
      // 77K
      nterms_r=NTERMS_Wc100T77K_R;
      for(i=0;i<nterms_r;i++) {
	for(j=0;j<nterms_r;j++) {
	  bath_At_params[i][j]=ohmart_At_wc100T77K_params[i][j];
	}
      }
    } else if(fabs(1.4387/beta-298.0)<2.0) {
      // 298K
      nterms_r=NTERMS_Wc100T298K_R;
      for(i=0;i<nterms_r;i++) {
	for(j=0;j<nterms_r;j++) {
	  bath_At_params[i][j]=ohmart_At_wc100T298K_params[i][j];
	}
      }
    }
    /* imaginary part is temperature independent */
    nterms_i=NTERMS_Wc100_I;
    for(i=0;i<nterms_i;i++) {
      for(j=0;j<nterms_i;j++) {
	bath_Bt_params[i][j]=ohmart_Bt_wc100_params[i][j];
      }
    }
    // end the wc=100 cases
  } else if((
	     fabs(1.4387/beta-77.0)<2.0 // 77K
	     ) && fabs(wc-120.0)<2.0) {
    // Wc = 120
    if(fabs(1.4387/beta-77.0)<2.0) {
      // 77K
      nterms_r=NTERMS_Wc120T77K_R;
      for(i=0;i<nterms_r;i++) {
	for(j=0;j<nterms_r;j++) {
	  bath_At_params[i][j]=ohmart_At_wc120T77K_params[i][j];
	}
      }
    }
    /* imaginary part is temperature independent */
    nterms_i=NTERMS_Wc120_I;
    for(i=0;i<nterms_i;i++) {
      for(j=0;j<nterms_i;j++) {
	bath_Bt_params[i][j]=ohmart_Bt_wc120_params[i][j];
      }
    }
    // end the wc=120 cases
  } else if((
	     fabs(1.4387/beta-77.0)<2.0 || // 77K
	     fabs(1.4387/beta-100.0)<2.0 || // 100K
	     fabs(1.4387/beta-150.0)<2.0 || // 150K
	     fabs(1.4387/beta-200.0)<2.0 || // 200K
	     fabs(1.4387/beta-250.0)<2.0 || // 250K
	     fabs(1.4387/beta-300.0)<2.0 // 300K
	     ) && fabs(wc-150.0)<2.0) {
    // Wc = 150
    if(fabs(1.4387/beta-77.0)<2.0) {
      // 77K
      nterms_r=NTERMS_Wc150T77K_R;
      for(i=0;i<nterms_r;i++) {
	for(j=0;j<nterms_r;j++) {
	  bath_At_params[i][j]=ohmart_At_wc150T77K_params[i][j];
	}
      }
    } else if(fabs(1.4387/beta-100.0)<2.0) {
      // 100K
      nterms_r=NTERMS_Wc150T100K_R;
      for(i=0;i<nterms_r;i++) {
	for(j=0;j<nterms_r;j++) {
	  bath_At_params[i][j]=ohmart_At_wc150T100K_params[i][j];
	}
      }
    } else if(fabs(1.4387/beta-150.0)<2.0) {
      // 150K
      nterms_r=NTERMS_Wc150T150K_R;
      for(i=0;i<nterms_r;i++) {
	for(j=0;j<nterms_r;j++) {
	  bath_At_params[i][j]=ohmart_At_wc150T150K_params[i][j];
	}
      }
    } else if(fabs(1.4387/beta-200.0)<2.0) {
      // 200K
      nterms_r=NTERMS_Wc150T200K_R;
      for(i=0;i<nterms_r;i++) {
	for(j=0;j<nterms_r;j++) {
	  bath_At_params[i][j]=ohmart_At_wc150T200K_params[i][j];
	}
      }
    } else if(fabs(1.4387/beta-250.0)<2.0) {
      // 250K
      nterms_r=NTERMS_Wc150T250K_R;
      for(i=0;i<nterms_r;i++) {
	for(j=0;j<nterms_r;j++) {
	  bath_At_params[i][j]=ohmart_At_wc150T250K_params[i][j];
	}
      }
    } else if(fabs(1.4387/beta-300.0)<2.0) {
      // 300K
      nterms_r=NTERMS_Wc150T300K_R;
      for(i=0;i<nterms_r;i++) {
	for(j=0;j<nterms_r;j++) {
	  bath_At_params[i][j]=ohmart_At_wc150T300K_params[i][j];
	}
      }
    }
    /* imaginary part is temperature independent */
    nterms_i=NTERMS_Wc150_I;
    for(i=0;i<nterms_i;i++) {
      for(j=0;j<nterms_i;j++) {
	bath_Bt_params[i][j]=ohmart_Bt_wc150_params[i][j];
      }
    }
    // end the wc=150 cases
  } else if((
	     fabs(1.4387/beta-77.0)<2.0 // 77K
	     ) && fabs(wc-200.0)<2.0) {
    // Wc = 200
    if(fabs(1.4387/beta-77.0)<2.0) {
      // 77K
      nterms_r=NTERMS_Wc200T77K_R;
      for(i=0;i<nterms_r;i++) {
	for(j=0;j<nterms_r;j++) {
	  bath_At_params[i][j]=ohmart_At_wc200T77K_params[i][j];
	}
      }
    }
    /* imaginary part is temperature independent */
    nterms_i=NTERMS_Wc200_I;
    for(i=0;i<nterms_i;i++) {
      for(j=0;j<nterms_i;j++) {
	bath_Bt_params[i][j]=ohmart_Bt_wc200_params[i][j];
      }
    }
    // end the wc=200 cases
  } else {
    printf("Sorry, parameterizations at T=%.2fK and wc=%.2fcm-1 not supported yet!!\n",1.4387/beta,wc);
    exit(EXIT_FAILURE);
  }

  // assign gamma, beta, wc, S
  f->beta=beta;
  f->S=gsl_matrix_alloc(S->size1,S->size1);
  gsl_matrix_memcpy(f->S,S);

#ifdef BATH_DEBUG
  printf("\n");
  printf("Initializing bath correlation function ...\n");
  printf("  Beta = %f cm (%5.2f K)\n",beta,1.4387/beta);
  printf("  Ohmic spectral sensity:\n");
  printf("    Wc = %f\n",wc);
  printf("    Coupling Strength Gamma = %f\n",gamma);
  printf("\n");
#endif

  /* the real part, a(t) */
  // totally nterms terms from the parameterization
  for(i=0;i<nterms_r;i++) {
    f->Ar[i]=gsl_complex_rect(gamma*bath_At_params[i][0],0.0);
    f->Gr[i]=gsl_complex_rect(bath_At_params[i][1],bath_At_params[i][2]);
  }
  /* the negative imaginary part, b(t) */
  for(i=0;i<nterms_i;i++) {
    f->Ai[i]=gsl_complex_rect(gamma*bath_Bt_params[i][0],0.0);
    f->Gi[i]=gsl_complex_rect(bath_Bt_params[i][1],bath_Bt_params[i][2]);
  }

  if(gamma>1e-12) {
    f->Nr=nterms_r;
    f->Ni=nterms_i;
  } else {
    f->Nr=0;
    f->Ni=0;
  }

  /* for consistency, we compute the renormalization 
     energy from the imaginary parameterization. The results 
     should be close to mu=2.16*gamma*wc/M_PI for JS03 spectral density */
  /* mu = 2 * \sum_k Ai[k]/Gi[k] (remember the extra minus 
     sign in the definition of Ai[k] here.) */
  vren=gsl_complex_rect(0.0,0.0);
  for(i=0;i<nterms_i;i++) {
    vren=gsl_complex_add(vren,gsl_complex_div(f->Ai[i],f->Gi[i]));
  }
  vren=gsl_complex_mul(vren,gsl_complex_rect(2.0,0.0));
  f->Vren=GSL_REAL(vren);

#ifdef BATH_DEBUG
  printf("  M&T's parameterization:\n");
  printf("\n");
  printf("Vren = %f\n",f->Vren);
  printf("\n");
  printf("  a(t) = \\sum_k alpha_k*exp(gamma_k*t)\n");
  printf("\n");
  printf("    %30s %30s\n","alpha_k","gamma_k");
  for(i=0;i<f->Nr;i++) {
    printf("    %14f%+14f   %14f%+14f\n",
	   GSL_REAL(f->Ar[i]),GSL_IMAG(f->Ar[i]),GSL_REAL(f->Gr[i]),GSL_IMAG(f->Gr[i]));
  }
  printf("\n");
  printf("\n");
  printf("  b(t) = \\sum_k alpha_k*exp(gamma_k*t)\n");
  printf("\n");
  printf("    %30s %30s\n","alpha_k","gamma_k");
  for(i=0;i<f->Ni;i++) {
    printf("    %14f%+14f   %14f%+14f\n",
	   GSL_REAL(f->Ai[i]),GSL_IMAG(f->Ai[i]),GSL_REAL(f->Gi[i]),GSL_IMAG(f->Gi[i]));
  }
  printf("\n");
  printf("\n");
#endif

  // done!
}

/* interface functions */
int bath_ohmart_init_params(const size_t nsize, const double beta, 
			     const size_t bath_nparams, const double *bath_params)
{
  int i,idx;
  int nterms;
  size_t ndim;

  double gamma[256],wc[256];

  int opbra[256];
  int opket[256];
  gsl_matrix *S;

  ndim=nsize;
  S=gsl_matrix_alloc(ndim,ndim);

  // input in the mod_params array should be in a table form
  // n m gamma wc
  // parse them to the global variables
  nterms=bath_nparams/4; // number of simple system-bath coupling terms
  for(i=0;i<nterms;i++) {
    int na,nb;
    idx=i*4;
    // system OP represented by ket state |n> and bra state <m|, i.e. |n><m|.
    // We make sure n<=m
    na=(int)bath_params[idx];
    nb=(int)bath_params[idx+1];
    opket[i]=((na < nb) ? na : nb);
    opbra[i]=((na < nb) ? nb : na);
    gamma[i]=bath_params[idx+2];
    wc[i]=bath_params[idx+3];
    if(opbra[i]<=0 || opket[i]<=0) {
      printf("Error: can't parse system-bath interaction terms.\n");
      exit(EXIT_FAILURE);
    }
    /* construct the S operator */
    gsl_matrix_set_zero(S);
    gsl_matrix_set(S,opket[i]-1,opbra[i]-1,1.0);
    gsl_matrix_set(S,opbra[i]-1,opket[i]-1,1.0);
    bath_ohmart_ct_init(BATH_MT99ARTBathFunc+BATH_MT99ART_NOps+i,beta,gamma[i],wc[i],S);
  }

  /* show what we have just done */
  printf("Use artificial C(t) parameterization of Ohmic J(w):\n");
  printf("J(w)=gamma*exp(-w/wc)\n");
  printf("\n");
  printf("System-Op and J(w) for each site:\n");
  printf("\n");
  printf("  %16s %12s %12s %12s\n","System-Op","gamma","wc (cm^-1)","(Nr,Ni)");
  for(i=0;i<nterms;i++) {
    // print out the parameters...
    if(opbra[i] == opket[i]) {
      // diagonal e-ph couplings
      printf("  |%d><%d|            %12.4f %12.4f (%4lu,%4lu)\n",
	     opket[i],opbra[i],gamma[i],wc[i],
	     BATH_MT99ARTBathFunc[BATH_MT99ART_NOps+i].Nr,BATH_MT99ARTBathFunc[BATH_MT99ART_NOps+i].Ni);
    } else {
      // off-diagonal terms, note that we will automatically 
      // symmetrize the system part in all evaluations later.
      printf("  |%d><%d|+|%d><%d|   %12.4f %12.4f (%4lu,%4lu)\n",
	     opket[i],opbra[i],
	     opbra[i],opket[i],
	     gamma[i],wc[i],
	     BATH_MT99ARTBathFunc[BATH_MT99ART_NOps+i].Nr,BATH_MT99ARTBathFunc[BATH_MT99ART_NOps+i].Ni);
    }
  }
  printf("\n");

  /* add to the number of MT99 bath terms */
  BATH_MT99ART_NOps+=nterms;

  gsl_matrix_free(S);
  return 0;
}

#ifdef MAIN
int main()
{
  // unit in cm^-1
  double beta=0.018685; // 77K
  //  double beta=1.4387/300.0;
  //  double beta=0.0048281; // 298K
  double gamma=1.0;
  double wc=120.0;

  double t;
  gsl_complex ct;
  mt99art_bath_func bf;
  int i=0;
  gsl_matrix *S;

  S=gsl_matrix_alloc(2,2);

  bath_ohmart_ct_init(&bf,beta,gamma,wc,S);
  t=0.0; // t in fs
  for(i=0;i<=2000;i++) {
    ct=bath_mt99art_ct(&bf,t/TIME_CM2FS);
    printf("t=%f, Ct=%f + %fI\n",t/TIME_CM2FS,GSL_REAL(ct),GSL_IMAG(ct));
    printf("t=%f, Ht=%f + %fI\n",t/TIME_CM2FS,
	   bath_mt99art_ht_r(&bf,t/TIME_CM2FS),bath_mt99art_ht_i(&bf,t/TIME_CM2FS));
    printf("t=%f, Gt=%f + %fI\n",t/TIME_CM2FS,
	   bath_mt99art_gt_r(&bf,t/TIME_CM2FS),bath_mt99art_gt_i(&bf,t/TIME_CM2FS));

    t=t+1.0;
  }

  return 1;

}
#endif

/*
 * $Log$
 * Revision 1.8  2007/07/18 23:59:23  platin
 *   - support Ohmic bat at wc=200 cm-1, 77K.
 *
 * Revision 1.7  2007-07-18 23:03:02  platin
 *
 *   - support wc=100, 77K.
 *   - clean up the checking code...
 *
 * Revision 1.6  2007/06/21 00:35:32  platin
 *
 *   - support functions for finding maximum value of the real (or imag) part
 *     of a gsl_matrix_complex.
 *   - support wc=40 and T=77K ohmic bath.
 *
 * Revision 1.5  2007/05/22 21:13:26  platin
 *
 *   - more temperature fitting for wc=150cm^-1 ohmic spectral density.
 *
 * Revision 1.4  2007/05/07 23:50:31  platin
 *
 *   - add parameter for ohmart, 77K wc=150.
 *
 * Revision 1.3  2007/03/12 23:08:27  platin
 *
 *   - fixed MT99 bath modules to allow more general bath functions,
 *     now baths not defined by gamma and wc can be used.
 *
 * Revision 1.2  2006/11/28 01:24:16  platin
 *
 *   - add parameter for ohmart bath at wc=100 and T=298K.
 *
 * Revision 1.1  2006/11/08 06:23:50  platin
 *
 *   - minor changes and adding a artificial ohmic bath module that
 *     uses more efficient artificial fit to Ohmic correlation
 *     functions.
 *
 *
 */
