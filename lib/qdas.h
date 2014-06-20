/***************************************************
 * qdas.h
 * 
 * By Yuan-Chung Cheng <yuanchung@ntu.edu.tw>
 *
 * Header file for qdas.c, some global constants...
 *
 ***************************************************/

#ifndef _QDAS_H
#define _QDAS_H 1

#include "gsl/gsl_matrix.h"

/* Global constants */
#define QDAS_VERSION "0.9rc4"

/* buffer size should be bigger than many intrinsic 
   dimensions such as the size of the Hamiltonian */
#define BUFFER_SIZE (4096)

#ifndef TRUE
#define TRUE (1)
#endif
#ifndef FALSE
#define FALSE (0)
#endif

/* Some physical unit conversion constants */
// convert time from cm unit to fs
#define TIME_CM2FS (5309.1)

/* a list of methods for static disorder sampling. */
#define QDAS_SDMETHOD_MC (0)
#define QDAS_SDMETHOD_GH (1)

/* a list of methods for field polarization sampling. */
#define QDAS_POLMETHOD_MC (0)

/* this constant is related to the ad-hoc population-only 
   dynamics for TNL dynamics */
#define QDAS_SET_DIAGONAL_DYNAMICS 12

/* type def */

/* basic 3D coordinate structure */
typedef struct {
  double x,y,z;
} coord3d;

#define VIBMODE_MAX 512
/* structure for vibration modes assigned to each excitations */
typedef struct {
  size_t n; /* number of modes */
  double omega[VIBMODE_MAX]; /* mode frequency */
  double S[VIBMODE_MAX]; /* Huang-Rhys factor */
} vibmodes;

#define PULSES_MAX 16
/* Assignments of Gaussian pulses, note the labeling
   does not correspond to time ordering, instead, the
   labeling represents the k-vector numbering */
typedef struct {
  double E0; /* pulse field strength */
  double tau0; /* center of pulse in time domain */
  double fwhm_tau; /* FWHM duration of the pulse */
  double w0; /* center of pulse in cm^-1 */
  /* note tha pulse polarization vectors are in the "polar" structure */
} pulse;

/* linear polarization angles of pulses and detector in a three-pulse experiment */
typedef struct {
  double angle[4]; /* polarization angles of pulse 1-4 */
  /* note that pulse 4 is actualy the detector  */
} polar_3ppe;

/* detector properties */
typedef struct {
  double upper;
  double lower; /* upper and lower energy limits for detector; in cm^-1 */
  /* note that detector polarization vector is in the "polar" structure */
} detector;

#define DIPOLES_MAX 1024
/* dipole vector for transition n->m */
typedef struct {
  size_t n,m;
  double vec[3]; /* dipole vector */
} dipole;

#define SINMODUL_MAX 256
/* Jm*sin(t/Tm+Phi) modulation on Hamiltonian element n,m */
typedef struct {
  size_t n,m;
  double Jm; /* amplitude */
  double Tm; /* sin(2*pi*t/Tm+Phi) */
  double Phi; /* phase factor */
} sinmodul;

/* two-exciton states (TES); note we only support
   well defined doubly excited states in the localized
   basis; so a two-exciton state is defined by two
   state numbers |nm> means excitation on chromophore n and m */
#define TE_STATE_MAX 1024
typedef struct {
  /* read: two-exciton state labeled as "label" is
     made up by excitations on site1 and site2 */
  size_t label;
  size_t site1,site2;
} te_state;

/* maximal number of extra population transfer terms */
#define POP_TRANS_MAX 16
typedef struct {
  size_t src,dest; /* population transfer from exciton state src to dest */
  double k; /* rate constant (in 1/fs), the 
	       rate input in the keyword file is the time constant in fs */
} poptrans;

/* maximal number of extra site-representation population transfer terms */
#define SITEPOP_TRANS_MAX 16
typedef struct {
  size_t src,dest; /* population transfer from site state src to dest */
  double k; /* rate constant (in 1/fs), the 
	       rate input in the keyword file is the time constant in fs */
  gsl_matrix_complex *Op1; /* this will be used for the Lindblad operator Ak
			      in the eigenbasis */
  gsl_matrix_complex *Op2; /* this will be used for Ak^\dagger*Ak
			      in the eigenbasis */
} sitepoptrans;

/* Non-linear population transfer terms; these are pairwise "annihilation" rxns
   A + A -> B
   na(t) = - k*na(t)*na(t)
   nb(t) = 1/2 * k*na(t)*na(t)
 */
#define NLPOP_TRANS_MAX 16
typedef struct {
  size_t src,dest; /* population transfer from exciton state src to dest */
  double k; /* rate constant (in 1/fs), the 
	       rate input in the keyword file is the time constant in fs */
} nlpoptrans;

/* bath terms; note that up to 8 groups of "BATH" keywords are
   accepted.. */
#define BATHGROUP_MAX 8
typedef struct {
  int bath_mod; /* bath module, this also depends on the 
		   lineshape or dynamics module selected */
  size_t bath_nparams;
  double bath_params[BUFFER_SIZE]; /* parameters defining the bath spectral function,
				      will be pass to the bath module */
} bathgroup_func;

/* non-local bath terms requires specific assignment of 
   system-bath operators */
#define NLBATH_MAX 256
typedef struct {
  double gamma; /* gamma is a scaling factor */
  double wc; /* cut-off frequency */
  gsl_matrix *S; /* system operator for system-bath interacting term */
} nlbath_func;

/* general data in the .key file; this structure defines
   everything, each individual module has to examine/use
   required elements */
typedef struct {

  /* Parameters that defines the multilevel electronic system */ 
  size_t nsize; /* size of sites/system, N */
  gsl_matrix *He; /* Hamiltonian of electronic excited-states, NxN real matrix */

  /* transition dipole moments of site localized states */
  size_t ndipoles;
  dipole mu[DIPOLES_MAX]; 

  /* sin() modulations; will modulate Hamiltonian matrix element
     using Jm*sin(t/Tm) */
  size_t nsinmoduls;
  sinmodul smodul[SINMODUL_MAX];

  /* labeling two-exciton states */
  size_t ntes;
  te_state tes_list[TE_STATE_MAX];

  /* automatically construct two-exciton states */
  size_t nctes_elems; /* number of elements (1 exciton states) that will
			 be used to construct two-exciton states */
  size_t ctes_elems_list[BUFFER_SIZE]; /* list of 1 exciton states */

  /* Parameters that defines the harmonic bath and vibrations */ 
  double beta; /* inverse temperature */

  /* bath groups; each represents a BATH keyword */
  size_t bathgroup_nterms;
  bathgroup_func bathgroup[BATHGROUP_MAX];

  /* additional non-local baths */
  size_t nlbath_nterms;
  nlbath_func nlbath[NLBATH_MAX];

  vibmodes *vibrations; /* vibrational modes */

  /* Parameters that defines the static disorders */ 
  gsl_matrix *disorder; /* sdev of gaussian static disorder for each He 
			   elements, NxN real matrix */
  size_t sdmethod; /* method for static disorder sampling; see the list defined above */
  size_t sdm_mc_niter; /* number of iterations for the MC average method */
  size_t sdm_mc_rseed; /* seed for random number generator. */
  size_t sdm_gh_order; /* number of points for Gauss-Hermite integration in each dim. */
  double sdm_gh_wmin; /* minimum weights accepted for Gauss-Hermite average */

  /* lineshape module */
  //  char *lineshape_mod; /* method/module that will be used to compute time-domain lineshape */
  double spec_start; /* how to print the simulated spectrum */
  double spec_end;
  double spec_step;

  /* dynamics and spec. module */
  size_t npulses; /* number of pulses */
  pulse pulse_seq[PULSES_MAX]; /* Pulse sequence */
  detector detect; /* detector properties */

  polar_3ppe *polar; /* pulse polarization angles; NULL pointer means
			  pulse polarizations will be ignored */
  size_t polmethod; /* method for polarization sampling; see the list defined above */
  size_t polm_mc_niter; /* number of iterations for the MC polarization average method */
  size_t polm_mc_rseed; /* seed for polarization random number generator. */
  double align; /* alignment factor along the z-axis; for anisotropic samples */

  size_t nptrans; /* number of population transfer terms */
  poptrans extra_ptrans[POP_TRANS_MAX]; /* extra population transfer terms */
  size_t nsiteptrans; /* number of  site-representation population transfer terms */
  sitepoptrans extra_siteptrans[SITEPOP_TRANS_MAX]; /* extra population transfer terms */
  size_t nnlptrans; /* number of NL population transfer terms */
  nlpoptrans extra_nlptrans[POP_TRANS_MAX]; /* extra non-linear population transfer terms */
  double tstop; /* how to propagate the system in time */
  double tstep;
  double tprint; /* print output in each tprint fs */

  int printlv; /* output detail level, the higher the more details got printed */

  double cnconv; /* convergence threshold for Crank-Nicholson iterations */

  int diagonaldynamics; /* trun on the ad-hoc population-only dynamics for TNL module
			   if(keys->diagonaldynamics == QDAS_SET_DIAGONAL_DYNAMICS) */
} qdas_keys;

/* Global variables */

#endif /* qdas.h */

/*
 * $Log$
 * Revision 1.18  2007/08/09 23:50:35  platin
 *   - implemented site localized extra population dynamics in TNL.
 *   - implemented the DIAGONALDYNAMICS method for population-only
 *     incoherent dynamics.
 *
 * Revision 1.17  2007/06/23 00:54:01  platin
 *
 *   - bump to 0.8p3.
 *   - add support for index adjustment so that poptrans terms can be
 *     consistent when level energy crossover happens at disordered jobs.
 *
 * Revision 1.16  2007/06/21 22:11:29  platin
 *
 *   - code clean up.
 *
 * Revision 1.15  2007/06/14 18:07:42  platin
 *
 *   - bump to 0.8p2 for various minor changes.
 *
 * Revision 1.14  2007/06/01 17:58:39  platin
 *
 *   - bump to 0.8p1
 *   - improved the Gaussian-Hermite module; now users can select between
 *     1-20 point rules can be used.
 *   - added keyword to select between MC and GH samplings.
 *
 * Revision 1.13  2007/03/26 19:15:42  platin
 *
 *   - bump to v0.8.
 *
 * Revision 1.12  2007/03/09 08:57:58  platin
 *
 *   - bump version.
 *
 * Revision 1.11  2006/12/22 00:29:52  platin
 *
 *   - support nonlinear pairwise rxn via Linblad formalism.
 *
 * Revision 1.10  2006/10/30 06:33:01  platin
 *
 *   - change source of parser lib.
 *   - bump version to 0.5.
 *
 * Revision 1.9  2006/10/27 23:02:46  platin
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
 * Revision 1.6  2006/08/15 23:05:06  platin
 *
 *   - remove ad-hoc fix in tnl-* dynamics; now both programs correctly
 *     handle renormalization term (Hren) and eliminate erroneous dynamics,
 *     such as nonpositive dynamics and non-zero long time coherence terms...
 *
 * Revision 1.5  2006/07/15 07:52:12  platin
 *
 *   - minor changes.
 *
 * Revision 1.4  2006/07/13 23:37:29  platin
 *
 *   - add keyword TPRINT that adjusts time step of output.
 *   - minor message changes.
 *
 * Revision 1.3  2006/06/28 17:14:11  platin
 *
 *   - add the dm3pes module that computes the three-pulse photon-echo
 *     signals using Domcke's density-matrix based method.
 *   - aux functions used to implement dm3pes.
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
