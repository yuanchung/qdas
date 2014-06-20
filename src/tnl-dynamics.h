/***************************************************
 * tnl-dynamics.h
 * 
 * By Yuan-Chung Cheng <yuanchung@ntu.edu.tw>
 *
 * Header file for tnl-dynamics.c, some global constants...
 *
 ***************************************************/

#ifndef _TNL_DYNAMICS_H
#define _TNL_DYNAMICS_H 1

#include "gsl/gsl_matrix.h"
#include "qdas.h"

/* Global constants */

/* type def */


/* Global variables */

#endif /* tnl-dynamics.h */

/*
 * $Log$
 * Revision 1.5  2007/06/23 00:54:49  platin
 *   - bump to 0.8p3.
 *   - add support for index adjustment so that poptrans terms can be
 *     consistent when level energy crossover happens at disordered jobs.
 *
 * Revision 1.4  2007/06/14 18:08:46  platin
 *
 *   - support for RSEED.
 *   - bump to 0.8p2
 *
 * Revision 1.3  2007/06/01 18:00:02  platin
 *
 *   - bump to 0.8p1
 *   - improved the Gaussian-Hermite module; now users can select between
 *     1-20 point rules can be used.
 *   - added keyword to select between MC and GH samplings.
 *   - static disorder average loop in both tnl-dm3pes and tnl-dynamics
 *     is fixed.
 *
 * Revision 1.2  2007/03/26 19:15:25  platin
 *
 *   - bump to v0.8.
 *   - fix the pulse duration factor; the fwhm value is the field fwhm.
 *
 * Revision 1.1  2006/10/27 22:24:54  platin
 *
 *   - commit re-organized tnl-dm3pes and tnl-dynamics code.
 *
 * Revision 1.6  2006/08/26 05:01:17  platin
 *
 *   - include primitive two-exciton state support for dynamics and
 *     optical transitions.
 *
 * Revision 1.5  2006/08/15 23:05:06  platin
 *
 *   - remove ad-hoc fix in tnl-* dynamics; now both programs correctly
 *     handle renormalization term (Hren) and eliminate erroneous dynamics,
 *     such as nonpositive dynamics and non-zero long time coherence terms...
 *
 * Revision 1.4  2006/08/11 06:00:06  platin
 *
 *   - bump to 0.5p2
 *   - now propagates in S-picture, remove interacting picture code.
 *   - add a "Markovian-like" term in bath functions to preserve positivity.
 *
 * Revision 1.3  2006/07/31 03:27:33  platin
 *
 *   - pre-compute Sn(t) to avoid repeating evaluations; this gives
 *     a slight performance increase.
 *
 * Revision 1.2  2006/07/28 23:31:08  platin
 *
 *   - major change. Use short-time Crank-Nicholson propagator.
 *
 * Revision 1.1  2006/07/11 16:46:32  platin
 *
 *   - import the tnl-dynamics bit; still requires modification to take
 *     initial states etc.
 *
 *
 */
