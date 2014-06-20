/***************************************************
 * tnl-dm3pes.h
 * 
 * By Yuan-Chung Cheng <yuanchung@ntu.edu.tw>
 *
 * Header file for tnl-dm3pes.c, some global constants...
 *
 ***************************************************/

#ifndef _TNL_DM3PES_H
#define _TNL_DM3PES_H 1

#include "gsl/gsl_matrix.h"
#include "qdas.h"

/* Global constants */

/* type def */

/* Global variables */

#endif /* tnl-dm3pes.h */

/*
 * $Log$
 * Revision 1.6  2007/06/23 00:54:49  platin
 *   - bump to 0.8p3.
 *   - add support for index adjustment so that poptrans terms can be
 *     consistent when level energy crossover happens at disordered jobs.
 *
 * Revision 1.5  2007/06/14 18:08:46  platin
 *
 *   - support for RSEED.
 *   - bump to 0.8p2
 *
 * Revision 1.4  2007/05/31 01:29:46  platin
 *
 *   - implemented the more efficient verison in Egorova et al.
 *   - bump to 0.8p1.
 *
 * Revision 1.3  2007/03/26 19:15:25  platin
 *
 *   - bump to v0.8.
 *   - fix the pulse duration factor; the fwhm value is the field fwhm.
 *
 * Revision 1.2  2007/03/09 08:59:46  platin
 *
 *   - include Gauss-Hermite Quadrature code for static disorder of two DOF.
 *   - bump version.
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
 * Revision 1.4  2006/08/11 07:05:36  platin
 *
 *   - bump to 0.5p1
 *   - not using interaction picture anymore; add ad-hoc correction to
 *     the positivity problem.
 *
 * Revision 1.3  2006/07/31 08:36:43  platin
 *
 *   - align laser pulse phase as 0 at center of pulse.
 *   - collect total signal, instead of t>t3, to reflect expr. condition.
 *   - use Crank-Nicholson split-operator short-time propagation scheme.
 *
 * Revision 1.2  2006/07/20 18:41:46  platin
 *
 *   - sign problem should be fixed.
 *   - bump version number for tnl-dm3pes to 0.4.
 *
 * Revision 1.1  2006/07/10 20:25:41  platin
 *
 *   - a possible sign problem in dm3pes.c,
 *   - fist import of the tnl-dm2pes module, which calculates 3P PE signal
 *     using non-Markovian TNL method from Meier and Tannor.
 *
 *
 */
