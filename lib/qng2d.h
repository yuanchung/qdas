/***************************************************
 * qng2d.h
 * 
 * By Yuan-Chung Cheng
 * 
 * header file for qng2d package, wrapper for 2D qng integration.
 *
 ***************************************************/

#ifndef _QNG2D_H
#define _QNG2D_H 1

#include <stdio.h>
#include <math.h>
#include <errno.h>
#include <gsl/gsl_integration.h>

/* 
   typedef 
*/

/* 2D numerical function */
typedef double (*DoubleFunc2d)(double,double);

/* function prototypes */
extern size_t qng2d(DoubleFunc2d func, double xlower, double xupper,
		    double ylower, double yupper,
		    double abseps,double releps, double *result, double *abserr);

#endif /* qng2d.h */

/*
 * $Log$
 * Revision 1.1  2007/01/03 07:30:00  platin
 *   - add a new package in which pulse duration integral is carried out
 *     to include pulse overlap effect using impulsive response
 *     function formalism.
 *
 * 
 */
