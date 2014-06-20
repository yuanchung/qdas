/***************************************************
 * params.h
 * 
 * By Yuan-Chung Cheng <yuanchung@ntu.edu.tw>
 *
 * Header file for params.c, some global constants
 * related to the parameters.
 *
 ***************************************************/

#ifndef _PARAMS_H
#define _PARAMS_H 1

#include "qdas.h"

/* shared useful functions */
/* lookup table is always the fastest !! */

/* exported functions */

/* initialize the parameters using a key file. */
void params_init(char *filename, qdas_keys *key);
void params_close(void);

#endif /* params.h */

/*
 * $Log$
 * Revision 1.1  2006/05/24 00:42:19  platin
 * Initial revision
 *
 *
 */
