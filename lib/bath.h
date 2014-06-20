/***************************************************
 * bath.h
 * 
 * Header file for bath.c.
 *
 * By Yuan-Chung Cheng <yuanchung@ntu.edu.tw>
 *
 *
 ***************************************************/

#ifndef _BATH_H
#define _BATH_H 1

#include "qdas.h"
#include "spectrum.h"

/* exported functions */
void bath_mod_init();
int bath_lookup_id(char *keyword);
void bath_init_params(qdas_keys *keys);
void bath_free_params(int id);
char *bath_get_keyword(int id);
char *bath_get_description(int id);

#endif /* bath.h */

/*
 * $Log$
 * Revision 1.1  2006/05/24 00:42:18  platin
 * Initial revision
 *
 *
 */
