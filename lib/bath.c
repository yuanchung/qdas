/***************************************************
 * bath.c
 *
 * Part of the qdas package;
 * Bath module interface for the package.
 *
 * By Yuan-Chung Cheng <yuanchung@ntu.edu.tw>
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
#include <sys/types.h>
#include <sys/stat.h>

#include "bath.h"
#include "qdas.h"
#include "spectrum.h"

/* modules to load */
#include "bath_htgau.h"
#include "bath_js03.h"
#include "bathtype_mt99art.h"
#include "bath_mt99ohm.h"
#include "bath_js03art.h"
#include "bath_ohmart.h"
#include "bath_mboart.h"
#include "bath_mt99art.h"

/* Constants that define the capacity of the program */
#define MAX_BATH_NUMBERS (16)
#define STR_BUFFER_SIZE (256)

/* data structure to hold a "bath module" */
struct bath_module_struct 
{
  int id;
  char *keyword;
  char *description;
  int (* f_init_params)(const size_t, const double, const size_t, const double *);
  int (* f_free_params)();
};

typedef struct bath_module_struct bath_module;

/* Global variables */
bath_module BathModulesList[MAX_BATH_NUMBERS]; // no more than 16 modules now.

/* initial all modules, i.e. fill in and define all modules */
void bath_mod_init()
{
  int i;

  for(i=0;i<MAX_BATH_NUMBERS;i++) {
    BathModulesList[i].id = -1; // id <0 indicates the end of the list
    BathModulesList[i].keyword = NULL;
    BathModulesList[i].description = NULL;
    BathModulesList[i].f_init_params = NULL;
    BathModulesList[i].f_free_params = NULL;
  }

  BathModulesList[0].id = 0; // id should match index number
  BathModulesList[0].keyword = strdup("HTGAU");
  BathModulesList[0].description = strdup("High-temperature Gaussian bath model of JSF.");
  BathModulesList[0].f_init_params = bath_htgau_init_params;
  BathModulesList[0].f_free_params = bath_htgau_free_params;

  BathModulesList[1].id = 1; // id should match index number
  BathModulesList[1].keyword = strdup("JS03");
  BathModulesList[1].description = strdup("Jang and Silbey's ILT and BChl J(w).");
  BathModulesList[1].f_init_params = bath_js03_init_params;
  BathModulesList[1].f_free_params = bath_js03_free_params;

  BathModulesList[2].id = 2; // id should match index number
  BathModulesList[2].keyword = strdup("MT99OHM");
  BathModulesList[2].description = strdup("Meier and Tannor's parameterization of Ohmic J(w).");
  BathModulesList[2].f_init_params = bath_mt99ohm_init_params;
  BathModulesList[2].f_free_params = bath_mt99art_free_params;

  BathModulesList[3].id = 3; // id should match index number
  BathModulesList[3].keyword = strdup("JS03ART");
  BathModulesList[3].description = strdup("Meier and Tannor's parameterization of JS03 BChl J(w).");
  BathModulesList[3].f_init_params = bath_js03art_init_params;
  BathModulesList[3].f_free_params = bath_mt99art_free_params;

  BathModulesList[4].id = 4; // id should match index number
  BathModulesList[4].keyword = strdup("OHMART");
  BathModulesList[4].description = strdup("Artificial parameterization of Ohmic J(w).");
  BathModulesList[4].f_init_params = bath_ohmart_init_params;
  BathModulesList[4].f_free_params = bath_mt99art_free_params;

  BathModulesList[5].id = 5; // id should match index number
  BathModulesList[5].keyword = strdup("MBOART");
  BathModulesList[5].description = strdup("Artificial Multimode Brownian Oscillator model.");
  BathModulesList[5].f_init_params = bath_mboart_init_params;
  BathModulesList[5].f_free_params = bath_mt99art_free_params;

  BathModulesList[6].id = 6; // id should match index number
  BathModulesList[6].keyword = strdup("MT99ART");
  BathModulesList[6].description = strdup("Meier and Tannor's artificial bath model.");
  BathModulesList[6].f_init_params = bath_mt99art_init_params;
  BathModulesList[6].f_free_params = bath_mt99art_free_params;

}

int bath_lookup_id(char *keyword)
{
  int i;
  int ret;

  ret=-1;
  for(i=0;i<MAX_BATH_NUMBERS && BathModulesList[i].id >=0;i++) {
    if(!strcasecmp(keyword,BathModulesList[i].keyword)) {
      ret = BathModulesList[i].id;
      break;
    }
  }

  return ret;
}

// initialize the bath module as well as printing the descriptions
void bath_init_params(qdas_keys *keys)
{
  int i,id;
  // handle multiple bath groups
  for(i=0;i<keys->bathgroup_nterms;i++) {
    id=keys->bathgroup[i].bath_mod;
    printf("\n");
    printf("%s\n",bath_get_description(id));
    printf("\n");
    if(BathModulesList[id].f_init_params != NULL) {
      BathModulesList[id].f_init_params(keys->nsize,keys->beta,
					keys->bathgroup[i].bath_nparams,keys->bathgroup[i].bath_params);
    } else {
      printf("Bath module ID out of range!\n");
      exit(EXIT_FAILURE);
    }
  }
}

void bath_free_params(int id)
{
  if(BathModulesList[id].f_free_params != NULL) {
    BathModulesList[id].f_free_params();
  } else {
    printf("Bath module ID out of range!\n");
    exit(EXIT_FAILURE);
  }
}

char *bath_get_keyword(int id) {
  if(BathModulesList[id].keyword != NULL) {
    return BathModulesList[id].keyword;
  } else {
    printf("Bath module ID out of range!\n");
    exit(EXIT_FAILURE);
  }
}

char *bath_get_description(int id) {
  if(BathModulesList[id].keyword != NULL) {
    return BathModulesList[id].description;
  } else {
    printf("Bath module ID out of range!\n");
    exit(EXIT_FAILURE);
  }
}

