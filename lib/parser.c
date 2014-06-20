/***************************************************
 * parser.c
 * 
 * By Yuan-Chung Cheng <yccheng@mit.edu>
 *
 * General input file parser...
 *
 ***************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <ctype.h>
#include "parser.h"

#define TRUE 1
#define FALSE 0

char * activated_key=NULL;
parser_klist * activated_klist=NULL;

parser_klist* parser_klist_alloc(int n)
{
  parser_klist* w;

  if (!(w=(parser_klist*)malloc(sizeof(parser_klist)))) {
    fprintf(stderr,"error in parser.c (parser_klist_alloc): %s \n", 
	    strerror(errno));
    exit(errno);
  }

  w->size=0;
  strncpy(w->delim_str,PARSER_DEFAULT_DELIM_CHARS,PARSER_STRING_BUFFER_SIZE);
  strncpy(w->rem_str,PARSER_DEFAULT_REM_CHARS,PARSER_STRING_BUFFER_SIZE);

  if (!(w->member=(parser_keyword*)malloc(n*sizeof(parser_keyword)))) {
    fprintf(stderr,"error in parser.c (parser_klist_alloc): %s \n",
	    strerror(errno));
    exit(errno);
  }
  return w;
}

void parser_klist_free(parser_klist* w)
{
  free(w->member);
  free(w);
}

int parser_keyword_add(parser_klist* w, char *str, 
		       int status, parser_callback_function func)
{
  int i;

  i= w->size;
  w->member[i].name=str;
  w->member[i].status=status;
  w->member[i].function=func;
  w->size++;
}

int parser_set_delim(parser_klist* w, char *str)
{
  strncpy(w->delim_str,str,PARSER_STRING_BUFFER_SIZE);
}

int parser_set_rem(parser_klist* w, char *str)
{
  strncpy(w->rem_str,str,PARSER_STRING_BUFFER_SIZE);
}

int parser_default_func(char * str)
{
  printf("ERROR: Unknown keyword \"%s\" !\n",str);
  return 0;
}

parser_default_function parser_f=parser_default_func;

int parser_set_default_func(parser_default_function func)
{
  parser_f=func;
}

int isinstr(char c,char *str)
{
  char *p;

  p=str;

  while(*p != '\0')
    if( c == *p++)
      return TRUE;

  return FALSE;
}


char *parser_strtok(char *s,char *delim)
{
  static char *pstr=NULL;
  static char *p=NULL;

  if (s != NULL) {
    pstr=s;
  } else {
    pstr=p;
  }

  /* skip a delim sections to the beginning of the next token */
  while(isinstr(*pstr,delim))
    pstr++;

  if (*pstr == '\0')
    return NULL;

  /* put the ending NULL char */
  p=pstr;
  while(!isinstr(*p,delim))
    p++;

  *p='\0';

  p++;

#ifdef DEBUG_PARSER
  printf("STRTOK: %s\n",pstr);
#endif

  return (char *) pstr;
}

int prefix_strcasecmp(const char *a, const char *b)
{
  int i;
  i=0;
  
  while(a[i] != '\0' && b[i] != '\0') {
    if(toupper(a[i]) != toupper(b[i]))
      return 1;
    i++;
  }

  return 0;
}

int parser_getline(FILE *ifile, size_t buffersize, char *buffer, 
		   char *delim, char *rem)
{
  int i,j;
  int len;
  int delim_state;
  
  buffer[0]='\\';
  buffer[1]='\0';
  len=1;

  while(len < buffersize && buffer[len-1] == '\\') {

    fgets(buffer+len-1,buffersize-len,ifile);
    
    /* rephase the string */

    /* strip off comment chars */
    i=0;
    while(buffer[i] != '\0' && !isinstr(buffer[i],rem) )
      i++;
    buffer[i]='\0';
    
    len=strlen(buffer);

    if(buffer[len-1] == '\n') {
      buffer[len-1]='\0';
      len--;
    }

    if(len == 0) { /* remaining part is an empty string */
      buffer[0]='\\';
      buffer[1]='\0';
      len=1;
    }

    if (feof(ifile)) {
      if(len == 1 && buffer[0]=='\\') {
	len=0;
	buffer[0]='\0';
      }
      return len;
    }
  }

  return len;
}

char *parser_next_token()
{
  return parser_strtok(NULL,activated_klist->delim_str);
}

char *parser_next_line()
{
  char *pstr;

  if (!(pstr=parser_strtok(NULL,"\n"))) {
    fprintf(stderr,"Error while parsing keyword \"%s\"!\n",activated_key);
    exit(EXIT_FAILURE);
  }
  
  return pstr;
}

double parser_next_double()
{
  char *pstr;

  if (!(pstr=parser_strtok(NULL,activated_klist->delim_str))) {
    fprintf(stderr,"Error while parsing keyword \"%s\"!\n",activated_key);
    exit(EXIT_FAILURE);
  }
  
  return atof(pstr);  
}

int parser_next_darray(double *vec,int length)
{
  char *pstr;
  int i;

  for (i=0;i<length;i++) {
    if (!(pstr=parser_strtok(NULL,activated_klist->delim_str))) {
      fprintf(stderr,"Error while parsing keyword \"%s\"!\n",activated_key);
      exit(EXIT_FAILURE);
    }
    vec[i]=atof(pstr);
  }
  return length;
}

int parser_next_int()
{
  char *pstr;

  if (!(pstr=parser_strtok(NULL,activated_klist->delim_str))) {
    fprintf(stderr,"Error while parsing keyword \"%s\"!\n",activated_key);
    exit(EXIT_FAILURE);
  }
  
  return atoi(pstr);  
}

long parser_next_long()
{
  char *pstr;

  if (!(pstr=parser_strtok(NULL,activated_klist->delim_str))) {
    fprintf(stderr,"Error while parsing keyword \"%s\"!\n",activated_key);
    exit(EXIT_FAILURE);
  }
  
  return atol(pstr);  
}

int parser_dispatch(char *key, parser_klist *w)
{
  int wsize;
  int i;
  
  wsize=w->size;
  for(i=0;i<wsize;i++) {

#ifdef DEBUG_PARSER
    printf("key = \"%s\", name = \"%s\"\n",key,w->member[i].name);
#endif

    if(!prefix_strcasecmp(w->member[i].name,key)) {
      if(w->member[i].status==PARSER_DONE) {
	printf("ERROR: Mutiple definition of keyword \"%s\" !\n",key);
	return 0;
      }
      activated_key=w->member[i].name; /* set the activated keyword */
      w->member[i].function(); /* run the call back function */
      if(w->member[i].status != PARSER_MULTIPLE)
	w->member[i].status=PARSER_DONE;
      return 1; /* sucess */
    }
  }
  
  return parser_f(key);
}

int parser_parse_input(FILE *ifile,parser_klist *w)
{
  char cbuffer[PARSER_FILE_BUFFER_SIZE];
  char *key=NULL;
  int wsize;
  int retvalue=1;
  int i,j;
  
  activated_klist=w; /* set the activated keylist */
  wsize=w->size;

  /* put everything into cbuffer */
  j=0;
  while(i=parser_getline(ifile, PARSER_FILE_BUFFER_SIZE-j, cbuffer+j, 
			 activated_klist->delim_str, activated_klist->rem_str)) {
    j+=i;
    cbuffer[j]='\n'; /* space with a newline */
    j++;
    cbuffer[j]='\0';
  }

#ifdef DEBUG_PARSER
  printf("READ:\"\n%s\n\"\n",cbuffer);
#endif

  /* init the token pool */
  key=parser_strtok(cbuffer,activated_klist->delim_str);
  if( !parser_dispatch(key,w))
    return 0;

  /* process all tokens */
  while(key=parser_next_token()) {

#ifdef DEBUG_PARSER
    printf("Processing \"%s\"\n",key);
#endif

    if( !parser_dispatch(key,w))
      return 0;
  }
  
  /* check if we got al required parameters */
  for(i=0;i<wsize;i++){
    if(w->member[i].status == PARSER_REQUIRED) {
      printf("ERROR: REQUIRED keyword \"%s\" not found!\n",w->member[i].name);
      retvalue=0;
    }
  }

  return retvalue;
}

#ifdef TEST

int f_cell()
{
  double tok;

  tok=parser_next_double();
  printf("cell = \"%lf\"\n",tok);
}

int f_unit()
{
  char *tok;

  tok=parser_next_token();
  printf("unit = \"%s\"\n",tok);
}

int f_par()
{
  char *tok;

  tok=parser_next_token();
  printf("par = \"%s\"\n",tok);
}

int f_title()
{
  char *tok;

  tok=parser_next_line();
  printf("title = \"%s\"\n",tok);
}

main ()
{
  int num,i;
  parser_klist *list;
  
  list=parser_klist_alloc(300);

  parser_keyword_add(list,"Cell",PARSER_REQUIRED,f_cell);
  parser_keyword_add(list,"Unit",PARSER_OPTIONAL,f_unit);
  parser_keyword_add(list,"par",PARSER_OPTIONAL,f_par);
  parser_keyword_add(list,"title",PARSER_OPTIONAL,f_title);

  if(!parser_parse_input(stdin,list))
    printf("Input file parsing error!\n");

  parser_klist_free(list);
}

#endif

/*
 * $Log$
 * Revision 1.1  2006/10/30 06:31:31  platin
 *   - add parser library to avoid external dependence.
 *
 * Revision 1.1.1.1  2001/09/08 04:23:49  platin
 *
 * 1st time to collect & put useful functions into CVS, 
 * including:
 *
 * setup.sh    ->  Shell script to help user (me) to setup the librarys
 *
 * mathconst.h ->  Header file providing useful math constants (UMC_*)
 * physconst.h ->  Header file providing useful physical constants (FPC_*)
 * unitconv.h  ->  Header file providing unit conversion factors (UCF_*)
 *
 * libzind     ->  Library used to parser ZINDO output files.
 * parser      ->  Platin's input-file parser.
 * prot        ->  planar rotation functions
 * quad        ->  wrappers for doing numerical integration by using
 *                  quadrature methods in GSL.
 * trapezoid   ->  functions to do numerical integration by using
 *                  trapezoid method.
 * vector      ->  vector related functions,, simply simulate vector by 
 *                  array of double...
 *
 *
 * Revision 1.1.1.1  2001/07/23 20:29:31  platin
 *
 * Program to calculate Dimer NonBonded Interaction with
 * Force Field methods.
 *
 * First workable version, with WS1,WS2,SJ96 force fields.
 *
 *
 * Revision 1.2  2001/06/20 06:34:57  platin
 *
 * minor bug fixed.
 *
 * Revision 1.1  2001/06/20 00:07:19  platin
 *
 * finished the parser functions,, ready to rewrite the input file parser
 * for bandm. :)
 *
 *
 */
