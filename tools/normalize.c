
/***************************************************
 * normalize.c
 *
 * Simply take x,y pairs of data and normalize the set
 * so that the maximum of y is 1.0.
 *
 * By Yuan-Chung Cheng <yccheng@berkeley.edu>
 *
 *
 ***************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <ctype.h>
#include <string.h>
#include <math.h>

/* Maxinum number of points, size of cache. */
#define MAX_NPOINTS  10000

int main()
{
  char buffer[256];

  double x[MAX_NPOINTS],y[MAX_NPOINTS];
  double max;

  int i,N;
  int kmax;

  // read data straight from stdin
  N=0;
  while(1) {
    fscanf(stdin,"%s",buffer);
    if(feof(stdin)) break;
    x[N]=strtod(buffer,NULL);
    buffer[0]='\0'; /* clear buffer in case the next read fails */
    fscanf(stdin,"%s",buffer);
    y[N]=strtod(buffer,NULL);
    N++;
  }

  kmax=0;
  max=y[0];
  for(i=0;i<N;i++) {
    if(y[i]>max) {
      kmax=i;
      max=y[i];
    }
  }

  /* output normalized data */
  for(i=0;i<N;i++) {
    printf("%20.6f %20.16f\n",x[i],y[i]/max);
  }

  return 1;

}

/*
 * $Log$
 * Revision 1.2  2006/07/12 17:44:20  platin
 *   - minor change to leave space in front of numbers.
 *
 * Revision 1.1  2006/07/12 17:32:42  platin
 *
 *   - cspline_max: use cspline interpolation to estimate maximum points.
 *   - normalize: normalize data set so that maximum of y is 1.0.
 *
 *
 */
