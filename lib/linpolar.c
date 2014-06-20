/***************************************************
 * linpolar.c
 * 
 * By Yuan-Chung Cheng <yuanchung@ntu.edu.tw>
 *
 * Part of the qdas package;
 * Operations related to the handeling of 
 * linear polarization in a three-pulse experiment.
 *
 ***************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "linpolar.h"
#include "aux.h"

/* allocate a linpolar vector */
linpolar_vec * linpolar_vec_alloc()
{
  return (linpolar_vec *)gsl_vector_alloc(3);
}

void linpolar_vec_free(linpolar_vec *vec)
{
  gsl_vector_free((gsl_vector *)vec);
}

/* allocate a linpolar matrix */
linpolar_4vecs * linpolar_4vecs_alloc()
{
  return (linpolar_4vecs *)gsl_matrix_alloc(3,4);
}

void linpolar_4vecs_free(linpolar_4vecs *mat)
{
  gsl_matrix_free((gsl_matrix *)mat);
}

/* allocate a 3x3 rotation matrix */
linpolar_rotmat * linpolar_rotmat_alloc()
{
  return (linpolar_rotmat *)gsl_matrix_alloc(3,3);
}

void linpolar_rotmat_free(linpolar_rotmat *mat)
{
  gsl_matrix_free((gsl_matrix *)mat);
}

/* construct a linear polarization vector from a polarization angle.
   We adopt the lab frame in which the x direction is the pulse propagation
   direction and the z direction is the electric component direction. So,
   the angle given specifies a vector on the y-z plane and has a angle from
   the z axis. */
void linpolar_vec_set(linpolar_vec *vec, double angle)
{
  gsl_vector_set(vec,0,0.0); /* x is the pulse propagation direction */
  gsl_vector_set(vec,1,sin(angle));
  gsl_vector_set(vec,2,cos(angle));
}

/* this set four vectors; the ang[] should contain 4 polarization angles */
void linpolar_4vecs_set(linpolar_4vecs *mat, double angles[])
{
  int i;
  for(i=0;i<4;i++) {
    gsl_matrix_set(mat,0,i,0.0);
    gsl_matrix_set(mat,1,i,sin(angles[i]));
    gsl_matrix_set(mat,2,i,cos(angles[i]));
  }
}

/* print a vector */
void linpolar_vec_print(linpolar_vec *vec)
{
  printf("(%f, %f, %f)\n",
	 gsl_vector_get(vec,0),
	 gsl_vector_get(vec,1),
	 gsl_vector_get(vec,2));
}

void linpolar_4vecs_print(linpolar_4vecs *vecs)
{
  int i;
  for(i=0;i<4;i++) {
#ifdef NEVER
    /* this is useful for gnuplot */
    printf("0 0 0\n");
    printf("%f %f %f\n",
	   gsl_matrix_get(vecs,0,i),
	   gsl_matrix_get(vecs,1,i),
	   gsl_matrix_get(vecs,2,i));
    printf("\n");
    printf("\n");
#endif
    printf("v%d = (%f, %f, %f)\n",
	   i+1,
	   gsl_matrix_get(vecs,0,i),
	   gsl_matrix_get(vecs,1,i),
	   gsl_matrix_get(vecs,2,i));
  }
}

inline double * normalize_vec(double *vec)
{
  double vabs;
  vabs=sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]);
  vec[0]=vec[0]/vabs;
  vec[1]=vec[1]/vabs;
  vec[2]=vec[2]/vabs;
  return vec;
}

/* randomly rotate a set of 4vecs in place; this is 
   effectly a change of coordinate to a random molecular frame */
void linpolar_4vecs_ranrot(linpolar_4vecs *mat, double align, gsl_rng *r)
{
  int i;
  double a[3];
  double b[3];
  double c[3];
  double product;
  linpolar_rotmat *rotmat=linpolar_rotmat_alloc();
  linpolar_4vecs *mtmp=linpolar_4vecs_alloc();

  /* generate the molecular frame */
  /* first, generate a random vector that will be the molecular z-axis */
  gsl_ran_dir_3d(r,c,c+1,c+2);
  /* apply the alignment factor on the z-axis for orienated samples */
  c[2]=c[2]*align;
  normalize_vec(c);

  /* second, generate another random vector that will be the molecular x-axis */
  gsl_ran_dir_3d(r,a,a+1,a+2);
  /* orthogonize them */
  product=a[0]*c[0]+a[1]*c[1]+a[2]*c[2];
  a[0]=a[0]-product*c[0];
  a[1]=a[1]-product*c[1];
  a[2]=a[2]-product*c[2];
  normalize_vec(a);
  /* third axis (molecular y-axis) b = c x a */
  b[0]=c[1]*a[2]-c[2]*a[1];
  b[1]=c[2]*a[0]-c[0]*a[2];
  b[2]=c[0]*a[1]-c[1]*a[0];
  /* rotation matrix; a, b, c in each column */
  for(i=0;i<3;i++) {
    //    gsl_matrix_set(rotmat,i,0,a[i]);
    //    gsl_matrix_set(rotmat,i,1,b[i]);
    //    gsl_matrix_set(rotmat,i,2,c[i]);
    gsl_matrix_set(rotmat,0,i,a[i]);
    gsl_matrix_set(rotmat,1,i,b[i]);
    gsl_matrix_set(rotmat,2,i,c[i]);
  }
  /* rotate the 4vecs */
  gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,rotmat,mat,0.0,mtmp);
  gsl_matrix_memcpy(mat,mtmp);

  /* DONE */
  linpolar_rotmat_free(rotmat);
  linpolar_4vecs_free(mtmp);
}

/* FIXME: need to find a better way to do this */
/* Compute effective interaction matrices for four pulses;
   this corresponds to a three-pulse experiment (the fourth pulse
   corresponds to the emission/detection). Linear polarization of
   each pulse and the molecular transition dipole matrix are considered to
   calculate these interaction matrices.

   Input:
          mX, mY, mZ: x,y, and z components of the molecular transition 
          dipole matrix, a mattrix of vectors.

          polar: 3x4 matrix containing the 4-vector polarization of the pulses (column vectors).

  Output: mIp1 ... mIp4: the four interaction matrices
*/

void compute_pulse_interactions(gsl_matrix *mIp1, gsl_matrix *mIp2,gsl_matrix *mIp3,gsl_matrix *mIp4,
			        gsl_matrix *mX, gsl_matrix *mY, gsl_matrix *mZ,
				const linpolar_4vecs *polar)
{
  size_t ndim;
  int i,j;
  double E;
  gsl_matrix *mtmp;
  gsl_matrix *mI[4];
  gsl_matrix *mol[3];
  ndim=mX->size1;

  mI[0]=mIp1;
  mI[1]=mIp2;
  mI[2]=mIp3;
  mI[3]=mIp4;
  mol[0]=mX;
  mol[1]=mY;
  mol[2]=mZ;

  mtmp=gsl_matrix_alloc(ndim,ndim);
  for(i=0;i<4;i++){
    gsl_matrix_set_zero(mI[i]);
    for(j=0;j<3;j++) {
      E=gsl_matrix_get(polar,j,i);
      gsl_matrix_memcpy(mtmp,mol[j]);
      gsl_matrix_scale(mtmp,E);
      gsl_matrix_add(mI[i],mtmp);
    }
  }
  gsl_matrix_free(mtmp);
}

#ifdef MAIN
int main()
{
  linpolar_4vecs *mm;
  //  double ang[]={M_PI/3.0,-1.0*M_PI/3.0,0.0,0.0};
  double ang[]={0.0,0.0,0.0,0.0};
  gsl_rng *r=gsl_rng_alloc(gsl_rng_ranlxd2);
  gsl_rng_set (r,(unsigned)time(NULL));

  gsl_matrix *xx=gsl_matrix_alloc(3,3); // 3x3, groundstate and two excited states...
  gsl_matrix *yy=gsl_matrix_alloc(3,3);
  gsl_matrix *zz=gsl_matrix_alloc(3,3);
  gsl_matrix *p1=gsl_matrix_alloc(3,3);
  gsl_matrix *p2=gsl_matrix_alloc(3,3);
  gsl_matrix *p3=gsl_matrix_alloc(3,3);
  gsl_matrix *p4=gsl_matrix_alloc(3,3);

  int i,niter=100000;
  double aabb,abab,aaaa,bbbb;
  double MaMaMbMb,MaMbMaMb,MaMaMaMa,MbMbMbMb; /* we will calculate these orientational averages */

  /* alignment factor for anisotropic sample */
  double falign=2.0;

  /* a dimer with orth dipole moments */
  gsl_matrix_set_zero(xx);
  gsl_matrix_set_zero(yy);
  gsl_matrix_set_zero(zz);
  gsl_matrix_set(yy,0,1,1.0); // |g> -> |e1> along y
  gsl_matrix_set(zz,0,2,1.0); // |g> -> |e2> along z

  mm=linpolar_4vecs_alloc();

  MaMaMaMa=MaMaMbMb=MbMbMbMb=MaMbMaMb=0.0;
  for(i=0;i<niter;i++) {
    /* set pulse polarizations and bring them to a random molecular frame */
    linpolar_4vecs_set(mm,ang);
    //    linpolar_4vecs_print(mm);
    linpolar_4vecs_ranrot(mm,falign,r);
    //    linpolar_4vecs_print(mm);
    compute_pulse_interactions(p1,p2,p3,p4,xx,yy,zz,mm);

    /*
    printf("\n"); 
    printf("mX =\n");
    gsl_matrix_print(xx);
    printf("mY =\n");
    gsl_matrix_print(yy);
    printf("mZ =\n");
    gsl_matrix_print(zz);
    
    printf("\n");
    printf("p1 =\n");
    gsl_matrix_print(p1);
    printf("p2 =\n");
    gsl_matrix_print(p2);
    printf("p3 =\n");
    gsl_matrix_print(p3);
    printf("p4 =\n");
    gsl_matrix_print(p4);
    */
    aaaa=gsl_matrix_get(p1,0,1)*gsl_matrix_get(p2,0,1)*
      gsl_matrix_get(p3,0,1)*gsl_matrix_get(p3,0,1);
    abab=gsl_matrix_get(p1,0,1)*gsl_matrix_get(p2,0,2)*
      gsl_matrix_get(p3,0,1)*gsl_matrix_get(p3,0,2);
    aabb=gsl_matrix_get(p1,0,1)*gsl_matrix_get(p2,0,1)*
      gsl_matrix_get(p3,0,2)*gsl_matrix_get(p3,0,2);
    bbbb=gsl_matrix_get(p1,0,2)*gsl_matrix_get(p2,0,2)*
      gsl_matrix_get(p3,0,2)*gsl_matrix_get(p3,0,2);

    //    printf(" i = %d, aaaa = %f, aabb = %f, bbbb = %f\n",i,aaaa,aabb,bbbb);
    MaMaMaMa+=aaaa;
    MaMbMaMb+=abab;
    MaMaMbMb+=aabb;
    MbMbMbMb+=bbbb;
  }
   
  printf("\n");
  printf("Average (%d):\n",niter);
  printf("MaMaMaMa (%d) = %f\n",niter,MaMaMaMa/(double)niter);
  printf("MaMbMaMb (%d) = %f\n",niter,MaMbMaMb/(double)niter);
  printf("MaMaMbMb (%d) = %f\n",niter,MaMaMbMb/(double)niter);
  printf("MbMbMbMb (%d) = %f\n",niter,MbMbMbMb/(double)niter);

  linpolar_4vecs_free(mm);

  return 0;
}
#endif

/*
 * $Log$
 */
