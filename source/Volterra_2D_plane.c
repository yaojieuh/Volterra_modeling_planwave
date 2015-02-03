//!########################################################################
//!                                                                       #
//! Copyright (C) University .                     #
//! This file is a modeling for 2d acoustic media  #
//!                                                                       #
//! 
//! Feb. 2014                                                      # 
//!                                                                       #
//!########################################################################
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include "init.h"
#include "P0an.h"
//#include "Reflectiondata.h"

int main( int argc, char *argv[] )
{
  // grids dimensions 
  int nz    = 101; 
  int nx    = 101;
  int ndims = nx*nz;       // total grid points number



  // physical grid size
  double dz = 5.0;     // modeling grid spacial step in z
  double dx = 5.0;    // modeling grid spacial step in x

 
  FILE* file1;
  char fname1[100];
  int ret;
  sprintf(fname1,"test_2D_Volterra_acoustic.out");
  file1 = fopen(fname1,"w");


  fprintf(file1, "--!\t                                     \n");
  fprintf(file1, "--!\t2D acoustic Volterra modeling research code\n");
  fprintf(file1, "--!\t                                     \n");

  ret = fflush(file1);

  // velocity field
  double c0=1500; //reference velocity
  double* vel  = malloc( ndims*sizeof(double) );
  double* vpe  = malloc( ndims*sizeof(double) );
  //init_acou_homog( file1, dim_w, vel);
  init_acou_layer(file1,  dx,dz,nx,nz, c0, vel,vpe);
                                                                                                                                                                            

 
   // time and frequency discretization parameter
  int   nt = 1000;        // number of steps
  double dt = 0.001;      // time step
  int nw=51;
  double dw=2;
  double* sourcet = malloc(nt*sizeof(double) );
  double* sourcefren  = malloc( (2*nw)*sizeof(double) );
  double* fren  = malloc( nw*sizeof(double) );
  init_source_ricker_fwps( file1, nt, dt,nw,dw, sourcet,fren, sourcefren); 

 

   int nz2    = 2*nz+1; 
   int nx2    = 2*nx+1;
   double* G0r = malloc(nz2*nx2*nw*sizeof(double) );
   double* G0i = malloc(nz2*nx2*nw*sizeof(double) );

   //note: Volterra Green's function is real, however, we treat it as complex here.
   double* VG0r = malloc(nz2*nx2*nw*sizeof(double) );
   double* VG0i = malloc(nz2*nx2*nw*sizeof(double) );
   greenfunctiontable(file1,nx2, nz2, dx,dz, nw,  fren, c0,G0r, G0i,VG0r,VG0i);

   double pi=3.1415926;
   int nks=nx,nks2=(nks-1)/2;
   double dks=2*pi/(nx-1)/dx;
   double ks;
   int ik;
   double* P0r = malloc(nz*nx*nw*sizeof(double) );
   double* P0i = malloc(nz*nx*nw*sizeof(double) );
   double* P1r = malloc(nz*nx*nw*sizeof(double) );
   double* P1i = malloc(nz*nx*nw*sizeof(double) );
   for(ik=nks2;ik<nks2+1;ik++){
	ks=(ik-nks2)*dks;
   	
	P0num(file1, nw, nx, nz,dx,dz, c0, ks, fren,  P0r,P0i);     
        P1num(file1,  nw, nx,  nz,dx, dz,c0,  sourcefren, fren, VG0r, VG0i, vpe, P0r,P0i, P1r,P1i);
   
   }

   fprintf(file1, "--!\t        Done                             \n");
   
}
