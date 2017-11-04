#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <f2c.h>
#include <clapack.h>

//nvcc Open.cu -lm -lcuda -lcudart -llapack -lblas -lf2c -lcublas -arch=sm_20 -o box
//I.V. Ionova, E.A. Carter, "Error vector choice in direct inversion in the iterative subspace method, J. Compt. Chem. 17, 1836-1847, 1996.

int N,N1,N2,n_seg,n_box;
flat CHI,n,GN;
int i,j,k,k_diis,m,s,Mx,My,Mz,M,jx,jy,it,bx1,by1,bz1,bxm,bym,bzm,iv;
float  error = 1, tolerance = 1e-7, eta = 0.1;
float *Aij,*Ci,*Apij,*H_phi,*H_u,*phi,*phi,*phitot,*G1,*alpha,*Gg_f,*Gg_b,*phi_side,*x,*x0,*g,*xR,*x_x0,*mask1,*mask2;
int *px1,*py1,*pz1,*px2,*py2,*pz2;
float *u; //alleen een pointer, voor ons gemak.
boolean SCF_iteration; 

