#include "tools.h"
#ifdef HAS_BLAS
#include <cblas.h>
#endif

void YisCtimesX(Vector Y,Vector X, double C, int last) {double *pY=&Y[1]; double *pX=&X[1]; YisCtimesX(pY,pX,C,last); }
void YisCtimesX(double *Y,double *X, double C, int last) {
#ifdef _OPENMP
	#pragma omp parallel for
#endif
	for (int i=0; i<last; i++) Y[i]=X[i] * C;
}
void YplusisCtimesX(Vector Y,Vector X, double C, int last) {double *pY=&Y[1]; double *pX=&X[1]; YplusisCtimesX(pY,pX,C,last); }
void YplusisCtimesX(double *Y,double *X, double C, int last) {
#ifdef _OPENMP
	#pragma omp parallel for
#endif
	for (int i=0; i<last; i++) Y[i]+=X[i] * C;
}

void times(Vector P,  const Vector A,  const Vector B, const int last) {double *pP=&P[1]; const double *pA=&A[1]; const double *pB=&B[1]; times(pP,pA,pB,last);}
void times(double *P, const double *A,  const double *B, const int last) {
#ifdef _OPENMP
	#pragma omp parallel for
#endif
	for (int i=0; i<last; i++) P[i]=A[i] * B[i];}
void times(Vector P,  const Vector A,  const Vector B, const Vector C, const int last) {double *pP=&P[1]; const double *pA=&A[1]; const double *pB=&B[1]; const double *pC=&C[1]; times(pP,pA,pB,pC,last);}
void times(double *P, const double *A,  const double *B, const double *C, const int last) {
#ifdef _OPENMP
	#pragma omp parallel for
#endif
	for (int i=0; i<last; i++) P[i]=A[i] * B[i]* C[i];}
void addTimes(Vector P, const Vector A, const Vector B, const int last) {double *pP=&P[1]; const double *pA=&A[1]; const double *pB=&B[1]; addTimes(pP,pA,pB,last);}

void addTimes(double *P, const double *A, const double *B, const int last) {
#ifdef _OPENMP
    #pragma omp parallel for
#endif
	for (int i=0; i<last; i++) P[i] +=A[i] * B[i];
}

void addTimesF(Vector P, const Vector A, const Vector B, const double f, const int last) {double *pP=&P[1]; const double *pA=&A[1]; const double *pB=&B[1];  addTimesF(pP,pA,pB,f,last);}
void addTimesF(double *P, const double *A, const double *B, const double f, const int last) {
#ifdef _OPENMP
	#pragma omp parallel for
#endif
	for (int i=0; i<last; i++) P[i] +=A[i] * B[i] *f;}

void addTimesC(Vector P, const Vector A, const Vector B, const double C, const int last) {double *pP=&P[1]; const double *pA=&A[1]; const double *pB=&B[1];  addTimesC(pP,pA,pB,C,last);}
void addTimesC(double *P, const double *A, const double *B, const double C, const int last) {
#ifdef _OPENMP
	#pragma omp parallel for
#endif
	for (int i=0; i<last; i++) P[i] +=A[i] * B[i] *C;}



void addTimes(Vector P, const Vector A, const Vector B, const Vector C, const int last) {double *pP=&P[1]; const double *pA=&A[1]; const double *pB=&B[1]; const double *pC=&C[1]; addTimes(pP,pA,pB,pC,last);}
void addTimes(double *P, const double *A, const double *B, const double *C, const int last) {
#ifdef _OPENMP
	#pragma omp parallel for
#endif
	for (int i=0; i<last; i++) P[i] +=A[i] * B[i] * C[i];}


void norm(Vector P, const double C, const int last)  {double *pP=&P[1]; norm(pP,C,last);};
void norm(double *P, const double C, const int last) {
#ifndef HAS_BLAS
#ifdef _OPENMP
	#pragma omp parallel for
#endif
	for (int i=0; i<last; i++) P[i] *=C;
#endif
#ifdef HAS_BLAS
	cblas_dscal(last,C,P,1);
#endif
}

void timesC(Vector P, const Vector A, double C, const int last) {double *pP=&P[1]; const double *pA=&A[1]; timesC(pP,pA,C,last);}
void timesC(double *P, const double *A, double C, const int last) {
#ifdef _OPENMP
	#pragma omp parallel for
#endif
	for (int i=0; i<last; i++) P[i] = A[i] * C;}
void div(Vector P, const Vector A, const int last) {double *pP=&P[1]; const double *pA=&A[1]; div(pP,pA,last);}
void div(double *P, const double *A, const int last) {
#ifdef _OPENMP
	#pragma omp parallel for
#endif
	for (int i=0; i<last; i++) if (A[i] >0) {P[i] /=A[i] ;} else P[i]=0;}
void cp(Vector P, const Vector A, const int last) {double *pP=&P[1]; const double *pA=&A[1]; cp(pP,pA,last);}
void cp(double *P, const double *A, const int last) {
#ifndef HAS_BLAS
#ifdef _OPENMP
	#pragma omp parallel for
#endif
	for (int i=0; i<last; i++) P[i]=A[i];
#else
	cblas_dcopy(last,A,1,P,1);
#endif
}
void set(Vector P, const Vector A, const double C,  const int last) {double *pP=&P[1]; const double *pA=&A[1]; set(pP,pA,C,last);}
void set(double *P, const double *A, const double C, const int last) {
#ifdef _OPENMP
	#pragma omp parallel for
#endif
	for (int i=0; i<last; i++) P[i]=A[i]*C; }
void addset(Vector P, const Vector A, const double C, const int last) {double *pP=&P[1]; const double *pA=&A[1]; addset(pP,pA,C,last);}
void addset(double *P,const  double *A,const  double C,const  int last) {
#ifndef HAS_BLAS
#ifdef _OPENMP
	#pragma omp parallel for
#endif
	for (int i=0; i<last; i++) P[i]+=A[i]*C;
#else
	cblas_daxpy(last,C,A,1,P,1);
#endif
}
void add(Vector P, const Vector A, const Vector B, const int last) {double *pP=&P[1]; const double *pA=&A[1]; const double *pB=&B[1];  add(pP,pA,pB,last); }
void add(double *P, const double *A, const double *B, const int last) {
#ifdef _OPENMP
	#pragma omp parallel for
#endif
	for (int i=0; i<last; i++) P[i]=A[i]+B[i]; }
void add2(Vector P, const Vector A, const int last) {double *pP=&P[1];const  double *pA=&A[1];  add2(pP,pA,last);}
void add2(double *P, const double *A, const int last) {
#ifndef HAS_BLAS
#ifdef _OPENMP
	#pragma omp parallel for
#endif
	for (int i=0; i<last; i++) P[i]+=A[i];
#else
	cblas_daxpy(last,1.0,A,1,P,1);
#endif
}
void min2(Vector P, const Vector A, const Vector B, const int last) {double *pP=&P[1];const  double *pA=&A[1];const  double *pB=&B[1]; min2(pP,pA,pB,last);}
void min2(double *P,const  double *A, const double *B, const int last) {
#ifdef _OPENMP
	#pragma omp parallel for
#endif
	for (int i=0; i<last; i++) {P[i]=A[i]-B[i]; P[i]*=P[i];}}
void addmin2(Vector P,const  Vector A, const Vector B, const int last) {double *pP=&P[1]; const double *pA=&A[1];const  double *pB=&B[1]; addmin2(pP,pA,pB,last);}
void addmin2(double *P, const double *A,const  double *B, const int last) {double Pi;
#ifdef _OPENMP
	#pragma omp parallel for
#endif
	for (int i=0; i<last; i++) {Pi=A[i]-B[i]; P[i]+=Pi*Pi;}}
void Zero(Vector P, const int last) {double *pP=&P[1]; Zero(pP,last);}
void Zero(double *P,const int last) {
//	#pragma omp parallel for
//	for (int i=0; i<last; i++) P[i]=0;
memset(P,0,last*sizeof(*P));}

