#ifndef TOOLSxH
#define TOOLSxH
#include "fenk.h"
#ifdef _OPENMP
#include <omp.h>
#endif

void times(Vector, const  Vector, const  Vector, const int) ;
void times(double*,  const double*,  const double*, const int) ;
void YisCtimesX(Vector ,Vector , double , int ) ;
void YisCtimesX(double* ,double* , double , int ) ;
void YplusisCtimesX(Vector, Vector , double , int ) ;
void YplusisCtimesX(double* ,double* , double , int ) ;

void times(Vector, const  Vector, const  Vector, const Vector, const int) ;
void times(double*,  const double*,  const double*, const double*, const int) ;
void addTimes(Vector, const Vector, const Vector,const Vector, const int)  ;
void addTimes(double*, const double*, const double *, const double *, const int)  ;
void addTimes(Vector, const Vector, const Vector, const int)  ;
void addTimes(double*, const double*, const double *, const int)  ;

void addTimesF(Vector, const Vector, const Vector, const double, const int)  ;
void addTimesF(double*, const double*, const double *, const double,  const int)  ;
void addTimesC(Vector, const Vector, const Vector, const double, const int)  ;
void addTimesC(double*, const double*, const double *, const double,  const int)  ;


void BLAS_test(double*, double*, double*, int);
void timesC(Vector, const Vector, const double, const int)  ;
void timesC(double*, const double*, const double, const int)  ;
void norm(Vector, const double, const int) ;
void norm(double*, const double, const int) ;
void div(Vector,  const Vector, const int)  ;
void div(double*, const  double*,const int)  ;
void cp(Vector, const Vector, const int) ;
void cp(double*, const double*, const int) ;
void set(Vector, const Vector, const double,  const int) ;
void set(double*, const double*, const double,  const int) ;
void addset(Vector, const Vector, const double, const int) ;
void addset(double*, const double*, const double, const int) ;
void add(Vector, const Vector, const Vector, const int) ;
void add(double*, const double*, const double*, const int) ;
void add2(Vector, const Vector, const int) ;
void add2(double*, const double*, const int) ;
void min2(Vector, const Vector, const Vector, const int) ;
void min2(double*, const double*, const double*, const int) ;
void addmin2(Vector, const Vector, const Vector, const int) ;
void addmin2(double*, const double*, const double*, const int) ;
void Zero(Vector,const int);
void Zero(double*, const int);


#endif
