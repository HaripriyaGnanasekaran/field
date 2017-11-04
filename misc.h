#ifndef MISCxH
#define MISCxH

#include <cmath>
#include <fenk/array.h>

#ifdef _WIN32
//#using <mscorlib.dll>
//#define M_PI System::Math::PI
#define M_PI       3.14159265358979323846

double
log1p(double);
#endif

///
int
NumExtrema(const Vector,int,int);

///
Vector
PlaceExtrema(const Vector,int,int,int);

///
Vector
ValueExtrema(const Vector,int,int,int);

#endif
