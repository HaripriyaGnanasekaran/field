#include "Lat1DCyl2ndO.h"
#include "tools.h"
#include <iostream>

// static Vectors for the Propagate functions. This avoid allocating memory
// at each call to a Propagate function, which is very time intensive
Vector Lat1DCyl2ndO::Gx0;
Vector Lat1DCyl2ndO::Hv;


Lat1DCyl2ndO::Lat1DCyl2ndO(Input* MyInput_, Text name_) :
Lat1DFlat(MyInput_,name_), Lat1DCylinder(MyInput_,name_)  {

	M=GetTotalNumLayers();
	Gx0.Dim(1,3*M);
	Hv.Dim(1,M);
	H=&Hv[1];
	l1 = new double[M];
	l_1 = new double[M];
	l11 = new double[M];
	l_11 = new double[M];
	for (int i=0; i<M; i++) {
		l1[i]=6*lambda1[i+1];  l11[i]=1.0-l1[i];
		l_1[i]=6*lambda_1[i+1]; l_11[i]=1.0-l_1[i];
	}

}

Lat1DCyl2ndO::~Lat1DCyl2ndO() {
	delete [] l11;
	delete [] l_11;
	delete [] l1;
	delete [] l_1;
}
Boolean
Lat1DCyl2ndO::OverflowProtection(void) const {
	return false;
}
void
Lat1DCyl2ndO::MakeSafe(Vector) const {
}
void
Lat1DCyl2ndO::RestoreFromSafe(Vector) const {
}

void
Lat1DCyl2ndO::GetLatticeInfo(int *Info) const {
	Message(fatal,"No getlatticeinfo in 1d");
}

void
Lat1DCyl2ndO::LReflect(double *Pout, double *Pin,int pos) const {
	times(Pout,l_1+1,Pin,M-1);
	addTimes(Pout,l_11+1,Pin+(1-pos)*2*M+1,M-1);
}

void
Lat1DCyl2ndO::UReflect(double *Pout, double *Pin,int pos) const {
	times(Pout+1,l1,Pin+1,M-1);
	addTimes(Pout+1,l11,Pin+(1-pos)*2*M,M-1);
}

void
Lat1DCyl2ndO::PropagateF(Matrix Gi, Vector G, const int s, const double S) const {
	//int M =GetTotalNumLayers();

	double L=(1.0-S)/4.0;
	double *gs0 = &Gi[1][s];

	double *gz0 = &Gi[1][s-1];
	double *gz1 = &Gi[1+M][s-1];
	double *gz2 = &Gi[1+2*M][s-1];

	SetBx1(gz0,gz2);SetBx1(gz1);
	SetBxm(gz0,gz2);SetBxm(gz1);

	double *gx0 = &Gx0[1];
	double *gx1 = &Gx0[1+M];
	double *gx2 = &Gx0[1+2*M];

	double *g = &G[1];

	LReflect(H,gz0,0); 	YisCtimesX(gx0+1,H,S,M-1);
	LReflect(H,gz1,1); 	YplusisCtimesX(gx0+1,H,4*L,M-1);

	    				YisCtimesX(gx1,gz0,L,M);
						YplusisCtimesX(gx1,gz1,2*L+S,M);
						YplusisCtimesX(gx1,gz2,L,M);

	UReflect(H,gz1,1); 	YisCtimesX(gx2,H+1,4*L,M-1);
	UReflect(H,gz2,2);  YplusisCtimesX(gx2,H+1,S,M-1);

	for (int k=0; k<3; k++) times(gs0+k*M,gx0+k*M,g,M);
}

void
Lat1DCyl2ndO::PropagateF(Matrix Gi, Vector G, const int s, const double S, const bool stiff_range) const {
Message(fatal,"PropagateF with stiff-range not implemented in 2d");
}

void
Lat1DCyl2ndO::PropagateF(Matrix Gi, Vector G, const int s) const {
Message(fatal,"PropagateF without S not implemented in 2d");
}

void
Lat1DCyl2ndO::PropagateF(Vector Gi, Vector G, const double S) const {
	//int M=GetTotalNumLayers();

	double L=(1.0-S)/4.0;
	double *gs0 = &Gi[1];

	double *gz0 = &Gi[1];
	double *gz1 = &Gi[1+M];
	double *gz2 = &Gi[1+2*M];

	SetBx1(gz0,gz2);SetBx1(gz1);
	SetBxm(gz0,gz2);SetBxm(gz1);

	double *gx0 = &Gx0[1];
	double *gx1 = gx0+M;
	double *gx2 = gx0+2*M;
	double *g = &G[1];


	LReflect(H,gz0,0); 	YisCtimesX(gx0+1,H,S,M-1);
	LReflect(H,gz1,1); 	YplusisCtimesX(gx0+1,H,4*L,M-1);

	   					YisCtimesX(gx1,gz0,L,M);
						YplusisCtimesX(gx1,gz1,2*L+S,M);
						YplusisCtimesX(gx1,gz2,L,M);

	UReflect(H,gz1,1); 	YisCtimesX(gx2,H+1,4*L,M-1);
	UReflect(H,gz2,2); 	YplusisCtimesX(gx2,H+1,S,M-1);

	for (int k=0; k<3; k++) times(gs0+k*M,gx0+k*M,g,M);
}

void
Lat1DCyl2ndO::PropagateF(Vector Gi, Vector G, const double S, const bool stiffness) const {
	Message(fatal,"PropagateF vector with stiff-range not implemented in 2d");
}

void
Lat1DCyl2ndO::PropagateF(Vector Gi, Vector G) const {
Message(fatal,"PropagateF without S not implemented in 2d");
}
void
Lat1DCyl2ndO::PropagateB(Vector Gi, Vector G,const double S) const {
	//int M=GetTotalNumLayers();

	double L=(1.0-S)/4.0;
	double *gs0 = &Gi[1];

	double *gz0 = &Gi[1];
	double *gz1 = &Gi[1+M];
	double *gz2 = &Gi[1+2*M];

	SetBx1(gz0,gz2);SetBx1(gz1);
	SetBxm(gz0,gz2);SetBxm(gz1);

	double *gx0 = &Gx0[1];
	double *gx1 = gx0+M;
	double *gx2 = gx0+2*M;

	double *g = &G[1];

	UReflect(H,gz0,0);
	    YisCtimesX(gx0,H+1,S,M-1);
		YplusisCtimesX(gx0,gz1,4*L,M);

	    YisCtimesX(gx1,H+1,L,M-1);
		YplusisCtimesX(gx1,gz1,2*L+S,M);
	LReflect(H,gz2,2);
		YplusisCtimesX(gx1+1,H,L,M);

	    YisCtimesX(gx2+1,H,S,M-1);
		YplusisCtimesX(gx2,gz1,4*L,M);

	for (int k=0; k<3; k++) times(gs0+k*M,gx0+k*M,g,M);
}
void
Lat1DCyl2ndO::PropagateB(Vector Gi, Vector G,const double S, const bool stiffness) const {
Message(fatal,"PropagateB vector with stiff-range not implemented in 2d");
}
void
Lat1DCyl2ndO::PropagateB(Vector Gi, Vector G) const {
Message(fatal,"PropagateB vector without S not implemented in 2d");
}

void
Lat1DCyl2ndO::Init2G(Vector Gi1, Vector Gi2, const Vector G, const LatticeRange* LatRange) const {
	Message(fatal,"Init2G not implemented in 3d");
}
void
Lat1DCyl2ndO::Propagate2G(Matrix Gi1, Matrix Gi2, const Vector G, const int s,const LatticeRange* LatRange) const {
Message(fatal,"Propagate 2G not implemented in 3d");
}
void
Lat1DCyl2ndO::Propagate2G(Matrix Gi1, Matrix Gi2, const Vector G, const int s,const LatticeRange* LatRange,const double f) const {
Message(fatal,"Propagate 2G not implemented in 3d");
}

void
Lat1DCyl2ndO::Propagate2G(Vector Gi1, Vector Gi2, const Vector G, const LatticeRange* LatRange) const {
Message(fatal,"Propagate 2G not implemented in 3d");
}
void
Lat1DCyl2ndO::Propagate2G(Vector Gi1, Vector Gi2, const Vector G, const LatticeRange* LatRange,const double f) const {
Message(fatal,"Propagate 2G not implemented in 3d");
}

void
Lat1DCyl2ndO::PropagateG(Matrix Gi1, const Vector G1, const int s) const {}

void
Lat1DCyl2ndO::PropagateG(Vector G1, const Vector G2) const {}

void
Lat1DCyl2ndO::PropagateG(const Vector Gi, const Vector G, Vector Gout) const {}

void
Lat1DCyl2ndO::PropagateG(Matrix Gi1, const Vector G1, const int s,const double f) const {}

void
Lat1DCyl2ndO::PropagateG(Vector G1, const Vector G2,const double f) const {}

void
Lat1DCyl2ndO::PropagateG(const Vector Gi, const Vector G, Vector Gout, const double f) const {}

void
Lat1DCyl2ndO::ConnectG(Vector GiA, Matrix GiB, int s, Vector out) const {
	//int M=GetTotalNumLayers();
	double *p_out = &out[1];
	double *p_giB = &GiB[1][s];
	double *p_giA = &GiA[1];
	for (int k=0; k<3; k++) {
		removeboundaries(p_giB+k*M);
		removeboundaries(p_giA+k*M);
	}
	addTimes(p_out,p_giA,p_giB,M);
	addTimesC(p_out,p_giA+M,p_giB+M,4.0,M);
	addTimes(p_out,p_giA+2*M,p_giB+2*M,M);
}
void
Lat1DCyl2ndO::ConnectG(Vector GiA, Matrix GiB, int s, Vector out, double *OUT) const {
	//int M=GetTotalNumLayers();
	double *p_out = &out[1];
	double *p_giB = &GiB[1][s];
	double *p_giA = &GiA[1];
	for (int k=0; k<3; k++) {
		removeboundaries(p_giB+k*M);
		removeboundaries(p_giA+k*M);
	}
	addTimes(p_out,p_giA,p_giB,M);         addTimes(OUT,p_giA,p_giB,M);
	addTimesC(p_out,p_giA+M,p_giB+M,4.0,M);addTimesC(OUT+M,p_giA+M,p_giB+M,4.0,M);
	addTimes(p_out,p_giA+2*M,p_giB+2*M,M); addTimes(OUT,p_giA+2*M,p_giB+2*M,M);
}

void
Lat1DCyl2ndO::ConnectG(Vector GiA, Vector GiB, Vector out) const {
	//int M=GetTotalNumLayers();
	double *p_out = &out[1];
	double *p_giB = &GiB[1];
	double *p_giA = &GiA[1];
	for (int k=0; k<3; k++) {
		removeboundaries(p_giB+k*M);
		removeboundaries(p_giA+k*M);
	}
	addTimes(p_out,p_giA,p_giB,M);
	addTimesC(p_out,p_giA+M,p_giB+M,4.0,M);
	addTimes(p_out,p_giA+2*M,p_giB+2*M,M);
}

void
Lat1DCyl2ndO::ConnectG(Vector GiA, Vector GiB, Vector out, double *OUT) const {
	//int M=GetTotalNumLayers();
	double *p_out = &out[1];
	double *p_giB = &GiB[1];
	double *p_giA = &GiA[1];
	for (int k=0; k<3; k++) {
		removeboundaries(p_giB+k*M);
		removeboundaries(p_giA+k*M);
	}
	addTimes(p_out,p_giA,p_giB,M);         addTimes(OUT,p_giA,p_giB,M);
	addTimesC(p_out,p_giA+M,p_giB+M,4.0,M);addTimesC(OUT+M,p_giA+M,p_giB+M,4.0,M);
	addTimes(p_out,p_giA+2*M,p_giB+2*M,M); addTimes(OUT,p_giA+2*M,p_giB+2*M,M);
}


void
Lat1DCyl2ndO::Connect2G(const Vector GiA1, const Matrix GiB1, const int s1, const Vector GiA2,
						 const Matrix GiB2, const int s2, Vector out) const {
	Message(fatal,"Connect 2G not implemented in 3d");
}

void
Lat1DCyl2ndO::Connect2G(const Vector GiA1, const Vector GiB1, const Vector GiA2,  const Vector GiB2, Vector out) const {
	Message(fatal,"Connect 2G not implemented in 3d");
}

void
Lat1DCyl2ndO::CorrectDoubleCountG(Vector in, const Vector G) const {
	double *p_in = &in[1];
	const double *p_G = &G[1];
    //Message(debug,"CorrectDoubleCountG(Vector in, const Vector G)");
	div(p_in,p_G,M);
}

void
Lat1DCyl2ndO::CorrectDoubleCountG(Vector in, double *PHI, const Vector G) const {
	double *p_in = &in[1];
	const double *p_G = &G[1];
    //Message(debug,"CorrectDoubleCountG(Vector in, const Vector G)");
	div(p_in,p_G,M); div(PHI,p_G,M); div(PHI+M,p_G,M);
}

double
Lat1DCyl2ndO::ComputeLnGN(Vector Gi) const {
	int z;
	double value = (Gi[2]+4.0*Gi[M+2]+Gi[2*M+2])*L[2];
	if (bound1 == 3) value/=2;
	for (z=3; z<=M-2; z++) {
		value += (Gi[z]+4.0*Gi[M+z]+Gi[2*M+z])*L[z];
	}
	if (bound2 == M-2) value += (Gi[M-1]+4.0*Gi[M+M-1]+Gi[2*M+M-1])*L[M-1]/2;
	else if (M > 3) value += (Gi[M-1]+4.0*Gi[M+M-1]+Gi[2*M+M-1])*L[M-1];
	value += numExtraLatticeSites; // usually zero.
  	return log(value/6.0);
}
void
Lat1DCyl2ndO::NormPhiFree(Vector phi, const double C) const {
	//Message(debug,"NormPhiFree(Vector phi, const double C)");
	norm(phi,C,M);
}
void
Lat1DCyl2ndO::NormPhiFree(Vector phi,double *PHI, const double C) const {
	//Message(debug,"NormPhiFree(Vector phi, const double C)");
	norm(phi,C,M); norm(PHI,C,M); norm(PHI+M,C,M);
}
void
Lat1DCyl2ndO::NormPhiRestr(Vector phi, const Vector Gi, double C) const {
	C /= 6.0*exp(ComputeLnGN(Gi));
	norm(phi,C,M);
}
void
Lat1DCyl2ndO::NormPhiRestr(Vector phi, double *PHI, const Vector Gi, double C) const {
	C /= 6.0*exp(ComputeLnGN(Gi));
	norm(phi,C,M); norm(PHI,C,M); norm(PHI+M,C,M);
}
void
Lat1DCyl2ndO::UpdateBoundaries(Vector A) const {
	//Message(debug,"UpdateBoundaries(Vector A)");
	SetBoundaries(A);
}
void
Lat1DCyl2ndO::UpdateBoundaries(Matrix A, const int s) const {
	Message(fatal,"report error: UpdateBoundaries not implemented in lat1d2nd0");
}
