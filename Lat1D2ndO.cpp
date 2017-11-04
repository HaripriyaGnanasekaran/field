#include "Lat1D2ndO.h"
#include "tools.h"
#include <iostream>

// static Vectors for the Propagate functions. This avoid allocating memory
// at each call to a Propagate function, which is very time intensive
Vector Lat1D2ndO::Gx0;


Lat1D2ndO::Lat1D2ndO(Input* MyInput_, Text name_) : Lat1DFlat(MyInput_,name_) {
	Gx0.Dim(1,3*GetTotalNumLayers());
}

Lat1D2ndO::~Lat1D2ndO() {
}
Boolean
Lat1D2ndO::OverflowProtection(void) const {
	return false;
}
void
Lat1D2ndO::MakeSafe(Vector) const {
}
void
Lat1D2ndO::RestoreFromSafe(Vector) const {
}

void
Lat1D2ndO::GetLatticeInfo(int *Info) const {
	Message(fatal,"No getlatticeinfo in 1d");
}

void
Lat1D2ndO::PropagateF(Matrix Gi, Vector G, const int s, const double S) const {
	int M =GetTotalNumLayers();

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

	    YisCtimesX(gx0+1,gz0,S,M-1);
	YplusisCtimesX(gx0+1,gz1,4*L,M-1);

	    YisCtimesX(gx1,gz0,L,M);
	YplusisCtimesX(gx1,gz1,2*L+S,M);
	YplusisCtimesX(gx1,gz2,L,M);

	    YisCtimesX(gx2,gz1+1,4*L,M-1);
	YplusisCtimesX(gx2,gz2+1,S,M-1);

	for (int k=0; k<3; k++) times(gs0+k*M,gx0+k*M,g,M);
}

void
Lat1D2ndO::PropagateF(Matrix Gi, Vector G, const int s, const double S, const bool stiff_range) const {
Message(fatal,"PropagateF with stiff-range not implemented in 2d");
}

void
Lat1D2ndO::PropagateF(Matrix Gi, Vector G, const int s) const {
Message(fatal,"PropagateF without S not implemented in 2d");
}

void
Lat1D2ndO::PropagateF(Vector Gi, Vector G, const double S) const {
	int M=GetTotalNumLayers();

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


	    YisCtimesX(gx0+1,gz0,S,M-1);
	YplusisCtimesX(gx0+1,gz1,4*L,M-1);

	    YisCtimesX(gx1,gz0,L,M);
	YplusisCtimesX(gx1,gz1,2*L+S,M);
	YplusisCtimesX(gx1,gz2,L,M);

	    YisCtimesX(gx2,gz1+1,4*L,M-1);
	YplusisCtimesX(gx2,gz2+1,S,M-1);

	for (int k=0; k<3; k++) times(gs0+k*M,gx0+k*M,g,M);
}

void
Lat1D2ndO::PropagateF(Vector Gi, Vector G, const double S, const bool stiffness) const {
	Message(fatal,"PropagateF vector with stiff-range not implemented in 2d");
}

void
Lat1D2ndO::PropagateF(Vector Gi, Vector G) const {
Message(fatal,"PropagateF without S not implemented in 2d");
}
void
Lat1D2ndO::PropagateB(Vector Gi, Vector G,const double S) const {
	int M=GetTotalNumLayers();

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

	    YisCtimesX(gx0,gz0+1,S,M-1);
	YplusisCtimesX(gx0,gz1,4*L,M);

	    YisCtimesX(gx1,gz0+1,L,M-1);
	YplusisCtimesX(gx1,gz1,2*L+S,M);
	YplusisCtimesX(gx1+1,gz2,L,M);

	    YisCtimesX(gx2+1,gz2,S,M-1);
	YplusisCtimesX(gx2,gz1,4*L,M);

	for (int k=0; k<3; k++) times(gs0+k*M,gx0+k*M,g,M);
}
void
Lat1D2ndO::PropagateB(Vector Gi, Vector G,const double S, const bool stiffness) const {
Message(fatal,"PropagateB vector with stiff-range not implemented in 2d");
}
void
Lat1D2ndO::PropagateB(Vector Gi, Vector G) const {
Message(fatal,"PropagateB vector without S not implemented in 2d");
}

void
Lat1D2ndO::Init2G(Vector Gi1, Vector Gi2, const Vector G, const LatticeRange* LatRange) const {
	Message(fatal,"Init2G not implemented in 3d");
}
void
Lat1D2ndO::Propagate2G(Matrix Gi1, Matrix Gi2, const Vector G, const int s,const LatticeRange* LatRange) const {
Message(fatal,"Propagate 2G not implemented in 3d");
}
void
Lat1D2ndO::Propagate2G(Matrix Gi1, Matrix Gi2, const Vector G, const int s,const LatticeRange* LatRange,const double f) const {
Message(fatal,"Propagate 2G not implemented in 3d");
}

void
Lat1D2ndO::Propagate2G(Vector Gi1, Vector Gi2, const Vector G, const LatticeRange* LatRange) const {
Message(fatal,"Propagate 2G not implemented in 3d");
}
void
Lat1D2ndO::Propagate2G(Vector Gi1, Vector Gi2, const Vector G, const LatticeRange* LatRange,const double f) const {
Message(fatal,"Propagate 2G not implemented in 3d");
}

void
Lat1D2ndO::PropagateG(Matrix Gi1, const Vector G1, const int s) const {}

void
Lat1D2ndO::PropagateG(Vector G1, const Vector G2) const {}

void
Lat1D2ndO::PropagateG(const Vector Gi, const Vector G, Vector Gout) const {}

void
Lat1D2ndO::PropagateG(Matrix Gi1, const Vector G1, const int s,const double f) const {}

void
Lat1D2ndO::PropagateG(Vector G1, const Vector G2,const double f) const {}

void
Lat1D2ndO::PropagateG(const Vector Gi, const Vector G, Vector Gout, const double f) const {}


void
Lat1D2ndO::ConnectG(Vector GiA, Matrix GiB, int s, Vector out) const {
	int M=GetTotalNumLayers();
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
Lat1D2ndO::ConnectG(Vector GiA, Matrix GiB, int s, Vector out, double *PHI) const {
	int M=GetTotalNumLayers();
	double *p_out = &out[1];
	double *p_giB = &GiB[1][s];
	double *p_giA = &GiA[1];
	for (int k=0; k<3; k++) {
		removeboundaries(p_giB+k*M);
		removeboundaries(p_giA+k*M);
	}
	addTimes(p_out,p_giA,p_giB,M);         addTimes(PHI,p_giA,p_giB,M);
	addTimesC(p_out,p_giA+M,p_giB+M,4.0,M);addTimesC(PHI+M,p_giA+M,p_giB+M,4.0,M);
	addTimes(p_out,p_giA+2*M,p_giB+2*M,M); addTimes(PHI,p_giA+2*M,p_giB+2*M,M);
}

void
Lat1D2ndO::ConnectG(Vector GiA, Vector GiB, Vector out) const {
	int M=GetTotalNumLayers();
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
Lat1D2ndO::ConnectG(Vector GiA, Vector GiB, Vector out, double* PHI) const {
	int M=GetTotalNumLayers();
	double *p_out = &out[1];
	double *p_giB = &GiB[1];
	double *p_giA = &GiA[1];
	for (int k=0; k<3; k++) {
		removeboundaries(p_giB+k*M);
		removeboundaries(p_giA+k*M);
	}
	addTimes(p_out,p_giA,p_giB,M);          addTimes(PHI,p_giA,p_giB,M);
	addTimesC(p_out,p_giA+M,p_giB+M,4.0,M); addTimesC(PHI+M,p_giA+M,p_giB+M,4.0,M);
	addTimes(p_out,p_giA+2*M,p_giB+2*M,M);  addTimes(PHI,p_giA+2*M,p_giB+2*M,M);
}


void
Lat1D2ndO::Connect2G(const Vector GiA1, const Matrix GiB1, const int s1, const Vector GiA2,
						 const Matrix GiB2, const int s2, Vector out) const {
	Message(fatal,"Connect 2G not implemented in 3d");
}

void
Lat1D2ndO::Connect2G(const Vector GiA1, const Vector GiB1, const Vector GiA2,  const Vector GiB2, Vector out) const {
	Message(fatal,"Connect 2G not implemented in 3d");
}

void
Lat1D2ndO::CorrectDoubleCountG(Vector in, const Vector G) const {

	double *p_in = &in[1];
	const double *p_G = &G[1];
    //Message(debug,"CorrectDoubleCountG(Vector in, const Vector G)");
	div(p_in,p_G,GetTotalNumLayers());
}

void
Lat1D2ndO::CorrectDoubleCountG(Vector in, double *IN, const Vector G) const {
	int M=GetTotalNumLayers();
	double *p_in = &in[1];
	const double *p_G = &G[1];
    //Message(debug,"CorrectDoubleCountG(Vector in, const Vector G)");
	div(p_in,p_G,M);
	div(IN,p_G,M);div(IN+M,p_G,M);
}

double
Lat1D2ndO::ComputeLnGN(Vector Gi) const {
	double value = 0;
	int M=GetTotalNumLayers();
	//SubtractBoundaries(Gi);
 //SubtractBoundaries(pGi);
	double *p_gi = &Gi[1];
	for (int k=0; k<3; k++) removeboundaries(p_gi+k*M);
	for (int i=0; i<M; i++) value +=p_gi[i]+4.0*p_gi[M+i]+p_gi[2*M+i];
  	return log(value/6.0);
}
void
Lat1D2ndO::NormPhiFree(Vector phi, const double C) const {
	//Message(debug,"NormPhiFree(Vector phi, const double C)");
	int M=GetTotalNumLayers();
	norm(phi,C,M);
}

void
Lat1D2ndO::NormPhiFree(Vector phi, double *PHI, const double C) const {
	int M=GetTotalNumLayers();
	//Message(debug,"NormPhiFree(Vector phi, const double C)");
	norm(phi,C,M); norm(PHI,C,M); norm(PHI+M,C,M);
}

void
Lat1D2ndO::NormPhiRestr(Vector phi, const Vector Gi, double C) const {
	int M=GetTotalNumLayers();
	C /= 6.0*exp(ComputeLnGN(Gi));
	norm(phi,C,M);
}
void
Lat1D2ndO::NormPhiRestr(Vector phi, double *PHI, const Vector Gi, double C) const {
	int M=GetTotalNumLayers();
	C /= 6.0*exp(ComputeLnGN(Gi));
	norm(phi,C,M); norm(PHI,C,M); norm(PHI+M,C,M);
}
void
Lat1D2ndO::UpdateBoundaries(Vector A) const {
	//Message(debug,"UpdateBoundaries(Vector A)");
	SetBoundaries(A);
}
void
Lat1D2ndO::UpdateBoundaries(Matrix A, const int s) const {
	Message(fatal,"report error: UpdateBoundaries not implemented in lat1d2nd0");
}
