#include "Lat2D2ndO.h"
#include "tools.h"
#include <iostream>

// static Vectors for the Propagate functions. This avoid allocating memory
// at each call to a Propagate function, which is very time intensive
Vector Lat2D2ndO::Gx0;

Lat2D2ndO::Lat2D2ndO(Input* MyInput_, Text name_) : Lat2DFlat(MyInput_,name_) {
	Gx0.Dim(1,5*M);
}

Lat2D2ndO::~Lat2D2ndO() {
}
Boolean
Lat2D2ndO::OverflowProtection(void) const {
	return false;
}
void
Lat2D2ndO::MakeSafe(Vector) const {
}
void
Lat2D2ndO::RestoreFromSafe(Vector) const {
}

void
Lat2D2ndO::GetLatticeInfo(int *Info) const {
	Message(fatal,"No getlatticeinfo in 2d");
}

void
Lat2D2ndO::PropagateF(Matrix Gi, Vector G, const int s, const double S) const {
	int jx;

	double L=(1.0-S)/4.0;
	double *gs0 = &Gi[1][s];
	double *gs1 = &Gi[1+M][s];
	double *gs2 = &Gi[1+2*M][s];
	double *gs3 = &Gi[1+3*M][s];
	double *gs4 = &Gi[1+4*M][s];

	double *gz0 = &Gi[1][s-1];
	double *gz1 = &Gi[1+M][s-1];
	double *gz2 = &Gi[1+2*M][s-1];
	double *gz3 = &Gi[1+3*M][s-1];
	double *gz4 = &Gi[1+4*M][s-1];

	SetBx1(gz0,gz4);SetBx1(gz1);SetBx1(gz2);SetBx1(gz3);
	SetBy1(gz1,gz3);SetBy1(gz0);SetBy1(gz2);SetBy1(gz4);
	SetBxm(gz0,gz4);SetBxm(gz1);SetBxm(gz2);SetBxm(gz3);
	SetBym(gz1,gz3);SetBym(gz0);SetBym(gz2);SetBym(gz4);

	double *gx0 = &Gx0[1];
	double *gx1 = &Gx0[1+M];
	double *gx2 = &Gx0[1+2*M];
	double *gx3 = &Gx0[1+3*M];;
	double *gx4 = &Gx0[1+4*M];

	double *g = &G[1];

	jx = GetNumLayers(2);

	    YisCtimesX(gx0+jx,gz0,S,M-jx);
	YplusisCtimesX(gx0+jx,gz1,L,M-jx);
	YplusisCtimesX(gx0+jx,gz2,2*L,M-jx);
	YplusisCtimesX(gx0+jx,gz3,L,M-jx);

	    YisCtimesX(gx1+1,gz0,L,M-1);
	YplusisCtimesX(gx1+1,gz1,S,M-1);
	YplusisCtimesX(gx1+1,gz2,2*L,M-1);
	YplusisCtimesX(gx1+1,gz4,L,M-1);

	    YisCtimesX(gx2,gz0,L,M);
	YplusisCtimesX(gx2,gz1,L,M);
	YplusisCtimesX(gx2,gz2,S,M);
	YplusisCtimesX(gx2,gz3,L,M);
	YplusisCtimesX(gx2,gz4,L,M);

	    YisCtimesX(gx3,gz0+1,L,M-1);
	YplusisCtimesX(gx3,gz2+1,2*L,M-1);
	YplusisCtimesX(gx3,gz3+1,S,M-1);
	YplusisCtimesX(gx3,gz4+1,L,M-1);

	    YisCtimesX(gx4,gz1+jx,L,M-jx);
	YplusisCtimesX(gx4,gz2+jx,2*L,M-jx);
	YplusisCtimesX(gx4,gz3+jx,L,M-jx);
	YplusisCtimesX(gx4,gz4+jx,S,M-jx);

	times(gs0,gx0,g,M);times(gs1,gx1,g,M); times(gs2,gx2,g,M); times(gs3,gx3,g,M);times(gs4,gx4,g,M);
}

void
Lat2D2ndO::PropagateF(Matrix Gi, Vector G, const int s, const double S, const bool stiff_range) const {
Message(fatal,"PropagateF with stiff-range not implemented in 2d");
}

void
Lat2D2ndO::PropagateF(Matrix Gi, Vector G, const int s) const {
Message(fatal,"PropagateF without S not implemented in 2d");
}

void
Lat2D2ndO::PropagateF(Vector Gi, Vector G, const double S) const {
	int jx;

	double L=(1.0-S)/4.0;
	double *gz0 = &Gi[1];
	double *gz1 = &Gi[1+M];
	double *gz2 = &Gi[1+2*M];
	double *gz3 = &Gi[1+3*M];
	double *gz4 = &Gi[1+4*M];

	SetBx1(gz0,gz4);SetBx1(gz1);SetBx1(gz2);SetBx1(gz3);
	SetBy1(gz1,gz3);SetBy1(gz0);SetBy1(gz2);SetBy1(gz4);
	SetBxm(gz0,gz4);SetBxm(gz1);SetBxm(gz2);SetBxm(gz3);
	SetBym(gz1,gz3);SetBym(gz0);SetBym(gz2);SetBym(gz4);

	double *gx0 = &Gx0[1];
	double *gx1 = gx0+M;
	double *gx2 = gx0+2*M;
	double *gx3 = gx0+3*M;
	double *gx4 = gx0+4*M;

	double *g = &G[1];

	jx=GetNumLayers(2);

	    YisCtimesX(gx0+jx,gz0,S,M-jx);
	YplusisCtimesX(gx0+jx,gz1,L,M-jx);
	YplusisCtimesX(gx0+jx,gz2,2*L,M-jx);
	YplusisCtimesX(gx0+jx,gz3,L,M-jx);

        YisCtimesX(gx1+1,gz0,L,M-1);
	YplusisCtimesX(gx1+1,gz1,S,M-1);
	YplusisCtimesX(gx1+1,gz2,2*L,M-1);
	YplusisCtimesX(gx1+1,gz4,L,M-1);

        YisCtimesX(gx2,gz0,L,M);
	YplusisCtimesX(gx2,gz1,L,M);
	YplusisCtimesX(gx2,gz2,S,M);
	YplusisCtimesX(gx2,gz3,L,M);
	YplusisCtimesX(gx2,gz4,L,M);

	    YisCtimesX(gx3,gz0+1,L,M-1);
	YplusisCtimesX(gx3,gz2+1,2*L,M-1);
	YplusisCtimesX(gx3,gz3+1,S,M-1);
	YplusisCtimesX(gx3,gz4+1,L,M-1);

	    YisCtimesX(gx4,gz1+jx,L,M-jx);
	YplusisCtimesX(gx4,gz2+jx,2*L,M-jx);
	YplusisCtimesX(gx4,gz3+jx,L,M-jx);
	YplusisCtimesX(gx4,gz4+jx,S,M-jx);

	times(gz0,gx0,g,M);times(gz1,gx1,g,M); times(gz2,gx2,g,M); times(gz3,gx3,g,M);times(gz4,gx4,g,M);
}

void
Lat2D2ndO::PropagateF(Vector Gi, Vector G, const double S, const bool stiffness) const {
	Message(fatal,"PropagateF vector with stiff-range not implemented in 2d");
}

void
Lat2D2ndO::PropagateF(Vector Gi, Vector G) const {
Message(fatal,"PropagateF without S not implemented in 2d");
}
void
Lat2D2ndO::PropagateB(Vector Gi, Vector G,const double S) const {
	int jx;

	double L=(1.0-S)/4.0;
	double *gz0 = &Gi[1];
	double *gz1 = &Gi[1+M];
	double *gz2 = &Gi[1+2*M];
	double *gz3 = &Gi[1+3*M];
	double *gz4 = &Gi[1+4*M];

	SetBx1(gz0,gz4);SetBx1(gz1);SetBx1(gz2);SetBx1(gz3);
	SetBy1(gz1,gz3);SetBy1(gz0);SetBy1(gz2);SetBy1(gz4);
	SetBxm(gz0,gz4);SetBxm(gz1);SetBxm(gz2);SetBxm(gz3);
	SetBym(gz1,gz3);SetBym(gz0);SetBym(gz2);SetBym(gz4);

	double *gx0 = &Gx0[1];
	double *gx1 = gx0+M;
	double *gx2 = gx0+2*M;
	double *gx3 = gx0+3*M;;
	double *gx4 = gx0+4*M;

	double *g = &G[1];

	jx=GetNumLayers(2);

	    YisCtimesX(gx0,gz0+jx,S,M-jx);
	YplusisCtimesX(gx0,gz1+1,L,M-1);
	YplusisCtimesX(gx0,gz2,2*L,M);
	YplusisCtimesX(gx0+1,gz3,L,M-1);

	    YisCtimesX(gx1,gz0+jx,L,M-jx);
	YplusisCtimesX(gx1,gz1+1,S,M-1);
	YplusisCtimesX(gx1,gz2,2*L,M);
	YplusisCtimesX(gx1+jx,gz4,L,M-jx);

	    YisCtimesX(gx2,gz0+jx,L,M-jx);
	YplusisCtimesX(gx2,gz1+1,L,M-1);
	YplusisCtimesX(gx2,gz2,S,M);
	YplusisCtimesX(gx2+1,gz3,L,M-1);
	YplusisCtimesX(gx2+jx,gz4,L,M-jx);

	    YisCtimesX(gx3,gz0+jx,L,M-jx);
	YplusisCtimesX(gx3,gz2,2*L,M);
	YplusisCtimesX(gx3+1,gz3,S,M-1);
	YplusisCtimesX(gx3+jx,gz4,L,M-jx);

	    YisCtimesX(gx4,gz1+1,L,M-1);
	YplusisCtimesX(gx4,gz2,2*L,M);
	YplusisCtimesX(gx4+1,gz3,L,M-1);
	YplusisCtimesX(gx4+jx,gz4,S,M-jx);

	times(gz0,gx0,g,M);times(gz1,gx1,g,M); times(gz2,gx2,g,M); times(gz3,gx3,g,M);times(gz4,gx4,g,M);
}
void
Lat2D2ndO::PropagateB(Vector Gi, Vector G,const double S, const bool stiffness) const {
Message(fatal,"PropagateB vector with stiff-range not implemented in 2d");
}
void
Lat2D2ndO::PropagateB(Vector Gi, Vector G) const {
Message(fatal,"PropagateB vector without S not implemented in 2d");
}

void
Lat2D2ndO::Init2G(Vector Gi1, Vector Gi2, const Vector G, const LatticeRange* LatRange) const {
	Message(fatal,"Init2G not implemented in 3d");
}
void
Lat2D2ndO::Propagate2G(Matrix Gi1, Matrix Gi2, const Vector G, const int s,const LatticeRange* LatRange) const {
Message(fatal,"Propagate 2G not implemented in 3d");
}
void
Lat2D2ndO::Propagate2G(Matrix Gi1, Matrix Gi2, const Vector G, const int s,const LatticeRange* LatRange,const double f) const {
Message(fatal,"Propagate 2G not implemented in 3d");
}

void
Lat2D2ndO::Propagate2G(Vector Gi1, Vector Gi2, const Vector G, const LatticeRange* LatRange) const {
Message(fatal,"Propagate 2G not implemented in 3d");
}
void
Lat2D2ndO::Propagate2G(Vector Gi1, Vector Gi2, const Vector G, const LatticeRange* LatRange,const double f) const {
Message(fatal,"Propagate 2G not implemented in 3d");
}

void
Lat2D2ndO::PropagateG(Matrix Gi1, const Vector G1, const int s) const {}

void
Lat2D2ndO::PropagateG(Vector G1, const Vector G2) const {}

void
Lat2D2ndO::PropagateG(const Vector Gi, const Vector G, Vector Gout) const {}

void
Lat2D2ndO::PropagateG(Matrix Gi1, const Vector G1, const int s,const double f) const {}

void
Lat2D2ndO::PropagateG(Vector G1, const Vector G2,const double f) const {}

void
Lat2D2ndO::PropagateG(const Vector Gi, const Vector G, Vector Gout, const double f) const {}

void
Lat2D2ndO::ConnectG(const Vector GiA, const Matrix GiB, const int s, Vector out) const {
	double *p_out = &out[1];
	const double *p_giB;
	const double *p_giA;
	for (int k=0; k<5; k++) {
		p_giB=&GiB[1+k*M][s];
		p_giA = &GiA[1+k*M];
		addTimes(p_out,p_giA,p_giB,M);
	}
	p_giB=&GiB[1+2*M][s];
	p_giA = &GiA[1+2*M];
	addTimes(p_out,p_giA,p_giB,M);
}

void
Lat2D2ndO::ConnectG(const Vector GiA, const Matrix GiB, const int s, Vector out, double *OUT) const {
	int M=GetTotalNumLayers();
	int kk;
	double *p_out = &out[1];
	const double *p_giB;
	const double *p_giA;
	for (int k=0; k<5; k++) {
		if (k==0 || k==5) kk=0;
		if (k==1 || k==4) kk=1;
		if (k==2) kk=2;
		p_giB=&GiB[1+k*M][s];
		p_giA = &GiA[1+k*M];
		addTimes(p_out,p_giA,p_giB,M);
		addTimes(OUT+kk*M,p_giA,p_giB,M);
	}
	p_giB=&GiB[1+2*M][s];
	p_giA = &GiA[1+2*M];
	addTimes(p_out,p_giA,p_giB,M);
	addTimes(OUT+2*M,p_giA,p_giB,M);
}

void
Lat2D2ndO::ConnectG(const Vector GiA, const Vector GiB, Vector out) const {
	cout << "ConnectG 2" << endl;
	double *p_out = &out[1];
	const double *p_giA;
	const double *p_giB;
	for (int k=0; k<5; k++) {
		p_giA = &GiA[1+k*M];
		p_giB = &GiB[1+k*M];
		addTimes(p_out,p_giA,p_giB,M);
	}
	p_giA = &GiA[1+2*M];
	p_giB = &GiB[1+2*M];
	addTimes(p_out,p_giA,p_giB,M);
}

void
Lat2D2ndO::ConnectG(const Vector GiA, const Vector GiB, Vector out, double* OUT) const {
	int M=GetTotalNumLayers();
	int kk;
	double *p_out = &out[1];
	const double *p_giA;
	const double *p_giB;
	for (int k=0; k<5; k++) {
		if (k==0 || k==5) kk=0;
		if (k==1 || k==4) kk=1;
		if (k==2) kk=2;
		p_giA = &GiA[1+k*M];
		p_giB = &GiB[1+k*M];
		addTimes(p_out,p_giA,p_giB,M);
		addTimes(OUT+kk*M,p_giA,p_giB,M);
	}
	p_giA = &GiA[1+2*M];
	p_giB = &GiB[1+2*M];
	addTimes(p_out,p_giA,p_giB,M);
	addTimes(OUT+2*M,p_giA,p_giB,M);
}

void
Lat2D2ndO::Connect2G(const Vector GiA1, const Matrix GiB1, const int s1, const Vector GiA2,
						 const Matrix GiB2, const int s2, Vector out) const {
	Message(fatal,"Connect 2G not implemented in 3d");
}

void
Lat2D2ndO::Connect2G(const Vector GiA1, const Vector GiB1, const Vector GiA2,  const Vector GiB2, Vector out) const {
	Message(fatal,"Connect 2G not implemented in 3d");
}

void
Lat2D2ndO::CorrectDoubleCountG(Vector in, const Vector G) const {
	double *p_in = &in[1];
	const double *p_G = &G[1];
    //Message(debug,"CorrectDoubleCountG(Vector in, const Vector G)");
	div(p_in,p_G,GetTotalNumLayers());
}

void
Lat2D2ndO::CorrectDoubleCountG(Vector in, double *PHI, const Vector G) const {
	double *p_in = &in[1];
	const double *p_G = &G[1];
    //Message(debug,"CorrectDoubleCountG(Vector in, const Vector G)");
	div(p_in,p_G,M);
	for (int k=0; k<3; k++) div(PHI+k*M,p_G,M);
}
//double
//Lat2D2ndO::ComputeLnGN(Vector Gi) const {
//	double value = 0;
//	int M=GetTotalNumLayers();
//	//SubtractBoundaries(Gi);
//	double *p_gi;
//	for (int k=0; k<5; k++) {
//		p_gi = &Gi[1+k*M];
//		removeboundaries(p_gi);
//		for (int i=0; i<M; i++) value +=p_gi[i];
//	}
//	p_gi = &Gi[1+2*M];
//	for (int i=0; i<M; i++) value +=p_gi[i];
//
//  	return log(value/6.0);
//}
double
Lat2D2ndO::ComputeLnGN(Vector Gi) const {
		int x,y,z;
		double value = 0;
		int jx=MY;
		for (x=1; x<MX-1; x++) {
			for (y=1; y<MY-1; y++) {
				z=jx*x+y+1;
				value += (Gi[z]+Gi[M+z]+2.0*Gi[2*M+z]+Gi[3*M+z]+Gi[4*M+z]);
			}
		}
		if (boundX1 == 3) {
			for (y=1; y<MY-1; y++) {
				x=1; z=jx*x+y+1;
				value -= 0.5*(Gi[z]+Gi[M+z]+2.0*Gi[2*M+z]+Gi[3*M+z]+Gi[4*M+z]);;
			}
		}
		if (boundX2 == MX-2) {
			for (y=1; y<MY-1; y++) {
				x=MX-2; z=jx*x+y+1;
				value -= 0.5*(Gi[z]+Gi[M+z]+2.0*Gi[2*M+z]+Gi[3*M+z]+Gi[4*M+z]);
			}
		}
		if (boundY1 == 3) {
			for (x=1; x<MX-1; x++) {
				y=1; z=jx*x+y+1;
				value -= 0.5*(Gi[z]+Gi[M+z]+2.0*Gi[2*M+z]+Gi[3*M+z]+Gi[4*M+z]);
			}
		}
		if (boundY2 == MY-2) {
			for (x=1; x<MX-1; x++) {
				y=MY-2; z=jx*x+y+1;
				value -= 0.5*(Gi[z]+Gi[M+z]+2.0*Gi[2*M+z]+Gi[3*M+z]+Gi[4*M+z]);
			}
		}
		if (boundX1 == 3 && boundY1 == 3) {
			x=1; y=1;  z=jx*x+y+1;
			value += 0.25*(Gi[z]+Gi[M+z]+2.0*Gi[2*M+z]+Gi[3*M+z]+Gi[4*M+z]);
		}
		if (boundX1 == 3 && boundY2 == MY-2) {
			x=1; y=MY-2;  z=jx*x+y+1;
			value += 0.25*(Gi[z]+Gi[M+z]+2.0*Gi[2*M+z]+Gi[3*M+z]+Gi[4*M+z]);
		}
		if (boundX2 == MX-2 && boundY1 == 3) {
			x=MX-2; y=1;  z=jx*x+y+1;
			value += 0.25*(Gi[z]+Gi[M+z]+2.0*Gi[2*M+z]+Gi[3*M+z]+Gi[4*M+z]);
		}
		if (boundX2 == MX-2 && boundY2 == MY-2) {
			x=MX-2; y=MY-2;  z=jx*x+y+1;
			value += 0.25*(Gi[z]+Gi[M+z]+2.0*Gi[2*M+z]+Gi[3*M+z]+Gi[4*M+z]);
	}
  	return log(value/6.0);
}

void
Lat2D2ndO::NormPhiFree(Vector phi, const double C) const {
	//Message(debug,"NormPhiFree(Vector phi, const double C)");
	norm(phi,C,M);
}

void
Lat2D2ndO::NormPhiFree(Vector phi, double *PHI, const double C) const {
	//Message(debug,"NormPhiFree(Vector phi, const double C)");
	norm(phi,C,M);
	for (int k=0; k<3; k++) norm(PHI+k*M,C,M);
}
void
Lat2D2ndO::NormPhiRestr(Vector phi, const Vector Gi, double C) const {
	C /= 6.0*exp(ComputeLnGN(Gi));
	norm(phi,C,M);
}

void
Lat2D2ndO::NormPhiRestr(Vector phi, double *PHI, const Vector Gi, double C) const {
	C /= 6.0*exp(ComputeLnGN(Gi));
	norm(phi,C,M);
	for (int k=0; k<3; k++) norm(PHI+k*M,C,M);
}

void
Lat2D2ndO::UpdateBoundaries(Vector A) const {
		SetBoundaries(A);
}
void
Lat2D2ndO::UpdateBoundaries(Matrix A, const int s) const {
		Message(fatal,"report error: UpdateBoundaries not implemented in lat2d2nd0");
}
