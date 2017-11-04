#include "Lat1DFlat1stO.h"

Lat1DFlat1stO::Lat1DFlat1stO(Input* MyInput_, Text name_)
	: Lat1DFlat(MyInput_,name_) {
}
Lat1DFlat1stO::~Lat1DFlat1stO() {
}
Boolean 
Lat1DFlat1stO::OverflowProtection(void) const {
	return false;
}
void 
Lat1DFlat1stO::MakeSafe(Vector) const {
}

void
Lat1DFlat1stO::GetLatticeInfo(int* Info) const{ 
	Info[0]=1; //Gradients
	Info[1]=1; //simpel cubic;
	//if (N_comp_ranges>1) {Message(fatal,"GPU not implemented for more than one comp-range. ");} 
	Info[2]=M;
	Info[3]=0;
	Info[4]=0; 
	//assumption that compranges =1;	
	Info[5]=Get_BL(1); //x lower bound
	Info[6]=Get_BL(2); //x upper bound
	Info[7]=0; //y lower bound
	Info[8]=0; //y upper bound
	Info[9]=0; //z lower bound
	Info[10]=0;//z upper bound
}

void 
Lat1DFlat1stO::RestoreFromSafe(Vector) const {
}
//! propagator for the classical matrix
void 
Lat1DFlat1stO::PropagateG(Matrix Gi, const Vector G, const int s) const {
	int z;
	double a,b,c;
	int s_1 = s-1;
	b = Gi[1][s_1];
	c = Gi[2][s_1];
	for (z=2; z<M; z++) {
		a = b; b = c; c = Gi[z+1][s_1];
		Gi[z][s] = G[z] * (b + l1*(a - b - b + c));
	}
	Gi[1][s] = Gi[bound1][s];
	Gi[M][s] = Gi[bound2][s];
}

void 
Lat1DFlat1stO::PropagateG(Matrix Gi, const Vector G, const int s, const double f) const {
	int z;
	double a,b,c;
	double norm = l1*(exp(-f)+exp(f)-2)+1;
	double f_1= l1*exp(-f)/norm;
	double f1= l1*exp(f)/norm;
	int s_1 = s-1;
	b = Gi[1][s_1];
	c = Gi[2][s_1];
	for (z=2; z<M; z++) {
		a = b; b = c; c = Gi[z+1][s_1];
		Gi[z][s] = G[z] * (b + f1*(a - b) + f_1*(c - b));
	}
	Gi[1][s] = Gi[bound1][s];
	Gi[M][s] = Gi[bound2][s];
}

void 
Lat1DFlat1stO::PropagateG(Vector Gi, const Vector G) const {
	int z;
	double a,b,c;
	b = Gi[1];
	c = Gi[2];
	for (z=2; z<M; z++) {
		a = b; b = c; c = Gi[z+1];
		Gi[z] = G[z] * (b + l1*(a - b - b + c));
	}
	Gi[1] = Gi[bound1];
	Gi[M] = Gi[bound2];	
}

void 
Lat1DFlat1stO::PropagateG(Vector Gi, const Vector G, const double f) const {
	int z;
	double a,b,c;
	double norm = l1*(exp(-f)+exp(f)-2)+1;
	double f_1= l1*exp(-f)/norm;
	double f1= l1*exp(f)/norm;
	b = Gi[1];
	c = Gi[2];
	for (z=2; z<M; z++) {
		a = b; b = c; c = Gi[z+1];
		Gi[z] = G[z] * (b + f1*(a - b) + f_1*(c - b));
	}
	Gi[1] = Gi[bound1];
	Gi[M] = Gi[bound2];	
}

void 
Lat1DFlat1stO::PropagateG(const Vector Gi, const Vector G, Vector Gout) const {
	int z;
	double a,b,c;
	b = Gi[1];
	c = Gi[2];
	for (z=2; z<M; z++) {
		a = b; b = c; c = Gi[z+1];
		Gout[z] = G[z] * (b + l1*(a - b - b + c));
	}
	Gout[1] = Gout[bound1];
	Gout[M] = Gout[bound2];	
}

void 
Lat1DFlat1stO::PropagateG(const Vector Gi, const Vector G, Vector Gout, const double f) const {
	int z;
	double a,b,c;
	double norm = l1*(exp(-f)+exp(f)-2)+1;
	double f_1= l1*exp(-f)/norm;
	double f1= l1*exp(f)/norm;
	b = Gi[1];
	c = Gi[2];
	for (z=2; z<M; z++) {
		a = b; b = c; c = Gi[z+1];
		Gout[z] = G[z] * (b + f1*(a - b) + f_1*(c - b));
	}
	Gout[1] = Gout[bound1];
	Gout[M] = Gout[bound2];	
}

void 
Lat1DFlat1stO::Init2G(Vector Gi1,Vector Gi2,const Vector G, const LatticeRange* LatRange) const {
	int z;
	for (z=2; z<M; z++) {
		if (LatRange->InRange(z)) {
			Gi1[z] = G[z];
			Gi2[z] = 0;
		} else {
			Gi1[z] = 0;
			Gi2[z] = G[z];
		}
	}
	Gi1[1] = Gi1[bound1];
	Gi2[1] = Gi2[bound1];
	Gi1[M] = Gi1[bound2];
	Gi2[M] = Gi2[bound2];
}


void 
Lat1DFlat1stO::Propagate2G(Matrix Gi1,  Matrix Gi2,  const Vector G, const int s,  const LatticeRange* LatRange) const {
	int z;
	double a1,b1,c1,a2,b2,c2;
	int s_1 = s-1;
	b1 = Gi1[bound1][s_1];
	b2 = Gi2[bound1][s_1];
	c1 = Gi1[2][s_1];
	c2 = Gi2[2][s_1];
	for (z=2; z<M; z++) {
		a1 = b1; b1 = c1; c1 = Gi1[z+1][s_1];
		a2 = b2; b2 = c2; c2 = Gi2[z+1][s_1];
		if (LatRange->InRange(z)) {
			Gi1[z][s] = G[z] * (b1 + l1*(a1 - b1 - b1 + c1 + a2 + c2));
			Gi2[z][s] = 0;
		} else {
			Gi1[z][s] = G[z] * (b1 + l1*(a1 - b1 - b1 + c1));
			Gi2[z][s] = G[z] * (b2 + l1*(a2 - b2 - b2 + c2));
		}
	}
	Gi1[1][s] = Gi1[bound1][s];
	Gi2[1][s] = Gi2[bound1][s];
	Gi1[M][s] = Gi1[bound2][s];
	Gi2[M][s] = Gi2[bound2][s];
}

void 
Lat1DFlat1stO::Propagate2G(Matrix Gi1,   Matrix Gi2, const Vector G,  const int s,const LatticeRange* LatRange, const double f) const {
	int z;
	double a1,b1,c1,a2,b2,c2;
	double norm = l1*(exp(-f)+exp(f)-2)+1;
	double f_1= l1*exp(-f)/norm;
	double f1= l1*exp(f)/norm;
	int s_1 = s-1;
	b1 = Gi1[bound1][s_1];
	b2 = Gi2[bound1][s_1];
	c1 = Gi1[2][s_1];
	c2 = Gi2[2][s_1];
	for (z=2; z<M; z++) {
		a1 = b1; b1 = c1; c1 = Gi1[z+1][s_1];
		a2 = b2; b2 = c2; c2 = Gi2[z+1][s_1];
		if (LatRange->InRange(z)) {
			Gi1[z][s] = G[z] * (b1 + f1*(a1 + a2 - b1)  + f_1*(c1 + c2 - b1));
			Gi2[z][s] = 0;
		} else {
			Gi1[z][s] = G[z] * (b1 + f1*(a1 - b1) + f_1*(c1 - b1));
			Gi2[z][s] = G[z] * (b2 + f1*(a2 - b2) + f_1*(c2 - b2));
		}
	}
	Gi1[1][s] = Gi1[bound1][s];
	Gi2[1][s] = Gi2[bound1][s];
	Gi1[M][s] = Gi1[bound2][s];
	Gi2[M][s] = Gi2[bound2][s];
}

void 
Lat1DFlat1stO::Propagate2G(Vector Gi1, Vector Gi2, const Vector G, const LatticeRange* LatRange) const {
	int z;
	double a1,b1,c1,a2,b2,c2;
	b1 = Gi1[bound1];
	b2 = Gi2[bound1];
	c1 = Gi1[2];
	c2 = Gi2[2];
	for (z=2; z<M; z++) {
		a1 = b1; b1 = c1; c1 = Gi1[z+1];
		a2 = b2; b2 = c2; c2 = Gi2[z+1];
		if (LatRange->InRange(z)) {
			Gi1[z] = G[z] * (b1 + l1*(a1 - b1 - b1 + c1 + a2 + c2));
			Gi2[z] = 0;
		} else {
			Gi1[z] = G[z] * (b1 + l1*(a1 - b1 - b1 + c1));
			Gi2[z] = G[z] * (b2 + l1*(a2 - b2 - b2 + c2));
		}
	}
	Gi1[1] = Gi1[bound1];
	Gi2[1] = Gi2[bound1];
	Gi1[M] = Gi1[bound2];
	Gi2[M] = Gi2[bound2];
}

void 
Lat1DFlat1stO::Propagate2G(Vector Gi1,  Vector Gi2,  const Vector G,   const LatticeRange* LatRange, const double f) const {
	int z;
	double a1,b1,c1,a2,b2,c2;
	double norm = l1*(exp(-f)+exp(f)-2)+1;
	double f_1= l1*exp(-f)/norm;
	double f1= l1*exp(f)/norm;
	b1 = Gi1[bound1];
	b2 = Gi2[bound1];
	c1 = Gi1[2];
	c2 = Gi2[2];
	for (z=2; z<M; z++) {
		a1 = b1; b1 = c1; c1 = Gi1[z+1];
		a2 = b2; b2 = c2; c2 = Gi2[z+1];
		if (LatRange->InRange(z)) {
			Gi1[z] = G[z] * (b1 + f1*(a1 + a2 - b1)  + f_1*(c1 + c2 - b1));
			Gi2[z] = 0;
		} else {
			Gi1[z] = G[z] * (b1 + f1*(a1 - b1) + f_1*(c1 - b1));
			Gi2[z] = G[z] * (b2 + f1*(a2 - b2) + f_1*(c2 - b2));
		}
	}
	Gi1[1] = Gi1[bound1];
	Gi2[1] = Gi2[bound1];
	Gi1[M] = Gi1[bound2];
	Gi2[M] = Gi2[bound2];
}

void 
Lat1DFlat1stO::ConnectG(const Vector GiA,const Matrix GiB, const int s, Vector out) const {
	int z;
	for (z=1; z<=M; z++) {
		out[z] += GiA[z]*GiB[z][s];
	}
}

//! Multiply GiA and GiB on all sites of the lattice 
void 
Lat1DFlat1stO::ConnectG(const Vector GiA, const Vector GiB, Vector out) const {
	int z;
	for (z=1; z<=M; z++) {
		out[z] += GiA[z]*GiB[z];
	}
}
void 
Lat1DFlat1stO::Connect2G(const Vector GiA1,  const Matrix GiB1,  const int s1, const Vector GiA2, const Matrix GiB2,const int s2, Vector out) const {
	int z;
	for (z=1; z<=M; z++) {
		out[z] += GiA1[z]*GiB2[z][s1] + GiA2[z]*GiB1[z][s2];
	}
}

//! ConnectG for the 2nd generation
void 
Lat1DFlat1stO::Connect2G(const Vector GiA1, 
						 const Vector GiB1, 
						 const Vector GiA2, 
						 const Vector GiB2, 
						 Vector out) const {
	for (int z=1; z<=M; z++) {
		out[z] += GiA1[z]*GiB2[z] + GiA2[z]*GiB1[z];
	}
}
void 
Lat1DFlat1stO::CorrectDoubleCountG(Vector in, const Vector G) const {
	for (int z=1; z<=M; z++) {
		if (G[z] > 0) in[z] /= G[z];
		else in[z] = 0;
	}
}
double
Lat1DFlat1stO::ComputeLnGN(const Vector Gi) const {
	double value = Gi[2]*layerAdjustment;
	if (bound1 == 3) value/=2;
	for (int z=3; z<=M-2; z++) {
		value += Gi[z];
	}
	if (bound2 == numLayers) value += Gi[M-1]/2; // mirror2
	else if (M > 3) {
		value += Gi[M-1]*(2-layerAdjustment);
	}
	if (numExtraLatticeSites != 0) {
		value += numExtraLatticeSites; // usually zero.
	}
	return log(value);
}
void 
Lat1DFlat1stO::NormPhiFree(Vector phi, const double norm) const {
	for (int z=1; z<=M; z++) {
		phi[z] *= norm;
	}
}
void 
Lat1DFlat1stO::NormPhiRestr(Vector phi, 
					   const Vector Gi, 
					   double norm) const {
	norm /= exp(ComputeLnGN(Gi));
	for (int z=1; z<=M; z++) {
		phi[z] *= norm;
	}
}

	
