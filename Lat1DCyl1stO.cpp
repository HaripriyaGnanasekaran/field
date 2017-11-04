#include "Lat1DCyl1stO.h"

Lat1DCyl1stO::Lat1DCyl1stO(Input* MyInput_, Text name_)
	: Lat1DFlat(MyInput_,name_), Lat1DCylinder(MyInput_,name_) {
}
Lat1DCyl1stO::~Lat1DCyl1stO() {
}

void
Lat1DCyl1stO::GetLatticeInfo(int* Info) const{ 
}

void 
Lat1DCyl1stO::PropagateG(Matrix Gi, const Vector G, const int s) const {
	int z;
	double a,b,c;
	int s_1 = s-1;
	b = Gi[1][s_1];
	c = Gi[2][s_1];
	for (z=2; z<M; z++) {
		a = b; b = c; c = Gi[z+1][s_1];
		Gi[z][s] = G[z] * (b + lambda_1[z]*(a - b) + lambda1[z]*(c - b));
	}
	Gi[1][s] = Gi[bound1][s];
	Gi[M][s] = Gi[bound2][s];
}

void 
Lat1DCyl1stO::PropagateG(Matrix Gi, const Vector G, const int s, const double f) const {
	//cout << "Implementation_error for force ensemble in non-sperical geometry" << endl; 
	double f_1 = exp(f);
	double f1 = exp(-f);
	double fnorm=(4.0+f1+f_1)/6.0;
	f1=f1/fnorm;
	f_1=f_1/fnorm; 
	double f0z,f1z,f_1z;
	int z;
	double a,b,c;
	int s_1 = s-1;
	b = Gi[1][s_1];
	c = Gi[2][s_1];
	for (z=2; z<M; z++) {
		a = b; b = c; c = Gi[z+1][s_1];
		
		f1z=lambda1[z]*f1;
		f_1z=lambda_1[z]*f_1;
		f0z=(1-lambda1[z]-lambda_1[z])/fnorm;
		
		//Gi[z][s] = G[z] * (b + lambda_1[z]*(a - b) + lambda1[z]*(c - b));
		Gi[z][s] = G[z]*(f_1z*a+f0z*b+f1z*c);
	}
	Gi[1][s] = Gi[bound1][s];
	Gi[M][s] = Gi[bound2][s];
}

void 
Lat1DCyl1stO::PropagateG(Vector Gi, const Vector G) const {
	int z;
	double a,b,c;
	b = Gi[1];
	c = Gi[2];
	for (z=2; z<M; z++) {
		a = b; b = c; c = Gi[z+1];
		Gi[z] = G[z] * (b + lambda_1[z]*(a - b) + lambda1[z]*(c - b));
	}
	Gi[1] = Gi[bound1];
	Gi[M] = Gi[bound2];	
}

void 
Lat1DCyl1stO::PropagateG(Vector Gi, const Vector G, const double f) const {
//cout << "Implementation_error for force ensemble in non-sperical geometry" << endl; 
    double f_1 = exp(f);
	double f1 = exp(-f);
	double fnorm=(4.0+f1+f_1)/6.0;
	f1=f1/fnorm;
	f_1=f_1/fnorm; 
	double f0z,f1z,f_1z;
	int z;
	double a,b,c;
	b = Gi[1];
	c = Gi[2];
	for (z=2; z<M; z++) {
		a = b; b = c; c = Gi[z+1];
		f1z=lambda1[z]*f1;
		f_1z=lambda_1[z]*f_1;
		f0z=(1-lambda1[z]-lambda_1[z])/fnorm;
		
		//Gi[z] = G[z] * (b + lambda_1[z]*(a - b) + lambda1[z]*(c - b));
		Gi[z] = G[z]*(f_1z*a+f0z*b+f1z*c);
	}
	Gi[1] = Gi[bound1];
	Gi[M] = Gi[bound2];	
}


void 
Lat1DCyl1stO::PropagateG(const Vector Gi, const Vector G, Vector Gout) const {
	int z;
	double a,b,c;
	b = Gi[1];
	c = Gi[2];
	for (z=2; z<M; z++) {
		a = b; b = c; c = Gi[z+1];
		Gout[z] = G[z] * (b + lambda_1[z]*(a - b) + lambda1[z]*(c - b));
	}
	Gout[1] = Gout[bound1];
	Gout[M] = Gout[bound2];	
}

void 
Lat1DCyl1stO::PropagateG(const Vector Gi, const Vector G, Vector Gout, const double f) const {
	double f_1 = exp(f);
	double f1 = exp(-f);
	double fnorm=(4.0+f1+f_1)/6.0;
	f1=f1/fnorm;
	f_1=f_1/fnorm; 
	double f0z,f1z,f_1z;
	int z;
	double a,b,c;
	b = Gi[1];
	c = Gi[2];
	for (z=2; z<M; z++) {
		a = b; b = c; c = Gi[z+1];
		f1z=lambda1[z]*f1;
		f_1z=lambda_1[z]*f_1;
		f0z=(1-lambda1[z]-lambda_1[z])/fnorm;
		//Gout[z] = G[z] * (b + lambda_1[z]*(a - b) + lambda1[z]*(c - b));
		Gout[z] = G[z]*(f_1z*a+f0z*b+f1z*c);
	}
	Gout[1] = Gout[bound1];
	Gout[M] = Gout[bound2];	
}

void 
Lat1DCyl1stO::Propagate2G(Matrix Gi1, Matrix Gi2, const Vector G, const int s, const LatticeRange* LatRange) const {
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
			Gi1[z][s] = G[z] * (b1 + lambda_1[z]*(a1 + a2 - b1) + lambda1[z]*(c1 + c2 - b1));
			Gi2[z][s] = 0;
		} else {
			Gi1[z][s] = G[z] * (b1 + lambda_1[z]*(a1 - b1) + lambda1[z]*(c1 - b1));
			Gi2[z][s] = G[z] * (b2 + lambda_1[z]*(a2 - b2) + lambda1[z]*(c2 - b2));
		}
	}
	Gi1[1][s] = Gi1[bound1][s];
	Gi2[1][s] = Gi2[bound1][s];
	Gi1[M][s] = Gi1[bound2][s];
	Gi2[M][s] = Gi2[bound2][s];
}

void 
Lat1DCyl1stO::Propagate2G(Matrix Gi1, Matrix Gi2, const Vector G, const int s, const LatticeRange* LatRange, const double f) const {
	double f_1 = exp(f);
	double f1 = exp(-f);
	double fnorm=(4.0+f1+f_1)/6.0;
	f1=f1/fnorm;
	f_1=f_1/fnorm; 
	double f0z,f1z,f_1z;
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
		f1z=lambda1[z]*f1;
		f_1z=lambda_1[z]*f_1;
		f0z=(1-lambda1[z]-lambda_1[z])/fnorm;
		if (LatRange->InRange(z)) {
			//Gi1[z][s] = G[z] * (b1 + lambda_1[z]*(a1 + a2 - b1) + lambda1[z]*(c1 + c2 - b1));
			Gi1[z][s] = G[z]*(f_1z*(a1+a2)+f0z*b1+f1z*(c1+c2));
			Gi2[z][s] = 0;
		} else {
			//Gi1[z][s] = G[z] * (b1 + lambda_1[z]*(a1 - b1) + lambda1[z]*(c1 - b1));
			//Gi2[z][s] = G[z] * (b2 + lambda_1[z]*(a2 - b2) + lambda1[z]*(c2 - b2));
			Gi1[z][s] = G[z]*(f_1z*(a1)+f0z*b1+f1z*(c1));
			Gi2[z][s] = G[z]*(f_1z*(a2)+f0z*b2+f1z*(c2));
		}
	}
	Gi1[1][s] = Gi1[bound1][s];
	Gi2[1][s] = Gi2[bound1][s];
	Gi1[M][s] = Gi1[bound2][s];
	Gi2[M][s] = Gi2[bound2][s];
}


void
Lat1DCyl1stO::Propagate2G(Vector Gi1, Vector Gi2, const Vector G, const LatticeRange* LatRange) const {
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
			Gi1[z] = G[z] * (b1 + lambda_1[z]*(a1 + a2 - b1) + lambda1[z]*(c1 + c2 - b1));
			Gi2[z] = 0;
		} else {
			Gi1[z] = G[z] * (b1 + lambda_1[z]*(a1 - b1) + lambda1[z]*(c1 - b1));
			Gi2[z] = G[z] * (b2 +  lambda_1[z]*(a2 - b2) + lambda1[z]*(c2 - b2));
		}
	}
	Gi1[1] = Gi1[bound1];
	Gi2[1] = Gi2[bound1];
	Gi1[M] = Gi1[bound2];
	Gi2[M] = Gi2[bound2];
}

void
Lat1DCyl1stO::Propagate2G(Vector Gi1, Vector Gi2, const Vector G, const LatticeRange* LatRange,const double f) const {
	double f_1 = exp(f);
	double f1 = exp(-f);
	double fnorm=(4.0+f1+f_1)/6.0;
	f1=f1/fnorm;
	f_1=f_1/fnorm; 
	double f0z,f1z,f_1z;
	int z;
	double a1,b1,c1,a2,b2,c2;
	b1 = Gi1[bound1];
	b2 = Gi2[bound1];
	c1 = Gi1[2];
	c2 = Gi2[2];
	for (z=2; z<M; z++) {
		a1 = b1; b1 = c1; c1 = Gi1[z+1];
		a2 = b2; b2 = c2; c2 = Gi2[z+1];
		f1z=lambda1[z]*f1;
		f_1z=lambda_1[z]*f_1;
		f0z=(1-lambda1[z]-lambda_1[z])/fnorm;
		if (LatRange->InRange(z)) {
			//Gi1[z] = G[z] * (b1 + lambda_1[z]*(a1 + a2 - b1) + lambda1[z]*(c1 + c2 - b1));
			Gi1[z] = G[z]*(f_1z*(a1+a2)+f0z*b1+f1z*(c1+c2));
			Gi2[z] = 0;
		} else {
			//Gi1[z] = G[z] * (b1 + lambda_1[z]*(a1 - b1) + lambda1[z]*(c1 - b1));
			//Gi2[z] = G[z] * (b2 +  lambda_1[z]*(a2 - b2) + lambda1[z]*(c2 - b2));
			Gi1[z] = G[z]*(f_1z*(a1)+f0z*b1+f1z*(c1));
			Gi2[z] = G[z]*(f_1z*(a2)+f0z*b2+f1z*(c2));
		}
	}
	Gi1[1] = Gi1[bound1];
	Gi2[1] = Gi2[bound1];
	Gi1[M] = Gi1[bound2];
	Gi2[M] = Gi2[bound2];
}


double
Lat1DCyl1stO::ComputeLnGN(const Vector Gi) const {
	int z;
	double value = Gi[2]*L[2];
	if (bound1 == 3) value/=2;
	for (z=3; z<=M-2; z++) {
		value += Gi[z]*L[z];
	}
	if (bound2 == M-2) value += Gi[M-1]*L[M-1]/2;
	else if (M > 3) value += Gi[M-1]*L[M-1];
	value += numExtraLatticeSites; // usually zero.
	return log(value);
}
void 
Lat1DCyl1stO::NormPhiRestr(Vector phi, 
					  const Vector Gi, 
					  double norm) const {
	int z;
	norm /= exp(ComputeLnGN(Gi));
	for (z=1; z<=M; z++) {
		phi[z] *= norm;
	}
}
