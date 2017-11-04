#include "Lat1DCyl1stOS.h"

Lat1DCyl1stOS::Lat1DCyl1stOS(Input* MyInput_, Text name_)
	: Lat1DFlat(MyInput_,name_), Lat1DCylinder(MyInput_,name_) {
}
Lat1DCyl1stOS::~Lat1DCyl1stOS() {
}
void
Lat1DCyl1stOS::GetLatticeInfo(int* Info) const{ 
}

void
Lat1DCyl1stOS::PropagateG(Matrix Gi, const Vector G, const int s) const {
	int z;
	double a,b,c;
	int s_1 = s-1;
	b = Gi[1][s_1];
	c = Gi[2][s_1];
	for (z=2; z<M; z++) {
		a = b; b = c; c = Gi[z+1][s_1];
		if (G[z] == LOGZERO) {
			Gi[z][s] = LOGZERO;
			continue;
		}
		double lambda0 = 1-lambda1[z]-lambda_1[z];
		if (b > LOGZERO && lambda0 > 0) {
			Gi[z][s] = G[z] + log(lambda0) 
			+ b + log1p((lambda_1[z]/lambda0)*exp(a-b) 
			+ (lambda1[z]/lambda0)*exp(c-b));
		} else if (a > LOGZERO) {
			Gi[z][s] = G[z] + log(lambda_1[z]) + a
			+ log1p((lambda1[z]/lambda_1[z])*exp(c-a));
		} else if (c > LOGZERO) {
			Gi[z][s] = G[z] + log(lambda1[z]) + c;
		} else {
			Gi[z][s] = LOGZERO;
		}		
	}
	Gi[1][s] = LOGZERO;
	Gi[M][s] = LOGZERO;
	Gi[1][s] = Gi[bound1][s];
	Gi[M][s] = Gi[bound2][s];
}
void
Lat1DCyl1stOS::PropagateG(const Vector Gi, const Vector G, Vector Gout) const {
	int z;
	double a,b,c;
	b = Gi[1];
	c = Gi[2];
	for (z=2; z<M; z++) {
		a = b; b = c; c = Gi[z+1];
		if (G[z] == LOGZERO) {
			Gout[z] = LOGZERO;
			continue;
		}
		double lambda0 = 1-lambda1[z]-lambda_1[z];
		if (b > LOGZERO && lambda0 > 0) {
			Gout[z] = G[z] + log(lambda0) 
			+ b + log1p((lambda_1[z]/lambda0)*exp(a-b) 
			+ (lambda1[z]/lambda0)*exp(c-b));
		} else if (a > LOGZERO) {
			Gout[z] = G[z] + log(lambda_1[z]) + a
			+ log1p((lambda1[z]/lambda_1[z])*exp(c-a));
		} else if (c > LOGZERO) {
			Gout[z] = G[z] + log(lambda1[z]) + c;
		} else {
			Gout[z] = LOGZERO;
		}		
	}
	Gout[1] = LOGZERO;
	Gout[M] = LOGZERO;
	Gout[1] = Gout[bound1];
	Gout[M] = Gout[bound2];
}
void
Lat1DCyl1stOS::PropagateG(Vector Gi, const Vector G) const {
	int z;
	double a,b,c;
	b = Gi[1];
	c = Gi[2];
	for (z=2; z<M; z++) {
		a = b; b = c; c = Gi[z+1];
		if (G[z] == LOGZERO) {
			Gi[z] = LOGZERO;
			continue;
		}
		double lambda0 = 1-lambda1[z]-lambda_1[z];
		if (b > LOGZERO && lambda0 > 0) {
			Gi[z] = G[z] + log(lambda0) 
			+ b + log1p((lambda_1[z]/lambda0)*exp(a-b) 
			+ (lambda1[z]/lambda0)*exp(c-b));
		} else if (a > LOGZERO) {
			Gi[z] = G[z] + log(lambda_1[z]) + a
			+ log1p((lambda1[z]/lambda_1[z])*exp(c-a));
		} else if (c > LOGZERO) {
			Gi[z] = G[z] + log(lambda1[z]) + c;
		} else {
			Gi[z] = LOGZERO;
		}		
	}
	Gi[1] = LOGZERO;
	Gi[M] = LOGZERO;
	Gi[1] = Gi[bound1];
	Gi[M] = Gi[bound2];
}
void
Lat1DCyl1stOS::Propagate2G(Matrix Gi1, Matrix Gi2, const Vector G, const int s, const LatticeRange* LatRange) const {
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
		double L_1L0 = (lambda_1[z]/(1-lambda1[z]-lambda_1[z]));
		double L1L0 = (lambda1[z]/(1-lambda1[z]-lambda_1[z]));
		double L1L_1 = (lambda1[z]/lambda_1[z]);
		if (LatRange->InRange(z)) {
 			Gi2[z][s] = LOGZERO;
 			if (b1 > LOGZERO && l1 < 0.5) {
				Gi1[z][s] = G[z] + log(1-lambda1[z]-lambda_1[z]) + b1 
 				+ log1p(L_1L0*(exp(a1-b1)+exp(a2-b1)) + L1L0*(exp(c1-b1)+exp(c2-b1)));
			} else if (a1 > LOGZERO) {
				Gi1[z][s] = G[z] + log(lambda_1[z]) + a1 
 				+ log1p(L1L_1*(exp(c1-a1)+exp(c2-a1)) + exp(a2-a1));
			} else if (a2 > LOGZERO) {
				Gi1[z][s] = G[z] + log(lambda_1[z]) + a2 
 				+ log1p(L1L_1*(exp(c1-a2)+exp(c2-a2)));
			} else if (c1 > LOGZERO) {
				Gi1[z][s] = G[z] + log(lambda1[z]) + c1 
 				+ log1p(exp(c2-c1));
			} else if (c2 > LOGZERO) {
				Gi1[z][s] = G[z] + log(lambda1[z]) + c2;
			} else {
				Gi1[z][s] = LOGZERO;
			}
		} else {
			if (b1 > LOGZERO && l1 < 0.5) {
				Gi1[z][s] = G[z] + log(1-lambda1[z]-lambda_1[z]) + b1 
				+ log1p(L_1L0*exp(a1-b1) + L1L0*exp(c1-b1));
			} else if (a1 > LOGZERO) {
				Gi1[z][s] = G[z] + log(lambda_1[z]) + a1 
 				+ log1p(L1L_1*exp(c1-a1));
			} else if (c1 > LOGZERO) {
				Gi1[z][s] = G[z] + log(lambda1[z]) + c1;
			} else {
				Gi1[z][s] = LOGZERO;
			}
			if (b2 > LOGZERO && l1 < 0.5) {
				Gi2[z][s] = G[z] + log(1-lambda1[z]-lambda_1[z]) + b2 
				+ log1p(L_1L0*exp(a2-b2) + L1L0*exp(c2-b2));
			} else if (a2 > LOGZERO) {
				Gi2[z][s] = G[z] + log(lambda_1[z]) + a2 
 				+ log1p(L1L_1*exp(c2-a2));
			} else if (c2 > LOGZERO) {
				Gi2[z][s] = G[z] + log(lambda1[z]) + c2; 
			} else {
				Gi2[z][s] = LOGZERO;
			}
		}
	}
	Gi1[1][s] = LOGZERO;
	Gi1[M][s] = LOGZERO;
	Gi2[1][s] = LOGZERO;
	Gi2[M][s] = LOGZERO;
	Gi1[1][s] = Gi1[bound1][s];
	Gi2[1][s] = Gi2[bound1][s];
	Gi1[M][s] = Gi1[bound2][s];
	Gi2[M][s] = Gi2[bound2][s];
}
void 
Lat1DCyl1stOS::Propagate2G(Vector Gi1, Vector Gi2, const Vector G, const LatticeRange* LatRange) const {
	int z;
	double a1,b1,c1,a2,b2,c2;
	b1 = Gi1[bound1];
	b2 = Gi2[bound1];
	c1 = Gi1[2];
	c2 = Gi2[2];
	for (z=2; z<M; z++) {
		a1 = b1; b1 = c1; c1 = Gi1[z+1];
		a2 = b2; b2 = c2; c2 = Gi2[z+1];
		double L_1L0 = (lambda_1[z]/(1-lambda1[z]-lambda_1[z]));
		double L1L0 = (lambda1[z]/(1-lambda1[z]-lambda_1[z]));
		double L1L_1 = (lambda1[z]/lambda_1[z]);
		if (LatRange->InRange(z)) {
 			Gi2[z] = LOGZERO;
 			if (b1 > LOGZERO && l1 < 0.5) {
				Gi1[z] = G[z] + log(1-lambda1[z]-lambda_1[z]) + b1 
 				+ log1p(L_1L0*(exp(a1-b1)+exp(a2-b1)) + L1L0*(exp(c1-b1)+exp(c2-b1)));
			} else if (a1 > LOGZERO) {
				Gi1[z] = G[z] + log(lambda_1[z]) + a1 
 				+ log1p(L1L_1*(exp(c1-a1)+exp(c2-a1)) + exp(a2-a1));
			} else if (a2 > LOGZERO) {
				Gi1[z] = G[z] + log(lambda_1[z]) + a2 
 				+ log1p(L1L_1*(exp(c1-a2)+exp(c2-a2)));
			} else if (c1 > LOGZERO) {
				Gi1[z] = G[z] + log(lambda1[z]) + c1 
 				+ log1p(exp(c2-c1));
			} else if (c2 > LOGZERO) {
				Gi1[z] = G[z] + log(lambda1[z]) + c2;
			} else {
				Gi1[z] = LOGZERO;
			}
		} else {
			if (b1 > LOGZERO && l1 < 0.5) {
				Gi1[z] = G[z] + log(1-lambda1[z]-lambda_1[z]) + b1 
				+ log1p(L_1L0*exp(a1-b1) + L1L0*exp(c1-b1));
			} else if (a1 > LOGZERO) {
				Gi1[z] = G[z] + log(lambda_1[z]) + a1 
 				+ log1p(L1L_1*exp(c1-a1));
			} else if (c1 > LOGZERO) {
				Gi1[z] = G[z] + log(lambda1[z]) + c1;
			} else {
				Gi1[z] = LOGZERO;
			}
			if (b2 > LOGZERO && l1 < 0.5) {
				Gi2[z] = G[z] + log(1-lambda1[z]-lambda_1[z]) + b2 
				+ log1p(L_1L0*exp(a2-b2) + L1L0*exp(c2-b2));
			} else if (a2 > LOGZERO) {
				Gi2[z] = G[z] + log(lambda_1[z]) + a2 
 				+ log1p(L1L_1*exp(c2-a2));
			} else if (c2 > LOGZERO) {
				Gi2[z] = G[z] + log(lambda1[z]) + c2; 
			} else {
				Gi2[z] = LOGZERO;
			}
		}
	}
	Gi1[1] = LOGZERO;
	Gi1[M] = LOGZERO;
	Gi2[1] = LOGZERO;
	Gi2[M] = LOGZERO;
	Gi1[1] = Gi1[bound1];
	Gi2[1] = Gi2[bound1];
	Gi1[M] = Gi1[bound2];
	Gi2[M] = Gi2[bound2];
}
double
Lat1DCyl1stOS::ComputeLnGN(const Vector Gi) const {
	int z;
	double value = LOGZERO;
	if (Gi[2] > LOGZERO) {
		value = log(L[2])+Gi[2];
		if (bound1 == 3) value-=log(2.0);
	}
	for (z=3; z<=M-2; z++) {
		if (Gi[z] > LOGZERO) {
			double x = log(L[z])+Gi[z]-value;
			if (x > LOGMAXDOUBLE-1 || value == LOGZERO) {		
				value = log(L[z])+Gi[z];
			} else {
				value += log1p(exp(x));
			}
		}
	}
	if (bound2 == M-2 && Gi[M-1] > LOGZERO) {
		double x = log(L[M-1])+Gi[M-1]-value-log(2.0);
		if (x > LOGMAXDOUBLE-1 || value == LOGZERO) {
			value = log(L[M-1])+Gi[M-1]-log(2.0);
		} else {
			value += log1p(exp(x));
		}
	} else if (M > 3) {
		if (Gi[M-1] > LOGZERO) {
			double x = log(L[M-1])+Gi[M-1]-value;
			if (x > LOGMAXDOUBLE-1 || value == LOGZERO) {
				value = log(L[M-1])+Gi[M-1];
			} else {
				value += log1p(exp(x));
			}
		}
	}
	if (numExtraLatticeSites != 0) {
		double x = log(numExtraLatticeSites) - value;
		value += log1p(exp(x));
	}
	return value;
}
void 
Lat1DCyl1stOS::NormPhiRestr(Vector phi, 
					   const Vector Gi, 
					   double norm) const {
	double lnGN = ComputeLnGN(Gi);
	double logNorm = log(norm)-lnGN;
	for (int z=1; z<=M; z++) {
		if (phi[z] > LOGZERO) {
			phi[z] += logNorm;
		}
	}
}
