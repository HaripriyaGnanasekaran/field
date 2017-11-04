#include "Lat1DFlat1stOS.h"

Lat1DFlat1stOS::Lat1DFlat1stOS(Input* MyInput_, Text name_)
	: Lat1DFlat(MyInput_,name_) {
	LogL1L0 = log(l1/(1-2*l1));
	LogL0 = log(1-2*l1);
	LogL1 = log(l1);

}
Lat1DFlat1stOS::~Lat1DFlat1stOS() {
}
Boolean 
Lat1DFlat1stOS::OverflowProtection(void) const {
	return true;
}

void
Lat1DFlat1stOS::GetLatticeInfo(int* Info) const{ 
}

void
Lat1DFlat1stOS::MakeSafe(Vector A) const {
	int z;
	double y;

	for (z=1; z<=M; z++) {
		if (A[z] > 0) {
			if (A[z] > CUTOFF) {
				y = (log(A[z]) - LOGCUTOFF)/LOGMAXPREFACTOR;
				A[z] = LOGCUTOFF
					+ 0.5*LOGMAXPREFACTOR*log((1+y)/(1-y));
			} else {
				A[z] = log(A[z]);
			}
		} else {
			A[z] = LOGZERO;
		}
	}
}
void
Lat1DFlat1stOS::RestoreFromSafe(Vector A) const {
	int z;

	for (z=1; z<=M; z++) {
		if (A[z] > LOGZERO) {
			if (A[z] > LOGCUTOFF) {
				A[z] = exp(LOGCUTOFF+LOGMAXPREFACTOR
				*tanh((A[z]-LOGCUTOFF)/LOGMAXPREFACTOR));;
			} else {
				A[z] = exp(A[z]);
			}
		} else {
			A[z] = 0;
		}
	}
}
void
Lat1DFlat1stOS::PropagateG(Matrix Gi, const Vector G, const int s) const {
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
		if (b > LOGZERO && l1 < 0.5) {
			Gi[z][s] = G[z] + LogL0 + b + log1p(exp(LogL1L0+a-b) + exp(LogL1L0+c-b)); 
		} else if (a > LOGZERO) {
			Gi[z][s] = G[z] + LogL1 + a + log1p(exp(c-a));
		} else if (c > LOGZERO) {
			Gi[z][s] = G[z] + LogL1 + c;
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
Lat1DFlat1stOS::PropagateG(Matrix Gi, const Vector G, const int s, const double f) const {
	Message(fatal,MyInput,"Force ensemble not implemented in Lat1DFlat1stOS");
}

void
Lat1DFlat1stOS::PropagateG(Vector Gi, const Vector G) const {
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
		if (b > LOGZERO && l1 < 0.5) {
			Gi[z] = G[z] + LogL0 + b + log1p(exp(LogL1L0+a-b) + exp(LogL1L0+c-b)); 
		} else if (a > LOGZERO) {
			Gi[z] = G[z] + LogL1 + a + log1p(exp(c-a));
		} else if (c > LOGZERO) {
			Gi[z] = G[z] + LogL1 + c;
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
Lat1DFlat1stOS::PropagateG(Vector Gi, const Vector G, const double f) const {
	Message(fatal,MyInput,"Force ensemble not implemented in Lat1DFlat1stOS");
}

void
Lat1DFlat1stOS::PropagateG(const Vector Gi, const Vector G, Vector Gout) const {
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
		if (b > LOGZERO && l1 < 0.5) {
			Gout[z] = G[z] + LogL0 + b + log1p(exp(LogL1L0+a-b) + exp(LogL1L0+c-b)); 
		} else if (a > LOGZERO) {
			Gout[z] = G[z] + LogL1 + a + log1p(exp(c-a));
		} else if (c > LOGZERO) {
			Gout[z] = G[z] + LogL1 + c;
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
Lat1DFlat1stOS::PropagateG(const Vector Gi, const Vector G, Vector Gout, const double f) const {
	Message(fatal,MyInput,"Force ensemble not implemented in Lat1DFlat1stOS");
}

void
Lat1DFlat1stOS::Init2G(Vector Gi1, 
					   Vector Gi2, 
					   const Vector G, 
					   const LatticeRange* LatRange) const {
	int z;
	for (z=2; z<M; z++) {
		if (LatRange->InRange(z)) {
			Gi1[z] = G[z];
			Gi2[z] = LOGZERO;
		} else {
			Gi1[z] = LOGZERO;
			Gi2[z] = G[z];
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
void
Lat1DFlat1stOS::Propagate2G(Matrix Gi1, 
							Matrix Gi2, 
							const Vector G, 
							const int s, 
							const LatticeRange* LatRange) const {
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
		if (G[z] == LOGZERO) {
			Gi1[z][s] = LOGZERO;
			Gi2[z][s] = LOGZERO;
			continue;
		}
		if (LatRange->InRange(z)) {
			Gi2[z][s] = LOGZERO;
			if (b1 > LOGZERO && l1 < 0.5) {
				Gi1[z][s] = G[z] + LogL0 + b1 + log1p(exp(LogL1L0+a1-b1) 
				+ exp(LogL1L0+c1-b1) + exp(LogL1L0+a2-b1) + exp(LogL1L0+c2-b1));
			} else if (a1 > LOGZERO) {
				Gi1[z][s] = G[z] + LogL1 + a1 + log1p(exp(c1-a1) 
				+ exp(a2-a1) + exp(c2-a1));
			} else if (c1 > LOGZERO) {
				Gi1[z][s] = G[z] + LogL1 + c1 + log1p(exp(a2-c1) + exp(c2-c1));
			} else if (a2 > LOGZERO) {
				Gi1[z][s] = G[z] + LogL1 + a2 + log1p(exp(c2-a2));
			} else if (c2 > LOGZERO) {
				Gi1[z][s] = G[z] + LogL1 + c2;
			} else {
				Gi1[z][s] = LOGZERO;
			}
		} else {
			if (b1 > LOGZERO && l1 < 0.5) {
				Gi1[z][s] = G[z] + LogL0 + b1 + log1p(exp(LogL1L0+a1-b1) 
				+ exp(LogL1L0+c1-b1));
			} else if (a1 > LOGZERO) {
				Gi1[z][s] = G[z] + LogL1 + a1 + log1p(exp(c1-a1));
			} else if (c1 > LOGZERO) {
				Gi1[z][s] = G[z] + LogL1 + c1;
			} else {
				Gi1[z][s] = LOGZERO;
			}
			if (b2 > LOGZERO && l1 < 0.5) {
				Gi2[z][s] = G[z] + LogL0 + b2+ log1p(exp(LogL1L0+a2-b2) 
				+ exp(LogL1L0+c2-b2));
			} else if (a2 > LOGZERO) {
				Gi2[z][s] = G[z] + LogL1 + a2 + log1p(exp(c2-a2));
			} else if (c2 > LOGZERO) {
				Gi2[z][s] = G[z] + LogL1 + c2;
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
Lat1DFlat1stOS::Propagate2G(Matrix Gi1, 
							Matrix Gi2, 
							const Vector G, 
							const int s, 
							const LatticeRange* LatRange,
							const double f) const {
	Message(fatal,MyInput,"Force ensemble not implemented in Lat1DFlat1stOS");
}


void
Lat1DFlat1stOS::Propagate2G(Vector Gi1, 
							Vector Gi2, 
							const Vector G, 
							const LatticeRange* LatRange) const {
	int z;
	double a1,b1,c1,a2,b2,c2;
	b1 = Gi1[bound1];
	b2 = Gi2[bound1];
	c1 = Gi1[2];
	c2 = Gi2[2];
	for (z=2; z<M; z++) {
		a1 = b1; b1 = c1; c1 = Gi1[z+1];
		a2 = b2; b2 = c2; c2 = Gi2[z+1];
		if (G[z] == LOGZERO) {
			Gi1[z] = LOGZERO;
			Gi2[z] = LOGZERO;
			continue;
		}
		if (LatRange->InRange(z)) {
			Gi2[z] = LOGZERO;
			if (b1 > LOGZERO && l1 < 0.5) {
				Gi1[z] = G[z] + LogL0 + b1 + log1p(exp(LogL1L0+a1-b1) 
				+ exp(LogL1L0+c1-b1) + exp(LogL1L0+a2-b1) + exp(LogL1L0+c2-b1));
			} else if (a1 > LOGZERO) {
				Gi1[z] = G[z] + LogL1 + a1 + log1p(exp(c1-a1) 
				+ exp(a2-a1) + exp(c2-a1));
			} else if (c1 > LOGZERO) {
				Gi1[z] = G[z] + LogL1 + c1 + log1p(exp(a2-c1) + exp(c2-c1));
			} else if (a2 > LOGZERO) {
				Gi1[z] = G[z] + LogL1 + a2 + log1p(exp(c2-a2));
			} else if (c2 > LOGZERO) {
				Gi1[z] = G[z] + LogL1 + c2;
			} else {
				Gi1[z] = LOGZERO;
			}
		} else {
			if (b1 > LOGZERO && l1 < 0.5) {
				Gi1[z] = G[z] + LogL0 + b1 + log1p(exp(LogL1L0+a1-b1) 
				+ exp(LogL1L0+c1-b1));
			} else if (a1 > LOGZERO) {
				Gi1[z] = G[z] + LogL1 + a1 + log1p(exp(c1-a1));
			} else if (c1 > LOGZERO) {
				Gi1[z] = G[z] + LogL1 + c1;
			} else {
				Gi1[z] = LOGZERO;
			}
			if (b2 > LOGZERO && l1 < 0.5) {
				Gi2[z] = G[z] + LogL0 + b2+ log1p(exp(LogL1L0+a2-b2) 
				+ exp(LogL1L0+c2-b2));
			} else if (a2 > LOGZERO) {
				Gi2[z] = G[z] + LogL1 + a2 + log1p(exp(c2-a2));
			} else if (c2 > LOGZERO) {
				Gi2[z] = G[z] + LogL1 + c2;
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

void
Lat1DFlat1stOS::Propagate2G(Vector Gi1, 
							Vector Gi2, 
							const Vector G, 
							const LatticeRange* LatRange,
							const double f) const {
	Message(fatal,MyInput,"Force ensemble not implemented in Lat1DFlat1stOS");
}

void
Lat1DFlat1stOS::ConnectG(const Vector GiA, 
						 const Matrix GiB, 
						 const int s, 
						 Vector out) const {
	for (int z=1; z<=M; z++) {
		if (GiA[z] > LOGZERO && GiB[z][s] > LOGZERO) {
			double x = GiA[z]+GiB[z][s]-out[z];
			if (x > LOGMAXDOUBLE-1 || out[z] == LOGZERO) {
				out[z] = GiA[z]+GiB[z][s];
			} else {
				out[z] += log1p(exp(x));
			}
		}
	}
}


void
Lat1DFlat1stOS::ConnectG(const Vector GiA, 
						 const Vector GiB, 
						 Vector out) const {
	for (int z=1; z<=M; z++) {
		if (GiA[z] > LOGZERO && GiB[z] > LOGZERO) {
			double x = GiA[z]+GiB[z]-out[z];
			if (x > LOGMAXDOUBLE-1 || out[z] == LOGZERO) {
				out[z] = GiA[z]+GiB[z];
			} else {
				out[z] += log1p(exp(x));
			}
		}
	}
}
void
Lat1DFlat1stOS::Connect2G(const Vector GiA1, 
						  const Matrix GiB1, 
						  const int s1, 
						  const Vector GiA2, 
						  const Matrix GiB2, 
						  const int s2, 
						  Vector out) const {
	double x;
	for (int z=1; z<=M; z++) {
		if (GiA1[z] > LOGZERO && GiB2[z][s1] > LOGZERO) {
			x = GiA1[z]+GiB2[z][s1]-out[z];
			if (x > LOGMAXDOUBLE-1 || out[z] == LOGZERO) {
				out[z] = GiA1[z]+GiB2[z][s1];
			} else {
				out[z] += log1p(exp(x));
			}
		}
		if (GiA2[z] > LOGZERO && GiB1[z][s2] > LOGZERO) {
			x = GiA2[z]+GiB1[z][s2]-out[z];
			if (x > LOGMAXDOUBLE-1 || out[z] == LOGZERO) {
				out[z] = GiA2[z]+GiB1[z][s2];
			} else {
				out[z] += log1p(exp(x));
			}
		}
	}
}
void
Lat1DFlat1stOS::Connect2G(const Vector GiA1, 
						  const Vector GiB1, 
						  const Vector GiA2, 
						  const Vector GiB2, 
						  Vector out) const {
	double x;
	for (int z=1; z<=M; z++) {
		if (GiA1[z] > LOGZERO && GiB2[z] > LOGZERO) {
			x = GiA1[z]+GiB2[z]-out[z];
			if (x > LOGMAXDOUBLE-1 || out[z] == LOGZERO) {
				out[z] = GiA1[z]+GiB2[z];
			} else {
				out[z] += log1p(exp(x));
			}
			x = GiA2[z]+GiB1[z]-out[z];
			if (x > LOGMAXDOUBLE-1 || out[z] == LOGZERO) {
				out[z] = GiA2[z]+GiB1[z];
			} else {
				out[z] += log1p(exp(x));
			}
		}
	}
}
void
Lat1DFlat1stOS::CorrectDoubleCountG(Vector in, const Vector G) const {
	for (int z=1; z<=M; z++) {
		if (G[z] > LOGZERO && in[z] > LOGZERO) {
			in[z] -= G[z];
		} else {
			in[z] = LOGZERO;
		}
	}
}
double
Lat1DFlat1stOS::ComputeLnGN(const Vector Gi) const {
	int z;
	double value = LOGZERO;
	if (Gi[2] > LOGZERO) {
		value = Gi[2];
		if (bound1 == 3) value-=log(2.0);
	}
	for (z=3; z<=M-2; z++) {
		if (Gi[z] > LOGZERO) {
			double x = Gi[z]-value;
			if (x > LOGMAXDOUBLE-1 || value == LOGZERO) {		
				value = Gi[z];
			} else {
				value += log1p(exp(x));
			}
		}
	}
	if (bound2 == M-2 && Gi[M-1] > LOGZERO) {
		double x = Gi[M-1]-log(2.0)-value;
		if (x > LOGMAXDOUBLE-1 || value == LOGZERO) {
			value = Gi[M-1]-log(2.0);
		} else {
			value += log1p(exp(x));
		}
	} else if (M > 3) {
		if (Gi[M-1] > LOGZERO) {
			double x = Gi[M-1]-value;
			if (x > LOGMAXDOUBLE-1 || value == LOGZERO) {
				value = Gi[M-1];
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
Lat1DFlat1stOS::NormPhiFree(Vector phi, const double norm) const {
	double logNorm = log(norm);
	for (int z=1; z<=M; z++) {
		if (phi[z] > LOGZERO) {
			phi[z] += logNorm;
		}
	}
}
void
Lat1DFlat1stOS::NormPhiRestr(Vector phi, 
						const Vector Gi, 
						double norm) const {
	int z;
	double logNorm = log(norm)-ComputeLnGN(Gi);
	for (z=1; z<=M; z++) {
		if (phi[z] > LOGZERO) {
			phi[z] += logNorm;
		}
	}
}

