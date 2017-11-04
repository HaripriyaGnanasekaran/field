#include "Lat2DFlat1stOS.h"

// static Vectors for the Propagate functions. This avoid allocating memory
// at each call to a Propagate function, which is very time intensive
Vector Lat2DFlat1stOS::Gi1Prime;
Vector Lat2DFlat1stOS::Gi2Prime;

Lat2DFlat1stOS::Lat2DFlat1stOS(Input* MyInput_, Text name_)
	: Lat2DFlat(MyInput_,name_) {
	Log_L1L0 = log(l1/(1-4*l1));
	Log_L0 = log(1-4*l1);
	Log_L1 = log(l1);
	Gi1Prime.Dim(1,MY);
	Gi2Prime.Dim(1,MY);
}
Lat2DFlat1stOS::~Lat2DFlat1stOS() {
}
Boolean 
Lat2DFlat1stOS::OverflowProtection(void) const {
	return true;
}

void
Lat2DFlat1stOS::GetLatticeInfo(int* Info) const{ 
}

void 
Lat2DFlat1stOS::MakeSafe(Vector A) const {
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
Lat2DFlat1stOS::RestoreFromSafe(Vector A) const {
	int z;

	for (z=1; z<=M; z++) {
		if (A[z] > LOGZERO) {
			if (A[z] > LOGCUTOFF) {
				A[z] = exp(LOGCUTOFF+LOGMAXPREFACTOR
				*tanh((A[z]-LOGCUTOFF)/LOGMAXPREFACTOR));
			} else {
				A[z] = exp(A[z]);
			}
		} else {
			A[z] = 0;
		}
	}
}

void 
Lat2DFlat1stOS::PropagateG(Matrix Gi, const Vector G, const int s) const {
	int x, y, x_1, x1, z;
	double b;
	double n[4];
	int s_1;

	s_1 = s-1;
	x_1 = 2;
	z = MY+2;
	x1 = MY*2 + 2;
	for (x=2; x<MX; x++) {
		b = Gi[z-1][s_1];
		n[1] = Gi[z][s_1];
		for (y=2; y<MY; y++) {
			n[0] = b; b = n[1]; n[1] = Gi[z+1][s_1];
			if (G[z] == LOGZERO) {
				Gi[z][s] = LOGZERO;
				z++;
				x_1++;
				x1++;
				continue;
			}
			n[2] = Gi[x_1++][s_1]; n[3] = Gi[x1++][s_1];
			Gi[z][s] = getNewGi(G[z], b, n, 4);
			z++;
		}
		x1+=2;
		x_1+=2;
		z+=2;
	}
	UpdateBoundaries(Gi,s);
}
void
Lat2DFlat1stOS::PropagateG(Matrix Gi, const Vector G, const int s, const double f) const {
	Message(fatal,MyInput,"Force ensemble not implemented in Lat2DFlat1stOS");
}

void 
Lat2DFlat1stOS::PropagateG(Vector Gi, const Vector G) const {
	int x,y,x1,z;
	double b;
	double n[4];

	z = MY+2;
	x1 = MY*2 + 2;
	for (y=2; y<MY; y++) {
		Gi1Prime[y] = Gi[y];
	}
	for (x=2; x<MX; x++) {
		b = Gi[z-1];
		n[1] = Gi[z];
		for (y=2; y<MY; y++) {
			n[0] = b; b = n[1]; n[1] = Gi[z+1];
			n[2] = Gi1Prime[y]; n[3] = Gi[x1++];
			Gi1Prime[y] = b;
			if (G[z] == LOGZERO) {
				Gi[z] = LOGZERO;
				z++;
				continue;
			}
			Gi[z] = getNewGi(G[z], b, n, 4);
			z++;
		}
		x1+=2;
		z+=2;
	}
	UpdateBoundaries(Gi);
}

void 
Lat2DFlat1stOS::PropagateG(Vector Gi, const Vector G,const double f) const {
	Message(fatal,MyInput,"Force ensemble not implemented in Lat2DFlat1stOS");
}

void 
Lat2DFlat1stOS::PropagateG(const Vector Gi, const Vector G, Vector Gout) const {
	int x,y,x1,z;
	double b;
	double n[4];

	z = MY+2;
	x1 = MY*2 + 2;
	for (y=2; y<MY; y++) {
		Gi1Prime[y] = Gi[y];
	}
	for (x=2; x<MX; x++) {
		b = Gi[z-1];
		n[1] = Gi[z];
		for (y=2; y<MY; y++) {
			n[0] = b; b = n[1]; n[1] = Gi[z+1];
			n[2] = Gi1Prime[y]; n[3] = Gi[x1++];
			Gi1Prime[y] = b;
			if (G[z] == LOGZERO) {
				Gout[z] = LOGZERO;
				z++;
				continue;
			}
			Gout[z] = getNewGi(G[z], b, n, 4);
			z++;
		}
		x1+=2;
		z+=2;
	}
	UpdateBoundaries(Gout);
}

void 
Lat2DFlat1stOS::PropagateG(const Vector Gi, const Vector G, Vector Gout,const double f) const {
	Message(fatal,MyInput,"Force ensemble not implemented in Lat2DFlat1stOS");
}

void 
Lat2DFlat1stOS::Init2G(Vector Gi1, Vector Gi2, const Vector G, const LatticeRange* LatRange) const {
	int x,y,z;

	z = MY+2;
	for (x=2; x<MX; x++) {
		for (y=2; y<MY; y++) {
			if (LatRange->InRange(z)) {
				Gi1[z] = G[z];
				Gi2[z] = LOGZERO;
			} else {
				Gi1[z] = LOGZERO;
				Gi2[z] = G[z];
			}
			z++;
		}
		z+=2;
	}
	// define boundaries where no mirror is
	if (boundX1 == 1) {
		for (z=1; z<=MY; z++) {
			Gi1[z] = Gi2[z] = LOGZERO;
		}
	}
	if (boundX2 == MX) {
		for (z=M-MY+1; z<=M; z++) {
			Gi1[z] = Gi2[z] = LOGZERO;
		}
	}
	if (boundY1 == 1) {
		for (z=1; z<M; z+=MY) {
			Gi1[z] = Gi2[z] = LOGZERO;
		}
	}
	if (boundY2 == MY) {
		for (z=MY; z<=M; z+=MY) {
			Gi1[z] = Gi2[z] = LOGZERO;
		}
	}
//	Gi1[1]      = Gi2[1]      = LOGZERO;
//	Gi1[MY]     = Gi2[MY]     = LOGZERO;
//	Gi1[M-MY+1] = Gi2[M-MY+1] = LOGZERO;
//	Gi1[M]      = Gi2[M]      = LOGZERO;
	UpdateBoundaries(Gi1);
	UpdateBoundaries(Gi2);
}

void
Lat2DFlat1stOS::Propagate2G(Matrix Gi1, Matrix Gi2, const Vector G, const int s,
				const LatticeRange* LatRange) const {
	int x, y, x_1, x1, z;
	int s_1;
	double b1,b2;
	double n[8];

	s_1 = s-1;
	x_1 = 2;
	z = MY+2;
	x1 = MY*2 + 2;
	for (x=2; x<MX; x++) {
		b1 = Gi1[z-1][s_1];
		b2 = Gi2[z-1][s_1];
		n[1] = Gi1[z][s_1];
		n[5] = Gi2[z][s_1];
		for (y=2; y<MY; y++) {
			n[0] = b1; b1 = n[1]; n[1] = Gi1[z+1][s_1];
			n[4] = b2; b2 = n[5]; n[5] = Gi2[z+1][s_1];
			if (G[z] == LOGZERO) {
				Gi1[z][s] = LOGZERO;
				Gi2[z][s] = LOGZERO;
				x_1++;
				x1++;
				z++;
				continue;
			}
			n[2] = Gi1[x_1][s_1]; n[3] = Gi1[x1][s_1];
			n[6] = Gi2[x_1][s_1]; n[7] = Gi2[x1][s_1];
			if (LatRange->InRange(z)) {
				Gi1[z][s] = getNewGi(G[z], b1, n, 8);
				Gi2[z][s] = LOGZERO;
			} else {
				Gi1[z][s] = getNewGi(G[z], b1, n, 4);
				Gi2[z][s] = getNewGi(G[z], b2, n+4, 4);
			}
			x_1++;
			x1++;
			z++;
		}
		x1+=2;
		x_1+=2;
		z+=2;
	}
	UpdateBoundaries(Gi1,s);
	UpdateBoundaries(Gi2,s);
}

void
Lat2DFlat1stOS::Propagate2G(Matrix Gi1, Matrix Gi2, const Vector G, const int s,
				const LatticeRange* LatRange,const double f) const {
	Message(fatal,MyInput,"Force ensemble not implemented in Lat2DFlat1stOS");
}

void
Lat2DFlat1stOS::Propagate2G(Vector Gi1, Vector Gi2, const Vector G,
			const LatticeRange* LatRange) const {
	int x,y,x1,z;
	double b1,b2;
	double n[8];

	for (y=2; y<MY; y++) {
		Gi1Prime[y] = Gi1[y];
		Gi2Prime[y] = Gi2[y];
	}
	z = MY+2;
	x1 = MY*2 + 2;
	for (x=2; x<MX; x++) {
		b1 = Gi1[z-1];
		b2 = Gi2[z-1];
		n[1] = Gi1[z];
		n[5] = Gi2[z];
		for (y=2; y<MY; y++) {
			n[0] = b1; b1 = n[1]; n[1] = Gi1[z+1];
			n[4] = b2; b2 = n[5]; n[5] = Gi2[z+1];
			if (G[z] == LOGZERO) {
				Gi1[z] = LOGZERO;
				Gi2[z] = LOGZERO;
				Gi1Prime[y] = b1;
				Gi2Prime[y] = b2;
				x1++;
				z++;
				continue;
			}
			n[2] = Gi1Prime[y]; n[3] = Gi1[x1];
			n[6] = Gi2Prime[y]; n[7] = Gi2[x1];
			Gi1Prime[y] = b1;
			Gi2Prime[y] = b2;
			if (LatRange->InRange(z)) {
				Gi1[z] = getNewGi(G[z], b1, n, 8);
				Gi2[z] = LOGZERO;
			} else {
				Gi1[z] = getNewGi(G[z], b1, n, 4);
				Gi2[z] = getNewGi(G[z], b2, n+4, 4);
			}
			x1++;
			z++;
		}
		x1+=2;
		z+=2;
	}
	UpdateBoundaries(Gi1);
	UpdateBoundaries(Gi2);
}

void
Lat2DFlat1stOS::Propagate2G(Vector Gi1, Vector Gi2, const Vector G,
			const LatticeRange* LatRange,const double f) const {
	Message(fatal,MyInput,"Force ensemble not implemented in Lat2DFlat1stOS");
}
double
Lat2DFlat1stOS::getNewGi(const double G, const double b,
			const double *n, short nn) const {
	short i,j;
	double lnterm;

	lnterm = 0;
	if (b != LOGZERO && l1 < 0.25) {
		for (i=0; i<nn; i++) {
			if (n[i] != LOGZERO) {
				lnterm += exp(Log_L1L0 + n[i] - b);
			}
		}
		return G + Log_L0 + b + log1p(lnterm);
	} else {
		i = 0;
		while (i < nn && n[i] == LOGZERO) i++;
		if (i == nn) {
			return LOGZERO;
		} else {
			for (j=i+1; j<nn; j++) {
				if (n[j] != LOGZERO) {
					lnterm += exp(n[j] - n[i]);
				}
			}
			return G + Log_L1 + n[i] + log1p(lnterm); 
		}
	}
}
void
Lat2DFlat1stOS::ConnectG(const Vector GiA, const Matrix GiB, const int s,
						 Vector out) const {
	int z;

	for (z=1; z<=M; z++) {
		addLntoLn(GiA[z], GiB[z][s], out[z], true);
	}
}
void
Lat2DFlat1stOS::ConnectG(const Vector GiA, const Vector GiB, Vector out) const {
	int z;

	for (z=1; z<=M; z++) {
		addLntoLn(GiA[z], GiB[z], out[z], true);
	}
}
void 
Lat2DFlat1stOS::Connect2G(const Vector GiA1, const Matrix GiB1, const int s1,
			  const Vector GiA2, const Matrix GiB2, const int s2,
						  Vector out) const {
	int z;

	for (z=1; z<=M; z++) {
		addLntoLn(GiA1[z], GiB2[z][s1], out[z], true);
		addLntoLn(GiA2[z], GiB1[z][s2], out[z], true);
	}
}
void 
Lat2DFlat1stOS::Connect2G(const Vector GiA1, const Vector GiB1,
			  const Vector GiA2, const Vector GiB2,
						  Vector out) const {
	int z;

	for (z=1; z<=M; z++) {
		addLntoLn(GiA1[z], GiB2[z], out[z], true);
		addLntoLn(GiA2[z], GiB1[z], out[z], true);
	}
}
void 
Lat2DFlat1stOS::CorrectDoubleCountG(Vector in, const Vector G) const {
	int z;

	for (z=1; z<=MX*MY; z++) {
		if (G[z] != LOGZERO && in[z] != LOGZERO) {
			in[z] -= G[z];
		} else {
			in[z] = LOGZERO;
		}
	}
}
/*
	We want to calculate: V += V1*V2,
	but we have only the logarithms of these values:
	v = ln(V), v1 = ln(V1), v2 = ln(V2)

	We can derive a new formula
	Vnew = Vold + V1*V2
	ln(Vnew) = ln(Vold + V1*V2) -> ln(Vnew) = ln(Vold) + ln(1+V1*V2/Vold)
	-> ln(Vnew) = ln(Vold) + ln( 1+exp(ln(V1) + ln(V2) - ln(Vold) ))
	
	->	vnew = vold + ln(1 + exp( v1 + v2 - vold )

	If we want to calculate V -= V1*V2, the formula becomes:
	->	vnew = vold + ln(1 + exp(-(v1 + v2 - vold))
*/
void
Lat2DFlat1stOS::addLntoLn(const double v1, const double v2, double &v,
		const bool add) const {
	double sum;

	if (v1 != LOGZERO && v2 != LOGZERO) {
		if (v == LOGZERO) {
			v = v1 + v2;
		} else {
			sum = v1 + v2 - v;
			if (sum > LOGMAXDOUBLE-1) {
				v = v1 + v2;
			} else {
				sum = (add) ?exp(sum) :-exp(sum);
				v += log1p(sum);
			}
		}
	}
}
double
Lat2DFlat1stOS::ComputeLnGN(const Vector Gi) const {
	int x,y,z;
	double value, log2, log4;

	value = LOGZERO;
	log2 = log(2.0);
	log4 = log(4.0);
	z = MY+2;
	for (x=2; x<MX; x++) {
		for (y=2; y<MY; y++) {
			addLntoLn(Gi[z++], 0, value, true);
		}
		z+=2;
	}
	if (boundX1 == 3) {
		z = MY+2;
		for (y=2; y<MY; y++) {
			addLntoLn(Gi[z++], -log2, value, false);
		}
	}		
	if (boundX2 == MX-2) {
		z = MY*(MX-2)+2;
		for (y=2; y<MY; y++) {
			addLntoLn(Gi[z++], -log2, value, false);
		}
	}
	if (boundY1 == 3) {
		z = MY+2;
		for (x=2; x<MX; x++) {
			addLntoLn(Gi[z], -log2, value, false);
			z+=MY;
		}
	}
	if (boundY2 == MY-2) {
		z = 2*MY-2;
		for (x=2; x<MX; x++) {
			addLntoLn(Gi[z], -log2, value, false);
			z+=MY;
		}
	}
	if (boundX1 == 3 && boundY1 == 3) {
		z = MY + 2;
		addLntoLn(Gi[z], -log4, value, true);
	}
	if (boundX1 == 3 && boundY2 == MY-2) {
		z = 2*MY-2;
		addLntoLn(Gi[z], -log4, value, true);
	}
	if (boundX2 == MX-2 && boundY1 == 3) {
		z = MY*(MX-2)+2;
		addLntoLn(Gi[z], -log4, value, true);
	}
	if (boundX2 == MX-2 && boundY2 == MY-2) {
		z =(MX-1)*MY-2;
		addLntoLn(Gi[z], -log4, value, true);
	}
	return value;
}
void 
Lat2DFlat1stOS::NormPhiFree(Vector phi, const double norm) const {
	int z;
	double logNorm;
	
	if (norm == 0) {
		for (z=1; z<=M; z++) {
			phi[z] = LOGZERO;
		}
		return;
	}
	logNorm = log(norm);
	for (z=1; z<=M; z++) {
		if (phi[z] > LOGZERO) {
			phi[z] += logNorm;
		}
	}
}
void 
Lat2DFlat1stOS::NormPhiRestr(Vector phi, const Vector Gi, double norm) const {
	int z;
	double logNorm;

	logNorm = log(norm)-ComputeLnGN(Gi);
	for (z=1; z<=MX*MY; z++) {
		if (phi[z] > LOGZERO) {
			phi[z] += logNorm;
		}
	}
}
void
Lat2DFlat1stOS::UpdateBoundaries(Vector A) const {
	int y,z;
	int i1, i2, i3;

	for (z=MY; z<M-MY; z+=MY) {
		A[z+1]  = LOGZERO;
		A[z+1]  = A[z+boundY1];
		A[z+MY] = LOGZERO;
		A[z+MY] = A[z+boundY2];
	}
	i1 = (boundX1-1)*MY;
	i2 = M-MY;
	i3 = (boundX2-1)*MY;
	for (y=2; y<MY; y++) {
		A[y] = LOGZERO;
		A[y] = A[i1+y];
		A[i2+y] = LOGZERO;
		A[i2+y] = A[i3+y];
	}
	// corner x=1,y=1,z=1
	A[1] = LOGZERO;
	A[1] = A[boundY1];
	A[1] = A[i1+1];
	// corner x=1,y=MY,z=MY
	A[MY] = LOGZERO;
	A[MY] = A[boundY2];
	A[MY] = A[i1+MY];
	// corner x=MX,y=1,z=M-MY+1=i2+1
	A[i2+1] = LOGZERO;
	A[i2+1] = A[i2+boundY1];
	A[i2+1] = A[i3+1];
	// corner x=MX,y=MY,z=M
	A[M] = LOGZERO;
	A[M] = A[i2+boundY2];
	A[M] = A[i3+MY];
}
void
Lat2DFlat1stOS::UpdateBoundaries(Matrix A, const int s) const {
	int y,z;
	int i1, i2, i3;

	for (z=MY; z<M-MY; z+=MY) {
		A[z+1][s]  = LOGZERO;
		A[z+1][s]  = A[z+boundY1][s];
		A[z+MY][s] = LOGZERO;
		A[z+MY][s] = A[z+boundY2][s];
	}
	i1 = (boundX1-1)*MY;
	i2 = M-MY;
	i3 = (boundX2-1)*MY;
	for (y=2; y<MY; y++) {
		A[y][s] = LOGZERO;
		A[y][s] = A[i1+y][s];
		A[i2+y][s] = LOGZERO;
		A[i2+y][s] = A[i3+y][s];
	}
	// corner x=1,y=1,z=1
	A[1][s] = LOGZERO;
	A[1][s] = A[boundY1][s];
	A[1][s] = A[i1+1][s];
	// corner x=1,y=MY,z=MY
	A[MY][s] = LOGZERO;
	A[MY][s] = A[boundY2][s];
	A[MY][s] = A[i1+MY][s];
	// corner x=MX,y=1,z=M-MY+1=i2+1
	A[i2+1][s] = LOGZERO;
	A[i2+1][s] = A[i2+boundY1][s];
	A[i2+1][s] = A[i3+1][s];
	// corner x=MX,y=MY,z=M
	A[M][s] = LOGZERO;
	A[M][s] = A[i2+boundY2][s];
	A[M][s] = A[i3+MY][s];
}
