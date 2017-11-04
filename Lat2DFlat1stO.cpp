#include "Lat2DFlat1stO.h"

// static Vectors for the Propagate functions. This avoid allocating memory
// at each call to a Propagate function, which is very time intensive
Vector Lat2DFlat1stO::Gi1Prime;
Vector Lat2DFlat1stO::Gi2Prime;

Lat2DFlat1stO::Lat2DFlat1stO(Input* MyInput_, Text name_)
	: Lat2DFlat(MyInput_,name_) {
	Gi1Prime.Dim(1,MY);
	Gi2Prime.Dim(1,MY);
}
Lat2DFlat1stO::~Lat2DFlat1stO() {
}
Boolean 
Lat2DFlat1stO::OverflowProtection(void) const {
	return false;
}
void 
Lat2DFlat1stO::MakeSafe(Vector) const {
}

void
Lat2DFlat1stO::GetLatticeInfo(int* Info) const{ 
}


void
Lat2DFlat1stO::RestoreFromSafe(Vector) const {
}
void 
Lat2DFlat1stO::PropagateG(Matrix Gi, const Vector G, const int s) const {
	int x,y,x_1,x1,z;
	double a,b,c;
	int s_1 = s-1;
	x_1 = 2;
	z = MY+2;
	x1 = MY*2 + 2;
	for (x=2; x<MX; x++) {
		b = Gi[z-1][s_1];
		c = Gi[z][s_1];
		for (y=2; y<MY; y++) {
			a = b; b = c; c = Gi[z+1][s_1];
			Gi[z][s] = G[z] * (l0*b + l1y*(a + c )+ l1x*(Gi[x_1++][s_1] + Gi[x1++][s_1]));
			z++;
		}
		x1+=2;
		x_1+=2;
		z+=2;
	}
	UpdateBoundaries(Gi,s);
}

void 
Lat2DFlat1stO::PropagateG(Matrix Gi, const Vector G, const int s,const double f) const {
	Message(fatal,MyInput,"Force ensemble not implemented in Lat2DFlat1stO");
}

void 
Lat2DFlat1stO::PropagateG(Vector Gi, const Vector G) const {
	int x,y,x1,z;
	double a,b,c;

	z = MY+2;
	x1 = MY*2 + 2;
	for (y=2; y<MY; y++) {
		Gi1Prime[y] = Gi[y];
	}
	for (x=2; x<MX; x++) {
		b = Gi[z-1];
		c = Gi[z];
		for (y=2; y<MY; y++) {
			a = b; b = c; c = Gi[z+1];
			Gi[z] = G[z] * (l0*b + l1y*(a + c) + l1x*(Gi1Prime[y] + Gi[x1++]));
			Gi1Prime[y] = b;
			z++;
		}
		x1+=2;
		z+=2;
	}
	UpdateBoundaries(Gi);
}

void 
Lat2DFlat1stO::PropagateG(Vector Gi, const Vector G,const double f) const {
	Message(fatal,MyInput,"Force ensemble not implemented in Lat2DFlat1stO");
}

void 
Lat2DFlat1stO::PropagateG(const Vector Gi, const Vector G, Vector Gout) const {
	int x,y,x1,z;
	double a,b,c;

	z = MY+2;
	x1 = MY*2 + 2;
	for (y=2; y<MY; y++) {
		Gi1Prime[y] = Gi[y];
	}
	for (x=2; x<MX; x++) {
		b = Gi[z-1];
		c = Gi[z];
		for (y=2; y<MY; y++) {
			a = b; b = c; c = Gi[z+1];
			Gout[z] = G[z] * (l0*b + l1y*(a + c) + l1x*(Gi1Prime[y] + Gi[x1++]));
			Gi1Prime[y] = b;
			z++;
		}
		x1+=2;
		z+=2;
	}
	UpdateBoundaries(Gout);
}

void 
Lat2DFlat1stO::PropagateG(const Vector Gi, const Vector G, Vector Gout,const double f) const {
	Message(fatal,MyInput,"Force ensemble not implemented in Lat2DFlat1stO");
}

void
Lat2DFlat1stO::Init2G(Vector Gi1, 
					  Vector Gi2, 
					  const Vector G, 
					  const LatticeRange* LatRange) const {
	int x,y,z;

	z = MY+2;
	for (x=2; x<MX; x++) {
		for (y=2; y<MY; y++) {
			if (LatRange->InRange(z)) {
				Gi1[z] = G[z];
				Gi2[z] = 0;
			} else {
				Gi1[z] = 0;
				Gi2[z] = G[z];
			}
			z++;
		}
		z+=2;
	}
	UpdateBoundaries(Gi1);
	UpdateBoundaries(Gi2);
}
void
Lat2DFlat1stO::Propagate2G(Matrix Gi1, 
						   Matrix Gi2, 
						   const Vector G, const int s,
						   const LatticeRange* LatRange) const {
	int x,y,x_1,x1,z;
	double a1,b1,c1,a2,b2,c2;
	int s_1 = s-1;
	x_1 = 2;
	z = MY+2;
	x1 = MY*2 + 2;
	for (x=2; x<MX; x++) {
		b1 = Gi1[z-1][s_1];
		b2 = Gi2[z-1][s_1];
		c1 = Gi1[z][s_1];
		c2 = Gi2[z][s_1];
		for (y=2; y<MY; y++) {
			a1 = b1; b1 = c1; c1 = Gi1[z+1][s_1];
			a2 = b2; b2 = c2; c2 = Gi2[z+1][s_1];
			if (LatRange->InRange(z)) {
				Gi1[z][s] = G[z] * (l0*b1 + l1y*(a1 + c1 + a2 + c2) + l1x*(Gi1[x_1][s_1] + Gi1[x1][s_1] + Gi2[x_1][s_1] + Gi2[x1][s_1]));
				Gi2[z][s] = 0;
			} else {
				Gi1[z][s] = G[z] * (l0*b1 + l1y*(a1 + c1) + l1x*(Gi1[x_1][s_1] + Gi1[x1][s_1]));
				Gi2[z][s] = G[z] * (l0*b2 + l1y*(a2 + c2) + l1x*(Gi2[x_1][s_1] + Gi2[x1][s_1]));
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
Lat2DFlat1stO::Propagate2G(Matrix Gi1, 
						   Matrix Gi2, 
						   const Vector G, const int s,
						   const LatticeRange* LatRange, 
						   const double f) const {
	Message(fatal,MyInput,"Force ensemble not implemented in Lat2DFlat1stO");
}


void
Lat2DFlat1stO::Propagate2G(Vector Gi1, 
						   Vector Gi2, 
						   const Vector G, 
						   const LatticeRange* LatRange) const {
	int x,y,x1,z;
	double a1,b1,c1,a2,b2,c2;

	for (y=2; y<MY; y++) {
		Gi1Prime[y] = Gi1[y];
		Gi2Prime[y] = Gi2[y];
	}
	z = MY+2;
	x1 = MY*2 + 2;
	for (x=2; x<MX; x++) {
		b1 = Gi1[z-1];
		b2 = Gi2[z-1];
		c1 = Gi1[z];
		c2 = Gi2[z];
		for (y=2; y<MY; y++) {
			a1 = b1; b1 = c1; c1 = Gi1[z+1];
			a2 = b2; b2 = c2; c2 = Gi2[z+1];
			if (LatRange->InRange(z)) {
				Gi1[z] = G[z] * (l0*b1 + l1y*(a1 + c1 + a2 + c2)+ l1x*(Gi1Prime[y] + Gi1[x1]  + Gi2Prime[y] + Gi2[x1]));
				Gi2[z] = 0;
			} else {
				Gi1[z] = G[z] * (l0*b1 + l1y*(a1 + c1) + l1x*(Gi1Prime[y] + Gi1[x1]));
				Gi2[z] = G[z] * (l0*b2 + l1y*(a2 + c2) + l1x*(Gi2Prime[y] + Gi2[x1]));
			}
			Gi1Prime[y] = b1;
			Gi2Prime[y] = b2;
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
Lat2DFlat1stO::Propagate2G(Vector Gi1, 
						   Vector Gi2, 
						   const Vector G, 
						   const LatticeRange* LatRange,
						   const double f) const {
	Message(fatal,MyInput,"Force ensemble not implemented in Lat2DFlat1stO");
}

void
Lat2DFlat1stO::ConnectG(const Vector GiA, 
						const Matrix GiB, 
						const int s, 
						Vector out) const {
	int z;
	for (z=1; z<=MX*MY; z++) {
		out[z] += GiA[z]*GiB[z][s];
	}
}
void
Lat2DFlat1stO::ConnectG(const Vector GiA, 
						const Vector GiB, 
						Vector out) const {
	int z;
	for (z=1; z<=MX*MY; z++) {
		out[z] += GiA[z]*GiB[z];
	}
}
void 
Lat2DFlat1stO::Connect2G(const Vector GiA1, 
						 const Matrix GiB1, 
						 const int s1, 
						 const Vector GiA2, 
						 const Matrix GiB2, 
						 const int s2, 
						 Vector out) const {
	int z;
	for (z=1; z<=MX*MY; z++) {
		out[z] += GiA1[z]*GiB2[z][s1] + GiA2[z]*GiB1[z][s2];
	}
}
void 
Lat2DFlat1stO::Connect2G(const Vector GiA1, 
						 const Vector GiB1, 
						 const Vector GiA2, 
						 const Vector GiB2, 
						 Vector out) const {
	int z;
	for (z=1; z<=MX*MY; z++) {
		out[z] += GiA1[z]*GiB2[z] + GiA2[z]*GiB1[z];
	}
}
void 
Lat2DFlat1stO::CorrectDoubleCountG(Vector in, const Vector G) const {
	int z;
	for (z=1; z<=MX*MY; z++) {
		if (G[z] > 0) in[z] /= G[z];
		else in[z] = 0;
	}
}
double
Lat2DFlat1stO::ComputeLnGN(const Vector Gi) const {
	int x,y,z;
	double value = 0;
	z = MY+2;
	for (x=2; x<MX; x++) {
		for (y=2; y<MY; y++) {
			value += Gi[z++];
		}
		z+=2;
	}
	if (boundX1 == 3) {
		z = MY+2;
		for (y=2; y<MY; y++) {
			value -= 0.5*Gi[z++];
		}
	}
	if (boundX2 == MX-2) {
		z = MY*(MX-2)+2;
		for (y=2; y<MY; y++) {
			value -= 0.5*Gi[z++];
		}
	}
	if (boundY1 == 3) {
		z = MY+2;
		for (x=2; x<MX; x++) {
			value -= 0.5*Gi[z];
			z+=MY;
		}
	}
	if (boundY2 == MY-2) {
		z = 2*MY-2;
		for (x=2; x<MX; x++) {
			value -= 0.5*Gi[z];
			z+=MY;
		}
	}
	if (boundX1 == 3 && boundY1 == 3) {
		value += 0.25*Gi[MY+2];
	}
	if (boundX1 == 3 && boundY2 == MY-2) {
		value += 0.25*Gi[2*MY-2];
	}
	if (boundX2 == MX-2 && boundY1 == 3) {
		value += 0.25*Gi[MY*(MX-2)+2];
	}
	if (boundX2 == MX-2 && boundY2 == MY-2) {
		value += 0.25*Gi[(MX-1)*MY-2];
	}
	return log(value);
}
void 
Lat2DFlat1stO::NormPhiFree(Vector phi, const double norm) const {
	int z;
	for (z=1; z<=MX*MY; z++) {
		phi[z] *= norm;
	}
}
void 
Lat2DFlat1stO::NormPhiRestr(Vector phi, 
					   const Vector Gi, 
					   double norm) const {
	int z;
	norm /= exp(ComputeLnGN(Gi));
	for (z=1; z<=MX*MY; z++) {
		phi[z] *= norm;
	}
}
void
Lat2DFlat1stO::UpdateBoundaries(Vector A) const {
	SetBoundaries(A);
}
void
Lat2DFlat1stO::UpdateBoundaries(Matrix A, const int s) const {
	int z;
	// set boundaries without corners
	// set lower boundary
	if (boundX1 != 1) {
		for (z=2; z<MY; z++) {
			A[z][s] = A[(boundX1-1)*MY+z][s];
		}
	}
	// set upper boundary
	if (boundX2 != MX) {
		for (z=2; z<MY; z++) {
			A[M-MY+z][s] = A[(boundX2-1)*MY+z][s];
		}
	}
	// set left boundary
	if (boundY1 != 1) {
		for (z=MY; z<M-MY; z+=MY) {
			A[z+1][s] = A[z+boundY1][s];
		}
	}
	// set right boundary
	if (boundY2 != MY) {
		for (z=MY+MY; z<M; z+=MY) {
			A[z][s] = A[z-MY+boundY2][s];
		}
	}

	// Set corners.
	// If one boundary is periodic or mirror, it's easy.
	// Else the corner is defined already.
	if (boundY1 != 1 && !bulkBoundY1) {
		A[1][s] = A[boundY1][s];
		A[M-MY+1][s] = A[M-MY+boundY1][s];
	}
	if (boundY2 != MY && !bulkBoundY2) {
		A[MY][s] = A[boundY2][s];
		A[M][s] = A[M-MY+boundY2][s];
	}
	if (boundX1 != 1 && !bulkBoundX1) {
		A[1][s] = A[(boundX1-1)*MY+1][s];
		A[MY][s] = A[boundX1*MY][s];
	}
	if (boundX2 != MX && !bulkBoundX2) {
		A[M-MY+1][s] = A[(boundX2-1)*MY+1][s];
		A[M][s] = A[boundX2*MY][s];
	}
}
