#include "Lat2DFlat1stOSten.h"

Lat2DFlat1stOSten::Lat2DFlat1stOSten(Input* MyInput_, Text name_)
	: Lat2DFlat(MyInput_,name_) {
}

Lat2DFlat1stOSten::~Lat2DFlat1stOSten() {
}

void
Lat2DFlat1stOSten::GetLatticeInfo(int* Info) const{ 
}


void 
Lat2DFlat1stOSten::PropagateG(Matrix Gi, const Vector G, const int s) const {
	int x,y,z,zc,zi,s_1;
	double a,b,c,d,e,f,g,h,i;

	zc = 1;
	z = MY+1;
	zi = MY+MY+1;
	s_1 = s-1;
	for (x=2; x<MX; x++) {
		b = Gi[zc++][s_1];
		c = Gi[zc++][s_1];
		e = Gi[z++][s_1];
		f = Gi[z][s_1];
		h = Gi[zi++][s_1];
		i = Gi[zi++][s_1];
		for (y=2; y<MY; y++) {
			a=b; b=c; c=Gi[zc][s_1];
			d=e; e=f; f=Gi[z+1][s_1];
			g=h; h=i; i=Gi[zi][s_1];
			Gi[z][s] = G[z] * (l0*e + l1*(b+d+f+h) + l2*(a+c+g+i));
			zc++;
			z++;
			zi++;
		}
		z++;
	}
	UpdateBoundaries(Gi,s);
}

void 
Lat2DFlat1stOSten::PropagateG(Matrix Gi,const Vector G,const int s,const double f) const {
	Message(fatal,MyInput,"Force ensemble not implemented in Lat2DFlat1stOSten");
}


void 
Lat2DFlat1stOSten::PropagateG(const Vector Gi, const Vector G, Vector Gout) const {
	int x,y,z,zc,zi;
	double a,b,c,d,e,f,g,h,i;
	Vector GiPrime(2,MY);

	for (y=2; y<=MY; y++) {
		GiPrime[y] = Gi[y];
	}
	zc = 2;
	z = MY+1;
	zi = MY+MY+1;
	for (x=2; x<MX; x++) {
		b = Gi[z-MY];
		c = GiPrime[zc];
		e = Gi[z++];
		f = Gi[z];
		h = Gi[zi++];
		i = Gi[zi++];
		GiPrime[zc++] = f;
		for (y=2; y<MY; y++) {
			a=b; b=c; c=GiPrime[zc];
			d=e; e=f; f=Gi[z+1];
			g=h; h=i; i=Gi[zi];
			Gout[z] = G[z] * (l0*e + l1*(b+d+f+h) + l2*(a+c+g+i));
			GiPrime[zc]=f;
			zc++;
			z++;
			zi++;
		}
		z++;
		zc-=MY-1;
	}
	UpdateBoundaries(Gout);
}

void 
Lat2DFlat1stOSten::PropagateG(const Vector Gi, const Vector G, Vector Gout,const double f) const {
	Message(fatal,MyInput,"Force ensemble not implemented in Lat2DFlat1stOSten");
}

void 
Lat2DFlat1stOSten::PropagateG(Vector Gi, const Vector G) const {
	int x,y,z,zc,zi;
	double a,b,c,d,e,f,g,h,i;
	Vector GiPrime(2,MY);

	for (y=2; y<=MY; y++) {
		GiPrime[y] = Gi[y];
	}
	zc = 2;
	z = MY+1;
	zi = MY+MY+1;
	for (x=2; x<MX; x++) {
		b = Gi[z-MY];
		c = GiPrime[zc];
		e = Gi[z++];
		f = Gi[z];
		h = Gi[zi++];
		i = Gi[zi++];
		GiPrime[zc++] = f;
		for (y=2; y<MY; y++) {
			a=b; b=c; c=GiPrime[zc];
			d=e; e=f; f=Gi[z+1];
			g=h; h=i; i=Gi[zi];
			Gi[z] = G[z] * (l0*e + l1*(b+d+f+h) + l2*(a+c+g+i));
			GiPrime[zc]=f;
			zc++;
			z++;
			zi++;
		}
		z++;
		zc-=MY-1;
	}
	UpdateBoundaries(Gi);
}

void 
Lat2DFlat1stOSten::PropagateG(Vector Gi, const Vector G, const double f) const {
	Message(fatal,MyInput,"Force ensemble not implemented in Lat2DFlat1stOSten");
}

void
Lat2DFlat1stOSten::Propagate2G(Matrix Gi1, 
						   Matrix Gi2, 
						   const Vector G, const int s,
						   const LatticeRange* LatRange) const {
	int x,y,z,zc,zi,s_1;
	double a1,a2,b1,b2,c1,c2,d1,d2,e1,e2,f1,f2,g1,g2,h1,h2,i1,i2;

	zc = 1;
	z = MY+1;
	zi = MY+MY+1;
	s_1 = s-1;
	for (x=2; x<MX; x++) {
		b1 = Gi1[zc][s_1];
		b2 = Gi2[zc++][s_1];
		c1 = Gi1[zc][s_1];
		c2 = Gi2[zc++][s_1];
		e1 = Gi1[z][s_1];
		e2 = Gi2[z++][s_1];
		f1 = Gi1[z][s_1];
		f2 = Gi2[z][s_1];
		h1 = Gi1[zi][s_1];
		h2 = Gi2[zi++][s_1];
		i1 = Gi1[zi][s_1];
		i2 = Gi2[zi++][s_1];
		for (y=2; y<MY; y++) {
			a1=b1; b1=c1; c1=Gi1[zc][s_1];
			a2=b2; b2=c2; c2=Gi2[zc][s_1];
			d1=e1; e1=f1; f1=Gi1[z+1][s_1];
			d2=e2; e2=f2; f2=Gi2[z+1][s_1];
			g1=h1; h1=i1; i1=Gi1[zi][s_1];
			g2=h2; h2=i2; i2=Gi2[zi][s_1];
			if (LatRange->InRange(z)) {
				Gi1[z][s] = G[z] * (l0*e1 // e2 = 0
									+ l1*(b1+d1+f1+h1+b2+d2+f2+h2)
									+ l2*(a1+c1+g1+i1+a2+c2+g2+i2));
				Gi2[z][s] = 0;
			} else {
				Gi1[z][s] = G[z] * (l0*e1 + l1*(b1+d1+f1+h1) + l2*(a1+c1+g1+i1));
				Gi2[z][s] = G[z] * (l0*e2 + l1*(b2+d2+f2+h2) + l2*(a2+c2+g2+i2));
			}
			zc++;
			z++;
			zi++;
		}
		z++;
	}
	UpdateBoundaries(Gi1,s);
	UpdateBoundaries(Gi2,s);
}

void
Lat2DFlat1stOSten::Propagate2G(Matrix Gi1, 
						   Matrix Gi2, 
						   const Vector G, const int s,
						   const LatticeRange* LatRange,
						   const double f) const {
	Message(fatal,MyInput,"Force ensemble not implemented in Lat2DFlat1stOSten");
}

void
Lat2DFlat1stOSten::Propagate2G(Vector Gi1, 
						   Vector Gi2, 
						   const Vector G, 
						   const LatticeRange* LatRange) const {
	int x,y,z,zc,zi;
	double a1,a2,b1,b2,c1,c2,d1,d2,e1,e2,f1,f2,g1,g2,h1,h2,i1,i2;
	Vector GiPrime1(2,MY), GiPrime2(2,MY);

	for (y=2; y<=MY; y++) {
		GiPrime1[y] = Gi1[y];
		GiPrime2[y] = Gi2[y];
	}
	zc = 2;
	z = MY+1;
	zi = MY+MY+1;
	for (x=2; x<MX; x++) {
		b1 = Gi1[z-MY];
		b2 = Gi2[z-MY];
		c1 = GiPrime1[zc];
		c2 = GiPrime2[zc];
		e1 = Gi1[z];
		e2 = Gi2[z++];
		f1 = Gi1[z];
		f2 = Gi2[z];
		h1 = Gi1[zi];
		h2 = Gi2[zi++];
		i1 = Gi1[zi];
		i2 = Gi2[zi++];
		GiPrime1[zc] = f1;
		GiPrime2[zc++] = f2;
		for (y=2; y<MY; y++) {
			a1=b1; b1=c1; c1=GiPrime1[zc];
			a2=b2; b2=c2; c2=GiPrime2[zc];
			d1=e1; e1=f1; f1=Gi1[z+1];
			d2=e2; e2=f2; f2=Gi2[z+1];
			g1=h1; h1=i1; i1=Gi1[zi];
			g2=h2; h2=i2; i2=Gi2[zi];
			if (LatRange->InRange(z)) {
				Gi1[z] = G[z] * (l0*e1
								+ l1*(b1+d1+f1+h1+b2+d2+f2+h2)
								+ l2*(a1+c1+g1+i1+a2+c2+g2+i2));
				Gi2[z] = 0;
			} else {
				Gi1[z] = G[z] * (l0*e1 + l1*(b1+d1+f1+h1) + l2*(a1+c1+g1+i1));
				Gi2[z] = G[z] * (l0*e2 + l1*(b2+d2+f2+h2) + l2*(a2+c2+g2+i2));
			}
			GiPrime1[zc]=f1;
			GiPrime2[zc]=f2;
			zc++;
			z++;
			zi++;
		}
		z++;
		zc-=MY-1;
	}
	UpdateBoundaries(Gi1);
	UpdateBoundaries(Gi2);
}

void
Lat2DFlat1stOSten::Propagate2G(Vector Gi1, 
						   Vector Gi2, 
						   const Vector G, 
						   const LatticeRange* LatRange,
						   const double f) const {
	Message(fatal,MyInput,"Force ensemble not implemented in Lat2DFlat1stOSten");
}

void
Lat2DFlat1stOSten::UpdateBoundaries(Vector A) const {
	SetBoundaries(A);
}

void
Lat2DFlat1stOSten::UpdateBoundaries(Matrix A, const int s) const {
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
