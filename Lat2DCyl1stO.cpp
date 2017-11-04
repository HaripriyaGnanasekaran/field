#include "Lat2DCyl1stO.h"

Lat2DCyl1stO::Lat2DCyl1stO(Input* MyInput_, Text name_)
	: Lat2DFlat(MyInput_,name_), Lat2DCylinder(MyInput_,name_) {
	Gi1Prime.Dim(1,MY);
	Gi2Prime.Dim(1,MY);
}
Lat2DCyl1stO::~Lat2DCyl1stO() {
}

void
Lat2DCyl1stO::GetLatticeInfo(int* Info) const{ 
}


void 
Lat2DCyl1stO::PropagateG(Matrix Gi, const Vector G, const int s) const {
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
			Gi[z][s] = G[z] * (b + lambda1_1[x]*(Gi[x_1++][s_1]-b) + 
					                lambda11[x]*(Gi[x1++][s_1]-b) + 
					                l1*(a + c - b - b));
			z++;
		}
		x_1+=2;
		x1+=2;
		z+=2;
	}
	UpdateBoundaries(Gi,s);
}
// FIXME ppp - I am not sure if it is correct

void
Lat2DCyl1stO::PropagateG(Matrix Gi, const Vector G, const int s, const double f) const {
	int x,y,x_1,x1,z;
	double a,b,c;
	
	double f_1= exp(-f);
	double f1= exp(f);
	double norm = l1*(f_1+f1-2)+1;
	int s_1 = s-1;
	x_1 = 2;
	z = MY+2;   
	x1 = MY*2 + 2;
	for (x=2; x<MX; x++) {
		b = Gi[z-1][s_1];
		c = Gi[z][s_1];
		for (y=2; y<MY; y++) {
			a = b; b = c; c = Gi[z+1][s_1];
			Gi[z][s] = G[z] * (b + lambda1_1[x]*(Gi[x_1++][s_1]-b) + 
					                lambda11[x]*(Gi[x1++][s_1]-b) + 
					                l1*(a*f1 + c*f_1 - b - b))/norm;
			z++;
		}
		x_1+=2;
		x1+=2;
		z+=2;
	}
	UpdateBoundaries(Gi,s);
}


void
Lat2DCyl1stO::PropagateG(const Vector Gi, const Vector G, Vector Gout) const {
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
			Gout[z] = G[z] * (b + lambda1_1[x]*(Gi1Prime[y]-b) + 
					              lambda11[x]*(Gi[x1++]-b) + 
					              l1*(a + c - b - b));
			Gi1Prime[y] = b;
			z++;
		}
		x1+=2;
		z+=2;
	}
	UpdateBoundaries(Gout);
}

void
Lat2DCyl1stO::PropagateG(const Vector Gi, const Vector G, Vector Gout, const double f) const {
	int x,y,x1,z;
	double a,b,c;
	
	double f_1= exp(-f);
	double f1= exp(f);
	double norm = l1*(f_1+f1-2)+1;
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
			Gout[z] = G[z] * (b + lambda1_1[x]*(Gi1Prime[y]-b) + 
					              lambda11[x]*(Gi[x1++]-b) + 
					              l1*(a*f1 + c*f_1 - b - b))/norm;
			Gi1Prime[y] = b;
			z++;
		}
		x1+=2;
		z+=2;
	}
	UpdateBoundaries(Gout);
}



void
Lat2DCyl1stO::PropagateG(Vector Gi, const Vector G) const {
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
			Gi[z] = G[z] * (b + lambda1_1[x]*(Gi1Prime[y]-b) + 
					            lambda11[x]*(Gi[x1++]-b) + 
					            l1*(a + c - b - b));
			Gi1Prime[y] = b;
			z++;
		}
		x1+=2;
		z+=2;
	}
	UpdateBoundaries(Gi);
}

void
Lat2DCyl1stO::PropagateG(Vector Gi, const Vector G, const double f) const {
	int x,y,x1,z;
	double a,b,c;

	double f_1= exp(-f);
	double f1= exp(f);
	double norm = l1*(f_1+f1-2)+1;
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
			Gi[z] = G[z] * (b + lambda1_1[x]*(Gi1Prime[y]-b) + 
					            lambda11[x]*(Gi[x1++]-b) + 
					            l1*(a*f1 + c*f_1 - b - b))/norm;
			Gi1Prime[y] = b;
			z++;
		}
		x1+=2;
		z+=2;
	}
	UpdateBoundaries(Gi);
}




void
Lat2DCyl1stO::Propagate2G(Matrix Gi1, 
						  Matrix Gi2, 
						  const Vector G, 
						  const int s, 
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
				Gi1[z][s] = G[z] * (b1 + lambda1_1[x]*(Gi1[x_1][s_1] + Gi2[x_1][s_1] - b1) +
					                      lambda11[x]*(Gi1[x1][s_1] + Gi2[x1][s_1] - b1) + 
					                      l1*(a1 + c1 + a2 + c2 - b1 -b1));
				Gi2[z][s] = 0;
			} else {
				Gi1[z][s] = G[z] * (b1 + lambda1_1[x]*(Gi1[x_1][s_1]-b1) + 
						                 lambda11[x]*(Gi1[x1][s_1]-b1) + 
						                 l1*(a1+c1-b1-b1));
				Gi2[z][s] = G[z] * (b2 + lambda1_1[x]*(Gi2[x_1][s_1]-b2) + 
						                 lambda11[x]*(Gi2[x1][s_1]-b2) + 
						                 l1*(a2+c2-b2-b2));
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
Lat2DCyl1stO::Propagate2G(Matrix Gi1, 
						  Matrix Gi2, 
						  const Vector G, 
						  const int s, 
						  const LatticeRange* LatRange,
						  const double f) const {
	int x,y,x_1,x1,z;
	double a1,b1,c1,a2,b2,c2;

	double f_1= exp(-f);
	double f1= exp(f);
	double norm = l1*(f_1+f1-2)+1;
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
				Gi1[z][s] = G[z] * (b1 + lambda1_1[x]*(Gi1[x_1][s_1] + Gi2[x_1][s_1] - b1) +
					                      lambda11[x]*(Gi1[x1][s_1] + Gi2[x1][s_1] - b1) + 
					                      l1*((a1 + a2)*f1 + (c1 + c2)*f_1 - b1  - b1))/norm;
				Gi2[z][s] = 0;
			} else {
				Gi1[z][s] = G[z] * (b1 + lambda1_1[x]*(Gi1[x_1][s_1]-b1) + 
						                 lambda11[x]*(Gi1[x1][s_1]-b1) + 
						                 l1*(a1*f1 + c1*f_1 - b1 - b1))/norm;
				Gi2[z][s] = G[z] * (b2 + lambda1_1[x]*(Gi2[x_1][s_1]-b2) + 
						                 lambda11[x]*(Gi2[x1][s_1]-b2) + 
						                 l1*(a2*f1 + c2*f_1 - b2 - b2))/norm;
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
Lat2DCyl1stO::Propagate2G(Vector Gi1, 
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
				Gi1[z] = G[z] * (b1 + lambda1_1[x]*(Gi1Prime[y] + Gi2Prime[y] - b1) +
					                  lambda11[x]*(Gi1[x1] + Gi2[x1] - b1) + 
					                  l1*(a1 + c1 + a2 + c2 - b1 - b1));
				Gi2[z] = 0;
			} else {
				Gi1[z] = G[z] * (b1 + lambda1_1[x]*(Gi1Prime[y]-b1) + 
						              lambda11[x]*(Gi1[x1]-b1) + 
						              l1*(a1+c1-b1-b1));
				Gi2[z] = G[z] * (b2 + lambda1_1[x]*(Gi2Prime[y]-b2) + 
						              lambda11[x]*(Gi2[x1]-b2) + 
						              l1*(a2+c2-b2-b2));
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
Lat2DCyl1stO::Propagate2G(Vector Gi1, 
						  Vector Gi2, 
						  const Vector G, 
						  const LatticeRange* LatRange,
						  const double f) const {
	int x,y,x1,z;
	double a1,b1,c1,a2,b2,c2;

	double f_1= exp(-f);
	double f1= exp(f);
	double norm = l1*(f_1+f1-2)+1;
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
				Gi1[z] = G[z] * (b1 + lambda1_1[x]*(Gi1Prime[y] + Gi2Prime[y] - b1) +
					                  lambda11[x]*(Gi1[x1] + Gi2[x1] - b1) + 
					                  l1*((a1 + a2)*f1 + (c1 + c2)*f_1 - b1 - b1))/norm;
				Gi2[z] = 0;
			} else {
				Gi1[z] = G[z] * (b1 + lambda1_1[x]*(Gi1Prime[y]-b1) + 
						              lambda11[x]*(Gi1[x1]-b1) + 
						              l1*(a1*f1 + c1*f_1- b1  - b1))/norm;
				Gi2[z] = G[z] * (b2 + lambda1_1[x]*(Gi2Prime[y]-b2) + 
						              lambda11[x]*(Gi2[x1]-b2) + 
						              l1*(a2*f1 + c2*f_1 - b2 - b2))/norm;
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



double
Lat2DCyl1stO::ComputeLnGN(const Vector Gi) const {
	int x,y,z;
	double value = 0;
	z = MY+2;
	for (x=2; x<MX; x++) {
		for (y=2; y<MY; y++) {
			value += L[x]*Gi[z++];
		}
		z+=2;
	}
	if (boundX1 == 3) {
		z = MY+2;
		for (y=2; y<MY; y++) {
			value -= 0.5*L[2]*Gi[z++];
		}
	}
	if (boundX2 == MX-2) {
		z = MY*(MX-2)+2;
		for (y=2; y<MY; y++) {
			value -= 0.5*L[MX-1]*Gi[z++];
		}
	}
	if (boundY1 == 3) {
		z = MY+2;
		for (x=2; x<MX; x++) {
			value -= 0.5*L[x]*Gi[z];
			z+=MY;
		}
	}
	if (boundY2 == MY-2) {
		z = 2*MY-2;
		for (x=2; x<MX; x++) {
			value -= 0.5*L[x]*Gi[z];
			z+=MY;
		}
	}
	if (boundX1 == 3 && boundY1 == 3) {
		value += 0.25*L[2]*Gi[MY+2];
	}
	if (boundX1 == 3 && boundY2 == MY-2) {
		value += 0.25*L[2]*Gi[2*MY-2];
	}
	if (boundX2 == MX-2 && boundY1 == 3) {
		value += 0.25*L[MX-1]*Gi[MY*(MX-2)+2];
	}
	if (boundX2 == MX-2 && boundY2 == MY-2) {
		value += 0.25*L[MX-1]*Gi[(MX-1)*MY-2];
	}
	return log(value);
}

