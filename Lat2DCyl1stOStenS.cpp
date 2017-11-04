#include <iostream>
#include "Lat2DCyl1stOStenS.h"

Lat2DCyl1stOStenS::Lat2DCyl1stOStenS(Input* MyInput_, Text name_)
	: Lat2DFlat(MyInput_,name_),
	Lat2DCylinder(MyInput_,name_),
	Lat2DCylinderSten(MyInput_,name_) {
	int x;
	double lmbda0;

	// a log() is about 8x as expensive as a multiplication, so we cache it
	// if lambda <= 0 the log is not taken an a value of 0 is written
	// to the array. Of course, this value may not be used in the
	// calculations, so always check if lambda >= 0 before using loglambda.
	logl0.Dim(2,MX-1);
	logl11.Dim(2,MX-1);
	logl1_1.Dim(2,MX-1);
	logl21.Dim(2,MX-1);
	logl2_1.Dim(2,MX-1);
	for (x=2; x<MX; x++) {
		lmbda0 = 1.0 - lambda1_1[x]-lambda11[x]-2*l1;
		logl0[x] = (lmbda0>0) ?log(lmbda0) :0;
		logl11[x] = (lambda11[x]>0) ?log(lambda11[x]) :0;
		logl1_1[x] = (lambda1_1[x]>0) ?log(lambda1_1[x]) :0;
		logl21[x] = (lambda21[x]>0) ?log(lambda21[x]) :0;
		logl2_1[x] = (lambda2_1[x]>0) ?log(lambda2_1[x]) :0;
	}
	logl1 = (l1>0) ?log(l1) :0;

	Gi1Prime.Dim(1,MY);
	Gi2Prime.Dim(1,MY);
}
Lat2DCyl1stOStenS::~Lat2DCyl1stOStenS() {
}
/* Below, if we use an array for values associated with an array, the following 
 * scheme is used:
 *
 *   |  4|3|5
 *  x|  0|8|1
 *  \/  6|2|7
 *       --->
 *        y
 * Here, x, is the direction perpendicular to the cylinder axis and y the
 * direction along the axis.
 * l0 and l2 values can often be equal or smaller than 0, so placing them at
 * the end of the array speeds up, finding a lambda > 0.
 */

void
Lat2DCyl1stOStenS::GetLatticeInfo(int* Info) const{ 
}


void 
Lat2DCyl1stOStenS::PropagateG(Matrix Gi, const Vector G, const int s) const {
	short i, j;
	int x, y, z, y_1, y1, s_1;
	double lnterm;
	double neigh[9]; // neighbour values for G
	double nl[9];    // neighbour lambdas
	double lognl[9]; // neighbour log(lambda)

 

	s_1 = s-1;
	y_1 = 3;
	z = MY + 2;
	y1 = MY + MY + 3;
	nl[0] = nl[1] = l1;
	lognl[0] = lognl[1] = logl1;
	for (x=2; x<MX; x++) {
		neigh[3] = Gi[y_1-2][s_1];
		neigh[5] = Gi[y_1-1][s_1];
		neigh[8] = Gi[z-1][s_1];
		neigh[1] = Gi[z][s_1];
		neigh[2] = Gi[y1-2][s_1];
		neigh[7] = Gi[y1-1][s_1];
		nl[4] = nl[5] = lambda2_1[x];
		nl[3] = lambda1_1[x];
		nl[6] = nl[7] = lambda21[x];
		nl[2] = lambda11[x];
		nl[8] = 1.0-l1-l1-nl[2]-nl[3]-nl[4]-nl[5]-nl[6]-nl[7];
		lognl[4] = lognl[5] = logl2_1[x];
		lognl[3] = logl1_1[x];
		lognl[6] = lognl[7] = logl21[x];
		lognl[2] = logl11[x];
		lognl[8] = logl0[x];
		for (y=2; y<MY; y++) {
			neigh[4] = neigh[3];
			neigh[3] = neigh[5];
			neigh[5] = Gi[y_1++][s_1];
			neigh[0] = neigh[8];
			neigh[8] = neigh[1];
			neigh[1] = Gi[z+1][s_1];
			neigh[6] = neigh[2];
			neigh[2] = neigh[7];
			neigh[7] = Gi[y1++][s_1];
			if (G[z] == LOGZERO) {
				Gi[z++][s] = LOGZERO;
				continue;
			}
			i = 0;
			while (i<9 && (neigh[i] == LOGZERO || nl[i] <= 0)) i++;
			if (i == 9) {
				Gi[z++][s] = LOGZERO;
				continue;
			}
			lnterm = 0;
			for (j=0; j<9; j++) {
				if (j != i && neigh[j] != LOGZERO) {
					lnterm += nl[j]*exp(neigh[j]-neigh[i]);
				}
			}
			Gi[z][s] = G[z] + lognl[i] + neigh[i]
					+ log1p(lnterm/nl[i]);
			z++;
		}
		y_1+=2;
		y1+=2;
		z+=2;
	}
	UpdateBoundaries(Gi,s);
}

void 
Lat2DCyl1stOStenS::PropagateG(Matrix Gi, const Vector G, const int s,const double f) const {
	Message(fatal,MyInput,"Force ensemble not implemented in Lat2DCyl1stOStenS");
}

void
Lat2DCyl1stOStenS::PropagateG(const Vector Gi, const Vector G, Vector Gout) const {
	short i, j;
	int x, y, z, y1;
	double lnterm;
	double neigh[9]; // neighbour values for G
	double nl[9];    // neighbour lambdas
	double lognl[9]; // neighbour log(lambda)



	for (y=1; y<=MY; y++) {
		Gi1Prime[y] = Gi[y];
	}
	z = MY + 2;
	y1 = MY + MY + 3;
	nl[0] = nl[1] = l1;
	lognl[0] = lognl[1] = logl1;
	for (x=2; x<MX; x++) {
		neigh[3] = Gi1Prime[1];
		neigh[5] = Gi1Prime[2];
		neigh[8] = Gi1Prime[1] = Gi[z-1];
		neigh[1] = Gi1Prime[2] = Gi[z];
		neigh[2] = Gi[y1-2];
		neigh[7] = Gi[y1-1];
		nl[4] = nl[5] = lambda2_1[x];
		nl[3] = lambda1_1[x];
		nl[6] = nl[7] = lambda21[x];
		nl[2] = lambda11[x];
		nl[8] = 1.0-l1-l1-nl[2]-nl[3]-nl[4]-nl[5]-nl[6]-nl[7];
		lognl[4] = lognl[5] = logl2_1[x];
		lognl[3] = logl1_1[x];
		lognl[6] = lognl[7] = logl21[x];
		lognl[2] = logl11[x];
		lognl[8] = logl0[x];
		for (y=3; y<=MY; y++) {
			neigh[4] = neigh[3];
			neigh[3] = neigh[5];
			neigh[5] = Gi1Prime[y];
			neigh[0] = neigh[8];
			neigh[8] = neigh[1];
			neigh[1] = Gi1Prime[y] = Gi[z+1];
			neigh[6] = neigh[2];
			neigh[2] = neigh[7];
			neigh[7] = Gi[y1++];
			if (G[z] == LOGZERO) {
				Gout[z++] = LOGZERO;
				continue;
			}
			i = 0;
			while (i<9 && (neigh[i] == LOGZERO || nl[i] <= 0)) i++;
			if (i == 9) {
				Gout[z++] = LOGZERO;
				continue;
			}
			lnterm = 0;
			for (j=0; j<9; j++) {
				if (j != i && neigh[j] != LOGZERO) {
					lnterm += nl[j]*exp(neigh[j]-neigh[i]);
				}
			}
			Gout[z] = G[z] + lognl[i] + neigh[i]
					+ log1p(lnterm/nl[i]);
			z++;
		}
		y1+=2;
		z+=2;
	}
	UpdateBoundaries(Gout);
}

void
Lat2DCyl1stOStenS::PropagateG(const Vector Gi, const Vector G, Vector Gout,const double f) const {
	Message(fatal,MyInput,"Force ensemble not implemented in Lat2DCyl1stOStenS");
}


void
Lat2DCyl1stOStenS::PropagateG(Vector Gi, const Vector G) const {
	short i, j;
	int x, y, z, y1;
	double lnterm;
	double neigh[9]; // neighbour values for G
	double nl[9];    // neighbour lambdas
	double lognl[9]; // neighbour log(lambda)


	for (y=1; y<=MY; y++) {
		Gi1Prime[y] = Gi[y];
	}
	z = MY + 2;
	y1 = MY + MY + 3;
	nl[0] = nl[1] = l1;
	lognl[0] = lognl[1] = logl1;
	for (x=2; x<MX; x++) {
		neigh[3] = Gi1Prime[1];
		neigh[5] = Gi1Prime[2];
		neigh[8] = Gi1Prime[1] = Gi[z-1];
		neigh[1] = Gi1Prime[2] = Gi[z];
		neigh[2] = Gi[y1-2];
		neigh[7] = Gi[y1-1];
		nl[4] = nl[5] = lambda2_1[x];
		nl[3] = lambda1_1[x];
		nl[6] = nl[7] = lambda21[x];
		nl[2] = lambda11[x];
		nl[8] = 1.0-l1-l1-nl[2]-nl[3]-nl[4]-nl[5]-nl[6]-nl[7];
		lognl[4] = lognl[5] = logl2_1[x];
		lognl[3] = logl1_1[x];
		lognl[6] = lognl[7] = logl21[x];
		lognl[2] = logl11[x];
		lognl[8] = logl0[x];
		for (y=3; y<=MY; y++) {
			neigh[4] = neigh[3];
			neigh[3] = neigh[5];
			neigh[5] = Gi1Prime[y];
			neigh[0] = neigh[8];
			neigh[8] = neigh[1];
			neigh[1] = Gi1Prime[y] = Gi[z+1];
			neigh[6] = neigh[2];
			neigh[2] = neigh[7];
			neigh[7] = Gi[y1++];
			if (G[z] == LOGZERO) {
				Gi[z++] = LOGZERO;
				continue;
			}
			i = 0;
			while (i<9 && (neigh[i] == LOGZERO || nl[i] <= 0)) i++;
			if (i == 9) {
				Gi[z++] = LOGZERO;
				continue;
			}
			lnterm = 0;
			for (j=0; j<9; j++) {
				if (j != i && neigh[j] != LOGZERO) {
					lnterm += nl[j]*exp(neigh[j]-neigh[i]);
				}
			}
			Gi[z] = G[z] + lognl[i] + neigh[i]
					+ log1p(lnterm/nl[i]);
			z++;
		}
		y1+=2;
		z+=2;
	}
	UpdateBoundaries(Gi);
}

void
Lat2DCyl1stOStenS::PropagateG(Vector Gi, const Vector G, const double f) const {
	Message(fatal,MyInput,"Force ensemble not implemented in Lat2DCyl1stOStenS");
}

void
Lat2DCyl1stOStenS::Propagate2G(Matrix Gi1, Matrix Gi2, const Vector G, const int s,
			   const LatticeRange* LatRange) const {
	bool inrange;
	short i1, i2, j;
	int x, y, z, y_1, y1, s_1;
	double lnterm;
	double n1[9];    // neighbour values for Gi1
	double n2[9];    // neighbour values for Gi2
	double nl[9];    // neighbour lambdas
	double lognl[9]; // neighbour log(lambda)


	s_1 = s-1;
	y_1 = 3;
	z = MY + 2;
	y1 = MY + MY + 3;
	nl[0] = nl[1] = l1;
	lognl[0] = lognl[1] = logl1;
	for (x=2; x<MX; x++) {
		n1[3] = Gi1[y_1-2][s_1]; n1[5] = Gi1[y_1-1][s_1];
		n1[8] = Gi1[z-1][s_1];   n1[1] = Gi1[z][s_1];
		n1[2] = Gi1[y1-2][s_1];  n1[7] = Gi1[y1-1][s_1];
		n2[3] = Gi2[y_1-2][s_1]; n2[5] = Gi2[y_1-1][s_1];
		n2[8] = Gi2[z-1][s_1];   n2[1] = Gi2[z][s_1];
		n2[2] = Gi2[y1-2][s_1];  n2[7] = Gi2[y1-1][s_1];
		nl[4] = nl[5] = lambda2_1[x];
		nl[3] = lambda1_1[x];
		nl[6] = nl[7] = lambda21[x];
		nl[2] = lambda11[x];
		nl[8] = 1.0-l1-l1-nl[2]-nl[3]-nl[4]-nl[5]-nl[6]-nl[7];
		lognl[4] = lognl[5] = logl2_1[x];
		lognl[3] = logl1_1[x];
		lognl[6] = lognl[7] = logl21[x];
		lognl[2] = logl11[x];
		lognl[8] = logl0[x];
		for (y=2; y<MY; y++) {
			n1[4] = n1[3]; n1[3] = n1[5]; n1[5] = Gi1[y_1][s_1];
			n1[0] = n1[8]; n1[8] = n1[1]; n1[1] = Gi1[z+1][s_1];
			n1[6] = n1[2]; n1[2] = n1[7]; n1[7] = Gi1[y1][s_1];
			n2[4] = n2[3]; n2[3] = n2[5]; n2[5] = Gi2[y_1++][s_1];
			n2[0] = n2[8]; n2[8] = n2[1]; n2[1] = Gi2[z+1][s_1];
			n2[6] = n2[2]; n2[2] = n2[7]; n2[7] = Gi2[y1++][s_1];
			if (G[z] == LOGZERO) {
				Gi1[z][s] = LOGZERO;
				Gi2[z++][s] = LOGZERO;
				continue;
			}
			inrange = LatRange->InRange(z);
			i1 = 0;
			while (i1<9 && (n1[i1] == LOGZERO || nl[i1] <= 0)) i1++;
			i2 = 9;
			if (!inrange || i1 == 9) { // we need i2 too
				i2 = 0;
				while (i2<9
					&& (n2[i2] == LOGZERO || nl[i2]<=0)) {
					i2++;
				}
			}
			lnterm = 0;
			if (i1 != 9) {
				for (j=0; j<9; j++) {
					if (j != i1 && n1[j] != LOGZERO) {
						lnterm += nl[j]
							*exp(n1[j]-n1[i1]);
					}
					if (inrange && n2[j] != LOGZERO) {
						lnterm += nl[j]
							*exp(n2[j]-n1[i1]);
					}
				}
				Gi1[z][s] = G[z] + lognl[i1] + n1[i1]
					+ log1p(lnterm/nl[i1]);
			} else if (inrange && i2 != 9) {
				for (j=0; j<9; j++) {
					if (n1[j] != LOGZERO) {
						lnterm += nl[j]
							*exp(n1[j]-n2[i2]);
					}
					if (j != i2 && n2[j] != LOGZERO) {
						lnterm += nl[j]
							*exp(n2[j]-n2[i2]);
					}
				}
				Gi1[z][s] = G[z] + lognl[i2] + n1[i2]
					+ log1p(lnterm/nl[i2]);
			} else {
				Gi1[z][s] = LOGZERO;
			}
			if (inrange || i2 == 9) {
				Gi2[z][s] = LOGZERO;
			} else {
				lnterm = 0;
				for (j=0; j<9; j++) {
					if (j != i2 && n2[j] != LOGZERO) {
						lnterm += nl[j]
							* exp(n2[j]-n2[i2]);
					}
				}
				Gi2[z][s] = G[z] + lognl[i2] + n2[i2]
					+ log1p(lnterm/nl[i2]);
			}
			z++;
		}
		y_1+=2;
		y1+=2;
		z+=2;
	}
	UpdateBoundaries(Gi1,s);
	UpdateBoundaries(Gi2,s);
}


void
Lat2DCyl1stOStenS::Propagate2G(Matrix Gi1, Matrix Gi2, const Vector G, const int s,
			   const LatticeRange* LatRange,const double f) const {
	Message(fatal,MyInput,"Force ensemble not implemented in Lat2DCyl1stOStenS");
}

void
Lat2DCyl1stOStenS::Propagate2G(Vector Gi1, Vector Gi2, const Vector G,
		const LatticeRange* LatRange) const {
	bool inrange;
	short i1, i2, j;
	int x, y, z, y1;
	double lnterm;
	double n1[9];    // neighbour values for Gi1
	double n2[9];    // neighbour values for Gi2
	double nl[9];    // neighbour lambdas
	double lognl[9]; // neighbour log(lambda)


	for (y=1; y<=MY; y++) {
		Gi1Prime[y] = Gi1[y];
		Gi2Prime[y] = Gi2[y];
	}
	z = MY + 2;
	y1 = MY + MY + 3;
	nl[0] = nl[1] = l1;
	lognl[0] = lognl[1] = logl1;
	for (x=2; x<MX; x++) {
		n1[3] = Gi1Prime[1]; n1[5] = Gi1Prime[2];
		n1[8] = Gi1[z-1];    n1[1] = Gi1[z];
		n1[2] = Gi1[y1-2];   n1[7] = Gi1[y1-1];
		Gi1Prime[1] = Gi1[z-1]; Gi1Prime[2] = Gi1[z];
		n2[3] = Gi2Prime[1]; n2[5] = Gi2Prime[2];
		n2[8] = Gi2[z-1];    n2[1] = Gi2[z];
		n2[2] = Gi2[y1-2];   n2[7] = Gi2[y1-1];
		Gi2Prime[1] = Gi2[z-1]; Gi2Prime[2] = Gi2[z];
		nl[4] = nl[5] = lambda2_1[x];
		nl[3] = lambda1_1[x];
		nl[6] = nl[7] = lambda21[x];
		nl[2] = lambda11[x];
		nl[8] = 1.0-l1-l1-nl[2]-nl[3]-nl[4]-nl[5]-nl[6]-nl[7];
		lognl[4] = lognl[5] = logl2_1[x];
		lognl[3] = logl1_1[x];
		lognl[6] = lognl[7] = logl21[x];
		lognl[2] = logl11[x];
		lognl[8] = logl0[x];
		for (y=3; y<=MY; y++) {
			n1[4] = n1[3]; n1[3] = n1[5]; n1[5] = Gi1Prime[y];
			n1[0] = n1[8]; n1[8] = n1[1]; n1[1] = Gi1[z+1];
			n1[6] = n1[2]; n1[2] = n1[7]; n1[7] = Gi1[y1];
			Gi1Prime[y] = n1[1];
			n2[4] = n2[3]; n2[3] = n2[5]; n2[5] = Gi2Prime[y];
			n2[0] = n2[8]; n2[8] = n2[1]; n2[1] = Gi2[z+1];
			n2[6] = n2[2]; n2[2] = n2[7]; n2[7] = Gi2[y1++];
			Gi2Prime[y] = n2[1];
			if (G[z] == LOGZERO) {
				Gi1[z] = LOGZERO;
				Gi2[z++] = LOGZERO;
				continue;
			}
			inrange = LatRange->InRange(z);
			i1 = 0;
			while (i1<9 && (n1[i1] == LOGZERO || nl[i1] <= 0)) i1++;
			i2 = 9;
			if (!inrange || i1 == 9) { // we need i2 too
				i2 = 0;
				while (i2<9
					&& (n2[i2] == LOGZERO || nl[i2]<=0)) {
					i2++;
				}
			}
			lnterm = 0;
			if (i1 != 9) {
				for (j=0; j<9; j++) {
					if (j != i1 && n1[j] != LOGZERO) {
						lnterm += nl[j]
							*exp(n1[j]-n1[i1]);
					}
					if (inrange && n2[j] != LOGZERO) {
						lnterm += nl[j]
							*exp(n2[j]-n1[i1]);
					}
				}
				Gi1[z] = G[z] + lognl[i1] + n1[i1]
					+ log1p(lnterm/nl[i1]);
			} else if (inrange && i2 != 9) {
				for (j=0; j<9; j++) {
					if (n1[j] != LOGZERO) {
						lnterm += nl[j]
							*exp(n1[j]-n2[i2]);
					}
					if (j != i2 && n2[j] != LOGZERO) {
						lnterm += nl[j]
							*exp(n2[j]-n2[i2]);
					}
				}
				Gi1[z] = G[z] + lognl[i2] + n1[i2]
					+ log1p(lnterm/nl[i2]);
			} else {
				Gi1[z] = LOGZERO;
			}
			if (inrange || i2 == 9) {
				Gi2[z] = LOGZERO;
			} else {
				lnterm = 0;
				for (j=0; j<9; j++) {
					if (j != i2 && n2[j] != LOGZERO) {
						lnterm += nl[j]
							* exp(n2[j]-n2[i2]);
					}
				}
				Gi2[z] = G[z] + lognl[i2] + n2[i2]
					+ log1p(lnterm/nl[i2]);
			}
			z++;
		}
		y1+=2;
		z+=2;
	}
	UpdateBoundaries(Gi1);
	UpdateBoundaries(Gi2);
}
void
Lat2DCyl1stOStenS::Propagate2G(Vector Gi1, Vector Gi2, const Vector G,
		const LatticeRange* LatRange,const double f) const {
	Message(fatal,MyInput,"Force ensemble not implemented in Lat2DCyl1stOStenS");
}

