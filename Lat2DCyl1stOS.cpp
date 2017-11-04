#include "Lat2DCyl1stOS.h"

Lat2DCyl1stOS::Lat2DCyl1stOS(Input* MyInput_, Text name_)
	: Lat2DFlat(MyInput_,name_), Lat2DCylinder(MyInput_,name_) {
	int x;
	double lmbda0;

	// a log() is about 8x as expensive as a multiplication, so we cache it
	// if lambda <= 0 the log is not taken an a value of 0 is written
	// to the array. Of course, this value may not be used in the
	// calculations, so always check if lambda >= 0 before using loglambda.
	logl11.Dim(2,MX-1);
	logl1_1.Dim(2,MX-1);
	logl0.Dim(2,MX-1);
	for (x=2; x<MX; x++) {
		logl11[x] = (lambda11[x]>0) ?log(lambda11[x]) :0;
		logl1_1[x] = (lambda1_1[x]>0) ?log(lambda1_1[x]) :0;
		lmbda0 = 1.0 - lambda1_1[x]-lambda11[x]-2*l1;
		logl0[x] = (lmbda0>0) ?log(lmbda0) :0;
	}
	logl1 = (l1>0) ?log(l1) :0;

	Gi1Prime.Dim(1,MY);
	Gi2Prime.Dim(1,MY);
}
Lat2DCyl1stOS::~Lat2DCyl1stOS() {
}
/* Below, if we use an array for values associated with an array, the following 
 * scheme is used:
 *
 *   |  _|3|_
 *  x|  0|1|2
 *  \/   |4| 
 *       --->
 *        y
 * Here, x, is the direction perpendicular to the cylinder axis and y the
 * direction along the axis.
 */

void
Lat2DCyl1stOS::GetLatticeInfo(int* Info) const{ 
}


void 
Lat2DCyl1stOS::PropagateG(Matrix Gi, const Vector G, const int s) const {
	short i, j;
	int x, y, z, y_1, y1, s_1;
	double lnterm;
	double neigh[5]; // neighbour values for G
	double nl[5];    // neighbour lambdas
	double lognl[5]; // neighbour log(lambda)

	s_1 = s-1;
	y_1 = 2;
	z = MY + 2;
	y1 = MY + MY + 2;
	nl[0] = nl[2] = l1;
	lognl[0] = lognl[2] = logl1;
	for (x=2; x<MX; x++) {
		neigh[1] = Gi[z-1][s_1];
		neigh[2] = Gi[z][s_1];
		nl[3] = lambda1_1[x];
		nl[4] = lambda11[x];
		nl[1] = 1.0 - nl[0] - nl[2] - l1 - l1;
		lognl[3] = logl1_1[x];
		lognl[1] = logl0[x];
		lognl[4] = logl11[x];
		for (y=2; y<MY; y++) {
			neigh[0] = neigh[1];
			neigh[1] = neigh[2];
			neigh[2] = Gi[z+1][s_1];
			neigh[3] = Gi[y_1++][s_1];
			neigh[4] = Gi[y1++][s_1];
			if (G[z] == LOGZERO) {
				Gi[z++][s] = LOGZERO;
				continue;
			}
			i = 0;
			while (i<5 && (neigh[i] == LOGZERO || nl[i] <= 0)) i++;
			if (i == 5) {
				Gi[z++][s] = LOGZERO;
				continue;
			}
			lnterm = 0;
			for (j=0; j<5; j++) {
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
Lat2DCyl1stOS::PropagateG(Matrix Gi, const Vector G, const int s,const double f) const {
	Message(fatal,MyInput,"Force ensemble not implemented in Lat2DCyl1stOS");
}


void
Lat2DCyl1stOS::PropagateG(Vector Gi, const Vector G) const {
	short i, j;
	int x, y, z, y1;
	double lnterm;
	double neigh[5]; // neighbour values for G
	double nl[5];    // neighbour lambdas
	double lognl[5]; // neighbour log(lambda)

	for (y=2; y<MY; y++) {
		Gi1Prime[y] = Gi[y];
	}
	z = MY + 2;
	y1 = MY + MY + 2;
	nl[0] = nl[2] = l1;
	lognl[0] = lognl[2] = logl1;
	for (x=2; x<MX; x++) {
		neigh[1] = Gi[z-1];
		neigh[2] = Gi[z];
		nl[3] = lambda1_1[x];
		nl[4] = lambda11[x];
		nl[1] = 1.0 - nl[0] - nl[2] - l1 - l1;
		lognl[3] = logl1_1[x];
		lognl[1] = logl0[x];
		lognl[4] = logl11[x];
		for (y=2; y<MY; y++) {
			neigh[0] = neigh[1];
			neigh[1] = neigh[2];
			neigh[2] = Gi[z+1];
			neigh[3] = Gi1Prime[y];
			neigh[4] = Gi[y1++];
			Gi1Prime[y] = neigh[1];
			if (G[z] == LOGZERO) {
				Gi[z++] = LOGZERO;
				continue;
			}
			i = 0;
			while (i<5 && (neigh[i] == LOGZERO || nl[i] <= 0)) i++;
			if (i == 5) {
				Gi[z++] = LOGZERO;
				continue;
			}
			lnterm = 0;
			for (j=0; j<5; j++) {
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
Lat2DCyl1stOS::PropagateG(Vector Gi, const Vector G, const double f) const {
	Message(fatal,MyInput,"Force ensemble not implemented in Lat2DCyl1stOS");
}


void
Lat2DCyl1stOS::PropagateG(const Vector Gi, const Vector G, Vector Gout) const {
	short i, j;
	int x, y, z, y1;
	double lnterm;
	double neigh[5]; // neighbour values for G
	double nl[5];    // neighbour lambdas
	double lognl[5]; // neighbour log(lambda)

	for (y=2; y<MY; y++) {
		Gi1Prime[y] = Gi[y];
	}
	z = MY + 2;
	y1 = MY + MY + 2;
	nl[0] = nl[2] = l1;
	lognl[0] = lognl[2] = logl1;
	for (x=2; x<MX; x++) {
		neigh[1] = Gi[z-1];
		neigh[2] = Gi[z];
		nl[3] = lambda1_1[x];
		nl[4] = lambda11[x];
		nl[1] = 1.0 - nl[0] - nl[2] - l1 - l1;
		lognl[3] = logl1_1[x];
		lognl[1] = logl0[x];
		lognl[4] = logl11[x];
		for (y=2; y<MY; y++) {
			neigh[0] = neigh[1];
			neigh[1] = neigh[2];
			neigh[2] = Gi[z+1];
			neigh[3] = Gi1Prime[y];
			neigh[4] = Gi[y1++];
			Gi1Prime[y] = neigh[1];
			if (G[z] == LOGZERO) {
				Gout[z++] = LOGZERO;
				continue;
			}
			i = 0;
			while (i<5 && (neigh[i] == LOGZERO || nl[i] <= 0)) i++;
			if (i == 5) {
				Gout[z++] = LOGZERO;
				continue;
			}
			lnterm = 0;
			for (j=0; j<5; j++) {
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
Lat2DCyl1stOS::PropagateG(const Vector Gi, const Vector G, Vector Gout,const double f) const {
	Message(fatal,MyInput,"Force ensemble not implemented in Lat2DCyl1stOS");
}


void
Lat2DCyl1stOS::Propagate2G(Matrix Gi1, Matrix Gi2, const Vector G, const int s,
			   const LatticeRange* LatRange) const {
	bool inrange;
	short i1, i2, j;
	int x, y, z, y_1, y1, s_1;
	double lnterm;
	double n1[5];    // neighbour values for Gi1
	double n2[5];    // neighbour values for Gi2
	double nl[5];    // neighbour lambdas
	double lognl[5]; // neighbour log(lambda)

	s_1 = s-1;
	y_1 = 2;
	z = MY + 2;
	y1 = MY + MY + 2;
	nl[0] = nl[2] = l1;
	lognl[0] = lognl[2] = logl1;
	for (x=2; x<MX; x++) {
		n1[1] = Gi1[z-1][s_1];
		n1[2] = Gi1[z][s_1];
		n2[1] = Gi2[z-1][s_1];
		n2[2] = Gi2[z][s_1];
		nl[3] = lambda1_1[x];
		nl[4] = lambda11[x];
		nl[1] = 1.0 - nl[0] - nl[2] - l1 - l1;
		lognl[3] = logl1_1[x];
		lognl[1] = logl0[x];
		lognl[4] = logl11[x];
		for (y=2; y<MY; y++) {
			n1[0] = n1[1];
			n1[1] = n1[2];
			n1[2] = Gi1[z+1][s_1];
			n1[3] = Gi1[y_1][s_1];
			n1[4] = Gi1[y1][s_1];
			n2[0] = n2[1];
			n2[1] = n2[2];
			n2[2] = Gi2[z+1][s_1];
			n2[3] = Gi2[y_1++][s_1];
			n2[4] = Gi2[y1++][s_1];
			if (G[z] == LOGZERO) {
				Gi1[z][s] = LOGZERO;
				Gi2[z++][s] = LOGZERO;
				continue;
			}
			inrange = LatRange->InRange(z);
			i1 = 0;
			while (i1<5 && (n1[i1] == LOGZERO || nl[i1] <= 0)) i1++;
			i2 = 5;
			if (!inrange || i1 == 5) { // we need i2 too
				i2 = 0;
				while (i2<5
					&& (n2[i2] == LOGZERO || nl[i2]<=0)) {
					i2++;
				}
			}
			lnterm = 0;
			if (i1 != 5) {
				for (j=0; j<5; j++) {
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
			} else if (inrange && i2 != 5) {
				for (j=0; j<5; j++) {
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
			if (inrange || i2 == 5) {
				Gi2[z][s] = LOGZERO;
			} else {
				lnterm = 0;
				for (j=0; j<5; j++) {
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
Lat2DCyl1stOS::Propagate2G(Matrix Gi1, Matrix Gi2, const Vector G, const int s,
			   const LatticeRange* LatRange,const double f) const {
	Message(fatal,MyInput,"Force ensemble not implemented in Lat2DCyl1stOS");
}

void
Lat2DCyl1stOS::Propagate2G(Vector Gi1, Vector Gi2, const Vector G,
		const LatticeRange* LatRange) const {
	bool inrange;
	short i1, i2, j;
	int x, y, z, y1;
	double lnterm;
	double n1[5];    // neighbour values for Gi1
	double n2[5];    // neighbour values for Gi2
	double nl[5];    // neighbour lambdas
	double lognl[5]; // neighbour log(lambda)

	for (y=2; y<MY; y++) {
		Gi1Prime[y] = Gi1[y];
		Gi2Prime[y] = Gi2[y];
	}
	z = MY + 2;
	y1 = MY + MY + 2;
	nl[0] = nl[2] = l1;
	lognl[0] = lognl[2] = logl1;
	for (x=2; x<MX; x++) {
		n1[1] = Gi1[z-1];
		n1[2] = Gi1[z];
		n2[1] = Gi2[z-1];
		n2[2] = Gi2[z];
		nl[3] = lambda1_1[x];
		nl[4] = lambda11[x];
		nl[1] = 1.0 - nl[0] - nl[2] - l1 - l1;
		lognl[3] = logl1_1[x];
		lognl[1] = logl0[x];
		lognl[4] = logl11[x];
		for (y=2; y<MY; y++) {
			n1[0] = n1[1];
			n1[1] = n1[2];
			n1[2] = Gi1[z+1];
			n1[3] = Gi1Prime[y];
			n1[4] = Gi1[y1];
			Gi1Prime[y] = n1[1];
			n2[0] = n2[1];
			n2[1] = n2[2];
			n2[2] = Gi2[z+1];
			n2[3] = Gi2Prime[y];
			n2[4] = Gi2[y1++];
			Gi2Prime[y] = n2[1];
			if (G[z] == LOGZERO) {
				Gi1[z] = LOGZERO;
				Gi2[z++] = LOGZERO;
				continue;
			}
			inrange = LatRange->InRange(z);
			i1 = 0;
			while (i1<5 && (n1[i1] == LOGZERO || nl[i1] <= 0)) i1++;
			i2 = 5;
			if (!inrange || i1 == 5) { // we need i2 too
				i2 = 0;
				while (i2<5
					&& (n2[i2] == LOGZERO || nl[i2]<=0)) {
					i2++;
				}
			}
			lnterm = 0;
			if (i1 != 5) {
				for (j=0; j<5; j++) {
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
			} else if (inrange && i2 != 5) {
				for (j=0; j<5; j++) {
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
			if (inrange || i2 == 5) {
				Gi2[z] = LOGZERO;
			} else {
				lnterm = 0;
				for (j=0; j<5; j++) {
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
Lat2DCyl1stOS::Propagate2G(Vector Gi1, Vector Gi2, const Vector G,
		const LatticeRange* LatRange,const double f) const {
	Message(fatal,MyInput,"Force ensemble not implemented in Lat2DCyl1stOS");
}

double
Lat2DCyl1stOS::ComputeLnGN(const Vector Gi) const {
	int x,y,z;
	double value;
	double LogLx;
	double log2, log4;

	value = LOGZERO;
	log2 = log(2.0);
	log4 = log(4.0);
	z = MY+2;
	for (x=2; x<MX; x++) {
		LogLx = log(L[x]);
		for (y=2; y<MY; y++) {
			addLntoLn(Gi[z++], LogLx, value, true);
		}
		z+=2;
	}
	if (boundX1 == 3) {
		z = MY+2;
		LogLx = log(L[2]);
		for (y=2; y<MY; y++) {
			addLntoLn(Gi[z++], LogLx-log2, value, false);
		}
	}		
	if (boundX2 == MX-2) {
		z = MY*(MX-2)+2;
		LogLx = log(L[MX-1]);
		for (y=2; y<MY; y++) {
			addLntoLn(Gi[z++], LogLx-log2, value, false);
		}
	}
	if (boundY1 == 3) {
		z = MY+2;
		for (x=2; x<MX; x++) {
			LogLx = log(L[x]);
			addLntoLn(Gi[z], LogLx-log2, value, false);
			z+=MY;
		}
	}
	if (boundY2 == MY-2) {
		z = 2*MY-2;
		for (x=2; x<MX; x++) {
			LogLx = log(L[x]);
			addLntoLn(Gi[z], LogLx-log2, value, false);
			z+=MY;
		}
	}
	LogLx = log(L[2]);
	if (boundX1 == 3 && boundY1 == 3) {
		z = MY + 2;
		addLntoLn(Gi[z], LogLx+log4, value, true);
	}
	if (boundX1 == 3 && boundY2 == MY-2) {
		z = 2*MY-2;
		addLntoLn(Gi[z], LogLx+log4, value, true);
	}
	LogLx = log(L[MX-1]);
	if (boundX2 == MX-2 && boundY1 == 3) {
		z = MY*(MX-2)+2;
		addLntoLn(Gi[z], LogLx+log4, value, true);
	}
	if (boundX2 == MX-2 && boundY2 == MY-2) {
		z =(MX-1)*MY-2;
		addLntoLn(Gi[z], LogLx+log4, value, true);
	}
	return value;
}
