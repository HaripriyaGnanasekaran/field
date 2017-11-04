#include "Lat2DCyl1stOSten.h"

Lat2DCyl1stOSten::Lat2DCyl1stOSten(Input* MyInput_, Text name_)
	: Lat2DFlat(MyInput_,name_), Lat2DCylinder(MyInput_,name_),
	Lat2DCylinderSten(MyInput_,name_) {
	Gi1Prime.Dim(1,MY);
	Gi2Prime.Dim(1,MY);
}
Lat2DCyl1stOSten::~Lat2DCyl1stOSten() {
}

void
Lat2DCyl1stOSten::GetLatticeInfo(int* Info) const{ 
}


void
Lat2DCyl1stOSten::PropagateG(Matrix Gi, const Vector G, const int s) const {
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
			Gi[z][s] = G[z] * (e + l1*(d+f-e-e) + lambda1_1[x]*(b-e)
				+ lambda11[x]*(h-e) + lambda2_1[x]*(a+c-e-e)
				+ lambda21[x]*(g+i-e-e));
			zc++;
			z++;
			zi++;
		}
		z++;
	}
	UpdateBoundaries(Gi,s);
}

void
Lat2DCyl1stOSten::PropagateG(Matrix Gi, const Vector G, const int s, const double ff) const {
	int x,y,z,zc,zi,s_1;
	double a,b,c,d,e,f,g,h,i;
	double f_1= exp(-ff);
	double f1= exp(ff);
	double norm = 1+(l1+2*l2)*(f_1+f1-2);
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
			Gi[z][s] = G[z] * (e + l1*(d*f1 - e + f*f_1 - e) + 
				     lambda1_1[x]*(b - e)
				+    lambda11[x]*(h - e) 
				+ lambda2_1[x]*(a*f1 - e + c*f_1 - e)
				+ lambda21[x]*(g*f1 - e + i*f_1 - e))/norm;
			zc++;
			z++;
			zi++;
		}
		z++;
	}
	UpdateBoundaries(Gi,s);
}

void 
Lat2DCyl1stOSten::PropagateG(Vector Gi, const Vector G) const {
	int x,y,z,zc,zi;
	double a,b,c,d,e,f,g,h,i;

	for (y=2; y<=MY; y++) {
		Gi1Prime[y] = Gi[y];
	}
	zc = 2;
	z = MY+1;
	zi = MY+MY+1;
	for (x=2; x<MX; x++) {
		b = Gi[z-MY];
		c = Gi1Prime[zc];
		e = Gi[z++];
		f = Gi[z];
		h = Gi[zi++];
		i = Gi[zi++];
		Gi1Prime[zc++] = f;
		for (y=2; y<MY; y++) {
			a=b; b=c; c=Gi1Prime[zc];
			d=e; e=f; f=Gi[z+1];
			g=h; h=i; i=Gi[zi];
			Gi[z] = G[z] * (e + l1*(d+f-e-e) + lambda1_1[x]*(b-e)
				+ lambda11[x]*(h-e) + lambda2_1[x]*(a+c-e-e)
				+ lambda21[x]*(g+i-e-e));
			Gi1Prime[zc]=f;
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
Lat2DCyl1stOSten::PropagateG(Vector Gi, const Vector G, const double ff) const {
	int x,y,z,zc,zi;
	double a,b,c,d,e,f,g,h,i;
	double f_1= exp(-ff);
	double f1= exp(ff);
	double norm = 1+(l1+2*l2)*(f_1+f1-2);
	for (y=2; y<=MY; y++) {
		Gi1Prime[y] = Gi[y];
	}
	zc = 2;
	z = MY+1;
	zi = MY+MY+1;
	for (x=2; x<MX; x++) {
		b = Gi[z-MY];
		c = Gi1Prime[zc];
		e = Gi[z++];
		f = Gi[z];
		h = Gi[zi++];
		i = Gi[zi++];
		Gi1Prime[zc++] = f;
		for (y=2; y<MY; y++) {
			a=b; b=c; c=Gi1Prime[zc];
			d=e; e=f; f=Gi[z+1];
			g=h; h=i; i=Gi[zi];
			Gi[z] = G[z] * (e + l1*(d*f1 + f*f_1 - e - e) 
					+ lambda1_1[x]*(b - e)
					+ lambda11[x]*(h - e) 
					+ lambda2_1[x]*(a*f1 + c*f_1 - e - e)
					+ lambda21[x]*(g*f1 + i*f_1 - e - e))/norm;
			Gi1Prime[zc]=f;
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
Lat2DCyl1stOSten::PropagateG(const Vector Gi, const Vector G, Vector Gout) const {
	int x,y,z,zc,zi;
	double a,b,c,d,e,f,g,h,i;

	for (y=2; y<=MY; y++) {
		Gi1Prime[y] = Gi[y];
	}
	zc = 2;
	z = MY+1;
	zi = MY+MY+1;
	for (x=2; x<MX; x++) {
		b = Gi[z-MY];
		c = Gi1Prime[zc];
		e = Gi[z++];
		f = Gi[z];
		h = Gi[zi++];
		i = Gi[zi++];
		Gi1Prime[zc++] = f;
		for (y=2; y<MY; y++) {
			a=b; b=c; c=Gi1Prime[zc];
			d=e; e=f; f=Gi[z+1];
			g=h; h=i; i=Gi[zi];
			Gout[z] = G[z] * (e + l1*(d+f-e-e) + lambda1_1[x]*(b-e)
				+ lambda11[x]*(h-e) + lambda2_1[x]*(a+c-e-e)
				+ lambda21[x]*(g+i-e-e));
			Gi1Prime[zc]=f;
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
Lat2DCyl1stOSten::PropagateG(const Vector Gi, const Vector G, Vector Gout,const double ff) const {
	int x,y,z,zc,zi;
	double a,b,c,d,e,f,g,h,i;
	double f_1= exp(-ff);
	double f1= exp(ff);
	double norm = 1+(l1+2*l2)*(f_1+f1-2);
	for (y=2; y<=MY; y++) {
		Gi1Prime[y] = Gi[y];
	}
	zc = 2;
	z = MY+1;
	zi = MY+MY+1;
	for (x=2; x<MX; x++) {
		b = Gi[z-MY];
		c = Gi1Prime[zc];
		e = Gi[z++];
		f = Gi[z];
		h = Gi[zi++];
		i = Gi[zi++];
		Gi1Prime[zc++] = f;
		for (y=2; y<MY; y++) {
			a=b; b=c; c=Gi1Prime[zc];
			d=e; e=f; f=Gi[z+1];
			g=h; h=i; i=Gi[zi];
			Gout[z] = G[z] * (e + l1*(d*f1 + f*f_1 - e - e) 
					+ lambda1_1[x]*(b - e)
					+ lambda11[x]*(h - e) 
					+ lambda2_1[x]*(a*f1 + c*f_1 - e - e)
					+ lambda21[x]*(g*f1 + i*f_1 - e - e))/norm;
			Gi1Prime[zc]=f;
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
Lat2DCyl1stOSten::Propagate2G(Matrix Gi1, 
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
				Gi1[z][s] = G[z] * (e1 + l1*(d1+f1-e1-e1+d2+f2)
					+ lambda1_1[x]*(b1-e1+b2) + lambda11[x]*(h1-e1+h2)
					+ lambda2_1[x]*(a1+c1-e1-e1+a2+c2)
					+ lambda21[x]*(g1+i1-e1-e1+g2+i2));
				Gi2[z][s] = 0;
			} else {
				Gi1[z][s] = G[z] * (e1 + l1*(d1+f1-e1-e1) + lambda1_1[x]*(b1-e1)
					+ lambda11[x]*(h1-e1) + lambda2_1[x]*(a1+c1-e1-e1)
					+ lambda21[x]*(g1+i1-e1-e1));
				Gi2[z][s] = G[z] * (e2 + l1*(d2+f2-e2-e2) + lambda1_1[x]*(b2-e2)
					+ lambda11[x]*(h2-e2) + lambda2_1[x]*(a2+c2-e2-e2)
					+ lambda21[x]*(g2+i2-e2-e2));
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
Lat2DCyl1stOSten::Propagate2G(Matrix Gi1, 
						   Matrix Gi2, 
						   const Vector G, const int s,
						   const LatticeRange* LatRange,
						   const double f) const {
	int x,y,z,zc,zi,s_1;
	double a1,a2,b1,b2,c1,c2,d1,d2,e1,e2,f1,f2,g1,g2,h1,h2,i1,i2;
	double F_1= exp(-f);
	double F1= exp(f);
	double norm = 1+(l1+2*l2)*(F_1+F1-2);
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
				Gi1[z][s] = G[z] * (e1 + 
						l1*((d1 + d2)*F1  - e1 - e1 + (f2+ f1)*F_1)
						+ lambda1_1[x]*(b1-e1+b2) 
						+ lambda11[x]*(h1-e1+h2)
						+ lambda2_1[x]*((a1 + a2)*F1 - e1 - e1 + (c1 + c2)*F_1)
						+ lambda21[x]*((g1 + g2)*F1 - e1 - e1+ (i1 + i2)*F_1))/norm;
				Gi2[z][s] = 0;
			} else {
				Gi1[z][s] = G[z] * (e1 
						+ l1*(d1*F1 + f1*F_1 - e1 - e1) 
						+ lambda1_1[x]*(b1 - e1)
						+ lambda11[x]*(h1 - e1) 
						+ lambda2_1[x]*(a1*F1 + c1*F_1 - e1 - e1)
						+ lambda21[x]*(g1*F1 +i1*F_1 - e1 - e1))/norm;
				Gi2[z][s] = G[z] * (e2 
						+ l1*(d2*F1 + f2*F_1 - e2 - e2) 
						+ lambda1_1[x]*(b2 - e2)
						+ lambda11[x]*(h2 - e2) 
						+ lambda2_1[x]*(a2*F1 + c2*F_1 - e2 - e2)
						+ lambda21[x]*(g2*F1 + i2*F_1 - e2 - e2))/norm;
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
Lat2DCyl1stOSten::Propagate2G(Vector Gi1, 
						   Vector Gi2, 
						   const Vector G, 
						   const LatticeRange* LatRange) const {
	int x,y,z,zc,zi;
	double a1,a2,b1,b2,c1,c2,d1,d2,e1,e2,f1,f2,g1,g2,h1,h2,i1,i2;

	for (y=2; y<=MY; y++) {
		Gi1Prime[y] = Gi1[y];
		Gi2Prime[y] = Gi2[y];
	}
	zc = 2;
	z = MY+1;
	zi = MY+MY+1;
	for (x=2; x<MX; x++) {
		b1 = Gi1[z-MY];
		b2 = Gi2[z-MY];
		c1 = Gi1Prime[zc];
		c2 = Gi2Prime[zc];
		e1 = Gi1[z];
		e2 = Gi2[z++];
		f1 = Gi1[z];
		f2 = Gi2[z];
		h1 = Gi1[zi];
		h2 = Gi2[zi++];
		i1 = Gi1[zi];
		i2 = Gi2[zi++];
		Gi1Prime[zc] = f1;
		Gi2Prime[zc++] = f2;
		for (y=2; y<MY; y++) {
			a1=b1; b1=c1; c1=Gi1Prime[zc];
			a2=b2; b2=c2; c2=Gi2Prime[zc];
			d1=e1; e1=f1; f1=Gi1[z+1];
			d2=e2; e2=f2; f2=Gi2[z+1];
			g1=h1; h1=i1; i1=Gi1[zi];
			g2=h2; h2=i2; i2=Gi2[zi];
			if (LatRange->InRange(z)) {
				Gi1[z] = G[z] * (e1 + l1*(d1+f1-e1-e1+d2+f2)
					+ lambda1_1[x]*(b1-e1+b2) + lambda11[x]*(h1-e1+h2)
					+ lambda2_1[x]*(a1+c1-e1-e1+a2+c2)
					+ lambda21[x]*(g1+i1-e1-e1+g2+i2));
				Gi2[z] = 0;
			} else {
				Gi1[z] = G[z] * (e1 + l1*(d1+f1-e1-e1) + lambda1_1[x]*(b1-e1)
					+ lambda11[x]*(h1-e1) + lambda2_1[x]*(a1+c1-e1-e1)
					+ lambda21[x]*(g1+i1-e1-e1));
				Gi2[z] = G[z] * (e2 + l1*(d2+f2-e2-e2) + lambda1_1[x]*(b2-e2)
					+ lambda11[x]*(h2-e2) + lambda2_1[x]*(a2+c2-e2-e2)
					+ lambda21[x]*(g2+i2-e2-e2));
			}
			Gi1Prime[zc]=f1;
			Gi2Prime[zc]=f2;
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
Lat2DCyl1stOSten::Propagate2G(Vector Gi1, 
						   Vector Gi2, 
						   const Vector G, 
						   const LatticeRange* LatRange,
						   const double f) const {
	int x,y,z,zc,zi;
	double a1,a2,b1,b2,c1,c2,d1,d2,e1,e2,f1,f2,g1,g2,h1,h2,i1,i2;
	double F_1= exp(-f);
	double F1= exp(f);
	double norm = 1+(l1+2*l2)*(F_1+F1-2);
	for (y=2; y<=MY; y++) {
		Gi1Prime[y] = Gi1[y];
		Gi2Prime[y] = Gi2[y];
	}
	zc = 2;
	z = MY+1;
	zi = MY+MY+1;
	for (x=2; x<MX; x++) {
		b1 = Gi1[z-MY];
		b2 = Gi2[z-MY];
		c1 = Gi1Prime[zc];
		c2 = Gi2Prime[zc];
		e1 = Gi1[z];
		e2 = Gi2[z++];
		f1 = Gi1[z];
		f2 = Gi2[z];
		h1 = Gi1[zi];
		h2 = Gi2[zi++];
		i1 = Gi1[zi];
		i2 = Gi2[zi++];
		Gi1Prime[zc] = f1;
		Gi2Prime[zc++] = f2;
		for (y=2; y<MY; y++) {
			a1=b1; b1=c1; c1=Gi1Prime[zc];
			a2=b2; b2=c2; c2=Gi2Prime[zc];
			d1=e1; e1=f1; f1=Gi1[z+1];
			d2=e2; e2=f2; f2=Gi2[z+1];
			g1=h1; h1=i1; i1=Gi1[zi];
			g2=h2; h2=i2; i2=Gi2[zi];
			if (LatRange->InRange(z)) {
				Gi1[z] = G[z] * (e1 + 
						l1*((d1 + d2)*F1- e1 - e1  + (f1 + f2)*F_1)
						+ lambda1_1[x]*(b1 - e1 +b2) 
						+ lambda11[x]*(h1 - e1 + h2)
						+ lambda2_1[x]*((a1 + a2)*F1 - e1 - e1 + (c1 + c2)*F_1)
						+ lambda21[x]*((g1 + g2)*F1 - e1 - e1 + (i1 + i2)*F_1))/norm;
				Gi2[z] = 0;
			} else {
				Gi1[z] = G[z] * (e1 
						+ l1*(d1*F1 + f1*F_1 - e1 - e1) 
						+ lambda1_1[x]*(b1 - e1)
						+ lambda11[x]*(h1 - e1) 
						+ lambda2_1[x]*(a1*F1 + c1*F_1 - e1 - e1)
						+ lambda21[x]*(g1*F1 + i1*F_1 - e1 - e1))/norm;
				Gi2[z] = G[z] * (e2 
						+ l1*(d2*F1 + f2*F_1 - e2 - e2) 
						+ lambda1_1[x]*(b2 - e2)
						+ lambda11[x]*(h2 - e2) 
						+ lambda2_1[x]*(a2*F1 + c2*F_1 - e2 - e2)
						+ lambda21[x]*(g2*F1 + i2*F_1 - e2 - e2))/norm;
			}
			Gi1Prime[zc]=f1;
			Gi2Prime[zc]=f2;
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
