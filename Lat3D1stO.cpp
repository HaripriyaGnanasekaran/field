#include "Lat3D1stO.h"
#include "tools.h"
#include <iostream>


// static Vectors for the Propagate functions. This avoid allocating memory
// at each call to a Propagate function, which is very time intensive
Vector Lat3D1stO::Gx;
Vector Lat3D1stO::Gt;
Vector Lat3D1stO::Gp;
Vector Lat3D1stO::G2;


Lat3D1stO::Lat3D1stO(Input* MyInput_, Text name_)
	: Lat3D(MyInput_,name_) {
	Gx.Dim(1,Mmax);
	Gt.Dim(1,Mmax);
	Gp.Dim(1,Mmax);
	G2.Dim(1,Mmax);
}

Lat3D1stO::~Lat3D1stO() {
}
Boolean
Lat3D1stO::OverflowProtection(void) const {
	return false;
}
void
Lat3D1stO::MakeSafe(Vector) const {
}

void
Lat3D1stO::RestoreFromSafe(Vector) const {
}

void
Lat3D1stO::GetLatticeInfo(int* Info) const{ 
	Info[0]=3; //Dimensions
	Info[1]=1; //simpel cubic;
	if (N_comp_ranges>1) {Message(fatal,"GPU not implemented for more than one comp-range. ");} 
	Info[2]=CR_Q[1]->GetNumLayers_x();
	Info[3]=CR_Q[1]->GetNumLayers_y();
	Info[4]=CR_Q[1]->GetNumLayers_z(); 
	//assumption that compranges =1;	
	Info[5]=Get_BL(1); //x lower bound
	Info[6]=Get_BL(2); //x upper bound
	Info[7]=Get_BL(3); //y lower bound
	Info[8]=Get_BL(4); //y upper bound
	Info[9]=Get_BL(5); //z lower bound
	Info[10]=Get_BL(6);//z upper bound
}

void
Lat3D1stO::PropagateG(const Vector Gi, const Vector G, Vector Gout) const{

	int i,jx,jy,M,Mprev=0;
	double *gout = &Gout[1];
	double *gs   = &Gx[1];
	const double *pgi=&Gi[1];
	const double *pg = &G[1];
	double *gs_1 = const_cast<double*>(pgi);
	double *g = const_cast<double*>(pg);
	//const double *gs_1 = &Gi[1];
	//const double *g = &G[1];
	SetBx1(gs_1);SetBxm(gs_1);SetBy1(gs_1);SetBym(gs_1);SetBz1(gs_1);SetBzm(gs_1);

	for (i=1; i<= N_comp_ranges; i++){
		M=CR_Q[i]->Get_M();
		if (i>1) Mprev +=CR_Q[i-1]->Get_M();
		jx=CR_Q[i]->Get_jx();
		jy=CR_Q[i]->Get_jy();
		gs_1 +=Mprev;
		g    +=Mprev;
		gs   +=Mprev;
		gout +=Mprev;

           	times(gs+jx,gs_1,   g+jx,M-jx);
		addTimes(gs,   gs_1+jx,g,   M-jx);
		addTimes(gs+jy,gs_1,   g+jy,M-jy);
        	addTimes(gs,   gs_1+jy,g,   M-jy);
        	addTimes(gs+jz,gs_1,   g+jz,M-jz);
        	addTimes(gs,   gs_1+jz,g,   M-jz);
       	timesC(gout,gs,1.0/6.0,M);

		gs_1-=Mprev;
		g   -=Mprev;
		gs  -=Mprev;
		gout-=Mprev;
	}
}


void
Lat3D1stO::PropagateG(const Vector Gi, const Vector G, Vector Gout, const double f) const{
	const double f1 = exp(f);
	const double f_1 = exp(-f);
	double fnorm = 1.0/(4.0 + exp(f) + exp(-f));
	int i,jx,jy,M,Mprev=0;
	double *gout = &Gout[1];
	double *gs   = &Gx[1];
	const double *pgi=&Gi[1];
	const double *pg = &G[1];
	double *gs_1 = const_cast<double*>(pgi);
	double *g = const_cast<double*>(pg);
	//const double *gs_1 = &Gi[1];
	//const double *g = &G[1];
	SetBx1(gs_1);SetBxm(gs_1);SetBy1(gs_1);SetBym(gs_1);SetBz1(gs_1);SetBzm(gs_1);

	for (i=1; i<= N_comp_ranges; i++){
		M=CR_Q[i]->Get_M();
		if (i>1) Mprev +=CR_Q[i-1]->Get_M();
		jx=CR_Q[i]->Get_jx();
		jy=CR_Q[i]->Get_jy();
		gs_1 +=Mprev;
		g    +=Mprev;
		gs   +=Mprev;
		gout +=Mprev;

           	times(gs+jx,gs_1,   g+jx,M-jx);
		addTimes(gs,   gs_1+jx,g,   M-jx);
		addTimes(gs+jy,gs_1,   g+jy,M-jy);
        	addTimes(gs,   gs_1+jy,g,   M-jy);
       	addTimesF(gs+jz,gs_1,   g+jz,f1,M-jz);
       	addTimesF(gs,   gs_1+jz,g,   f_1, M-jz);
       	timesC(gout,gs,fnorm,M);

		gs_1-=Mprev;
		g   -=Mprev;
		gs  -=Mprev;
		gout-=Mprev;
	}
}

void
Lat3D1stO::PropagateG(Matrix Gi, const Vector G, const int s) const {
	int i,jx,jy,M,Mprev=0;
	double *gs   = &Gi[1][s];
	double *gs_1 = &Gi[1][s-1];
	const double *g = &G[1];
	SetBx1(gs_1);SetBxm(gs_1);SetBy1(gs_1);SetBym(gs_1);SetBz1(gs_1);SetBzm(gs_1);

	for (i=1; i<= N_comp_ranges; i++){
		M=CR_Q[i]->Get_M();
		if (i>1) Mprev +=CR_Q[i-1]->Get_M();
		jx=CR_Q[i]->Get_jx();
		jy=CR_Q[i]->Get_jy();
		gs_1+=Mprev;
		g +=Mprev;
		gs+=Mprev;

		   times(gs+jx,gs_1,   g+jx,M-jx);
        	addTimes(gs   ,gs_1+jx,g,   M-jx);
		addTimes(gs+jy,gs_1,   g+jy,M-jy);
		addTimes(gs,   gs_1+jy,g,   M-jy);
        	addTimes(gs+jz,gs_1,   g+jz,M-jz);
       	addTimes(gs,   gs_1+jz,g,   M-jz);
       	norm(gs,1.0/6.0,M);

        	gs_1-=Mprev;
       	g   -=Mprev;
        	gs  -=Mprev;
	}
}

void
Lat3D1stO::PropagateG(Matrix Gi, const Vector G, const int s, const double f) const {
	const double f1 = exp(f);
	const double f_1 = exp(-f);
	double fnorm = 1.0/(4.0 + exp(f) + exp(-f));
	int i,jx,jy,M,Mprev=0;
	double *gs   = &Gi[1][s];
	double *gs_1 = &Gi[1][s-1];
	const double *g = &G[1];
	SetBx1(gs_1);SetBxm(gs_1);SetBy1(gs_1);SetBym(gs_1);SetBz1(gs_1);SetBzm(gs_1);

	for (i=1; i<= N_comp_ranges; i++){
		M=CR_Q[i]->Get_M();
		if (i>1) Mprev +=CR_Q[i-1]->Get_M();
		jx=CR_Q[i]->Get_jx();
		jy=CR_Q[i]->Get_jy();
		gs_1+=Mprev;
		g +=Mprev;
		gs+=Mprev;

		   times(gs+jx,gs_1,   g+jx,M-jx);
        	addTimes(gs   ,gs_1+jx,g,   M-jx);
		addTimes(gs+jy,gs_1,   g+jy,M-jy);
		addTimes(gs,   gs_1+jy,g,   M-jy);
        	addTimesF(gs+jz,gs_1,   g+jz,f1, M-jz);
        	addTimesF(gs,   gs_1+jz,g,   f_1,  M-jz);
        	norm(gs,fnorm,M);

        	gs_1-=Mprev;
        	g   -=Mprev;
        	gs  -=Mprev;
	}

}

void
Lat3D1stO::PropagateG(Vector Gi, const Vector G) const {
	int i,jx,jy,M,Mprev=0;
	double *gs   = &Gx[1];
	double *gs_1 = &Gi[1];
	const double *g = &G[1];
	SetBx1(gs_1);SetBxm(gs_1);SetBy1(gs_1);SetBym(gs_1);SetBz1(gs_1);SetBzm(gs_1);

	for (i=1; i<= N_comp_ranges; i++){
		M=CR_Q[i]->Get_M();
		if (i>1) Mprev +=CR_Q[i-1]->Get_M();
		jx=CR_Q[i]->Get_jx();
		jy=CR_Q[i]->Get_jy();
		gs_1+=Mprev;
		g +=Mprev;

        	times(gs+jx,gs_1,   g+jx,M-jx);
		addTimes(gs,   gs_1+jx,g,   M-jx);
		addTimes(gs+jy,gs_1,   g+jy,M-jy);
       	addTimes(gs,   gs_1+jy,g,   M-jy);
        	addTimes(gs+jz,gs_1,   g+jz,M-jz);
       	addTimes(gs,   gs_1+jz,g,   M-jz);
       	timesC(gs_1,gs,1.0/6.0,M);

		gs_1-=Mprev;
		g -=Mprev;
	}
}

void
Lat3D1stO::PropagateG(Vector Gi, const Vector G, const double f) const {
	const double f1 = exp(f);
	const double f_1 = exp(-f);
	double fnorm = 1.0/(4.0 + exp(f) + exp(-f));
	int i,jx,jy,M,Mprev=0;
	double *gs   = &Gx[1];
	double *gs_1 = &Gi[1];
	const double *g = &G[1];
	SetBx1(gs_1);SetBxm(gs_1);SetBy1(gs_1);SetBym(gs_1);SetBz1(gs_1);SetBzm(gs_1);

	for (i=1; i<= N_comp_ranges; i++){
		M=CR_Q[i]->Get_M();
		if (i>1) Mprev +=CR_Q[i-1]->Get_M();
		jx=CR_Q[i]->Get_jx();
		jy=CR_Q[i]->Get_jy();
		gs_1+=Mprev;
		g +=Mprev;

        	times(gs+jx,gs_1,   g+jx,M-jx);
		addTimes(gs,   gs_1+jx,g,   M-jx);
		addTimes(gs+jy,gs_1,   g+jy,M-jy);
        	addTimes(gs,   gs_1+jy,g,   M-jy);
        	addTimesF(gs+jz,gs_1,   g+jz,f1,M-jz);
        	addTimesF(gs,   gs_1+jz,g,   f_1, M-jz);
        	timesC(gs_1,gs,fnorm,M);

		gs_1-=Mprev;
		g -=Mprev;
	}

}

void
Lat3D1stO::Init2G(Vector Gi1, Vector Gi2, const Vector G, const LatticeRange* LatRange) const {
	int i,z,M,Mprev=0;
	for (i=1; i<=N_comp_ranges; i++) {
		M=CR_Q[i]->Get_M();
		if (i>1) Mprev +=CR_Q[i-1]->Get_M();
		for (z=1; z<=M; z++) {
			if (LatRange->InRange(Mprev+z)) {
				Gi1[Mprev+z] = G[Mprev+z];
				Gi2[Mprev+z] = 0;
			} else {
				Gi1[Mprev+z] = 0;
				Gi2[Mprev+z] = G[Mprev+z];
			}
		}
	}
}

void
Lat3D1stO::Propagate2G(Matrix Gi1, Matrix Gi2, const Vector G, const int s,const LatticeRange* LatRange) const {
	Message(fatal,"Propagate2G not implemented in 3d");
}
void
Lat3D1stO::Propagate2G(Matrix Gi1, Matrix Gi2, const Vector G, const int s,const LatticeRange* LatRange,const double f) const {
	Message(fatal,MyInput,"Force ensemble not implemented in Lat3D1stO");
}

void
Lat3D1stO::Propagate2G(Vector Gi1, Vector Gi2, const Vector G, const LatticeRange* LatRange) const {
	Message(fatal,"Propagate2G not implemented in 3d");
}

void
Lat3D1stO::Propagate2G(Vector Gi1, Vector Gi2, const Vector G, const LatticeRange* LatRange, const double f) const {
	Message(fatal,MyInput,"Force ensemble not implemented in Lat3D1stO");
}

void
Lat3D1stO::ConnectG(const Vector GiA, const Matrix GiB, const int s, Vector out) const {
	const double *p_giB=&GiB[1][s];
	const double *p_giA = &GiA[1];
	double *p_out = &out[1];
	addTimes(p_out,p_giA,p_giB,GetTotalNumLayers());
}
void
Lat3D1stO::ConnectG(const Vector GiA, const Vector GiB, Vector out) const {
	const double *p_giA = &GiA[1];
	const double *p_giB = &GiB[1];
	double *p_out = &out[1];
	addTimes(p_out,p_giA,p_giB,GetTotalNumLayers());
}
void
Lat3D1stO::Connect2G(const Vector GiA1, const Matrix GiB1, const int s1, const Vector GiA2,
						 const Matrix GiB2, const int s2, Vector Out) const {
	const double *giA1 = &GiA1[1];
	const double *giB1 = &GiB1[1][s2];
	const double *giA2 = &GiA2[1];
	const double *giB2 = &GiB2[1][s1];
	double *out =&Out[1];

	addTimes(out,giA1,giB2,GetTotalNumLayers());
	addTimes(out,giA2,giB1,GetTotalNumLayers());
}

void
Lat3D1stO::Connect2G(const Vector GiA1, const Vector GiB1, const Vector GiA2,  const Vector GiB2, Vector Out) const {
	const double *giA1 = &GiA1[1];
	const double *giB1 = &GiB1[1];
	const double *giA2 = &GiA2[1];
	const double *giB2 = &GiB2[1];
	double *out =&Out[1];

	addTimes(out,giA1,giB2,GetTotalNumLayers());
	addTimes(out,giA2,giB1,GetTotalNumLayers());
}

void
Lat3D1stO::CorrectDoubleCountG(Vector in, const Vector G) const {
	double *p_in = &in[1];
	const double *p_G = &G[1];
	div(p_in,p_G,GetTotalNumLayers());
}
double
Lat3D1stO::ComputeLnGN(Vector Gi) const {
	double *pGi = &Gi[1];
	//Message(debug,"ComputeLnGN(Vector Gi)");
	SubtractBoundaries(pGi);
	int i,x,y,z,px,py,NlayersX,NlayersY,NlayersZ;
	//int Previous=0;
	int Current=0,Next=0,jx,jy;
	double value = 0;
	for (i=1; i<=N_comp_ranges; i++) {
		//Previous=Current; 
		Current=Next; Next +=CR_Q[i]->Get_M();
		NlayersX=CR_Q[i]->GetNumLayers_x();
		NlayersY=CR_Q[i]->GetNumLayers_y();
		NlayersZ=CR_Q[i]->GetNumLayers_z();
		jx=CR_Q[i]->Get_jx();
		jy=CR_Q[i]->Get_jy();
		px=0;
		for (x=1; x<=NlayersX; x++) {
			px+=jx; py=0; for (y=1; y<=NlayersY; y++) {
				py+=jy; for (z=1; z<=NlayersZ; z++) value+=pGi[Current+px+py+z];
			}
		}
	}
	RestoreBoundaries(pGi);
	return log(value);
}
void
Lat3D1stO::NormPhiFree(Vector phi, const double C) const {
	//Message(debug,"NormPhiFree(Vector phi, const double C)");
	norm(phi,C,GetTotalNumLayers());
}
void
Lat3D1stO::NormPhiRestr(Vector phi, const Vector Gi, double C) const {
	//Message(debug,"NormPhiRestr(Vector phi)");
	C /= exp(ComputeLnGN(Gi));
	norm(phi,C,GetTotalNumLayers());
}
void
Lat3D1stO::UpdateBoundaries(Vector A) const {
	//Message(debug,"UpdateBoundaries(Vector A)");
	double *pA=&A[1];
	SetBoundaries(pA);
}
void
Lat3D1stO::UpdateBoundaries(Matrix A, const int s) const {
	double *pA=&A[1][s];
	SetBoundaries(pA);
}
