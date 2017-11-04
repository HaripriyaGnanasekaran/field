#include "Lat3D2ndO.h"
#include "tools.h"
#include <iostream>



// static Vectors for the Propagate functions. This avoid allocating memory
// at each call to a Propagate function, which is very time intensive
Vector Lat3D2ndO::Gx0;
Vector Lat3D2ndO::Gx1;
Vector Lat3D2ndO::Gx2;
Vector Lat3D2ndO::Gx3;
Vector Lat3D2ndO::Gx4;
Vector Lat3D2ndO::Gx5;
Vector Lat3D2ndO::Pl;
Vector Lat3D2ndO::Ps;
Vector Lat3D2ndO::GL;
Vector Lat3D2ndO::GS;


Lat3D2ndO::Lat3D2ndO(Input* MyInput_, Text name_) : Lat3D(MyInput_,name_) {
	Gx0.Dim(1,Mmax);
	Gx1.Dim(1,Mmax);
	Gx2.Dim(1,Mmax);
	Gx3.Dim(1,Mmax);
	Gx4.Dim(1,Mmax);
	Gx5.Dim(1,Mmax);
	Pl.Dim(1,GetTotalNumLayers());
	Ps.Dim(1,GetTotalNumLayers());
	GL.Dim(1,GetTotalNumLayers());
	GS.Dim(1,GetTotalNumLayers());
}

Lat3D2ndO::~Lat3D2ndO() {
}
Boolean
Lat3D2ndO::OverflowProtection(void) const {
	return false;
}
void
Lat3D2ndO::MakeSafe(Vector) const {
}
void
Lat3D2ndO::RestoreFromSafe(Vector) const {
}

void
Lat3D2ndO::GetLatticeInfo(int *Info) const {
	Info[0]=3; //Dimensions
	Info[1]=1; //simpel cubic;
	if (N_comp_ranges>1) {Message(fatal,"GPU not implemented for more than one comp-range. Deactivate GPU for molecule. ");}
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
Lat3D2ndO::PropagateF(Matrix Gi, Vector G, const int s, const double S) const {
	int i,jx,jy, M,Mprev=0;
	M=GetTotalNumLayers();

	double L=(1.0-S)/4.0;
	double *gs0 = &Gi[1][s];
	double *gs1 = &Gi[1+M][s];
	double *gs2 = &Gi[1+2*M][s];
	double *gs3 = &Gi[1+3*M][s];
	double *gs4 = &Gi[1+4*M][s];
	double *gs5 = &Gi[1+5*M][s];

	double *gz0 = &Gi[1][s-1];
	double *gz1 = &Gi[1+M][s-1];
	double *gz2 = &Gi[1+2*M][s-1];
	double *gz3 = &Gi[1+3*M][s-1];
	double *gz4 = &Gi[1+4*M][s-1];
	double *gz5 = &Gi[1+5*M][s-1];

	SetBx1(gz0,gz5);SetBx1(gz1);SetBx1(gz2);SetBx1(gz3);SetBx1(gz4);
	SetBy1(gz1,gz4);SetBy1(gz0);SetBy1(gz2);SetBy1(gz3);SetBy1(gz5);
	SetBz1(gz2,gz3);SetBz1(gz0);SetBz1(gz1);SetBz1(gz4);SetBz1(gz5);
	SetBxm(gz0,gz5);SetBxm(gz1);SetBxm(gz2);SetBxm(gz3);SetBxm(gz4);
	SetBym(gz1,gz4);SetBym(gz0);SetBym(gz2);SetBym(gz3);SetBym(gz5);
	SetBzm(gz2,gz3);SetBzm(gz0);SetBzm(gz1);SetBzm(gz4);SetBzm(gz5);

	double *gx0 = &Gx0[1];
	double *gx1 = &Gx1[1];
	double *gx2 = &Gx2[1];
	double *gx3 = &Gx3[1];
	double *gx4 = &Gx4[1];
	double *gx5 = &Gx5[1];

	double *g = &G[1];
	double *gL = &GL[1]; timesC(gL,g,L,M);
	double *gS = &GS[1]; timesC(gS,g,S,M);
	for (i=1; i<= N_comp_ranges; i++){
		M=CR_Q[i]->Get_M();
		if (i>1) Mprev +=CR_Q[i-1]->Get_M();
		jx=CR_Q[i]->Get_jx();
		jy=CR_Q[i]->Get_jy();
		gz0+=Mprev; gz1+=Mprev; gz2+=Mprev; gz3+=Mprev; gz4+=Mprev; gz5+=Mprev;
		gs0+=Mprev; gs1+=Mprev; gs2+=Mprev; gs3+=Mprev; gs4+=Mprev; gs5+=Mprev;
		gS+=Mprev; gL+=Mprev;
		times(gx0+jx,gz0,gS+jx,M-jx);
		addTimes(gx0+jx,gz1,gL+jx,M-jx);
		addTimes(gx0+jx,gz2,gL+jx,M-jx);
		addTimes(gx0+jx,gz3,gL+jx,M-jx);
		addTimes(gx0+jx,gz4,gL+jx,M-jx);

        times(gx1+jy,gz0,gL+jy,M-jy);
		addTimes(gx1+jy,gz1,gS+jy,M-jy);
		addTimes(gx1+jy,gz2,gL+jy,M-jy);
		addTimes(gx1+jy,gz3,gL+jy,M-jy);
		addTimes(gx1+jy,gz5,gL+jy,M-jy);

        times(gx2+jz,gz0,gL+jz,M-jz);
		addTimes(gx2+jz,gz1,gL+jz,M-jz);
		addTimes(gx2+jz,gz2,gS+jz,M-jz);
		addTimes(gx2+jz,gz4,gL+jz,M-jz);
		addTimes(gx2+jz,gz5,gL+jz,M-jz);

        times(gx3,gz0+jz,gL,M-jz);
		addTimes(gx3,gz1+jz,gL,M-jz);
		addTimes(gx3,gz3+jz,gS,M-jz);
		addTimes(gx3,gz4+jz,gL,M-jz);
		addTimes(gx3,gz5+jz,gL,M-jz);

        times(gx4,gz0+jy,gL,M-jy);
		addTimes(gx4,gz2+jy,gL,M-jy);
		addTimes(gx4,gz3+jy,gL,M-jy);
		addTimes(gx4,gz4+jy,gS,M-jy);
		addTimes(gx4,gz5+jy,gL,M-jy);

        times(gx5,gz1+jx,gL,M-jx);
		addTimes(gx5,gz2+jx,gL,M-jx);
		addTimes(gx5,gz3+jx,gL,M-jx);
		addTimes(gx5,gz4+jx,gL,M-jx);
		addTimes(gx5,gz5+jx,gS,M-jx);
		cp(gs0,gx0,M);cp(gs1,gx1,M);cp(gs2,gx2,M);cp(gs3,gx3,M);cp(gs4,gx4,M);cp(gs5,gx5,M);

		gz0-=Mprev; gz1-=Mprev; gz2-=Mprev; gz3-=Mprev; gz4-=Mprev; gz5-=Mprev;
		gs0-=Mprev; gs1-=Mprev; gs2-=Mprev; gs3-=Mprev; gs4-=Mprev; gs5-=Mprev;
		gS-=Mprev;  gL-=Mprev;
	}
	//double sums_1=0, sums=0;
	//for (int i = 0; i<6*M; i++) { sums_1+=Gi[i+1][s-1]; sums +=Gi[i+1][s];} cout<< s << " FgiGSs " << sums_1 << "  " << sums << endl;
}

void
Lat3D2ndO::PropagateF(Matrix Gi, Vector G, const int s, const double S, const bool stiff_range) const {
	int i,jx,jy, M,Mprev=0;
	M=GetTotalNumLayers();
	double *pl=&Pl[1];
	double *ps=&Ps[1];
	for (int i=0; i<M; i++) {
		if (stiff->InRange(i+1)) {ps[i]=S;} else {ps[i]=0.2; }  pl[i]=(1.0-ps[i])/4.0;
	}

	double *gs0 = &Gi[1][s];
	double *gs1 = &Gi[1+M][s];
	double *gs2 = &Gi[1+2*M][s];
	double *gs3 = &Gi[1+3*M][s];
	double *gs4 = &Gi[1+4*M][s];
	double *gs5 = &Gi[1+5*M][s];

	double *gz0 = &Gi[1][s-1];
	double *gz1 = &Gi[1+M][s-1];
	double *gz2 = &Gi[1+2*M][s-1];
	double *gz3 = &Gi[1+3*M][s-1];
	double *gz4 = &Gi[1+4*M][s-1];
	double *gz5 = &Gi[1+5*M][s-1];

	SetBx1(gz0,gz5);SetBx1(gz1);SetBx1(gz2);SetBx1(gz3);SetBx1(gz4);
	SetBy1(gz1,gz4);SetBy1(gz0);SetBy1(gz2);SetBy1(gz3);SetBy1(gz5);
	SetBz1(gz2,gz3);SetBz1(gz0);SetBz1(gz1);SetBz1(gz4);SetBz1(gz5);
	SetBxm(gz0,gz5);SetBxm(gz1);SetBxm(gz2);SetBxm(gz3);SetBxm(gz4);
	SetBym(gz1,gz4);SetBym(gz0);SetBym(gz2);SetBym(gz3);SetBym(gz5);
	SetBzm(gz2,gz3);SetBzm(gz0);SetBzm(gz1);SetBzm(gz4);SetBzm(gz5);

	double *gx0 = &Gx0[1];
	double *gx1 = &Gx1[1];
	double *gx2 = &Gx2[1];
	double *gx3 = &Gx3[1];
	double *gx4 = &Gx4[1];
	double *gx5 = &Gx5[1];

	double *g = &G[1];

	for (i=1; i<= N_comp_ranges; i++){
		M=CR_Q[i]->Get_M();
		if (i>1) Mprev +=CR_Q[i-1]->Get_M();
		jx=CR_Q[i]->Get_jx();
		jy=CR_Q[i]->Get_jy();
		gz0+=Mprev; gz1+=Mprev; gz2+=Mprev; gz3+=Mprev; gz4+=Mprev; gz5+=Mprev;
		gs0+=Mprev; gs1+=Mprev; gs2+=Mprev; gs3+=Mprev; gs4+=Mprev; gs5+=Mprev;
		g+=Mprev; ps+=Mprev; pl+=Mprev;


		times(gx0+jx,gz0,g+jx,ps,M-jx);
		addTimes(gx0+jx,gz1,g+jx,pl,M-jx);
		addTimes(gx0+jx,gz2,g+jx,pl,M-jx);
		addTimes(gx0+jx,gz3,g+jx,pl,M-jx);
		addTimes(gx0+jx,gz4,g+jx,pl,M-jx);

        times(gx1+jy,gz0,g+jy,pl,M-jy);
		addTimes(gx1+jy,gz1,g+jy,ps,M-jy);
		addTimes(gx1+jy,gz2,g+jy,pl,M-jy);
		addTimes(gx1+jy,gz3,g+jy,pl,M-jy);
		addTimes(gx1+jy,gz5,g+jy,pl,M-jy);

        times(gx2+jz,gz0,g+jz,pl,M-jz);
		addTimes(gx2+jz,gz1,g+jz,pl,M-jz);
		addTimes(gx2+jz,gz2,g+jz,ps,M-jz);
		addTimes(gx2+jz,gz4,g+jz,pl,M-jz);
		addTimes(gx2+jz,gz5,g+jz,pl,M-jz);

        times(gx3,gz0+jz,g,pl,M-jz);
		addTimes(gx3,gz1+jz,g,pl,M-jz);
		addTimes(gx3,gz3+jz,g,ps,M-jz);
		addTimes(gx3,gz4+jz,g,pl,M-jz);
		addTimes(gx3,gz5+jz,g,pl,M-jz);

        times(gx4,gz0+jy,g,pl,M-jy);
		addTimes(gx4,gz2+jy,g,pl,M-jy);
		addTimes(gx4,gz3+jy,g,pl,M-jy);
		addTimes(gx4,gz4+jy,g,ps,M-jy);
		addTimes(gx4,gz5+jy,g,pl,M-jy);

        times(gx5,gz1+jx,g,pl,M-jx);
		addTimes(gx5,gz2+jx,g,pl,M-jx);
		addTimes(gx5,gz3+jx,g,pl,M-jx);
		addTimes(gx5,gz4+jx,g,pl,M-jx);
		addTimes(gx5,gz5+jx,g,ps,M-jx);


		cp(gs0,gx0,M);cp(gs1,gx1,M);cp(gs2,gx2,M);cp(gs3,gx3,M);cp(gs4,gx4,M);cp(gs5,gx5,M);

		gz0-=Mprev; gz1-=Mprev; gz2-=Mprev; gz3-=Mprev; gz4-=Mprev; gz5-=Mprev;
		gs0-=Mprev; gs1-=Mprev; gs2-=Mprev; gs3-=Mprev; gs4-=Mprev; gs5-=Mprev;
		g-=Mprev; ps-=Mprev; pl-=Mprev;
	}
}

void
Lat3D2ndO::PropagateF(Matrix Gi, Vector G, const int s) const {
	int i,jx,jy, M,Mprev=0;
	M=GetTotalNumLayers();

	double *gs0 = &Gi[1][s];
	double *gs1 = &Gi[1+M][s];
	double *gs2 = &Gi[1+2*M][s];
	double *gs3 = &Gi[1+3*M][s];
	double *gs4 = &Gi[1+4*M][s];
	double *gs5 = &Gi[1+5*M][s];

	double *gz0 = &Gi[1][s-1];
	double *gz1 = &Gi[1+M][s-1];
	double *gz2 = &Gi[1+2*M][s-1];
	double *gz3 = &Gi[1+3*M][s-1];
	double *gz4 = &Gi[1+4*M][s-1];
	double *gz5 = &Gi[1+5*M][s-1];

	SetBx1(gz0,gz5);SetBx1(gz1);SetBx1(gz2);SetBx1(gz3);SetBx1(gz4);
	SetBy1(gz1,gz4);SetBy1(gz0);SetBy1(gz2);SetBy1(gz3);SetBy1(gz5);
	SetBz1(gz2,gz3);SetBz1(gz0);SetBz1(gz1);SetBz1(gz4);SetBz1(gz5);
	SetBxm(gz0,gz5);SetBxm(gz1);SetBxm(gz2);SetBxm(gz3);SetBxm(gz4);
	SetBym(gz1,gz4);SetBym(gz0);SetBym(gz2);SetBym(gz3);SetBym(gz5);
	SetBzm(gz2,gz3);SetBzm(gz0);SetBzm(gz1);SetBzm(gz4);SetBzm(gz5);

	double *gx0   = &Gx0[1];
	double *gx1   = &Gx1[1];
	double *gx2   = &Gx2[1];
	double *gx3   = &Gx3[1];
	double *gx4   = &Gx4[1];
	double *gx5   = &Gx5[1];

	double *g = &G[1];
	double *gI = &GL[1]; timesC(gI,g,1.0/6.0,M); //isotropic GL used because this vector was already generated...

	for (i=1; i<= N_comp_ranges; i++){
		M=CR_Q[i]->Get_M();
		if (i>1) Mprev +=CR_Q[i-1]->Get_M();
		jx=CR_Q[i]->Get_jx();
		jy=CR_Q[i]->Get_jy();
		gz0+=Mprev; gz1+=Mprev; gz2+=Mprev; gz3+=Mprev; gz4+=Mprev; gz5+=Mprev;
		gs0+=Mprev; gs1+=Mprev; gs2+=Mprev; gs3+=Mprev; gs4+=Mprev; gs5+=Mprev;
		gI+=Mprev;

		   times(gx0+jx,gz0,gI+jx,M-jx);
		addTimes(gx0+jx,gz1,gI+jx,M-jx);
		addTimes(gx0+jx,gz2,gI+jx,M-jx);
		addTimes(gx0+jx,gz3,gI+jx,M-jx);
		addTimes(gx0+jx,gz4,gI+jx,M-jx);
		addTimes(gx0+jx,gz5,gI+jx,M-jx);

           times(gx1+jy,gz0,gI+jy,M-jy);
		addTimes(gx1+jy,gz1,gI+jy,M-jy);
		addTimes(gx1+jy,gz2,gI+jy,M-jy);
		addTimes(gx1+jy,gz3,gI+jy,M-jy);
		addTimes(gx1+jy,gz4,gI+jy,M-jy);
		addTimes(gx1+jy,gz5,gI+jy,M-jy);

           times(gx2+jz,gz0,gI+jz,M-jz);
		addTimes(gx2+jz,gz1,gI+jz,M-jz);
		addTimes(gx2+jz,gz2,gI+jz,M-jz);
		addTimes(gx2+jz,gz3,gI+jz,M-jz);
		addTimes(gx2+jz,gz4,gI+jz,M-jz);
		addTimes(gx2+jz,gz5,gI+jz,M-jz);

        times(gx3,gz0+jz,gI,M-jz);
		addTimes(gx3,gz1+jz,gI,M-jz);
		addTimes(gx3,gz2+jz,gI,M-jz);
		addTimes(gx3,gz3+jz,gI,M-jz);
		addTimes(gx3,gz4+jz,gI,M-jz);
		addTimes(gx3,gz5+jz,gI,M-jz);

           times(gx4,gz0+jy,gI,M-jy);
		addTimes(gx4,gz2+jy,gI,M-jy);
		addTimes(gx4,gz2+jy,gI,M-jy);
		addTimes(gx4,gz3+jy,gI,M-jy);
		addTimes(gx4,gz4+jy,gI,M-jy);
		addTimes(gx4,gz5+jy,gI,M-jy);

        times(gx5,gz0+jx,gI,M-jx);
		addTimes(gx5,gz1+jx,gI,M-jx);
		addTimes(gx5,gz2+jx,gI,M-jx);
		addTimes(gx5,gz3+jx,gI,M-jx);
		addTimes(gx5,gz4+jx,gI,M-jx);
		addTimes(gx5,gz5+jx,gI,M-jx);


		cp(gs0,gx0,M);cp(gs1,gx1,M);cp(gs2,gx2,M);cp(gs3,gx3,M);cp(gs4,gx4,M);cp(gs5,gx5,M);
		gz0-=Mprev; gz1-=Mprev; gz2-=Mprev; gz3-=Mprev; gz4-=Mprev; gz5-=Mprev;
		gs0-=Mprev; gs1-=Mprev; gs2-=Mprev; gs3-=Mprev; gs4-=Mprev; gs5-=Mprev;
		gI-=Mprev;
	}
	//double sums_1=0, sums=0;
	//for (int i = 0; i<6*M; i++) { sums_1+=Gi[i+1][s-1]; sums +=Gi[i+1][s];} cout<< s << " FgiGs " << sums_1 << "  " << sums << endl;
}

void
Lat3D2ndO::PropagateF(Vector Gi, Vector G, const double S) const {
	int i,jx,jy, M,Mprev=0;
	M=GetTotalNumLayers();
	double L=(1.0-S)/4.0;
	double *gz0 = &Gi[1];
	double *gz1 = &Gi[1+M];
	double *gz2 = &Gi[1+2*M];
	double *gz3 = &Gi[1+3*M];
	double *gz4 = &Gi[1+4*M];
	double *gz5 = &Gi[1+5*M];


	SetBx1(gz0,gz5);SetBx1(gz1);SetBx1(gz2);SetBx1(gz3);SetBx1(gz4);
	SetBy1(gz1,gz4);SetBy1(gz0);SetBy1(gz2);SetBy1(gz3);SetBy1(gz5);
	SetBz1(gz2,gz3);SetBz1(gz0);SetBz1(gz1);SetBz1(gz4);SetBz1(gz5);
	SetBxm(gz0,gz5);SetBxm(gz1);SetBxm(gz2);SetBxm(gz3);SetBxm(gz4);
	SetBym(gz1,gz4);SetBym(gz0);SetBym(gz2);SetBym(gz3);SetBym(gz5);
	SetBzm(gz2,gz3);SetBzm(gz0);SetBzm(gz1);SetBzm(gz4);SetBzm(gz5);

	double *gx0 = &Gx0[1];
	double *gx1 = &Gx1[1];
	double *gx2 = &Gx2[1];
	double *gx3 = &Gx3[1];
	double *gx4 = &Gx4[1];
	double *gx5 = &Gx5[1];

	double *g = &G[1];
	double *gL = &GL[1]; timesC(gL,g,L,M);
	double *gS = &GS[1]; timesC(gS,g,S,M);

	for (i=1; i<= N_comp_ranges; i++){
		M=CR_Q[i]->Get_M();
		if (i>1) Mprev +=CR_Q[i-1]->Get_M();
		jx=CR_Q[i]->Get_jx();
		jy=CR_Q[i]->Get_jy();
		gz0+=Mprev; gz1+=Mprev; gz2+=Mprev; gz3+=Mprev; gz4+=Mprev; gz5+=Mprev;
		gS+=Mprev; gL+=Mprev;

		times(gx0+jx,gz0,gS+jx,M-jx);
		addTimes(gx0+jx,gz1,gL+jx,M-jx);
		addTimes(gx0+jx,gz2,gL+jx,M-jx);
		addTimes(gx0+jx,gz3,gL+jx,M-jx);
		addTimes(gx0+jx,gz4,gL+jx,M-jx);

        times(gx1+jy,gz0,gL+jy,M-jy);
		addTimes(gx1+jy,gz1,gS+jy,M-jy);
		addTimes(gx1+jy,gz2,gL+jy,M-jy);
		addTimes(gx1+jy,gz3,gL+jy,M-jy);
		addTimes(gx1+jy,gz5,gL+jy,M-jy);

        times(gx2+jz,gz0,gL+jz,M-jz);
		addTimes(gx2+jz,gz1,gL+jz,M-jz);
		addTimes(gx2+jz,gz2,gS+jz,M-jz);
		addTimes(gx2+jz,gz4,gL+jz,M-jz);
		addTimes(gx2+jz,gz5,gL+jz,M-jz);

        times(gx3,gz0+jz,gL,M-jz);
		addTimes(gx3,gz1+jz,gL,M-jz);
		addTimes(gx3,gz3+jz,gS,M-jz);
		addTimes(gx3,gz4+jz,gL,M-jz);
		addTimes(gx3,gz5+jz,gL,M-jz);

        times(gx4,gz0+jy,gL,M-jy);
		addTimes(gx4,gz2+jy,gL,M-jy);
		addTimes(gx4,gz3+jy,gL,M-jy);
		addTimes(gx4,gz4+jy,gS,M-jy);
		addTimes(gx4,gz5+jy,gL,M-jy);

        times(gx5,gz1+jx,gL,M-jx);
		addTimes(gx5,gz2+jx,gL,M-jx);
		addTimes(gx5,gz3+jx,gL,M-jx);
		addTimes(gx5,gz4+jx,gL,M-jx);
		addTimes(gx5,gz5+jx,gS,M-jx);

		cp(gz0,gx0,M);cp(gz1,gx1,M);cp(gz2,gx2,M);cp(gz3,gx3,M);cp(gz4,gx4,M);cp(gz5,gx5,M);

		gz0-=Mprev; gz1-=Mprev; gz2-=Mprev; gz3-=Mprev; gz4-=Mprev; gz5-=Mprev;
		gS-=Mprev;  gL-=Mprev;
	}
}

void
Lat3D2ndO::PropagateF(Vector Gi, Vector G, const double S, const bool stiffness) const {
	int i,jx,jy,M,Mprev=0;
	M=GetTotalNumLayers();
	double *pl=&Pl[1];
	double *ps=&Ps[1];
	for (int i=0; i<M; i++) {
		if (stiff->InRange(i+1)) {ps[i]=S;} else {ps[i]=0.2;}  pl[i]=(1.0-ps[i])/4.0;
	}
	//cout << "with stiff range S = " << S << endl;
	double *gz0 = &Gi[1];
	double *gz1 = &Gi[1+M];
	double *gz2 = &Gi[1+2*M];
	double *gz3 = &Gi[1+3*M];
	double *gz4 = &Gi[1+4*M];
	double *gz5 = &Gi[1+5*M];

	SetBx1(gz0,gz5);SetBx1(gz1);SetBx1(gz2);SetBx1(gz3);SetBx1(gz4);
	SetBy1(gz1,gz4);SetBy1(gz0);SetBy1(gz2);SetBy1(gz3);SetBy1(gz5);
	SetBz1(gz2,gz3);SetBz1(gz0);SetBz1(gz1);SetBz1(gz4);SetBz1(gz5);
	SetBxm(gz0,gz5);SetBxm(gz1);SetBxm(gz2);SetBxm(gz3);SetBxm(gz4);
	SetBym(gz1,gz4);SetBym(gz0);SetBym(gz2);SetBym(gz3);SetBym(gz5);
	SetBzm(gz2,gz3);SetBzm(gz0);SetBzm(gz1);SetBzm(gz4);SetBzm(gz5);

	double *gx0 = &Gx0[1];
	double *gx1 = &Gx1[1];
	double *gx2 = &Gx2[1];
	double *gx3 = &Gx3[1];
	double *gx4 = &Gx4[1];
	double *gx5 = &Gx5[1];

	double *g = &G[1];
	for (i=1; i<= N_comp_ranges; i++){
		M=CR_Q[i]->Get_M();
		if (i>1) Mprev +=CR_Q[i-1]->Get_M();
		jx=CR_Q[i]->Get_jx();
		jy=CR_Q[i]->Get_jy();
		gz0+=Mprev; gz1+=Mprev; gz2+=Mprev; gz3+=Mprev; gz4+=Mprev; gz5+=Mprev;
		g+=Mprev; ps+=Mprev; pl+=Mprev;

		times(gx0+jx,gz0,g+jx,ps,M-jx);
		addTimes(gx0+jx,gz1,g+jx,pl,M-jx);
		addTimes(gx0+jx,gz2,g+jx,pl,M-jx);
		addTimes(gx0+jx,gz3,g+jx,pl,M-jx);
		addTimes(gx0+jx,gz4,g+jx,pl,M-jx);

        times(gx1+jy,gz0,g+jy,pl,M-jy);
		addTimes(gx1+jy,gz1,g+jy,ps,M-jy);
		addTimes(gx1+jy,gz2,g+jy,pl,M-jy);
		addTimes(gx1+jy,gz3,g+jy,pl,M-jy);
		addTimes(gx1+jy,gz5,g+jy,pl,M-jy);

        times(gx2+jz,gz0,g+jz,pl,M-jz);
		addTimes(gx2+jz,gz1,g+jz,pl,M-jz);
		addTimes(gx2+jz,gz2,g+jz,ps,M-jz);
		addTimes(gx2+jz,gz4,g+jz,pl,M-jz);
		addTimes(gx2+jz,gz5,g+jz,pl,M-jz);

        times(gx3,gz0+jz,g,pl,M-jz);
		addTimes(gx3,gz1+jz,g,pl,M-jz);
		addTimes(gx3,gz3+jz,g,ps,M-jz);
		addTimes(gx3,gz4+jz,g,pl,M-jz);
		addTimes(gx3,gz5+jz,g,pl,M-jz);

        times(gx4,gz0+jy,g,pl,M-jy);
		addTimes(gx4,gz2+jy,g,pl,M-jy);
		addTimes(gx4,gz3+jy,g,pl,M-jy);
		addTimes(gx4,gz4+jy,g,ps,M-jy);
		addTimes(gx4,gz5+jy,g,pl,M-jy);

        times(gx5,gz1+jx,g,pl,M-jx);
		addTimes(gx5,gz2+jx,g,pl,M-jx);
		addTimes(gx5,gz3+jx,g,pl,M-jx);
		addTimes(gx5,gz4+jx,g,pl,M-jx);
		addTimes(gx5,gz5+jx,g,ps,M-jx);

		cp(gz0,gx0,M);cp(gz1,gx1,M);cp(gz2,gx2,M);cp(gz3,gx3,M);cp(gz4,gx4,M);cp(gz5,gx5,M);

		gz0-=Mprev; gz1-=Mprev; gz2-=Mprev; gz3-=Mprev; gz4-=Mprev; gz5-=Mprev;
		g-=Mprev; ps-=Mprev; pl-=Mprev;
	}
}

void
Lat3D2ndO::PropagateF(Vector Gi, Vector G) const {
	int i,jx,jy,M,Mprev=0;
	M=GetTotalNumLayers();

	double *gz0 = &Gi[1];
	double *gz1 = &Gi[1+M];
	double *gz2 = &Gi[1+2*M];
	double *gz3 = &Gi[1+3*M];
	double *gz4 = &Gi[1+4*M];
	double *gz5 = &Gi[1+5*M];

	SetBx1(gz0,gz5);SetBx1(gz1);SetBx1(gz2);SetBx1(gz3);SetBx1(gz4);
	SetBy1(gz1,gz4);SetBy1(gz0);SetBy1(gz2);SetBy1(gz3);SetBy1(gz5);
	SetBz1(gz2,gz3);SetBz1(gz0);SetBz1(gz1);SetBz1(gz4);SetBz1(gz5);
	SetBxm(gz0,gz5);SetBxm(gz1);SetBxm(gz2);SetBxm(gz3);SetBxm(gz4);
	SetBym(gz1,gz4);SetBym(gz0);SetBym(gz2);SetBym(gz3);SetBym(gz5);
	SetBzm(gz2,gz3);SetBzm(gz0);SetBzm(gz1);SetBzm(gz4);SetBzm(gz5);

	double *gx0 = &Gx0[1];
	double *gx1 = &Gx1[1];
	double *gx2 = &Gx2[1];
	double *gx3 = &Gx3[1];
	double *gx4 = &Gx4[1];
	double *gx5 = &Gx5[1];

	double *g = &G[1];
	double *gI = &GL[1]; timesC(gI,g,1.0/6.0,M); //to be called as the last forward propagator
	for (i=1; i<= N_comp_ranges; i++){
		M=CR_Q[i]->Get_M();
		if (i>1) Mprev +=CR_Q[i-1]->Get_M();
		jx=CR_Q[i]->Get_jx();
		jy=CR_Q[i]->Get_jy();
		gz0+=Mprev; gz1+=Mprev; gz2+=Mprev; gz3+=Mprev; gz4+=Mprev; gz5+=Mprev;
		gI+=Mprev;

		times(gx0+jx,gz0,gI+jx,M-jx);
		addTimes(gx0+jx,gz1,gI+jx,M-jx);
		addTimes(gx0+jx,gz2,gI+jx,M-jx);
		addTimes(gx0+jx,gz3,gI+jx,M-jx);
		addTimes(gx0+jx,gz4,gI+jx,M-jx);
		addTimes(gx0+jx,gz5,gI+jx,M-jx);

        times(gx1+jy,gz0,gI+jy,M-jy);
		addTimes(gx1+jy,gz1,gI+jy,M-jy);
		addTimes(gx1+jy,gz2,gI+jy,M-jy);
		addTimes(gx1+jy,gz3,gI+jy,M-jy);
		addTimes(gx1+jy,gz4,gI+jy,M-jy);
		addTimes(gx1+jy,gz5,gI+jy,M-jy);

        times(gx2+jz,gz0,gI+jz,M-jz);
		addTimes(gx2+jz,gz1,gI+jz,M-jz);
		addTimes(gx2+jz,gz2,gI+jz,M-jz);
		addTimes(gx2+jz,gz3,gI+jz,M-jz);
		addTimes(gx2+jz,gz4,gI+jz,M-jz);
		addTimes(gx2+jz,gz5,gI+jz,M-jz);

        times(gx3,gz0+jz,gI,M-jz);
		addTimes(gx3,gz1+jz,gI,M-jz);
		addTimes(gx3,gz2+jz,gI,M-jz);
		addTimes(gx3,gz3+jz,gI,M-jz);
		addTimes(gx3,gz4+jz,gI,M-jz);
		addTimes(gx3,gz5+jz,gI,M-jz);

        times(gx4,gz0+jy,gI,M-jy);
		addTimes(gx4,gz1+jy,gI,M-jy);
		addTimes(gx4,gz2+jy,gI,M-jy);
		addTimes(gx4,gz3+jy,gI,M-jy);
		addTimes(gx4,gz4+jy,gI,M-jy);
		addTimes(gx4,gz5+jy,gI,M-jy);

        times(gx5,gz0+jx,gI,M-jx);
		addTimes(gx5,gz1+jx,gI,M-jx);
		addTimes(gx5,gz2+jx,gI,M-jx);
		addTimes(gx5,gz3+jx,gI,M-jx);
		addTimes(gx5,gz4+jx,gI,M-jx);
		addTimes(gx5,gz5+jx,gI,M-jx);

		cp(gz0,gx0,M);cp(gz1,gx1,M);cp(gz2,gx2,M);cp(gz3,gx3,M);cp(gz4,gx4,M);cp(gz5,gx5,M);


		gz0-=Mprev; gz1-=Mprev; gz2-=Mprev; gz3-=Mprev; gz4-=Mprev; gz5-=Mprev;
		gI-=Mprev;
	}
}
void
Lat3D2ndO::PropagateB(Vector Gi, Vector G,const double S) const {
	int i,jx,jy, M,Mprev=0;
	M=GetTotalNumLayers();

	double L=(1.0-S)/4.0;
	double *gz0 = &Gi[1];
	double *gz1 = &Gi[1+M];
	double *gz2 = &Gi[1+2*M];
	double *gz3 = &Gi[1+3*M];
	double *gz4 = &Gi[1+4*M];
	double *gz5 = &Gi[1+5*M];

	SetBx1(gz0,gz5);SetBx1(gz1);SetBx1(gz2);SetBx1(gz3);SetBx1(gz4);
	SetBy1(gz1,gz4);SetBy1(gz0);SetBy1(gz2);SetBy1(gz3);SetBy1(gz5);
	SetBz1(gz2,gz3);SetBz1(gz0);SetBz1(gz1);SetBz1(gz4);SetBz1(gz5);
	SetBxm(gz0,gz5);SetBxm(gz1);SetBxm(gz2);SetBxm(gz3);SetBxm(gz4);
	SetBym(gz1,gz4);SetBym(gz0);SetBym(gz2);SetBym(gz3);SetBym(gz5);
	SetBzm(gz2,gz3);SetBzm(gz0);SetBzm(gz1);SetBzm(gz4);SetBzm(gz5);

	double *gx0 = &Gx0[1];
	double *gx1 = &Gx1[1];
	double *gx2 = &Gx2[1];
	double *gx3 = &Gx3[1];
	double *gx4 = &Gx4[1];
	double *gx5 = &Gx5[1];

	double *g = &G[1];
	double *gL = &GL[1]; timesC(gL,g,L,M);
	double *gS = &GS[1]; timesC(gS,g,S,M);
	for (i=1; i<= N_comp_ranges; i++){
		M=CR_Q[i]->Get_M();
		if (i>1) Mprev +=CR_Q[i-1]->Get_M();
		jx=CR_Q[i]->Get_jx();
		jy=CR_Q[i]->Get_jy();
		gz0+=Mprev; gz1+=Mprev; gz2+=Mprev; gz3+=Mprev; gz4+=Mprev; gz5+=Mprev;
		gS+=Mprev; gL+=Mprev;

        times(gx1+jx,gz5,gL+jx,M-jx);
        times(gx2+jx,gz5,gL+jx,M-jx);
        times(gx3+jx,gz5,gL+jx,M-jx);
		times(gx4+jx,gz5,gL+jx,M-jx);
		times(gx5+jx,gz5,gS+jx,M-jx);

        times(gx0+jy,gz4,gL+jy,M-jy);
		addTimes(gx2+jy,gz4,gL+jy,M-jy);
		addTimes(gx3+jy,gz4,gL+jy,M-jy);
		addTimes(gx4+jy,gz4,gS+jy,M-jy);
		addTimes(gx5+jy,gz4,gL+jy,M-jy);

        addTimes(gx0+jz,gz3,gL+jz,M-jz);
		addTimes(gx1+jz,gz3,gL+jz,M-jz);
		addTimes(gx3+jz,gz3,gS+jz,M-jz);
		addTimes(gx4+jz,gz3,gL+jz,M-jz);
		addTimes(gx5+jz,gz3,gL+jz,M-jz);

        addTimes(gx0,gz2+jz,gL,M-jz);
		addTimes(gx1,gz2+jz,gL,M-jz);
		addTimes(gx2,gz2+jz,gS,M-jz);
		addTimes(gx4,gz2+jz,gL,M-jz);
		addTimes(gx5,gz2+jz,gL,M-jz);

        addTimes(gx0,gz1+jy,gL,M-jy);
		addTimes(gx1,gz1+jy,gS,M-jy);
		addTimes(gx2,gz1+jy,gL,M-jy);
		addTimes(gx3,gz1+jy,gL,M-jy);
		addTimes(gx5,gz1+jy,gL,M-jy);

        addTimes(gx0,gz0+jx,gS,M-jx);
		addTimes(gx1,gz0+jx,gL,M-jx);
		addTimes(gx2,gz0+jx,gL,M-jx);
		addTimes(gx3,gz0+jx,gL,M-jx);
		addTimes(gx4,gz0+jx,gL,M-jx);


        cp(gz0,gx0,M);cp(gz1,gx1,M);cp(gz2,gx2,M);cp(gz3,gx3,M);cp(gz4,gx4,M);cp(gz5,gx5,M);


		gz0-=Mprev; gz1-=Mprev; gz2-=Mprev; gz3-=Mprev; gz4-=Mprev; gz5-=Mprev;
		gS-=Mprev; gL-=Mprev;
	}

}
void
Lat3D2ndO::PropagateB(Vector Gi, Vector G,const double S, const bool stiffness) const {
	int i,jx,jy, M,Mprev=0;
	M=GetTotalNumLayers();

	double *pl=&Pl[1];
	double *ps=&Ps[1];
	for (int i=0; i<M; i++) {
		if (stiff->InRange(i+1)) {ps[i]=S;} else {ps[i]=0.2;}  pl[i]=(1.0-ps[i])/4.0;
	}
	double *gz0 = &Gi[1];
	double *gz1 = &Gi[1+M];
	double *gz2 = &Gi[1+2*M];
	double *gz3 = &Gi[1+3*M];
	double *gz4 = &Gi[1+4*M];
	double *gz5 = &Gi[1+5*M];

	SetBx1(gz0,gz5);SetBx1(gz1);SetBx1(gz2);SetBx1(gz3);SetBx1(gz4);
	SetBy1(gz1,gz4);SetBy1(gz0);SetBy1(gz2);SetBy1(gz3);SetBy1(gz5);
	SetBz1(gz2,gz3);SetBz1(gz0);SetBz1(gz1);SetBz1(gz4);SetBz1(gz5);
	SetBxm(gz0,gz5);SetBxm(gz1);SetBxm(gz2);SetBxm(gz3);SetBxm(gz4);
	SetBym(gz1,gz4);SetBym(gz0);SetBym(gz2);SetBym(gz3);SetBym(gz5);
	SetBzm(gz2,gz3);SetBzm(gz0);SetBzm(gz1);SetBzm(gz4);SetBzm(gz5);

	double *gx0 = &Gx0[1];
	double *gx1 = &Gx1[1];
	double *gx2 = &Gx2[1];
	double *gx3 = &Gx3[1];
	double *gx4 = &Gx4[1];
	double *gx5 = &Gx5[1];

	double *g = &G[1]; //SetBoundaries(g);

	for (i=1; i<= N_comp_ranges; i++){
		M=CR_Q[i]->Get_M();
		if (i>1) Mprev +=CR_Q[i-1]->Get_M();
		jx=CR_Q[i]->Get_jx();
		jy=CR_Q[i]->Get_jy();
		gz0+=Mprev; gz1+=Mprev; gz2+=Mprev; gz3+=Mprev; gz4+=Mprev; gz5+=Mprev;
		g+=Mprev; ps+=Mprev; pl+=Mprev;

        times(gx1+jx,gz5,g+jx,pl,M-jx);
        times(gx2+jx,gz5,g+jx,pl,M-jx);
		times(gx3+jx,gz5,g+jx,pl,M-jx);
		times(gx4+jx,gz5,g+jx,pl,M-jx);
		times(gx5+jx,gz5,g+jx,ps,M-jx);

        times(gx0+jy,gz4,g+jy,pl,M-jy);
		addTimes(gx2+jy,gz4,g+jy,pl,M-jy);
		addTimes(gx3+jy,gz4,g+jy,pl,M-jy);
		addTimes(gx4+jy,gz4,g+jy,ps,M-jy);
		addTimes(gx5+jy,gz4,g+jy,pl,M-jy);

        addTimes(gx0+jz,gz3,g+jz,pl,M-jz);
		addTimes(gx1+jz,gz3,g+jz,pl,M-jz);
		addTimes(gx3+jz,gz3,g+jz,ps,M-jz);
		addTimes(gx4+jz,gz3,g+jz,pl,M-jz);
		addTimes(gx5+jz,gz3,g+jz,pl,M-jz);

        addTimes(gx0,gz2+jz,g,pl,M-jz);
		addTimes(gx1,gz2+jz,g,pl,M-jz);
		addTimes(gx2,gz2+jz,g,ps,M-jz);
		addTimes(gx4,gz2+jz,g,pl,M-jz);
		addTimes(gx5,gz2+jz,g,pl,M-jz);

        addTimes(gx0,gz1+jy,g,pl,M-jy);
		addTimes(gx1,gz1+jy,g,ps,M-jy);
		addTimes(gx2,gz1+jy,g,pl,M-jy);
		addTimes(gx3,gz1+jy,g,pl,M-jy);
		addTimes(gx5,gz1+jy,g,pl,M-jy);

        addTimes(gx0,gz0+jx,g,ps,M-jx);
		addTimes(gx1,gz0+jx,g,pl,M-jx);
		addTimes(gx2,gz0+jx,g,pl,M-jx);
		addTimes(gx3,gz0+jx,g,pl,M-jx);
		addTimes(gx4,gz0+jx,g,pl,M-jx);


        cp(gz0,gx0,M);cp(gz1,gx1,M);cp(gz2,gx2,M);cp(gz3,gx3,M);cp(gz4,gx4,M);cp(gz5,gx5,M);


		gz0-=Mprev; gz1-=Mprev; gz2-=Mprev; gz3-=Mprev; gz4-=Mprev; gz5-=Mprev;
		g-=Mprev; ps-=Mprev; pl-=Mprev;
	}
}
void
Lat3D2ndO::PropagateB(Vector Gi, Vector G) const {
	int i,jx,jy, M,Mprev=0;
	M=GetTotalNumLayers();

	double *gz0 = &Gi[1];
	double *gz1 = &Gi[1+M];
	double *gz2 = &Gi[1+2*M];
	double *gz3 = &Gi[1+3*M];
	double *gz4 = &Gi[1+4*M];
	double *gz5 = &Gi[1+5*M];

	SetBx1(gz0,gz5);SetBx1(gz1);SetBx1(gz2);SetBx1(gz3);SetBx1(gz4);
	SetBy1(gz1,gz4);SetBy1(gz0);SetBy1(gz2);SetBy1(gz3);SetBy1(gz5);
	SetBz1(gz2,gz3);SetBz1(gz0);SetBz1(gz1);SetBz1(gz4);SetBz1(gz5);
	SetBxm(gz0,gz5);SetBxm(gz1);SetBxm(gz2);SetBxm(gz3);SetBxm(gz4);
	SetBym(gz1,gz4);SetBym(gz0);SetBym(gz2);SetBym(gz3);SetBym(gz5);
	SetBzm(gz2,gz3);SetBzm(gz0);SetBzm(gz1);SetBzm(gz4);SetBzm(gz5);

	double *gx0   = &Gx0[1];
	double *gx1   = &Gx1[1];
	double *gx2   = &Gx2[1];
	double *gx3   = &Gx3[1];
	double *gx4   = &Gx4[1];
	double *gx5   = &Gx5[1];

	double *g = &G[1];
	double *gI = &GL[1]; timesC(gI,g,1.0/6.0,M); //to be called as the first backward step.

	for (i=1; i<= N_comp_ranges; i++){
		M=CR_Q[i]->Get_M();
		if (i>1) Mprev +=CR_Q[i-1]->Get_M();
		jx=CR_Q[i]->Get_jx();
		jy=CR_Q[i]->Get_jy();
		gz0+=Mprev; gz1+=Mprev; gz2+=Mprev; gz3+=Mprev; gz4+=Mprev; gz5+=Mprev;
		gI+=Mprev;

		times(gx0+jx,gz5,gI+jx,M-jx);
		times(gx1+jx,gz5,gI+jx,M-jx);
		times(gx2+jx,gz5,gI+jx,M-jx);
		times(gx3+jx,gz5,gI+jx,M-jx);
		times(gx4+jx,gz5,gI+jx,M-jx);
		times(gx5+jx,gz5,gI+jx,M-jx);

		addTimes(gx0+jy,gz4,gI+jy,M-jy);
		addTimes(gx1+jy,gz4,gI+jy,M-jy);
		addTimes(gx2+jy,gz4,gI+jy,M-jy);
		addTimes(gx3+jy,gz4,gI+jy,M-jy);
		addTimes(gx4+jy,gz4,gI+jy,M-jy);
		addTimes(gx5+jy,gz4,gI+jy,M-jy);

        addTimes(gx0+jz,gz3,gI+jz,M-jz);
		addTimes(gx1+jz,gz3,gI+jz,M-jz);
		addTimes(gx2+jz,gz3,gI+jz,M-jz);
		addTimes(gx3+jz,gz3,gI+jz,M-jz);
		addTimes(gx4+jz,gz3,gI+jz,M-jz);
		addTimes(gx5+jz,gz3,gI+jz,M-jz);

        addTimes(gx0,gz2+jz,gI,M-jz);
		addTimes(gx1,gz2+jz,gI,M-jz);
		addTimes(gx2,gz2+jz,gI,M-jz);
		addTimes(gx3,gz2+jz,gI,M-jz);
		addTimes(gx4,gz2+jz,gI,M-jz);
		addTimes(gx5,gz2+jz,gI,M-jz);

              addTimes(gx0,gz1+jy,gI,M-jy);
		addTimes(gx1,gz1+jy,gI,M-jy);
		addTimes(gx2,gz1+jy,gI,M-jy);
		addTimes(gx3,gz1+jy,gI,M-jy);
		addTimes(gx4,gz1+jy,gI,M-jy);
		addTimes(gx5,gz1+jy,gI,M-jy);

              addTimes(gx0,gz0+jx,gI,M-jx);
		addTimes(gx1,gz0+jx,gI,M-jx);
		addTimes(gx2,gz0+jx,gI,M-jx);
		addTimes(gx3,gz0+jx,gI,M-jx);
		addTimes(gx4,gz0+jx,gI,M-jx);
		addTimes(gx5,gz0+jx,gI,M-jx);


        cp(gz0,gx0,M);cp(gz1,gx1,M);cp(gz2,gx2,M);cp(gz3,gx3,M);cp(gz4,gx4,M);cp(gz5,gx5,M);


		gz0-=Mprev; gz1-=Mprev; gz2-=Mprev; gz3-=Mprev; gz4-=Mprev; gz5-=Mprev;
		gI-=Mprev;
	}
}


void
Lat3D2ndO::Init2G(Vector Gi1, Vector Gi2, const Vector G, const LatticeRange* LatRange) const {
	Message(fatal,"Init2G not implemented in 3d");
}
void
Lat3D2ndO::Propagate2G(Matrix Gi1, Matrix Gi2, const Vector G, const int s,const LatticeRange* LatRange) const {
Message(fatal,"Propagate 2G not implemented in 3d");
}
void
Lat3D2ndO::Propagate2G(Matrix Gi1, Matrix Gi2, const Vector G, const int s,const LatticeRange* LatRange,const double f) const {
Message(fatal,"Propagate 2G not implemented in 3d");
}

void
Lat3D2ndO::Propagate2G(Vector Gi1, Vector Gi2, const Vector G, const LatticeRange* LatRange) const {
Message(fatal,"Propagate 2G not implemented in 3d");
}
void
Lat3D2ndO::Propagate2G(Vector Gi1, Vector Gi2, const Vector G, const LatticeRange* LatRange,const double f) const {
Message(fatal,"Propagate 2G not implemented in 3d");
}

void
Lat3D2ndO::PropagateG(Matrix Gi1, const Vector G1, const int s) const {}

void
Lat3D2ndO::PropagateG(Vector G1, const Vector G2) const {}

void
Lat3D2ndO::PropagateG(const Vector Gi, const Vector G, Vector Gout) const {}

void
Lat3D2ndO::PropagateG(Matrix Gi1, const Vector G1, const int s,const double f) const {}

void
Lat3D2ndO::PropagateG(Vector G1, const Vector G2,const double f) const {}

void
Lat3D2ndO::PropagateG(const Vector Gi, const Vector G, Vector Gout, const double f) const {}


void
Lat3D2ndO::ConnectG(const Vector GiA, const Matrix GiB, const int s, Vector out) const {
	double *p_out = &out[1];
	for (int k=0; k<6; k++) {
		const double *p_giB=&GiB[1+k*GetTotalNumLayers()][s];
		const double *p_giA = &GiA[1+k*GetTotalNumLayers()];
		addTimes(p_out,p_giA,p_giB,GetTotalNumLayers());
	}
}

void
Lat3D2ndO::ConnectG(const Vector GiA, const Matrix GiB, const int s, Vector out, double *OUT) const {
	int kk;
	int M=GetTotalNumLayers();
	double *p_out = &out[1];
	for (int k=0; k<6; k++) {
		if (k==0 || k==5) kk=0;
		if (k==1 || k==4) kk=1;
		if (k==2 || k==3) kk=2;
		const double *p_giB=&GiB[1+k*M][s];
		const double *p_giA = &GiA[1+k*M];
		addTimes(p_out,p_giA,p_giB,M);
		addTimes(OUT+kk*M,p_giA,p_giB,M);
	}
}


//void
//Lat3D2ndO::ConnectG(const Vector GiA, const Matrix GiB, const int s, Vector phi, Vector phi_p) const {
//	double *p_phi = &phi[1];
//	double *p_phi_p = &phi_p[1];
//	for (int k=0; k<6; k++) {
//		const double *p_giB=&GiB[1+k*GetTotalNumLayers()][s];
//		const double *p_giA = &GiA[1+k*GetTotalNumLayers()];
//		addTimes(p_phi,p_giA,p_giB,GetTotalNumLayers());
//		if (k==2 || k==3) addTimes(p_phi_p,p_giA,p_giB,GetTotalNumLayers());
//	}
//}
void
Lat3D2ndO::ConnectG(const Vector GiA, const Vector GiB, Vector out) const {
	double *p_out = &out[1];
	for (int k=0; k<6; k++) {
		const double *p_giA = &GiA[1+k*GetTotalNumLayers()];
		const double *p_giB = &GiB[1+k*GetTotalNumLayers()];
		addTimes(p_out,p_giA,p_giB,GetTotalNumLayers());
	}
}
void
Lat3D2ndO::ConnectG(const Vector GiA, const Vector GiB, Vector out, double *OUT) const {
	int kk;
	int M=GetTotalNumLayers();
	double *p_out = &out[1];
	for (int k=0; k<6; k++) {
		if (k==0 || k==5) kk=0;
		if (k==1 || k==4) kk=1;
		if (k==2 || k==3) kk=2;
		const double *p_giA = &GiA[1+k*M];
		const double *p_giB = &GiB[1+k*M];
		addTimes(p_out,p_giA,p_giB,M);
		addTimes(OUT+kk*M,p_giA,p_giB,M);
	}
}


//void
//Lat3D2ndO::ConnectG(const Vector GiA, const Vector GiB, Vector phi, Vector phi_p) const {
//	double *p_phi = &phi[1];
//	double *p_phi_p = &phi_p[1];
//	for (int k=0; k<6; k++) {
//		const double *p_giA = &GiA[1+k*GetTotalNumLayers()];
//		const double *p_giB = &GiB[1+k*GetTotalNumLayers()];
//		addTimes(p_phi,p_giA,p_giB,GetTotalNumLayers());
//		if (k==2 || k==3) addTimes(p_phi_p,p_giA,p_giB,GetTotalNumLayers());
//	}
//}

void
Lat3D2ndO::Connect2G(const Vector GiA1, const Matrix GiB1, const int s1, const Vector GiA2,
						 const Matrix GiB2, const int s2, Vector out) const {
	Message(fatal,"Connect 2G not implemented in 3d");
}

void
Lat3D2ndO::Connect2G(const Vector GiA1, const Vector GiB1, const Vector GiA2,  const Vector GiB2, Vector out) const {
	Message(fatal,"Connect 2G not implemented in 3d");
}
void
Lat3D2ndO::CorrectDoubleCountG(Vector in, const Vector G) const {
	double *p_in = &in[1];
	const double *p_G = &G[1];
    //Message(debug,"CorrectDoubleCountG(Vector in, const Vector G)");
	div(p_in,p_G,GetTotalNumLayers());
}
void
Lat3D2ndO::CorrectDoubleCountG(Vector in, double *IN, const Vector G) const {
	int M=GetTotalNumLayers();
	double *p_in = &in[1];
	const double *p_G = &G[1];
    //Message(debug,"CorrectDoubleCountG(Vector in, const Vector G)");
	div(p_in,p_G,M); for (int k=0; k<3; k++) div(IN+k*M,p_G,M);
}
double
Lat3D2ndO::ComputeLnGN(Vector Gi) const {
	double value = 0;
	int i,x,y,z,px,py,NlayersX,NlayersY,NlayersZ;
	//int Previous=0;
	int Current=0,Next=0,jx,jy;
	for (int k=0; k<6; k++) {
		double *pGi = &Gi[1+k*GetTotalNumLayers()];
		//Message(debug,"ComputeLnGN(Vector Gi)");
		SubtractBoundaries(pGi);
		//Previous=0;
		Current=0;Next=0;
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
	}
	return log(value/6.0);
}
void
Lat3D2ndO::NormPhiFree(Vector phi, const double C) const {
	norm(phi,C,GetTotalNumLayers());
}
void
Lat3D2ndO::NormPhiFree(Vector phi, double *PHI, const double C) const {
	int  M=GetTotalNumLayers();
	norm(phi,C,M); for (int k=0; k<3; k++) norm(PHI+k*M,C,M);
}
void
Lat3D2ndO::NormPhiRestr(Vector phi, const Vector Gi, double C) const {
	C /= 6.0*exp(ComputeLnGN(Gi));
	norm(phi,C,GetTotalNumLayers());
}

void
Lat3D2ndO::NormPhiRestr(Vector phi,double *PHI, const Vector Gi, double C) const {
	int M=GetTotalNumLayers();
	C /= 6.0*exp(ComputeLnGN(Gi));
	norm(phi,C,M); for (int k=0; k<3; k++) norm(PHI+k*M,C,M);
}


void
Lat3D2ndO::UpdateBoundaries(Vector A) const {
	//Message(debug,"UpdateBoundaries(Vector A)");
	SetBoundaries(A);
}
void
Lat3D2ndO::UpdateBoundaries(Matrix A, const int s) const {
	for (int k=0; k<6; k++) {
		double *pA=&A[1+k*GetTotalNumLayers()][s];
		SetBoundaries(pA);
	}
}
