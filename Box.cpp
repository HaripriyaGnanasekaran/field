#include "Box.h"
#include <iostream>
#include <string>
#include <fstream>
#include <algorithm>
#include <f2c.h>
#include <clapack.h> 
using namespace std;


void bx(double *P, int mmx, int My, int Mz, int bx1, int bxm, int jx, int jy){
	int i;
	int jx_mmx=jx*mmx;
	int jx_bxm=jx*bxm;
	int bx1_jx=bx1*jx;
	for (int y=0; y<My; y++)
	for (int z=0; z<Mz; z++){
		i=jy*y+z;
		P[i]=P[bx1_jx+i];
		P[jx_mmx+i]=P[jx_bxm+i];
	}
}

void b_x(double *P, int mmx, int My, int Mz, int bx1, int bxm, int jx, int jy){
	int i, jx_mmx=jx*mmx;// jx_bxm=jx*bxm, bx1_jx=bx1*jx;
	for (int y=0; y<My; y++)
	for (int z=0; z<Mz; z++){
		i=jy*y+z;
		P[i]=0;
		P[jx_mmx+i]=0;
	}
}
void by(double *P, int Mx, int mmy, int Mz, int by1, int bym, int jx, int jy){
	int i, jy_mmy=jy*mmy, jy_bym=jy*bym, jy_by1=jy*by1;
	for (int x=0; x<Mx; x++)
	for (int z=0; z<Mz; z++) {
		i=jx*x+z;
		P[i]=P[jy_by1+i];
		P[jy_mmy+i]=P[jy_bym+i];
	}
}
void b_y(double *P, int Mx, int mmy, int Mz, int by1, int bym, int jx, int jy){
	int i, jy_mmy=jy*mmy;// jy_bym=jy*bym, jy_by1=jy*by1;
	for (int x=0; x<Mx; x++)
	for (int z=0; z<Mz; z++) {
		i=jx*x+z;
		P[i]=0;
		P[jy_mmy+i]=0;
	}
}
void bz(double *P, int Mx, int My, int mmz, int bz1, int bzm, int jx, int jy){
	int i;
	for (int x=0; x<Mx; x++)
	for (int y=0; y<My; y++) {
		i=jx*x+jy*y;
		P[i]=P[i+bz1];
		P[i+mmz]=P[i+bzm];
	}
}
void b_z(double *P, int Mx, int My, int mmz, int bz1, int bzm, int jx, int jy){
	int i;
	for (int x=0; x<Mx; x++)
	for (int y=0; y<My; y++) {
		i=jx*x+jy*y;
		P[i]=0;
		P[i+mmz]=0;
	}
}

double Dot(double *x,double *y,int M){
	double result=0;
 	for (int i=0; i<M; i++) result +=x[i]*y[i];
	return result;
}

double Sum(double *x,int M){
	double result=0;
 	for (int i=0; i<M; i++) result +=x[i];
	return result;
}

void AddTimes(double *P, double *A, double *B, int M){
	for (int i=0; i<M; i++) P[i]+=A[i]*B[i];
}
void Times(double *P, double *A, double *B, int M){
	for (int i=0; i<M; i++) P[i]=A[i]*B[i];
}
void Norm(double *P, double C, int M){
	for (int i=0; i<M; i++) P[i] *= C;

}
void Zero(double* P, int M){
	for (int i=0; i<M; i++) P[i] =0;
}
void Cp(double *P,double *A, int M){
	for (int i=0; i<M; i++) P[i] = A[i];
}
void YisAminB(double *Y, double *A, double *B, int M){
	for (int i=0; i<M; i++) Y[i] = A[i]-B[i];
}
void YplusisCtimesX(double *Y, double *X, double C, int M) {
	for (int i=0; i<M; i++) Y[i] += C*X[i];
}
void Add(double *P, double *A, int M){
	for (int i=0; i<M; i++) P[i]+=A[i];
}
void Dubble(double *P, double *A, double norm){
	for (int i=0; i<M; i++) P[i]*=norm/A[i];

}
void Boltzmann(double *P, double *A){
	for (int i=0; i<M; i++) P[i]=exp(-A[i]);
}
void PutAlpha(double *g, double *phitot, double *phi_side, double chi, double phibulk){
	for (int i=0; i<M; i++) g[i] = g[i] - chi*(phi_side[i]/phitot[i]-phibulk);
}
void AddG(double *g, double *phitot, double *alpha){
	for (int i=0; i<M; i++) g[i]= g[i] -alpha[i] +1/phitot[i]-1;
}

void SetBoundaries(double *P) {
	bx(P,Mx+1,My+2,Mz+2,bx1,bxm,jx,jy);
	by(P,Mx+2,My+1,Mz+2,by1,bym,jx,jy);
	bz(P,Mx+2,My+2,Mz+1,bz1,bzm,jx,jy);
}
void RemoveBoundaries(double *g) {
	b_x(g,Mx+1,My+2,Mz+2,bx1,bxm,jx,jy);
	b_y(g,Mx+2,My+1,Mz+2,by1,bym,jx,jy);
	b_z(g,Mx+2,My+2,Mz+1,bz1,bzm,jx,jy);
}
void Side(double *X_side, double *X) {
	Zero(X_side,M); SetBoundaries(X);
	Add(X_side+jx,X,M-jx); Add(X_side,X+jx,M-jx);
	Add(X_side+jy,X,M-jy); Add(X_side,X+jy,M-jy);
	Add(X_side+1,X,M-1);  Add(X_side,X+1, M-1);
	Norm(X_side,1.0/6.0,M);
}
void Propagate(double *G, double *G1, int s_from, int s_to) {
	double *gs = G+M*s_to, *gs_1 = G+M*s_from, *g = G1;
	Side(gs,gs_1);
	Times(gs,gs,g,M);
}
void ComputePhi(){
	Boltzmann(G1,u); SetBoundaries(G1);
	Cp(phi,G1,M);
	Cp(Gg_f,mask1,M); 
	for (s=1; s<=N; s++) Propagate(Gg_f,G1,s-1,s); Propagate(Gg_f,mask2,N,N+1);
	Cp(Gg_b+((N+1)%2)*M,mask2,M);
	Zero(phi+M,M); 
	for (int s=N; s>0; s--) {
		Propagate(Gg_b,G1,((s+1)%2),(s%2));
		AddTimes(phi+M,Gg_f+s*M,Gg_b+(s%2)*M,M);
	}
	GN=Sum(Gg_f+(N+1)*M,M);
	Dubble(phi+M,G1,1/GN);
}

void ComputeG(){
	ComputePhi(); 
	Cp(phitot,phi,M); Add(phitot,phi+M,M);
	for (int i=0; i<n_seg; i++) Side(phi_side+i*M,phi+i*M);
	Cp(g,u,iv);
	PutAlpha(g,phitot,phi_side+M,CHI,1.0);
	PutAlpha(g+M,phitot,phi_side,CHI,0.0);

	Cp(alpha,g,M); Add(alpha,g+M,M);  Norm(alpha,1.0/n_seg,M); //alpha contains average alpha
	for (int i=0; i<n_seg; i++) {AddG(g+i*M,phitot,alpha); RemoveBoundaries(g+i*M);}
}

void Ax(double* A, double* X, int N){//From Ax_B; below B is not used: it is assumed to contain a row of unities.
	double *U = new double[N*N];
	double *S = new double[N];
	double *VT = new double[N*N];
	integer MM = (integer)N, NN = (integer)N;
	integer LDA=MM, LDU=MM, LDVT=NN, INFO, LWORK;
	int lwork;
	double WKOPT;
	double* WORK;
	char JOBU='S'; //'S' is nodig om alleen de eerste N colommen in U te schrijven.
	char JOBVT='A';

	LWORK = -1; //grootte hulpgeheugen aanvragen
	dgesvd_( &JOBU, &JOBVT, &MM, &NN, A, &LDA, S, U, &LDU, VT, &LDVT, &WKOPT, &LWORK, &INFO );
	lwork = (int)WKOPT;
	WORK = (double*)malloc( lwork*sizeof(double) );
	LWORK = (integer)lwork; //nu uitrekenen.
	dgesvd_( &JOBU, &JOBVT, &MM, &NN, A, &LDA, S, U, &LDU, VT, &LDVT,WORK, &LWORK, &INFO );
	if (INFO >0) { //error message genereren
	};
	free(WORK);
	for (int i=0; i<N; i++) X[i]=0;
	for (int i=0; i<N; i++) for (int j=0; j<N; j++) X[i] += U[i*N + j];//*B[j];
	for (int i=0; i<N; i++) {S[i] = X[i]/S[i]; X[i]=0;} //S is use decause it is no longer needed.
	for (int i=0; i<N; i++) for (int j=0; j<N; j++) X[i] += VT[i*N + j]*S[j];
	delete U;
	delete S;
	delete VT;
}
void DIIS(double *x, double *x_x0, double *xR, double *Aij, double *Apij,double *Ci, int k, int m, int iv) {
	double normC=0; int posi;
	if (k_diis>m) { k_diis =m;
		for (int i=1; i<m; i++) for (int j=1; j<m; j++)
		Aij[m*(i-1)+j-1]=Aij[m*i+j]; //remove oldest elements
	}
	for (int i=0; i<k_diis; i++) {posi = k-k_diis+1+i; if (posi<0) posi +=m;
		Aij[i+m*(k_diis-1)] = Aij[k_diis-1+m*i] = Dot(x_x0+posi*iv, x_x0+k*iv,iv);	}
		// write to (compressed) matrix Apij
	for (int i=0; i<k_diis; i++) for (int j=0; j<k_diis; j++) {
		Apij[j+k_diis*i] = Aij[j+m*i];
	}
	Ax(Apij,Ci,k_diis);
	for (int i=0; i<k_diis; i++) normC +=Ci[i];
	for (int i=0; i<k_diis; i++) {Ci[i] =Ci[i]/normC; }
	Zero(x,iv);
	posi = k-k_diis+1; if (posi<0) posi +=m;

	YplusisCtimesX(x,xR+posi*iv,Ci[0],iv); //pv = Ci[0]*xR[0];
	for (int i=1; i<k_diis; i++) {
		posi = k-k_diis+1+i; if (posi<0) posi +=m;
		YplusisCtimesX(x,xR+posi*iv,Ci[i],iv);
	}
}

double Helmholtz(){
	double F_Helmholtz=0, Uint=0, entropy=0;
	RemoveBoundaries(phi); RemoveBoundaries(phi+M);
	Uint = CHI*Dot(phi_side,phi+M,M);
	entropy =  Dot(phi+M,u+M,M)-n*log(n*N/GN);
	entropy += Dot(phi,u,M);
	F_Helmholtz = Uint - entropy;
	return F_Helmholtz;
}

double SCF() {
	Zero(x,iv); Zero(x0,iv);
	it=0; k_diis=1; k=0; 
	ComputeG();
	YplusisCtimesX(x,g,-eta,iv);
	YisAminB(x_x0,x,x0,iv);
	Cp(xR,x,iv);
	error = sqrt(Dot(g,g,iv));
	printf("DIIS has been notified\n");
	printf("Your guess = %1e \n",error);
	while (error > tolerance && it < 1000) {
		it++;
		Cp(x0,x,iv); ComputeG();
		k=it % m; k_diis++; //plek voor laatste opslag
		YplusisCtimesX(x,g,-eta,iv);
		Cp(xR+k*iv,x,iv); YisAminB(x_x0+k*iv,x,x0,iv);
		DIIS(x,x_x0,xR,Aij,Apij,Ci,k,m,iv);
		error = sqrt(Dot(g,g,iv));
		printf("it = %i error = %1e \n",it,error);
	}
	return Helmholtz();
};


int main(int argc, char *argv[]) {
	string fname;
	string filename;
	string line;

	if (argc == 2) fname = argv[1]; else {printf("Use: box filename \n"); return 1;}
	filename = fname + ".in";

	ifstream in_file;
	string calculation_type;
	string SCF_type;

	in_file.open(filename.c_str());
		getline(in_file,calculation_type,':'); 
		getline(in_file,SCF_type); 

    		SCF_type.erase(std::remove(SCF_type.begin(), SCF_type.end(), ' '), SCF_type.end());
		if (SCF_type == "FULL") cout << "SCF calculations" << endl; else cout << "PHI calculations" << endl;
	in_file.close(); 
	
//	filename = fname + ".pot";
//	in_file.open(filename.c_str());
//		in_file >> Mx >> My >> Mz; 
//	        jx = (My+2)*(Mz+2); jy = Mz+2; M=(Mx+2)*(My+2)*(Mz+2);
//		for (int i=0; i<Mx; i++) 
//		for (int j=0; j<My; j++)
//		for (int k=0; k<Mz; k++) 
//		in_file >> u[jx*i+jy*j+k];
//
//	in_file.close();

	Mx=My=Mz=50; 	//box size
	N=100; 	//chain length
	n_box =1;	//number of boxex

	n_seg=2; //number of segment types.
 	n=n_box; //number of chains

//	printf("Tolerance (1e-7): "); scanf("%lf", &tolerance);
// 	printf("Regularisation parameter, eta (0.1): "); scanf("%lf", &eta);
//	printf("Memory depth (m) : "); scanf("%d", &m);

	m=10;

	jx = (My+2)*(Mz+2); jy = Mz+2; M=(Mx+2)*(My+2)*(Mz+2);
	//bx1=Mx; bxm=1; by1=My; bym=1; bz1=Mz; bzm=1; // periodic
	bx1=1; bxm=Mx; by1=1; bym=My; bz1=1; bzm=Mz; //reflecting

	iv=n_seg*M;

	Aij = new double[m*m]; for (int i=0; i<m*m; i++) Aij[i]=0;
	Ci = new double[m]; for (int i=0; i<m; i++) Ci[i]=0;
	Apij = new double[m*m]; for (int i=0; i<m*m; i++) Apij[i]=0;

	phi = new double[iv];
	phitot = new double[M]; G1 = new double[iv]; alpha = new double[M];
	Gg_f = new double[M*(N+2)]; Gg_b = new double[M*2]; phi_side = new double[iv]; phi = new double [iv]; x = new double[iv];
	x0 = new double[iv]; g = new double[iv]; xR = new double[m*iv]; x_x0 =new double[m*iv];
	mask1 = new double[M]; 
	mask2 = new double[M];
	px1=new int[n_box]; py1=new int[n_box]; pz1=new int[n_box]; px2=new int[n_box]; py2=new int[n_box]; pz2=new int[n_box];

	u = x;

	Zero(mask1,M); Zero(mask2,M);
	for (int p=0; p<n_box; p++) {
		px1[p]=5;
		py1[p]=5;
		pz1[p]=5;
		mask1[px1[p]*jx+py1[p]*jy+pz1[p]]=1;
		px2[p]=10;
		py2[p]=10;
		pz2[p]=10;
		mask2[px2[p]*jx+py2[p]*jy+pz2[p]]=1;
	}

	double Free_energy = SCF();
	cout << " phitot = " << Sum(phi+M,M) <<  endl ;
 

	//printf("z \t phi  \n"); for (int zz=1; zz<=Mz; zz++) printf(" %i \t %1f \n", zz, phi[jx+jy+zz]);
	printf("Free energy : %1f \n", Free_energy);

	//fprintf (vtkfile,"# vtk DataFile Version 3.0\nvtk output\nASCII\nDATASET STRUCTURED_POINTS\nDIMENSIONS %i %i %i\n",Mx,My,Mz);
	//fprintf (vtkfile,"SPACING 1 1 1\nORIGIN 0 0 0\nPOINT_DATA %i\n",Mx*My*Mz);
	//fprintf (vtkfile,"SCALARS SFBox_profile double\nLOOKUP_TABLE default\n");
	//for (int i=1;i<=Mx;i++){
	//	for (int j=1;j<=My;j++){
	//		for (int l=1;l<=Mz;l++){
	//			fprintf (vtkfile,"%e\n",phi[i*jx+j*jy+l]);
	//		}
	//	}
	//}
	//fclose (vtkfile);

	return 0;
};
