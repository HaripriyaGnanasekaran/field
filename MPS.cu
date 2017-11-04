#include <math.h>
#include <stdio.h> 
#include <cuda.h>
#include <cublas_v2.h> 
#include <cuda_runtime.h> 
#include <f2c.h>  
#include <clapack.h> 
#include <cstdlib> 
#include <ctime>
 
//nvcc Open.cu -lm -lcuda -lcudart -llapack -lblas -lf2c -lcublas -arch=sm_20 -o open
//I.V. Ionova, E.A. Carter, "Error vector choice in direct inversion in the iterative subspace method, J. Compt. Chem. 17, 1836-1847, 1996. 

double *BlasResult;
cublasStatus_t stat; 
cublasHandle_t handle;
int N,N1,N2,n_seg;
double CHI,n,GN;
int block_size=256,i,j,k,k_diis,m,s,Mx,My,Mz,M,jx,jy,it,bx1,by1,bz1,bxm,bym,bzm,iv;
double  error = 1, tolerance = 1e-7, eta = 0.1;
double *Aij,*Ci,*Apij,*H_phi,*H_u,*phi,*phi_pinned,*phitot,*G1,*alpha,*Gg_f,*Gg_b,*phi_side,*x,*x0,*g,*xR,*x_x0;
double *u; //alleen een pointer, voor ons gemak.

__global__ void times(double *P, double *A, double *B, int M){int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) P[idx]=A[idx]*B[idx];
}
__global__ void addtimes(double *P, double *A, double *B, int M){int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) P[idx]+=A[idx]*B[idx];
}
__global__ void norm(double *P, double C, int M){
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) P[idx] *= C;
}
__global__ void zero(double *P, int M){int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) P[idx] = 0.0;
}
__global__ void cp (double *P, double *A, int M){int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) P[idx] = A[idx];
}
__global__ void yisaminb(double *Y, double *A,double *B, int M){int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) Y[idx] = A[idx]-B[idx];
}
__global__ void yplusisctimesx(double *Y, double *X, double C, int M){int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) Y[idx] += C*X[idx];
}
__global__ void add(double *P, double *A, int M){int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) P[idx]+=A[idx];
}
__global__ void dubble(double *P, double *A, double norm, int M){
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) P[idx]*=norm/A[idx];
}
__global__ void boltzmann(double *P, double *A, int M){int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) P[idx]=exp(-A[idx]);
}
__global__ void putalpha(double *g,double *phitot,double *phi_side,double chi,double phibulk,int M){
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) g[idx] = g[idx] - chi*(phi_side[idx]/phitot[idx]-phibulk);
}
__global__ void addg(double *g, double *phitot, double *alpha, int M) { 
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) {
		g[idx]= g[idx] -alpha[idx] +1/phitot[idx]-1;
	}
}
__global__ void bx(double *P, int mmx, int My, int Mz, int bx1, int bxm, int jx, int jy){
	int idx, jx_mmx=jx*mmx, jx_bxm=jx*bxm, bx1_jx=bx1*jx;
	int yi =blockIdx.x*blockDim.x+threadIdx.x, zi =blockIdx.y*blockDim.y+threadIdx.y;
	if (yi<My && zi<Mz) {
		idx=jy*yi+zi;
		P[idx]=P[bx1_jx+idx];
		P[jx_mmx+idx]=P[jx_bxm+idx];
	}
}
__global__ void b_x(double *P, int mmx, int My, int Mz, int bx1, int bxm, int jx, int jy){
	int idx, jx_mmx=jx*mmx;
	int yi =blockIdx.x*blockDim.x+threadIdx.x, zi =blockIdx.y*blockDim.y+threadIdx.y;
	if (yi<My && zi<Mz) {
		idx=jy*yi+zi;
		P[idx]=0;
		P[jx_mmx+idx]=0;
	}
}
__global__ void by(double *P, int Mx, int mmy, int Mz, int by1, int bym, int jx, int jy){
	int idx, jy_mmy=jy*mmy, jy_bym=jy*bym, jy_by1=jy*by1;
	int xi =blockIdx.x*blockDim.x+threadIdx.x, zi =blockIdx.y*blockDim.y+threadIdx.y;
	if (xi<Mx && zi<Mz) {
		idx=jx*xi+zi;
		P[idx]=P[jy_by1+idx];
		P[jy_mmy+idx]=P[jy_bym+idx];
	}
}
__global__ void b_y(double *P, int Mx, int mmy, int Mz, int by1, int bym, int jx, int jy){
	int idx, jy_mmy=jy*mmy;
	int xi =blockIdx.x*blockDim.x+threadIdx.x, zi =blockIdx.y*blockDim.y+threadIdx.y;
	if (xi<Mx && zi<Mz) {
		idx=jx*xi+zi;
		P[idx]=0;
		P[jy_mmy+idx]=0;
	}
}
__global__ void bz(double *P, int Mx, int My, int mmz, int bz1, int bzm, int jx, int jy){
	int idx, xi =blockIdx.x*blockDim.x+threadIdx.x, yi =blockIdx.y*blockDim.y+threadIdx.y;
	if (xi<Mx && yi<My) {
		idx=jx*xi+jy*yi;
		P[idx]=P[idx+bz1];
		P[idx+mmz]=P[idx+bzm];
	}
}
__global__ void b_z(double *P, int Mx, int My, int mmz, int bz1, int bzm, int jx, int jy){
	int idx, xi =blockIdx.x*blockDim.x+threadIdx.x, yi =blockIdx.y*blockDim.y+threadIdx.y;
	if (xi<Mx && yi<My) {
		idx=jx*xi+jy*yi;
		P[idx]=0; 
		P[idx+mmz]=0;
	}
}
double Dot(double *x,double *y,int M){ 
	double result;
 	cublasDdot(handle,M,x,1,y,1,BlasResult);
	cudaMemcpy(&result,BlasResult,sizeof(double),cudaMemcpyDeviceToHost);
	return result;
}
double Sum(double *x,int M){
	double result;
 	cublasDasum(handle,M,x,1,BlasResult);
	cudaMemcpy(&result,BlasResult,sizeof(double),cudaMemcpyDeviceToHost);
	return result;
}
bool GPU_present() {
    int deviceCount =0; cuDeviceGetCount(&deviceCount);
    if (deviceCount ==0) printf("There is no device supporting Cuda.\n");
    else cudaDeviceReset();
	return deviceCount > 0;
}
double *AllOnDev(int N) {
	double *X;
	if (cudaSuccess != cudaMalloc((void **) &X, sizeof(double)*N))
	printf("Memory allocation on GPU failed.\n Please reduce size of lattice and/or chain length(s)\n");
	return X;
}
void TransferDataToHost(double *H, double *D) {
	cudaMemcpy(H, D, sizeof(double)*M,cudaMemcpyDeviceToHost);
}
void TransferDataToDevice(int M, double *H, double *D ) { 
	cudaMemcpy(D, H, sizeof(double)*M,cudaMemcpyHostToDevice);
}
void AddTimes(double *P, double *A, double *B, int M){int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	addtimes<<<n_blocks,block_size>>>(P,A,B,M);
}
void Times(double *P, double *A, double *B, int M){ int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	times<<<n_blocks,block_size>>>(P,A,B,M);
}
void Norm(double *P, double C, int M){
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	norm<<<n_blocks,block_size>>>(P,C,M);
}
void Zero(double* P, int M){int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	zero<<<n_blocks,block_size>>>(P,M);
}
void Cp(double *P,double *A, int M){int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	cp<<<n_blocks,block_size>>>(P,A,M);
}
void YisAminB(double *Y, double *A, double *B, int M){ int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	yisaminb<<<n_blocks,block_size>>>(Y,A,B,M);
}
void Add(double *P, double *A, int M){int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	add<<<n_blocks,block_size>>>(P,A,M);
}
void Dubble(double *P, double *A, double norm){ 
       int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	dubble<<<n_blocks,block_size>>>(P,A,norm,M);
}
void Boltzmann(double *P, double *A){int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	boltzmann<<<n_blocks,block_size>>>(P,A,M);
}
void PutAlpha(double *g, double *phitot, double *phi_side, double chi, double phibulk){
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	putalpha<<<n_blocks,block_size>>>(g,phitot,phi_side,chi,phibulk,M);
}
void SetBoundaries(double *P) {
	dim3 dimBlock(16,16);
	dim3 dimGridz((Mx+dimBlock.x+1)/dimBlock.x,(My+dimBlock.y+1)/dimBlock.y);
	dim3 dimGridy((Mx+dimBlock.x+1)/dimBlock.x,(Mz+dimBlock.y+1)/dimBlock.y);
	dim3 dimGridx((My+dimBlock.x+1)/dimBlock.x,(Mz+dimBlock.y+1)/dimBlock.y);
	bx<<<dimGridx,dimBlock>>>(P,Mx+1,My+2,Mz+2,bx1,bxm,jx,jy);
	by<<<dimGridy,dimBlock>>>(P,Mx+2,My+1,Mz+2,by1,bym,jx,jy);
	bz<<<dimGridz,dimBlock>>>(P,Mx+2,My+2,Mz+1,bz1,bzm,jx,jy);
}
void RemoveBoundaries(double *g) {
	dim3 dimBlock(16,16);
	dim3 dimGridz((Mx+dimBlock.x+1)/dimBlock.x,(My+dimBlock.y+1)/dimBlock.y);
	dim3 dimGridy((Mx+dimBlock.x+1)/dimBlock.x,(Mz+dimBlock.y+1)/dimBlock.y);
	dim3 dimGridx((My+dimBlock.x+1)/dimBlock.x,(Mz+dimBlock.y+1)/dimBlock.y);
	b_x<<<dimGridx,dimBlock>>>(g,Mx+1,My+2,Mz+2,bx1,bxm,jx,jy);
	b_y<<<dimGridy,dimBlock>>>(g,Mx+2,My+1,Mz+2,by1,bym,jx,jy);
	b_z<<<dimGridz,dimBlock>>>(g,Mx+2,My+2,Mz+1,bz1,bzm,jx,jy);
}
void Side(double *X_side, double *X) {
	Zero(X_side,M); SetBoundaries(X);
	Add(X_side+jx,X,M-jx); Add(X_side,X+jx,M-jx);
	Add(X_side+jy,X,M-jy); Add(X_side,X+jy,M-jy);
	Add(X_side+1,X,M-1);  Add(X_side,X+1, M-1);
	Norm(X_side,1.0/6.0,M);
}
void Propagate(double *G, double *G1, int s_from, int s_to) {
	double *gs = G+M*(s_to-1), *gs_1 = G+M*(s_from-1), *g = G1;
	Side(gs,gs_1); 
	Times(gs,gs,g,M);
};
void ComputePhi(){
	for (int i=0; i<n_seg; i++) { Boltzmann(G1+i*M,u+i*M); SetBoundaries(G1+i*M); }
	Cp(Gg_f,G1,M); 
	for (s=2; s<=N1; s++) Propagate(Gg_f,G1,s-1,s);
	for (s=N1+1; s<=N; s++) Propagate(Gg_f,G1+M,s-1,s);
 
	Cp(Gg_b+(N%2)*M,G1+M,M); 	
	Times(phi+M,Gg_f+(N-1)*M,Gg_b+(N%2)*M,M);	
	
	for (int s=N-1; s>N1; s--) {
		Propagate(Gg_b,G1+M,((s+1)%2)+1,(s%2)+1); 
		AddTimes(phi+M,Gg_f+(s-1)*M,Gg_b+(s%2)*M,M);
	}
	Zero(phi,M);
	for (int s=N1; s>=1; s--) {
		Propagate(Gg_b,G1,((s+1)%2)+1,(s%2)+1); 
		AddTimes(phi,Gg_f+(s-1)*M,Gg_b+(s%2)*M,M); 
	}
	RemoveBoundaries(Gg_f+(N-1)*M);
	GN=Sum(Gg_f+(N-1)*M,M); 
	for (int i=0; i<n_seg; i++) Dubble(phi+i*M,G1+i*M,n/GN); 
}
void AddG(float *g, float *phitot, float *alpha){
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	addg<<<n_blocks,block_size>>>(g,phitot,alpha,M);
}
void ComputeG(){
	ComputePhi();
	Cp(phitot,phi,M); Add(phitot,phi+M,M);
	for (int i=0; i<n_seg; i++) Side(phi_side+i*M,phi+i*M);

	Cp(g,u,iv); 
	PutAlpha(g,phitot,phi_side+M,CHI,N2/N);
	PutAlpha(g+M,phitot,phi_side,CHI,N1/N);

	Cp(alpha,g,M); Add(alpha,g+M,M);  Norm(alpha,1.0/n_seg,M); //alpha contains average alpha
	for (int i=0; i<n_seg; i++) {AddG(g+i*M,phitot,alpha); RemoveBoundaries(g+i*M);}
}
void YplusisCtimesX(double *Y, double *X, double C, int M) {
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	yplusisctimesx<<<n_blocks,block_size>>>(Y,X,C,M);
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
	dgesvd_( &JOBU, &JOBVT, &MM, &NN, A, &LDA, S, U, &LDU, VT, &LDVT, WORK, &LWORK, &INFO );
	if (INFO >0) { //error message genereren
	};
	delete WORK;
	for (int i=0; i<N; i++) X[i]=0;
	for (int i=0; i<N; i++) for (int j=0; j<N; j++) X[i] += U[i*N + j];//*B[j];
	for (int i=0; i<N; i++) {S[i] = X[i]/S[i]; X[i]=0;} //S is use decause it is no longer needed.
	for (int i=0; i<N; i++) for (int j=0; j<N; j++) X[i] += VT[i*N + j]*S[j];
	delete U,S,VT;
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
	entropy = Dot(phi+M,u+M,M)-n*log(n*N/GN);	
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
}

main() {  
	cudaDeviceReset();
	char fname[150];
	stat = cublasCreate(&handle); if (stat !=CUBLAS_STATUS_SUCCESS) {printf("CUBLAS failed \n");}
	cublasSetPointerMode(handle, CUBLAS_POINTER_MODE_DEVICE);
	printf("n_layers x (10): "); scanf("%d", &Mx);
	printf("n_layers y (10): "); scanf("%d", &My);
	printf("n_layers z (10): "); scanf("%d", &Mz);		

	printf("(A)N_1-(B)N_2: Give N_1:"); scanf("%d", &N1);
	printf("(A)N_1-(B)N_2: Give N_2:"); scanf("%d", &N2);
	N=N1+N2;
	printf("chi(A,B) "); scanf("%lf", &CHI); 
       printf("Enter file name for vtk-output: ");  scanf("%s", &fname);
	FILE *vtkfile = fopen(fname,"w");
	n_seg=2;
	n=Mx*My*Mz/N; //number of chains in system.
	
//	printf("Tolerance (1e-7): "); scanf("%lf", &tolerance);
// 	printf("Regularisation parameter, eta (0.1): "); scanf("%lf", &eta);
//	printf("Memory depth (m) : "); scanf("%d", &m);
m=10;
	jx = (My+2)*(Mz+2); jy = Mz+2; M=(Mx+2)*(My+2)*(Mz+2);
	bx1=Mx; bxm=1; by1=My; bym=1; bz1=Mz; bzm=1; // periodic
	iv=n_seg*M;

	Aij = new double[m*m]; for (int i=0; i<m*m; i++) Aij[i]=0;
	Ci = new double[m]; for (int i=0; i<m; i++) Ci[i]=0;
	Apij = new double[m*m]; for (int i=0; i<m*m; i++) Apij[i]=0;

	H_phi = new double[M];  H_u = new double[M]; phi = (double*)AllOnDev(iv);
	phi_pinned = (double*)AllOnDev(M); phitot = (double*)AllOnDev(M); G1 = (double*)AllOnDev(iv); alpha = (double*)AllOnDev(M);
	Gg_f = (double*)AllOnDev(M*N); Gg_b = (double*)AllOnDev(M*2); phi_side = (double*)AllOnDev(iv); x = (double*)AllOnDev(iv);
	x0 = (double*)AllOnDev(iv); g = (double*)AllOnDev(iv); xR = (double*)AllOnDev(m*iv); x_x0 = (double*)AllOnDev(m*iv);
	BlasResult = (double*)AllOnDev(1);

	u = x; 
	double Free_energy = SCF();
//-------------
	TransferDataToHost(H_phi,phi);
	printf("z \t phi  \n"); for (int zz=1; zz<=Mz; zz++) printf(" %i \t %1f \n", zz, H_phi[jx+jy+zz]);
	printf("Free energy : %1f \n", Free_energy); 

	fprintf (vtkfile,"# vtk DataFile Version 3.0\nvtk output\nASCII\nDATASET STRUCTURED_POINTS\nDIMENSIONS %i %i %i\n",Mx,My,Mz);
	fprintf (vtkfile,"SPACING 1 1 1\nORIGIN 0 0 0\nPOINT_DATA %i\n",Mx*My*Mz);
	fprintf (vtkfile,"SCALARS SFBox_profile double\nLOOKUP_TABLE default\n");
	for (int i=1;i<=Mx;i++){
		for (int j=1;j<=My;j++){
			for (int l=1;l<=Mz;l++){
				fprintf (vtkfile,"%e\n",H_phi[i*jx+j*jy+l]);
			}		
		}		
	}		
	fclose (vtkfile);

	return(0);
};
