#include <math.h>
#include <stdio.h> 
#include <cuda.h>
#include <cublas_v2.h> 
#include <cuda_runtime.h> 
#include <f2c.h>  
#include <clapack.h> 
#include <cstdlib>
#include <ctime>

//nvcc LL.cu -lm -lcuda -lcudart -llapack -lblas -lf2c -lcublas -arch=sm_20 -o LL
//I.V. Ionova, E.A. Carter, "Error vector choice in direct inversion in the iterative subspace method, J. Compt. Chem. 17, 1836-1847, 1996. 

double *BlasResult;
cublasStatus_t stat; 
cublasHandle_t handle;
int n_mol, n_seg;
int *N;
double *GN, *phibulk, *theta, *n;
double chi; 
int block_size=256,i,j,k,kk,k_diis,kk_diis,m,mm,s,M,Mz,it,b1,bm,iv,iiv;
double sigma, error = 1, tolerance = 1e-7, eta = 0.1, normC, pi=4.0*atan(1.0);
double *Aij,*Ci,*Apij,*H_phi,*H_mask,*H_u,*mask,*phi,*phi_pinned,*phitot,*G1,*alpha,*Gg_f,*Gg_b,*phi_side,*x,*x0,*g,*xR,*x_x0;
double *u,*f,*r,*r0;
double *AAij,*CCi,*AApij,*r_r0;

__global__ void times(double *P, double *A, double *B, int M){int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) P[idx]=A[idx]*B[idx];
}
__global__ void addtimes(double *P, double *A, double *B, int M){int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) P[idx]+=A[idx]*B[idx];
}
__global__ void norm(double *P, double C, int M){int idx = blockIdx.x*blockDim.x+threadIdx.x;
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
__global__ void addg(double *g, double *phitot, double *alpha, int M) {int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) {g[idx]= g[idx] -alpha[idx] +1/phitot[idx]-1;
	}
}
__global__ void bz(double *P, int mmz, int b1, int bm ){
	P[0]=P[b1];
	P[mmz]=P[bm];
}
__global__ void b_z(double *P, int Mx, int My, int mmz, int bz1, int bzm, int jx, int jy){
	P[0]=0; 
	P[mmz]=0;
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
	bz<<<n_blocks,block_size>>>(P,Mz+1,b1,bm);
}
void RemoveBoundaries(double *g) {
	b_z<<<n_blocks,block_size>>>(g,Mz+1,b1,bm);
}
void Side(double *phi_side, double *phi) {
	Zero(phi_side,M); SetBoundaries(phi);
       yplusisctimesx(phi_side,phi,4.0,M);
	Add(phi_side+1, phi,M-1);  Add(phi_side,phi+1, M-1);
	Norm(phi_side,1.0/6.0,M);
}
void Propagate(double *G, double *G1, int s_from, int s_to) {
	double *gs = G+M*(s_to-1), *gs_1 = G+M*(s_from-1), *g = G1;
	Side(gs,gs_1);
	times(gs,gs,g);
};
void ComputePhi(){
	int N1=N[1];
	int Ndiv2=N1/2;
	for (int i=0; i<n_seg; i++) { Boltzmann(G1+i*M,u+i*M); SetBoundaries(G1+i*M); }
	Cp(phi,G1,M);//phi for solvent;	
	Cp(Gg_f,G1+M,M); for (s=2; s<=Ndiv2; s++) Propagate(Gg_f,G1+M,s-1,s);
	Cp(Gg_b+(Ndiv2%2)*M,Gf_f+(Ndiv2-1)*M,M); 
	Zero(phi+M,M);
	for (int s=Ndiv2+1; s<=N1; s++) {
		Propagate(Gg_b,G1+M,((s-1)%2)+1,(s%2)+1); 
		AddTimes(phi+M,Gg_f+(N1-s)*M,Gg_b+(s%2)*M,M);
	}
	RemoveBoundaries(Gg_b+(N1%2)*M);
	GN[1]=Sum(Gg_b+(N1%2)*M,M); 
 	Dubble(phi+M,G1+M,2*n[1]/GN[1]); 
	       	
	phibulk[1]=N1*n[1]/GN[1];
	phibulk[0]=1-phibulk[1];
	Norm(phi,phibulk[0]);
}
void AddG(double *g, double *phitot, double *alpha){
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	addg<<<n_blocks,block_size>>>(g,phitot,alpha,M);
}
void ComputeG(){
	ComputePhi();
	Cp(phitot,phi,M); for (int i=1; i<n_seg; i++) Add(phitot,phi+i*M,M);
	for (int i=0; i<n_seg; i++) Side(phi_side+i*M,phi+i*M); 	}
	Cp(g,u,iv); 
	for (int i=0; i<n_seg; i++) for (int j=0; j<n_seg; j++) {
		if ((i != j) && (chi[i+n_seg*j] !=0) ) {
			PutAlpha(g+i*M,phitot,phi_side+j*M,chi[i+n_seg*j],phibulk[j]);
		}
	} //now g contains segment type dependent alpha.
	Cp(alpha,g,M); for (int i=1; i<n_seg; i++) Add(alpha,g+i*M,M);  Norm(alpha,1.0/n_seg,M); //alpha contains average alpha
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
void DIIS(double *x, double *x_x0, double *xR, double *Aij, double *Apij, double *Ci, int k, int m, int iv) {
	double normC=0; int posi;  
	if (k_diis>m) { k_diis =m;
		for (int i=1; i<m; i++) for (int j=1; j<m; j++) 
		Aij[m*(i-1)+j-1]=Aij[m*i+j]; //remove oldest elements
	}
	for (int i=0; i<k_diis; i++) {
		posi = k-k_diis+1+i; if (posi<0) posi +=m; 
		Aij[i+m*(k_diis-1)] = Aij[k_diis-1+m*i] = Dot(x_x0+posi*iv, x_x0+k*iv,iv);	
	}
	for (int i=0; i<k_diis; i++) for (int j=0; j<k_diis; j++) {
		Apij[j+k_diis*i] = Aij[j+m*i];
	}
	Ax(Apij,Ci,k_diis);		
	for (int i=0; i<k_diis; i++) normC +=Ci[i];
	for (int i=0; i<k_diis; i++) Ci[i] =Ci[i]/normC; 
	Zero(x,iv);
	posi = k-k_diis+1; if (posi<0) posi +=m; 
		
	YplusisCtimesX(x,xR+posi*iv,Ci[0],iv); //pv = Ci[0]*xR[0];
	for (int i=1; i<k_diis; i++) { 
		posi = k-k_diis+1+i; if (posi<0) posi +=m; 
		YplusisCtimesX(x,xR+posi*iv,Ci[i],iv); 
	}
}
double Helmholtz(){
	double F_Helmholtz=0;
	F_Helmholtz = -n[1]*log(GN[1]/n[1]/N[1]);
	RemoveBoundaries(alpha);
	F_Helmholtz -= Sum(alpha,M);
	return F_Helmholtz;
}
double SCF() {
	TransferDataToDevice(M, H_mask, mask); 
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

	stat = cublasCreate(&handle); if (stat !=CUBLAS_STATUS_SUCCESS) {printf("CUBLAS failed \n");}
	cublasSetPointerMode(handle, CUBLAS_POINTER_MODE_DEVICE);
	printf("n_layers z: "); scanf("%d", &Mz);	
n_mol=2;	
	n = new double[n_mol]; n[0]=0; 
	phibulk = new double[n_mol]; 
	GN = new double[n_mol];
	N = new int[n_mol];
n_seg=2;
	printf("mol 1 : solvent S \n"); N[0]=1;
	printf("mol 2 : (A)N;  Give N:"); scanf("%d", &N[1]); 
	printf("give n mol 2: "); scanf("%lf", &n[1]);
	
	printf("chi "); scanf("%lf", &chi); 
//	printf("Tolerance (1e-7): "); scanf("%lf", &tolerance);
// 	printf("Regularisation parameter, eta (0.1): "); scanf("%lf", &eta);
//	printf("Memory depth (m) : "); scanf("%d", &m);
m=10;
	M=Mz+2;
	iv=n_seg*M;
	b1=1; bm=Mz;

	Aij = new double[m*m]; for (int i=0; i<m*m; i++) Aij[i]=0;
	Ci = new double[m]; for (int i=0; i<m; i++) Ci[i]=0;
	Apij = new double[m*m]; for (int i=0; i<m*m; i++) Apij[i]=0;

	H_phi = new double[M]; 	H_mask = new double[M]; H_u = new double[M]; mask = (double*)AllOnDev(M); phi = (double*)AllOnDev(iv);
	phi_pinned = (double*)AllOnDev(M); phitot = (double*)AllOnDev(M); G1 = (double*)AllOnDev(iv); alpha = (double*)AllOnDev(M);
	Gg_f = (double*)AllOnDev(M*N[1]); Gg_b = (double*)AllOnDev(M*2); phi_side = (double*)AllOnDev(iv); x = (double*)AllOnDev(iv);
	x0 = (double*)AllOnDev(iv); g = (double*)AllOnDev(iv); xR = (double*)AllOnDev(m*iv); x_x0 = (double*)AllOnDev(m*iv);
	BlasResult = (double*)AllOnDev(1);

	u = x; 

	double Free_energy = SCF();
	

	
	//TransferDataToHost(H_phi,phi);
	//printf("z \t phi  \n"); for (int zz=1; zz<=Mz; zz++) printf(" %i \t %1f \n", zz, H_phi[jx+jy+zz]);
	//TransferDataToHost(H_phi,phi+1*M);
	//printf("z \t phi  \n"); for (int zz=1; zz<=Mz; zz++) printf(" %i \t %1f \n", zz, H_phi[jx+jy+zz]);
	//TransferDataToHost(H_phi,phi+2*M);
	//printf("z \t phi  \n"); for (int zz=1; zz<=Mz; zz++) printf(" %i \t %1f \n", zz, H_phi[jx+jy+zz]);

	free(H_phi);
	cudaFree(phi); cudaFree(x);
	cudaFree(G1);cudaFree(alpha);cudaFree(Gg_f);cudaFree(Gg_b);cudaFree(phi_side);
	cudaFree(x_x0);
	cublasDestroy(handle);
	return(0);
};
