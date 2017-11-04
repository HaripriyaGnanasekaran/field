#include <math.h>
#include <stdio.h> 
#include <stdlib.h> 
#include <cuda.h>
#include <cublas_v2.h> //versie 4?   
#include <cuda_runtime.h> 
#include <f2c.h>  
#include <clapack.h>
#include "mtrand.h"
//#include <cstdio> // may not be needed

//nvcc Openbox.cu mtrand.cpp -lm -lcuda -lcudart -llapack -lblas -lf2c -lcublas -arch=sm_20 -o open
//I.V. Ionova, E.A. Carter, "Error vector choice in direct inversion in the iterative subspace method, J. Compt. Chem. 17, 1836-1847, 1996. 

// Declare Global Variables

double *DotResult;
int *IntResult;
cublasStatus_t stat; 
cublasHandle_t handle;
const int block_size=256;
int i,j,l,k,m,k_diis,s,M,jx,jy,it,bx1,by1,bz1,bxm,bym,bzm,pol,sol,iv,Nchain,forwards=0,underflow,numlength,baselength,randstart=0;// The maximum value that can be set for block size512
bool solve;
double sigma, error = 1, normC, GN, theta, delta_u=0.1,sumamp=0;
double *Aij,*Ci,*Apij, *averageG; 
double *H_randvector,*H_vector,*H_mask,*mask,*phi,*phitot,*G1,*alpha,*Gg_f,*Gg_b,*phi_side,*x,*x0,*g,*xR,*x_x0,*temp;
double *u;
int *Md;
char digits[] = {'0', '1', '2', '3', '4', '5', '6', '7', '8', '9' };
int mcs =0; // the number of montecarlosteps taken
FILE * failedparticleconfigs;
#include "data.h" 

// initialise random number generator
MTRand_int32 irand(4); // initialize mersenne twister random number generater with value 4 it also declares irand which returns a random 32 bit integer 	
MTRand frand; // double in [0, 1) generator, already init downloaded from http://www.bedaux.net/mtrand/

// Cuda function Declarations

__global__ void powerinplace(double *P,double C, int M){int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) P[idx] = pow(P[idx],C);
	
}
__global__ void times(double *P, double *A, double *B, int M){int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) P[idx]=A[idx]*B[idx];
}
__global__ void timesinplace(double *P, double *A, int M){int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) P[idx] *= A[idx];
}
__global__ void addtimes(double *P, double *A, double *B, int M){int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) P[idx]+=A[idx]*B[idx];
}
__global__ void norm(double *P, double C, int M){
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) P[idx] *= C;
}
__global__ void addconstinplace(double *P, double C, int M){
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) P[idx] += C;
}
__global__ void zero(double *P, int M){int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) P[idx] = 0.0;
}

__global__ void add(double *P, double *A, int M){int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if(idx< M) P[idx] += A[idx];	
}
__global__ void yisaminb(double *Y, double *A,double *B, int M){int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) Y[idx] = A[idx]-B[idx];
}
__global__ void yplusisctimesx(double *Y, double *X, double C, int M){int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) Y[idx] += C*X[idx];
}
__global__ void cp(double *P, double *A, int M){int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) P[idx] =A[idx];
}
__global__ void dubble(double *P, double *A, double norm, int M){
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) P[idx]*=norm/A[idx];
}
__global__ void boltzmann(double *P, double *A, int M){int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) P[idx]=exp(-A[idx]);
}
__global__ void putalpha(double *g,double *u,double *phitot,double *phi_side,double chi,double phibulk,int M){
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) g[idx]=u[idx]-chi*(phi_side[idx]/phitot[idx]-phibulk);
}
__global__ void addg(double *g, double *phitot, double *alpha, int M) { 
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) {
		g[idx]=g[idx]-alpha[idx] +(1/phitot[idx]-phitot[idx]);
		g[idx+M]=g[idx+M]-alpha[idx] +(1/phitot[idx]-phitot[idx]);
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
__global__ void initialize (double *Gg_b, double *G1, double *mask, int N, int M){
	double *P=Gg_b+(N-1)*M;
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) P[idx] = G1[idx]*mask[idx];
}

__global__ void side(double *P,double *A,int M,int jx, int jy){
	//int idx, jx_mmx=jx*Mx+1, jx_bxm=jx*bxm, bx1_jx=bx1*jx;
	//int yi =blockIdx.x*blockDim.x+threadIdx.x, zi =blockIdx.y*blockDim.y+threadIdx.y;
	//if (yi<My+2 && zi<Mz+2) {
	//	idx=jy*yi+zi;
	//	A[idx]=A[bx1_jx+idx];
	//	A[jx_mmx+idx]=A[jx_bxm+idx];
	//}
	//int idx, jy_mmy=jy*My+1, jy_bym=jy*bym, jy_by1=jy*by1;
	//int xi =blockIdx.x*blockDim.x+threadIdx.x, zi =blockIdx.y*blockDim.y+threadIdx.y;
	//if (xi<Mx+2 && zi<Mz+2) {
	//	idx=jx*xi+zi;
	//	A[idx]=A[jy_by1+idx];
	//	A[jy_mmy+idx]=A[jy_bym+idx];
	//}
	//int idx, xi =blockIdx.x*blockDim.x+threadIdx.x, yi =blockIdx.y*blockDim.y+threadIdx.y;
	//if (xi<Mx+2 && yi<My+2) {
	//	idx=jx*xi+jy*yi;
	//	A[idx]=A[idx+bz1];
	//	A[idx+mmz]=A[idx+bzm];
	//}
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if(idx<M-jx && idx>=jx) {P[idx] = (A[idx-jx]+A[idx-jy]+A[idx-1]+A[idx+1]+A[idx+jy]+A[idx+jx])/6.0;}
	//else{ P[idx]=0;}
}

double Dot(double *x,double *y,int M){ 
	double result;
 	cublasDdot(handle,M,x,1,y,1,DotResult);
	cudaMemcpy(&result,DotResult,sizeof(double),cudaMemcpyDeviceToHost);
	return result;
}

double Sum(double *x,int M){
	double result;
 	cublasDasum(handle,M,x,1,DotResult);
	cudaMemcpy(&result,DotResult,sizeof(double),cudaMemcpyDeviceToHost);
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

// Function declarations

void TransferDataToHost(double *H, double *D) {
	cudaMemcpy(H, D, sizeof(double)*M,cudaMemcpyDeviceToHost);
}
void TransferDataToDevice(int M, double *H, double *D ) { 
	cudaMemcpy(D, H, sizeof(double)*M,cudaMemcpyHostToDevice);
}
void AddTimes(double *P, double *A, double *B, int M){int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	addtimes<<<n_blocks,block_size>>>(P,A,B,M);
}
void Powerinplace(double *P,double C,int M){ int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	powerinplace<<<n_blocks,block_size>>>(P,C,M);
}
void Times(double *P, double *A, double *B, int M){ int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	times<<<n_blocks,block_size>>>(P,A,B,M);
}
void Timesinplace(double *P, double *A, int M){ int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	timesinplace<<<n_blocks,block_size>>>(P,A,M);
}
void Norm(double *P, double C, int M){
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	norm<<<n_blocks,block_size>>>(P,C,M);
}
void Addconstinplace(double *P, double C, int M){
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	addconstinplace<<<n_blocks,block_size>>>(P,C,M);
}
void Zero(double* P, int M){//int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	//zero<<<n_blocks,block_size>>>(P,M);
	cudaMemset(P,0.0,M*sizeof(double));
}

void Cp(double *P,double *A, int M){int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	cp<<<n_blocks,block_size>>>(P,A,M);
}
void YisAminB(double *Y, double *A, double *B, int M){ int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	yisaminb<<<n_blocks,block_size>>>(Y,A,B,M);
}
void Add(double *P, double *A, int M){int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	//cudaEvent_t start,stop;
	//float time;
	//cudaEventCreate(&start);
	//cudaEventCreate(&stop);
	//cudaEventRecord(start,0);
	add<<<n_blocks,block_size>>>(P,A,M);
	//cudaEventRecord(stop,0);
	//cudaEventSynchronize(stop);
	//cudaEventElapsedTime(&time,start,stop);//givestimeinmilliseconds
	//cudaEventDestroy(start);
	//cudaEventDestroy(stop);
	//printf("time to do addition %e nblock %d blocksize %d \n",time,n_blocks,block_size);
}
void Dubble(double *P, double *A, double norm){ 
      	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	dubble<<<n_blocks,block_size>>>(P,A,norm,M);
}
void Boltzmann(double *P, double *A){int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	boltzmann<<<n_blocks,block_size>>>(P,A,M);
}
void PutAlpha(double *g, double *u, double *phitot, double *phi_side, double chi, double phibulk){
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	putalpha<<<n_blocks,block_size>>>(g,u,phitot,phi_side,chi,phibulk,M);
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
void Initialize(double *Gg_b, double *G1, double *mask){
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	initialize<<<n_blocks,block_size>>>(Gg_b, G1, mask, N_g, M);
}

//_global__ void side(double *P,double *A,int M,int Mx,int My,int Mz,int bx1,int by1,int bz1, int bxm,int bym,int bzm, int jx, int jy){



void Side(double *phi_side, double *phi) {int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	
	SetBoundaries(phi);
	side<<<n_blocks,block_size>>>(phi_side,phi,M,jx,jy);
	//Add(phi_side+jx,phi,M-jx); Add(phi_side,phi+jx,M-jx);
	//Add(phi_side+jy,phi,M-jy); Add(phi_side,phi+jy,M-jy);
	//Add(phi_side+1, phi,M-1);  Add(phi_side,phi+1, M-1);
	//Norm(phi_side,1.0/6.0,M);
}
void Propagate(double *G, double *G1, int s_from, int s_to) {
	double *gs = G+M*(s_to-1), *gs_1 = G+M*(s_from-1), *g = G1; 
	SetBoundaries(gs_1);
	Side(gs,gs_1);
  	//Cp(gs,gs_1+1,M-1); Add(gs+1,gs_1,M-1); // summing over the z direction (I hope) // The last point of gs does not need to be set to zero because it only effects the boundaries which are over written with the set boundaries command.
  	//Add(gs,gs_1+jy,M-jy);Add(gs+jy,gs_1,M-jy); // summing over the y direction (I hope)
  	//Add(gs,gs_1+jx,M-jx);Add(gs+jx,gs_1,M-jx); // summing over the x direction (I hope)
  	//Norm(gs,1.0/6.0,M);
	if (dendrimer==1 && forwards ==1 &&(s_to-1)%N_g ==0){		// if the new segment is a node 		
  			Powerinplace(gs,dendfunc-1,M);  // add propagators of the different arms of the dendrimer
  	}
  	Timesinplace(gs,g,M);
}

void ComputePhi(double *phi, double *G1, double *u, double *Gg_f, double *Gg_b, int sol, int pol){
	sumamp=0;
	int namp = 0;
	// calculate boltzman weights for single segments of type pol and sol
	Boltzmann(phi+sol,u+sol); Boltzmann(G1,u+pol); SetBoundaries(G1); 
	
	// check if the values stored in G1 are reasonable.
	double SumG1 = Sum(G1,M);
	//printf("check-5 %e \n",SumG1);
	if (!isfinite(SumG1)){
		solve= false;						
		return;
	}
	if (SumG1> ((double)100000*M)){
		printf("Sum G1 is suspiciously big, lowering u's by %e \n",0.2*log(SumG1/M)); 
		//restart calculation of phi's 
		//lower u values so G1 will be smaller next time
		double delta = 0.2*log(SumG1/M);
		Addconstinplace(u+pol,delta,M);//subtract value from u to prevent overflow
		ComputePhi(phi,G1,u,Gg_f,Gg_b,sol,pol);
		Addconstinplace(u+pol,-delta,M);//restore values of u just incase
		return;	
	}
	// start calculation propegator
	if(dendrimer ==1){ 
		forwards = 1;
		if (endnode==1){Times(Gg_f,G1,mask+generations*M,M);} // if the endpoints of the dendrimer should be located on MC segments a mask should be applied to the terminal segment.
		else {Cp(Gg_f,G1,M);}
		for (s=2; s<=N_g*generations+1; s++){ 	//loop over all segments of the dendrimer to calculate the propagator.
			Propagate(Gg_f,G1,s-1,s);
			RemoveBoundaries(Gg_f+M*(s-1));
			//printf("check-4 %d %e \n",s,Sum(Gg_f+M*(s-1),M));
			if ((s-1)%N_g ==0){  									// when the segment is a node
				Timesinplace(Gg_f+M*(s-1),mask+M*(generations-(s-1)/N_g),M);		// apply the appropriate mask
			}
			//printf("check-3 %e \n",(Sum(Gg_f+M*(s-1),M)/M));
			if (underflowprotect == 1 && s%underflowinterval==0){	// to prevent over and underflows we adjust the value of the propagator 
				i = generations-ceil((double)s/N_g); // calculate how many segmenst s there are.
				namp = (s/underflowinterval)-1;
				//printf("check-3 %e \n",(Sum(Gg_f+M*(s-1),M)/(Mx*My*Mz)));
				averageG[namp] =1/(Sum(Gg_f+M*(s-1),M)/(Mx*My*Mz));
				//printf("check-2 %d %f \n",s,averageG[namp]);
				if (s<N_g*generations+1){
					sumamp += log(averageG[namp])*dendfunc*pow((dendfunc-1.0),i); 
				}else {sumamp +=log(averageG[namp]);}
				Norm(Gg_f+M*(s-1),averageG[namp],M);
				//printf("check-1.5 %e \n",Sum(Gg_f+M*(s-1),M));
			}
		}
		// Ok done with the forward propagation.starting with backwards propagation and calculation of volume fractions.

		// first calculate the weight for the central segment. As there is only one central segment 
		//SetBoundaries(Gg_f+(N_g*generations-1)*M); // boundaries zijn al in propagator gezet
		Cp(Gg_b,Gg_f+(N_g*generations-1)*M+1,M-1);Add(Gg_b+1,Gg_f+(N_g*generations-1)*M,M-1); // summing over the z direction (I hope) // temporarely use Gg_b for storing sum Gg_f
		Add(Gg_b,Gg_f+(N_g*generations-1)*M+jy,M-jy);Add(Gg_b+jy,Gg_f+(N_g*generations-1)*M,M-jy); // summing over the y direction (I hope)
  		Add(Gg_b,Gg_f+(N_g*generations-1)*M+jx,M-jx);Add(Gg_b+jx,Gg_f+(N_g*generations-1)*M,M-jx); // summing over the x direction (I hope)
  		Norm(Gg_b,1.0/6.0,M);
		Times(phi+pol,Gg_f+N_g*generations*M,Gg_b,M);
		
		RemoveBoundaries(phi+pol);
		GN=Sum(phi+pol,M); // calculate the Zustandsumme based on the central segment.
		printf("GN %e \n",GN);
		// check if GN has a finite value if not modify segment potentials and try again.
		if (isinf(GN)|| GN > 1e300){ //check if GN is Inf or so big it may cause over/underflows later
			printf("Overflow detected \n"); // probably one of the numbers in phi 
			//restart calculation of phi's 
			//lower u values so GN won't be inf next time but order one
			Addconstinplace(u+pol,0.1,M);//subtract value from u to prevent overflow
			ComputePhi(phi,G1,u,Gg_f,Gg_b,sol,pol);
			Addconstinplace(u+pol,-0.1,M);//restore values of u just incase
			return;
		}
		if (isnan(GN)){ //check if GN is Inf or NAN
			printf("Overflow detected \n"); // probably one of the numbers in phi 
			//restart calculation of phi's 
			//lower u values so GN won't be inf next time but order one
			//double delta = -log(1/(Sum(G1,M)/M));
			Addconstinplace(u+pol,1,M);//subtract value from u to prevent overflow
			ComputePhi(phi,G1,u,Gg_f,Gg_b,sol,pol);
			Addconstinplace(u+pol,-1,M);//restore values of u just incase
			return;
		}
		if (GN < 1e-300){ //check for underflow
			printf("Underflow detected \n"); // probably one of the numbers in phi 
			//restart calculation of phi's 
			//raise u values so GN won't be inf next time but order one
			//double delta = -log(1/(Sum(G1,M)/M));
			underflow++;
			if (underflow > 11){
				solve= false;						
				return;
			}else{
				Addconstinplace(u+pol,-0.09123456789,M);//subtract value from u to prevent overflow
				ComputePhi(phi,G1,u,Gg_f,Gg_b,sol,pol);
				Addconstinplace(u+pol,0.09123456789,M);//restore values of u just incase
				return;
			}
		}
		Timesinplace(phi+pol,G1,M);// in the end when we calculate the total volume fraction of polymer phi the weights are devided by G1 so we have to multiply by G1 here.
		Cp(Gg_b,Gg_f+N_g*generations*M,M);
		int switch2=1; forwards = 0;
		for (s=N_g*generations;s>=1;s--){ // go back from the centre over all segments.
			//printf("check-1 %d \n",s);
			switch2++;
			switch2=switch2%2;
			Propagate(Gg_b,G1,switch2+1,2-switch2);
			
			if ((s-1)%N_g ==0 && s !=1-endnode){
				Timesinplace(Gg_b+(1-switch2)*M,mask+M*(generations-(s-1)/N_g),M);// apply mask
			}
			//printf("check1 \n");
			
			//printf("check2 \n");	
			Norm(Gg_f+(s-1)*M,dendfunc*pow((dendfunc-1.0),((double)generations-ceil((double)s/N_g))),M); // to calculate the volume fraction distribution we need to weigh the propegator with the number of segments.
			if (s==1){
				Times(temp,Gg_b+(1-switch2)*M,Gg_f+(s-1)*M,M); // store endpoint distribution in tenmp.
				Add(phi+pol,temp,M);
			}else{
				AddTimes(phi+pol,Gg_b+(1-switch2)*M,Gg_f+(s-1)*M,M);
			}
			if (underflowprotect == 1 && s%underflowinterval==0){	// to prevent over and underflows we adjust the value of the propagator 
				//i = generations-ceil((double)s/N_g); // calculate how many segments s there are.
				namp = (s/underflowinterval)-1;
				Norm(Gg_b+(1-switch2)*M,averageG[namp],M);
				}
			//printf("check3 \n");
			if ((s-1)%N_g ==0 && s !=1){			// add the propagator of the side chain attached to 
				//printf("check3.1 \n");
				Cp(Gg_b+switch2*M,Gg_f+(s-2)*M+1,M-1);Add(Gg_b+switch2*M+1,Gg_f+(s-2)*M,M-1); // summing over the z direction (I hope) // temporarely use Gg_b for storing sum Gg_f the boundaries have already been set in the propagator so we do not need to do it here.`
				Add(Gg_b+switch2*M,Gg_f+(s-2)*M+jy,M-jy);Add(Gg_b+switch2*M+jy,Gg_f+(s-2)*M,M-jy); // summing over the y direction (I hope)
  				Add(Gg_b+switch2*M,Gg_f+(s-2)*M+jx,M-jx);Add(Gg_b+switch2*M+jx,Gg_f+(s-2)*M,M-jx); // summing over the x direction (I hope)
   				Norm(Gg_b+switch2*M,1.0/6.0,M);
  				Powerinplace(Gg_b+switch2*M,dendfunc-2,M);
				Timesinplace(Gg_b+(1-switch2)*M,Gg_b+switch2*M,M);	
			}				
		}
		//printf("check4 \n");
		Dubble(phi+pol,G1,f/GN);
		RemoveBoundaries(phi+pol);
		//printf("check5 sumphi %e \n",Sum(phi+pol,M));
	}else{
		Times(Gg_f,G1,mask,M); //apply mask to the boltmann weights of first polymer segment. 
	
		for (s=2; s<=(N_g)/2; s++) Propagate(Gg_f,G1,s-1,s); // we willen gebruik maken van de symmetrie voor de lus vormende polymeren dus rekenen we propagator maar half uit	
		s = N_g/2;  Zero(phi+pol,M);   
		if (N_g%2 == 1) {  					// if the number of segment is even.
			Cp(Gg_b+(s%2)*M,Gg_f+(s-1)*M,M); 
			Propagate(Gg_b,G1,(s%2)+1,((s+1)%2)+1); 
			AddTimes(phi+pol,Gg_f+(N_g-s-2)*M,Gg_b+((s+1)%2)*M,M);
			Norm(phi+pol,0.5,M);
		}  else Cp(Gg_b+(s%2)*M,Gg_f+(s-1)*M,M);	
		for (s=(N_g+3)/2; s<=N_g; s++) {
			Propagate(Gg_b,G1,((s-1)%2)+1,(s%2)+1);
			AddTimes(phi+pol,Gg_f+(N_g-s)*M,Gg_b+(s%2)*M,M);
		} 
		RemoveBoundaries(Gg_b+(N_g%2)*M);
		GN=Dot(Gg_b+(N_g%2)*M,mask,M); 
		Dubble(phi+pol,G1,2.0*theta/N_g/GN); 
	}
	RemoveBoundaries(phi+pol);
	//printf("check6 \n");
}

void AddG(double *g,double *phitot, double *alpha){int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	addg<<<n_blocks,block_size>>>(g,phitot,alpha,M);
}
void ComputeG(double *g,double *phi, double *u, double *phitot, double *phi_side, double *alpha, int sol, int pol){
	Cp(phitot,phi+sol,M);  Add(phitot,phi+pol,M); Zero(phi_side,M); 
	Side(phi_side,phi+pol); PutAlpha(g,u+sol,phitot,phi_side,chi,0.0);		//could be adefaster using copy and only filling the remaining value with 0.
	Zero(phi_side,M);
	Side(phi_side,phi+sol); PutAlpha(g+M,u+pol,phitot,phi_side,chi,1.0);
	Cp(alpha,g,M); Add(alpha,g+M,M);Norm(alpha,0.5,M);
	AddG(g,phitot,alpha); RemoveBoundaries(g); RemoveBoundaries(g+M);
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
	RemoveBoundaries(phi+sol);	RemoveBoundaries(phi+pol); // remove boundaries so we do not count the free energy of the boundary double
	Uint = chi*Dot(phi_side,phi+pol,M); // calculate the interaction energy
	//Entropy polymer
	entropy = Dot(phi+pol,u+pol,M)-f*log(theta/GN)-sumamp;// the last term sumamp comes from the underflow protection scheme where 
	printf("sumamp %f \n",sumamp);
	printf("entropy molecul %f \n",entropy);
	// add Entropy solvent	
	entropy += Dot(phi+sol,u+sol,M);
	printf("entropy total %f \n",entropy);	
	F_Helmholtz = Uint - entropy; 
	return F_Helmholtz;
}


double SCF() {
	
	
	for (i=0; i<m*m; i++) Aij[i]=0;  // waarom niet met grafische kaart?
	for (i=0; i<m; i++) Ci[i]=0;
	for (i=0; i<m*m; i++) Apij[i]=0;
	if (randstart==false){Zero(x,iv); Zero(x0,iv);
	}else{ // try different starting situation if problem failed to converge
		H_randvector = new double[iv];
		for(i=0;i<iv;i++){
			H_randvector[i] = frand()-0.5; // hope that a different start will avoid the iterations from getting stuck.
		}
		TransferDataToDevice(iv, H_randvector, x);
		for(i=0;i<iv;i++){
			H_randvector[i] = frand()-0.5; // hope that a different start will avoid the iterations from getting stuck.
		}
		TransferDataToDevice(iv, H_randvector, x0);
	}
	
	it=0; j=1; k_diis=1; k=0;
	ComputePhi(phi,G1,u,Gg_f,Gg_b,sol,pol);
	ComputeG(g,phi,u,phitot,phi_side,alpha,sol,pol);
	YplusisCtimesX(x,g,-eta,iv);
	YisAminB(x_x0,x,x0,iv);
	Cp(xR,x,iv); 
	error = sqrt(Dot(g,g,iv));
	printf("DIIS has been notified\n");
	printf("Your guess = %1e \n",error);
	solve = true;
	while (error > tolerance && it < 1000) {
		it++;
		Cp(x0,x,iv); 
		underflow=0;
		ComputePhi(phi,G1,u,Gg_f,Gg_b,sol,pol);
		//printf("check7 \n");
		if (solve == false) {break;} // if the program failed to calculate the free energy there is no use in continuing to iterate therefore it is better to exit the loop and try again with a smaller m.
		ComputeG(g,phi,u,phitot,phi_side,alpha,sol,pol);
		//printf("check8 \n");
		k=it % m; k_diis++; //plek voor laatste opslag
		YplusisCtimesX(x,g,-eta,iv);
		Cp(xR+k*iv,x,iv); YisAminB(x_x0+k*iv,x,x0,iv); 	
		DIIS(x,x_x0,xR,Aij,Apij,Ci,k,m,iv); 
		if (it % j == 0) { 
			error = sqrt(Dot(g,g,iv));			
			if (isfinite(error)==0){solve = false;}
			printf("it = %i error = %1e \n",it,error);
			//j = log(error/tolerance); if (j<1) j=1;
		}
	}
	if (it == 1000) {solve = false;}
	if (solve == false){
		m--;
		printf("Failed to solve the problem Retry with smaller m %d \n",m);
		if (m < 5 || it == 1000) {
			for (i=0; i<Npart*3; i++){
				fprintf (failedparticleconfigs, "%f	",H_oldposlist[i]); //Write xyz coordinate of particle to a file
			}
			fprintf (failedparticleconfigs, "\n"); 
			printf("Failed to solve problem trying new random starting guess");
			randstart = true;
			m=m_user;
			SCF();
		}
		else {SCF();}
	}
	randstart = false;
	//RemoveBoundaries(phi+pol);
	//printf("total phi %e \n",Sum(phi+pol,M));
	return Helmholtz();
}

bool Setmask(double *poslist,double maskvalue){
	double distance;
	int coordinate;
	int maskno =0;
	for (i=0; i<Npart; i++){ 		// loop over all particles
		if (dendrimer == 1){
			if (i == 0){maskno =0;}
			else{maskno = 1+floor(log(ceil((double)i/(double)dendfunc))/log(dendfunc-1));} // calculate to which generation of nodes this node belongs and thus in which mask.
		}

		for (int xpos = ceil(poslist[i*3]-partrad);xpos <=floor(poslist[i*3]+partrad);xpos++){ // loop overall sites close to particle to check if they are within the particle
			for (int ypos = ceil(poslist[i*3+1]-partrad);ypos <=floor(poslist[i*3+1]+partrad);ypos++){
				for (int zpos = ceil(poslist[i*3+2]-partrad);zpos <=floor(poslist[i*3+2]+partrad);zpos++){
					distance = pow(pow(xpos-poslist[i*3],2.0)+pow(ypos-poslist[i*3+1],2.0)+pow(zpos-poslist[i*3+2],2.0),0.5); //calculate distance between the centre of the particle and the lattice site 
					if (distance <= partrad){										// if the distance between the particle and the lattice site is less than the radius of the particle
						coordinate = ((xpos+Mx)%Mx+1)*jx+((ypos+My)%My+1)*jy+(zpos+Mz)%Mz+1;				// determine the position on the lattice after applying periodic boundary
						if (H_mask[coordinate+maskno*M]!=maskvalue) {H_mask[coordinate+maskno*M]=maskvalue;}						// if the value is not yet set in the mask set value to one.// for now this check does not checke for overlap between different generations of the 
						else {
							printf("Mask value for Particle %d has already been set . \n",i);
							if (allowoverlap == false){
							return false;}									// if particles overlap return false
						}
					}
				}
			}
		}
	}
		
	return true;
}
int Randint(int i){						// Should return a random integer ranging from 0 to i
	uint Rndmax =  4294967295-(4294967295%(i+1)+1);		
	
	while (true){
		uint rndnum = irand();
		if (rndnum <= Rndmax) {return((int)(rndnum%(i+1)));}
	}
}
void VTKoutput(const char *basename,double *vector){
	FILE * vtkfile; // file to store density profiles
	baselength = strlen(basename); //determine the length of the base of the file name
	
	numlength = ceil(log10(NMCsteps/vtkstore+1));
	int namesize=baselength+4+numlength; // 4 characters are needed for the file extension. 
	char *filename;
	filename = new char[baselength+4+numlength]; 
	for (i =0;i<baselength;i++) {filename[i] = basename[i];} // copy filename
	i = mcs/vtkstore;
	for (j=1;j<=numlength;j++){
		l = floor(i/pow(10,(numlength-j)));
		filename[j+baselength-1] =digits[l];
		i=i-l*(pow(10,(numlength-j)));
	}	
	filename[namesize-4] ='.';filename[namesize-3] ='v';filename[namesize-2] ='t';filename[namesize-1] ='k';filename[namesize]='\0'; // declare file extension. for now in a rather ugly way.
	vtkfile = fopen(filename,"w");
	TransferDataToHost(H_vector,vector);
	fprintf (vtkfile,"# vtk DataFile Version 3.0\nvtk output\nASCII\nDATASET STRUCTURED_POINTS\nDIMENSIONS %i %i %i\n",Mx,My,Mz);
	fprintf (vtkfile,"SPACING 1 1 1\nORIGIN 0 0 0\nPOINT_DATA %i\n",Mx*My*Mz);
	fprintf (vtkfile,"SCALARS SFBox_profile double\nLOOKUP_TABLE default\n");
	for (i=1;i<=Mx;i++){
		for (j=1;j<=My;j++){
			for (l=1;l<=Mz;l++){
				fprintf (vtkfile,"%e\n",H_vector[i*jx+j*jy+l]);
			}		
		}		
	}		
	fclose (vtkfile);
}		
	
main() { 

	// initialize graphics card
	cudaDeviceReset();  
	stat = cublasCreate(&handle); if (stat !=CUBLAS_STATUS_SUCCESS) {printf("CUBLAS failed \n");}
	cublasSetPointerMode(handle, CUBLAS_POINTER_MODE_DEVICE);
	
	// calculate theta and incase dendrimer required number of nodes.
	if (dendrimer == 1){ // if dendrimers 
		int Nnode = 1;
		int Nnodeprev = dendfunc;
		Nchain = dendfunc*N_g;
		for (i = 2; i <=generations; i++){
			Nnode = Nnode + Nnodeprev;
			Nnodeprev = Nnodeprev *(dendfunc-1);
			Nchain += Nnodeprev*N_g;
		}
		if (endnode == 1){Nnode = Nnode + Nnodeprev;}
		theta = f*(Nchain+1);
		if (Nnode != Npart) {printf("The number of nodes in the dendrimer does not match the number of monte carlo sites. N MC sites %d NNodes %d \n",Npart,Nnode);// check if the number of nodes equals the number of branch points  + optionally the number of end points of the dendrimer
			return(0);
		}
	}else{				// polymer loops
		theta = N_g*f;
		//generations =1;
	}
	printf("Theta %f \n",theta);
	jx = (My+2)*(Mz+2); jy = Mz+2; M=(Mx+2)*(My+2)*(Mz+2);
	bx1=Mx;   bxm=1; //periodic
	by1=My;   bym=1; //periodic
	bz1=Mz;   bzm=1; // periodic
	pol = M; sol = 0*M;
	Md = new int[3];
	Md[0]= Mx; Md[1]= My; Md[2]= Mz;
	
	iv=2*M;
	int Ndiv2=N_g/2;
	m = m_user;
	Aij = new double[m*m];
	Ci = new double[m]; 
	Apij = new double[m*m]; 

	H_vector = new double[M];
	averageG = new double[(int)ceil((generations*N_g+1)/(double)underflowinterval)];
	
	phi = (double*)AllOnDev(2*M);
	temp = (double*)AllOnDev(M);
	phitot = (double*)AllOnDev(M);
	G1 = (double*)AllOnDev(M); 
	alpha = (double*)AllOnDev(M);
	if (dendrimer ==1){
		Gg_f = (double*)AllOnDev(M*(generations*N_g+2));
		H_mask = new double[M*(generations+endnode)];
		mask = (double*)AllOnDev(M*(generations+endnode)); 
	}
	else{
		Gg_f = (double*)AllOnDev(M*Ndiv2);
		H_mask = new double[M];
		mask = (double*)AllOnDev(M); 
	}
	Gg_b = (double*)AllOnDev(M*2);
	phi_side = (double*)AllOnDev(M);  
	x = (double*)AllOnDev(iv); //U_s + U 
	x0 = (double*)AllOnDev(iv);
	g = (double*)AllOnDev(iv);
	xR = (double*)AllOnDev(m*iv);
	x_x0 = (double*)AllOnDev(m*iv);
	DotResult = (double*)AllOnDev(1);

	u = x;	


//-------------
	//put here your MC code
	//int Nequisteps = 1000; // the number of equilibration steps during which the stepsize may be adjusted.

	FILE * poslistfile; // file to store the particle positions
	FILE * kalfile; // file to store other data such as free energy
	
	
	poslistfile = fopen ("poslist.G","w");
	kalfile = fopen ("kalfile.G","w");
	failedparticleconfigs = fopen ("poslistfailed.G","w");
	
	int index,direction;
	bool accept = true;
	double Free_energy, old_Free_energy, n_moves_new;
	i=0;
	while(i<M*(generations+endnode)){H_mask[i]=0.0;i++;}//initialize H_mask to zero
	if (Setmask(H_poslist,1.0)){		// if there is no overlap calculate free energy
		TransferDataToDevice(M*(generations+endnode), H_mask, mask); // copy mask to graphics card
		Free_energy= SCF();		// It would be nice if SCF would return a boolean so we can check if the calcultion was succesfull
	}
	else {
		printf("The initial particle positions are not allowed to overlap. The program will now exit");
		return(0);
	}
	Setmask(H_poslist,0.0);//set mask back to zero
	old_Free_energy = Free_energy;
	
	int Nnonoverlap=0; // counts the number of moves that were not rejected due to overlap between the nodes.
	n_moves_new = n_moves;
	while (mcs <= NMCsteps) {			
		// Generate new positions
		for (i=1; i<=n_moves;i++){
			direction =Randint(2);// randomly select direction 
			index = (nfixed + Randint(Npart-1-nfixed))*3; //randomly select particle			
			H_poslist[index+direction]=fmod(H_oldposlist[index+direction]-1+2*Randint(1)+Md[direction],Md[direction]); //move particle one step and check for periodic boundaries
		}
		// Update mask
		
		if (Setmask(H_poslist,1.0)){
			TransferDataToDevice(M*(generations+endnode), H_mask, mask); // copy mask to graphics card
			Nnonoverlap++; // count the number of moves that were not rejected due to overlap between the nodes.
			// calculate free energy and thus volume fraction distribution
			m = m_user;
			Free_energy = SCF();
			printf("Free energy itermediate %1f	MCS %d \n",Free_energy,mcs);
			if (Free_energy <= old_Free_energy){accept = true;}
			else{
				if (exp(old_Free_energy-Free_energy)> frand()) {accept = true;}
				else {accept = false;}
			}
			if (accept) {
				old_Free_energy = Free_energy;
				for(int i=0;i<Npart*3;i++){H_oldposlist[i]=H_poslist[i];}
			}
			if (Nequisteps > mcs){
				n_moves_new *= 1+(2*accept-2*desiredaccept)/(pow(Nnonoverlap,pow(0.5,0.5))+1);//pow(0.5,0.5) is thedamping factor, it sets the amount of damping for the adjustment of the size of the Monte-Carlo steps. The higher the value the stronger the damping and thus the smaller the fluctuations. It however also takes longer to reach the optimal value.
				n_moves = (int)n_moves_new;
				if (n_moves < 1){
					n_moves = 1 ;
				}
			}
		}
		else {accept = false;}
		Setmask(H_poslist,0.0);// set mask back to zero.
		// Output
		
		// Positions of the particles
		printf("Free energy %1f \n",Free_energy);
		for (i=0; i<Npart*3; i++){
			fprintf (poslistfile, "%f	",H_oldposlist[i]); //Write xyz coordinate of particle to a file
		}
		fprintf (poslistfile, "\n"); 
		// Kal file
		fprintf (kalfile, "%f	%f	%s	%d	%e	%d \n",old_Free_energy ,Free_energy, accept ? "true" : "false",n_moves,n_moves_new,it); 
		fflush(kalfile); // this line is for debugging file out put
		// vtk files
		if (mcs%vtkstore ==0){
			char name[] ="phi";
			VTKoutput(name,phi+pol);
			Dubble(temp,G1,f/GN); 
			char endname[] ="phiend";
			VTKoutput(endname,temp);
			 
//			i = mcs/vtkstore;
//			for (j=1;j<=numlength;j++){
//				l = floor(i/pow(10,(numlength-j)));
//				vtkfilename[j+baselength-1] =digits[l];
//				i=i-l*(pow(10,(numlength-j)));
//			}
//			
//			vtkfile = fopen (vtkfilename,"w");
//			
//			fprintf (vtkfile,"# vtk DataFile Version 3.0\nvtk output\nASCII\nDATASET STRUCTURED_POINTS\nDIMENSIONS %i %i %i\n",Mx,My,Mz);
//			fprintf (vtkfile,"SPACING 1 1 1\nORIGIN 0 0 0\nPOINT_DATA %i\n",Mx*My*Mz);
//			fprintf (vtkfile,"SCALARS SFBox_profile double\nLOOKUP_TABLE default\n");
//			for (i=1;i<=Mx;i++){
//				for (j=1;j<=My;j++){
//					for (l=1;l<=Mz;l++){
//						fprintf (vtkfile,"%e\n",H_phi[i*jx+j*jy+l]);
//					}		
//				}		
//			}		
//			fclose (vtkfile);
//
		}
		mcs++;
	}

	// clean up by freeing memory
	free(H_vector);free(H_mask);
	cudaFree(phi); cudaFree(x);cudaFree(phitot);cudaFree(phi_side);
	cudaFree(G1);cudaFree(alpha);cudaFree(Gg_f);cudaFree(Gg_b);cudaFree(phi_side);
	cudaFree(x_x0);cudaFree(mask);cudaFree(x0);cudaFree(g);cudaFree(xR);
	cublasDestroy(handle);
	fclose (poslistfile);
	fclose (kalfile);
	fclose (failedparticleconfigs);
	
	return(0);
};

