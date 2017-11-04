#include <stdio.h> 
#include <cuda.h> 
#include "cublas_v2.h"
#include <cuda_runtime.h>



//It works on Tesla 2050;
void GPU_present() {
    int deviceCount =0;
    cuDeviceGetCount(&deviceCount);
    if (deviceCount ==0) {
	printf("There is no device supporting Cuda.\n");
    } 
    else { 
	cudaDeviceReset();
	
	}
}

__global__ void Times(double *P, double *A, double *B, int M){
int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) P[idx]=A[idx]*B[idx];
// P =  A*B;
}

__global__ void add(double *P, double *A, int M){
int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) P[idx]+=A[idx];
// P =  P+B;
}

__global__ void AddTimes(double *P, double *A, double *B, int M){
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) P[idx]+=A[idx]*B[idx];
// P +=  A*B;
}

__global__ void AddTimesF(double *P, double *A, double *B, double F, int M){ 	
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) P[idx]+=F*A[idx]*B[idx];
// P+=(A*B)*F
}

__global__ void TimesF(double *P, double *A, double *B, double F, int M){ 	
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) P[idx]=F*A[idx]*B[idx];
// P=(A*B)*F
}


__global__ void TimesC(double *P,double *A, double C, int M){
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) P[idx]=A[idx]*C;
// P = A*C
}

__global__ void Div(double *P,double *A, int M){
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) {
		if (A[idx]>0) { P[idx]/=A[idx];} else {P[idx]=0.0;}
	}
}

__global__ void Norm(double *P, double C, int M){
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) P[idx] *= C;
}

__global__ void zero(double *P, int M){
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) P[idx] = 0.0;
}

__global__ void Cp (double *P, double *A, int M){
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) P[idx] = A[idx];
}

__global__ void SetBx(double *P,int mmx,int My,int Mz, int bx1, int bxm, int jx, int jy){ 
	int idx;
	int jx_mmx=jx*mmx;
	int jx_bxm=jx*bxm;
	int bx1_jx=bx1*jx;

	int yi =blockIdx.x*blockDim.x+threadIdx.x;
	int zi =blockIdx.y*blockDim.y+threadIdx.y;

	if (yi<My && zi<Mz) {
		idx=jy*yi+zi;
		P[idx]=P[bx1_jx+idx];
		P[jx_mmx+idx]=P[jx_bxm+idx];
	}
}
__global__ void SetBy(double *P,int Mx,int mmy,int Mz, int by1, int bym, int jx, int jy){ 
	int idx;
	int jy_mmy=jy*mmy;
	int jy_bym=jy*bym;
	int jy_by1=jy*by1;

	int xi =blockIdx.x*blockDim.x+threadIdx.x;
	int zi =blockIdx.y*blockDim.y+threadIdx.y;

	if (xi<Mx && zi<Mz) {
		idx=jx*xi+zi;
		P[idx]=P[jy_by1+idx];
		P[jy_mmy+idx]=P[jy_bym+idx];
	}
}
__global__ void SetBz(double *P,int Mx,int My,int mmz, int bz1, int bzm, int jx, int jy){  
	int idx;
	int xi =blockIdx.x*blockDim.x+threadIdx.x;
	int yi =blockIdx.y*blockDim.y+threadIdx.y;

	if (xi<Mx && yi<My) {
		idx=jx*xi+jy*yi;
		P[idx]=P[idx+bz1];
		P[idx+mmz]=P[idx+bzm];
	}
}

__global__ void SwitchBx(double *P,double *Q,int mmx,int My,int Mz, int bx1, int bxm, int jx, int jy){ 
	int idx;
	int jx_mmx=jx*mmx;
	int jx_bxm=jx*bxm;
	int bx1_jx=bx1*jx;

	int yi =blockIdx.x*blockDim.x+threadIdx.x;
	int zi =blockIdx.y*blockDim.y+threadIdx.y;

	if (yi<My && zi<Mz) {
		idx=jy*yi+zi;
		P[idx]=Q[bx1_jx+idx];
		P[jx_mmx+idx]=Q[jx_bxm+idx];
		Q[idx]=P[bx1_jx+idx];
		Q[jx_mmx+idx]=P[jx_bxm+idx];
	}
}
__global__ void SwitchBy(double *P,double *Q,int Mx,int mmy,int Mz, int by1, int bym, int jx, int jy){ 
	int idx;
	int jy_mmy=jy*mmy;
	int jy_bym=jy*bym;
	int jy_by1=jy*by1;

	int xi =blockIdx.x*blockDim.x+threadIdx.x;
	int zi =blockIdx.y*blockDim.y+threadIdx.y;

	if (xi<Mx && zi<Mz) {
		idx=jx*xi+zi;
		P[idx]=Q[jy_by1+idx];
		P[jy_mmy+idx]=Q[jy_bym+idx];
		Q[idx]=P[jy_by1+idx];
		Q[jy_mmy+idx]=P[jy_bym+idx];

	}
}
__global__ void SwitchBz(double *P,double *Q,int Mx,int My,int mmz, int bz1, int bzm, int jx, int jy){  
	int idx;
	int xi =blockIdx.x*blockDim.x+threadIdx.x;
	int yi =blockIdx.y*blockDim.y+threadIdx.y;

	if (xi<Mx && yi<My) {
		idx=jx*xi+jy*yi;
		P[idx]=Q[idx+bz1];
		P[idx+mmz]=Q[idx+bzm];
		Q[idx]=P[idx+bz1];
		Q[idx+mmz]=P[idx+bzm];
	}
}


double *  AllocateMemoryOnDevice(int N) {
	double *X;
	if (cudaSuccess != cudaMalloc((void **) &X, sizeof(double)*N))
	printf("Memory allocation on GPU failed.\n Please reduce size of system and/or chain length\n"); 	    
	return X;
} 

void TransferDataToDevice(int M, double *H, double *D ) { 
	cudaMemcpy(D, H, sizeof(double)*M,cudaMemcpyHostToDevice);
}

void InitializeForward(int s, int M, int Seg, double *Dgi, double *Dg) {
	int block_size=256;
	int n_blocks = M/block_size+(M%block_size == 0 ? 0:1);
	Cp <<< n_blocks,block_size >>> (Dgi+(s-1)*M,Dg+(Seg-1)*M,M);
}

void InitializeForward(int s, int M, int Seg, int Mem, double *Dgi, double *Dg) {
	int block_size=256;
	int n_blocks = M/block_size+(M%block_size == 0 ? 0:1);
	for (int k=0; k<Mem; k++) Cp <<< n_blocks,block_size >>> (Dgi+((s-1)*Mem+k)*M,Dg+(Seg-1)*M,M);
}

void InitializeBackward(int M, int Seg, double *Dgi_inv, double *Dg) {
	int block_size=256;
	int n_blocks = M/block_size+(M%block_size ==0 ? 0:1);
	Cp <<<n_blocks,block_size>>> (Dgi_inv,Dg+(Seg-1)*M,M);
}

void InitializeBackward(int M, int Seg, int Mem, double *Dgi_inv, double *Dg) {
	int block_size=256;
	int n_blocks = M/block_size+(M%block_size ==0 ? 0:1);
	for (int k=0; k<Mem; k++) Cp <<<n_blocks,block_size>>> (Dgi_inv+k*M,Dg+(Seg-1)*M,M);
}


void Composition(int s, int M, int Seg, double *Dgi_inv, double *Dgi, double *Dphi){
	int block_size=256;
	int n_blocks=M/block_size + (M%block_size == 0 ? 0:1); 
	AddTimes <<<n_blocks,block_size>>> (Dphi+(Seg-1)*M,Dgi_inv  ,Dgi+(s-1)*M,M); 
}

void Composition(int s, int M, int Seg, int Mem, double *Dgi_inv, double *Dgi, double *Dphi){
	int block_size=256;
	int n_blocks=M/block_size + (M%block_size == 0 ? 0:1); 
	for (int k=0; k<Mem; k++) AddTimes <<<n_blocks,block_size>>> (Dphi+(Seg-1)*M,Dgi_inv+k*M,Dgi+((s-1)*Mem+k)*M,M); 
}

void CorrectDoubleCounting(int M,int Seg,double *Dphi,double *Dg){
	int block_size=256;
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	Div <<<n_blocks,block_size>>>(Dphi +(Seg-1)*M,Dg +(Seg-1)*M,M);
}

void CorrectDoubleCounting(int M,int s,int Seg,double *Dphi,double *Dg){
	int block_size=256;
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	Div <<<n_blocks,block_size>>>(Dphi +(s-1)*M,Dg +(Seg-1)*M,M);
}

void NormPhi(int M,int Seg,double* Dphi,double C){
	int block_size=256;
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	Norm<<<n_blocks,block_size>>>(Dphi+(Seg-1)*M,C,M);
}

void Zero(int M,double* DX){
	int block_size=256;
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	zero<<<n_blocks,block_size>>>(DX,M);
}

void Add(int s, int M, int Seg, double* DX, double* Dphi){
	int block_size=256;
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	add<<<n_blocks,block_size>>>(Dphi+(Seg-1)*M,DX +(s-1)*M,M);
}

void TransferDataToHost(int M, double *H, double *D) {
	cudaMemcpy(H, D, sizeof(double)*M,cudaMemcpyDeviceToHost);
}

void SetBoundaries(double *P,int Mx,int My,int Mz,int bx1,int bxm,int by1,int bym,int bz1,int bzm) {
	int jx =(My+2)*(Mz+2);
	int jy = Mz+2;
	dim3 dimBlock(16,16);
	dim3 dimGridz((Mx+2+dimBlock.x-1)/dimBlock.x,(My+2+dimBlock.y-1)/dimBlock.y);
	dim3 dimGridy((Mx+2+dimBlock.x-1)/dimBlock.x,(Mz+2+dimBlock.y-1)/dimBlock.y);
	dim3 dimGridx((My+2+dimBlock.x-1)/dimBlock.x,(Mz+2+dimBlock.y-1)/dimBlock.y);
	SetBx<<<dimGridx,dimBlock>>>(P,Mx+1,My+2,Mz+2,bx1,bxm,jx,jy);
	SetBy<<<dimGridy,dimBlock>>>(P,Mx+2,My+1,Mz+2,by1,bym,jx,jy);
	SetBz<<<dimGridz,dimBlock>>>(P,Mx+2,My+2,Mz+1,bz1,bzm,jx,jy);
}

void SwitchBounds(double *P,double *Q,int Mx,int My,int Mz,int bx1,int bxm,int by1,int bym,int bz1,int bzm) {
	int jx =(My+2)*(Mz+2);
	int jy = Mz+2;
	dim3 dimBlock(16,16);
	dim3 dimGridz((Mx+2+dimBlock.x-1)/dimBlock.x,(My+2+dimBlock.y-1)/dimBlock.y);
	dim3 dimGridy((Mx+2+dimBlock.x-1)/dimBlock.x,(Mz+2+dimBlock.y-1)/dimBlock.y);
	dim3 dimGridx((My+2+dimBlock.x-1)/dimBlock.x,(Mz+2+dimBlock.y-1)/dimBlock.y);
	SwitchBx<<<dimGridx,dimBlock>>>(P,Q,Mx+1,My+2,Mz+2,bx1,bxm,jx,jy);
	SwitchBy<<<dimGridy,dimBlock>>>(P,Q,Mx+2,My+1,Mz+2,by1,bym,jx,jy);
	SwitchBz<<<dimGridz,dimBlock>>>(P,Q,Mx+2,My+2,Mz+1,bz1,bzm,jx,jy);
}

double Norm2(int M, int s, double *Da) {
	int s1=(s-1)*M;
	double *da = &Da[s1];
	cublasHandle_t handle;
	cublasCreate(&handle);
	double result[1]; 
	cublasDnrm2(handle,M,da,1,result); 
	cublasDestroy(handle);
	return result[0];
}

void Forward(int s, int M, int Seg, double *Dgi, double *Dg, double S, double f, int *Info) { 
	//second order....
	double f1=exp(f);
	double f_1=exp(-f);	
	//double fnorm = 1.0/(4.0 + f1 + f_1); //force not yet implemented
	
	int Mx = Info[2];
	int My = Info[3];
	int Mz = Info[4];
	int jx = (My+2)*(Mz+2);
	int jy = Mz+2;
	int jz = 1;
	int s2=(s-2)*M*6;
	int s1=(s-1)*M*6;
	int sm=(Seg-1)*M;
	int block_size=256; 
	int n_blocks;
	double *gs0 = &Dgi[s1];	//[1][s];
	double *gs1 = &Dgi[M+s1];	//[1+M][s];
	double *gs2 = &Dgi[2*M+s1];	//[1+2*M][s];
	double *gs3 = &Dgi[3*M+s1];	//[1+3*M][s];
	double *gs4 = &Dgi[4*M+s1];	//[1+4*M][s];
	double *gs5 = &Dgi[5*M+s1];	//[1+5*M][s];

	double *gz0 = &Dgi[s2];	//[1][s-1];
	double *gz1 = &Dgi[M+s2];	//[1+M][s-1];
	double *gz2 = &Dgi[2*M+s2];	//[1+2*M][s-1];
	double *gz3 = &Dgi[3*M+s2];	//[1+3*M][s-1];
	double *gz4 = &Dgi[4*M+s2];	//[1+4*M][s-1];
	double *gz5 = &Dgi[5*M+s2];	//[1+5*M][s-1];

	SetBoundaries(gz0,Mx,My,Mz,Info[5],Info[6],Info[7],Info[8],Info[9],Info[10]);
	SetBoundaries(gz1,Mx,My,Mz,Info[5],Info[6],Info[7],Info[8],Info[9],Info[10]);
	SetBoundaries(gz2,Mx,My,Mz,Info[5],Info[6],Info[7],Info[8],Info[9],Info[10]);
	SetBoundaries(gz3,Mx,My,Mz,Info[5],Info[6],Info[7],Info[8],Info[9],Info[10]);
	SetBoundaries(gz4,Mx,My,Mz,Info[5],Info[6],Info[7],Info[8],Info[9],Info[10]);
	SetBoundaries(gz5,Mx,My,Mz,Info[5],Info[6],Info[7],Info[8],Info[9],Info[10]);
	if (Info[5]<3) SwitchBounds(gz0,gz5,Mx,My,Mz,Info[5],Info[6],Info[7],Info[8],Info[9],Info[10]);
	if (Info[7]<3) SwitchBounds(gz1,gz4,Mx,My,Mz,Info[5],Info[6],Info[7],Info[8],Info[9],Info[10]);
	if (Info[9]<3) SwitchBounds(gz2,gz3,Mx,My,Mz,Info[5],Info[6],Info[7],Info[8],Info[9],Info[10]);

	switch (Info[1]) {
	case 1: 
		if (S<0) {
			n_blocks=(M-jx)/block_size + ((M-jx)%block_size == 0 ? 0:1);
			   Times<<<n_blocks,block_size>>>(gs0+jx,gz0,Dg+sm+jx,M-jx); 
			AddTimes<<<n_blocks,block_size>>>(gs0+jx,gz1,Dg+sm+jx,M-jx); 
			AddTimes<<<n_blocks,block_size>>>(gs0+jx,gz2,Dg+sm+jx,M-jx); 
			AddTimes<<<n_blocks,block_size>>>(gs0+jx,gz3,Dg+sm+jx,M-jx); 
			AddTimes<<<n_blocks,block_size>>>(gs0+jx,gz4,Dg+sm+jx,M-jx); 
			AddTimes<<<n_blocks,block_size>>>(gs0+jx,gz5,Dg+sm+jx,M-jx); 
 
        		   Times<<<n_blocks,block_size>>>(gs5,gz0+jx,Dg+sm,M-jx);
			AddTimes<<<n_blocks,block_size>>>(gs5,gz1+jx,Dg+sm,M-jx);
			AddTimes<<<n_blocks,block_size>>>(gs5,gz2+jx,Dg+sm,M-jx);    
			AddTimes<<<n_blocks,block_size>>>(gs5,gz3+jx,Dg+sm,M-jx);    
			AddTimes<<<n_blocks,block_size>>>(gs5,gz4+jx,Dg+sm,M-jx);    
			AddTimes<<<n_blocks,block_size>>>(gs5,gz5+jx,Dg+sm,M-jx);

			n_blocks=(M-jy)/block_size + ((M-jy)%block_size == 0 ? 0:1);
      		  	   Times<<<n_blocks,block_size>>>(gs1+jy,gz0,Dg+sm+jy,M-jy); 
			AddTimes<<<n_blocks,block_size>>>(gs1+jy,gz1,Dg+sm+jy,M-jy); 
			AddTimes<<<n_blocks,block_size>>>(gs1+jy,gz2,Dg+sm+jy,M-jy); 
			AddTimes<<<n_blocks,block_size>>>(gs1+jy,gz3,Dg+sm+jy,M-jy);
			AddTimes<<<n_blocks,block_size>>>(gs1+jy,gz4,Dg+sm+jy,M-jy); 
			AddTimes<<<n_blocks,block_size>>>(gs1+jy,gz5,Dg+sm+jy,M-jy); 
 
        		   Times<<<n_blocks,block_size>>>(gs4,gz0+jy,Dg+sm,M-jy); 
			AddTimes<<<n_blocks,block_size>>>(gs4,gz1+jy,Dg+sm,M-jy);
			AddTimes<<<n_blocks,block_size>>>(gs4,gz2+jy,Dg+sm,M-jy);    
			AddTimes<<<n_blocks,block_size>>>(gs4,gz3+jy,Dg+sm,M-jy);    
			AddTimes<<<n_blocks,block_size>>>(gs4,gz4+jy,Dg+sm,M-jy);    
			AddTimes<<<n_blocks,block_size>>>(gs4,gz5+jy,Dg+sm,M-jy); 

			n_blocks=(M-jz)/block_size + ((M-jz)%block_size == 0 ? 0:1);
        		   Times<<<n_blocks,block_size>>>(gs2+jz,gz0,Dg+sm+jz,M-jz); 
			AddTimes<<<n_blocks,block_size>>>(gs2+jz,gz1,Dg+sm+jz,M-jz); 
			AddTimes<<<n_blocks,block_size>>>(gs2+jz,gz2,Dg+sm+jz,M-jz);
			AddTimes<<<n_blocks,block_size>>>(gs2+jz,gz3,Dg+sm+jz,M-jz);
			AddTimes<<<n_blocks,block_size>>>(gs2+jz,gz4,Dg+sm+jz,M-jz); 
			AddTimes<<<n_blocks,block_size>>>(gs2+jz,gz5,Dg+sm+jz,M-jz); 

       	          	Times<<<n_blocks,block_size>>>(gs3,gz0+jz,Dg+sm,M-jz);    
			AddTimes<<<n_blocks,block_size>>>(gs3,gz1+jz,Dg+sm,M-jz);
			AddTimes<<<n_blocks,block_size>>>(gs3,gz2+jz,Dg+sm,M-jz);
			AddTimes<<<n_blocks,block_size>>>(gs3,gz3+jz,Dg+sm,M-jz);    
			AddTimes<<<n_blocks,block_size>>>(gs3,gz4+jz,Dg+sm,M-jz);    
			AddTimes<<<n_blocks,block_size>>>(gs3,gz5+jz,Dg+sm,M-jz);    

			n_blocks=6*M/block_size + ((6*M)%block_size == 0 ? 0:1);
  			Norm<<<n_blocks,block_size>>>(gs0,1.0/6.0,6*M);
		} else {
			double L=(1.0-S)/4.0;
			double SoverL=S/L;

			n_blocks=(M-jx)/block_size + ((M-jx)%block_size == 0 ? 0:1);
			  TimesF<<<n_blocks,block_size>>>(gs0+jx,gz0,Dg+sm+jx,SoverL,M-jx); 
			AddTimes<<<n_blocks,block_size>>>(gs0+jx,gz1,Dg+sm+jx,M-jx); 
			AddTimes<<<n_blocks,block_size>>>(gs0+jx,gz2,Dg+sm+jx,M-jx); 
			AddTimes<<<n_blocks,block_size>>>(gs0+jx,gz3,Dg+sm+jx,M-jx); 
			AddTimes<<<n_blocks,block_size>>>(gs0+jx,gz4,Dg+sm+jx,M-jx); 

        		    Times<<<n_blocks,block_size>>>(gs5,gz1+jx,Dg+sm,M-jx);    
			 AddTimes<<<n_blocks,block_size>>>(gs5,gz2+jx,Dg+sm,M-jx);    
			 AddTimes<<<n_blocks,block_size>>>(gs5,gz3+jx,Dg+sm,M-jx);    
			 AddTimes<<<n_blocks,block_size>>>(gs5,gz4+jx,Dg+sm,M-jx);    
	        	AddTimesF<<<n_blocks,block_size>>>(gs5,gz5+jx,Dg+sm,SoverL,M-jx);  

			n_blocks=(M-jy)/block_size + ((M-jy)%block_size == 0 ? 0:1);
       	 	    	Times<<<n_blocks,block_size>>>(gs1+jy,gz0,Dg+sm+jy,M-jy); 
			AddTimesF<<<n_blocks,block_size>>>(gs1+jy,gz1,Dg+sm+jy,SoverL,M-jy); 
			 AddTimes<<<n_blocks,block_size>>>(gs1+jy,gz2,Dg+sm+jy,M-jy); 
			 AddTimes<<<n_blocks,block_size>>>(gs1+jy,gz3,Dg+sm+jy,M-jy); 
			 AddTimes<<<n_blocks,block_size>>>(gs1+jy,gz5,Dg+sm+jy,M-jy);
 
        		    Times<<<n_blocks,block_size>>>(gs4,gz0+jy,Dg+sm,M-jy);    
			 AddTimes<<<n_blocks,block_size>>>(gs4,gz2+jy,Dg+sm,M-jy);    
			 AddTimes<<<n_blocks,block_size>>>(gs4,gz3+jy,Dg+sm,M-jy);    
			AddTimesF<<<n_blocks,block_size>>>(gs4,gz4+jy,Dg+sm,SoverL,M-jy);    
			 AddTimes<<<n_blocks,block_size>>>(gs4,gz5+jy,Dg+sm,M-jy);  

			n_blocks=(M-jz)/block_size + ((M-jz)%block_size == 0 ? 0:1);
       	 	    	Times<<<n_blocks,block_size>>>(gs2+jz,gz0,Dg+sm+jz,M-jz); 
			 AddTimes<<<n_blocks,block_size>>>(gs2+jz,gz1,Dg+sm+jz,M-jz); 
			AddTimesF<<<n_blocks,block_size>>>(gs2+jz,gz2,Dg+sm+jz,SoverL,M-jz); 
			 AddTimes<<<n_blocks,block_size>>>(gs2+jz,gz4,Dg+sm+jz,M-jz); 
			 AddTimes<<<n_blocks,block_size>>>(gs2+jz,gz5,Dg+sm+jz,M-jz); 

      			    Times<<<n_blocks,block_size>>>(gs3,gz0+jz,Dg+sm,M-jz);    
			 AddTimes<<<n_blocks,block_size>>>(gs3,gz1+jz,Dg+sm,M-jz);    
			AddTimesF<<<n_blocks,block_size>>>(gs3,gz3+jz,Dg+sm,SoverL,M-jz);    
			 AddTimes<<<n_blocks,block_size>>>(gs3,gz4+jz,Dg+sm,M-jz);    
			 AddTimes<<<n_blocks,block_size>>>(gs3,gz5+jz,Dg+sm,M-jz);    
 
			n_blocks=(6*M)/block_size + ((6*M)%block_size == 0 ? 0:1);
  			Norm<<<n_blocks,block_size>>>(gs0,L,6*M);
		}
		break;
	default: 
		printf("Unknown LatticeType in GPU module, 2nd order.\n");
		break; 
	}
}

void Forward(int s, int M, int Seg, double *Dgi, double *Dg, double f, int LB, int UB) { //LB = info 5 , UP = info 6	
	double fnorm = 1.0/6.0;
	int jz = 1;
	int s2=(s-2)*M;
	int s1=(s-1)*M;
	int sm=(Seg-1)*M;
	
	//SetBoundaries(Dgi+s2,M,LB,UB);
	int block_size=256; 
	int n_blocks; 
	n_blocks=M/block_size + (M%block_size == 0 ? 0:1);
       SetBz<<<n_blocks,block_size>>> (Dgi+s2,0,0,M, LB, UB, 0, 0);

 
	n_blocks=(M-1)/block_size + ((M-1)%block_size == 0 ? 0:1);
	Times <<<n_blocks,block_size>>>  (Dgi+s1+jz,Dgi+s2,   Dg+sm+jz,   M-jz);
	AddTimes <<<n_blocks,block_size>>>  (Dgi+s1,   Dgi+s2+jz,Dg+sm,   M-jz);
	n_blocks=M/block_size + (M%block_size == 0 ? 0:1);
	AddTimesF <<<n_blocks,block_size>>>  (Dgi+s1,   Dgi+s2,Dg+sm, 4.0,   M);
	Norm <<<n_blocks,block_size>>> (Dgi+s1,fnorm,M);	

}

void Forward(int s, int M, int Seg, double *Dgi, double *Dg, double f, int *Info) { 
	if (Info[0]==1) {Forward(s, M, Seg, Dgi, Dg, f, Info[5], Info[6]); } else {
	double f1=exp(f);
	double f_1=exp(-f);	
	double fnorm = 1.0/(4.0 + f1 + f_1);
	int Mx = Info[2];
	int My = Info[3];
	int Mz = Info[4];
	int jx = (My+2)*(Mz+2);
	int jy = Mz+2;
	int jz = 1;
	int s2=(s-2)*M;
	int s1=(s-1)*M;
	int sm=(Seg-1)*M;
	SetBoundaries(Dgi+s2,Mx,My,Mz,Info[5],Info[6],Info[7],Info[8],Info[9],Info[10]);

	int block_size=256; 
	int n_blocks;  

	switch (Info[1]) {
	case 1: //Simple cubic
		n_blocks=(M-jx)/block_size + ((M-jx)%block_size == 0 ? 0:1);
		Times    <<<n_blocks,block_size>>> (Dgi+s1+jx,Dgi+s2,   Dg+sm+jx,M-jx);
		AddTimes <<<n_blocks,block_size>>> (Dgi+s1   ,Dgi+s2+jx,Dg+sm   ,M-jx);
		n_blocks=(M-jy)/block_size + ((M-jy)%block_size == 0 ? 0:1);
		AddTimes <<<n_blocks,block_size>>> (Dgi+s1+jy,Dgi+s2,   Dg+sm+jy,M-jy);
		AddTimes <<<n_blocks,block_size>>> (Dgi+s1,   Dgi+s2+jy,Dg+sm   ,M-jy);
		n_blocks=(M-1)/block_size + ((M-1)%block_size == 0 ? 0:1);
        	if (f !=0)  {
			  AddTimesF <<<n_blocks,block_size>>> (Dgi+s1+jz,Dgi+s2,   Dg+sm+jz,f1, M-jz);
			  AddTimesF <<<n_blocks,block_size>>> (Dgi+s1,   Dgi+s2+jz,Dg+sm,   f_1,M-jz); 
		} else {
			  AddTimes <<<n_blocks,block_size>>>  (Dgi+s1+jz,Dgi+s2,   Dg+sm+jz,M-jz);
			  AddTimes <<<n_blocks,block_size>>>  (Dgi+s1,   Dgi+s2+jz,Dg+sm,   M-jz);
		}
				n_blocks=M/block_size + (M%block_size == 0 ? 0:1);
			  Norm <<<n_blocks,block_size>>> (Dgi+s1,fnorm,M);	
		break;	
	case 2: //FCC
		n_blocks=(M-jx-jy)/block_size + ((M-jx-jy)%block_size == 0 ? 0:1);
           	Times<<<n_blocks,block_size>>>(   Dgi+s1+jx+jy,Dgi+s2,       Dg+sm+jx+jy,M-jx-jy);
		AddTimes<<<n_blocks,block_size>>>(Dgi+s1,      Dgi+s2+jx+jy, Dg+sm,      M-jx-jy);
		AddTimes<<<n_blocks,block_size>>>(Dgi+s1+jx,   Dgi+s2+jy,    Dg+sm+jx,   M-jx-jy);
		AddTimes<<<n_blocks,block_size>>>(Dgi+s1+jy,   Dgi+s2+jx,    Dg+sm+jy,   M-jx-jy);
		n_blocks=(M-jx-jy)/block_size + ((M-jx-jy)%block_size == 0 ? 0:1);
		if (f !=0) {
			AddTimesF<<<n_blocks,block_size>>>(Dgi+s1+jy+jz,Dgi+s2,       Dg+sm+jy+jz,f1,M-jy-jz);
        		AddTimesF<<<n_blocks,block_size>>>(Dgi+s1,      Dgi+s2+jy+jz, Dg+sm,      f_1,M-jy-jz);
			AddTimesF<<<n_blocks,block_size>>>(Dgi+s1+jy,   Dgi+s2+jz,    Dg+sm+jy,   f_1,M-jy-jz);
			AddTimesF<<<n_blocks,block_size>>>(Dgi+s1+jz,   Dgi+s2+jy,    Dg+sm+jz,   f1,M-jy-jz);
			n_blocks=(M-jx-jz)/block_size + ((M-jx-jz)%block_size == 0 ? 0:1);
			AddTimesF<<<n_blocks,block_size>>>(Dgi+s1+jx+jz,Dgi+s2,       Dg+sm+jx+jz,f1,M-jx-jz);
        		AddTimesF<<<n_blocks,block_size>>>(Dgi+s1,      Dgi+s2+jx+jz, Dg,      f_1,M-jx-jz);
			AddTimesF<<<n_blocks,block_size>>>(Dgi+s1+jx,   Dgi+s2+jz,    Dg+sm+jx,   f_1,M-jx-jz);
			AddTimesF<<<n_blocks,block_size>>>(Dgi+s1+jz,   Dgi+s2+jx    ,Dg+sm+jz,   f1,M-jx-jz);
		} else {
			AddTimes<<<n_blocks,block_size>>>(Dgi+s1+jy+jz,Dgi+s2,       Dg+sm+jy+jz,M-jy-jz);
        		AddTimes<<<n_blocks,block_size>>>(Dgi+s1,      Dgi+s2+jy+jz, Dg+sm,      M-jy-jz);
			AddTimes<<<n_blocks,block_size>>>(Dgi+s1+jy,   Dgi+s2+jz,    Dg+sm+jy,   M-jy-jz);
			AddTimes<<<n_blocks,block_size>>>(Dgi+s1+jz,   Dgi+s2+jy,    Dg+sm+jz,   M-jy-jz);
			n_blocks=(M-jx-jz)/block_size + ((M-jx-jz)%block_size == 0 ? 0:1);
			AddTimes<<<n_blocks,block_size>>>(Dgi+s1+jx+jz,Dgi+s2,       Dg+sm+jx+jz,M-jx-jz);
        		AddTimes<<<n_blocks,block_size>>>(Dgi+s1,      Dgi+s2+jx+jz, Dg+sm,      M-jx-jz);
			AddTimes<<<n_blocks,block_size>>>(Dgi+s1+jx,   Dgi+s2+jz,    Dg+sm+jx,   M-jx-jz);
			AddTimes<<<n_blocks,block_size>>>(Dgi+s1+jz,   Dgi+s2+jx    ,Dg+sm+jz,   M-jx-jz);
		}

		n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
       	Norm<<<n_blocks,block_size>>>(Dgi+s1,1/12.0,M);

		break;
	case 3: //HEX
		n_blocks=(M-jx)/block_size + ((M-jx)%block_size == 0 ? 0:1);
		Times<<<n_blocks,block_size>>>(   Dgi+s1+jx,   Dgi+s2,       Dg+sm+jx,   M-jx);
		AddTimes<<<n_blocks,block_size>>>(Dgi+s1,      Dgi+s2+jx,    Dg+sm,      M-jx);
		n_blocks=(M-jx-jy)/block_size + ((M-jx-jy)%block_size == 0 ? 0:1);
		AddTimes<<<n_blocks,block_size>>>(Dgi+s1+jx,   Dgi+s2+jy,    Dg+sm+jx,   M-jx-jy);
		AddTimes<<<n_blocks,block_size>>>(Dgi+s1+jy,   Dgi+s2+jx,    Dg+sm+jy,   M-jy-jx);
		n_blocks=(M-jy)/block_size + ((M-jy)%block_size == 0 ? 0:1);
		AddTimes<<<n_blocks,block_size>>>(Dgi+s1,      Dgi+s2+jy,    Dg+sm,      M-jy);
		AddTimes<<<n_blocks,block_size>>>(Dgi+s1+jy,   Dgi+s2,       Dg+sm+jy ,  M-jy);
		if (f!=0) {
			n_blocks=(M-jx-jz)/block_size + ((M-jx-jz)%block_size == 0 ? 0:1);		
			AddTimesF<<<n_blocks,block_size>>>(Dgi+s1+jz,   Dgi+s2+jx,    Dg+sm+jz,   f1 ,M-jx-jz);
        		AddTimesF<<<n_blocks,block_size>>>(Dgi+s1+jx,   Dgi+s2+jz,    Dg+sm+jx,   f_1,M-jx-jz);
			n_blocks=(M-jy-jz)/block_size + ((M-jy-jz)%block_size == 0 ? 0:1);
			AddTimesF<<<n_blocks,block_size>>>(Dgi+s1+jy,   Dgi+s2+jz,    Dg+sm+jy,   f_1,M-jy-jz);
			AddTimesF<<<n_blocks,block_size>>>(Dgi+s1+jz,   Dgi+s2+jy,    Dg+sm+jz,   f1, M-jy-jz);
			n_blocks=(M-jz)/block_size + ((M-jz)%block_size == 0 ? 0:1);
			AddTimesF<<<n_blocks,block_size>>>(Dgi+s1+jz,   Dgi+s2,       Dg+sm+jz,   f1, M-jz);
			AddTimesF<<<n_blocks,block_size>>>(Dgi+s1,      Dgi+s2+jz,    Dg+sm,      f_1,M-jz);
			
		} else {
			n_blocks=(M-jx-jz)/block_size + ((M-jx-jz)%block_size == 0 ? 0:1);			
			AddTimesF<<<n_blocks,block_size>>>(Dgi+s1+jz,   Dgi+s2+jx,    Dg+sm+jz,   f1, M-jx-jz);
        		AddTimesF<<<n_blocks,block_size>>>(Dgi+s1+jx,   Dgi+s2+jz,    Dg+sm+jx,   f_1,M-jx-jz);
			n_blocks=(M-jy-jz)/block_size + ((M-jy-jz)%block_size == 0 ? 0:1);
			AddTimesF<<<n_blocks,block_size>>>(Dgi+s1+jy,   Dgi+s2+jz,    Dg+sm+jy,   f_1,M-jy-jz);
			AddTimesF<<<n_blocks,block_size>>>(Dgi+s1+jz,   Dgi+s2+jy,    Dg+sm+jz,   f1, M-jy-jz);
			n_blocks=(M-jz)/block_size + ((M-jz)%block_size == 0 ? 0:1);
			AddTimesF<<<n_blocks,block_size>>>(Dgi+s1+jz,   Dgi+s2,       Dg+sm+jz,   f1, M-jz);
			AddTimesF<<<n_blocks,block_size>>>(Dgi+s1,      Dgi+s2+jz,    Dg+sm,      f_1,M-jz);
			
		}
		n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
       		Norm<<<n_blocks,block_size>>>(Dgi+s1,1/12.0,M);

		break;
	default: 
		printf("Unknown LatticeType in GPU module.\n");
		break;
	}} //extra bracket to escape from the call to first order matrix forward;
}

void Backward(int M, int Seg, double *Dgi_inv, double *Dg, double *Dgx, double f, int LB, int UB) { 
	
	double fnorm = 1.0/(6.0);
	int jz = 1;
	int sm=(Seg-1)*M;
	//SetBoundaries(Dgi_inv,M,LB,UB);

       int block_size=256; 
	int n_blocks; 

	n_blocks=M/block_size + (M%block_size == 0 ? 0:1);
	SetBz<<<n_blocks,block_size>>>(Dgi_inv,0,0,M, LB, UB, 0, 0);

	n_blocks=(M-1)/block_size+((M-1)%block_size == 0 ? 0:1);
	Times<<<n_blocks,block_size>>>  (Dgx+jz,Dgi_inv,   Dg+sm+jz,M-jz);
	AddTimes<<<n_blocks,block_size>>>  (Dgx,   Dgi_inv+jz,Dg+sm,   M-jz);
	n_blocks=M/block_size+(M%block_size ==0 ? 0:1);
	AddTimesF<<<n_blocks,block_size>>> (Dgx,   Dgi_inv,Dg+sm,   4.0, M);
	TimesC<<<n_blocks,block_size>>> (Dgi_inv,Dgx,fnorm,M);
}

void Backward(int M, int Seg, double *Dgi_inv, double *Dg, double *Dgx, double f, int *Info) { 
	if (Info[0]==1) {Backward(M, Seg, Dgi_inv, Dg, Dgx, f, Info[5], Info[6]);} else {
	const double f1 = exp(f);
	const double f_1 =  exp(-f);
	double fnorm = 1.0/(4.0 + f1 + f_1);
	int Mx = Info[2];
	int My = Info[3];
	int Mz = Info[4];
	int jx = (My+2)*(Mz+2);
	int jy = (Mz+2);
	int jz = 1;
	int sm=(Seg-1)*M;
	SetBoundaries(Dgi_inv,Mx,My,Mz,Info[5],Info[6],Info[7],Info[8],Info[9],Info[10]);

       int block_size=256; 
	int n_blocks; 

	switch (Info[1]){
	case 1: //simple cubic
	         n_blocks=(M-jx)/block_size+ ((M-jx)%block_size == 0 ? 0:1);	
		  Times   <<<n_blocks,block_size>>> (Dgx+jx,Dgi_inv,   Dg+sm+jx,M-jx);
		  AddTimes<<<n_blocks,block_size>>> (Dgx,   Dgi_inv+jx,Dg+sm,   M-jx);

		  n_blocks=(M-jy)/block_size+ ((M-jy)%block_size == 0 ? 0:1);
		  AddTimes<<<n_blocks,block_size>>> (Dgx+jy,Dgi_inv,   Dg+sm+jy,M-jy);
		  AddTimes<<<n_blocks,block_size>>> (Dgx,   Dgi_inv+jy,Dg+sm,   M-jy);
		  n_blocks=(M-1)/block_size+((M-1)%block_size == 0 ? 0:1);
        	if (f !=0)  {
		  AddTimesF<<<n_blocks,block_size>>> (Dgx+jz,Dgi_inv,   Dg+sm+jz,f1, M-jz);
		  AddTimesF<<<n_blocks,block_size>>> (Dgx,   Dgi_inv+jz,Dg+sm,   f_1,M-jz);
		} else {
		  AddTimes<<<n_blocks,block_size>>>  (Dgx+jz,Dgi_inv,   Dg+sm+jz,M-jz);
		  AddTimes<<<n_blocks,block_size>>>  (Dgx,   Dgi_inv+jz,Dg+sm,   M-jz);
		}
		  n_blocks=M/block_size+(M%block_size ==0 ? 0:1);
		  TimesC<<<n_blocks,block_size>>> (Dgi_inv,Dgx,fnorm,M);	
		break;	
	case 2: //FCC
		n_blocks=(M-jx-jy)/block_size + ((M-jx-jy)%block_size == 0 ? 0:1);
           	Times<<<n_blocks,block_size>>>(   Dgx+jx+jy,Dgi_inv,       Dg+sm+jx+jy,M-jx-jy);
		AddTimes<<<n_blocks,block_size>>>(Dgx,      Dgi_inv+jx+jy, Dg+sm,      M-jx-jy);
		AddTimes<<<n_blocks,block_size>>>(Dgx+jx,   Dgi_inv+jy,    Dg+sm+jx,   M-jx-jy);
		AddTimes<<<n_blocks,block_size>>>(Dgx+jy,   Dgi_inv+jx,    Dg+sm+jy,   M-jx-jy);
		n_blocks=(M-jy-jz)/block_size + ((M-jy-jz)%block_size == 0 ? 0:1);
		if (f !=0)  {
			AddTimesF<<<n_blocks,block_size>>>(Dgx+jy+jz,Dgi_inv,       Dg+sm+jy+jz,f1, M-jy-jz);
        		AddTimesF<<<n_blocks,block_size>>>(Dgx,      Dgi_inv+jy+jz, Dg+sm,      f_1,M-jy-jz);
			AddTimesF<<<n_blocks,block_size>>>(Dgx+jy,   Dgi_inv+jz,    Dg+sm+jy,   f_1,M-jy-jz);
			AddTimesF<<<n_blocks,block_size>>>(Dgx+jz,   Dgi_inv+jy,    Dg+sm+jz,   f1, M-jy-jz);
			n_blocks=(M-jx-jz)/block_size + ((M-jx-jz)%block_size == 0 ? 0:1);
			AddTimesF<<<n_blocks,block_size>>>(Dgx+jx+jz,Dgi_inv,       Dg+sm+jx+jz,f1, M-jx-jz);
        		AddTimesF<<<n_blocks,block_size>>>(Dgx,      Dgi_inv+jx+jz, Dg+sm,      f_1,M-jx-jz);
			AddTimesF<<<n_blocks,block_size>>>(Dgx+jx,   Dgi_inv+jz,    Dg+sm+jx,   f_1,M-jx-jz);
			AddTimesF<<<n_blocks,block_size>>>(Dgx+jz,   Dgi_inv+jx    ,Dg+sm+jz,   f1, M-jx-jz);
		} else {
			AddTimes<<<n_blocks,block_size>>>(Dgx+jy+jz,Dgi_inv,       Dg+sm+jy+jz,M-jy-jz);
        		AddTimes<<<n_blocks,block_size>>>(Dgx,      Dgi_inv+jy+jz, Dg+sm,      M-jy-jz);
			AddTimes<<<n_blocks,block_size>>>(Dgx+jy,   Dgi_inv+jz,    Dg+sm+jy,   M-jy-jz);
			AddTimes<<<n_blocks,block_size>>>(Dgx+jz,   Dgi_inv+jy,    Dg+sm+jz,   M-jy-jz);
			n_blocks=(M-jx-jz)/block_size + ((M-jx-jz)%block_size == 0 ? 0:1);
			AddTimes<<<n_blocks,block_size>>>(Dgx+jx+jz,Dgi_inv,       Dg+sm+jx+jz,M-jx-jz);
        		AddTimes<<<n_blocks,block_size>>>(Dgx,      Dgi_inv+jx+jz, Dg+sm,      M-jx-jz);
			AddTimes<<<n_blocks,block_size>>>(Dgx+jx,   Dgi_inv+jz,    Dg+sm+jx,   M-jx-jz);
			AddTimes<<<n_blocks,block_size>>>(Dgx+jz,   Dgi_inv+jx    ,Dg+sm+jz,   M-jx-jz);
		}
		n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
       		TimesC<<<n_blocks,block_size>>>(Dgi_inv,Dgx,1/12.0,M);
		break;
	case 3: //HEX
		n_blocks=(M-jx)/block_size + ((M-jx)%block_size == 0 ? 0:1);
		Times<<<n_blocks,block_size>>>(   Dgx+jx,   Dgi_inv,       Dg+sm+jx,   M-jx);
		AddTimes<<<n_blocks,block_size>>>(Dgx,      Dgi_inv+jx,    Dg+sm,      M-jx);
		n_blocks=(M-jx-jy)/block_size + ((M-jx-jy)%block_size == 0 ? 0:1);
		AddTimes<<<n_blocks,block_size>>>(Dgx+jx,   Dgi_inv+jy,    Dg+sm+jx,   M-jx-jy);
		AddTimes<<<n_blocks,block_size>>>(Dgx+jy,   Dgi_inv+jx,    Dg+sm+jy,   M-jy-jx);
		n_blocks=(M-jy)/block_size + ((M-jy)%block_size == 0 ? 0:1);
		AddTimes<<<n_blocks,block_size>>>(Dgx,      Dgi_inv+jy,    Dg+sm,      M-jy);
		AddTimes<<<n_blocks,block_size>>>(Dgx+jy,   Dgi_inv,       Dg+sm+jy,   M-jy);
		if (f!=0) {
			n_blocks=(M-jx-jz)/block_size + ((M-jx-jz)%block_size == 0 ? 0:1);		
			AddTimesF<<<n_blocks,block_size>>>(Dgx+jz,   Dgi_inv+jx,    Dg+sm+jz,   f1 ,M-jx-jz);
        		AddTimesF<<<n_blocks,block_size>>>(Dgx+jx,   Dgi_inv+jz,    Dg+sm+jx,   f_1,M-jx-jz);
			n_blocks=(M-jy-jz)/block_size + ((M-jy-jz)%block_size == 0 ? 0:1);
			AddTimesF<<<n_blocks,block_size>>>(Dgx+jy,   Dgi_inv+jz,    Dg+sm+jy,   f_1,M-jy-jz);
			AddTimesF<<<n_blocks,block_size>>>(Dgx+jz,   Dgi_inv+jy,    Dg+sm+jz,   f1, M-jy-jz);
			n_blocks=(M-jz)/block_size + ((M-jz)%block_size == 0 ? 0:1);
			AddTimesF<<<n_blocks,block_size>>>(Dgx+jz,   Dgi_inv,       Dg+sm+jz,   f1, M-jz);
			AddTimesF<<<n_blocks,block_size>>>(Dgx,      Dgi_inv+jz,    Dg+sm,      f_1,M-jz);
			
		} else {
			n_blocks=(M-jx-jz)/block_size + ((M-jx-jz)%block_size == 0 ? 0:1);			
			AddTimesF<<<n_blocks,block_size>>>(Dgx+jz,   Dgi_inv+jx,    Dg+sm+jz,   f1, M-jx-jz);
        		AddTimesF<<<n_blocks,block_size>>>(Dgx+jx,   Dgi_inv+jz,    Dg+sm+jx,   f_1,M-jx-jz);
			n_blocks=(M-jy-jz)/block_size + ((M-jy-jz)%block_size == 0 ? 0:1);
			AddTimesF<<<n_blocks,block_size>>>(Dgx+jy,   Dgi_inv+jz,    Dg+sm+jy,   f_1,M-jy-jz);
			AddTimesF<<<n_blocks,block_size>>>(Dgx+jz,   Dgi_inv+jy,    Dg+sm+jz,   f1, M-jy-jz);
			n_blocks=(M-jz)/block_size + ((M-jz)%block_size == 0 ? 0:1);
			AddTimesF<<<n_blocks,block_size>>>(Dgx+jz,   Dgi_inv,       Dg+sm+jz,   f1, M-jz);
			AddTimesF<<<n_blocks,block_size>>>(Dgx,      Dgi_inv+jz,    Dg+sm,      f_1,M-jz);
			
		}
		n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
       		TimesC<<<n_blocks,block_size>>>(Dgi_inv,Dgx,1/12.0,M);
		break;
	default: 
		printf("Unknown LatticeType in GPU module.\n");
		break;
	}}

}

void Backward(int M, int Seg, double *Dgi_inv, double *Dg, double *Dgx, double S, double f, int *Info) { 
	//second order.
	const double f1 = exp(f);
	const double f_1 =  exp(-f);
	double fnorm = 1.0/(4.0 + f1 + f_1);
	int Mx = Info[2];
	int My = Info[3];
	int Mz = Info[4];
	int jx = (My+2)*(Mz+2);
	int jy = (Mz+2);
	int jz = 1;
	int sm=(Seg-1)*M;
	double *gz0 = &Dgi_inv[0]; 
	double *gz1 = &Dgi_inv[M]; 
	double *gz2 = &Dgi_inv[2*M]; 
	double *gz3 = &Dgi_inv[3*M]; 
	double *gz4 = &Dgi_inv[4*M]; 
	double *gz5 = &Dgi_inv[5*M]; 
			
	SetBoundaries(gz0,Mx,My,Mz,Info[5],Info[6],Info[7],Info[8],Info[9],Info[10]);
	SetBoundaries(gz1,Mx,My,Mz,Info[5],Info[6],Info[7],Info[8],Info[9],Info[10]);
	SetBoundaries(gz2,Mx,My,Mz,Info[5],Info[6],Info[7],Info[8],Info[9],Info[10]);
	SetBoundaries(gz3,Mx,My,Mz,Info[5],Info[6],Info[7],Info[8],Info[9],Info[10]);
	SetBoundaries(gz4,Mx,My,Mz,Info[5],Info[6],Info[7],Info[8],Info[9],Info[10]);
	SetBoundaries(gz5,Mx,My,Mz,Info[5],Info[6],Info[7],Info[8],Info[9],Info[10]);
	if (Info[5]<3) SwitchBounds(gz0,gz5,Mx,My,Mz,Info[5],Info[6],Info[7],Info[8],Info[9],Info[10]);
	if (Info[7]<3) SwitchBounds(gz1,gz4,Mx,My,Mz,Info[5],Info[6],Info[7],Info[8],Info[9],Info[10]);
	if (Info[9]<3) SwitchBounds(gz2,gz3,Mx,My,Mz,Info[5],Info[6],Info[7],Info[8],Info[9],Info[10]);
			
	double *gx0 = &Dgx[0];
	double *gx1 = &Dgx[M];
	double *gx2 = &Dgx[2*M];
	double *gx3 = &Dgx[3*M];
	double *gx4 = &Dgx[4*M];
	double *gx5 = &Dgx[5*M];

        int block_size=256; 
	int n_blocks; 

	switch (Info[1]) {
	case 1:
		if (S<0){
			n_blocks=(M-jx)/block_size + ((M-jx)%block_size == 0 ? 0:1);
			Times<<<n_blocks,block_size>>>(gx0+jx,gz5,Dg+sm+jx,M-jx);
			Times<<<n_blocks,block_size>>>(gx1+jx,gz5,Dg+sm+jx,M-jx);
			Times<<<n_blocks,block_size>>>(gx2+jx,gz5,Dg+sm+jx,M-jx);
			Times<<<n_blocks,block_size>>>(gx3+jx,gz5,Dg+sm+jx,M-jx);
			Times<<<n_blocks,block_size>>>(gx4+jx,gz5,Dg+sm+jx,M-jx);
			Times<<<n_blocks,block_size>>>(gx5+jx,gz5,Dg+sm+jx,M-jx);
 
			AddTimes<<<n_blocks,block_size>>>(gx0,gz0+jx,Dg+sm,M-jx);
			AddTimes<<<n_blocks,block_size>>>(gx1,gz0+jx,Dg+sm,M-jx);
			AddTimes<<<n_blocks,block_size>>>(gx2,gz0+jx,Dg+sm,M-jx);
			AddTimes<<<n_blocks,block_size>>>(gx3,gz0+jx,Dg+sm,M-jx);
			AddTimes<<<n_blocks,block_size>>>(gx4,gz0+jx,Dg+sm,M-jx);
			AddTimes<<<n_blocks,block_size>>>(gx5,gz0+jx,Dg+sm,M-jx); 

			n_blocks=(M-jy)/block_size + ((M-jy)%block_size == 0 ? 0:1);
			AddTimes<<<n_blocks,block_size>>>(gx0+jy,gz4,Dg+sm+jy,M-jy);
			AddTimes<<<n_blocks,block_size>>>(gx1+jy,gz4,Dg+sm+jy,M-jy);
			AddTimes<<<n_blocks,block_size>>>(gx2+jy,gz4,Dg+sm+jy,M-jy);
			AddTimes<<<n_blocks,block_size>>>(gx3+jy,gz4,Dg+sm+jy,M-jy);
			AddTimes<<<n_blocks,block_size>>>(gx4+jy,gz4,Dg+sm+jy,M-jy);
			AddTimes<<<n_blocks,block_size>>>(gx5+jy,gz4,Dg+sm+jy,M-jy);  

      		  	AddTimes<<<n_blocks,block_size>>>(gx0,gz1+jy,Dg+sm,M-jy);
			AddTimes<<<n_blocks,block_size>>>(gx1,gz1+jy,Dg+sm,M-jy);
			AddTimes<<<n_blocks,block_size>>>(gx2,gz1+jy,Dg+sm,M-jy);
			AddTimes<<<n_blocks,block_size>>>(gx3,gz1+jy,Dg+sm,M-jy);
			AddTimes<<<n_blocks,block_size>>>(gx4,gz1+jy,Dg+sm,M-jy);
			AddTimes<<<n_blocks,block_size>>>(gx5,gz1+jy,Dg+sm,M-jy);

			n_blocks=(M-jz)/block_size + ((M-jz)%block_size == 0 ? 0:1);
        		AddTimes<<<n_blocks,block_size>>>(gx0+jz,gz3,Dg+sm+jz,M-jz);
			AddTimes<<<n_blocks,block_size>>>(gx1+jz,gz3,Dg+sm+jz,M-jz);
			AddTimes<<<n_blocks,block_size>>>(gx2+jz,gz3,Dg+sm+jz,M-jz);
			AddTimes<<<n_blocks,block_size>>>(gx3+jz,gz3,Dg+sm+jz,M-jz);
			AddTimes<<<n_blocks,block_size>>>(gx4+jz,gz3,Dg+sm+jz,M-jz);
			AddTimes<<<n_blocks,block_size>>>(gx5+jz,gz3,Dg+sm+jz,M-jz); 

       			AddTimes<<<n_blocks,block_size>>>(gx0,gz2+jz,Dg+sm,M-jz);
			AddTimes<<<n_blocks,block_size>>>(gx1,gz2+jz,Dg+sm,M-jz);
			AddTimes<<<n_blocks,block_size>>>(gx2,gz2+jz,Dg+sm,M-jz);
			AddTimes<<<n_blocks,block_size>>>(gx3,gz2+jz,Dg+sm,M-jz);
			AddTimes<<<n_blocks,block_size>>>(gx4,gz2+jz,Dg+sm,M-jz);
			AddTimes<<<n_blocks,block_size>>>(gx5,gz2+jz,Dg+sm,M-jz);   			

			n_blocks=(6*M)/block_size + ((6*M)%block_size == 0 ? 0:1);
        		TimesC<<<n_blocks,block_size>>>(gz0,gx0,1.0/6.0,6*M);
		} else {
			double L=(1.0-S)/4.0;
			double SoverL=S/L;

			n_blocks=(M-jx)/block_size + ((M-jx)%block_size == 0 ? 0:1);
        		 Times<<<n_blocks,block_size>>>(gx1+jx,gz5,Dg+sm+jx,M-jx);
       		 	Times<<<n_blocks,block_size>>>(gx2+jx,gz5,Dg+sm+jx,M-jx);
        		 Times<<<n_blocks,block_size>>>(gx3+jx,gz5,Dg+sm+jx,M-jx);
			 Times<<<n_blocks,block_size>>>(gx4+jx,gz5,Dg+sm+jx,M-jx);
			TimesF<<<n_blocks,block_size>>>(gx5+jx,gz5,Dg+sm+jx,SoverL,M-jx); 

        		  TimesF<<<n_blocks,block_size>>>(gx0,gz0+jx,Dg+sm,SoverL,M-jx);
			AddTimes<<<n_blocks,block_size>>>(gx1,gz0+jx,Dg+sm,M-jx);
			AddTimes<<<n_blocks,block_size>>>(gx2,gz0+jx,Dg+sm,M-jx);
			AddTimes<<<n_blocks,block_size>>>(gx3,gz0+jx,Dg+sm,M-jx); 
			AddTimes<<<n_blocks,block_size>>>(gx4,gz0+jx,Dg+sm,M-jx); 

			n_blocks=(M-jy)/block_size + ((M-jy)%block_size == 0 ? 0:1);
        		 AddTimes<<<n_blocks,block_size>>>(gx0+jy,gz4,Dg+sm+jy,M-jy);
			 AddTimes<<<n_blocks,block_size>>>(gx2+jy,gz4,Dg+sm+jy,M-jy);
			 AddTimes<<<n_blocks,block_size>>>(gx3+jy,gz4,Dg+sm+jy,M-jy);
			AddTimesF<<<n_blocks,block_size>>>(gx4+jy,gz4,Dg+sm+jy,SoverL,M-jy);
			 AddTimes<<<n_blocks,block_size>>>(gx5+jy,gz4,Dg+sm+jy,M-jy); 

        		 AddTimes<<<n_blocks,block_size>>>(gx0,gz1+jy,Dg+sm,M-jy);
			AddTimesF<<<n_blocks,block_size>>>(gx1,gz1+jy,Dg+sm,SoverL,M-jy);
			 AddTimes<<<n_blocks,block_size>>>(gx2,gz1+jy,Dg+sm,M-jy);
			 AddTimes<<<n_blocks,block_size>>>(gx3,gz1+jy,Dg+sm,M-jy);
			 AddTimes<<<n_blocks,block_size>>>(gx5,gz1+jy,Dg+sm,M-jy);

			n_blocks=(M-jz)/block_size + ((M-jz)%block_size == 0 ? 0:1);
        		 AddTimes<<<n_blocks,block_size>>>(gx0+jz,gz3,Dg+sm+jz,M-jz);
			 AddTimes<<<n_blocks,block_size>>>(gx1+jz,gz3,Dg+sm+jz,M-jz);
			AddTimesF<<<n_blocks,block_size>>>(gx3+jz,gz3,Dg+sm+jz,SoverL,M-jz);
			 AddTimes<<<n_blocks,block_size>>>(gx4+jz,gz3,Dg+sm+jz,M-jz);
			 AddTimes<<<n_blocks,block_size>>>(gx5+jz,gz3,Dg+sm+jz,M-jz); 

        		 AddTimes<<<n_blocks,block_size>>>(gx0,gz2+jz,Dg+sm,M-jz);
			 AddTimes<<<n_blocks,block_size>>>(gx1,gz2+jz,Dg+sm,M-jz);
			AddTimesF<<<n_blocks,block_size>>>(gx2,gz2+jz,Dg+sm,SoverL,M-jz);
			 AddTimes<<<n_blocks,block_size>>>(gx4,gz2+jz,Dg+sm,M-jz);
			 AddTimes<<<n_blocks,block_size>>>(gx5,gz2+jz,Dg+sm,M-jz);

			n_blocks=(6*M)/block_size + ((6*M)%block_size == 0 ? 0:1);
        		TimesC<<<n_blocks,block_size>>>(gz0,gx0,L,6*M); 
		}
		break;
	default: 
		printf("Unknown LatticeType in GPU module.\n");
		break;
	}
}

void Propagate(int s_from, int s_to, int M, int Seg, double *Dgi, double *Dg, double f, int *Info) { 
	double f1=exp(f);
	double f_1=exp(-f);	
	double fnorm = 1.0/(4.0 + f1 + f_1);
	int Mx = Info[2];
	int My = Info[3];
	int Mz = Info[4];
	int jx = (My+2)*(Mz+2);
	int jy = Mz+2;
	int jz = 1;
	int s2=(s_from-1)*M;
	int s1=(s_to-1)*M;
	int sm=(Seg-1)*M;
	SetBoundaries(Dgi+s2,Mx,My,Mz,Info[5],Info[6],Info[7],Info[8],Info[9],Info[10]);

	int block_size=256; 
	int n_blocks;  

	switch (Info[1]) {
	case 1: //Simple cubic
		n_blocks=(M-jx)/block_size + ((M-jx)%block_size == 0 ? 0:1);
		Times    <<<n_blocks,block_size>>> (Dgi+s1+jx,Dgi+s2,   Dg+sm+jx,M-jx);
		AddTimes <<<n_blocks,block_size>>> (Dgi+s1   ,Dgi+s2+jx,Dg+sm   ,M-jx);
		n_blocks=(M-jy)/block_size + ((M-jy)%block_size == 0 ? 0:1);
		AddTimes <<<n_blocks,block_size>>> (Dgi+s1+jy,Dgi+s2,   Dg+sm+jy,M-jy);
		AddTimes <<<n_blocks,block_size>>> (Dgi+s1,   Dgi+s2+jy,Dg+sm   ,M-jy);
		n_blocks=(M-1)/block_size + ((M-1)%block_size == 0 ? 0:1);
        	if (f !=0)  {
			  AddTimesF <<<n_blocks,block_size>>> (Dgi+s1+jz,Dgi+s2,   Dg+sm+jz,f1, M-jz);
			  AddTimesF <<<n_blocks,block_size>>> (Dgi+s1,   Dgi+s2+jz,Dg+sm,   f_1,M-jz); 
		} else {
			  AddTimes <<<n_blocks,block_size>>>  (Dgi+s1+jz,Dgi+s2,   Dg+sm+jz,M-jz);
			  AddTimes <<<n_blocks,block_size>>>  (Dgi+s1,   Dgi+s2+jz,Dg+sm,   M-jz);
		}
				n_blocks=M/block_size + (M%block_size == 0 ? 0:1);
			  Norm <<<n_blocks,block_size>>> (Dgi+s1,fnorm,M);	
		break;	
	case 2: //FCC
		n_blocks=(M-jx-jy)/block_size + ((M-jx-jy)%block_size == 0 ? 0:1);
           	Times<<<n_blocks,block_size>>>(   Dgi+s1+jx+jy,Dgi+s2,       Dg+sm+jx+jy,M-jx-jy);
		AddTimes<<<n_blocks,block_size>>>(Dgi+s1,      Dgi+s2+jx+jy, Dg+sm,      M-jx-jy);
		AddTimes<<<n_blocks,block_size>>>(Dgi+s1+jx,   Dgi+s2+jy,    Dg+sm+jx,   M-jx-jy);
		AddTimes<<<n_blocks,block_size>>>(Dgi+s1+jy,   Dgi+s2+jx,    Dg+sm+jy,   M-jx-jy);
		n_blocks=(M-jx-jy)/block_size + ((M-jx-jy)%block_size == 0 ? 0:1);
		if (f !=0) {
			AddTimesF<<<n_blocks,block_size>>>(Dgi+s1+jy+jz,Dgi+s2,       Dg+sm+jy+jz,f1,M-jy-jz);
        		AddTimesF<<<n_blocks,block_size>>>(Dgi+s1,      Dgi+s2+jy+jz, Dg+sm,      f_1,M-jy-jz);
			AddTimesF<<<n_blocks,block_size>>>(Dgi+s1+jy,   Dgi+s2+jz,    Dg+sm+jy,   f_1,M-jy-jz);
			AddTimesF<<<n_blocks,block_size>>>(Dgi+s1+jz,   Dgi+s2+jy,    Dg+sm+jz,   f1,M-jy-jz);
			n_blocks=(M-jx-jz)/block_size + ((M-jx-jz)%block_size == 0 ? 0:1);
			AddTimesF<<<n_blocks,block_size>>>(Dgi+s1+jx+jz,Dgi+s2,       Dg+sm+jx+jz,f1,M-jx-jz);
        		AddTimesF<<<n_blocks,block_size>>>(Dgi+s1,      Dgi+s2+jx+jz, Dg,      f_1,M-jx-jz);
			AddTimesF<<<n_blocks,block_size>>>(Dgi+s1+jx,   Dgi+s2+jz,    Dg+sm+jx,   f_1,M-jx-jz);
			AddTimesF<<<n_blocks,block_size>>>(Dgi+s1+jz,   Dgi+s2+jx    ,Dg+sm+jz,   f1,M-jx-jz);
		} else {
			AddTimes<<<n_blocks,block_size>>>(Dgi+s1+jy+jz,Dgi+s2,       Dg+sm+jy+jz,M-jy-jz);
        		AddTimes<<<n_blocks,block_size>>>(Dgi+s1,      Dgi+s2+jy+jz, Dg+sm,      M-jy-jz);
			AddTimes<<<n_blocks,block_size>>>(Dgi+s1+jy,   Dgi+s2+jz,    Dg+sm+jy,   M-jy-jz);
			AddTimes<<<n_blocks,block_size>>>(Dgi+s1+jz,   Dgi+s2+jy,    Dg+sm+jz,   M-jy-jz);
			n_blocks=(M-jx-jz)/block_size + ((M-jx-jz)%block_size == 0 ? 0:1);
			AddTimes<<<n_blocks,block_size>>>(Dgi+s1+jx+jz,Dgi+s2,       Dg+sm+jx+jz,M-jx-jz);
        		AddTimes<<<n_blocks,block_size>>>(Dgi+s1,      Dgi+s2+jx+jz, Dg+sm,      M-jx-jz);
			AddTimes<<<n_blocks,block_size>>>(Dgi+s1+jx,   Dgi+s2+jz,    Dg+sm+jx,   M-jx-jz);
			AddTimes<<<n_blocks,block_size>>>(Dgi+s1+jz,   Dgi+s2+jx    ,Dg+sm+jz,   M-jx-jz);
		}

		n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
       	Norm<<<n_blocks,block_size>>>(Dgi+s1,1/12.0,M);

		break;
	case 3: //HEX
		n_blocks=(M-jx)/block_size + ((M-jx)%block_size == 0 ? 0:1);
		Times<<<n_blocks,block_size>>>(   Dgi+s1+jx,   Dgi+s2,       Dg+sm+jx,   M-jx);
		AddTimes<<<n_blocks,block_size>>>(Dgi+s1,      Dgi+s2+jx,    Dg+sm,      M-jx);
		n_blocks=(M-jx-jy)/block_size + ((M-jx-jy)%block_size == 0 ? 0:1);
		AddTimes<<<n_blocks,block_size>>>(Dgi+s1+jx,   Dgi+s2+jy,    Dg+sm+jx,   M-jx-jy);
		AddTimes<<<n_blocks,block_size>>>(Dgi+s1+jy,   Dgi+s2+jx,    Dg+sm+jy,   M-jy-jx);
		n_blocks=(M-jy)/block_size + ((M-jy)%block_size == 0 ? 0:1);
		AddTimes<<<n_blocks,block_size>>>(Dgi+s1,      Dgi+s2+jy,    Dg+sm,      M-jy);
		AddTimes<<<n_blocks,block_size>>>(Dgi+s1+jy,   Dgi+s2,       Dg+sm+jy ,  M-jy);
		if (f!=0) {
			n_blocks=(M-jx-jz)/block_size + ((M-jx-jz)%block_size == 0 ? 0:1);		
			AddTimesF<<<n_blocks,block_size>>>(Dgi+s1+jz,   Dgi+s2+jx,    Dg+sm+jz,   f1 ,M-jx-jz);
        		AddTimesF<<<n_blocks,block_size>>>(Dgi+s1+jx,   Dgi+s2+jz,    Dg+sm+jx,   f_1,M-jx-jz);
			n_blocks=(M-jy-jz)/block_size + ((M-jy-jz)%block_size == 0 ? 0:1);
			AddTimesF<<<n_blocks,block_size>>>(Dgi+s1+jy,   Dgi+s2+jz,    Dg+sm+jy,   f_1,M-jy-jz);
			AddTimesF<<<n_blocks,block_size>>>(Dgi+s1+jz,   Dgi+s2+jy,    Dg+sm+jz,   f1, M-jy-jz);
			n_blocks=(M-jz)/block_size + ((M-jz)%block_size == 0 ? 0:1);
			AddTimesF<<<n_blocks,block_size>>>(Dgi+s1+jz,   Dgi+s2,       Dg+sm+jz,   f1, M-jz);
			AddTimesF<<<n_blocks,block_size>>>(Dgi+s1,      Dgi+s2+jz,    Dg+sm,      f_1,M-jz);
			
		} else {
			n_blocks=(M-jx-jz)/block_size + ((M-jx-jz)%block_size == 0 ? 0:1);			
			AddTimesF<<<n_blocks,block_size>>>(Dgi+s1+jz,   Dgi+s2+jx,    Dg+sm+jz,   f1, M-jx-jz);
        		AddTimesF<<<n_blocks,block_size>>>(Dgi+s1+jx,   Dgi+s2+jz,    Dg+sm+jx,   f_1,M-jx-jz);
			n_blocks=(M-jy-jz)/block_size + ((M-jy-jz)%block_size == 0 ? 0:1);
			AddTimesF<<<n_blocks,block_size>>>(Dgi+s1+jy,   Dgi+s2+jz,    Dg+sm+jy,   f_1,M-jy-jz);
			AddTimesF<<<n_blocks,block_size>>>(Dgi+s1+jz,   Dgi+s2+jy,    Dg+sm+jz,   f1, M-jy-jz);
			n_blocks=(M-jz)/block_size + ((M-jz)%block_size == 0 ? 0:1);
			AddTimesF<<<n_blocks,block_size>>>(Dgi+s1+jz,   Dgi+s2,       Dg+sm+jz,   f1, M-jz);
			AddTimesF<<<n_blocks,block_size>>>(Dgi+s1,      Dgi+s2+jz,    Dg+sm,      f_1,M-jz);
			
		}
		n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
       		Norm<<<n_blocks,block_size>>>(Dgi+s1,1/12.0,M);

		break;
	default: 
		printf("Unknown LatticeType in GPU module.\n");
		break;
	}
}


void FreeMemoryOnDevice(double *D) { 
    	cudaFree(D);  
} 

//	Info[0]=3; //Dimensions
//	Info[1]=3; //3=HEX, 2=FCC, 1=simple cubic;
//	Info[2]=n_layers_x; 
//	Info[3]=n_layers_y;
//	Info[4]=n_layers_z; 
	//assumption that compranges =1;	
//	Info[5]=Get_BL(1); //x lower bound
//	Info[6]=Get_BL(2); //x upper bound
//	Info[7]=Get_BL(3); //y lower bound
//	Info[8]=Get_BL(4); //y upper bound
//	Info[9]=Get_BL(5); //z lower bound
//	Info[10]=Get_BL(6);//z upper bound
