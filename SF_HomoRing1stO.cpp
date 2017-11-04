#include "SF_HomoRing1stO.h"
//#include <iostream>
#ifdef CUDA
#include "Cuda_tools.h"
#endif

SF_HomoRing1stO::SF_HomoRing1stO(Text name_,
								   SF_SegmentList* SegQ_,
								   Lat1stO* Lat_,
								   Input* MyInput_)
: SF_Molecule(name_, SegQ_, Lat_, MyInput_)
{
	Lat = Lat_;
	norm_ref=0;
	Seg = Chain->GetSegment(1);

	force_set = (!GetForce()==0);

	if (force_set) {
		Message(fatal,MyInput,
					"Force ensemble not implemented in Homopol1stO; do not know how to deal with sign of force when inversion symmetry is imposed.");
	}

#ifdef CUDA
	if (GetGPU()) {GPUactivated = true;
	}else GPUactivated = false;
#else
	GPUactivated = false;
	if (GetGPU()) {
		Message(literal,"Compile with nvcc (option CUDA=1) to activate GPU... Contact Frans Leermakers. Going classical");
	}
#endif

	if (freedom == thirdGeneration) {
		CreatePhi(constrained);
		CreatePhi(unconstrained);
		CreatePhi(renorm);
	} else if (freedom == secondGeneration) {
		CreatePhi(constrained);
		CreatePhi(unconstrained);
	} else {
		CreatePhi(total);
	}
	n = int(pow(N*3,1.0/3)+0.5);
	if (2*N < 60) n++;
	if (2*N < 30) n++;
	SetPhiBulk(phiBulk);
}

SF_HomoRing1stO::~SF_HomoRing1stO() {
	if (freedom == thirdGeneration) {
		DeletePhi(constrained);
		DeletePhi(unconstrained);
		DeletePhi(renorm);
	} else if (freedom == secondGeneration) {
		DeletePhi(constrained);
		DeletePhi(unconstrained);
	} else {
		DeletePhi(total);
	}
}

void
SF_HomoRing1stO::ReComputePhi() {
	Message(fatal,"ReComputePhi not implemented in SF_HomoRing1stO");
}

Vector
SF_HomoRing1stO::GetLong(){
	Message(fatal,"GetLong not implemented in SF_HomoRing1stO");
	Vector x;
	return x;
}
Vector
SF_HomoRing1stO::GetShort(){
	Message(fatal,"GetShort not implemented in SF_HomoRing1stO");
	Vector x;
	return x;
}
Vector
SF_HomoRing1stO::GetBondOrientation( const DensityPart DensPar) {
	Message(fatal,"GetBondOrientation not implemented in SF_HomoRing1stO");
	Vector x;
	return x;
};

MoleculeType
SF_HomoRing1stO::GetMoleculeType() const {
	return Chain->GetMoleculeType();
}

void
SF_HomoRing1stO::GetLayerAnalysis(Output* Out) {}

void
SF_HomoRing1stO::ComputePhi() {
	Vector G = Seg->GetSWF();
	Lat->SetBoundaries(G);

	if (!saveMemory) {
		if (GPUactivated) CudaMatrix1(total); else Matrix1(total);
	}
	else {
		if (GPUactivated) CudaMatrix1Long(total); else Matrix1Long(total);
	}

}

void
SF_HomoRing1stO::Matrix1(const DensityPart DensPart) {
	int k,z,s; double a,b,c;
	int M_max; if (N/2>M) M_max = 2*M; else M_max = M + N/2;
	int z_min,z_max;
	double gz;
	Vector NORM(0,M);
	Vector phi = Seg->GetPhi(DensPart);
	Vector G = Seg->GetSWF();
	Lat->MakeSafe(G);
	Vector X(-N/2,N/2);
	Vector Y(0,M_max);
	if (N%2==1) {Message(fatal,MyInput,"rings should have even contour lengths, for the time being....");}

	if (norm_ref==0) {
		for (z=-N/2; z<=N/2; z++) X[z]=0; X[0]=1;
		a=0; b=0; c=0;
    	for (s=1; s<N/2; s++) {
			b=X[-s-1]; c=X[-s];
	    	for (z=-s; z<s; z++) {
				a=b; b=c; c=X[z+1];
				X[z] = (a +4.0*b + c)/6.0;
			}
		}
		NORM[0]=0;
		for (z=-N/2; z<N/2; z++) NORM[0]+=X[z]*X[z]; norm_ref=NORM[0];
		} else {
			NORM[0]=norm_ref;
	}

	for (z=1; z<=M; z++) phi[z] = 0;

	for (k=1; k<=M; k++){
		for (z=1; z<=M_max; z++) {
			if (z==k) Y[z] = G[z]; else  Y[z]=0;
		}
		a=0; b=0; c=0;
		for (s=1; s<N/2; s++) {
			z_min = k-s; if (z_min< 1) z_min = 1;
			z_max = k+s; if (z_max >M_max) z_max = M_max;
			b=Y[z_min-1]; c=Y[z_min]; Y[M_max]=Y[M_max-1];
			for (z=z_min; z<z_max; z++) {
				a=b; b=c; c=Y[z+1]; if (z<=M) gz=G[z]/6.0; else gz=1.0/6.0;
				Y[z]=gz*(a+4*b+c);
			}
		}
		NORM[k]=0;
		for (z=1; z<=M_max; z++) NORM[k]+=Y[z]*Y[z];
	}

	for (z=1; z<=M; z++) phi[z] =NORM[z]/NORM[0]*phiBulk;
}

void
SF_HomoRing1stO::CudaMatrix1(const DensityPart DensPart) {
#ifdef CUDA
	int z,s,s_inv,Ndiv2;
	double Chalf=0.5;
	int Info[11];
	int segnr=0;
	Vector phi,G;
	Vector Gi_inv(1,M);
	Lat->GetLatticeInfo(Info);

	Ndiv2=N/2;

	double *Hphi, *Hg, *Hgi_inv;  //pointers to host memory;
	double *Dphi, *Dg, *Dgi, *Dgi_inv, *Dgx; //pointers to device memory;
	Dphi =    (double *)AllocateMemoryOnDevice(M);
	Dg =      (double *)AllocateMemoryOnDevice(M);
	Dgi =     (double *)AllocateMemoryOnDevice(M*Ndiv2);
	Dgi_inv = (double *)AllocateMemoryOnDevice(M);
	Dgx =     (double *)AllocateMemoryOnDevice(M);
	Hphi =    (double *)malloc(sizeof(double)*M);
	Hg =      (double *)malloc(sizeof(double)*M);
	Hgi_inv = (double *)malloc(sizeof(double)*M);


	Seg = Chain->GetSegment(1);
	G = Seg->GetSWF();
	phi = Seg->GetPhi(DensPart);
	for (z=1; z<=M; z++) {
		phi[z] = 0;
		Hphi[z-1] = 0;
		Hg[z-1] =G[z];
	}

	TransferDataToDevice(M,Hphi, Dphi);
	TransferDataToDevice(M,Hg, Dg);
	segnr=1;
	InitializeForward(1, M, segnr, Dgi, Dg);
cout.precision(15);
cout << "norm1=" << Norm2(M,1,Dgi) << endl;
	for (s=2; s<=N/2; s++) {
		Forward(s,M,segnr,Dgi,Dg,GetForce(),Info);


	}
	s = N/2;

	InitializeBackward(M,s,Dgi_inv,Dgi); //truck; s en segnr worden in .cu file op dezelfde manier behandeld.
	if (N%2 == 1) {
		s = N/2 + 1;
		Backward(M,segnr,Dgi_inv,Dg,Dgx,-GetForce(),Info);
		Composition(1,M,segnr,Dgi_inv,Dgi_inv,Dphi); //ook een truck. Hier vul ik voor eerste argument een 1 in om het programma om de tuin te leiden...
		NormPhi(M,segnr,Dphi,Chalf); // correction needed for optimization
	}
	for (s=(N+3)/2; s<=N; s++) {
		s_inv = N-s+1;
		Backward(M,segnr,Dgi_inv,Dg,Dgx,-GetForce(),Info);
		Composition(s_inv,M,segnr,Dgi_inv,Dgi,Dphi);
	}
cout << "normN=" << Norm2(M,1,Dgi_inv) << endl;
	TransferDataToHost(M,Hgi_inv,Dgi_inv);
	for (int z=1; z<=M; z++) Gi_inv[z]=Hgi_inv[z-1];
	lnGN = Lat->ComputeLnGN(Gi_inv);
	CorrectDoubleCounting(M,1,Dphi,Dg);
	TransferDataToHost(M,Hphi,Dphi);
	phi = Seg->GetPhi(DensPart);
	for (int z=1; z<=M; z++) {
		phi[z]=Hphi[z-1];
	}

	if (freedom == fixedTheta) {
		lnCt = log(theta)-lnGN-log(1.0*N);
		Lat->NormPhiRestr(phi,Gi_inv,2*theta/N); // factor 2 optimization correction
		phiBulk = theta/exp(lnGN);
		SetPhiBulk(phiBulk);
	} else if (freedom == rangeRestricted) {
		double phiPreNorm = 0;
		for (z=1; z<=M; z++) {
			if (restrictedRange->InRange(z)) {
				phiPreNorm += phi[z];
			}
		}
		lnCt = log(phiRange/(2*phiPreNorm));// factor 2 optimization correction
		Lat->NormPhiFree(phi,phiRange/phiPreNorm);
		theta = ComputeTheta();
		phiBulk = exp(lnCt)*N;
		SetPhiBulk(phiBulk);
	} else {
		Lat->NormPhiFree(phi,2*phiBulk/N);// factor 2 optimization correction
		lnCb = log(phiBulk/N);
	}
	Lat->RestoreFromSafe(G);
	Lat->RestoreFromSafe(phi);
	FreeMemoryOnDevice(Dphi);
	FreeMemoryOnDevice(Dg);
	FreeMemoryOnDevice(Dgi);
	FreeMemoryOnDevice(Dgi_inv);
	FreeMemoryOnDevice(Dgx);
	free(Hphi); free(Hg); free(Hgi_inv);
#else
	Message(fatal,MyInput,"Programming error: entered Cuda enabled routine but the compilation was not done accordingly. Compile SFBox sith Cuda=1");
#endif
}


void
SF_HomoRing1stO::Matrix1Long(const DensityPart DensPart) {
	int z,s,s0,t0,v0,t,rs1;
	Matrix Gi(1,M,1,n);
	Vector Gi_inv(1,M), Gs(1,M);
	// in constructor: Seg=Chain->GetSegment(1)
	Vector phi = Seg->GetPhi(DensPart);
	Vector G = Seg->GetSWF();
	Lat->MakeSafe(G);
	// G and phi for the first segment
	for (z=1; z<=M; z++) {
		Gi[z][1] = Gs[z] = G[z];
		phi[z] = 0;
	}
    Lat->MakeSafe(phi);
	t=1;
	v0=t0=s0 = 0;
	for (s=2; s<=N/2; s++) {
		t++;
		if (t>n) {    //(++t>n)
			t0++;
			if (t0 == n) t0 = ++v0;
			t = t0 + 1;
			s0 = s - t0 - 1;
		}
		if (force_set) {
			Lat->PropagateG(Gs,G,GetForce());
		} else
		Lat->PropagateG(Gs,G);
		if ((t == t0+1 && t0 == v0)
		   || (t == t0+1 && ((n-t0)*(n-t0+1) >= N-1-2*(s0+t0)))
		   || (2*(n-t+s) >= N-1)) {
			for (z=1; z<=M; z++) {
				Gi[z][t] = Gs[z];
			}
		}
	}
	for (z=1; z<=M; z++) {
		Gi_inv[z] = Gs[z];
	}
	if (N%2 == 1) {
		s = N/2 + 1;
		if (force_set) {
			Lat->PropagateG(Gi_inv,G,GetForce());
		} else
		Lat->PropagateG(Gi_inv,G);
		Lat->ConnectG(Gi_inv,Gi_inv,phi);
		Lat->NormPhiFree(phi,0.5);
	}
	for (s=(N+3)/2; s<=N; s++) {
		if (force_set) {
			Lat->PropagateG(Gi_inv,G,GetForce());
		} else
		Lat->PropagateG(Gi_inv,G);
		t = N - s + 1 - s0;
		if (t == t0) {
			s0 += - n + t0;
			if (t0 == v0)
				s0 -= ((n - t0)*(n - t0 + 1))/2;
			t0 --;
			if (t0 < v0)
				v0 = t0;
			for (z=1; z<=M; z++) {
				Gs[z] = Gi[z][t];
			}
			for (rs1=s0+t0+2; rs1<=(N-s+1); rs1++) {
				t++;
				if (force_set) {
					Lat->PropagateG(Gs,G,GetForce());
				} else
				Lat->PropagateG(Gs,G);
				if (t == t0+1 || s0+n == N-s+1) {
					for (z=1; z<=M; z++) {
						Gi[z][t] = Gs[z];
					}
				}
				if (t == n && s0+n < N-s+1) {
					t  = ++t0;
					s0 += n - t0;
				}
			}
			t = n;
		}
		Lat->ConnectG(Gi_inv,Gi,t,phi);
	}
	Lat->CorrectDoubleCountG(phi,G);
	lnGN = Lat->ComputeLnGN(Gi_inv);
	if (freedom == fixedTheta) {
		lnCt = log(theta)-lnGN-log(1.0*N);
		Lat->NormPhiRestr(phi,Gi_inv,2*theta/N);
		phiBulk = theta/exp(lnGN);
		SetPhiBulk(phiBulk);
	} else if (freedom == rangeRestricted) {
		double phiPreNorm = 0;
		for (z=1; z<=M; z++) {
			if (restrictedRange->InRange(z)) {
				phiPreNorm += phi[z];
			}
		}
		lnCt = log(phiRange/(2*phiPreNorm));
		Lat->NormPhiFree(phi,2*exp(lnCt));
		theta = ComputeTheta();
		phiBulk = theta/exp(lnGN);
		SetPhiBulk(phiBulk);
	} else {
		Lat->NormPhiFree(phi,2*phiBulk/N);
		lnCb = log(phiBulk/N);
	}
	Lat->RestoreFromSafe(G);
	Lat->RestoreFromSafe(phi);
}

void
SF_HomoRing1stO::CudaMatrix1Long(const DensityPart DensPart) {
#ifdef CUDA
	int z,s,s_inv,s0,t0,v0,t,rs1,k,segnr;
	double Chalf=0.5;
	int Info[11];
	Vector phi,G;
	Vector Gi_inv(1,M);

	double *Hphi, *Hg;
	double *Hgi_inv=&Gi_inv[1];
	double *Dphi, *Dg, *Dgi, *Dgi_inv, *Dgx, *Dgs;

	Lat ->GetLatticeInfo(Info);
	Dphi =    (double *)AllocateMemoryOnDevice(M);
	Dg =      (double *)AllocateMemoryOnDevice(M);
	Dgi =     (double *)AllocateMemoryOnDevice(M*n);
	Dgi_inv = (double *)AllocateMemoryOnDevice(M);
	Dgx =     (double *)AllocateMemoryOnDevice(M);
	Dgs =     (double *)AllocateMemoryOnDevice(M);
	Hphi =    (double *)malloc(sizeof(double)*M);
	Hg =      (double *)malloc(sizeof(double)*M);

	k=0;
	Seg->GetPhi(DensPart);
	G = Seg->GetSWF();
	phi = Seg->GetPhi(DensPart);
	for (z=1; z<=M; z++) {
		phi[z] = 0;
		Hphi[k] = 0;
		Hg[k] =G[z];
		k++;
	}

	TransferDataToDevice(M,Hphi, Dphi);
	TransferDataToDevice(M,Hg, Dg);

	segnr=1; while (Chain->GetSegment(1) != Chain->GetDiffSegment(segnr)) segnr++;
	InitializeForward(1, M, segnr, Dgi, Dg);
	InitializeForward(1, M, segnr, Dgs, Dg);

	t=1;
	v0=t0=s0 = 0;
	for (s=2; s<=N/2; s++) {
		t++;
		if (t>n) {
			t0++;
			if (t0 == n) t0 = ++v0;
			t = t0 + 1;
			s0 = s - t0 - 1;
		}
		Backward(M,segnr,Dgs,Dg,Dgx,GetForce(),Info);

		if ((t == t0+1 && t0 == v0)
		   || (t == t0+1 && ((n-t0)*(n-t0+1) >= N-1-2*(s0+t0)))
		   || (2*(n-t+s) >= N-1)) {
			InitializeForward(t, M, 1, Dgi, Dgs); //the way to store....sorry...
		}
	}

	InitializeBackward(M,1,Dgi_inv,Dgs);
	if (N%2 == 1) {
		s = N/2 + 1;
		Backward(M,segnr,Dgi_inv,Dg,Dgx,-GetForce(),Info);
		Composition(1,M,segnr,Dgi_inv,Dgi_inv,Dphi);
		NormPhi(M,segnr,Dphi,Chalf);
	}
	for (s=(N+3)/2; s<=N; s++) {
		Backward(M,segnr,Dgi_inv,Dg,Dgx,-GetForce(),Info);
		s_inv = N-s+1;
		t = N - s + 1 - s0;
		if (t == t0) {
			s0 += - n + t0;
			if (t0 == v0) s0 -= ((n - t0)*(n - t0 + 1))/2;
			t0 --;
			if (t0 < v0)	v0 = t0;
			InitializeBackward(M,t,Dgs,Dgi);
			for (rs1=s0+t0+2; rs1<=(N-s+1); rs1++) {
				t++;
				Backward(M,segnr,Dgs,Dg,Dgx,-GetForce(),Info);
				if (t == t0+1 || s0+n == N-s+1) {
					InitializeForward(t, M, 1, Dgi, Dgs);
				}
				if (t == n && s0+n < N-s+1) {
					t  = ++t0;
					s0 += n - t0;
				}
			}
			t = n;
		}
		Composition(t,M,segnr,Dgi_inv,Dgi,Dphi);
	}
	TransferDataToHost(M,Hgi_inv,Dgi_inv);
	lnGN = Lat->ComputeLnGN(Gi_inv);

	CorrectDoubleCounting(M,1,Dphi,Dg);
	TransferDataToHost(M,Hphi,Dphi);


	k=0;
	phi = Seg->GetPhi(DensPart);
	for (int z=1; z<=M; z++) {
		phi[z]=Hphi[k];
		k++;
	}

	if (freedom == fixedTheta) {
		lnCt = log(theta)-lnGN-log(1.0*N);
		Lat->NormPhiRestr(phi,Gi_inv,2*theta/N); // factor 2 optimization correction
		phiBulk = theta/exp(lnGN);
		SetPhiBulk(phiBulk);
	} else if (freedom == rangeRestricted) {
		double phiPreNorm = 0;
		for (z=1; z<=M; z++) {
			if (restrictedRange->InRange(z)) {
				phiPreNorm += phi[z];
			}
		}
		lnCt = log(phiRange/(2*phiPreNorm));// factor 2 optimization correction
		Lat->NormPhiFree(phi,2*exp(lnCt));// factor 2 optimization correction
		theta = ComputeTheta();
		phiBulk = theta/exp(lnGN);
		SetPhiBulk(phiBulk);
	} else {
		Lat->NormPhiFree(phi,2*phiBulk/N);// factor 2 optimization correction
		lnCb = log(phiBulk/N);
	}
	Lat->RestoreFromSafe(G);
	Lat->RestoreFromSafe(phi);


	FreeMemoryOnDevice(Dphi);
	FreeMemoryOnDevice(Dg);
	FreeMemoryOnDevice(Dgi);
	FreeMemoryOnDevice(Dgi_inv);
	FreeMemoryOnDevice(Dgx);
	FreeMemoryOnDevice(Dgs);
	free(Hphi); free(Hg);

#else
	Message(fatal,MyInput,"Programming error: entered Cuda enabled routine but the compilation was not done accordingly. Compile SFBox sith Cuda=1");
#endif

}




