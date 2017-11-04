#include "SF_Copol2ndO.h"
#include "tools.h"
#include <iostream>
#ifdef CUDA
#include "Cuda_tools.h"
#endif

SF_Copol2ndO::SF_Copol2ndO(Text name_,
						   SF_SegmentList* SegQ_,
						   Lat2ndO* Lat_,
						   Input* MyInput_)
: SF_Molecule(name_, SegQ_, Lat_, MyInput_)
{
	//if (Chain->GetMoleculeType() != copolymer) {
	//	Message(fatal,MyInput,
	//	"Programming inconsistency..., using copolymer module for molecule '" +	name + "'");
	//}
	Array<Text> sysNames = MyInput->GetNames("sys");
	Text sysName = sysNames[1];

	//MayerSaupe = MyInput->GetBoolean("sys",sysName,"MayerSaupe",false);
	symmetric = Chain->Symmetric();
	if (symmetric) {
		symmetric =false;
		Message(literal,"Symmetry rejected in second order calculations");
	}
	force_set = (!GetForce()==0);

	if (force_set) {
		Message(fatal,MyInput,"Force ensemble not implemented in Copol2ndO");
	}
	if (force_set) symmetric=false;
#ifdef CUDA
	if (GetGPU()) {GPUactivated = true;
		Message(literal, "Still under construction for GPU; Outcome uncertain. In case of trouble contact frans or do not invoke GPU" );
	}else GPUactivated = false;
#else
	GPUactivated = false;
	if (GetGPU()) {
		Message(literal,"Compile with nvcc (option CUDA=1) to activate GPU... Contact Frans Leermakers. Going classical");
	}
#endif

	numDiffSegments = Chain->GetNumDiffSegments();
	Lat = Lat_;
	if (freedom == secondGeneration) {
		CreatePhi(constrained);
		CreatePhi(unconstrained);
	} else {
		CreatePhi(total);
	}
	//if (MayerSaupe) {
	//	CreatePhi(parallel_to_director);
	//}

	if (gradients==3) {size=6; }//CreatePhi(x_dir);CreatePhi(y_dir);CreatePhi(z_dir);}
	if (gradients==2) {size=5; }//CreatePhi(x_dir);CreatePhi(y_dir);CreatePhi(z_dir);}
	if (gradients==1) {size=3; }//CreatePhi(x_dir); CreatePhi(yz_dir);}
	if (size==3) Mphi=2; else Mphi=3;
	PHI = new double[Mphi*numDiffSegments*M];

	if (symmetric) {
		n = int(pow(N*3,1.0/3)+0.5);
		if (2*N < 120) n++;
		if (2*N < 60) n++;
	} else {
		n = int(pow(N*6,1.0/3)+0.5);
		if (N < 120) n++;
		if (N < 60) n++;
	}
	SetPhiBulk(phiBulk);

}
SF_Copol2ndO::~SF_Copol2ndO() {
	if (freedom == secondGeneration) {
		DeletePhi(constrained);
		DeletePhi(unconstrained);
	} else {
		DeletePhi(total);
	}

	//if (gradients ==1) {  DeletePhi(x_dir); DeletePhi(yz_dir); } else{
	//	DeletePhi(x_dir); DeletePhi(y_dir); DeletePhi(z_dir);
	//}
	delete [] PHI;
	//if (MayerSaupe) {
	//	DeletePhi(parallel_to_director);
	//}
}
void
SF_Copol2ndO::ReComputePhi() {
	Message(fatal,"ReComputePhi not implemented in SF_Copol2ndO");
}

Vector
SF_Copol2ndO::GetLong(){
	Message(fatal,"GetLong not implemented in SF_Copol2ndO");
	Vector x;
	return x;
}
Vector
SF_Copol2ndO::GetShort(){
	Message(fatal,"GetShort not implemented in SF_Copol2ndO");
	Vector x;
	return x;
}
Vector
SF_Copol2ndO::GetBondOrientation( const DensityPart DIR) {
	Vector phi_dir;
	phi_dir.Dim(1,M);
	for (int i=1; i<=numDiffSegments; i++) {
		if (DIR==x_dir) for (int z=0; z<M; z++) phi_dir[z+1]+=               PHI[z+      M*Mphi*(i-1)];
		if (DIR==yz_dir)  for (int z=0; z<M; z++) phi_dir[z+1]+=             PHI[z+ M+   M*Mphi*(i-1)];
		if (DIR==y_dir) for (int z=0; z<M; z++) phi_dir[z+1]+=               PHI[z+ M+   M*Mphi*(i-1)];
		if (gradients>1 && DIR==z_dir) for (int z=0; z<M; z++) phi_dir[z+1]+=PHI[z+ 2*M+ M*Mphi*(i-1)];	
	}
	Lat->SetBoundaries(phi_dir);
	return phi_dir; 
};

MoleculeType
SF_Copol2ndO::GetMoleculeType() const {
	return Chain->GetMoleculeType();
}

void
SF_Copol2ndO::ComputePhi() {
	Vector phix,phiy,phiz,phiyz,phiTot;
	if (Lat->GetNamesBulkBoundaries().Upperbound() > 0) {
		CreatePhi(bulk);
		MatrixBulk(bulk);
	}
	int i;
	for (i=1; i<=numDiffSegments; i++) {
		Seg = Chain->GetDiffSegment(i);
		Vector G = Seg->GetSWF();
		Lat->SetBoundaries(G);
	}
	int z;
	if (symmetric) {
		if (freedom != secondGeneration) {
			if (!saveMemory) {
				if (GPUactivated) {CudaSymMatrix1(total);} else SymMatrix1(total);
			}
			else {
				if (GPUactivated) {Message(literal," SaveMemory not on GPU implemented. going classical instead." );}
				SymMatrix1Long(total);
			}
		} else {
			if (!saveMemory) SymMatrix2ndGen(LatRange,constrained,unconstrained);
			else SymMatrix2ndGenLong(LatRange,constrained,unconstrained);
		}
	} else {
		if (freedom != secondGeneration) {
			if (!saveMemory) {
				if (GPUactivated) {CudaMatrix1(total);} else Matrix1(total);
			}
			else {
				if (GPUactivated) {Message(literal," SaveMemory not on GPU implemented. going classical instead." );}
				Matrix1Long(total);
			}
		} else {
			if (GPUactivated) {Message(literal," Second Generation not on GPU implemented. going classical instead." );}
			if (!saveMemory) Matrix2ndGen(LatRange,constrained,unconstrained);
			else Matrix2ndGenLong(LatRange,constrained,unconstrained);
		}
	}
	if (Chain->SomeSegmentsPinned()) SetPhiBulk(0);
	if (Chain->SomeSegmentsGrafted()) SetPhiBulk(0);
	if (Lat->GetNamesBulkBoundaries().Upperbound() > 0) {
		if (freedom == thirdGeneration || freedom == secondGeneration) {
			CopyBulkBoundaries(unconstrained);
		} else {
			CopyBulkBoundaries(total);
		}
		DeletePhi(bulk);
	}
	if (freedom == secondGeneration) {
		Vector phi1Tot = GetPhi(constrained);;
		Vector phi2Tot = GetPhi(unconstrained);
		for (z=1; z<=M; z++) {
			phi1Tot[z] = 0;
			phi2Tot[z] = 0;
		}
		for (i=1; i<=numDiffSegments; i++) {
			Seg = Chain->GetDiffSegment(i);
			Vector phi1 = Seg->GetPhi(constrained);
			Vector phi2 = Seg->GetPhi(unconstrained);
			for (z=1; z<=M; z++) {
				phi1Tot[z] += phi1[z];
				phi2Tot[z] += phi2[z];
			}
		}
		if (muFixed) {
			Lat->SubtractBoundaries(phi1Tot);
			Lat->MultiplyWithLatticeSites(phi1Tot);
			theta = 0;
			for (int z=1; z<=M; z++) {
				theta += phi1Tot[z];
			}
			Lat->DivideByLatticeSites(phi1Tot);
			Lat->RestoreBoundaries(phi1Tot);
		}
	} else {

		Vector phiTot = GetPhi(total);
		Zero(phiTot,M);
		for (i=1; i<=numDiffSegments; i++) {
			Seg = Chain->GetDiffSegment(i);
			Vector phi = Seg->GetPhi(total);
			add2(phiTot,phi,M);
		}
	}
}
void
SF_Copol2ndO::GetLayerAnalysis(Output* Out) {
	//Message(literal,"getLayerAnalysis not implemented in 3d 2nd order");
}
void
SF_Copol2ndO::StateDistribution(Output* Out) const {
	// function to calculate the dissociation distribution function
 			Message(literal, "Unable to generate multistate "
 			"output in 3D 2nd order");
}

void
SF_Copol2ndO::ContactNumberDistr(Output* Out, Boolean partial) const {
	// function to calculate the endpoint distribution function as
	// a function of the number of surface contacts.
     Message(literal, "ContactNumberDistr output not generated in 3D 2nd order");
}


double SF_Copol2ndO::GetPf(int s1,int s2) {
	int NumPfs = SegQ->GetPf_count();
	Vector Pf = SegQ->GetPfVector();
	Array <Text> Pfx = SegQ->GetPfx();
	Array <Text> Pfy = SegQ->GetPfy();
	Text name1 = Chain->GetSegment(s1)->GetName();
	Text name2 = Chain->GetSegment(s2)->GetName();

	for (int i=1; i<=NumPfs; i++){
		if ((Pfx[i]==name1 && Pfy[i] == name2) || (Pfy[i]==name1 && Pfx[i] == name2)) {
		return Pf[i];}
	}
	return 1.0/5.0;
}

void
SF_Copol2ndO::Matrix1(const DensityPart DensPart) {//tested; should work now.
	int z,s,i,k;
	int I=0;
	SF_MolSegment* SEG;
	Vector phi,G;
	Matrix Gi(1,M*size,1,N);
	Vector Gi_inv(1,M*size);
	for (i=1; i<=numDiffSegments; i++) {
		Seg = Chain->GetDiffSegment(i);
		G = Seg->GetSWF();
		Lat->MakeSafe(G);
		phi = Seg->GetPhi(DensPart);
		Zero(phi,M);
		Lat->MakeSafe(phi);
	}
	Zero(PHI,Mphi*numDiffSegments*M);
	G = Chain->GetSegment(1)->GetSWF();
	double *gi;
	double *g;g=&G[1];
	for (k=0; k<size; k++) {
		gi=&Gi[1+k*M][1]; cp(gi,g,M);
	}
	for (s=2; s<=N; s++) {
		G = Chain->GetSegment(s)->GetSWF();
		//cout << "pf s= " << s << ": " << GetPf(s-1,s)<< endl;
		Lat->PropagateF(Gi,G,s,GetPf(s-1,s)); //this one...
	}
	Seg = Chain->GetSegment(N);
	G = Seg->GetSWF();
	double *gi_inv;
	g=&G[1];
	for (k=0; k<size; k++) {gi_inv=&Gi_inv[1+k*M]; cp(gi_inv,g,M); }

	phi = Seg->GetPhi(DensPart);
	for (i=1; i<=numDiffSegments; i++) {
		SEG = Chain->GetDiffSegment(i);
		if (SEG->GetName()==Seg->GetName()) I=i;

	}
	Lat->ConnectG(Gi_inv,Gi,N,phi,PHI+(I-1)*Mphi*M);


	for (s=N-1; s>=1; s--) {
		Seg = Chain->GetSegment(s);
		G = Seg->GetSWF();
		phi = Seg->GetPhi(DensPart);
		//cout << "pf s= " << s << ": " << GetPf(s,s+1)<< endl;
		Lat->PropagateB(Gi_inv,G,GetPf(s,s+1));
		for (i=1; i<=numDiffSegments; i++) {
			SEG = Chain->GetDiffSegment(i);
			if (SEG->GetName()==Seg->GetName()) I=i;
		}
		Lat->ConnectG(Gi_inv,Gi,s,phi,PHI+(I-1)*Mphi*M);

	}
	lnGN = Lat->ComputeLnGN(Gi_inv);
	for (i=1; i<=numDiffSegments; i++) {
		Seg = Chain->GetDiffSegment(i);
		G = Seg->GetSWF();
		phi = Seg->GetPhi(DensPart);
		Lat->CorrectDoubleCountG(phi,PHI+(i-1)*Mphi*M,G);

	}
	if (freedom == fixedTheta) {
		lnCt = log(theta) - lnGN - log(1.0*N);
		phiBulk = theta/exp(lnGN);
		SetPhiBulk(phiBulk);
	} else if (freedom == rangeRestricted) {
		double phiPreNorm = 0;
		for (i=1; i<=numDiffSegments; i++) {
			Seg = Chain->GetDiffSegment(i);
			phi = Seg->GetPhi(DensPart);
			for (z=1; z<=M; z++) {
				if (restrictedRange->InRange(z)) {
					phiPreNorm += phi[z];
				}
			}
		}
		lnCt = log(phiRange/phiPreNorm);
	} else {
		lnCb = log(phiBulk/N);
	}
	for (i=1; i<=numDiffSegments; i++) {
		Seg = Chain->GetDiffSegment(i);
		phi = Seg->GetPhi(DensPart);
		if (freedom == fixedTheta) {
			Lat->NormPhiRestr(phi,PHI+(i-1)*Mphi*M,Gi_inv,theta/N);

		} else if (freedom == rangeRestricted) {

			//Lat->NormPhiFree(phi,exp(lnCt));
			Lat->NormPhiFree(phi,PHI+(i-1)*Mphi*M,exp(lnCt));
		} else {
			Lat->NormPhiFree(phi,PHI+(i-1)*Mphi*M,phiBulk/(6*N));

		}
		G = Seg->GetSWF();
		Lat->RestoreFromSafe(G);
		Lat->RestoreFromSafe(phi);
	}
	if (freedom == rangeRestricted) {
		theta = ComputeTheta();
		phiBulk = theta/exp(lnGN);
		SetPhiBulk(phiBulk);
	}

}

void
SF_Copol2ndO::CudaMatrix1(const DensityPart DensPart) {//not tested.
#ifdef CUDA
	int Info[11];
	int z,s,i,k,segnr;
	double S;
	Vector phi,G;
	Vector Gi_inv(1,M*size);

	double *Hphi, *Hg;				//pointers to host memory;
	double *Hgi_inv = &Gi_inv[1];
	double *Dphi, *Dg, *Dgi, *Dgi_inv, *Dgx; //pointers to device memory;

	Lat ->GetLatticeInfo(Info);

	Dphi =    (double *)AllocateMemoryOnDevice(M*numDiffSegments);
	Dg =      (double *)AllocateMemoryOnDevice(M*numDiffSegments);
	Dgi =     (double *)AllocateMemoryOnDevice(M*N*size);
	Dgi_inv = (double *)AllocateMemoryOnDevice(M*size);
	Dgx =     (double *)AllocateMemoryOnDevice(M*size);
	Hphi =    (double *)malloc(sizeof(double)*M*numDiffSegments);
	Hg =      (double *)malloc(sizeof(double)*M*numDiffSegments);

	k=0;
	for (i=1; i<=numDiffSegments; i++) {
		Seg = Chain->GetDiffSegment(i);
		G = Seg->GetSWF();
		phi = Seg->GetPhi(DensPart);
		for (z=1; z<=M; z++) {
			phi[z] = 0;
			Hphi[k] = 0;
			Hg[k] =G[z];
			k++;
		}
	}

	TransferDataToDevice(M*numDiffSegments,Hphi, Dphi);
	TransferDataToDevice(M*numDiffSegments,Hg, Dg);

	segnr=1; while (Chain->GetSegment(1) != Chain->GetDiffSegment(segnr)) segnr++;
	InitializeForward(1,M,segnr,6,Dgi,Dg);

	//G = Chain->GetSegment(1)->GetSWF();
	//for (k=0; k<6; k++) {double *gi=&Gi[1+k*M][1], *g=&G[1]; cp(gi,g,M); }

	for (s=2; s<=N; s++) {
		segnr=1; while (Chain->GetSegment(s) != Chain->GetDiffSegment(segnr)) segnr++;
		S=GetPf(s-1,s);
		Forward(s,M,segnr,Dgi,Dg,S,GetForce(),Info);

		//G = Chain->GetSegment(s)->GetSWF();
		//Lat->PropagateF(Gi,G,s,GetPf(s-1,s));
	}

	InitializeBackward(M, segnr, 6, Dgi_inv, Dg);

	//Seg = Chain->GetSegment(N);
	//G = Seg->GetSWF();
	//for (k=0; k<6; k++) {double *gi_inv=&Gi_inv[1+k*M], *g=&G[1]; cp(gi_inv,g,M); }

	//phi = Seg->GetPhi(DensPart);
	//Lat->ConnectG(Gi_inv,Gi,N,phi);
	Composition(N,M,segnr,6,Dgi_inv,Dgi,Dphi);

	for (s=N-1; s>=1; s--) {
		//Seg = Chain->GetSegment(s);
		//G = Seg->GetSWF();
		//phi = Seg->GetPhi(DensPart);

		//Lat->PropagateB(Gi_inv,G,GetPf(s,s+1));
		//Lat->ConnectG(Gi_inv,Gi,s,phi);
		segnr=1; while (Chain->GetSegment(s) != Chain->GetDiffSegment(segnr)) segnr++;
		S=GetPf(s,s+1);  //Needs to be checked.
		Backward(M,segnr,Dgi_inv,Dg,Dgx,S,-GetForce(),Info);
		Composition(s,M,segnr,6,Dgi_inv,Dgi,Dphi);

	}
    TransferDataToHost(size*M,Hgi_inv,Dgi_inv);
	lnGN = Lat->ComputeLnGN(Gi_inv);

	CorrectDoubleCounting(M*numDiffSegments,1,Dphi,Dg);
	TransferDataToHost(M*numDiffSegments,Hphi,Dphi);

	k=0;
	for (i=1; i<=numDiffSegments; i++) {
		Seg = Chain->GetDiffSegment(i);
		phi = Seg->GetPhi(DensPart);
		for (int z=1; z<=M; z++) {
			phi[z]=Hphi[k];
			k++;
		}
	}

	if (freedom == fixedTheta) {
		lnCt = log(theta) - lnGN - log(1.0*N);
		phiBulk = theta/exp(lnGN);
		SetPhiBulk(phiBulk);
	} else if (freedom == rangeRestricted) {
		double phiPreNorm = 0;
		for (i=1; i<=numDiffSegments; i++) {
			Seg = Chain->GetDiffSegment(i);
			phi = Seg->GetPhi(DensPart);
			for (z=1; z<=M; z++) {
				if (restrictedRange->InRange(z)) {
					phiPreNorm += phi[z];
				}
			}
		}
		lnCt = log(phiRange/phiPreNorm);
	} else {
		lnCb = log(phiBulk/N);
	}
	for (i=1; i<=numDiffSegments; i++) {
		Seg = Chain->GetDiffSegment(i);
		phi = Seg->GetPhi(DensPart);
		if (freedom == fixedTheta) {
			Lat->NormPhiRestr(phi,Gi_inv,theta/N);
		} else if (freedom == rangeRestricted) {
			Lat->NormPhiFree(phi,exp(lnCt));
		} else {
			Lat->NormPhiFree(phi,phiBulk/(6*N));
		}
		G = Seg->GetSWF();
		Lat->RestoreFromSafe(G);
		Lat->RestoreFromSafe(phi);
	}
	if (freedom == rangeRestricted) {
		theta = ComputeTheta();
		phiBulk = theta/exp(lnGN);
		SetPhiBulk(phiBulk);
	}
	FreeMemoryOnDevice(Dphi);
	FreeMemoryOnDevice(Dg);
	FreeMemoryOnDevice(Dgi);
	FreeMemoryOnDevice(Dgi_inv);
	FreeMemoryOnDevice(Dgx);
	free(Hphi); free(Hg);
#else
	Message(fatal,MyInput,"Programming error: entered Cuda enabled routine but the compilation was not done accordingly. Compile SFBox sith Cuda=1");
#endif
}


void
SF_Copol2ndO::Matrix1Long(const DensityPart DensPart) {//tested should work.
	int z,s,s0,t0,v0,t,rs1,i,k,I=0;
	SF_MolSegment* SEG;
	Vector phi,G;
	Matrix Gi(1,M*size,1,n);
	Vector Gi_inv(1,M*size), Gs(1,M*size);

	for (i=1; i<=numDiffSegments; i++) {
		Seg = Chain->GetDiffSegment(i);
		G = Seg->GetSWF();
		Lat->MakeSafe(G);
		phi = Seg->GetPhi(DensPart);
		Zero(phi,M);
		Lat->MakeSafe(phi);
	}
	Zero(PHI,Mphi*numDiffSegments*M);
	G = Chain->GetSegment(1)->GetSWF();
	for (k=0; k<size; k++) {double *gi=&Gi[1+k*M][1], *g=&G[1], *gs=&Gs[1+k*M]; cp(gi,g,M); cp(gs,g,M); }
	t=1;
	v0=t0=s0 = 0;
	for (s=2; s<=N; s++) {
		t++;
		if (t>n) {
			t0++;
			if (t0 == n) t0 = ++v0;
			t = t0 + 1;
			s0 = s - t0 - 1;
		}
		G = Chain->GetSegment(s)->GetSWF();
		Lat->PropagateF(Gs,G,GetPf(s-1,s));
		if ((t == t0+1 && t0 == v0)
		   || (t == t0+1 && ((n-t0)*(n-t0+1) >= N-1-2*(s0+t0)))
		   || (2*(n-t+s) >= N-1))
		  {double *gi=&Gi[1][t], *gs=&Gs[1]; cp(gi,gs,size*M);}
	}
	for (s=N; s>=1; s--) {
		G = Chain->GetSegment(s)->GetSWF();
		if (s == N) {
			for (k=0; k<size; k++) {double *gi_inv=&Gi_inv[1+k*M],*g=&G[1]; cp(gi_inv,g,M);}
		} else {
			Lat->PropagateB(Gi_inv,G,GetPf(s,s+1));
		}
		t = s - s0;
		if (t == t0) {
			s0 += - n + t0;
			if (t0 == v0 ) {
				s0 -= ((n - t0)*(n - t0 + 1))/2;
			}
			t0 --;
			if (t0 < v0) {
				v0 = t0;
			}
			double *gi=&Gi[1][t], *gs=&Gs[1]; cp(gs,gi,size*M);

			for (rs1=s0+t0+2; rs1<=s; rs1++) {
				t++;
				G = Chain->GetSegment(rs1)->GetSWF();
				Lat->PropagateF(Gs,G,GetPf(s-1,s));
				if (t == t0+1 || s0+n == s) {
					double *gi=&Gi[1][t], *gs=&Gs[1]; cp(gi,gs,size*M);
				}
				if (t == n && s0+n < s) {
					t  = ++t0;
					s0 += n - t0;
				}
			}
			t = n;
		}
		Seg = Chain->GetSegment(s);
		phi = Seg->GetPhi(DensPart);
		for (i=1; i<=numDiffSegments; i++) {
			SEG = Chain->GetDiffSegment(i);
			if (SEG==Seg) I=i;
		}
		Lat->ConnectG(Gi_inv,Gi,t,phi,PHI+(I-1)*Mphi*M);
	}
	lnGN = Lat->ComputeLnGN(Gi_inv);
	for (i=1; i<=numDiffSegments; i++) {
		Seg = Chain->GetDiffSegment(i);
		G = Seg->GetSWF();
		phi = Seg->GetPhi(DensPart);
		Lat->CorrectDoubleCountG(phi,PHI+(i-1)*Mphi*M,G);
	}
	if (freedom == fixedTheta) {
		lnCt = log(theta) - lnGN - log(1.*N);
		phiBulk = theta/exp(lnGN);
		SetPhiBulk(phiBulk);
	} else if (freedom == rangeRestricted) {
		double phiPreNorm = 0;
		for (i=1; i<=numDiffSegments; i++) {
			Seg = Chain->GetDiffSegment(i);
			phi = Seg->GetPhi(DensPart);
			for (z=1; z<=M; z++) {
				if (restrictedRange->InRange(z)) {
					phiPreNorm += phi[z];
				}
			}
		}
		lnCt = log(phiRange/phiPreNorm);
	} else {
		lnCb = log(phiBulk/N);
	}
	for (i=1; i<=numDiffSegments; i++) {
		Seg = Chain->GetDiffSegment(i);
		phi = Seg->GetPhi(DensPart);
		if (freedom == fixedTheta) {
			Lat->NormPhiRestr(phi,PHI+(I-1)*Mphi*M,Gi_inv,theta/N);
		} else if (freedom == rangeRestricted) {
			Lat->NormPhiFree(phi,PHI+(I-1)*Mphi*M,exp(lnCt));
		} else {
			Lat->NormPhiFree(phi,PHI+(I-1)*Mphi*M,phiBulk/(6*N));
		}
		G = Seg->GetSWF();
		Lat->RestoreFromSafe(G);
		Lat->RestoreFromSafe(phi);
	}
	if (freedom == rangeRestricted) {
		theta = ComputeTheta();
		phiBulk = theta/exp(lnGN);
		SetPhiBulk(phiBulk);
	}
}
void SF_Copol2ndO::Matrix2ndGen(const LatticeRange* Range,
								const DensityPart DensPart1,
								const DensityPart DensPart2) {
	Message(literal, "Matrix2nd Gen not implemented in 3D 2nd order");
}
void
SF_Copol2ndO::Matrix2ndGenLong(const LatticeRange* Range,
							   const DensityPart DensPart1,
							   const DensityPart DensPart2) {
	Message(literal, "Matrix2ndGenLong not implemented in 3D 2nd order");
}
void
SF_Copol2ndO::SymMatrix1(const DensityPart DensPart) { //tested: should work.
 Message(fatal, "SymMatrix1 can not work in 2nd order");
}

void
SF_Copol2ndO::CudaSymMatrix1(const DensityPart DensPart) {//not tested possibly error by halfway point
#ifdef CUDA
	int Info[12];
	int z,s,i,k,segnr;
	double S;
	Message(fatal, "CudaSymMatrix1 can not work in 2nd order");
	Vector phi,G;
	//Matrix Gi(1,M*size,1,N/2);
	Vector Gi_inv(1,M*size);
	//Vector Gc(1,M*size); //


	double *Hphi, *Hg;
	double *Hgi_inv=&Gi_inv[1];  //pointers to host memory;
	double *Dphi, *Dg, *Dgi, *Dgi_inv, *Dgx, *Dgc; //pointers to device memory;
	Lat ->GetLatticeInfo(Info);

	Dphi =    (double *)AllocateMemoryOnDevice(M*numDiffSegments);
	Dg =      (double *)AllocateMemoryOnDevice(M*numDiffSegments);
	Dgi =     (double *)AllocateMemoryOnDevice(M*N*3);
	Dgi_inv = (double *)AllocateMemoryOnDevice(M*size);
	Dgc =     (double *)AllocateMemoryOnDevice(M*size);
	Dgx =     (double *)AllocateMemoryOnDevice(M*size);
	Hphi =    (double *)malloc(sizeof(double)*M*numDiffSegments);
	Hg =      (double *)malloc(sizeof(double)*M*numDiffSegments);

	k=0;
	for (i=1; i<=numDiffSegments; i++) {
		Seg = Chain->GetDiffSegment(i);
		G = Seg->GetSWF();
		phi = Seg->GetPhi(DensPart);
		for (z=1; z<=M; z++) {
			phi[z] = 0;
			Hphi[k] = 0;
			Hg[k] =G[z];
			k++;
		}
	}

	TransferDataToDevice(M*numDiffSegments,Hphi, Dphi);
	TransferDataToDevice(M*numDiffSegments,Hg, Dg);
	segnr=1; while (Chain->GetSegment(1) != Chain->GetDiffSegment(segnr)) segnr++;
	InitializeForward(1,M,segnr,6,Dgi,Dg);

	//G = Chain->GetSegment(1)->GetSWF();
	//for (k=0; k<6; k++) {double *gi=&Gi[1+k*M][1], *g=&G[1]; cp(gi,g,M); }

	for (s=2; s<=N/2; s++) {
		segnr=1; while (Chain->GetSegment(s) != Chain->GetDiffSegment(segnr)) segnr++;
		if (s==N) S=-1; else S=GetPf(s-1,s);//Negative S signals end of chain.
		Forward(s,M,segnr,Dgi,Dg,S,GetForce(),Info);

		//G = Chain->GetSegment(s)->GetSWF();
		//Lat->PropagateF(Gi,G,s,GetPf(s-1,s));
	}
	s = N/2;
	segnr=1; while (Chain->GetSegment(s) != Chain->GetDiffSegment(segnr)) segnr++;
	InitializeBackward(size*M, s, Dgi_inv, Dgi);
	InitializeBackward(size*M, s, Dgc, Dgi);
	//double *gi_inv=&Gi_inv[1], *gi=&Gi[1][s], *gc=&Gc[1]; cp(gi_inv,gi,size*M); cp(gc,gi,size*M);

	if (N%2 == 1) {
		s = N/2 + 1;
		segnr=1; while (Chain->GetSegment(s) != Chain->GetDiffSegment(segnr)) segnr++;
		Backward(M,segnr,Dgc,Dg,Dgx,GetPf(s-1,s),GetForce(),Info);
		Backward(M,segnr,Dgi_inv,Dg,Dgx,GetPf(s,s+1),GetForce(),Info);
		Composition(1,M,segnr,6,Dgi_inv,Dgc,Dphi);
		NormPhi(M,segnr,Dphi,0.5);
		//Seg = Chain->GetSegment(s);
		//G = Seg->GetSWF();
		//Lat->PropagateF(Gc,G,GetPf(s-1,s));
		//phi = Seg->GetPhi(DensPart);
		//Lat->PropagateB(Gi_inv,G,GetPf(s,s+1));
		//Lat->ConnectG(Gi_inv,Gc,phi);
		//Lat->ConnectG(Gi_inv,Gi_inv,phi);
		//Lat->NormPhiFree(phi,0.5); // correction needed for optimization
	}
	for (s=(N-1)/2; s>=1; s--) {
		segnr=1; while (Chain->GetSegment(s) != Chain->GetDiffSegment(segnr)) segnr++;
		Backward(M,segnr,Dgi_inv,Dg,Dgx,GetPf(s,s+1),GetForce(),Info);
		Composition(s,M,segnr,6,Dgi_inv,Dgi,Dphi);
		//Seg = Chain->GetSegment(s);
		//G = Seg->GetSWF();
		//Lat->PropagateB(Gi_inv,G,GetPf(s,s+1));
		//phi = Seg->GetPhi(total);
		//Lat->ConnectG(Gi_inv,Gi,s_inv,phi); // optimization, chain is symmetric
			// Connect only half the chain into phi
	}
       TransferDataToHost(size*M,Hgi_inv,Dgi_inv);
	lnGN = Lat->ComputeLnGN(Gi_inv);

	CorrectDoubleCounting(M*numDiffSegments,1,Dphi,Dg);
	TransferDataToHost(M*numDiffSegments,Hphi,Dphi);

	k=0;
	for (i=1; i<=numDiffSegments; i++) {
		Seg = Chain->GetDiffSegment(i);
		phi = Seg->GetPhi(DensPart);
		for (int z=1; z<=M; z++) {
			phi[z]=Hphi[k];
			k++;
		}
	}

	if (freedom == fixedTheta) {
		lnCt = log(theta) - lnGN -log(1.0*N);
		phiBulk = theta/exp(lnGN);
		Chain->SetPhiBulk(phiBulk);
	} else if (freedom == rangeRestricted) {
		double phiPreNorm = 0;
		for (i=1; i<=numDiffSegments; i++) {
			Seg = Chain->GetDiffSegment(i);
			phi = Seg->GetPhi(DensPart);
			for (z=1; z<=M; z++) {
				if (restrictedRange->InRange(z)) {
					phiPreNorm += phi[z];
				}
			}
		}
		lnCt = log(phiRange/(2*phiPreNorm)); // factor 2 is correction needed for optimization
	} else {
		lnCb = log(phiBulk/N);
	}
	for (i=1; i<=numDiffSegments; i++) {
		Seg = Chain->GetDiffSegment(i);
		phi = Seg->GetPhi(DensPart);
		if (freedom == fixedTheta) {
			Lat->NormPhiRestr(phi,Gi_inv,2*theta/N); // factor 2 is correction needed for optimization
		} else if (freedom == rangeRestricted) {
			Lat->NormPhiFree(phi,2*exp(lnCt));// factor 2 is correction needed for optimization
		} else {
			Lat->NormPhiFree(phi,2*phiBulk/(6*N));// factor 2 is correction needed for optimization
		}
		G = Seg->GetSWF();
		Lat->RestoreFromSafe(G);
		Lat->RestoreFromSafe(phi);
	}
	if (freedom == rangeRestricted) {
		theta = ComputeTheta();
		phiBulk = theta/exp(lnGN);
		SetPhiBulk(phiBulk);
	}
	FreeMemoryOnDevice(Dphi);
	FreeMemoryOnDevice(Dg);
	FreeMemoryOnDevice(Dgi);
	FreeMemoryOnDevice(Dgi_inv);
	FreeMemoryOnDevice(Dgx);
	FreeMemoryOnDevice(Dgc);
	free(Hphi); free(Hg);
#else
	Message(fatal,MyInput,"Programming error: entered Cuda enabled routine but the compilation was not done accordingly. Compile SFBox sith Cuda=1");
#endif
}

void
SF_Copol2ndO::SymMatrix1Long(const DensityPart DensPart) {//should work.
	int z,s;
	//int s_inv;
	int s0,t0,v0,t,rs1,i,k;
Message(fatal, "SymMatrix1Long can not work in 2nd order");
	Vector phi,G;
	Matrix Gi(1,M*size,1,n);
	Vector Gi_inv(1,M*size), Gs(1,M*size);
	Vector Gc(1,M*size);
	for (i=1; i<=numDiffSegments; i++) {
		Seg = Chain->GetDiffSegment(i);
		G = Seg->GetSWF();
		Lat->MakeSafe(G);
		phi = Seg->GetPhi(DensPart);
		Zero(phi,M);
		Lat->MakeSafe(phi);
	}
	G = Chain->GetSegment(1)->GetSWF();
	for (k=0; k<size; k++){ double *gi=&Gi[1+k*M][1], *gs=&Gs[1+k*M], *g=&G[1]; cp(gi,g,M); cp(gs,g,M);}

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
		G = Chain->GetSegment(s)->GetSWF();
		Lat->PropagateF(Gs,G,GetPf(s-1,s));
		if ((t == t0+1 && t0 == v0)
		   || (t == t0+1 && ((n-t0)*(n-t0+1) >= N-1-2*(s0+t0)))
		   || (2*(n-t+s) >= N-1)) {
		    double *gi=&Gi[1][t], *gs=&Gs[1]; cp(gi,gs,size*M);
		}
	}
	double *gi_inv=&Gi_inv[1], *gs=&Gs[1], *gc=&Gc[1]; cp(gi_inv,gs,size*M); cp(gc,gs,size*M);
	if (N%2 == 1) {
		s = N/2 + 1;
		Seg = Chain->GetSegment(s);
		G = Seg->GetSWF();
		Lat->PropagateF(Gc,G,GetPf(s-1,s));
		phi = Seg->GetPhi(DensPart);
		Lat->PropagateB(Gi_inv,G,GetPf(s,s+1));
		Lat->ConnectG(Gi_inv,Gc,phi);
		Lat->NormPhiFree(phi,0.5);
	}
	for (s=(N+3)/2; s<=N; s++) {
		G = Chain->GetSegment(s)->GetSWF();
		Lat->PropagateB(Gi_inv,G,GetPf(s,s+1));
		//s_inv = N-s+1;
		t = N - s + 1 - s0;
		if (t == t0) {
			s0 += - n + t0;
			if (t0 == v0)
				s0 -= ((n - t0)*(n - t0 + 1))/2;
			t0 --;
			if (t0 < v0)
				v0 = t0;
			double *gs=&Gs[1], *gi=&Gi[1][t]; cp(gs,gi,size*M);

			for (rs1=s0+t0+2; rs1<=(N-s+1); rs1++) {
				t++;
				G = Chain->GetSegment(rs1)->GetSWF();
				Lat->PropagateF(Gs,G,GetPf(s-1,s));
				if (t == t0+1 || s0+n == N-s+1) {
					double *gs=&Gs[1], *gi=&Gi[1][t]; cp(gi,gs,size*M);
				}
				if (t == n && s0+n < N-s+1) {
					t  = ++t0;
					s0 += n - t0;
				}
			}
			t = n;
		}
		Seg = Chain->GetSegment(s);
		phi = Seg->GetPhi(DensPart);
		Lat->ConnectG(Gi_inv,Gi,t,phi);

	}
	lnGN = Lat->ComputeLnGN(Gi_inv);
	for (i=1; i<=numDiffSegments; i++) {
		Seg = Chain->GetDiffSegment(i);
		G = Seg->GetSWF();
		phi = Seg->GetPhi(DensPart);
		Lat->CorrectDoubleCountG(phi,G);
	}
	if (freedom == fixedTheta) {
		lnCt = log(theta) - lnGN -log(1.*N);
		phiBulk = theta/exp(lnGN);
		Chain->SetPhiBulk(phiBulk);
	} else if (freedom == rangeRestricted) {
		double phiPreNorm = 0;
		for (i=1; i<=numDiffSegments; i++) {
			Seg = Chain->GetDiffSegment(i);
			phi = Seg->GetPhi(DensPart);
			for (z=1; z<=M; z++) {
				if (restrictedRange->InRange(z)) {
					phiPreNorm += phi[z];
				}
			}
		}
		lnCt = log(phiRange/(2*phiPreNorm)); // factor 2 is correction needed for optimization
	} else {
		lnCb = log(phiBulk/N);
	}
	for (i=1; i<=numDiffSegments; i++) {
		Seg = Chain->GetDiffSegment(i);
		phi = Seg->GetPhi(DensPart);
		if (freedom == fixedTheta) {
			Lat->NormPhiRestr(phi,Gi_inv,2*theta/N); // factor 2 is correction needed for optimization
		} else if (freedom == rangeRestricted) {
			Lat->NormPhiFree(phi,2*exp(lnCt));// factor 2 is correction needed for optimization
		} else {
			Lat->NormPhiFree(phi,2*phiBulk/(6*N));// factor 2 is correction needed for optimization
		}
		G = Seg->GetSWF();
		Lat->RestoreFromSafe(G);
		Lat->RestoreFromSafe(phi);
	}
	if (freedom == rangeRestricted) {
		theta = ComputeTheta();
		phiBulk = theta/exp(lnGN);
		SetPhiBulk(phiBulk);
	}
}
void SF_Copol2ndO::SymMatrix2ndGen(const LatticeRange* Range,
								   const DensityPart DensPart1,
								   const DensityPart DensPart2) {
	Message(literal,"SymMatrix2ndGen not implemented in 3d 2nd order");
}
void
SF_Copol2ndO::SymMatrix2ndGenLong(const LatticeRange* Range,
								  const DensityPart DensPart1,
								  const DensityPart DensPart2) {
	Message(literal,"SymMatrix2ndGenLong not implemented in 3d 2nd order");
}
void
SF_Copol2ndO::MatrixBulk(const DensityPart DensPart) {
	int i,z;
	Vector phiTot = GetPhi(DensPart);
	for (z=1; z<=M; z++) {
		if (Lat->BulkBoundary(z)) {
			phiTot[z] = 0;
		}
	}
	Seg = Chain->GetDiffSegment(1);
	Vector G = Seg->GetSWF();
	for (z=1; z<=M; z++) {
		if (Lat->BulkBoundary(z)) {
			phiTot[z] = pow(G[z],Chain->GetAvNumSegments(Seg));
		}
	}
	for (i=2; i<=numDiffSegments; i++) {
		Seg = Chain->GetDiffSegment(i);
		Vector G = Seg->GetSWF();
		for (z=1; z<=M; z++) {
			if (Lat->BulkBoundary(z)) {
				phiTot[z] *= pow(G[z],Chain->GetAvNumSegments(Seg));
			}
		}
	}
	for (z=1; z<=M; z++) {
		if (Lat->BulkBoundary(z)) {
			phiTot[z] *= phiBulk;
		}
	}
	for (i=1; i<=numDiffSegments; i++) {
		Seg = Chain->GetDiffSegment(i);
		Vector phi = Seg->GetPhi(DensPart);
		for (z=1; z<=M; z++) {
			if (Lat->BulkBoundary(z)) {
				phi[z] = phiTot[z]*Chain->GetAvNumSegments(Seg)/Chain->GetAvLength();
			}
		}
	}
}
void
SF_Copol2ndO::CopyBulkBoundaries(const DensityPart DensPart) {
	for (int i=1; i<=numDiffSegments; i++) {
		Seg = Chain->GetDiffSegment(i);
		Vector phi = Seg->GetPhi(DensPart);
		Vector PhiBulk = Seg->GetPhi(bulk);
		for (int z=1; z<=M; z++) {
			if (Lat->BulkBoundary(z)) {
				phi[z] = PhiBulk[z];
			}
		}
	}
}

