#include "SF_Homopol2ndO.h"
#include "tools.h"

SF_Homopol2ndO::SF_Homopol2ndO(Text name_,
								   SF_SegmentList* SegQ_,
								   Lat2ndO* Lat_,
								   Input* MyInput_)
: SF_Molecule(name_, SegQ_, Lat_, MyInput_)
{
	Lat = Lat_;
	Array<Text> sysNames = MyInput->GetNames("sys");
	Text sysName = sysNames[1];
	//MayerSaupe = MyInput->GetBoolean("sys",sysName,"MayerSaupe",false);
	force_set = (!GetForce()==0);
	if (force_set) {
		Message(fatal,MyInput,"Force ensemble not implemented in Homopol2stO");
	}
	if (GetGPU()) {
		Message(literal,MyInput,"For homopolymers 2nd order the GPU card is not yet activated. Contact Frans Leermakers. Going classical instead");
	}
	//SegQ = SegQ_;
	Seg = Chain->GetSegment(1);
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
	//if (MayerSaupe) {
	//	CreatePhi(parallel_to_director);
	//}

	n = int(pow(N*3,1.0/3)+0.5);

	if (gradients==3) {size=6;}
	if (gradients==2) {size=5;}
	if (gradients==1) {size=3;}
	if (size==3) Mphi=2; else Mphi=3;
	PHI = new double[Mphi*M];

	if (2*N < 60) n++;
	if (2*N < 30) n++;
	SetPhiBulk(phiBulk);


}
SF_Homopol2ndO::~SF_Homopol2ndO() {
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
	delete [] PHI;
}
void
SF_Homopol2ndO::ReComputePhi() {
	Message(fatal,"ReComputePhi not implemented in SF_Homopol2ndO");
}

Vector
SF_Homopol2ndO::GetLong(){
	Message(fatal,"GetLong not implemented in SF_Homopol2ndO");
	Vector x;
	return x;
}
Vector
SF_Homopol2ndO::GetShort(){
	Message(fatal,"GetShort not implemented in SF_Homopol2ndO");
	Vector x;
	return x;
}
Vector
SF_Homopol2ndO::GetBondOrientation( const DensityPart DIR) {
	Vector phi_dir;
	phi_dir.Dim(1,M);
	if (DIR==x_dir) for (int z=0; z<M; z++)  phi_dir[z+1]=              PHI[z    ];
	if (DIR==yz_dir) for (int z=0; z<M; z++) phi_dir[z+1]=              PHI[z+M  ];
	if (DIR==y_dir) for (int z=0; z<M; z++)  phi_dir[z+1]=              PHI[z+M  ];
	if (gradients>1 && DIR==z_dir) for (int z=0; z<M; z++) phi_dir[z+1]=PHI[z+2*M];	
	Lat->SetBoundaries(phi_dir);
	return phi_dir; 
}; 

MoleculeType
SF_Homopol2ndO::GetMoleculeType() const {
	return Chain->GetMoleculeType();
}


void
SF_Homopol2ndO::ComputePhi() {
	if (Lat->GetNamesBulkBoundaries().Upperbound() > 0) {
		CreatePhi(bulk);
		MatrixBulk(bulk);
	}
	Vector G = Seg->GetSWF();
	Lat->SetBoundaries(G);
	if (freedom == thirdGeneration) {
		if (!saveMemory) Matrix3rdGen(LatRange,constrained,unconstrained,renorm,
			LatRangeRenorm);
		else  Matrix3rdGenLong(LatRange,constrained,unconstrained,renorm,
			LatRangeRenorm);
	} else if (freedom == secondGeneration) {
		if (!saveMemory) Matrix2ndGen(LatRange,constrained,unconstrained);
		else Matrix2ndGenLong(LatRange,constrained,unconstrained);
	} else {
		if (!saveMemory) Matrix1(total);
		else Matrix1Long(total);
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

	Vector Phi;
}

void
SF_Homopol2ndO::GetLayerAnalysis(Output* Out) {
	//Message(debug,"GetLayerAnalysis not implemented in second order"); return;
}

double SF_Homopol2ndO::GetPf() {
	int NumPfs = SegQ->GetPf_count();
	Vector Pf = SegQ->GetPfVector();
	Array<Text> Pfx= SegQ->GetPfx();
	Array<Text> Pfy= SegQ->GetPfy();
	Text name = Seg->GetName();
	for (int i=1; i<=NumPfs; i++){
		if (Pfx[i]==name && Pfy[i] == name ) return Pf[i];
	}
	return 1.0/5.0;
}


void
SF_Homopol2ndO::Matrix1(const DensityPart DensPart) {
	int z,k,s;
	double S=GetPf();
	Matrix Gi(1,M*size,1,N);
	//double *gi; for (s=1; s<=N/2; s++) {*gi=&Gi[1][s]; Zero(gi,6*M);}
	Vector Gi_inv(1,M*size);
	Vector phi = Seg->GetPhi(DensPart);
	Vector G = Seg->GetSWF();

	Lat->MakeSafe(G);
	Zero(phi,M);
	Zero(PHI,Mphi*M);


	for (k=0; k<size; k++) {double *gi = &Gi[1+k*M][1], *g=&G[1]; cp(gi,g,M);	}

	Lat->MakeSafe(phi);
	for (s=2; s<=N; s++) {Lat->PropagateF(Gi,G,s,S);}


	for (k=0; k<size; k++) {double *gi_inv=&Gi_inv[1+k*M], *g=&G[1]; cp(gi_inv,g,M); }

	Lat->ConnectG(Gi_inv,Gi,N,phi,PHI);

	for (s=N-1; s>=1; s--) {
		Lat->PropagateB(Gi_inv,G,S);
		Lat->ConnectG(Gi_inv,Gi,s,phi,PHI);
	}

	Lat->CorrectDoubleCountG(phi,PHI,G);
	lnGN = Lat->ComputeLnGN(Gi_inv);
	if (freedom == fixedTheta) {
		lnCt = log(theta)-lnGN-log(1.0*N);
		Lat->NormPhiRestr(phi,PHI,Gi_inv,theta/N); // factor 2 optimization correction
		phiBulk = theta/exp(lnGN);
		SetPhiBulk(phiBulk);
	} else if (freedom == rangeRestricted) {
		double phiPreNorm = 0;
		for (z=1; z<=M; z++) {
			if (restrictedRange->InRange(z)) {
				phiPreNorm += phi[z];
			}
		}
		lnCt = log(phiRange/(phiPreNorm));// factor 2 optimization correction
		Lat->NormPhiFree(phi,PHI,phiRange/phiPreNorm);
		theta = ComputeTheta();
		phiBulk = 6.0*exp(lnCt)*N;
		SetPhiBulk(phiBulk);
	} else {
		Lat->NormPhiFree(phi,PHI,phiBulk/(6*N));// factor 2 optimization correction
		lnCb = log(phiBulk/N);
	}
	Lat->RestoreFromSafe(G);
	Lat->RestoreFromSafe(phi);
}
void
SF_Homopol2ndO::Matrix1Long(const DensityPart DensPart) {//not tested yet; possibly works
	int k,z,s,s0,t0,v0,t,rs1;
	double S=GetPf();
	Matrix Gi(1,M*size,1,n);
	Vector Gi_inv(1,M*size), Gs(1,M*size);
	Vector phi = Seg->GetPhi(DensPart);
	Vector G = Seg->GetSWF();
	Lat->MakeSafe(G);
	Zero(phi,M);
	Zero(PHI,Mphi*M);
	for (k=0; k<size; k++) {double *gi = &Gi[1+k*M][1],*gs = &Gs[1+k*M], *g=&G[1]; cp(gi,g,M); cp(gs,g,M);}

    Lat->MakeSafe(phi);
	t=1;
	v0=t0=s0 = 0;
	for (s=2; s<=N; s++) {
		t++;
		if (t>n) {    //(++t>n)
			t0++;
			if (t0 == n) t0 = ++v0;
			t = t0 + 1;
			s0 = s - t0 - 1;
		}
		Lat->PropagateF(Gs,G,S);
		if ((t == t0+1 && t0 == v0)
		   || (t == t0+1 && ((n-t0)*(n-t0+1) >= N-1-2*(s0+t0)))
		   || (2*(n-t+s) >= N-1)) {
			   double *gi = &Gi[1][t],*gs = &Gs[1]; cp(gi,gs,size*M);
		}
	}

	for (s=N; s>=1; s++) {
		if (s == N) {
			for (k=0; k<size; k++) {double *gi_inv=&Gi_inv[1+k*M],*g=&G[1]; cp(gi_inv,g,M);}
		} else {
			Lat->PropagateB(Gi_inv,G,S);
		}
		t = s - s0;
		if (t == t0) {
			s0 += - n + t0;
			if (t0 == v0)
				s0 -= ((n - t0)*(n - t0 + 1))/2;
			t0 --;
			if (t0 < v0)
				v0 = t0;
			double *gs = &Gs[1],*gi = &Gi[1][t]; cp(gs,gi,size*M);
			for (rs1=s0+t0+2; rs1<=s; rs1++) {
				t++;
				Lat->PropagateF(Gs,G,S);
				if (t == t0+1 || s0+n == s) {
					double *gi = &Gi[1][t],*gs = &Gs[1]; cp(gi,gs,size*M);
				}
				if (t == n && s0+n < s) {
					t  = ++t0;
					s0 += n - t0;
				}
			}
			t = n;
		}
		Lat->ConnectG(Gi_inv,Gi,t,phi,PHI); // optimization, chain is symmetric
	}
	Lat->CorrectDoubleCountG(phi,PHI,G);
	lnGN = Lat->ComputeLnGN(Gi_inv);
	if (freedom == fixedTheta) {
		lnCt = log(theta)-lnGN-log(1.0*N);
		Lat->NormPhiRestr(phi,PHI,Gi_inv,2*theta/N); // factor 2 optimization correction
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
		Lat->NormPhiFree(phi,PHI,2*exp(lnCt));// factor 2 optimization correction
		theta = ComputeTheta();
		phiBulk = theta/exp(lnGN);
		SetPhiBulk(phiBulk);
	} else {
		Lat->NormPhiFree(phi,PHI,2*phiBulk/(6*N));// factor 2 optimization correction
		lnCb = log(phiBulk/N);
	}
	Lat->RestoreFromSafe(G);
	Lat->RestoreFromSafe(phi);
}
void
SF_Homopol2ndO::Matrix2ndGen(const LatticeRange* Range,
							 const DensityPart DensPart1,
							 const DensityPart DensPart2) {
	Message(debug,"Matrix2ndGen not implemented in second order");
}
void
SF_Homopol2ndO::Matrix2ndGenLong(const LatticeRange* Range,
								 const DensityPart DensPart1,
								 const DensityPart DensPart2) {
	Message(debug,"Matrix2ndGenLong  not implemented in second order");
}
void
SF_Homopol2ndO::Matrix3rdGen(const LatticeRange* LatRangeConstr,
							 const DensityPart DensPart1,
							 const DensityPart DensPart2,
							 const DensityPart DensPart3,
							 const LatticeRange* LatRangeRenormed) {
    Message(debug, "Matrix3rdGen not implemented in second order");
}
void
SF_Homopol2ndO::Matrix3rdGenLong(const LatticeRange* LatRangeConstr,
								 const DensityPart constr,
								 const DensityPart unconstr,
								 const DensityPart renormed,
								 const LatticeRange* LatRangeRenormed) {
	Message(debug, "Matrix3rdGenLong not implemented in second order");
}
void
SF_Homopol2ndO::MatrixBulk(const DensityPart DensPart) {
	Vector G = Seg->GetSWF();
	Vector phi = GetPhi(DensPart);
	for (int z=1; z<=M; z++) {
		if (Lat->BulkBoundary(z)) {
			phi[z] = phiBulk*pow(G[z],N);
		}
	}
}
void
SF_Homopol2ndO::CopyBulkBoundaries(const DensityPart DensPart) {
	Vector phi = GetPhi(DensPart);
	Vector PhiBulk = GetPhi(bulk);
	for (int z=1; z<=M; z++) {
		if (Lat->BulkBoundary(z)) {
			phi[z] = PhiBulk[z];
		}
	}
}


