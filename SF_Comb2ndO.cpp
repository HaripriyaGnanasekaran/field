#include "SF_Comb2ndO.h"
#include "tools.h"

SF_Comb2ndO::SF_Comb2ndO(Text name_,
						   SF_SegmentList* SegQ_,
						   Lat2ndO* Lat_,
						   Input* MyInput_)
: SF_Molecule(name_, SegQ_, Lat_, MyInput_)
{
	Lat = Lat_;
	M=Lat->GetTotalNumLayers();
	//if (MayerSaupeSet()) {
	//	Message(fatal,MyInput,"In SF_Comb2ndO: Mayer_Saupe not implemented for first order propagators");
	//}
	force_set = (!GetForce()==0);
	if (force_set) {
		Message(fatal,MyInput,
				"Force ensemble not implemented in SF_Comb2ndO");
	}
	if (GetGPU()) {
		Message(literal,MyInput,"For comb-polymer 2nd order the GPU card is not yet activated. Contact Frans Leermakers. Going classical instead");
	}

	if (Chain->GetMoleculeType() != comb) {
		Message(fatal,MyInput,
		"Programming error, trying to create a comb-chain for molecule '" +
		name + "'");
	}

	numDiffSegments = Chain->GetNumDiffSegments(); NDS=numDiffSegments;
	if (freedom == secondGeneration) {
		CreatePhi(constrained);
		CreatePhi(unconstrained);
	} else {
		CreatePhi(total);
	}

	SetPhiBulk(phiBulk);
	if (gradients==3) {size=6;}
	if (gradients==2) {size=5;}
	if (gradients==1) {size=3;}
	if (size==3) Mphi=2; else Mphi=3;
	PHI = new double[Mphi*NDS*M];
	labda.Dim(0,5);
	labda[0]=labda[1]=labda[2]=labda[3]=labda[4]=labda[5]=1.0/6.0;
	if (size==3) labda[1] *=4.0;
	if (size==5) labda[2] *=2.0;

	//numGenerations=Chain->GetNumGenerations(); NG=numGenerations;
	numArms=Chain->GetNumArms(); NG=numArms;
	//side_type = 1 : linear sides
	//side_type = 2 : dendrimers
	//side_type = 3 : asymmetric dendrimers
	//side_type = 4 : branched;
	side_type = Chain->GetSideType();
	if (side_type>1) Message(fatal,MyInput,"Only linear sides implemented in 2nd order comb calculations...contact FL when this needs to be changed.");

	SegBlock1 = Chain->GetSegmentBlockComb(1); N1=SegBlock1->GetMaxLength();
	SegBlock2 = Chain->GetSegmentBlockComb(2); N2=SegBlock2->GetMaxLength();
	SegBlock3 = Chain->GetSegmentBlockComb(3); N3=SegBlock3->GetMaxLength();
	SegBlock4 = Chain->GetSegmentBlockComb(4); N4=SegBlock4->GetMaxLength();
	phi_long.Dim(1,M); //backbone
	phi_short.Dim(1,M); //sides
}
SF_Comb2ndO::~SF_Comb2ndO() {
	if (freedom == secondGeneration) {
		DeletePhi(constrained);
		DeletePhi(unconstrained);
	} else {
		DeletePhi(total);
	}

	delete [] PHI;
}
void
SF_Comb2ndO::ReComputePhi() {
	Message(fatal,"ReComputePhi not implemented in SF_Comb2ndO");
}
void
SF_Comb2ndO::ComputePhi() {
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
	if (freedom == thirdGeneration) {
		Message(implementation, "thirdGeneration not implemented for a dendrimer in SF_Comb2ndO");
	}
	if (freedom != secondGeneration) {
		if (!saveMemory) {
			Matrix1(total);
		} else {
			Message(implementation, "cannot save memory for a dendrimer in SF_Comb2ndO");
		}
	} else {
		if (!saveMemory) {
			Matrix2ndGen(LatRange,constrained,unconstrained);
		} else {
			Message(implementation, "cannot save memory for a dendrimer in SF_Comb2ndO");
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
	int z;
	if (freedom == secondGeneration && numDiffSegments > 1) {
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
	} else if (numDiffSegments > 1) {
		Vector phiTot = GetPhi(total);
		for (z=1; z<=M; z++) {
			phiTot[z] = 0;
		}
		for (i=1; i<=numDiffSegments; i++) {
			Seg = Chain->GetDiffSegment(i);
			Vector phi = Seg->GetPhi(total);
			for (z=1; z<=M; z++) {
				phiTot[z] += phi[z];
			}
		}
	}
	if (freedom == secondGeneration) {
		if (muFixed) {
			Vector phi1Tot = GetPhi(constrained);;
			Lat->SubtractBoundaries(phi1Tot);
			Lat->MultiplyWithLatticeSites(phi1Tot);
			theta = 0;
			for (int z=1; z<=M; z++) {
				theta += phi1Tot[z];
			}
			Lat->DivideByLatticeSites(phi1Tot);
			Lat->RestoreBoundaries(phi1Tot);
		}
	}
}
void
SF_Comb2ndO::GetLayerAnalysis(Output*) {
	if (LatRangeStartLoops == NULL || LatRangeTrains == NULL) {
		return;
	}
	Message(implementation,"SF_Comb2ndO::GetLayerAnalysis");
}
Vector
SF_Comb2ndO::GetLong() {
	return phi_long;
}
Vector
SF_Comb2ndO::GetShort() {
	return phi_short;
}

Vector
SF_Comb2ndO::GetBondOrientation( const DensityPart DensPar) {
	Message(fatal,"GetBondOrientation not implemented in SF_Comb2ndO");
	Vector x;
	return x;
};

MoleculeType
SF_Comb2ndO::GetMoleculeType() const {
	return Chain->GetMoleculeType();
}

double SF_Comb2ndO::GetPf(SF_SegmentBlock* SegBlock, int s1,int s2) {
	int NumPfs = SegQ->GetPf_count();
	Vector Pf = SegQ->GetPfVector();
	Array <Text> Pfx = SegQ->GetPfx();
	Array <Text> Pfy = SegQ->GetPfy();
	Text name1 = SegBlock->GetSegment(s1)->GetName();
	Text name2 = SegBlock->GetSegment(s2)->GetName();

	for (int i=1; i<=NumPfs; i++){
		if ((Pfx[i]==name1 && Pfy[i] == name2) || (Pfy[i]==name1 && Pfx[i] == name2)) {
		return Pf[i];}
	}
	return 1.0/5.0;
}

double SF_Comb2ndO::GetPf(SF_SegmentBlock* SegBlock1,SF_SegmentBlock* SegBlock2, int s1,int s2) {
	int NumPfs = SegQ->GetPf_count();
	Vector Pf = SegQ->GetPfVector();
	Array <Text> Pfx = SegQ->GetPfx();
	Array <Text> Pfy = SegQ->GetPfy();
	Text name1 = SegBlock1->GetSegment(s1)->GetName();
	Text name2 = SegBlock2->GetSegment(s2)->GetName();

	for (int i=1; i<=NumPfs; i++){
		if ((Pfx[i]==name1 && Pfy[i] == name2) || (Pfy[i]==name1 && Pfx[i] == name2)) {
		return Pf[i];}
	}
	return 1.0/5.0;
}


double SF_Comb2ndO::GetPf(SF_MolSegment* Seg1,SF_MolSegment* Seg2) {
	int NumPfs = SegQ->GetPf_count();
	Vector Pf = SegQ->GetPfVector();
	Array <Text> Pfx = SegQ->GetPfx();
	Array <Text> Pfy = SegQ->GetPfy();
	Text name1 = Seg1->GetName();
	Text name2 = Seg2->GetName();

	for (int i=1; i<=NumPfs; i++){
		if ((Pfx[i]==name1 && Pfy[i] == name2) || (Pfy[i]==name1 && Pfx[i] == name2)) {
		return Pf[i];}
	}
	return 1.0/5.0;
}


void
SF_Comb2ndO::Matrix1(const DensityPart DensPart) {
	int z,s,ss,i,k;
	Vector phi,G;
	Matrix Gi(1,M*size,1,N1+N4+N2*NG);
	Matrix gi(1,M*size,1,N3+1);
	Vector Gi_inv(1,M*size);
	Vector Gi_temp(1,M*size);
	Vector Gb(1,M);

	for (i=1; i<=numDiffSegments; i++) {
		Seg = Chain->GetDiffSegment(i);
		G = Seg->GetSWF();
		Lat->MakeSafe(G);
		phi = Seg->GetPhi(DensPart);
		for (z=1; z<=M; z++) {phi[z] = 0;}
		Lat->MakeSafe(phi);
	}
	//Zero(PHI,Mphi*NDS*M);
	Zero(Gb,M);
	Zero(phi_long,M);
	Zero(phi_short,M);
	double *ig;
	double *g;
	double *iG;
	s = 0;
	ss = 0;
	for (int jj=N3; jj>=1; jj--) { //first do the side chain propagation from free end to branchpoint
		ss++;
		Seg_old=Seg;
		Seg=SegBlock3->GetSegment(jj);
		G=Seg->GetSWF();
		if (jj==N3) {
			g=&G[1];
			for (k=0; k<size; k++) {ig=&gi[1+k*M][ss]; cp(ig,g,M); }
		} else {
			Lat->PropagateF(gi,G,ss,GetPf(Seg_old,Seg));
		}
	}

	for (int j=1; j<=N1; j++) {
		s++;
		Seg_old=Seg;
		Seg=SegBlock1->GetSegment(j);
		G=Seg->GetSWF(); g=&G[1];
		if (s==1) {
			for (k=0; k<size; k++) {iG=&Gi[1+k*M][s]; cp(iG,g,M); }
		} else {
			Lat->PropagateF(Gi,G,s,GetPf(Seg_old,Seg));
		}
	}
	for (i=1; i<=NG; i++) {
		ss=0;
		for (int j=1; j<=N2; j++){
			s++; ss++;
			Seg_old=Seg;
			Seg=SegBlock2->GetSegment(ss);
			G=Seg->GetSWF();
			if (ss==1 && i>1) {
				Lat->PropagateF(Gi_temp,G,GetPf(Seg_old,Seg));
				iG=&Gi[1][s]; g=&Gi_temp[1]; cp(iG,g,M*size);
			} else {
				Lat->PropagateF(Gi,G,s,GetPf(Seg_old,Seg));
			}
		}
		if (i==1) {
			Lat->PropagateF(gi,G,N3+1,0.2); //compute one extra step for the side chain
			for (z=1; z<=M; z++) {
				for (k=0; k<size; k++) Gb[z] += gi[z+k*M][N3+1]*labda[k];
				if (G[z]>0) Gb[z] /=G[z];
			}
		}
		for (z=1; z<=M; z++) for (k=0; k<size; k++) Gi_temp[z+k*M]=Gb[z]*Gi[z+k*M][s];
	}
	ss=0;
	for (int j=1; j<=N4; j++) {
	 	ss++; s++;
	 	Seg_old=Seg;
	 	Seg=SegBlock4->GetSegment(j);
	 	G=Seg->GetSWF();
	 	if (ss==1) {
			Lat->PropagateF(Gi_temp,G,GetPf(Seg_old,Seg));
			iG= &Gi[1][s]; g=&Gi_temp[1]; cp(iG,g,M*size);
		}
		else {
			Lat->PropagateF(Gi,G,s,GetPf(Seg_old,Seg));
		}
	}

	s=N1+N4+N2*NG+1;
	for (int j=N4; j>=1; j--) {
		s--;
		Seg_old=Seg;
		Seg=SegBlock4->GetSegment(j);
		G=Seg->GetSWF();
		phi=Seg ->GetPhi(DensPart);
		if (j==N4) {
			g=&G[1];
			for (k=0; k<size; k++) {iG=&Gi_inv[1+k*M]; cp(iG,g,M); }
		} else {
			Lat->PropagateB(Gi_inv,G,GetPf(Seg_old,Seg));
		}
		Lat->ConnectG(Gi_inv,Gi,s,phi);
	}
	for (i=NG; i>=1; i--) {
		for (int j=N2; j>=1; j--) {
			s--;
			Seg_old=Seg;
			Seg = SegBlock2->GetSegment(j);
			G=Seg->GetSWF();
			phi = Seg ->GetPhi(DensPart);
			Lat->PropagateB(Gi_inv,G,GetPf(Seg_old,Seg));
			if (j==N2) {
				for (z=1; z<=M; z++) for (k=0; k<size; k++) Gi_temp[z+k*M]=Gi_inv[z+k*M]*Gb[z];
				Lat->ConnectG(Gi_temp,Gi,s,phi);
				for (z=1; z<=M; z++) for (k=0; k<size; k++)
				if (G[z]>0) Gi[z+k*M][s] =Gi[z+k*M][s]*Gi_inv[z+k*M]/G[z]; //prepare for phi side
				for (z=1; z<=M*size; z++) Gi_inv[z] = Gi_temp[z]; //prepare for further propagation
			} else {
				Lat->ConnectG(Gi_inv,Gi,s,phi);
			}
		}
	}

	for (int j=N1; j>=1; j--) {
		s--;
		Seg_old=Seg;
		Seg=SegBlock1->GetSegment(j);
		G=Seg->GetSWF();
		phi=Seg ->GetPhi(DensPart);
		Lat->PropagateB(Gi_inv,G,GetPf(Seg_old,Seg));
		Lat->ConnectG(Gi_inv,Gi,s,phi);
	}
	lnGN = Lat->ComputeLnGN(Gi_inv);

	for (i=1; i<=numDiffSegments; i++) {
		Seg_old=Seg;
		Seg = Chain->GetDiffSegment(i);
		G = Seg->GetSWF();
		phi = Seg->GetPhi(DensPart);
	}

	for (z=1; z<=M; z++) {
		if (G[z]>0) phi_long[z] +=phi[z]/G[z]; //backbone
	}

	for (i=1; i<=NG; i++) {
		Zero(Gb,M);
		for (z=1; z<=M; z++) {
			for (k=0; k<size; k++) Gb[z] += Gi[z+k*M][N1+i*N2]*labda[k];
		}
		for (z=1; z<=M; z++) for (k=0; k<size; k++) Gi_inv[z+k*M]=Gb[z];
		for (int j=1; j<=N3; j++){
			Seg_old=Seg;
			Seg =SegBlock3->GetSegment(j);
			G=Seg->GetSWF();
			phi = Seg->GetPhi(DensPart);
			double Pf;
			if (j==1) Pf =1.0/5.0; else Pf = GetPf(Seg_old,Seg);
			Lat->PropagateB(Gi_inv,G,Pf);
			Lat->ConnectG(Gi_inv,gi,N3-j+1,phi);
		}
	}

	for (i=1; i<=numDiffSegments; i++) {
		Seg_old=Seg;
		Seg = Chain->GetDiffSegment(i);
		G = Seg->GetSWF();
		phi = Seg->GetPhi(DensPart);
		//Lat->CorrectDoubleCountG(phi,G);
		for (int z=1; z<=M; z++) {
			if (G[z] > 0) {
				phi[z] /=G[z];
				phi_short[z] +=phi[z]; //teeth + backbone
			}
		}
	}
	for (int z=1; z<=M; z++) {
		phi_short[z] -= phi_long[z]; //teeth only
	}
	// and do the various renormalizations
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
			//Lat->NormPhiRestr(phi,Gi_inv,theta/N);
			double C=theta/(6*N*exp(lnGN));
			for (int z=1; z<=M; z++) {
				phi[z] *= C;
				phi_long[z] *=C;
				phi_short[z]*=C;
			}
		} else if (freedom == rangeRestricted) {
			//Lat->NormPhiFree(phi,exp(lnCt));
			double C=exp(lnCt);
			for (int z=1; z<=M; z++) {
				phi[z] *= C;
				phi_long[z] *=C;
				phi_short[z]*=C;
			}
		} else {
			//Lat->NormPhiFree(phi,phiBulk/N);
			double C=phiBulk/(6*N);
			for (int z=1; z<=M; z++) {
				phi[z] *= C;
				phi_long[z] *=C;
				phi_short[z]*=C;
			}

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
SF_Comb2ndO::Matrix2ndGen(const LatticeRange* Range,
							const DensityPart DensPart1,
							const DensityPart DensPart2) {

}
void
SF_Comb2ndO::MatrixBulk(const DensityPart DensPart) {

}
void
SF_Comb2ndO::CopyBulkBoundaries(const DensityPart DensPart) {
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

