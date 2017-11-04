#include "SF_AsymDend1stO.h"
#include "tools.h"

SF_AsymDend1stO::SF_AsymDend1stO(Text name_,
						   SF_SegmentList* SegQ_,
						   Lat1stO* Lat_,
						   Input* MyInput_)
: SF_Molecule(name_, SegQ_, Lat_, MyInput_)
{
	Lat = Lat_;
	//if (MayerSaupeSet()) {
	//	Message(fatal,MyInput,"In SF_AsymDend1stO: Mayer_Saupe not implemented for first order propagators");
	//}
	force_set = (!GetForce()==0);
	if (force_set) {
		Message(fatal,MyInput,
				"Force ensemble not implemented in SF_AsymDend1stO");
	}
	if (GetGPU()) {
		Message(literal,MyInput,"For Asymdendrimer 1st order the GPU card is not yet activated. Contact Frans Leermakers. Going classical instead");
	}

	if (Chain->GetMoleculeType() != asymmetric_dendrimer) {
		Message(fatal,MyInput,
		"Programming error, trying to create a Asymmetric dendrimer for molecule '" +
		name + "'");
	}

	numDiffSegments = Chain->GetNumDiffSegments();
	if (freedom == secondGeneration) {
		CreatePhi(constrained);
		CreatePhi(unconstrained);
	} else {
		CreatePhi(total);
	}


	SetPhiBulk(phiBulk);

	numGenerations=Chain->GetNumGenerations();
	cumlengths_short = new int [numGenerations+1];cumlengths_short[0]=0;
	cumlengths_long = new int [numGenerations+1]; cumlengths_long[0]=0;

	SF_SegmentBlock* SegBlock1=NULL;
	SF_SegmentBlock* SegBlock2=NULL;
	shortest=1; longest=2;
	for (int i=numGenerations; i>=1; i--) {
		SegBlock1 = Chain->GetSegmentBlock(i,shortest);
		SegBlock2 = Chain->GetSegmentBlock(i,longest);
		cumlengths_short[numGenerations-i+1] = cumlengths_short[numGenerations-i]+SegBlock1->GetMaxLength();
		cumlengths_long[numGenerations-i+1]  = cumlengths_long[numGenerations-i]+ SegBlock2->GetMaxLength();
	}

	//for (int i=1; i<=numGenerations; i++){
	//	cout << "i=" << i << "short" << cumlengths_short[i] << "long = " << cumlengths_long[i] << endl;
	///}
	M=Lat->GetTotalNumLayers();
	NDS=numDiffSegments;
	NG=numGenerations;
	phi_long.Dim(1,(NG+1)*M);
	phi_short.Dim(1,2*M);
	PHI = new double[2*NDS*M];
	PHI_ = new double [(NG+1)*NDS*M];


}
SF_AsymDend1stO::~SF_AsymDend1stO() {
	if (freedom == secondGeneration) {
		DeletePhi(constrained);
		DeletePhi(unconstrained);
	} else {
		DeletePhi(total);
	}
	delete [] cumlengths_short;
	delete [] cumlengths_long;
	delete [] PHI;
	delete [] PHI_;

}
void
SF_AsymDend1stO::ReComputePhi() {
	Message(fatal,"ReComputePhi not implemented in SF_AsymDend1stO");
}
void
SF_AsymDend1stO::ComputePhi() {
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
		Message(implementation, "thirdGeneration not implemented for a dendrimer in SF_AsymDend1stO");
	}
	if (freedom != secondGeneration) {
		if (!saveMemory) {
			Matrix1(total);
		} else {
			Message(implementation, "cannot save memory for a dendrimer in SF_AsymDend1stO");
		}
	} else {
		if (!saveMemory) {
			Matrix2ndGen(LatRange,constrained,unconstrained);
		} else {
			Message(implementation, "cannot save memory for a dendrimer in SF_AsymDend1stO");
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
	Zero(phi_short,2*M);
	Zero(phi_long,(NG+1)*M);

	for (i=1; i<=numDiffSegments; i++) {

		Seg = Chain->GetDiffSegment(i);
		for (int z=0; z<M; z++) {
			phi_short[z+1]+=PHI[M*2*(i-1)+z];
			phi_short[M+z+1]+=PHI[M*2*(i-1)+M+z];
			for (int p=0; p<=NG; p++) phi_long[p*M+z+1] += PHI_[(i-1)*M*(NG+1)+p*M+z];
		}
	}
}
void
SF_AsymDend1stO::GetLayerAnalysis(Output*) {
	if (LatRangeStartLoops == NULL || LatRangeTrains == NULL) {
		return;
	}
	Message(implementation,"SF_AsymDend1stO::GetLayerAnalysis");
}
Vector
SF_AsymDend1stO::GetLong() {
	return phi_long;
}
Vector
SF_AsymDend1stO::GetShort() {
	return phi_short;
}

Vector
SF_AsymDend1stO::GetBondOrientation( const DensityPart DensPar) {
	Message(fatal,"GetBondOrientation not implemented in SF_AsymDend1stO");
	Vector x;
	return x; 
};

MoleculeType
SF_AsymDend1stO::GetMoleculeType() const {
	return Chain->GetMoleculeType();
}

void
SF_AsymDend1stO::ComputePhi(int i, Matrix Gi, Matrix gi, Vector Gb, Vector Gi_inv, bool short_, bool long_, int nshort, int nlong, const DensityPart DensPart) {
	if (i>numGenerations) return;
	int I=0;
	double phiz;
	int length_short,length_long;
	SF_MolSegment* SEG;
	bool long__=long_;
	bool short__=short_;
	int Nshort =nshort;
	int Nlong = nlong;
	SF_SegmentBlock* SegBlock;
	Vector phi,G;
	Vector Gi_temp(1,M); for (int z=1; z<=M; z++) Gi_temp[z]=Gb[z];
	SegBlock2 = Chain->GetSegmentBlock(i,longest); length_long=SegBlock2->GetMaxLength();
	SegBlock1 = Chain->GetSegmentBlock(i,shortest); length_short=SegBlock1->GetMaxLength();
	for (int s=1; s<=length_long; s++) {
		Seg = SegBlock2->GetSegment(s);
		for (int q=1; q<=numDiffSegments; q++) {
			SEG = Chain->GetDiffSegment(q);
			if (SEG->GetName()==Seg->GetName()) I=q;
		}
		phi = Seg->GetPhi(DensPart);
		G = Seg->GetSWF();
		if (s==1) {
			for (int z=1; z<=M; z++) Gi_inv[z]=Gi_temp[z]*gi[z][cumlengths_short[numGenerations-i+1]]; //omdat Gb al dedeeld is door gz hoeft dit niet meer.
		} else {
			Lat->PropagateG(Gi_inv,G);
			for (int z=1; z<=M; z++) {
				phiz=Gi_inv[z]*Gi[z][cumlengths_long[numGenerations-i+1]-s+1];
				phi[z]+=phiz;
				if (long__) PHI[(I-1)*2*M+M+z-1] += phiz;
				for (int p=Nshort+1; p<=NG-Nlong; p++) PHI_[(I-1)*(NG+1)*M + p*M + z - 1] +=phiz;
			}
		}
	}
	if (i<numGenerations) {
		SegBlock = Chain->GetSegmentBlock(i+1,longest);
		Seg = SegBlock->GetSegment(1);
		for (int q=1; q<=numDiffSegments; q++) {
			SEG = Chain->GetDiffSegment(q);
			if (SEG->GetName()==Seg->GetName()) I=q;
		}
		G = Seg->GetSWF();
		phi = Seg->GetPhi(DensPart);
		Lat->PropagateG(Gi_inv,G);
		Lat->CorrectDoubleCountG(Gi_inv,G);

		for (int z=1; z<=M; z++) {
			Gb[z]=Gi_inv[z];
			phiz=Gi[z][cumlengths_long[numGenerations-i]]*gi[z][cumlengths_short[numGenerations-i]]*Gi_inv[z];
			phi[z]+=phiz;
			if (long__) PHI[(I-1)*2*M+M+z-1]+=phiz;
			for (int p=Nshort+1; p<=NG-Nlong; p++) PHI_[(I-1)*(NG+1)*M + p*M + z - 1] +=phiz;
		}
		ComputePhi(i+1,Gi,gi,Gb,Gi_inv, false, long__, Nshort+1, Nlong, DensPart);
	}
	for (int s=1; s<=length_short; s++) {
		Seg = SegBlock1->GetSegment(s);
		for (int q=1; q<=numDiffSegments; q++) {
			SEG = Chain->GetDiffSegment(q);
			if (SEG->GetName()==Seg->GetName()) I=q;
		}
		phi = Seg->GetPhi(DensPart);
		G = Seg->GetSWF();
		if (s==1) {
			for (int z=1; z<=M; z++) Gi_inv[z]=Gi_temp[z]*Gi[z][cumlengths_long[numGenerations-i+1]];
		} else {
			Lat->PropagateG(Gi_inv,G);
			for (int z=1; z<=M; z++) {
				phiz=Gi_inv[z]*gi[z][cumlengths_short[numGenerations-i+1]-s+1];
				phi[z] +=phiz;
				if (short__) PHI[(I-1)*2*M+z-1] +=phiz;
				for (int p=Nshort; p<=NG-Nlong-1; p++) PHI_[(I-1)*(NG+1)*M + p*M + z - 1] +=phiz;
			}
		}
	}
	if (i<numGenerations) {
		SegBlock = Chain->GetSegmentBlock(i+1,shortest);
		Seg = SegBlock->GetSegment(1);
		for (int q=1; q<=numDiffSegments; q++) {
			SEG = Chain->GetDiffSegment(q);
			if (SEG->GetName()==Seg->GetName()) I=q;
		}
		G = Seg->GetSWF();
		phi = Seg->GetPhi(DensPart);
		Lat->PropagateG(Gi_inv,G);
		Lat->CorrectDoubleCountG(Gi_inv,G);

		for (int z=1; z<=M; z++) {
			Gb[z]=Gi_inv[z];
			phiz=Gi[z][cumlengths_long[numGenerations-i]]*gi[z][cumlengths_short[numGenerations-i]]*Gi_inv[z];
			phi[z]+=phiz;
			if (short__) PHI[(I-1)*2*M + z-1]+=phiz;
			for (int p=Nshort; p<=NG-Nlong-1; p++) PHI_[(I-1)*(NG+1)*M + p*M + z - 1] +=phiz;
		}
		ComputePhi(i+1,Gi,gi,Gb,Gi_inv, short__,false, Nshort, Nlong+1, DensPart);
	}
}

void
SF_AsymDend1stO::Matrix1(const DensityPart DensPart) {
	int z,s,ss,i, length_short=0, length_long=0, I=0;
	Vector phi,G;
	Matrix Gi(1,M,1,cumlengths_long[numGenerations]);
	Matrix gi(1,M,1,cumlengths_short[numGenerations]);
	Vector Gb(1,M);
	Vector Gi_inv(1,M);
	Vector Gi_temp(1,M);
	SF_MolSegment* SEG;
	for (i=1; i<=numDiffSegments; i++) {
		Seg = Chain->GetDiffSegment(i);
		G = Seg->GetSWF();
		Lat->MakeSafe(G);
		phi = Seg->GetPhi(DensPart);
		for (z=1; z<=M; z++) {phi[z] = 0;}
		Lat->MakeSafe(phi);
	}
	Zero(PHI,2*numDiffSegments*M);
	Zero(PHI_,(NG+1)*numDiffSegments*M);
	s = 0;
	ss = 0;
	Vector Gi_add;
		Gi_add.Dim(1,M);
		Lat->MakeSafe(Gi_add);
	Vector Gold;
	for (i=numGenerations; i>=1; i--) {
		SegBlock2 = Chain->GetSegmentBlock(i,longest); length_long=SegBlock2->GetMaxLength();
		SegBlock1 = Chain->GetSegmentBlock(i,shortest); length_short=SegBlock1->GetMaxLength();

		for (int j=length_long; j>=1; j--) {
			Gold = G;
			G = SegBlock2->GetSegment(j)->GetSWF();
			s++;
			if (s == 1) {
				for (z=1; z<=M; z++) Gi[z][1] = G[z];
			} else if (j==length_long) {

				for (z=1; z<=M; z++) {
					Gi_inv[z] = Gi[z][s-1]; // Gi_inv used locally here
					Gi_temp[z] = gi[z][ss]; //ss van vorige loop
				}
				for (z=1; z<=M; z++) Gi_add[z]=0;
				Lat->ConnectG(Gi_inv,Gi_temp,Gi_add);
				Lat->CorrectDoubleCountG(Gi_add,Gold);
				for (z=1; z<=M; z++) {Gi_temp[z] = Gi_add[z];} //used in the other loop.

				Lat->PropagateG(Gi_add,G);
				for (z=1; z<=M; z++) {Gi[z][s] = Gi_add[z];	}
			} else {
				Lat->PropagateG(Gi,G,s);
			}
		}
		for (int k=length_short; k>=1; k--) {
			G = SegBlock1->GetSegment(k)->GetSWF();
			ss++;
			if (ss ==1) {for (z=1; z<=M; z++) gi[z][1] = G[z];}
			else if (k == length_short) {
				Lat->PropagateG(Gi_temp,G);
				for (z=1; z<=M; z++) {gi[z][ss] = Gi_temp[z];}
			} else {
				Lat->PropagateG(gi,G,ss);
			}
		}
	}
	Seg = SegBlock1->GetSegment(1);
	phi = Seg ->GetPhi(DensPart);
	G = Seg->GetSWF();
	for (int q=1; q<=numDiffSegments; q++) {
		SEG = Chain->GetDiffSegment(q);
		if (SEG->GetName()==Seg->GetName()) I=q;
	}
	int shorts = cumlengths_short[numGenerations];
	int longs = cumlengths_long[numGenerations];
	for (int z=1; z<=M; z++) Gb[z]=1;
	for (int z=1; z<=M; z++) {
		phi[z]+=gi[z][shorts]*Gi[z][longs];
		Gi_inv[z]=phi[z];
		PHI[2*M*(I-1)+z-1]=phi[z];
		PHI[2*M*(I-1)+M+z-1]=phi[z];
		for (int p=0; p<=NG; p++) PHI_[(I-1)*(NG+1)*M + p*M + z - 1]=phi[z];
	}
	ComputePhi(1,Gi,gi,Gb,Gi_temp,true,true, 0, 0, DensPart);
	Lat->CorrectDoubleCountG(Gi_inv,G);
	lnGN = Lat->ComputeLnGN(Gi_inv);
	for (i=1; i<=numDiffSegments; i++) {
		Seg = Chain->GetDiffSegment(i);
		G = Seg->GetSWF();
		phi = Seg->GetPhi(DensPart);
		//Lat->CorrectDoubleCountG(phi,G);
		for (int z=1; z<=M; z++) {
			if (G[z] > 0) {
				phi[z] /=G[z];
				PHI[(i-1)*2*M+z-1] /= G[z];
				PHI[(i-1)*2*M+M+z-1] /= G[z];
				for (int p=0; p<=NG; p++) PHI_[(i-1)*(NG+1)*M + p*M + z - 1]/=G[z];
			}
		}
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
			double C=theta/(N*exp(lnGN));
			for (int z=1; z<=M; z++) {
				phi[z] *= C;
				PHI[(i-1)*2*M+z-1] *=C;
				PHI[(i-1)*2*M+M+z-1] *=C;
				for (int p=0; p<=NG; p++) PHI_[(i-1)*(NG+1)*M + p*M + z - 1]*=C;
			}
		} else if (freedom == rangeRestricted) {
			//Lat->NormPhiFree(phi,exp(lnCt));
			double C=exp(lnCt);
			for (int z=1; z<=M; z++) {
				phi[z] *= C;
				PHI[(i-1)*2*M+z-1] *=C;
				PHI[(i-1)*2*M+M+z-1] *=C;
				for (int p=0; p<=NG; p++) PHI_[(i-1)*(NG+1)*M + p*M + z - 1]*=C;
			}
		} else {
			//Lat->NormPhiFree(phi,phiBulk/N);
			double C=phiBulk/N;
			for (int z=1; z<=M; z++) {
				phi[z] *= C;
				PHI[(i-1)*2*M+z-1] *=C;
				PHI[(i-1)*2*M+M+z-1] *=C;
				for (int p=0; p<=NG; p++) PHI_[(i-1)*(NG+1)*M + p*M + z - 1]*=C;
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
SF_AsymDend1stO::Matrix2ndGen(const LatticeRange* Range,
							const DensityPart DensPart1,
							const DensityPart DensPart2) {

}
void
SF_AsymDend1stO::MatrixBulk(const DensityPart DensPart) {

}
void
SF_AsymDend1stO::CopyBulkBoundaries(const DensityPart DensPart) {
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

