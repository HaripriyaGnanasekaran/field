#include "SF_Dend1stO.h"

SF_Dend1stO::SF_Dend1stO(Text name_,
						   SF_SegmentList* SegQ_,
						   Lat1stO* Lat_,
						   Input* MyInput_)
: SF_Molecule(name_, SegQ_, Lat_, MyInput_)
{
	Lat = Lat_;

	force_set = (!GetForce()==0);
	if (force_set) {
		Message(fatal,MyInput,
				"Force ensemble not implemented in SF_Dend1stO");
	}
	if (GetGPU()) {
		Message(literal,MyInput,"For dendrimer 1st order the GPU card is not yet activated. Contact Frans Leermakers. Going classical instead");
	}

	if (Chain->GetMoleculeType() != dendrimer) {
		Message(fatal,MyInput,
		"Programming error, trying to create a dendrimer for molecule '" +
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
}
SF_Dend1stO::~SF_Dend1stO() {
	if (freedom == secondGeneration) {
		DeletePhi(constrained);
		DeletePhi(unconstrained);
	} else {
		DeletePhi(total);
	}
}
void
SF_Dend1stO::ReComputePhi() {
	Message(fatal,"ReComputePhi not implemented in SF_Dend1stO");
}
Vector
SF_Dend1stO::GetLong(){
	Message(fatal,"GetLong not implemented in SF_Homopol2ndO");
	Vector x;
	return x;
}
Vector
SF_Dend1stO::GetShort(){
	Message(fatal,"GetShort not implemented in SF_Homopol2ndO");
	Vector x;
	return x;
}
Vector
SF_Dend1stO::GetBondOrientation( const DensityPart DensPar) {
	Message(fatal,"GetBondOrientation not implemented in SF_Dend1stO");
	Vector x;
	return x; 
};
MoleculeType
SF_Dend1stO::GetMoleculeType() const {
	return Chain->GetMoleculeType();
}



void
SF_Dend1stO::ComputePhi() {
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
		Message(implementation, "thirdGeneration not implemented for a dendrimer in SF_Dend1stO");
	}
	if (freedom != secondGeneration) {
		if (!saveMemory) {
			Matrix1(total);
		} else {
			Message(implementation, "cannot save memory for a dendrimer in SF_Dend1stO");
		}
	} else {
		if (!saveMemory) {
			Matrix2ndGen(LatRange,constrained,unconstrained);
		} else {
			Message(implementation, "cannot save memory for a dendrimer in SF_Dend1stO");
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
SF_Dend1stO::GetLayerAnalysis(Output*) {
	if (LatRangeStartLoops == NULL || LatRangeTrains == NULL) {
		return;
	}
	Message(implementation,"SF_Dend1stO::GetLayerAnalysis");
}

// the dendrimer version of Matrix, using the dendrimer symmetry
void
SF_Dend1stO::Matrix1(const DensityPart DensPart) {
	int z,s,i,repeat;
	Vector phi,G;
	int Nshortcut = 0; // Nshortcut is distance from the origin to the end
	SF_SegmentBlock* SegBlock;
	for (i=1; i<=Chain->GetNumGenerations(); i++) {
		SegBlock = Chain->GetSegmentBlock(i);
		Nshortcut += SegBlock->GetMaxLength();
	}
	// SegBlock is now the last block of the last generation
	// the matrix size is Mx(distance from the centre to the end)
	Matrix Gi(1,M,1,Nshortcut);
	Vector Gi_inv(1,M);
	for (i=1; i<=numDiffSegments; i++) {
		Seg = Chain->GetDiffSegment(i);
		G = Seg->GetSWF();
		Lat->MakeSafe(G);
		phi = Seg->GetPhi(DensPart);
		for (z=1; z<=M; z++) {
			phi[z] = 0;
		}
		Lat->MakeSafe(phi);
	}
	s = 0;
	Vector Gi_add;
	Vector Gold;
	// SegBlock is now the last block of the last generation
	// descend back from the end towards the origin
	for (i=Chain->GetNumGenerations(); i>=1; i--) {
		SegBlock = Chain->GetSegmentBlock(i);
		// descend through different homo-blocks within each generation
		for (int j=SegBlock->GetMaxLength(); j>=1; j--) {
			Gold = G;
			G = SegBlock->GetSegment(j)->GetSWF();
			s++;
			if (s == 1) {
				for (z=1; z<=M; z++) {
					Gi[z][1] = G[z];
				}
			} else if (j==SegBlock->GetMaxLength() && Chain->GetNumRepeat(i+1) > 1) {
				// if we get to the node, repeat the block
				repeat = Chain->GetNumRepeat(i+1);
				Vector Gi_temp(1,M);
				for (z=1; z<=M; z++) {
					Gi_inv[z] = Gi[z][s-1]; // Gi_inv used locally here
					Gi_temp[z] = Gi[z][s-1];
				}
				for (int k=1; k<=repeat-1; k++) {
					Gi_add.Dim(1,M);
					Lat->MakeSafe(Gi_add);
					Lat->ConnectG(Gi_inv,Gi_temp,Gi_add);
					Lat->CorrectDoubleCountG(Gi_add,Gold);
					for (z=1; z<=M; z++) {
						Gi_temp[z] = Gi_add[z];
					}
				}
				Lat->PropagateG(Gi_add,G);
				for (z=1; z<=M; z++) {
					Gi[z][s] = Gi_add[z];
				}
			} else {
				Lat->PropagateG(Gi,G,s);
			}
		}
	}
	s = Nshortcut;
	double NormRepeat = 1;
	// the same from inside out with Gi_inv and computing Phi at the same time
	for (i=1; i<=Chain->GetNumGenerations(); i++) {
		SegBlock = Chain->GetSegmentBlock(i);
		repeat = Chain->GetNumRepeat(i);
		for (int j=1; j<=SegBlock->GetMaxLength(); j++) {
			Seg = SegBlock->GetSegment(j);
			G = Seg->GetSWF();
			if (s == Nshortcut) {
				for (z=1; z<=M; z++) {
					Gi_inv[z] = G[z];
				}
			} else {
				Lat->PropagateG(Gi_inv,G);
			}
			phi = Seg->GetPhi(DensPart);
			if (j == 1) { // branchpoint
				for (int k=1; k<=repeat-1; k++) {
					Vector Gi_add(1,M);
					Lat->MakeSafe(Gi_add);
					Lat->ConnectG(Gi_inv,Gi,s,Gi_add);
					Lat->CorrectDoubleCountG(Gi_add,G);
					for (z=1; z<=M; z++) {
						Gi_inv[z] = Gi_add[z];
					}
				}
				Lat->NormPhiFree(Gi_inv,NormRepeat);
				Lat->ConnectG(Gi_inv,Gi,s,phi);
				Lat->NormPhiFree(Gi_inv,1.0/NormRepeat);
				NormRepeat *= repeat; // NormRepeat = NormRepeat * repeat
			} else {
				Lat->NormPhiFree(Gi_inv,NormRepeat);
				Lat->ConnectG(Gi_inv,Gi,s,phi);
				Lat->NormPhiFree(Gi_inv,1.0/NormRepeat);
			}
			s--;
		}
	}
	// we have the G, we can compute phi
	lnGN = Lat->ComputeLnGN(Gi_inv);
	for (i=1; i<=numDiffSegments; i++) {
		Seg = Chain->GetDiffSegment(i);
		G = Seg->GetSWF();
		phi = Seg->GetPhi(DensPart);
		Lat->CorrectDoubleCountG(phi,G);
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
			Lat->NormPhiRestr(phi,Gi_inv,theta/N);
		} else if (freedom == rangeRestricted) {
			Lat->NormPhiFree(phi,exp(lnCt));
		} else {
			Lat->NormPhiFree(phi,phiBulk/N);
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
SF_Dend1stO::Matrix2ndGen(const LatticeRange* Range,
							const DensityPart DensPart1,
							const DensityPart DensPart2) {
	int z,s,i,repeat;
	Vector phi1,phi2,G;
	int Nshortcut = 0;
	SF_SegmentBlock* SegBlock;
	for (i=1; i<=Chain->GetNumGenerations(); i++) {
		SegBlock = Chain->GetSegmentBlock(i);
		Nshortcut += SegBlock->GetMaxLength();
	}
	Matrix Gi1(1,M,1,Nshortcut);
	Vector Gi_inv1(1,M);
	Matrix Gi2(1,M,1,Nshortcut);
	Vector Gi_inv2(1,M);
	Vector E(1,M);
	for (i=1; i<=numDiffSegments; i++) {
		Seg = Chain->GetDiffSegment(i);
		G = Seg->GetSWF();
		Lat->MakeSafe(G);
		phi1 = Seg->GetPhi(DensPart1);
		phi2 = Seg->GetPhi(DensPart2);
		for (z=1; z<=M; z++) {
			phi1[z] = 0;
			phi2[z] = 0;
			E[z] = 1;
		}
		Lat->MakeSafe(phi1);
		Lat->MakeSafe(phi2);
		Lat->MakeSafe(E);
	}
	s = 0;
	Vector Gi_add1;
	Vector Gi_add2;
	Vector Gold;
	for (i=Chain->GetNumGenerations(); i>=1; i--) {
		SegBlock = Chain->GetSegmentBlock(i);
		for (int j=SegBlock->GetMaxLength(); j>=1; j--) {
			Gold = G;
			G = SegBlock->GetSegment(j)->GetSWF();
			Lat->Init2G(Gi_inv1,Gi_inv2,G,Range); // Gi_inv used here to get values
			s++;
			if (s == 1) {
				for (z=1; z<=M; z++) {
					Gi1[z][1] = Gi_inv1[z];
					Gi2[z][1] = Gi_inv2[z];
				}
			} else if (j==SegBlock->GetMaxLength() && Chain->GetNumRepeat(i+1) > 1) {
				repeat = Chain->GetNumRepeat(i+1);
				Vector Gi_temp1(1,M);
				Vector Gi_temp2(1,M);
				for (z=1; z<=M; z++) {
					Gi_inv1[z] = Gi1[z][s-1]; // Gi_inv1 used locally here
					Gi_temp1[z] = Gi1[z][s-1];
					Gi_inv2[z] = Gi2[z][s-1]; // Gi_inv2 used locally here
					Gi_temp2[z] = Gi2[z][s-1];
				}
				for (int k=1; k<=repeat-1; k++) {
					Gi_add1.Dim(1,M);
					Lat->MakeSafe(Gi_add1);
					Lat->ConnectG(Gi_inv1,Gi_temp1,Gi_add1);
					Lat->Connect2G(Gi_inv1,Gi_temp1,Gi_inv2,Gi_temp2,Gi_add1);
					Lat->CorrectDoubleCountG(Gi_add1,Gold);
					Gi_add2.Dim(1,M);
					Lat->MakeSafe(Gi_add2);
					Lat->ConnectG(Gi_inv2,Gi_temp2,Gi_add2);
					Lat->CorrectDoubleCountG(Gi_add2,Gold);
					for (z=1; z<=M; z++) {
						Gi_temp1[z] = Gi_add1[z];
						Gi_temp2[z] = Gi_add2[z];
					}
				}
				Lat->Propagate2G(Gi_add1,Gi_add2,G,Range);
				for (z=1; z<=M; z++) {
					Gi1[z][s] = Gi_add1[z];
					Gi2[z][s] = Gi_add2[z];
				}
			} else {
				Lat->Propagate2G(Gi1,Gi2,G,s,Range);
			}
		}
	}
	s = Nshortcut;
	double NormRepeat = 1;
	for (i=1; i<=Chain->GetNumGenerations(); i++) {
		SegBlock = Chain->GetSegmentBlock(i);
		repeat = Chain->GetNumRepeat(i);
		for (int j=1; j<=SegBlock->GetMaxLength(); j++) {
			Seg = SegBlock->GetSegment(j);
			G = Seg->GetSWF();
			phi1 = Seg->GetPhi(DensPart1);
			phi2 = Seg->GetPhi(DensPart2);
			if (s == Nshortcut) {
				Lat->Init2G(Gi_inv1,Gi_inv2,G,Range);
			} else {
				Lat->Propagate2G(Gi_inv1,Gi_inv2,G,Range);
			}
			if (j == 1) { // branchpoint
				for (int k=1; k<=repeat-1; k++) {
					Vector Gi_add1(1,M);
					Lat->MakeSafe(Gi_add1);
					Lat->ConnectG(Gi_inv1,Gi1,s,Gi_add1);
					Lat->Connect2G(Gi_inv1,Gi1,s,Gi_inv2,Gi2,s,Gi_add1);
					Lat->CorrectDoubleCountG(Gi_add1,G);
					Vector Gi_add2(1,M);
					Lat->MakeSafe(Gi_add2);
					Lat->ConnectG(Gi_inv2,Gi2,s,Gi_add2);
					Lat->CorrectDoubleCountG(Gi_add2,G);
					for (z=1; z<=M; z++) {
						Gi_inv1[z] = Gi_add1[z];
						Gi_inv2[z] = Gi_add2[z];
					}
				}
				Lat->NormPhiFree(Gi_inv1,NormRepeat);
				Lat->NormPhiFree(Gi_inv2,NormRepeat);
				Lat->ConnectG(Gi_inv1,Gi1,s,phi1);
				Lat->Connect2G(Gi_inv1,Gi1,s,Gi_inv2,Gi2,s,phi1); //'tails'
				Lat->ConnectG(Gi_inv2,Gi2,s,phi2);
				Lat->NormPhiFree(Gi_inv1,1.0/NormRepeat);
				Lat->NormPhiFree(Gi_inv2,1.0/NormRepeat);
				NormRepeat *= repeat; // NormRepeat = NormRepeat * repeat
			} else {
				Lat->NormPhiFree(Gi_inv1,NormRepeat);
				Lat->NormPhiFree(Gi_inv2,NormRepeat);
				Lat->ConnectG(Gi_inv1,Gi1,s,phi1);
				Lat->Connect2G(Gi_inv1,Gi1,s,Gi_inv2,Gi2,s,phi1); //'tails'
				Lat->ConnectG(Gi_inv2,Gi2,s,phi2);
				Lat->NormPhiFree(Gi_inv1,1.0/NormRepeat);
				Lat->NormPhiFree(Gi_inv2,1.0/NormRepeat);
			}
			s--;
		}
	}
	lnGN = Lat->ComputeLnGN(Gi_inv1);
	if (!muFixed) {
		lnCt = log(theta) - lnGN - log(1.*N);
	}
	lnCb = log(phiBulk/N);
	for (i=1; i<=numDiffSegments; i++) {
		Seg = Chain->GetDiffSegment(i);
		G = Seg->GetSWF();
		phi1 = Seg->GetPhi(DensPart1);
		phi2 = Seg->GetPhi(DensPart2);
		Lat->CorrectDoubleCountG(phi1,G);
		Lat->CorrectDoubleCountG(phi2,G);
		if (muFixed) {
			Lat->NormPhiFree(phi1,exp(lnCt));
		} else {
			Lat->NormPhiRestr(phi1,Gi_inv1,theta/N);
		}
		Lat->NormPhiFree(phi2,phiBulk/N);
		Lat->RestoreFromSafe(G);
		Lat->RestoreFromSafe(phi1);
		Lat->RestoreFromSafe(phi2);
	}
}
void
SF_Dend1stO::MatrixBulk(const DensityPart DensPart) {
	int i,z;
	Vector phiTot = GetPhi(DensPart);
	for (z=1; z<=M; z++) {
		if (Lat->BulkBoundary(z)) {
			phiTot[z] = 0;
		}
	}
	for (i=1; i<=numDiffSegments; i++) {
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
SF_Dend1stO::CopyBulkBoundaries(const DensityPart DensPart) {
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

