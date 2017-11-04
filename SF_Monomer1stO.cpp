#include "SF_Monomer1stO.h"


SF_Monomer1stO::SF_Monomer1stO(Text name_,
								   SF_SegmentList* SegQ_,
								   Lat1stO* Lat_,
								   Input* MyInput_)
: SF_Molecule(name_, SegQ_, Lat_, MyInput_)
{
	if (Chain->GetMoleculeType() != monomer) {
		Message(fatal,MyInput,
		"Programming error, trying to create a monomer for molecule '" +
		name + "'");
	}
	Lat = Lat_;

	Seg = Chain->GetSegment(1);
	if (freedom == secondGeneration) {
		CreatePhi(constrained);
		CreatePhi(unconstrained);
	} else if (freedom == thirdGeneration) {
		CreatePhi(renorm);
		CreatePhi(constrained);
		CreatePhi(unconstrained);
	} else {
		CreatePhi(total);
	}
	SetPhiBulk(phiBulk);
	RenormQ.Dim(1,1);
}
SF_Monomer1stO::~SF_Monomer1stO() {
	if (freedom == secondGeneration) {
		DeletePhi(constrained);
		DeletePhi(unconstrained);
	} else if (freedom == thirdGeneration) {
		DeletePhi(renorm);
		DeletePhi(constrained);
		DeletePhi(unconstrained);
	} else {
		DeletePhi(total);
	}
}

void
SF_Monomer1stO::ReComputePhi() {
	Message(fatal,"ReComputePhi not implemented in SF_Monomer1stO");
}
Vector
SF_Monomer1stO::GetLong(){
	Message(fatal,"GetLong not implemented in SF_Branched1stO");
	Vector x;
	return x;
}
Vector
SF_Monomer1stO::GetShort(){
	Message(fatal,"GetShort not implemented in SF_Branched1stO");
	Vector x;
	return x;
}

Vector
SF_Monomer1stO::GetBondOrientation( const DensityPart DensPar) {
	Message(fatal,"GetBondOrientation not implemented in SF_Branched1stO");
	Vector x;
	return x; 
}; 

MoleculeType
SF_Monomer1stO::GetMoleculeType() const {
	return Chain->GetMoleculeType();
}

void
SF_Monomer1stO::ComputePhi() {
	if (Lat->GetNamesBulkBoundaries().Upperbound() > 0) {
		CreatePhi(bulk);
		MatrixBulk(bulk);
	}
	Vector G = Seg->GetSWF();
	Lat->SetBoundaries(G);
	if (freedom == thirdGeneration) { // broken!!
//		Matrix3thGen(LatRange,constrained, unconstrained,renorm,LatRange);
	} else if (freedom == secondGeneration) {
		Matrix2ndGen(LatRange,constrained, unconstrained);
	} else Matrix1();
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
}
void
SF_Monomer1stO::GetLayerAnalysis(Output*) {
}
void
SF_Monomer1stO::Matrix1() {
	int z;
	Vector G = Seg->GetSWF();
	Vector phi = GetPhi(total);
	Lat->MakeSafe(G);
	for (z=1; z<=M; z++) {
		phi[z] = G[z];
	}
	lnGN = Lat->ComputeLnGN(G);
	if (freedom == fixedTheta) {
		lnCt = log(theta) - lnGN;
		Lat->NormPhiRestr(phi,G,theta);
		phiBulk = theta/exp(lnGN);
		SetPhiBulk(phiBulk);
	} else if (freedom == rangeRestricted) {
		double phiPreNorm = 0;
		for (z=1; z<=M; z++) {
			if (restrictedRange->InRange(z)) {
				phiPreNorm += phi[z];
			}
		}
		lnCt = log(phiRange/phiPreNorm);
		Lat->NormPhiFree(phi,exp(lnCt));
		theta = ComputeTheta();
		phiBulk = theta/exp(lnGN);
		SetPhiBulk(phiBulk);
	} else {
		lnCb = log(phiBulk);
		Lat->NormPhiFree(phi,phiBulk);
	}
	Lat->RestoreFromSafe(G);
	Lat->RestoreFromSafe(phi);
}
void
SF_Monomer1stO::Matrix2ndGen (const LatticeRange* Range,
							  const DensityPart DensPart1,
							  const DensityPart DensPart2) {
	Vector G1(1,M);
	Vector G2(1,M);
	Vector phi1 = GetPhi(DensPart1);
	Vector phi2 = GetPhi(DensPart2);
	Vector G = Seg->GetSWF();
	Lat->MakeSafe(G);
	Lat->Init2G(G1,G2,G,Range);
	for (int z=1; z<=M; z++) {
		phi1[z] = G1[z];
		phi2[z] = G2[z];
	}
	lnGN = Lat->ComputeLnGN(G1);
	if (muFixed) {
		Lat->NormPhiFree(phi1,exp(lnCt));
		Lat->SubtractBoundaries(phi1);
		Lat->MultiplyWithLatticeSites(phi1);
		theta = 0;
		for (int z=1; z<=M; z++) {
			theta += phi1[z];
		}
		Lat->DivideByLatticeSites(phi1);
		Lat->RestoreBoundaries(phi1);
	} else {
		lnCt = log(theta)-lnGN-log(1.*N);
		Lat->NormPhiRestr(phi1,G1,theta);
	}
	lnCb = log(phiBulk);
	Lat->NormPhiFree(phi2,phiBulk);
	Lat->RestoreFromSafe(G);
	Lat->RestoreFromSafe(phi1);
	Lat->RestoreFromSafe(phi2);
}
/*
void
SF_Monomer1stO::Matrix3thGen (const LatticeRange* Range,
							  const DensityPart DensPart1,
							  const DensityPart DensPart2,
							  const DensityPart DensPart3,
							  const LatticeRange* RangeRenorm,
							  const DensityPart DensPartRenorm) { //broken!!
	int z;
	Vector G1(1,M);
	Vector G2(1,M);
	Vector G = Seg->GetSWF();
	Lat->MakeSafe(G);
	Lat->Init2G(G1,G2,G,Range);
	Vector phi1 = GetPhi(DensPart1);
	Vector phi2 = GetPhi(DensPart2);
	for (z=1; z<=M; z++) {
		phi1[z] = G1[z];
		phi2[z] = G2[z];
	}
	GN = Lat->ComputeGN(G1);
	Ct = GN/(theta);
	Cb = phiBulk;
	Lat->NormPhiRestr(phi1,G1,theta);
	Lat->NormPhiFree(phi2,Cb);
	Lat->Init2G(G1,G2,G,RangeRenorm);
	Vector phi3 = GetPhi(DensPartRenorm);
	for (z=1; z<=M; z++) {
		phi3[z] = G1[z];
	}
	RenormQ[1] = Lat->RenormPhi(phi1,G1,phi3,thetaRenorm);
	Lat->RestoreFromSafe(G);
	Lat->RestoreFromSafe(phi3);
}*/
void
SF_Monomer1stO::MatrixBulk(const DensityPart DensPart) {
	Vector G = Seg->GetSWF();
	Vector phi = GetPhi(DensPart);
	for (int z=1; z<=M; z++) {
		if (Lat->BulkBoundary(z)) {
			phi[z] = phiBulk*G[z];
		}
	}
}
void
SF_Monomer1stO::CopyBulkBoundaries(const DensityPart DensPart) {
	Vector phi = GetPhi(DensPart);
	Vector PhiBulk = GetPhi(bulk);
	for (int z=1; z<=M; z++) {
		if (Lat->BulkBoundary(z)) {
			phi[z] = PhiBulk[z];
		}
	}
}
