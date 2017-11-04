#include "SF_Monomer2ndO.h"
#include "tools.h"

SF_Monomer2ndO::SF_Monomer2ndO(Text name_,
								   SF_SegmentList* SegQ_,
								   Lat2ndO* Lat_,
								   Input* MyInput_)
: SF_Molecule(name_, SegQ_, Lat_, MyInput_)
{
	if (Chain->GetMoleculeType() != monomer) {
		Message(fatal,MyInput,
		"Programming error, trying to create a monomer for molecule '" +
		name + "'");
	}
	Lat = Lat_;
	//force_set = (!GetForce()==0);
	Array<Text> sysNames = MyInput->GetNames("sys");
	Text sysName = sysNames[1];
	//MayerSaupe = MyInput->GetBoolean("sys",sysName,"MayerSaupe",false);
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
	gradients = Lat->GetNumGradients();


	if (gradients==3) {size=6;}// CreatePhi(x_dir); CreatePhi(y_dir); CreatePhi(z_dir);}
	if (gradients==2) {size=5;}// CreatePhi(x_dir); CreatePhi(y_dir); CreatePhi(z_dir);}
	if (gradients==1) {size=3;}// CreatePhi(x_dir); CreatePhi(yz_dir);}

	//if (MayerSaupe) {
	//	CreatePhi(parallel_to_director);
	//}
	SetPhiBulk(phiBulk);
	RenormQ.Dim(1,1);
}
SF_Monomer2ndO::~SF_Monomer2ndO() {
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
	//if (gradients ==1) {DeletePhi(x_dir); DeletePhi(yz_dir); } else{
	//	DeletePhi(x_dir); DeletePhi(y_dir); DeletePhi(z_dir);
	//}
	//if (MayerSaupe) {
	//	DeletePhi(parallel_to_director);
	//}
}
void
SF_Monomer2ndO::ReComputePhi() {
	Message(fatal,"ReComputePhi not implemented in SF_Monomer2ndO");
}

Vector
SF_Monomer2ndO::GetLong(){
	Message(fatal,"GetLong not implemented in SF_Branched1stO");
	Vector x;
	return x;
}
Vector
SF_Monomer2ndO::GetShort(){
	Message(fatal,"GetShort not implemented in SF_Branched1stO");
	Vector x;
	return x;
}
Vector
SF_Monomer2ndO::GetBondOrientation( const DensityPart DensPar) {
	Message(fatal,"GetBondOrientation not implemented in SF_Monomer2ndO");
	Vector x;
	return x; 
};


MoleculeType
SF_Monomer2ndO::GetMoleculeType() const {
	return Chain->GetMoleculeType();
}


void
SF_Monomer2ndO::ComputePhi() {
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
SF_Monomer2ndO::GetLayerAnalysis(Output*) {
}
void
SF_Monomer2ndO::Matrix1() {
	int z,k;

	Vector Gi(1,M*size);
	Vector G = Seg->GetSWF();
	Vector phi = GetPhi(total);
	Lat->MakeSafe(G);

	cp(phi,G,M);
	for (k=0; k<size; k++) {double *gi = &Gi[1+k*M], *g = &G[1]; cp(gi,g,M);}
	lnGN = Lat->ComputeLnGN(Gi);
	if (freedom == fixedTheta) {
		lnCt = log(theta) - lnGN;
		Lat->NormPhiRestr(phi,Gi,theta*6);
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
SF_Monomer2ndO::Matrix2ndGen (const LatticeRange* Range,
							  const DensityPart DensPart1,
							  const DensityPart DensPart2) {
     Message(fatal,"Matrix2ndGen not implemented in SF_Monomer2ndO");
}

void
SF_Monomer2ndO::MatrixBulk(const DensityPart DensPart) {
	Vector G = Seg->GetSWF();
	Vector phi = GetPhi(DensPart);
	for (int z=1; z<=M; z++) {
		if (Lat->BulkBoundary(z)) {
			phi[z] = phiBulk*G[z];
		}
	}
}
void
SF_Monomer2ndO::CopyBulkBoundaries(const DensityPart DensPart) {
	Vector phi = GetPhi(DensPart);
	Vector PhiBulk = GetPhi(bulk);
	for (int z=1; z<=M; z++) {
		if (Lat->BulkBoundary(z)) {
			phi[z] = PhiBulk[z];
		}
	}
}
