#include "SF_MolState.h"

Text STATE = "state";
Text SEGMENT = "mon";

SF_MolState::SF_MolState(Text name_,
						Text molName_,
						double valence_,
						double epsilon_,
						double alphaBulk_,
						SegmentFreedom freedom_,
						LatticeRange* LatRange_,
						Lattice* Lat_) {
	name = name_;
	molName = molName_;
	valence = valence_;
	epsilon = epsilon_;
	alphaBulk = alphaBulk_;
	freedom = freedom_;
	LatRange = LatRange_;
	Lat = Lat_;
	alphaBulkBoundaries.Dim(1,2*Lat->GetNumGradients());
	for (int i=1; i<=2*Lat->GetNumGradients(); i++) {
		alphaBulkBoundaries[i] = alphaBulk;
	}
	numPhi = 0;
	phiBulk = 0;
	phiRef = 0;
	extPot = false;
	internalFreeEnergy = 0;
}
SF_MolState::~SF_MolState() {
	Out();
}
Text
SF_MolState::GetName() const {
	return name;
}
Text
SF_MolState::GetMolName() const {
	return molName;
}
double
SF_MolState::GetValence() const {
	return valence;
}
double
SF_MolState::GetEpsilon() const {
	return epsilon;
}
double
SF_MolState::GetInternalFreeEnergy() const {
	return internalFreeEnergy;
}
void
SF_MolState::SetInternalFreeEnergy(double value) {
	internalFreeEnergy = value;
}
double
SF_MolState::GetAlphaBulk() const {
	return alphaBulk;
}
void
SF_MolState::SetAlphaBulk(double value) {
	alphaBulk = value;
}
Vector
SF_MolState::GetAlphaBulkBoundaries() const {
	return alphaBulkBoundaries;
}
void
SF_MolState::SetAlphaBulkBoundaries(Vector value) {
	alphaBulkBoundaries = value;
}
double
SF_MolState::GetPhiBulk() const {
	return phiBulk;
}
void
SF_MolState::SetPhiBulk(double phiBulk_) {
	phiBulk = phiBulk_*alphaBulk;
}
Vector
SF_MolState::GetPhiBulkBoundaries() const {
	return phiBulkBoundaries;
}
void
SF_MolState::SetPhiBulkBoundaries(Vector phiBulkBoundaries_) {
	int numBounds = Lat->GetNumGradients()*2;
	phiBulkBoundaries.Dim(1,numBounds);
	for (int i=1; i<=numBounds; i++) {
		phiBulkBoundaries[i] = phiBulkBoundaries_[i]*alphaBulkBoundaries[i];
	}
}
double
SF_MolState::GetPhiRef() const {
	return phiRef;
}
void
SF_MolState::SetPhiRef(double phiRef_) {
	phiRef = phiRef_;
}
SegmentFreedom
SF_MolState::GetFreedom() const {
	return freedom;
}
LatticeRange*
SF_MolState::GetLatRange() const {
	return LatRange;
}
void
SF_MolState::SetSWF(Vector swf_) {
	swfFull = swf_;
	if (freedom != loose) {
		swf.Dim(1,Lat->GetTotalNumLayers());
	} else {
		swf = swfFull;
	}
}
Vector
SF_MolState::GetSWF() const {
	return swf;
}

void
SF_MolState::DelPos(int r,int SpotType, int SpotSize) {
	if (LatRange->InRange(r)) LatRange->SetPos(r,SpotType,SpotSize,0.0);
}

void
SF_MolState::ClearAllPos() {
	LatRange->ClearAllPos();
}

void
SF_MolState::UpdatePos(int r,int SpotType, int SpotSize) {
double waarde = 1.0;
	LatRange->SetPos(r,SpotType,SpotSize,waarde);
}

void
SF_MolState::UpdatePos(double x, double y, double z, double *submask) {
	LatRange->SetPos(x,y,z,submask);
}

void
SF_MolState::UpdateSWF() {
	if (freedom != loose) {
		LatRange->MapVector(swf,swfFull);
	}
}


void
SF_MolState::SetExternalPotential(Vector pot) {
	externalPotential = pot;
	extPot = true;
}
Vector
SF_MolState::GetExternalPotential(void) {
	return externalPotential;
}
Boolean
SF_MolState::ExternalPotential(void) {
	return extPot;
}
Vector
SF_MolState::GetPhi(DensityPart densityPart) const {
	int i;
	if (densityPart == total) {
		if (numPhi == 1) {
			//Lat->SubtractBoundaries(PhiQ[1]);
			//Lat->RestoreBoundaries(PhiQ[1]);
			return PhiQ[1];
		}
		int M = Lat->GetTotalNumLayers();
		Vector phiTotal(1,M);
		Vector phi;
		for (i=1; i<=numPhi; i++) {
			phi = PhiQ[i];
			for (int z=1; z<=M; z++) {
				phiTotal[z] += phi[z];

			}
		}
		//Lat->SubtractBoundaries(phiTotal);
		//Lat->RestoreBoundaries(phiTotal);
		return phiTotal;
	}
	for (i=1; i<=numPhi; i++) {
		if (DensityPartQ[i] == densityPart) {
			return PhiQ[i];
		}
	}
	Message(fatal,"Programming error in call to "
		"SF_MolState::GetPhi(DensityPart densityPart), "
		"Create phi before using it");
	return PhiQ[1]; // never get here
}
double
SF_MolState::GetTheta(DensityPart densityPart) const {
	double theta = 0;
	int M = Lat->GetTotalNumLayers();
	Vector phi = GetPhi(densityPart);
	Lat->SubtractBoundaries(phi);
	Lat->MultiplyWithLatticeSites(phi);
	for (int z=1; z<=M; z++) {
		theta += phi[z];
	}
	Lat->DivideByLatticeSites(phi);
	Lat->RestoreBoundaries(phi);
	return theta;
}

// allocate the Phi[z] vectors
void
SF_MolState::CreatePhi(DensityPart densityPart) {
	if (numPhi == 0) {
		PhiQ.Dim(1,1);
		DensityPartQ.Dim(1,1);
		Vector phi(1,Lat->GetTotalNumLayers());
		PhiQ[1] = phi;
		DensityPartQ[1] = densityPart;
		numPhi++;
	} else {
		Array<Vector> PhiQtemp(1,numPhi+1);
		Array<DensityPart> DensityPartQtemp(1,numPhi+1);
		for (int i=1; i<=numPhi; i++) {
			PhiQtemp[i] = PhiQ[i];
			DensityPartQtemp[i] = DensityPartQ[i];
		}
		Vector phi(1,Lat->GetTotalNumLayers());
		numPhi++;
		DensityPartQtemp[numPhi] = densityPart;
		PhiQtemp[numPhi] = phi;
		PhiQ = PhiQtemp;
		DensityPartQ = DensityPartQtemp;
	}
}
void
SF_MolState::DeletePhi(DensityPart densityPart) {
	int i;
	int number = 0;
	for (i=1; i<=numPhi; i++) {
		if (DensityPartQ[i] == densityPart) number = i;
	}
	if (number == 0) {
		Message(fatal,"Programming error in call "
		"to SF_MolState::DeletePhi(DensityPart densityPart), "
		"Create phi before deleting it");
	}

	Array<Vector> PhiQtemp(1,numPhi-1);
	Array<DensityPart> DensityPartQtemp(1,numPhi-1);
	for (i=1; i<=numPhi; i++) {
		if (i < number) {
			PhiQtemp[i] = PhiQ[i];
			DensityPartQtemp[i] = DensityPartQ[i];
		} else if (i > number) {
			PhiQtemp[i-1] = PhiQ[i];
			DensityPartQtemp[i-1] = DensityPartQ[i];
		}
	}
	numPhi--;
	PhiQ = PhiQtemp;
	DensityPartQ = DensityPartQtemp;
}
Array<DensityPart>
SF_MolState::GetDensityPartQ() const {
	return DensityPartQ;
}
