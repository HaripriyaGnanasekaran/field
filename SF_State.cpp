#include "SF_State.h"


SF_State::SF_State(Text name_,
				   double valence_,
				   double epsilon_,
				   double phiBulk_,
				   double alphaBulk_,
				   Boolean internalFreeEnergyGiven_,
				   double internalDegeneration_,
				   double internalEnergy_,
				   SegmentFreedom freedom_,
				   LatticeRange* LatRange_,
				   Lattice* Lat_) {
	name = name_;
	valence = valence_;
	epsilon = epsilon_;
	freedom = freedom_;
	LatRange = LatRange_;
	Lat = Lat_;
	phiBulkFixed = phiBulk_;
	if (freedom == frozen && phiBulkFixed > 0) {
		Message(fatal, "Cannot specify a phibulk for " + STATE + " : " + name
		+ " when freedom is set to frozen");
	}
	if (phiBulkFixed > 0) BoolPhiBulkFixed = true;
	else BoolPhiBulkFixed = false;
	alphaBulkFixed = alphaBulk_;
	if (alphaBulkFixed > 0) BoolAlphaBulkFixed = true;
	else BoolAlphaBulkFixed = false;
	alphaBulk = 1;
	internalFreeEnergyGiven = internalFreeEnergyGiven_;
	internalFreeEnergySet = false;
	if (!internalFreeEnergyGiven) {
		internalFreeEnergy = 0;
	} else {
		internalDegeneration = internalDegeneration_;
		internalEnergy = internalEnergy_;
	}
	phiBulk = 0;
	equilibrium = false;
	extPot = false;
	swf.Dim(1,Lat->GetTotalNumLayers());
	rho.Dim(1,Lat->GetTotalNumLayers());
	Rho_up_to_date = false;
	//LD=false;
	//MC=false;

}
SF_State::SF_State() {
}
SF_State::~SF_State() {
}
Text
SF_State::GetName() const {
	return name;
}
void
SF_State::GetOutput(Output* Out) const {
	Out->PutReal(STATE,name,"valence",valence);
	Out->PutReal(STATE,name,"epsilon",epsilon);
	Out->PutReal(STATE,name,"alphabulk",alphaBulk);
	Out->PutReal(STATE,name,"int.free energy",internalFreeEnergy);
	if (internalFreeEnergyGiven) {
		Out->PutReal(STATE,name,"int.energy",internalEnergy);
		Out->PutReal(STATE,name,"int.degeneration",internalDegeneration);
	}
	Out->PutReal(STATE,name,"phibulk",phiBulk);
	int M = Lat->GetTotalNumLayers();
	double theta = 0;
	Vector phi = GetPhi();
	Lat->SubtractBoundaries(phi);
	Lat->MultiplyWithLatticeSites(phi);
	int z;
	for (z=1; z<=M; z++) {
		theta += phi[z];
	}
	Out->PutReal(STATE,name,"theta",theta);
//  FL wrong theta excess calculation;
//	theta = 0;
//	phi = GetPhi();
//	Lat->SubtractBoundaries(phi);
//	Lat->MultiplyWithLatticeSites(phi);
//	for (z=1; z<=M; z++) {
//		if (phi[z] > 0) {
//			theta += phi[z] - phiBulk;
//		}
//	}
	theta -=Lat->GetTotalNumLatticeSites()*phiBulk;
	Out->PutReal(STATE,name,"theta excess",theta);
	Out->PutProfile(STATE,name,"phi",GetPhi());
	if (freedom != frozen && equilibrium) {
		Out->PutProfile(STATE,name,"G",GetSWF());
	}
	if (freedom != frozen && !equilibrium) {
		SF_MolState* MolState;
		for (int i=1; i<=MolStateQ.Cardinal(); i++) {
			MolState = (SF_MolState*) MolStateQ[i];
			Out->PutProfile("mol",
							MolState->GetMolName(),
							"G " + name,
							MolState->GetSWF());
		}
	}
}
double
SF_State::GetValence() const {
	return valence;
}
double
SF_State::GetEpsilon() const {
	return epsilon;
}
Boolean
SF_State::InternalFreeEnergyGiven() const {
	return internalFreeEnergyGiven;
}
void
SF_State::UpdateInternalFreeEnergy(double roomTemp, double temp) {
	if (!internalFreeEnergyGiven) return;
	double value = internalEnergy*roomTemp/temp-log(internalDegeneration);
	SetInternalFreeEnergy(value);
}
void
SF_State::SetInternalFreeEnergy(double value) {
	internalFreeEnergy = value;
	internalFreeEnergySet = true;
	SF_MolState* MolState;
	for (int i=1; i<=MolStateQ.Cardinal(); i++) {
		MolState = (SF_MolState*) MolStateQ[i];
		MolState->SetInternalFreeEnergy(value);
	}
}
Boolean
SF_State::InternalFreeEnergySet() const {
	return internalFreeEnergySet;
}
double
SF_State::GetInternalFreeEnergy() const {
	return internalFreeEnergy;
}
Boolean
SF_State::AlphaBulkFixed() const {
	return BoolAlphaBulkFixed;
}
double
SF_State::GetAlphaBulk() const {
	return alphaBulk;
}
double
SF_State::GetAlphaBulkFixed() const {
	return alphaBulkFixed;
}
void
SF_State::SetAlphaBulk(double value) {
	SF_MolState* MolState;
	alphaBulk = value;
	for (int i=1; i<=MolStateQ.Cardinal(); i++) {
		MolState = (SF_MolState*) MolStateQ[i];
		MolState->SetAlphaBulk(value);
	}
}
SegmentFreedom
SF_State::GetFreedom() const {
	return freedom;
}
Boolean
SF_State::PhiBulkFixed() const {
	return BoolPhiBulkFixed;
}
double
SF_State::GetPhiBulk() const {
	return phiBulk;
}
double
SF_State::GetPhiBulkFixed() const {
	return phiBulkFixed;
}
void
SF_State::UpdatePhiBulk() {
	int i;
	phiBulk = 0;
	SF_MolState* MolState;
	for (i=1; i<=MolStateQ.Cardinal(); i++) {
		MolState = (SF_MolState*) MolStateQ[i];
		phiBulk += MolState->GetPhiBulk();
	}
}
double
SF_State::GetTheta() const {
	int M = Lat->GetTotalNumLayers();
	Vector phi = GetPhi();
	Lat->SubtractBoundaries(phi);
	Lat->MultiplyWithLatticeSites(phi);
	double theta = 0;
	for (int z=1; z<=M; z++) {
		theta += phi[z];
	}
	return theta;
}
LatticeRange*
SF_State::GetLatRange() const {
	return LatRange;
}
Vector
SF_State::GetSWF() const {
	return swf;
}
Vector
SF_State::GetRHO() const {
	return rho;
}

double
SF_State::GetRHO(int z) const{
	return rho[z];
}

void
SF_State::PutRHO(){
	Rho_up_to_date =true;
	Lat->Sides(rho,GetPhi());
}

void
SF_State::PutRHO(double A, int z) {
	Rho_up_to_date =true;
	rho[z]=A;
}
Boolean
SF_State::RhoSet(){
	return Rho_up_to_date;
}

void
SF_State::ResetRhoSet(){
	Rho_up_to_date = false;
}

void
SF_State::SetSWF(Vector swf_) {
	swfFull = swf_;
	Rho_up_to_date = false;
	SF_MolState* MolState;
	for (int i=1; i<=MolStateQ.Cardinal(); i++) {
		MolState = (SF_MolState*) MolStateQ[i];
		MolState->SetSWF(swfFull);
	}
	if (freedom != loose) {
		swf.Dim(1,Lat->GetTotalNumLayers());
		rho.Dim(1,Lat->GetTotalNumLayers());
		Rho_up_to_date = false;
	} else {
		swf = swfFull;
	}
	equilibrium = true;
}

void
SF_State::DelPos(int r, int SpotType, int SpotSize) {
	if (LatRange->InRange(r)) LatRange->SetPos(r,SpotType,SpotSize,0.0);
}

void
SF_State::ClearAllPos() {
	LatRange->ClearAllPos();
}

void
SF_State::UpdatePos(double x, double y, double z, double *submask) {
	LatRange->SetPos(x,y,z,submask);
}

void
SF_State::UpdatePos(int r, int SpotType, int SpotSize) {
	double waarde = 1.0;
	LatRange->SetPos(r,SpotType,SpotSize,waarde);
}

void
SF_State::UpdateSWF() {
	if (!equilibrium) return;
	Rho_up_to_date = false;
	if (freedom != loose) {
		LatRange->MapVector(swf,swfFull);
	}
}
void
SF_State::SetExternalPotential(Vector pot) {
	externalPotential = pot;
	extPot = true;
}
Vector
SF_State::GetExternalPotential(void) {
	return externalPotential;
}
Boolean
SF_State::ExternalPotential(void) {
	return extPot;
}
Vector
SF_State::GetPhi() const {
	int M = Lat->GetTotalNumLayers();
	SF_MolState* MolState;
	Vector molPhi;
	Vector phi(1,M);
	Array<DensityPart> DensityPartQ;
	for (int i=1; i<=MolStateQ.Cardinal(); i++) {
		MolState = (SF_MolState*) MolStateQ[i];
		DensityPartQ = MolState->GetDensityPartQ();
		int numDensPart = DensityPartQ.Upperbound() - DensityPartQ.Lowerbound()+1;
		for (int j=1; j<=numDensPart; j++) {
			molPhi = MolState->GetPhi(DensityPartQ[j]);
			for (int z=1; z<=M; z++) {
				phi[z] += molPhi[z];
			}
		}
	}
	if (freedom == frozen) {
		for (int z=1; z<=M; z++) {
			// if (LatRange->InRange(z)) phi[z] = 1;  else phi[z] = 0;
			phi[z]=LatRange->GetRangeValue(z);
		}
		Lat->SetBoundaries(phi);
	}
	return phi;
}
void
SF_State::AddMolState(SF_MolState* MolState) {
	MolState->Into(MolStateQ);
}
int
SF_State::GetNumMolStates() const {
	return MolStateQ.Cardinal();
}
SF_MolState*
SF_State::GetMolState(int number) const {
	return (SF_MolState*) MolStateQ[number];
}
