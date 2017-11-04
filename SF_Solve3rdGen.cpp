#include "SF_Solve3rdGen.h"

SF_Solve3rdGen::SF_Solve3rdGen(Boolean compute_,
							   SF_ReactionList* ReactionQ_,
							   SF_MoleculeList* MolQ_,
							   SF_SegmentList* SegQ_,
							   Lattice* Lat_,
							   SF_ReactionList* ReactionQRN_,
							   SF_MoleculeList* MolQRN_,
							   SF_SegmentList* SegQRN_,
							   Lattice* LatRN_,
							   Input* MyInput_) :
		SF_Solve(compute_,ReactionQ_,MolQ_,SegQ_,Lat_,MyInput_) {
	
	ReactionQRN = ReactionQRN_;
	MolQRN = MolQRN_;
	SegQRN = SegQRN_;
	LatRN = LatRN_;
//	Boolean equal;
	StateQRN.Dim(1,numLiquidStates);
	int i;
	for (i=1; i<=numLiquidStates; i++) {
		StateQRN[i] = SegQRN->GetState((StateQ[i])->GetName());
	}
	SWFRN.Dim(1,numLiquidStates);
	MRN = LatRN->GetTotalNumLayers();
	for (i=1; i<=numLiquidStates; i++) {
		SWFRN[i].Dim(1,MRN);
		(StateQRN[i])->SetSWF(SWFRN[i]);
	}
	SegQRN->UpdateSWF();
	int z;
	int extraItVar = 0;
	for (z=1; z<=MRN; z++) {
		if (SegItVar(LatRN,z)) {
			extraItVar++;
		}
	}
	numItVar += extraItVar*numDiffItVar;
	if (charged) {
		psi0RN.Dim(1,MRN);
		for (z=1; z<=MRN; z++) {
			if (PotItVar(LatRN,z)) {
				numItVar++;
			}
		}
	}
	x.Dim(1,numItVar);
	SF_Molecule* Mol;
	MolNum3rdGen = 0;
	for (i=1; i<=MolQ->GetNumMolecules(); i++) {
		Mol = MolQ->GetMolecule(i);
		if (Mol->GetFreedom() == thirdGeneration) {
			if (MolNum3rdGen != 0) {
				Message(fatal, "define only one molecule with freedom"
					" 'third_generation'");
			} else {
				MolNum3rdGen = i;
			}
		}	
	}
}
SF_Solve3rdGen::~SF_Solve3rdGen() {
}
void
SF_Solve3rdGen::Iterate() {
	ProcessInput();
	iterate(&x[1],numItVar);
	if (getaccuracy() != fabs(getaccuracy())) { //true if accuracy == NaN
		WriteOverflowError();
	}
}
void
SF_Solve3rdGen::GetResiduals(Vector f) const {
	int i,z;
	Vector phiTotal = SegQ->GetPhiTotal();
	Vector uPrimeAv(1,M);
	double uPrime;
	Vector psi;
	if (charged) {
		if (potentialItVar) {
			psi = SegQ->ComputeElectricPotential(psi0);
		} else {
			psi = SegQ->ComputeElectricPotential();
		}
	}
	SF_State* State;
	for (z=1; z<=M; z++) {
		if (SegItVar(Lat,z)) {
			int count=0;
			for (i=1; i<=numDiffItVar; i++) {
				count+=NumStatesForItVar[i];
				State = StateQ[count];
				uPrime = -log(SWF[count][z]);
				uPrime -= SegQ->ChemInt(State,z)/phiTotal[z];
				uPrime += SegQ->ChemIntBulk(State);
				if (charged) {
					uPrime -= State->GetValence()*psi0[z]*
						ELEM_CHARGE/(BOLTZMANN*TEMPERATURE);
				}
				uPrimeAv[z] += uPrime;
			}
			uPrimeAv[z] /= numDiffItVar;
		}
	}
	int j=1;
	for (z=1; z<=M; z++) {
		if (PotItVar(Lat,z)) {
			f[j++] = psi0[z] - psi[z];
		}
		if (SegItVar(Lat,z)) {
			int count=0;
			for (i=1; i<=numDiffItVar; i++) {
				count+=NumStatesForItVar[i];
				State = StateQ[count];
				uPrime = -log(SWF[count][z]);
				uPrime -= SegQ->ChemInt(State,z)/phiTotal[z];
				uPrime += SegQ->ChemIntBulk(State);
				if (charged) {
					uPrime -= State->GetValence()*psi0[z]*
						ELEM_CHARGE/(BOLTZMANN*TEMPERATURE);
				}
				f[j++] = - 1 + 1/phiTotal[z] + uPrime - uPrimeAv[z];
			}
		}
	}
	phiTotal = SegQRN->GetPhiTotal();
	uPrimeAv.Dim(1,MRN);
	if (charged) {
		if (potentialItVar) {
			psi = SegQRN->ComputeElectricPotential(psi0RN);
		} else {
			psi = SegQRN->ComputeElectricPotential();
		}
	}
	for (z=1; z<=MRN; z++) {
		if (SegItVar(LatRN,z)) {
			int count=0;
			for (i=1; i<=numDiffItVar; i++) {
				count+=NumStatesForItVar[i];
				State = StateQRN[count];
				uPrime = -log(SWFRN[count][z]);
				uPrime -= SegQRN->ChemInt(State,z)/phiTotal[z];
				uPrime += SegQRN->ChemIntBulk(State);
				if (charged) {
					uPrime -= State->GetValence()*psi0RN[z]*
						ELEM_CHARGE/(BOLTZMANN*TEMPERATURE);
				}
				uPrimeAv[z] += uPrime;
			}
			uPrimeAv[z] /= numDiffItVar;
		}
	}
	for (z=1; z<=MRN; z++) {
		if (PotItVar(LatRN,z)) {
			f[j++] = psi0RN[z] - psi[z];
		}
		if (SegItVar(LatRN,z)) {
			int count=0;
			for (i=1; i<=numDiffItVar; i++) {
				count+=NumStatesForItVar[i];
				State = StateQRN[count];
				uPrime = -log(SWF[count][z]);
				uPrime -= SegQRN->ChemInt(State,z)/phiTotal[z];
				uPrime += SegQRN->ChemIntBulk(State);
				if (charged) {
					uPrime -= State->GetValence()*psi0RN[z]*
						ELEM_CHARGE/(BOLTZMANN*TEMPERATURE);
				}
				f[j++] = - 1 + 1/phiTotal[z] + uPrime - uPrimeAv[z];
			}
		}
	}
}
void 
SF_Solve3rdGen::residuals(double *const ff, double *const) {
	Vector f(ff, 1, numItVar);
	ComputePhi();
	GetResiduals(f);
}
void
SF_Solve3rdGen::inneriteration(double *const funct,
						 double *const itvar,
						 double acc) {
	SF_Solve::inneriteration(funct,itvar,acc);						
	if (multiState) {
		ReactionQRN->InnerIterate();
		if (ReactionQ->getiterations() > 0) {
			residuals(funct,itvar);
		}
	}
}
void
SF_Solve3rdGen::ComputePhi() {
	UpdateSWF();
	if (multiState && neutralizerPresent) {
		MolQ->ComputePhi(phiBulkSolvent,phiBulkNeutralizer);
		MolQRN->ComputePhi(phiBulkSolvent,phiBulkNeutralizer);
	} else if (multiState) {
		MolQ->ComputePhi(phiBulkSolvent);
		MolQRN->ComputePhi(phiBulkSolvent);
	} else {
		MolQ->ComputePhi();
		MolQRN->ComputePhi();
	}
	SF_Molecule* Mol = MolQ->GetMolecule(MolNum3rdGen);
	SF_Molecule* MolRN = MolQRN->GetMolecule(MolNum3rdGen);
	Vector Renorm = Lat->RenormPhi(Mol->GetPhi(renorm),
								   MolRN->GetPhi(constrained),
								   Mol->GetPhi(constrained),
								   Mol->GetThetaRenorm());
	Mol->SetRenorm(Renorm,1);
	MolRN->SetRenorm(Renorm,1);
}	
void
SF_Solve3rdGen::UpdateSWF() {
	int i,j,k,z;
	j=1;
	for (z=1; z<=M; z++) {
		if (PotItVar(Lat,z)) {
			psi0[z] = x[j++];
		}
		if (SegItVar(Lat,z)) {
			int count = 0;
			SF_State* State;
			for (i=1; i<=numDiffItVar; i++) {
				for (k=1; k<=NumStatesForItVar[i]; k++) {
					count++;
					SWF[count][z] = exp(-x[j]);
					if (charged) {
						State = StateQ[count];
						SWF[count][z] *= exp(-State->GetValence()*psi0[z]*
							ELEM_CHARGE/(BOLTZMANN*TEMPERATURE));
					}
				}
				j++;
			}
		}
	}
	int count = 0;
	for (i=1; i<=numDiffItVar; i++) {
		for (k=1; k<=NumStatesForItVar[i]; k++) {
			count++;
			Lat->SetBoundaries(SWF[count]);
		}
	}
	for (z=1; z<=MRN; z++) {
		if (PotItVar(LatRN,z)) {
			psi0RN[z] = x[j++];
		}
		if (SegItVar(LatRN,z)) {
			int count = 0;
			SF_State* State;
			for (i=1; i<=numDiffItVar; i++) {
				for (k=1; k<=NumStatesForItVar[i]; k++) {
					count++;
					SWFRN[count][z] = exp(-x[j]);
					if (charged) {
						State = StateQRN[count];
						SWFRN[count][z] *= exp(-State->GetValence()*psi0RN[z]*
							ELEM_CHARGE/(BOLTZMANN*TEMPERATURE));
					}
				}
				j++;
			}
		}
	}
	count = 0;
	for (i=1; i<=numDiffItVar; i++) {
		for (k=1; k<=NumStatesForItVar[i]; k++) {
			count++;
			LatRN->SetBoundaries(SWFRN[count]);
		}
	}
	if (neutralizerPresent) {
		phiBulkNeutralizer = (MolQ->GetNeutralizer())->GetPhiBulk();
	}
	phiBulkSolvent = (MolQ->GetSolvent())->GetPhiBulk();
	if (charged) {
		Lat->SetBoundaries(psi0);
		LatRN->SetBoundaries(psi0RN);
	}
	SegQ->UpdatePhiBulk();
	SegQ->UpdateSWF();
	SegQRN->UpdatePhiBulk();
	SegQRN->UpdateSWF();
}
void
SF_Solve3rdGen::UpdateItVar(Boolean offset) {
	int i,j,z;
	SF_State* State;
	j=1;
	for (z=1; z<=M; z++) {
		if (PotItVar(Lat,z)) {
			x[j++] = psi0[z];
		}
		if (SegItVar(Lat,z)) {
			int count=0;
			for (i=1; i<=numDiffItVar; i++) {
				count += NumStatesForItVar[i];
				State = StateQ[count];
				if (SWF[count][z] > 0) {
					x[j] = -log(SWF[count][z]);
					if (charged) {
						x[j] -= State->GetValence()*psi0[z]*
							ELEM_CHARGE/(BOLTZMANN*TEMPERATURE);
					}
				}
				j++;
			}
		}
	}
	for (z=1; z<=MRN; z++) {
		if (PotItVar(Lat,z)) {
			x[j++] = psi0RN[z];
		}
		if (SegItVar(LatRN,z)) {
			int count=0;
			for (i=1; i<=numDiffItVar; i++) {
				count += NumStatesForItVar[i];
				State = StateQRN[count];
				if (SWFRN[count][z] > 0) {
					x[j] = -log(SWFRN[count][z]);
					if (charged) {
						x[j] -= State->GetValence()*psi0RN[z]*
							ELEM_CHARGE/(BOLTZMANN*TEMPERATURE);
					}
				}
				j++;
			}
		}
	}
	if (neutralizerPresent) {
		phiBulkNeutralizer = (MolQ->GetNeutralizer())->GetPhiBulk();
	}
	phiBulkSolvent = (MolQ->GetSolvent())->GetPhiBulk();
}
Boolean
SF_Solve3rdGen::SegItVar(const Lattice* SomeLat, const int z) const {
	if (!SomeLat->WithinBoundaries(z)) return false;
	SF_State* State;
	int i;
	LatticeRange* LatRange;
	if (SomeLat->GetNumLayers(1) == 3) {
		int numStates = SegQRN->GetNumStates();
		Boolean itVar = true;
		for (i=1; i<=numStates; i++) {
			State = SegQRN->GetState(i);
			if (State->GetFreedom() == frozen) {
				LatRange = State->GetLatRange();
				if (LatRange->InRange(z)) itVar = false;
			}
		}
		return itVar;
	} else {
		int numStates = SegQ->GetNumStates();
		Boolean itVar = true;
		for (i=1; i<=numStates; i++) {
			State = SegQ->GetState(i);
			if (State->GetFreedom() == frozen) {
				LatRange = State->GetLatRange();
				if (LatRange->InRange(z)) itVar = false;
			}
		}
		return itVar;
	}
}
Boolean
SF_Solve3rdGen::PotItVar(const Lattice* SomeLat, const int z) const {
	if (!charged) return false;
	if (SomeLat->WithinBoundaries(z)) {
		return potentialItVar;
	}
	else return false;
}

