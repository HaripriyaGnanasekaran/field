#include "SF_SolveFraaije.h"

SF_SolveFraaije::SF_SolveFraaije(Boolean compute_,
				   SF_ReactionList* ReactionQ_,
				   SF_MoleculeList* MolQ_,
				   SF_SegmentList* SegQ_,
				   SF_ReactionList* ReactionNewQ_,
				   SF_MoleculeList* MolNewQ_,
				   SF_SegmentList* SegNewQ_,
				   Lattice* Lat_,
				   Input* MyInput_) {
	
	firstIterate = true;
	compute = compute_;
	ReactionQ = ReactionQ_;
	MolQ = MolQ_;
	SegQ = SegQ_;
	ReactionNewQ = ReactionNewQ_;
	MolNewQ = MolNewQ_;
	SegNewQ = SegNewQ_;
	Lat = Lat_;
	MyInput = MyInput_;
	ProcessInput();
	charged = SegQ->Charged();
	if (MyInput->GetNumNames("reaction") > 0) {
		multiState = true;
	} else {
		multiState = false;
	}
	smallAlphaCount = 0;
	reverseDirection.Dim(1,reverseDirectionRange);
	numLiquidStates = 0;
	SF_MolStructure* Chain;
	SF_MolSegment* Seg;
	int i,j;
	for (i=1; i<=MolQ->GetNumMolecules(); i++) {
		Chain = (MolQ->GetMolecule(i))->GetMolStructure();
		for (j=1; j<=Chain->GetNumDiffSegments(); j++) {
			Seg = Chain->GetDiffSegment(j);
			numLiquidStates+=Seg->GetNumStates();
		}
	}
	StateQ.Dim(1,numLiquidStates);
	StateNewQ.Dim(1,numLiquidStates);
	numLiquidStates = 0;
	SF_MolStructure* ChainNew;
	SF_MolSegment* SegNew;
	for (i=1; i<=MolQ->GetNumMolecules(); i++) {
		Chain = (MolQ->GetMolecule(i))->GetMolStructure();
		ChainNew = (MolNewQ->GetMolecule(i))->GetMolStructure();
		for (j=1; j<=Chain->GetNumDiffSegments(); j++) {
			Seg = Chain->GetDiffSegment(j);
			SegNew = ChainNew->GetDiffSegment(j);
			for (int k=1; k<=Seg->GetNumStates(); k++) {
				numLiquidStates++;
				StateQ[numLiquidStates] = Seg->GetState(k);
				StateNewQ[numLiquidStates] = SegNew->GetState(k);
			}
		}
	}
	int z;
	M = Lat->GetTotalNumLayers();
	SWF.Dim(1,numLiquidStates);
	SWFNew.Dim(1,numLiquidStates);
	dPhi.Dim(1,numLiquidStates);
	dPhiNew.Dim(1,numLiquidStates);
	for (i=1; i<=numLiquidStates; i++) {
		SWF[i].Dim(1,M);
		(StateQ[i])->SetSWF(SWF[i]);
		SWFNew[i].Dim(1,M);
		(StateNewQ[i])->SetSWF(SWFNew[i]);
		dPhi[i].Dim(1,M);
		dPhiNew[i].Dim(1,M);
		for (z=1; z<=M; z++) {
			SWFNew[i][z] = 1;
			SWF[i][z] = 1;
		}
	}
	SegQ->UpdateSWF();
	SegNewQ->UpdateSWF();
	if (MolQ->GetNeutralizer() == NULL) neutralizerPresent = false;
	else neutralizerPresent = true;
	if (multiState) {
		Message(fatal,"non equilibrium multistates are not implemented");
		if (!neutralizerPresent && ReactionQ->NeutralizerNeeded()) {
			Message(fatal,"Please define a molecule with freedom neutralizer");
		}
	}
	numItVar = 0;
	for (z=1; z<=M; z++) {
		if (SegItVar(z)) {
			numItVar++;
		}
		if (BulkItVar(z)) {
			numItVar++;
		}
	}
	numItVar*=numLiquidStates;
	if (charged) {
		psi0.Dim(1,M);
		potentialItVar = /*SegQ->WeakElectrolytes();*/ true;                                 
		for (z=1; z<=M; z++) {
			if (PotItVar(z)) {
				numItVar++;
			}
		}
	}
	x.Dim(1,numItVar);
	oldNumIterations = 0;
	oldNumFunctions = 0;
}
SF_SolveFraaije::~SF_SolveFraaije() {
}
void
SF_SolveFraaije::GetOutput(Output* Out) const {
	Out->PutText("newton",name,"iteration method","Crank Nicolson");
	Out->PutInt("newton",name,"iterations",getiterations());
	Out->PutInt("newton",name,"iterationlimit",getiterationlimit());
	Out->PutBoolean("newton",name,"error occurred",errorOccurred);
	if (errorOccurred) {
		Out->SetErrorOccurred();
	}
	Out->PutReal("newton",name,"accuracy",getaccuracy());
	Out->PutReal("newton",name,"tolerance",gettolerance());
	Out->PutReal("newton",name,"deltamin",getdeltamin());
	Out->PutReal("newton",name,"deltamax",getdeltamax());
	Out->PutBoolean("newton",name,"pseudohessian",pseudohessian);
	Out->PutInt("newton",name,"number iter var",numItVar);
	Out->PutVector("newton",name,"iter var",x,numItVar);
	Vector funct(1,numItVar);
	GetResiduals(funct);
	Out->PutVector("newton",name,"residuals",funct,numItVar);
}
void 
SF_SolveFraaije::SetInitialGuessFromPrevious(const Lattice* LatOld, 
											 const SF_SegmentList* SegQOld, 
											 const SF_MoleculeList* MolQOld) {
	SF_MolState* StateOld;
	SF_MolState* State;
	SF_MolState* StateNew;
	Array<SF_MolState*> StateQOld;
	Text stateName;
	int gradients = LatOld->GetNumGradients();
	int numLiquidStatesOld = 0;
	SF_MolStructure* Chain;
	SF_MolSegment* Seg;
	int i,j;
	for (i=1; i<=MolQOld->GetNumMolecules(); i++) {
		Chain = (MolQOld->GetMolecule(i))->GetMolStructure();
		for (j=1; j<=Chain->GetNumDiffSegments(); j++) {
			Seg = Chain->GetDiffSegment(j);
			numLiquidStatesOld+=Seg->GetNumStates();
		}
	}
	StateQOld.Dim(1,numLiquidStatesOld);
	numLiquidStatesOld = 0;
	for (i=1; i<=MolQOld->GetNumMolecules(); i++) {
		Chain = (MolQOld->GetMolecule(i))->GetMolStructure();
		for (j=1; j<=Chain->GetNumDiffSegments(); j++) {
			Seg = Chain->GetDiffSegment(j);
			for (int k=1; k<=Seg->GetNumStates(); k++) {
				StateQOld[++numLiquidStatesOld] = Seg->GetState(k);
			}
		}
	}
	for (i=1; i<=numLiquidStates; i++) {
		StateNew = StateNewQ[i];
		State = StateQ[i];
		for (int j=1; j<=numLiquidStatesOld; j++) {
			StateOld = StateQOld[j];
			Boolean found = false;
			if ( *(StateNew->GetName()) == 	*(StateOld->GetName()) && 
				*(StateNew->GetMolName()) == *(StateOld->GetMolName()) ) {
					found = true;
			}
			if (found) {
				if (StateNew->GetFreedom() != frozen 
				&& StateOld->GetFreedom() != frozen) {
					if (gradients == 1) {
						int Min = LatOld->GetTotalNumLayers();
						CopyInitialGuess(StateOld->GetSWF(),StateNew->GetSWF(),Min);
						CopyInitialGuess(StateOld->GetSWF(),State->GetSWF(),Min);
					} else if (gradients == 2) {
						int MXin = LatOld->GetNumLayers(1);
						int MYin = LatOld->GetNumLayers(2);
						CopyInitialGuess(StateOld->GetSWF(),StateNew->GetSWF(),MXin,MYin);
						CopyInitialGuess(StateOld->GetSWF(),State->GetSWF(),MXin,MYin);
					}
				}
			}
		}
	}
	if (gradients == 1 && charged) {
		int Min = LatOld->GetTotalNumLayers();
		CopyInitialGuess(SegQOld->GetElectricPotential(),psi0,Min);
	}
	if (gradients == 2 && charged) {
		int MXin = LatOld->GetNumLayers(1);
		int MYin = LatOld->GetNumLayers(2);
		CopyInitialGuess(SegQOld->GetElectricPotential(),psi0,MXin,MYin);
	}
	(MolQ->GetSolvent())->SetPhiBulk((MolQOld->GetSolvent())->GetPhiBulk());
	phiBulkSolvent = (MolQ->GetSolvent())->GetPhiBulk();
	SegQ->UpdateSWF();
	SegQ->UpdatePhiBulk();
	UpdateItVar(false);
	UpdateSWF();
	if (neutralizerPresent) {
		MolQ->ComputePhi(phiBulkSolvent,phiBulkNeutralizer);
	} else {
		MolQ->ComputePhi(phiBulkSolvent);
	}		
	if (multiState) {
		SegQ->UpdatePhiBulk();
		SegQ->UpdatePhiStates();
	}
}
void
SF_SolveFraaije::ComputePhi() {
	UpdateSWF();
	if (multiState && neutralizerPresent) {
		MolNewQ->ComputePhi(phiBulkSolvent,phiBulkNeutralizer);
	} else if (multiState) {
		MolNewQ->ComputePhi(phiBulkSolvent);
	} else {
	 	MolNewQ->ComputePhi();
	}
	UpdateSWF();
	if (!firstIterate)  {
		CalculateDPhiNew();
	}
}
void
SF_SolveFraaije::GetResiduals(Vector f) const {
	Vector psi;
	if (charged) {
		if (potentialItVar) {
			psi = SegNewQ->ComputeElectricPotential(psi0);
		} else {
			psi = SegNewQ->ComputeElectricPotential();
		}
	}
	int j=1;
	Vector phi,phiNew;
	if (firstIterate)  {
		for (int z=1; z<=M; z++) {
			if (PotItVar(z)) {
				f[j++] = psi0[z] - psi[z];
			}
			if (SegItVar(z) || BulkItVar(z)) {
				for (int i=1; i<=numLiquidStates; i++) {
					phi = (StateQ[i])->GetPhi(total);
					phiNew = (StateNewQ[i])->GetPhi(total);
					f[j++] = phi[z] - phiNew[z];
				}
			}
		}
	} else {
		Vector phiTotal = SegNewQ->GetPhiTotal();
		for (int z=1; z<=M; z++) {
			if (PotItVar(z)) {
				f[j++] = psi0[z] - psi[z];
			/*		Sysout().Outtext("psi0[z] : ");
					Sysout().Outreal(psi0[z],4,0);
					Sysout().Outtext(" psi[z] : ");
					Sysout().Outreal(psi[z],4,0);
					Sysout().Outimage();*/
			}
			if (SegItVar(z)) {
				for (int i=1; i<=numLiquidStates-1; i++) {
					phi = (StateQ[i])->GetPhi(total);
					phiNew = (StateNewQ[i])->GetPhi(total);
					f[j++] = (phi[z] + (dPhi[i][z] + dPhiNew[i][z])*timeStep-phiNew[z])/phi[z];
				/*	Sysout().Outtext("phi[z] : ");
					Sysout().Outreal(phi[z],4,0);
					Sysout().Outtext(" dPhi[i][z] : ");
					Sysout().Outreal(dPhi[i][z],4,0);
					Sysout().Outtext(" dPhiNew[i][z] : ");
					Sysout().Outreal(dPhiNew[i][z],4,0);
					Sysout().Outtext(" phiNew[z] : ");
					Sysout().Outreal(phiNew[z],4,0);
					Sysout().Outimage();*/
				}
				f[j++] = 1/phiTotal[z] - 1;
			} else if (BulkItVar(z)) {
				int m=1;
				for (int i=1; i<=MolNewQ->GetNumMolecules(); i++) {
					SF_MolStructure* Chain = MolNewQ->GetMolecule(i)->GetMolStructure();
					double dPhiNewAv = 0;
					int k;
					for (k=1; k<=Chain->GetNumDiffSegments(); k++) {
						dPhiNewAv += dPhiNew[m++][z]/Chain->GetAvNumSegments(Chain->GetSegment(k));
					}
					m -= Chain->GetNumDiffSegments();
					for (k=1; k<=Chain->GetNumDiffSegments(); k++) {
						phiNew = (StateNewQ[m])->GetPhi(total);
						double phiSet = (StateQ[m]->GetPhiBulkBoundaries())[Lat->NumBulkBoundary(z)];
						f[j++] = phiNew[z] - phiSet + dPhiNewAv*timeStep - dPhiNew[m][z]*timeStep/Chain->GetAvNumSegments(Chain->GetSegment(k));
				/*	Sysout().Outtext("phibulk[z] : ");
					Sysout().Outreal(phiNew[z],4,0);
					Sysout().Outtext(" phiSet : ");
					Sysout().Outreal(phiSet,4,0);
					Sysout().Outtext(" dPhiNewAv : ");
					Sysout().Outreal(dPhiNewAv,4,0);
					Sysout().Outtext(" phiNew[z] : ");
					Sysout().Outreal(phiNew[z],4,0);
					Sysout().Outimage();*/
						m++;
					}
				}
			}
		}
	}
}
void 
SF_SolveFraaije::residuals(double *const ff, double *const) {
	Vector f(ff, 1, numItVar);
	ComputePhi();
//	GetResiduals(f);
}
void
SF_SolveFraaije::CalculateDPhiNew() {
	int i,j,k,z;
	int count=0;
	SF_MolStructure* Chain;
	SF_MolSegment* Seg;
	SF_MolState* State;
	for (i=1; i<=numLiquidStates; i++) {
		for (z=1; z<=M; z++) {
			dPhiNew[i][z] = 0;
		}
	}
	int N;
	for (i=1; i<=MolNewQ->GetNumMolecules(); i++) {
		Chain = (MolNewQ->GetMolecule(i))->GetMolStructure();
		N = Chain->GetMaxLength();
		for (j=1; j<=Chain->GetNumDiffSegments(); j++) {
			Seg = Chain->GetDiffSegment(j);
			for (k=1; k<=Seg->GetNumStates(); k++) {
				State = Seg->GetState(k);
				count++;
				Vector flux = CalcFluxState(State,N);
				for (z=1; z<=M; z++) {
					dPhiNew[count][z] += flux[z];
				}
			}
		}
	}
}
Vector
SF_SolveFraaije::CalcFluxState(const SF_MolState* State1, const int N1) {
	int i,j,k,z;
	Vector flux(1,M);
	if (State1->GetFreedom() != loose) {
		return flux;
	}
	Vector uPrime = CalcUPrime(State1);
	int N2;
	Vector partFlux;
	SF_MolStructure* Chain;
	SF_MolSegment* Seg;
	SF_MolState* State2;
	for (i=1; i<=MolNewQ->GetNumMolecules(); i++) {
		Chain = (MolNewQ->GetMolecule(i))->GetMolStructure();
		N2 = Chain->GetMaxLength();
		if (N1 > 1 && N2 > 1) {
			continue;
		}
		for (j=1; j<=Chain->GetNumDiffSegments(); j++) {
			Seg = Chain->GetDiffSegment(j);
			for (k=1; k<=Seg->GetNumStates(); k++) {
				State2 = Seg->GetState(k);
				if (State2->GetFreedom() != loose) {
					continue;
				}
				if (State1 == State2) {
					continue;
				}
				partFlux = CalcFluxTwoStates(uPrime,State1->GetPhi(total),State2);
				for (z=1; z<=M; z++) {
					flux[z] += partFlux[z];
				}
			}
		}
	}
	return flux;
}
Vector
SF_SolveFraaije::CalcFluxTwoStates(const Vector uPrime1, 
								   const Vector phi1, 
								   const SF_MolState* State2) {
	Vector phi2 = State2->GetPhi(total);
	Vector uPrime2 = CalcUPrime(State2);
	Vector L(1,M);
	Vector mu(1,M);
	int z;
	for (z=1; z<=M; z++) {
		L[z] = phi1[z]*phi2[z];
		mu[z] = -uPrime1[z] + uPrime2[z];
	}
	Vector flux = Lat->Div1Grad2(L,mu);
	return flux;
}
Vector
SF_SolveFraaije::CalcUPrime(const SF_MolState* State) {
	Vector uPrime(1,M);
	for (int z=1; z<=M; z++) {
		Vector G = State->GetSWF();
		uPrime[z] = -log(G[z]);
		uPrime[z] -= SegNewQ->ChemInt(State,z);
		uPrime[z] += SegNewQ->ChemIntBulk(State);
		if (charged /* && PotItVar(z)*/) {
			uPrime[z] -= State->GetValence()*psi0[z]*ELEM_CHARGE/(BOLTZMANN*TEMPERATURE);
		}
	}
	return uPrime;
}	
void
SF_SolveFraaije::UpdateDensities() {
	for (int i=1; i<=numLiquidStates; i++) {
		for (int z=1; z<=M; z++) {
			dPhi[i][z] = dPhiNew[i][z];
			SWF[i][z] = SWFNew[i][z];
		}
	}
	SegQ->UpdateSWF();
	if (multiState && neutralizerPresent) {
		MolQ->ComputePhi(phiBulkSolvent,phiBulkNeutralizer);
	} else if (multiState) {
		MolQ->ComputePhi(phiBulkSolvent);
	} else {
	 	MolQ->ComputePhi();
	}
}
double
SF_SolveFraaije::CalcStop() const {
	double value = 0;
	for (int i=1; i<=numLiquidStates; i++) {
		Vector phi = (StateQ[i])->GetPhi(total);
		Vector phiStepNew = dPhiNew[i];
		for (int z=1; z<=M; z++) {
			if (fabs(phiStepNew[z]/phi[z]) > value && phi[z] > 0) {
				value = fabs(phiStepNew[z]/phi[z]);
			}
		}
	}
	return value;
}
double
SF_SolveFraaije::CalcError() const {
	double value = 0;
	for (int i=1; i<=numLiquidStates; i++) {
		Vector phi = (StateQ[i])->GetPhi(total);
		Vector phiStep = dPhi[i];
		Vector phiStepNew = dPhiNew[i];
		for (int z=1; z<=M; z++) {
			if (fabs(phiStep[z] - phiStepNew[z])*timeStep/phi[z] > value && phi[z] > 0) {
				value = fabs(phiStep[z] - phiStepNew[z])*timeStep/phi[z];
			}
		}
	}
	return value;	
}
void
SF_SolveFraaije::UpdateSWF() {
	int i,j,z;
	j=1;
	for (z=1; z<=M; z++) {
		if (PotItVar(z)) {
			psi0[z] = x[j++];
		}
		if (SegItVar(z) || BulkItVar(z)) {
			for (i=1; i<=numLiquidStates; i++) {
				SWFNew[i][z] = exp(-x[j++]);
			}
		}
	}
	if (neutralizerPresent) {
		phiBulkNeutralizer = (MolNewQ->GetNeutralizer())->GetPhiBulk();
	}
	phiBulkSolvent = (MolNewQ->GetSolvent())->GetPhiBulk();
	if (charged) {
		Lat->SetBoundaries(psi0);
	}
	SegNewQ->UpdatePhiBulk();
	SegNewQ->UpdateSWF();
}
void
SF_SolveFraaije::UpdateItVar(Boolean offset) {
	int i,j,z;
	j=1;
	for (z=1; z<=M; z++) {
		if (PotItVar(z)) {
			x[j++] = psi0[z];
		}
		if (SegItVar(z) || BulkItVar(z)) {
			for (i=1; i<=numLiquidStates; i++) {
				if (SWFNew[i][z] > 0) {
					x[j] = -log(SWFNew[i][z]);
				}
				j++;
			}
		}
	}
	if (neutralizerPresent) {
		phiBulkNeutralizer = (MolNewQ->GetNeutralizer())->GetPhiBulk();
	}
	phiBulkSolvent = (MolNewQ->GetSolvent())->GetPhiBulk();
}
void
SF_SolveFraaije::Iterate(Boolean firstIterate_, double timeStep_) {
	firstIterate = firstIterate_;
	ProcessInput();
	ComputePhi();
	timeStep = timeStep_;
	iterate(&x[1],numItVar);
	CheckSolution();
	if (firstIterate) {
		CalculateDPhiNew();
	}
}
void
SF_SolveFraaije::SetFrozenSWFZero(void) {
	// to be made
}

Boolean
SF_SolveFraaije::BulkItVar(const int z) const {
	return Lat->BulkBoundary(z);
}
