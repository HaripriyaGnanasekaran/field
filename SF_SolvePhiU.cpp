#include "SF_SolvePhiU.h"
#include <iostream> 

SF_SolvePhiU::SF_SolvePhiU(Boolean compute_,SF_ReactionList* ReactionQ_, SF_MoleculeList* MolQ_, SF_SegmentList* SegQ_, Lattice* Lat_, Input* MyInput_) 
	: SF_Solve(compute_,ReactionQ_,MolQ_,SegQ_,Lat_,MyInput_) {

	use_vector_iterations = false; LM_BFGS=false; pseudohessian=true; 
	compute = compute_;
	ReactionQ = ReactionQ_;
	MolQ = MolQ_;
	SegQ = SegQ_;
	Lat = Lat_;
	MyInput = MyInput_;
	ProcessInput();
	if (extPotSet) {
		SetExternalPotentials();
	}
	reverseDirection.Dim(1,reverseDirectionRange);
	charged = SegQ->Charged();
	if (MyInput->GetNumNames("reaction") > 0) {
		multiState = true;
	} else {
		multiState = false;
	}
	smallAlphaCount = 0;
	M = Lat->GetTotalNumLayers();
	U_offset.Dim(1,M); for (int i=1; i<=M; i++) U_offset[i]=0.0;
	Offset = false;
	reset_hess = false;
	reset_pseudohessian = false;  
	int i,j;
	Boolean equal;
	SF_State* State1;
	SF_State* State2;
	numLiquidStates = 0;
	for (i=1; i<=SegQ->GetNumStates(); i++) {
		State1 = SegQ->GetState(i);
		if (State1->GetFreedom() != frozen) {
			numLiquidStates++;
		}
	}
	StateQ.Dim(1,numLiquidStates);
	numLiquidStates = 0;
	for (i=1; i<=SegQ->GetNumStates(); i++) {
		State1 = SegQ->GetState(i);
		if (State1->GetFreedom() != frozen) {
			StateQ[++numLiquidStates] = SegQ->GetState(i);
		}
	}
	numDiffItVar = numLiquidStates;
	SF_State* StateDummy;
	for (i=1; i<=numLiquidStates; i++) {
		State1 = StateQ[i];
		for (j=i+1; j<=numLiquidStates; j++) {
			State2 = StateQ[j];
			equal = SegQ->ChemIntStatesEqual(State1,State2);
			if (equal && !State1->ExternalPotential()
					&& !State2->ExternalPotential()) {
				if (j!=i+1) {
					StateDummy = StateQ[i+1];
					StateQ[i+1] = State2;
					StateQ[j] = StateDummy;
				}
				i++;
				numDiffItVar--;
			}
		}
	}
	NumStatesForItVar.Dim(1,numDiffItVar);
	numDiffItVar=0;
	for (i=1; i<=numLiquidStates; i++) {
		State1 = StateQ[i];
		NumStatesForItVar[++numDiffItVar] = 1;
		for (j=i+1; j<=numLiquidStates; j++) {
			State2 = StateQ[j];
			equal = SegQ->ChemIntStatesEqual(State1,State2);
			if (equal && !State1->ExternalPotential()
					&& !State2->ExternalPotential()) {
				(NumStatesForItVar[numDiffItVar])++;
				i++;
			}
		}
	}
	SWF.Dim(1,numLiquidStates);
	for (i=1; i<=numLiquidStates; i++) {
		SWF[i].Dim(1,M);
		(StateQ[i])->SetSWF(SWF[i]);
	}
	PHI.Dim(1,numLiquidStates);
	for (i=1; i<=numLiquidStates; i++) {
		PHI[i].Dim(1,M);
	}

	SegQ->UpdateSWF();
	if (MolQ->GetNeutralizer() == NULL) neutralizerPresent = false;
	else neutralizerPresent = true; 
	if (multiState) {
		if (!neutralizerPresent && ReactionQ->NeutralizerNeeded()) {
			Message(fatal,"Please define a molecule with freedom neutralizer");
		}
	} else {
		if (!neutralizerPresent && SegQ->Charged()) {
			SF_Molecule* Mol;
			for (i=1; i<=MolQ->GetNumMolecules(); i++) {
				Mol = MolQ->GetMolecule(i);
				if (Mol->GetBulkCharge() != 0 && Mol->GetPhiBulk() != 0) {
					Message(fatal,"Please define a molecule with freedom neutralizer");
				}
			}
			double totalCharge = 0;
			for (i=1; i<=MolQ->GetNumMolecules(); i++) {
				Mol = MolQ->GetMolecule(i);
				totalCharge += Mol->GetBulkCharge()*Mol->GetTheta()/Mol->GetChainLength();
			}
			SF_State* State;
			for (i=1; i<=SegQ->GetNumStates(); i++) {
				State = SegQ->GetState(i);
				if (State->GetFreedom() == frozen) {
					totalCharge += State->GetTheta()*State->GetValence();
				}
			}
			if (fabs(totalCharge) > gettolerance()) {
				Message(fatal,"There's no neutralizer defined, which is fine but "
					"the charge of the molecules and surfaces in the system does not add up to zero");
			}
		}
	}	
	int z;

	//mmm 
	ItVarArray.Dim(1,M);
	SF_State* State;
	LatticeRange* LatRange;
	int numStates = SegQ->GetNumStates();
	numItVar = 0;
	for (z=1; z<=M; z++) {
		ItVarArray[z] = (Lat->WithinBoundaries(z));
		for (int i=1; ItVarArray[z] && i<=numStates; i++) {
			State = SegQ->GetState(i);
			if (State->GetFreedom() == frozen) {
				LatRange = State->GetLatRange();
				if (LatRange->InRange(z)) {
					if (LatRange->GetRangeValue(z)>0.999)
					ItVarArray[z] = false;
				}
			}
		}
		if (ItVarArray[z]) {
			numItVar++;
		}
	}
//	ItVarRange = LatticeRangeFile(ItVarArray,M);
	//end mmm
	
	double epsilon = (StateQ[1])->GetEpsilon();
	diffEpsilon = false;
	if (charged) {
		for (i=2; i<=numLiquidStates; i++) {
			if ((StateQ[i])->GetEpsilon() != epsilon) {
				diffEpsilon = true;
			}
		}
	}
	numItVar*=(numDiffItVar*2 +1);
	//u_prime.Dim(1,M);

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
	ESquared.Dim(1,M);
	phiTotal_.Dim(1,M);
	phi_tmp.Dim(1,M);
}
SF_SolvePhiU::~SF_SolvePhiU() {
}
void
SF_SolvePhiU::GetOutput(Output* Out)  {
	Out->PutText("newton",name,"iteration solver","PhiU");
	Out->PutInt("newton",name,"iterations",getiterations());
	Out->PutInt("newton",name,"iterationlimit",getiterationlimit());
	Out->PutBoolean("newton",name,"error occurred",errorOccurred);
	if (errorOccurred) {Out->SetErrorOccurred();}
	Out->PutReal("newton",name,"accuracy",getaccuracy());
	Out->PutReal("newton",name,"tolerance",gettolerance());
	Out->PutReal("newton",name,"deltamin",getdeltamin());
	Out->PutReal("newton",name,"deltamax",getdeltamax());
	//Out->PutBoolean("newton",name,"pseudohessian",pseudohessian);
	Out->PutInt("newton",name,"number iter var",numItVar);
	Out->PutVector("newton",name,"iter var",x,numItVar);
	Vector funct(1,numItVar);
	GetResiduals(funct);
	Out->PutVector("newton",name,"residuals",funct,numItVar);
}
void
SF_SolvePhiU::GetResiduals(Vector f) {
	int i,j,k,z,count,count2;
	//int pitvar=0;
	//int jp;
	double Ex=0,phiz;

	//cout << "Get Residuals" << endl; 
	SF_State* State;
	Vector phiTotal = SegQ->GetPhiTotal();

	for (i=1; i<=numLiquidStates; i++) {
		State = StateQ[i];
		Vector phi=State->GetPhi();
		for (z=1; z<=M; z++) PHI[i][z]=phi[z];
	}
	Vector psi;
	if (charged) {
		if (potentialItVar) {
			psi = SegQ->ComputeElectricPotential(psi0);
			//pitvar=1;
		} else {
			psi = SegQ->ComputeElectricPotential();
		}
	}
  	j=1;
	for (z=1; z<=M; z++) {
		if (PotItVar(z)) {
			f[j++] = psi0[z] - psi[z]; //pot itvar - pot comp
		}
		
		if (SegItVar(z)) {
			//jp=j;
			f[j++]=(1.0-phiTotal[z])/phiTotal[z]; //for u-prime; note that this is -(phitot-1) (maximum for uprime itvar)
			count=0; count2=0; 
			for (i=1; i<=numDiffItVar; i++) {
				count+=NumStatesForItVar[i];
				State = StateQ[count];

				Ex= SegQ->ChemInt(State,z);//  /phiTotal_[z];
				Ex -= SegQ->ChemIntBulk(State);

				if (State->ExternalPotential()) {
					Ex += (State->GetExternalPotential())[z];
				}

				f[j] = (Ex-x[j+1])/phiTotal[z]; //PhiItvar: E computed - E itvar; (minimum for phi itvar)
				j++;
				phiz=0;
				for (k=1; k<=NumStatesForItVar[i]; k++) {
					count2++;
					phiz +=PHI[count2][z];
				} 
				f[j] =(x[j-1]-phiz)/phiTotal[z];///phiz; //EItvar: Phi itvar - phi computed  (maximum for E Itvar)
							  //(normalisation is to keep terms of equal size)
				
				j++; 
				
			}
		}
	}
}
void 
SF_SolvePhiU::residuals(double *const ff, double *const) {
	Vector f(ff, 1, numItVar);
	ComputePhi();
	GetResiduals(f);
}

void 
SF_SolvePhiU::inneriteration(double *const f, double *const itvar, double error) {
}

void
SF_SolvePhiU::UpdateSWF() { //order is psi -> u_prime -> [ ln phi -> U ]x
	int i,j,k,z, count,jp,count2;
	SF_State* State;
	j=1;
	for (z=1; z<=M; z++) {
		if (PotItVar(z)) {
			psi0[z] = x[j++];
		}
		phiTotal_[z]=0;
		if (SegItVar(z)) {			
			jp=j++; //location for u_prime(z)
			count = 0; count2=0;
			for (i=1; i<=numDiffItVar; i++) { 
				for (k=1; k<NumStatesForItVar[i]; k++) {
					count2++;
					PHI[count2][z]=0;	
				} count2++;
				PHI[count2][z] = x[j++]; //fill last one with itvar
				
				//cout << "i" << i << "  z=  " << z << " phi i z " << PHI[count2][z] << endl; 
				for (k=1; k<=NumStatesForItVar[i]; k++) {
					count++;
					//SWF[count][z] = exp(-x[j]);
					SWF[count][z] = exp(-(x[jp]+x[j]));
				}
				j++;
			}
		}
	}
	count=0;
	for (i=1; i<=numDiffItVar; i++) {	
		for (k=1; k<=NumStatesForItVar[i]; k++) { count++;
			State = StateQ[count];
			Vector phi=State->GetPhi();
			for (z=1; z<=M; z++) { phi[z]=PHI[count][z]; phiTotal_[z] +=phi[z];}
			SegQ->PutSides(State);
		}
	}
	count=0;
	//for (i=1; i<=numDiffItVar; i++) {
	//	count+=NumStatesForItVar[i];
	//	State = StateQ[count];
	//	SegQ->PutSides(State);
	//}


	if (charged) {
		Lat->SetBoundaries(psi0);
		double preFactor = (EPS0*Lat->GetSiteDistance())
			/(BOLTZMANN*TEMPERATURE);
		Lat->ElectricFieldSquared(ESquared,psi0,preFactor);
		for (z=1; z<=M; z++) {
			if (SegItVar(z)) {
				count = 0;
				SF_State* State;
				for (i=1; i<=numDiffItVar; i++) {
					for (k=1; k<=NumStatesForItVar[i]; k++) {
						count++;
						State = StateQ[count];
						SWF[count][z] *= exp(-State->GetValence()*psi0[z]*
							ELEM_CHARGE/(BOLTZMANN*TEMPERATURE));
						if (diffEpsilon) {
							SWF[count][z] *= exp(0.5*State->GetEpsilon()*ESquared[z]);
						}
					}
				}
			}
		}
	}
	count = 0;
	for (i=1; i<=numDiffItVar; i++) {
		for (k=1; k<=NumStatesForItVar[i]; k++) {
			count++;
			Lat->SetBoundaries(SWF[count]);
		}
	}
	if (neutralizerPresent) {
		phiBulkNeutralizer = (MolQ->GetNeutralizer())->GetPhiBulk();
	}
	phiBulkSolvent = (MolQ->GetSolvent())->GetPhiBulk();
	if (charged) {
		Lat->SetBoundaries(psi0);
	}
	SegQ->UpdatePhiBulk();
	SegQ->UpdateSWF();
}
void
SF_SolvePhiU::UpdateItVar(Boolean offset) {
	int i,j,k,z,count=0,count2=0,jp;
	double Ex,phiz; 
	SF_State* State;
	ComputePhi();
	for (i=1; i<=numLiquidStates; i++) {
		State = StateQ[i];
		SegQ->PutSides(State);
		Vector phi=State->GetPhi();
		for (z=1; z<=M; z++) PHI[i][z]=phi[z];
	}
	j=1;
	for (z=1; z<=M; z++) {
		if (PotItVar(z)) {
			x[j++] = psi0[z]; //itvar contains potential
		}
		if (SegItVar(z)) {
			count=0; count2=0;
			jp=j;
			x[j++] = 0; //itvar for u-prime. later more is added
			for (i=1; i<=numDiffItVar; i++) { 
				phiz=0;
				for (k=1; k<=NumStatesForItVar[i]; k++) { count2++;
					phiz +=PHI[count2][z];
				}
			        //cout << "z " << z<< "  phiz "<< phiz << endl;
				x[j++]=phiz;

				count+=NumStatesForItVar[i];
				State = StateQ[count];
				Ex= SegQ->ChemInt(State,z);
				Ex -= SegQ->ChemIntBulk(State);
				
				if (State->ExternalPotential()) {
					Ex += (State->GetExternalPotential())[z];
				}
				
				x[j++]=Ex; //it variable contains the interaction only
				
				State = StateQ[count];
				if (SWF[count][z] > 0) {
					x[jp] = -log(SWF[count][z]);

					if (charged) {
						x[jp] -= State->GetValence()*psi0[z]*
								ELEM_CHARGE/(BOLTZMANN*TEMPERATURE);
						if (diffEpsilon) {
							x[jp] += 0.5*State->GetEpsilon()*ESquared[z];
						}
					}
				}
				x[jp] -= Ex;
				//if (offset) x[j-1] -= U_offset[z];
			}
		}
	}

	if (neutralizerPresent) {
		phiBulkNeutralizer = (MolQ->GetNeutralizer())->GetPhiBulk();
	}
	phiBulkSolvent = (MolQ->GetSolvent())->GetPhiBulk();
}

