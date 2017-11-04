#include "SF_SolveCG.h"
#include <iostream>  
#include <math.h>

SF_SolveCG::SF_SolveCG(Boolean compute_,SF_ReactionList* ReactionQ_, SF_MoleculeList* MolQ_, SF_SegmentList* SegQ_, Lattice* Lat_, Input* MyInput_) 
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
	numItVar*=numDiffItVar;
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
}

SF_SolveCG::~SF_SolveCG() {
}

void
SF_SolveCG::GetOutput(Output* Out) const {
	Out->PutText("newton",name,"iteration solver","CG-F");
	Out->PutInt("newton",name,"iterations",getiterations());
	Out->PutInt("newton",name,"iterationlimit",getiterationlimit());
	Out->PutReal("newton",name,"accuracy",getaccuracy());
	Out->PutReal("newton",name,"tolerance",gettolerance());
	Out->PutInt("newton",name,"number iter var",numItVar);
	Out->PutVector("newton",name,"iter var",x,numItVar);
	Out->PutReal("newton",name,"deltamax",getdeltamax());
	Vector funct(1,numItVar);
	GetResiduals(funct);
	Out->PutVector("newton",name,"residuals",funct,numItVar);
}

void
SF_SolveCG::GetResiduals(Vector f) const {
	int i,j,z,count;
	//int count2;
	SF_State* State;
//	SF_State* State2;
	Vector phiTotal = SegQ->GetPhiTotal();
	Vector uPrime(1,numDiffItVar);
	double uPrimeAv;
	Vector psi;
	if (charged) {
		if (potentialItVar) {
			psi = SegQ->ComputeElectricPotential(psi0);
		} else {
			psi = SegQ->ComputeElectricPotential();
		}
	}
	count=0;
	//count2=0;
	
	for (i=1; i<=numLiquidStates; i++) {
	//for (i=1; i<=numDiffItVar; i++) {
		//count+=NumStatesForItVar[i];
		//State = StateQ[count];
		State = StateQ[i];
		SegQ->PutSides(State);
	}

        j=1;
	for (z=1; z<=M; z++) {
		if (PotItVar(z)) {
			f[j++] = psi0[z] - psi[z];
		}
		if (SegItVar(z)) {
			count=0; //count2=0;
			uPrimeAv=0;
			for (i=1; i<=numDiffItVar; i++) {
				count+=NumStatesForItVar[i];
				State = StateQ[count];
				uPrime[i] = -log(SWF[count][z]);
				uPrime[i] -= SegQ->ChemInt(State,z)/phiTotal[z];
				uPrime[i] += SegQ->ChemIntBulk(State);
				if (charged) {
					uPrime[i] -= State->GetValence()*psi0[z]*
						ELEM_CHARGE/(BOLTZMANN*TEMPERATURE);
					if (diffEpsilon) {
						uPrime[i] += 0.5*State->GetEpsilon()*ESquared[z];
					}
				}
				if (State->ExternalPotential()) {
					uPrime[i] -= (State->GetExternalPotential())[z];
				}
				uPrimeAv += uPrime[i];
				
			}
			uPrimeAv /= numDiffItVar;
			
			for (i=1; i<=numDiffItVar; i++){
				f[j++] = - 1.0 + 1.0/phiTotal[z]+ uPrime[i] - uPrimeAv;
				//f[j++] = 0.43*(1-phiTotal[z]) + (uPrime[i] - uPrimeAv); //perhaps better in terms of convergence
			}
			
		}
	}
	//for (z=1; z<=numItVar; z++) f[z]=exp(f[z]/10.0)-1.0;
}

void 
SF_SolveCG::inneriteration(double *const f, double *const itvar, double error) {
}

void 
SF_SolveCG::residuals(double *const ff, double *const) {
	Vector f(ff, 1, numItVar);
	ComputePhi();
	GetResiduals(f);
}

void
SF_SolveCG::UpdateSWF() { //order is psi -> u_prime -> [ E ]x
		int i,j,k,z, count;
	j=1;
	for (z=1; z<=M; z++) {
		if (PotItVar(z)) {
			//psi0[z] = 10*log(tanh(x[j++])+1.0);
			psi0[z] = x[j++];
		}
		if (SegItVar(z)) {
			count = 0;
			for (i=1; i<=numDiffItVar; i++) {
				for (k=1; k<=NumStatesForItVar[i]; k++) {
					count++;
					//SWF[count][z] = exp(-x[j]);
					//SWF[count][z] = exp(-10*log(tanh(x[j])+1.0)-U_offset[z]);
					SWF[count][z] = exp(-x[j]-U_offset[z]);
				}
				j++;
			}
		}
	}
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
SF_SolveCG::UpdateItVar(Boolean offset) {
	int i,j,z;
	SF_State* State;
	j=1;
	for (z=1; z<=M; z++) {
		if (PotItVar(z)) {
			x[j++] = psi0[z];
		}
		if (SegItVar(z)) {
			int count=0;
			for (i=1; i<=numDiffItVar; i++) {
				for (int k=1; k<=NumStatesForItVar[i]; k++) {
					count++;
					State = StateQ[count];
					if (SWF[count][z] > 0) {
						x[j] = -log(SWF[count][z]);

						if (offset) x[j] -= U_offset[z];
						if (charged) {
							x[j] -= State->GetValence()*psi0[z]*
								ELEM_CHARGE/(BOLTZMANN*TEMPERATURE);
							if (diffEpsilon) {
								x[j] += 0.5*State->GetEpsilon()*ESquared[z];
							}
						}
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
	
	//for (int z=1; z<=numItVar; z++) x[z]=exp(x[z]/10.0)-1.0;
}

void
SF_SolveCG::ProcessInput() {
	int numNames = MyInput->GetNumNames("newton",0,1);
	if (numNames == 1) {
		Array<Text> names(1,1);
		names = MyInput->GetNames("newton");
		name = names[1];
		Array<Text> param(1,15);
		param[1] = "solver";
		param[2] = "iterationlimit";
		param[3] = "tolerance";
		param[4] = "e_info";
		param[5] = "s_info";
		param[6] = "deltamax";
		param[7] = "print_iteration_info";
		param[8] = "initial_guess";
		param[9] = "initial_guess_input_file";
		param[10] = "initial_guess_output_file";
		param[11] = "write_initial_guess";
		param[12] = "initial_guess_symmetry";
		param[13] = "external_potential_file";
		param[14] = "number_of_stored_vectors";
		param[15] = "m"; 
 
		MyInput->CheckParameterNames("newton",name,param);

		e_info = MyInput->GetBoolean("newton",name,"e_info",true);
		s_info = MyInput->GetBoolean("newton",name,"s_info",false);
		setiterationlimit(MyInput->GetInt("newton",name,"iterationlimit",0,INT_MAX,5000));
		settolerance(MyInput->GetReal("newton",name,"tolerance",1e-15,1e-1,1e-7));
		setdeltamax(MyInput->GetReal("newton",name,"deltamax",1e-15,DBL_MAX,0.5));
		CG_F=true;
		amount_of_stored_vectors = MyInput->GetInt("newton",name,"m",0,64,0);
		if (amount_of_stored_vectors==0) amount_of_stored_vectors = MyInput->GetInt("newton",name,"number_of_stored_vectors",1,64,1);

		print_iteration_info = MyInput->GetInt("newton",name,"print_iteration_info",1,1000,1);
		Array<Text> initialGuess(1,6);
		initialGuess[1] = "file";
		initialGuess[2] = "previous_result";
		initialGuess[3] = "compute";
		initialGuess[4] = "bulk";
		initialGuess[5] = "none";
		initialGuess[6] = "move";
		int initGuess;
		if (MyInput->ValueSet("newton",name,"initial_guess_input_file")) {
			if (!MyInput->ValueSet("newton",name,"initial_guess")) {
				Message(fatal,"Set 'newton' : " + name + " : initial_guess' to 'file'"
					" when setting 'initial_guess_input_file'.");
			}
		}
		// default: use results from previous calculation
		initGuess = MyInput->GetChoice("newton",name,"initial_guess",initialGuess,2);
		switch(initGuess) {
		case(1):
			MyInput->AlwaysCombineParam("newton",name,"initial_guess","initial_guess_input_file");
			initialGuessInputFile = MyInput->GetText("newton",name,"initial_guess_input_file");
			initialGuessInputFile = parseFileName(initialGuessInputFile, MyInput);
			break;
		case(2):
			MyInput->DontCombineParam("newton",name,"initial_guess","initial_guess_input_file");
			break;
		case(3):
			MyInput->DontCombineParam("newton",name,"initial_guess","initial_guess_input_file");
			break;
		case(4):
			MyInput->DontCombineParam("newton",name,"initial_guess","initial_guess_input_file");
			break;
		case(5):
			MyInput->DontCombineParam("newton",name,"initial_guess","initial_guess_input_file");
			break;
		case(6):
			MyInput->DontCombineParam("newton",name,"initial_guess","initial_guess_input_file");
			break;
		}
		if (MyInput->ValueSet("newton",name,"initial_guess_output_file")) {
			initialGuessOutputFile = MyInput->GetText("newton",name,"initial_guess_output_file");
			initialGuessOutputFile = parseFileName(initialGuessOutputFile, MyInput);
			writeInitialGuess = MyInput->GetBoolean("newton",name,"write_initial_guess",true);;
		} else {
			writeInitialGuess = false;
		}
		Array<Text> initialGuessSym(1,2);
		initialGuessSym[1] = "first_to_last";
		initialGuessSym[2] = "bounds_to_middle";
		int choice = MyInput->GetChoice("newton",name,"initial_guess_symmetry",initialGuessSym,1);
		if (choice == 1) initGuessSym = firstToLast;
		if (choice == 2) initGuessSym = boundsToMiddle;
				
		if (MyInput->ValueSet("newton",name,"external_potential_file")) {
			externalPotentialFile = MyInput->GetText("newton",name,"external_potential_file");
			extPotSet = true;
		} else {
			extPotSet = false;
		}

	}
}
 



