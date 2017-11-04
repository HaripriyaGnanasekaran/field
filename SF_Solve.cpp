#include <iostream>
#include "SF_Solve.h"
SF_Solve::SF_Solve(Boolean compute_,
				   SF_ReactionList* ReactionQ_,
				   SF_MoleculeList* MolQ_,
				   SF_SegmentList* SegQ_,
				   Lattice* Lat_,
				   Input* MyInput_) {
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
		Text latName;
		MyInput ->GetNumNames("lat",1,1);
		Array<Text> latNames = MyInput->GetNames("lat");
		latName = latNames[1];
		if (!MyInput->ValueSet("lat",latName,"bondlength"))
		Message(fatal,"Property bondlength not found in 'lat': System contains charges and therefore you need to set the 'bondlength'parameter in 'lat'. A good choice would be to use the Bjerrumlength for this (approx 6e-10). In surfactant problems you meight consider a smaller value, e.g. 3e-10 to allow for a finer difinition of the molecules.");
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
SF_Solve::~SF_Solve() {

}
Text
SF_Solve::GetName() const {
	return name;
}
void
SF_Solve::GetOutput(Output* Out) const {
	Out->PutText("newton",name,"iteration solver","normal");
	if (MyInput->ValueSet("newton",name,"method")) {
		Out->PutText("newton",name,"iteration method",MyInput->GetText("newton",name,"method"));}
	else
		Out->PutText("newton",name,"iteration method","pseudohessian");
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
	Out->PutInt("newton",name,"number iter var",numItVar);
	Out->PutVector("newton",name,"iter var",x,numItVar);
	if (LM_BFGS) Out->PutInt("newton",name,"number of stored arrays",getNumStoredArrays());
	Vector funct(1,numItVar);
	GetResiduals(funct);
	Out->PutVector("newton",name,"residuals",funct,numItVar);
}

void
SF_Solve::Iterate() {
	ComputeAlphaBulk();
	numIterationsSinceHessian = 0;
	minAccuracySoFar = 1e30;
	if (reset_hess) resethessian();
	if (use_vector_iterations) {
		vector_iterate(&x[1],numItVar);
	}
	else {
		cout << "iterate " << numItVar << endl;
		iterate(&x[1],numItVar);
	}
	if (multiState) {
		ReactionQ->InnerIterate();
		while (ReactionQ->getiterations() > 0) {
			samehessian = true;
			numIterationsSinceHessian = 0;
			minAccuracySoFar = 1e30;
			if (e_info || s_info) cout << "NEW ITERATION with adjusted alpha bulk values" << endl;
			if (use_vector_iterations) {
				vector_iterate(&x[1],numItVar);
			}
			else {
				iterate(&x[1],numItVar);
			}
			ReactionQ->InnerIterate();
		}
	}
	CheckSolution();
}
Boolean
SF_Solve::ErrorOccurred(void) {
	return errorOccurred;
}
void
SF_Solve::SetErrorOccurred(Boolean a) {
	errorOccurred = a;
}


void
SF_Solve::MoveInitialGuess(int *newpos, int *oldpos, int numpos) {
	int i,rr;
	int xo,yo,zo,xn,yn,zn;
	double xyz;
	//double x1yz,xy1z,xyz1,x_1yz,xy_1z,xyz_1;

	//int n_layers_x = Lat->GetNumLayers(1)-2;
	int n_layers_y = Lat->GetNumLayers(2)-2;
	int n_layers_z = Lat->GetNumLayers(3)-2;
	int jx = (n_layers_y+2)*(n_layers_z+2);
	int jy = n_layers_z+2;
	int jz = 1;

	for (int p=1; p<=numpos; p++) {
		rr = oldpos[p]-1;
		zo = ((rr % jx) % jy)/jz;
		yo = ((rr % jx)-zo*jz)/jy;
		xo = (rr-yo*jy-zo*jz)/jx;

		rr = newpos[p]-1;
		zn = ((rr % jx) % jy)/jz;
		yn = ((rr % jx)-zn*jz)/jy;
		xn = (rr-yn*jy-zn*jz)/jx;
		//cout << "("<< xo <<","<<yo<<","<<zo <<")("<< xn <<","<<yn<<","<<zn <<")"<< endl;

		int count = 0;
		if (charged) {
			xyz=psi0[jx*xo+jy*yo+jz*zo+1];
			//x1yz=psi0[jx*(xo+1)+jy*yo+jz*zo+1];
			//xy1z=psi0[jx*xo+jy*(yo+1)+jz*zo+1];
			//xyz1=psi0[jx*xo+jy*yo+jz*(zo+1)+1];
			//x_1yz=psi0[jx*(xo-1)+jy*yo+jz*zo+1];
			//xy_1z=psi0[jx*xo+jy*(yo-1)+jz*zo+1];
			//xyz_1=psi0[jx*xo+jy*yo+jz*(zo-1)+1];

			psi0[jx*xo+jy*yo+jz*zo+1]=psi0[jx*xn+jy*yn+jz*zn+1];
			//psi0[jx*(xo+1)+jy*yo+jz*zo+1]=psi0[jx*(xn+1)+jy*yn+jz*zn+1];
			//psi0[jx*xo+jy*(yo+1)+jz*zo+1]=psi0[jx*xn+jy*(yn+1)+jz*zn+1];
			//psi0[jx*xo+jy*yo+jz*(zo+1)+1]=psi0[jx*xn+jy*yn+jz*(zn+1)+1];
			//psi0[jx*(xo-1)+jy*yo+jz*zo+1]=psi0[jx*(xn-1)+jy*yn+jz*zn+1];
			//psi0[jx*xo+jy*(yo-1)+jz*zo+1]=psi0[jx*xn+jy*(yn-1)+jz*zn+1];
			//psi0[jx*xo+jy*yo+jz*(zo-1)+1]=psi0[jx*xn+jy*yn+jz*(zn-1)+1];

			if (xyz !=0)psi0[jx*xn+jy*yn+jz*zn+1]=xyz;
			//if (x1yz !=0) psi0[jx*(xn+1)+jy*yn+jz*zn+1]=x1yz;
			//if (xy1z !=0) psi0[jx*xn+jy*(yn+1)+jz*zn+1]=xy1z;
			//if (xyz1 !=0) psi0[jx*xn+jy*yn+jz*(zn+1)+1]=xyz1;
			//if (x_1yz !=0) psi0[jx*(xn-1)+jy*yn+jz*zn+1]=x_1yz;
			//if (xy_1z !=0) psi0[jx*xn+jy*(yn-1)+jz*zn+1]=xy_1z;
			//if (xyz_1 !=0) psi0[jx*xn+jy*yn+jz*(zn-1)+1]=xyz_1;

		}
		if (diffEpsilon) {
			xyz=ESquared[jx*xo+jy*yo+jz*zo+1];
			//x1yz=ESquared[jx*(xo+1)+jy*yo+jz*zo+1];
			//xy1z=ESquared[jx*xo+jy*(yo+1)+jz*zo+1];
			//xyz1=ESquared[jx*xo+jy*yo+jz*(zo+1)+1];
			//x_1yz=ESquared[jx*(xo-1)+jy*yo+jz*zo+1];
			//xy_1z=ESquared[jx*xo+jy*(yo-1)+jz*zo+1];
			//xyz_1=ESquared[jx*xo+jy*yo+jz*(zo-1)+1];

			ESquared[jx*xo+jy*yo+jz*zo+1]=ESquared[jx*xn+jy*yn+jz*zn+1];
			//ESquared[jx*(xo+1)+jy*yo+jz*zo+1]=ESquared[jx*(xn+1)+jy*yn+jz*zn+1];
			//ESquared[jx*xo+jy*(yo+1)+jz*zo+1]=ESquared[jx*xn+jy*(yn+1)+jz*zn+1];
			//ESquared[jx*xo+jy*yo+jz*(zo+1)+1]=ESquared[jx*xn+jy*yn+jz*(zn+1)+1];
			//ESquared[jx*(xo-1)+jy*yo+jz*zo+1]=ESquared[jx*(xn-1)+jy*yn+jz*zn+1];
			//ESquared[jx*xo+jy*(yo-1)+jz*zo+1]=ESquared[jx*xn+jy*(yn-1)+jz*zn+1];
			//ESquared[jx*xo+jy*yo+jz*(zo-1)+1]=ESquared[jx*xn+jy*yn+jz*(zn-1)+1];

			if (xyz !=0)ESquared[jx*xn+jy*yn+jz*zn+1]=xyz;
			//if (x1yz !=0) ESquared[jx*(xn+1)+jy*yn+jz*zn+1]=x1yz;
			//if (xy1z !=0) ESquared[jx*xn+jy*(yn+1)+jz*zn+1]=xy1z;
			//if (xy1z1 !=0) ESquared[jx*xn+jy*yn+jz*(zn+1)+1]=xyz1;
			//if (x_1yz !=0) ESquared[jx*(xn-1)+jy*yn+jz*zn+1]=x_1yz;
			//if (xy_1z !=0) ESquared[jx*xn+jy*(yn-1)+jz*zn+1]=xy_1z;
			//if (xyz_1 !=0) ESquared[jx*xn+jy*yn+jz*(zn-1)+1]=xyz_1;
		}

		for (i=1; i<=numDiffItVar; i++){
			for (int k=1; k<=NumStatesForItVar[i]; k++) {
				count++;
				//cout << SWF[count][jx*xo+jy*yo+jz*zo+1] <<","<<SWF[count][jx*xn+jy*yn+jz*zn+1]<< endl;
				xyz=SWF[count][jx*xo+jy*yo+jz*zo+1];
				//x1yz=SWF[count][jx*(xo+1)+jy*yo+jz*zo+1];
				//xy1z=SWF[count][jx*xo+jy*(yo+1)+jz*zo+1];
				//xyz1=SWF[count][jx*xo+jy*yo+jz*(zo+1)+1];
				//x_1yz=SWF[count][jx*(xo-1)+jy*yo+jz*zo+1];
				//xy_1z=SWF[count][jx*xo+jy*(yo-1)+jz*zo+1];
				//xyz_1=SWF[count][jx*xo+jy*yo+jz*(zo-1)+1];

				if (SWF[count][jx*xn+jy*yn+jz*zn+1]!=0) SWF[count][jx*xo+jy*yo+jz*zo+1]=SWF[count][jx*xn+jy*yn+jz*zn+1];
				//if (SWF[count][jx*(xn+1)+jy*yn+jz*zn+1] !=0) SWF[count][jx*(xo+1)+jy*yo+jz*zo+1]=SWF[count][jx*(xn+1)+jy*yn+jz*zn+1];
				//if (SWF[count][jx*xn+jy*(yn+1)+jz*zn+1]!=0) SWF[count][jx*xo+jy*(yo+1)+jz*zo+1]=SWF[count][jx*xn+jy*(yn+1)+jz*zn+1];
				//if (SWF[count][jx*xn+jy*yn+jz*(zn+1)+1] !=0) SWF[count][jx*xo+jy*yo+jz*(zo+1)+1]=SWF[count][jx*xn+jy*yn+jz*(zn+1)+1];
				//if (SWF[count][jx*(xn-1)+jy*yn+jz*zn+1]!=0) SWF[count][jx*(xo-1)+jy*yo+jz*zo+1]=SWF[count][jx*(xn-1)+jy*yn+jz*zn+1];
				//if (SWF[count][jx*xn+jy*(yn-1)+jz*zn+1]!=0) SWF[count][jx*xo+jy*(yo-1)+jz*zo+1]=SWF[count][jx*xn+jy*(yn-1)+jz*zn+1];
				//if (SWF[count][jx*xn+jy*yn+jz*(zn-1)+1] !=0) SWF[count][jx*xo+jy*yo+jz*(zo-1)+1]=SWF[count][jx*xn+jy*yn+jz*(zn-1)+1];

				if (xyz !=0) SWF[count][jx*xn+jy*yn+jz*zn+1]=xyz;
				//if (x1yz !=0) SWF[count][jx*(xn+1)+jy*yn+jz*zn+1]=x1yz;
				//if (xy1z !=0) SWF[count][jx*xn+jy*(yn+1)+jz*zn+1]=xy1z;
				//if (xyz1 !=0) SWF[count][jx*xn+jy*yn+jz*(zn+1)+1]=xyz1;
				//if (x_1yz !=0) SWF[count][jx*(xn-1)+jy*yn+jz*zn+1]=x_1yz;
				//if (xy_1z !=0) SWF[count][jx*xn+jy*(yn-1)+jz*zn+1]=xy_1z;
				//if (xyz_1 !=0) SWF[count][jx*xn+jy*yn+jz*(zn-1)+1]=xyz_1;
			}
		}
	}
}

void
SF_Solve::SetInitialGuess(SF_Solve* SolveOld) {
	//SetUOffset();
	Array<Text> initialGuess(1,6);
	initialGuess[1] = "file";
	initialGuess[2] = "previous_result";
	initialGuess[3] = "compute";
	initialGuess[4] = "bulk";
	initialGuess[5] = "none";
	initialGuess[6] = "move";
	int initGuess = MyInput->GetChoice("newton",name,"initial_guess",initialGuess,2);
	switch(initGuess) {
		case(1):
			SetInitialGuessFromFile();SetFrozenSWFZero();UpdateItVar(Offset);
			break;
		case(2):
			if (SolveOld->MolQ != MolQ) { //when this is not the first calculation
				SetInitialGuessFromPrevious(SolveOld->Lat,SolveOld->SegQ, SolveOld->MolQ);
				SetFrozenSWFZero();UpdateItVar(Offset);
 			} else {
				ComputeAlphaBulk();SetFrozenSWFZero();UpdateItVar(false);
			}
			break;
		case(3):
			ComputeInitialGuess();
			ComputeAlphaBulk();SetFrozenSWFZero();UpdateItVar(false);
			break;
		case(4):
			ComputeAlphaBulk();SetFrozenSWFZero();UpdateItVar(false);
			break;
		case(5):
			NoInitialGuess(); SetFrozenSWFZero();UpdateItVar(false);
			break;
		case(6):SetFrozenSWFZero();UpdateItVar(false);
			break;
		default:SetFrozenSWFZero();UpdateItVar(false);
			break;
	}
}
Boolean
SF_Solve::NewSolveNeeded(SF_Solve* SolveOld) {
	if (MolQ == SolveOld->MolQ) {  // Old solve is not a calculation
		return true;
	}
	if (SegQ->GetNumStates() != SolveOld->SegQ->GetNumStates()) {
		return true;
	}
	if (numItVar != SolveOld->numItVar) {
		return true;
	}
	if (!transferHessian) {
		return true;
	}
	return false;
}
void
SF_Solve::UpdateSolve(SF_Solve* SolveNew) {
	MyInput = SolveNew->MyInput;
	Lat = SolveNew->Lat;
	SF_SegmentList* SegQOld = SegQ;
	SF_MoleculeList* MolQOld = MolQ;
	Text stateName;
	int numStates = SegQ->GetNumStates();
	SegQ = SolveNew->SegQ;
	MolQ = SolveNew->MolQ;
	SF_State* StateNew;
	SF_State* StateOld;
	int i;
	for (i=1; i<=numStates; i++) {
		StateOld = SegQOld->GetState(i);
		stateName = StateOld->GetName();
		if (SegQ->StateDefined(stateName)) {
			StateNew = SegQ->GetState(stateName);
			StateOld->SetAlphaBulk(StateNew->GetAlphaBulk());
		}
	}
	if (neutralizerPresent) {
		SF_Molecule* NeutralizerOld = MolQOld->GetNeutralizer();
		if (NeutralizerOld != NULL) {
			(MolQ->GetNeutralizer())->SetPhiBulk(NeutralizerOld->GetPhiBulk());
			phiBulkNeutralizer = (MolQ->GetNeutralizer())->GetPhiBulk();
		} else {
			phiBulkNeutralizer = 0;
		}
	}
	ReactionQ = SolveNew->ReactionQ;
	if (multiState) {
		ReactionQ->ComputeAlphaBulk();
	}
	for (i=1; i<=numLiquidStates; i++) {
		StateQ[i] = SolveNew->StateQ[i];
		(StateQ[i])->SetSWF(SWF[i]);
	}
	ProcessInput();
	SegQ->UpdateSWF();
	SegQ->UpdatePhiBulk();
	samehessian = true;
}
void
SF_Solve::WriteInitialGuessToFile() const{
	if (!writeInitialGuess) return;
	ofstream Out;
	Out.open(initialGuessOutputFile.MainC());
	Out.precision(20);
	Out.setf(ios::showpoint);
	Out.setf(ios::scientific);
	SF_State* State;
	Vector vectorOut;
	int gradients = Lat->GetNumGradients();
	Out << "gradients" << endl;
	Out << gradients << endl;
	int i;
	for (i=1; i<=gradients; i++) {
		Out << Lat->GetNumLayers(i) << endl;
	}
	for (i=1; i<=numLiquidStates; i++) {
		State = StateQ[i];
		Out << "molecule" << endl;
		Out << "all" << endl;
		Out << "state" << endl;
		Out << State->GetName().MainC() << endl;
		vectorOut = State->GetSWF();
		for (int z=1; z<=M; z++) {
			Out << vectorOut[z] << endl;
		}
	}
	if (charged) {
		Out << "electric potential" << endl;
		for (int z=1; z<=M; z++) {
			Out << ESquared[z] << endl;
		}
		if (diffEpsilon) {
			Out << "E^2" << endl;
			for (int z=1; z<=M; z++) {
				Out << ESquared[z] << endl;
			}
		}
	}
	Out << "phibulk solvent" << endl;
	Out << (MolQ->GetSolvent())->GetPhiBulk() << endl;
	if (MolQ->GetNeutralizer() != NULL) {
		Out << "phibulk neutralizer" << endl;
		Out << (MolQ->GetNeutralizer())->GetPhiBulk() << endl;
	}
	for (i=1; i<=numLiquidStates; i++) {
		State = StateQ[i];
		Out << "alphabulk" << endl;
		Out << State->GetName().MainC() << endl;
		Out << State->GetAlphaBulk() << endl;
	}
	Out.close();
}
int
SF_Solve::GetNumItVar() const {
	return numItVar;
}
Vector
SF_Solve::GetItVar() const {
	return x;
}
void
SF_Solve::GetResiduals(Vector f) const {
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
				f[j++] = - 1.0 + 1.0/phiTotal[z] + uPrime[i] - uPrimeAv;
				//f[j++] = 0.43*(1-phiTotal[z]) + (uPrime[i] - uPrimeAv); //perhaps better in terms of convergence
			}
		}
		//else {
			//if (getiterations()==10 ) cout << "z = " << z << "phit =" << phiTotal[z] << endl;
		//}
	}
}

double
SF_Solve::GetTolerance() const {
	return gettolerance();
}
void
SF_Solve::Jump(int jump) {
	for (int i=1; i<=SegQ->GetNumStates(); i++) {
		Vector swf = SegQ->GetState(i)->GetSWF();
		if (jump < 0) {
			for (int z = M; z>1; z--) {
				swf[z] = swf[z-1];
			}
			swf[1] = swf[2];
		} else {
			for (int z=2; z<=M; z++) {
				swf[z-1] = swf[z];
			}
			swf[M] = swf[M-1];
		}
	}
	if (charged) {
		if (jump < 0) {
			for (int z = M; z>1; z--) {
				psi0[z] = psi0[z-1];
				if (diffEpsilon) {
					ESquared[z] = ESquared[z-1];
				}
			}
			psi0[1] = psi0[2];
			if (diffEpsilon) {
				ESquared[1] = ESquared[2];
			}
		} else {
			for (int z=2; z<=M; z++) {
				psi0[z-1] = psi0[z];
				if (diffEpsilon) {
					ESquared[z-1] = ESquared[z];
				}
			}
			psi0[M] = psi0[M-1];
			if (diffEpsilon) {
				ESquared[M] = ESquared[M-1];
			}
		}
	}
	UpdateItVar(false);
	SegQ->UpdateSWF();
	ComputePhi();
}
//Boolean
//SF_Solve::GetOverflowProtectionNeeded(void) const {
//	return overflowProtectionNeeded;
//}
void
SF_Solve::residuals(double *const ff, double *const) {
	Vector f(ff, 1, numItVar);
	ComputePhi();
	GetResiduals(f);
}
void
SF_Solve::residuals(Vector f) {
	ComputePhi();
	GetResiduals(f);
}

void
SF_Solve::inneriteration(double *const f,
						 double *const itva,
						 double) {
	Vector funct(f, 1, numItVar);
	Vector itvar(itva, 1, numItVar);

	if (getiterations() > 0) {
		samehessian = false;
	}

	if (reset_pseudohessian) {
		reset_pseudohessian=false; pseudohessian = true;
	}

	if (getaccuracy() < minAccuracySoFar && getiterations() > 0 && getaccuracy() == fabs(getaccuracy()) ) {
		minAccuracySoFar = getaccuracy();
	}

	if (getaccuracy() > minAccuracySoFar*resetHessianCriterion && getaccuracy() == fabs(getaccuracy()) ) {
		if (s_info) {
			cout << getaccuracy() << '\t' << minAccuracySoFar << '\t' << resetHessianCriterion << endl;
			cout << "walking backwards: newton reset" << endl;
		}
		resethessian();
		minAccuracySoFar *=1.5;

		if (getdeltamax() >0.005) {setdeltamax(0.9*getdeltamax());}
		reverseDirection.Dim(1,reverseDirectionRange);
		numIterationsSinceHessian = 0;

	}

	if (historyTest) {
		int i;
		Vector f1(1,numItVar);
		GetResiduals(f1);
		Vector x1(1,numItVar);
		for (i=1; i<=numItVar; i++) {
			x1[i]=x[i];
			x[i]*=1.00001;
		}
		Vector f2(1,numItVar);
		GetResiduals(f2);
		for (i=1; i<=numItVar; i++) {
			x[i]=x1[i];
		}
		Vector f3(1,numItVar);
		GetResiduals(f3);
		for (i=1; i<=numItVar; i++) {
			if (f1[i] != f3[i]) Message(fatal,"History in residuals");
		}
	}

	if (getalpha() < smallAlpha) {
		smallAlphaCount++;
	} else {
		smallAlphaCount = 0;
	}

	if (smallAlphaCount == maxNumSmallAlpha) {
		smallAlphaCount = 0;
		reverseDirection.Dim(1,reverseDirectionRange);
		if (!s_info) {
			cout << "too many small alphas: newton reset" << endl;
		}
		resethessian();
		if (getdeltamax() >0.005) setdeltamax(0.9*getdeltamax());
		numIterationsSinceHessian = 0;
	}

	if (!getnewtondirection() && pseudohessian) {
		reverseDirection[getiterations()%reverseDirectionRange + 1] = 1;
	} else {
		reverseDirection[getiterations()%reverseDirectionRange + 1] = 0;
	}

	int numReverseDirection = 0;
	for (int i=1; i<=reverseDirectionRange; i++) {
		if (reverseDirection[i] > 0.5) {
			numReverseDirection++;
		}
	}

	numIterationsSinceHessian++;
	double frReverseDirection = double(numReverseDirection)/reverseDirectionRange;
	if ((frReverseDirection > maxFrReverseDirection && pseudohessian && getaccuracy() < minAccuracyForHessian) && !LM_BFGS) {
		Message(literal,"Bad convergence (reverse direction), computing full hessian...");
		pseudohessian = false; reset_pseudohessian =true;
		reverseDirection.Dim(1,reverseDirectionRange);
		numIterationsSinceHessian = 0;
	} else if ((numIterationsSinceHessian >= numIterationsForHessian &&
				getiterations() > 0 && getaccuracy() < minAccuracyForHessian && getminimum() < minAccuracyForHessian)&& !LM_BFGS) {
		Message(literal,"Still no solution, computing full hessian...");
		pseudohessian = false; reset_pseudohessian =true;
		numIterationsSinceHessian = 0;

	}
//	for (int i=1; i<=MolQ->GetNumMolecules(); i++) {
//		SF_Molecule* Mol = MolQ->GetMolecule(i);
//		cout << Mol->GetLnGN() << endl;
//		if (Mol->GetLnGN() > 100 || Mol->GetLnGN() < -100) {
//
//			if (e_info) {
//				cout << "Overflow protection enabled" << endl;
//			}
//		}
//
//	}
}
void
SF_Solve::ComputePhi() {

	UpdateSWF();
	if (multiState && neutralizerPresent) {
		MolQ->ComputePhi(phiBulkSolvent,phiBulkNeutralizer);
	} else if (multiState) {
		MolQ->ComputePhi(phiBulkSolvent);
	} else {
	 	MolQ->ComputePhi();
	}
}

void
SF_Solve::ReComputePhi(bool doit) {
	UpdateSWF();
	if (doit){
		MolQ->ReComputePhi();
	}
	else {
		MolQ->ComputePhi();
	}

}

void
SF_Solve::UpdateSWF() {
	int i,j,k,z, count;
	j=1;
	for (z=1; z<=M; z++) {
		if (PotItVar(z)) {
			psi0[z] = x[j++];
		}
		if (SegItVar(z)) {
			count = 0;
			for (i=1; i<=numDiffItVar; i++) {
				for (k=1; k<=NumStatesForItVar[i]; k++) {
					count++;
					//SWF[count][z] = exp(-x[j]);
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
SF_Solve::UpdateItVar(Boolean Compensate_offset) {
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

						if (Compensate_offset) x[j] -= U_offset[z];
						if (charged)
						{
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
}
void
SF_Solve::SetInitialGuessFromFile() {
	Infile in(initialGuessInputFile);
	if (!in.Open(Blanks(100))) {
		Message(fatal,"initial_guess_input_file : "
			+ initialGuessInputFile + " cannot be opened");
	}
	in.Inimage();
	if (in.Endfile()) {
		Message(fatal, "initial_guess_input_file : "
			+ initialGuessInputFile + " is empty");
	}
	Text line;
	line = Copy(in.Image.Strip());
	if (*line != *Copy("gradients")) {
		Message(fatal,"error reading '" +initialGuessInputFile
			+ "' unrecognized line: " + line);
	}
	in.Inimage();
	line = Copy(in.Image.Strip());
	int gradients = line.Getint();
	int Min=0;
	int MXin=0;
	int MYin=0;
	int MZin=0;
	if (gradients == 1) {
		in.Inimage();
		line = Copy(in.Image.Strip());
		Min = line.Getint();
	}
	if (gradients == 2) {
		in.Inimage();
		line = Copy(in.Image.Strip());
		MXin = line.Getint();
		in.Inimage();
		line = Copy(in.Image.Strip());
		MYin = line.Getint();
		Min = MXin*MYin;
	}
	if (gradients == 3) {
		in.Inimage();
		line = Copy(in.Image.Strip());
		MXin = line.Getint();
		in.Inimage();
		line = Copy(in.Image.Strip());
		MYin = line.Getint();
		in.Inimage();
		line = Copy(in.Image.Strip());
		MZin = line.Getint();
		Min = MXin*MYin*MZin;
	}
	Text stateName;
	Vector guess;
	guess.Dim(1,Min);
	if (!in.Endfile()) in.Inimage();
	line = Copy(in.Image.Strip());
	while (*line == *Copy("molecule")) {
		in.Inimage();
		line = Copy(in.Image.Strip());
		in.Inimage();
		line = Copy(in.Image.Strip());
		in.Inimage();
		stateName = Copy(in.Image.Strip());
		for (int z=1; z<=Min; z++) {
			in.Inimage();
			line = Copy(in.Image.Strip());
			guess[z] = line.Getreal();
		}
		if (SegQ->StateDefined(stateName)) {
			SF_State* State = SegQ->GetState(stateName);
			Vector vectorOut = State->GetSWF();
			if (gradients == 1) {
				CopyInitialGuess(guess,vectorOut,Min);
			}
			if (gradients == 2) {
				CopyInitialGuess(guess,vectorOut,MXin,MYin);
			}
			if (gradients == 3) {
				CopyInitialGuess(guess,vectorOut,MXin,MYin,MZin);
			}
		}
		if (!in.Endfile()) in.Inimage();
		line = Copy(in.Image.Strip());
	}
	if (*line == *Copy("electric potential") && charged) {
		for (int z=1; z<=Min; z++) {
			in.Inimage();
			line = Copy(in.Image.Strip());
			guess[z] = line.Getreal();
		}
		if (gradients == 1) {
			CopyInitialGuess(guess,psi0,Min);
		}
		if (gradients == 2) {
			CopyInitialGuess(guess,psi0,MXin,MYin);
		}
		if (gradients == 3) {
			CopyInitialGuess(guess,psi0,MXin,MYin,MZin);
		}
		if (!in.Endfile()) in.Inimage();
		line = Copy(in.Image.Strip());
	}
	if (*line == *Copy("E^2")) {
		for (int z=1; z<=Min; z++) {
			in.Inimage();
			line = Copy(in.Image.Strip());
			guess[z] = line.Getreal();
		}
		if (gradients == 1 && diffEpsilon) {
			CopyInitialGuess(guess,ESquared,Min);
		}
		if (gradients == 2 && diffEpsilon) {
			CopyInitialGuess(guess,ESquared,MXin,MYin);
		}
		if (gradients == 3 && diffEpsilon) {
			CopyInitialGuess(guess,ESquared,MXin,MYin,MZin);
		}
		if (!in.Endfile()) in.Inimage();
		line = Copy(in.Image.Strip());
	}
	if (*line == *Copy("phibulk solvent")) {
		if (!in.Endfile()) in.Inimage();
		line = Copy(in.Image.Strip());
		phiBulkSolvent = line.Getreal();
		if (!in.Endfile()) in.Inimage();
		line = Copy(in.Image.Strip());
	} else {
		phiBulkSolvent = 0;
	}
	(MolQ->GetSolvent())->SetPhiBulk(phiBulkSolvent);
	if (*line == *Copy("phibulk neutralizer")) {
		if (!in.Endfile()) in.Inimage();
		line = Copy(in.Image.Strip());
		if (MolQ->GetNeutralizer() != NULL) {
			phiBulkNeutralizer = line.Getreal();
		}
		if (!in.Endfile()) in.Inimage();
		line = Copy(in.Image.Strip());
	} else {
		phiBulkNeutralizer = 0;
	}
	if (neutralizerPresent) {
		(MolQ->GetNeutralizer())->SetPhiBulk(phiBulkNeutralizer);
	}
	while (*line == *Copy("alphabulk")) {
		if (!in.Endfile()) in.Inimage();
		line = Copy(in.Image.Strip());
		if (SegQ->StateDefined(line)) {
			SF_State* State = SegQ->GetState(line);
			if (!in.Endfile()) in.Inimage();
			line = Copy(in.Image.Strip());
			State->SetAlphaBulk(line.Getreal());
			if (!in.Endfile()) in.Inimage();
			line = Copy(in.Image.Strip());
		} else {
			if (!in.Endfile()) in.Inimage();
			line = Copy(in.Image.Strip());
			if (!in.Endfile()) in.Inimage();
			line = Copy(in.Image.Strip());
		}
	}
	if (multiState) {
		ReactionQ->ComputeAlphaBulk();
	}
	if (compute) {
		SegQ->UpdateSWF();
		SegQ->UpdatePhiBulk();
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
	in.Close();
}
void
SF_Solve::SetInitialGuessFromPrevious(const Lattice* LatOld,
									  const SF_SegmentList* SegQOld,
									  const SF_MoleculeList* MolQOld) {
	SF_State* StateOld;
	SF_State* StateNew;
	Text stateName;
	int gradients = LatOld->GetNumGradients();
	int numStates = SegQ->GetNumStates();
	int i;
	for (i=1; i<=numStates; i++) {
		StateNew = SegQ->GetState(i);
		stateName = StateNew->GetName();
		if (SegQOld->StateDefined(stateName)) {
			StateOld = SegQOld->GetState(stateName);
			if (StateNew->GetFreedom() != frozen
			&& StateOld->GetFreedom() != frozen) {
				if (gradients == 1) {
					int Min = LatOld->GetTotalNumLayers();
					CopyInitialGuess(StateOld->GetSWF(),StateNew->GetSWF(),Min);
				} else if (gradients == 2) {
					int MXin = LatOld->GetNumLayers(1);
					int MYin = LatOld->GetNumLayers(2);
					CopyInitialGuess(StateOld->GetSWF(),StateNew->GetSWF(),MXin,MYin);
				}else if (gradients == 3) {
					int MXin = LatOld->GetNumLayers(1);
					int MYin = LatOld->GetNumLayers(2);
					int MZin = LatOld->GetNumLayers(3);
					CopyInitialGuess(StateOld->GetSWF(),StateNew->GetSWF(),MXin,MYin,MZin);
				}
			}
		}
	}
	for (i=1; i<=numLiquidStates; i++) {
		StateNew = StateQ[i];
		Vector G = StateNew->GetSWF();
		for (int z=1; z<=M; z++) {
			SWF[i][z] = G[z];
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
	if (gradients == 3 && charged) {
		int MXin = LatOld->GetNumLayers(1);
		int MYin = LatOld->GetNumLayers(2);
		int MZin = LatOld->GetNumLayers(3);
		CopyInitialGuess(SegQOld->GetElectricPotential(),psi0,MXin,MYin,MZin);
	}

	if (charged && diffEpsilon ) {
		Lat->SetBoundaries(psi0);
		double preFactor = (EPS0*Lat->GetSiteDistance())
					/(BOLTZMANN*TEMPERATURE);
		Lat->ElectricFieldSquared(ESquared,psi0,preFactor);
	}

	(MolQ->GetSolvent())->SetPhiBulk((MolQOld->GetSolvent())->GetPhiBulk());
	phiBulkSolvent = (MolQ->GetSolvent())->GetPhiBulk();
	for (i=1; i<=numStates; i++) {
		StateNew = SegQ->GetState(i);
		stateName = StateNew->GetName();
		if (SegQOld->StateDefined(stateName)) {
			StateOld = SegQOld->GetState(stateName);
			StateNew->SetAlphaBulk(StateOld->GetAlphaBulk());
		}
	}
	if (neutralizerPresent) {
		SF_Molecule* NeutralizerOld = MolQOld->GetNeutralizer();
		if (NeutralizerOld != NULL) {
			(MolQ->GetNeutralizer())->SetPhiBulk(NeutralizerOld->GetPhiBulk());
			phiBulkNeutralizer = (MolQ->GetNeutralizer())->GetPhiBulk();
		} else {
			phiBulkNeutralizer = 0;
		}
	}
	if (multiState) {
		ReactionQ->ComputeAlphaBulk();
	}
	if (compute) {
		SegQ->UpdateSWF();
		SegQ->UpdatePhiBulk();
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
}

void
SF_Solve::NoInitialGuess() {
	int i,z;
	for (z=1; z<=M; z++) {
		if (SegItVar(z)) {
			for (i=1; i<=numLiquidStates; i++) {
				SWF[i][z] = 1;
			}
		}
		if (PotItVar(z)) {
			psi0[z] = 0;
		}
	}
	//SegQ->UpdateSWF();
	//UpdateItVar(Offset);
}

void
SF_Solve::ComputeInitialGuess() {
	int i,z;
	for (z=1; z<=M; z++) {
		if (SegItVar(z)) {
			for (i=1; i<=numLiquidStates; i++) {
				SF_State* State=StateQ[i];
				double u = SegQ->ChemInt(State,z);
				u -= SegQ->ChemIntBulk(State);
				if (State->ExternalPotential()) {
					u += (State->GetExternalPotential())[z];
				}
				SWF[i][z] = exp(-u);
			}
		}
		if (PotItVar(z)) {
			psi0[z] = 0;
		}
	}
	UpdateItVar(false);
	SegQ->UpdateSWF();
}
void
SF_Solve::ComputeAlphaBulk() {
	if (!multiState) return;
	int i,z;
	SF_Molecule* Mol;
	if (neutralizerPresent) {
		double bulkCharge = 0;
		for (i=1; i<=MolQ->GetNumMolecules(); i++) {
			Mol = MolQ->GetMolecule(i);
			if (Mol != MolQ->GetSolvent() && Mol!= MolQ->GetNeutralizer()) {
				bulkCharge += Mol->GetBulkCharge()*Mol->GetPhiBulk()/Mol->GetChainLength();
			}
		}
		Mol = MolQ->GetNeutralizer();
		phiBulkNeutralizer = -bulkCharge*Mol->GetChainLength()/Mol->GetBulkCharge();
		if (phiBulkNeutralizer <= 0) phiBulkNeutralizer = 1e-9;
		if (phiBulkNeutralizer >= 1) phiBulkNeutralizer = 1-1e-9;
		Mol->SetPhiBulk(phiBulkNeutralizer);
	}
	phiBulkSolvent = 1;
	for (i=1; i<=MolQ->GetNumMolecules(); i++) {
		Mol = MolQ->GetMolecule(i);
		if (Mol != MolQ->GetSolvent()) {
			phiBulkSolvent -= Mol->GetPhiBulk();
		}
	}
	if (phiBulkSolvent <= 0) phiBulkSolvent = 1e-9;
	if (phiBulkSolvent >= 1) phiBulkSolvent = 1-1e-9;
	Mol = MolQ->GetSolvent();
	Mol->SetPhiBulk(phiBulkSolvent);
	SegQ->UpdatePhiBulk();
	if (neutralizerPresent) {
		phiBulkNeutralizer = (MolQ->GetNeutralizer())->GetPhiBulk();
	}
	SF_Segment* Segment;
	double dummy;
	for (i=1; i<=SegQ->GetNumSegments(); i++) {
		Segment = SegQ->GetSegment(i);
		int numStates = Segment->GetNumStates();
		if (numStates > 1) {
			for (z=1; z<=numStates; z++) {
				dummy = 1.0/numStates + 1e-9;
				if (z == 2) dummy -= 5e-9;
				(Segment->GetState(z))->SetAlphaBulk(dummy);
			}
		}
	}
	UpdateSWF();
	if (neutralizerPresent) {
		MolQ->ComputePhi(phiBulkSolvent,phiBulkNeutralizer);
	} else {
		MolQ->ComputePhi(phiBulkSolvent);
	}
	SegQ->UpdatePhiBulk();
	SegQ->UpdatePhiStates();
	ReactionQ->ComputeAlphaBulk();
	phiBulkSolvent = 1;
	for (i=1; i<=MolQ->GetNumMolecules(); i++) {
		Mol = MolQ->GetMolecule(i);
		if (Mol != MolQ->GetSolvent()) {
			phiBulkSolvent -= Mol->GetPhiBulk();
		}
	}
	if (phiBulkSolvent <= 0) phiBulkSolvent = 1e-9;
	if (phiBulkSolvent >= 1) phiBulkSolvent = 1-1e-9;
	phiBulkSolvent += 1e-10;
	Mol = MolQ->GetSolvent();
	Mol->SetPhiBulk(phiBulkSolvent);
	SegQ->UpdatePhiBulk();
}
void
SF_Solve::SetExternalPotentials() {
	Infile in(externalPotentialFile);
	if (!in.Open(Blanks(100))) {
		Message(fatal,"external_potential_file : "
			+ externalPotentialFile + " cannot be opened");
	}
	in.Inimage();
	if (in.Endfile()) {
		Message(fatal, "external_potential_file : "
			+ externalPotentialFile + " is empty");
	}
	Text line;
	line = Copy(in.Image.Strip());
	if (*line != *Copy("gradients")) {
		Message(fatal,"error reading '" + externalPotentialFile
			+ "' unrecognized line: " + line);
	}
	in.Inimage();
	line = Copy(in.Image.Strip());
	int gradients = line.Getint();
	int Min=0;
	int MXin=0;
	int MYin=0;
	int MZin=0;
	if (gradients == 1) {
		in.Inimage();
		line = Copy(in.Image.Strip());
		Min = line.Getint();
	}
	if (gradients == 2) {
		in.Inimage();
		line = Copy(in.Image.Strip());
		MXin = line.Getint();
		in.Inimage();
		line = Copy(in.Image.Strip());
		MYin = line.Getint();
		Min = MXin*MYin;
	}
	if (gradients == 3) {
		in.Inimage();
		line = Copy(in.Image.Strip());
		MXin = line.Getint();
		in.Inimage();
		line = Copy(in.Image.Strip());
		MYin = line.Getint();
		in.Inimage();
		line = Copy(in.Image.Strip());
		MZin = line.Getint();
		Min = MXin*MYin*MZin;
	}
	Text stateName;
	if (!in.Endfile()) in.Inimage();
	line = Copy(in.Image.Strip());
	while (*line == *Copy("molecule")) {
		in.Inimage();
		line = Copy(in.Image.Strip());
		in.Inimage();
		line = Copy(in.Image.Strip());
		in.Inimage();
		stateName = Copy(in.Image.Strip());
		Vector guess;
		guess.Dim(1,Min);
		for (int z=1; z<=Min; z++) {
			in.Inimage();
			line = Copy(in.Image.Strip());
			guess[z] = line.Getreal();
		}
		if (SegQ->StateDefined(stateName)) {
			SF_State* State = SegQ->GetState(stateName);
			State->SetExternalPotential(guess);
		}
		if (!in.Endfile()) in.Inimage();
		line = Copy(in.Image.Strip());
	}
	in.Close();
}
void
SF_Solve::SetFrozenSWFZero() {
	for (int z=1; z<=M; z++) {
		if (!SegItVar(z)) {
			int count = 0;
			for (int i=1; i<=numDiffItVar; i++) {
				for (int k=1; k<=NumStatesForItVar[i]; k++) {
					count++;
					SWF[count][z] = 0;
				}
			}
		}
	}
}

void
SF_Solve::SetUOffset(){
	int numSeg = SegQ->GetNumSegments();
	int numMol = MolQ->GetNumMolecules();
	Vector phi;
	for (int i=1; i<=M; i++) U_offset[i]=0.0;
	for (int i=1; i<=numSeg; i++) {
		SF_Segment* Seg = SegQ->GetSegment(i);
		if (Seg->GetFreedom()==frozen) {
			phi = Seg->GetPhi();
			for (int z=1; z<=M; z++) U_offset[z]+=phi[z];
			Offset=true;
		}
	}
	for (int i=1; i<=numMol; i++) {
		SF_Molecule* Mol=MolQ->GetMolecule(i);
		if (Mol->GetMolStructure()->SomeSegmentsPinned() ) {
			double theta=Mol ->GetTheta();
			double N= 1.0*Mol->GetChainLength();
			SF_MolStructure* Chain = Mol->GetMolStructure();
			for (int k=1; k<=Chain->GetNumDiffSegments(); k++) {
				SF_MolSegment* MolSeg=Chain->GetDiffSegment(k);
				if (MolSeg->GetFreedom() == pinned) {
					LatticeRange* LatRange=MolSeg->GetLatRange();
					double Npos = 1.0*LatRange->GetNumPos();
					double avLengthPart = Chain->GetAvNumSegments(MolSeg);
					double phi_av=0;
					if (Npos != 0)  {
						phi_av=theta/N*avLengthPart/Npos;
						if (phi_av > 1) {
							    cout << " Number Pinned  = " << avLengthPart << endl;
								cout << " N = " << N << endl;
								cout << " Npos = " << Npos << endl;
								cout << " theta = "  << theta << endl;
								Message(literal,"pinned segments possibly have not enough space");}
						for (int z=1; z<=M; z++) if (LatRange -> InRange(z)) U_offset[z]+=phi_av;
						Offset=true;
					}
				}

			}
		}
	}
	if (Offset) {
		for (int z=1; z<=M; z++) {
			if (U_offset[z]>0.00001 && U_offset[z]<0.9999) U_offset[z]=-log(1-U_offset[z]);
		}
	}
}


void
SF_Solve::CopyInitialGuess(const Vector in,
						   Vector o,
						   const int Min) const {
	if (Lat->GetNumGradients() == 1) {
		switch(initGuessSym) {
			case firstToLast:
				Copy1Dto1DFirstToLast(in,o,Min);
				break;
			case boundsToMiddle:
				Copy1Dto1DBoundsToMiddle(in,o,Min);
				break;
			default:
				Message(fatal,"Programming error 1 in SF_Solve::CopyInitialGuess");
				break;
		}
	} else if (Lat->GetNumGradients() == 2) {
		switch(initGuessSym) {
			case firstToLast:
				Copy1Dto2DFirstToLast(in,o,Min);
				break;
			case boundsToMiddle:
				Copy1Dto2DBoundsToMiddle(in,o,Min);
				break;
			default:
				Message(fatal,"Programming error 2 in SF_Solve::CopyInitialGuess");
				break;
		}
	} else if (Lat->GetNumGradients() == 3) {
		switch(initGuessSym) {
			case firstToLast:
				Copy1Dto3DFirstToLast(in,o,Min);
			break;
			case boundsToMiddle:
				Message(fatal,"Programming error 2a in SF_Solve::CopyInitialGuess");
				break;
			default:
				Copy1Dto3DFirstToLast(in,o,Min);
				break;
		}
	}
}
void
SF_Solve::CopyInitialGuess(const Vector in,
						   Vector o,
						   const int MXin,
						   const int MYin
						   ) const {
	if (Lat->GetNumGradients() == 1) {
		switch(initGuessSym) {
			case firstToLast:
				Copy2Dto1DFirstToLast(in,o,MXin,MYin);
				break;
			case boundsToMiddle:
				Copy2Dto1DBoundsToMiddle(in,o,MXin,MYin);
				break;
			default:
				Message(fatal,"Programming error 3 in SF_Solve::CopyInitialGuess");
				break;
		}
	} else if (Lat->GetNumGradients() == 2) {
		switch(initGuessSym) {
			case firstToLast:
				Copy2Dto2DFirstToLast(in,o,MXin,MYin);
				break;
			case boundsToMiddle:
				Copy2Dto2DBoundsToMiddle(in,o,MXin,MYin);
				break;
			default:
				Message(fatal,"Programming error 4 in SF_Solve::CopyInitialGuess");
				break;
		}
	} else if (Lat->GetNumGradients() == 3) {
		switch(initGuessSym) {
			case firstToLast:
				Message(fatal,"Programming error 4a in SF_Solve::CopyInitialGuess");
				break;
			case boundsToMiddle:
				Message(fatal,"Programming error 4a in SF_Solve::CopyInitialGuess");
				break;
			default:
				Message(fatal,"Programming error 4a in SF_Solve::CopyInitialGuess");
				break;
		}
	}
}

void
SF_Solve::CopyInitialGuess(const Vector in,
						   Vector o,
						   const int MXin,
						   const int MYin,
						   const int MZin
						   ) const {
	if (Lat->GetNumGradients() == 1) {
		switch(initGuessSym) {
			case firstToLast:
				Message(fatal,"Programming error 5 in SF_Solve::CopyInitialGuess");
				break;
			case boundsToMiddle:
				Message(fatal,"Programming error 5 in SF_Solve::CopyInitialGuess");
				break;
			default:
				Message(fatal,"Programming error 5 in SF_Solve::CopyInitialGuess");
				break;
		}
	} else if (Lat->GetNumGradients() == 2) {
		switch(initGuessSym) {
			case firstToLast:
				Message(fatal,"Programming error 5 in SF_Solve::CopyInitialGuess");
				break;
			case boundsToMiddle:
				Message(fatal,"Programming error 5 in SF_Solve::CopyInitialGuess");
				break;
			default:
				Message(fatal,"Programming error 5 in SF_Solve::CopyInitialGuess");
				break;
		}
	}	else if (Lat->GetNumGradients() == 3) {
			switch(initGuessSym) {
				case firstToLast:
					Copy3Dto3DFirstToLast(in,o,MXin,MYin,MZin);
					break;
				case boundsToMiddle:
					Message(fatal,"Programming error 5 in SF_Solve::CopyInitialGuess");
					break;
				default:
					Message(fatal,"Programming error 5 in SF_Solve::CopyInitialGuess");
					break;
			}
		}
}

void
SF_Solve::Copy1Dto1DFirstToLast(const Vector in,
								Vector o,
								const int Min) const {
	int Mout = Lat->GetTotalNumLayers();
	for (int z=2; z<=Min; z++) {
		if (z < Mout) o[z] = in[z];
	}
}
void
SF_Solve::Copy2Dto1DFirstToLast(const Vector in,
								Vector o,
								const int MXin,
								const int MYin) const {
	int Mout = Lat->GetTotalNumLayers();
	int zIn;
	for (int xGrad=1; xGrad<=MXin; xGrad++) {
		for (int yGrad=1; yGrad<=MYin; yGrad++) {
			if (xGrad==2 && yGrad<Mout) {
				zIn = xGrad*(MYin-1) + yGrad;
				o[yGrad] = in[zIn];
			}
		}
	}
}
void
SF_Solve::Copy2Dto2DFirstToLast(const Vector in,
								Vector o,
								const int MXin,
								const int MYin) const {
	int MXout = Lat->GetNumLayers(1);
	int MYout = Lat->GetNumLayers(2);
	int zIn, zOut;
	for (int xGrad=1; xGrad<=MXin; xGrad++) {
		for (int yGrad=1; yGrad<=MYin; yGrad++) {
			if (xGrad<MXout && yGrad<MYout) {
				zIn = (xGrad-1)*MYin + yGrad;
				zOut = (xGrad-1)*MYout + yGrad;
				//zIn = xGrad*(MYin-1) + yGrad; // frans: I do not understand this, I expect it does not work
				//zOut = xGrad*(MYout-1) + yGrad;
				o[zOut] = in[zIn];
			}
		}
	}
}
void
SF_Solve::Copy3Dto3DFirstToLast(const Vector in,
								Vector o,
								const int MXin,
								const int MYin,
								const int MZin) const {
	int MXout = Lat->GetNumLayers(1);
	int MYout = Lat->GetNumLayers(2);
	int MZout = Lat->GetNumLayers(3);
	int zIn, zOut;

	for (int xGrad=1; xGrad<=MXin; xGrad++) {
		for (int yGrad=1; yGrad<=MYin; yGrad++) {
			for (int zGrad=1; zGrad<=MZin; zGrad++){
				if (xGrad<MXout && yGrad<MYout && zGrad<MZout) {
					zIn = (xGrad-1)*MYin*MZin + (yGrad-1)*MZin+ zGrad;
					zOut = (xGrad-1)*MYout*MZout + (yGrad-1)*MZout + zGrad;
					o[zOut] = in[zIn];
				}
			}
		}
	}
}
void
SF_Solve::Copy1Dto2DFirstToLast(const Vector in,
								Vector o,
								const int Min) const {
	int MXout = Lat->GetNumLayers(1);
	int MYout = Lat->GetNumLayers(2);
	int zOut;
	for (int yGrad=1; yGrad<=Min; yGrad++) {
		for (int xGrad=1; xGrad<=MXout; xGrad++) {
			if (yGrad<MYout) {
				zOut = (xGrad-1)*MYout + yGrad;
				o[zOut] = in[yGrad];
			}
		}
	}
}

void
SF_Solve::Copy1Dto3DFirstToLast(const Vector in, Vector o, const int Min) const {
	int MXout = Lat->GetNumLayers(1);
	int MYout = Lat->GetNumLayers(2);
	int MZout = Lat->GetNumLayers(3);
	int jx = MYout*MZout;
	int jy = MZout;
	//int jz = 1;

	int zOut;
	for (int zGrad=1; zGrad<=Min; zGrad++) {
		for (int yGrad=1; yGrad<=MYout; yGrad++) {
			for (int xGrad=1; xGrad<=MXout; xGrad++) {
				if (zGrad<=MZout) {
					zOut = jx*(xGrad-1)+jy*(yGrad-1) + zGrad;
					o[zOut] = in[zGrad];
				}
			}
		}
	}
}

void
SF_Solve::Copy1Dto1DBoundsToMiddle(const Vector in,
					 Vector o,
					 const int Min) const {
	int Mout = Lat->GetTotalNumLayers();
	int z;
	for (z=2; z<=Min/2; z++) {
		if (z < Mout/2) o[z] = in[z];
	}
	int zOut = Mout;
	for (z=Min-1; z>Min/2; z--) {
		zOut--;
		if (zOut > Mout/2) o[zOut] = in[z];
	}
}
void
SF_Solve::Copy2Dto1DBoundsToMiddle(const Vector in,
					 Vector o,
					 const int MXin,
					 const int MYin) const {
	int Mout = Lat->GetTotalNumLayers();
	int zIn;
	for (int xGrad=1; xGrad<=MXin; xGrad++) {
		for (int yGrad=1; yGrad<=MYin; yGrad++) {
			if (xGrad==2 && yGrad<Mout) {
				zIn = xGrad*(MYin-1) + yGrad;
				o[yGrad] = in[zIn];
			}
		}
	}
}
void
SF_Solve::Copy2Dto2DBoundsToMiddle(const Vector in,
					 Vector o,
					 const int MXin,
					 const int MYin) const {
	int MXout = Lat->GetNumLayers(1);
	int MYout = Lat->GetNumLayers(2);
	int zIn, zOut;
	for (int xGrad=1; xGrad<=MXin; xGrad++) {
		for (int yGrad=1; yGrad<=MYin; yGrad++) {
			if (xGrad<MXout && yGrad<MYout) {
				zIn = xGrad*(MYin-1) + yGrad;
				zOut = xGrad*(MYout-1) + yGrad;
				o[zOut] = in[zIn];
			}
		}
	}
}
void
SF_Solve::Copy1Dto2DBoundsToMiddle(const Vector in,
					 Vector o,
					 const int Min) const {
	int MXout = Lat->GetNumLayers(1);
	int MYout = Lat->GetNumLayers(2);
	int zOut;
	for (int yGrad=1; yGrad<=Min; yGrad++) {
		for (int xGrad=1; xGrad<=MXout; xGrad++) {
			if (yGrad<MYout) {
				zOut = xGrad*(MYout-1) + yGrad;
				o[zOut] = in[yGrad];
			}
		}
	}
}
void
SF_Solve::WriteOverflowError() const {
	Message(literal,"Some problem appears: \n"
		"Hmm, seems like the equations produced an overflow");
}
//Boolean
//SF_Solve::SegItVar(const int z) const {
//	if (!Lat->WithinBoundaries(z)) return false;
//	SF_State* State;
//	LatticeRange* LatRange;
//	int numStates = SegQ->GetNumStates();
//	for (int i=1; i<=numStates; i++) {
//		State = SegQ->GetState(i);
//		if (State->GetFreedom() == frozen) {
//			LatRange = State->GetLatRange();
//			if (LatRange->InRange(z)) {
//				return false;
//			}
//		}
//	}
//	return true;
//}

/* This function determines whether a certain position in the
   lattice should have an iteration variable for the electric potential.
   When the system is charged, all positions within boundaries or with
   frozen segments should have a potential iteration variable.
*/
Boolean
SF_Solve::PotItVar(const int z) const {
	if (!charged) return false;
	int i;
	SF_Segment *Segment;
	for (i=SegQ->GetNumSegments(); i>0; i--) {
		Segment = SegQ->GetSegment(i);
		if (Segment->GetFreedom() == frozen) {
			if (Segment->GetLatRange()->InRange(z)) {
				return potentialItVar;
			}
		}
	}
	if (Lat->WithinBoundaries(z)) {
		return potentialItVar;
	}
	else return false;
}
void
SF_Solve::ProcessInput() {
	int numNames = MyInput->GetNumNames("newton",0,1);
	if (numNames == 1) {
		Array<Text> names(1,1);
		names = MyInput->GetNames("newton");
		name = names[1];
		Array<Text> param(1,46);
		param[1] = "solver";
		param[2] = "iterationlimit";
		param[3] = "tolerance";
		param[4] = "deltamin";
		param[5] = "deltamax";
		param[6] = "linesearchlimit";
		param[7] = "d_info";
		param[8] = "e_info";
		param[9] = "g_info";
		param[10] = "h_info";
		param[11] = "s_info";
		param[12] = "x_info";
		param[13] = "i_info";
		param[14] = "history_check";
		param[15] = "initial_guess";
		param[16] = "initial_guess_symmetry";
		param[17] = "initial_guess_input_file";
		param[18] = "initial_guess_output_file";
		param[19] = "small_alpha";
		param[20] = "max_n_small_alpha";
		param[21] = "method";
		param[22] = "linetolerance";
		param[23] = "reverse_direction_range";
		param[24] = "max_fr_reverse_direction";
		param[25] = "min_accuracy_for_hessian";
		param[26] = "n_iterations_for_hessian";
		param[27] = "transfer_hessian";
		param[28] = "external_potential_file";
		param[29] = "reset_hessian_criterion";
		param[30] = "no_hessian_scaling";
		param[31] = "ignore_newton_direction";
		param[32] = "max_accuracy_for_hessian_scaling";
		param[33] = "write_initial_guess";
		param[34] = "number_of_stored_vectors";
		param[35] = "m";
		param[36] = "reset_limit_gradient_inversion";
		param[37] = "print_exportable_info";
		param[38] = "print_iteration_info";
		param[39] = "print_common_info";
		param[40] = "print_verbose_info";
		param[41] = "print_improvement_info";
		param[42] = "linesearch";
		//param[43] = "n_iterations_CG";
		param[43] = "reset_CG";
		param[44] = "DIIS";
		param[45] = "numchargeiter";
		param[46] = "no_stars";
		MyInput->CheckParameterNames("newton",name,param);

		setiterationlimit(MyInput->GetInt("newton",name,"iterationlimit",0,INT_MAX,5000));
		settolerance(MyInput->GetReal("newton",name,"tolerance",1e-15,1e-1,1e-7));

		d_info = MyInput->GetBoolean("newton",name,"d_info",false);
		e_info = MyInput->GetBoolean("newton",name,"e_info",!silentForced);
		g_info = MyInput->GetBoolean("newton",name,"g_info",false);
		h_info = MyInput->GetBoolean("newton",name,"h_info",false);
		s_info = MyInput->GetBoolean("newton",name,"s_info",false);
		x_info = MyInput->GetBoolean("newton",name,"x_info",false);
		i_info = MyInput->GetBoolean("newton",name,"i_info",false);
		no_stars = MyInput->GetBoolean("newton",name,"no_stars",false);

		print_iteration_info = MyInput->GetInt("newton",name,"print_iteration_info",0,INT_MAX,1);

		Array<Text> method(1,9);
		method[1] = "hessian";
		method[2] = "pseudohessian";
		method[3] = "LBFGS";
		method[4] = "CG";
		method[5] = "SD";
		method[6] = "TN";
		method[7] = "BRR";
		method[8] = "Pure_DIIS";
		method[9] = "experimental";

		Array<Text> line_method(1,4);
		line_method[1] = "Scheutjens";
		line_method[2] = "CG";
		line_method[3] = "BB";
		line_method[4] = "secant";


		int Newtonmethod=2;
		if (MyInput->ValueSet("newton",name,"method")){
			Newtonmethod = MyInput->GetChoice("newton",name,"method",method,2);
		};

		int Linemethod = 1;
		if (MyInput->ValueSet("newton",name,"linesearch")){
			Linemethod = MyInput->GetChoice("newton",name,"linesearch",line_method,1);
		};

		switch (Linemethod) {
		case(1) :
			Scheutjens=true;
			alpha_cg = false;
			secant = false;
			No_LS = false;
			setlinetolerance(MyInput->GetReal("newton",name,"linetolerance",0,1,0.9));
			setlinesearchlimit(MyInput->GetInt("newton",name,"linesearchlimit",1,1000,getlinesearchlimit()));
			setdeltamax(MyInput->GetReal("newton",name,"deltamax",1e-15,DBL_MAX,0.5));
			setdeltamin(MyInput->GetReal("newton",name,"deltamin",0,1,0));
			break;
		case(2) :
			Scheutjens=false;
			alpha_cg = true;
			secant = false;
			No_LS =false;
			setlinetolerance(MyInput->GetReal("newton",name,"linetolerance",0,1,0.9));
			setlinesearchlimit(MyInput->GetInt("newton",name,"linesearchlimit",1,1000,getlinesearchlimit()));
			setdeltamax(MyInput->GetReal("newton",name,"deltamax",1e-15,DBL_MAX,0.5));
			setdeltamin(MyInput->GetReal("newton",name,"deltamin",0,1,0));
			break;
		case(3) :
			Scheutjens=false;
			alpha_cg = false;
			secant = false;
			No_LS = true;
			setdeltamax(MyInput->GetReal("newton",name,"deltamax",1e-15,DBL_MAX,0.5));
			setdeltamin(MyInput->GetReal("newton",name,"deltamin",0,1,0));
			break;
		case(4) :
			Scheutjens=false;
			alpha_cg = false;
			No_LS = false;
			secant = true;
			setdeltamax(MyInput->GetReal("newton",name,"deltamax",1e-15,DBL_MAX,0.5));
			setdeltamin(MyInput->GetReal("newton",name,"deltamin",0,1,0));
			setlinetolerance(MyInput->GetReal("newton",name,"linetolerance",0,1,1e-5));
			setlinesearchlimit(MyInput->GetInt("newton",name,"linesearchlimit",1,10,getlinesearchlimit()));
			break;
		}

		switch (Newtonmethod) {
		case(1) :
			pseudohessian = false; SD=false; TN=false; BRR = false;Pure_DIIS = false;
			setDIIS(MyInput->GetInt("newton",name,"DIIS",0,64,0));
			MyInput->DontCombineParam("newton",name,"method","number_of_stored_vectors");
			MyInput->DontCombineParam("newton",name,"method","m");
			MyInput->DontCombineParam("newton",name,"method","print_common_info");
			MyInput->DontCombineParam("newton",name,"method","print_verbose_info");
			MyInput->DontCombineParam("newton",name,"method","print_improvement_info");
			MyInput->DontCombineParam("newton",name,"method","print_exportable_info");
			MyInput->DontCombineParam("newton",name,"method","reset_limit_gradient_inversion");
			reverseDirectionRange = MyInput->GetInt("newton",name,"reverse_direction_range",2,INT_MAX,50);
			maxFrReverseDirection = MyInput->GetReal("newton",name,"max_fr_reverse_direction",0,1,0.4);
			transferHessian = MyInput->GetBoolean("newton",name,"transfer_hessian",false);
			MyInput->DontCombineParam("newton",name,"method","n_iterations_for_hessian");
			MyInput->DontCombineParam("newton",name,"method","min_accuracy_for_hessian");
			setlinetolerance(MyInput->GetReal("newton",name,"linetolerance",0,1,0.9));
			setlinesearchlimit(MyInput->GetInt("newton",name,"linesearchlimit",1,1000,getlinesearchlimit()));
			//setdeltamax(MyInput->GetReal("newton",name,"deltamax",1e-15,DBL_MAX,0.5));
			//setdeltamin(MyInput->GetReal("newton",name,"deltamin",0,1,0));
			historyTest = MyInput->GetBoolean("newton",name,"history_check",false);
			break;
		case(2):
			pseudohessian = true; SD=false; TN=false; BRR = false;Pure_DIIS = false;
			setDIIS(MyInput->GetInt("newton",name,"DIIS",0,64,0));
			historyTest = MyInput->GetBoolean("newton",name,"history_check",false);
			setlinetolerance(MyInput->GetReal("newton",name,"linetolerance",0,1,0.9));
			//setdeltamax(MyInput->GetReal("newton",name,"deltamax",1e-15,DBL_MAX,0.5));
			//setdeltamin(MyInput->GetReal("newton",name,"deltamin",0,getdeltamax(),0));
			setlinesearchlimit(MyInput->GetInt("newton",name,"linesearchlimit",1,1000,getlinesearchlimit()));
			if (MyInput->ValueSet("newton",name, "n_iterations_for_hessian")) {
				minAccuracyForHessian = MyInput->GetReal("newton",name,"min_accuracy_for_hessian",0,1,0.1);
			} else {
				minAccuracyForHessian = MyInput->GetReal("newton",name,"min_accuracy_for_hessian",0,1,0);
			}
			MyInput->DontCombineParam("newton",name,"method","number_of_stored_vectors");
			MyInput->DontCombineParam("newton",name,"method","m");
			MyInput->DontCombineParam("newton",name,"method","print_common_info");
			MyInput->DontCombineParam("newton",name,"method","print_verbose_info");
			MyInput->DontCombineParam("newton",name,"method","print_improvement_info");
			MyInput->DontCombineParam("newton",name,"method","print_exportable_info");
			MyInput->DontCombineParam("newton",name,"method","reset_limit_gradient_inversion");
			reverseDirectionRange = MyInput->GetInt("newton",name,"reverse_direction_range",2,INT_MAX,50);
			maxFrReverseDirection = MyInput->GetReal("newton",name,"max_fr_reverse_direction",0,1,0.4);
			transferHessian = MyInput->GetBoolean("newton",name,"transfer_hessian",false);
			numIterationsForHessian = MyInput->GetInt("newton",name,"n_iterations_for_hessian",1,INT_MAX,getiterationlimit()+100);
			break;
		case(3) :
			LM_BFGS = true; pseudohessian = false;  CG=false; SD=false; TN = false; BRR = false; Pure_DIIS = false;
			setDIIS(MyInput->GetInt("newton",name,"DIIS",0,64,0));
			historyTest = MyInput->GetBoolean("newton",name,"history_check",false);
			amount_of_stored_vectors = MyInput->GetInt("newton",name,"m",0,64,0);
			if (amount_of_stored_vectors==0) amount_of_stored_vectors = MyInput->GetInt("newton",name,"number_of_stored_vectors",1,64,12);
		       reverseDirectionRange = MyInput->GetInt("newton",name,"reverse_direction_range",2,INT_MAX,50);
		       MyInput->DontCombineParam("newton",name,"method","min_accuracy_for_hessian");
			MyInput->DontCombineParam("newton",name,"method","transfer_hessian");
			MyInput->DontCombineParam("newton",name,"method","print_common_info");
			MyInput->DontCombineParam("newton",name,"method","print_verbose_info");
			MyInput->DontCombineParam("newton",name,"method","print_improvement_info");
			MyInput->DontCombineParam("newton",name,"method","print_exportable_info");
			MyInput->DontCombineParam("newton",name,"method","reset_limit_gradient_inversion");
			MyInput->DontCombineParam("newton",name,"method","max_fr_reverse_direction");
			MyInput->DontCombineParam("newton",name,"method","n_iterations_for_hessian");
			break;
		case(4) :
			LM_BFGS = true; pseudohessian = false;  CG=true;  SD =false; TN = false;BRR = false; Pure_DIIS = false;
			setDIIS(MyInput->GetInt("newton",name,"DIIS",0,64,0));
			historyTest = MyInput->GetBoolean("newton",name,"history_check",false);
			//setdeltamax(MyInput->GetReal("newton",name,"deltamax",1e-15,DBL_MAX,0.5));
			//setdeltamin(MyInput->GetReal("newton",name,"deltamin",0,getdeltamax(),0));
			//setlinetolerance(MyInput->GetReal("newton",name,"linetolerance",0,1,0.9));
			setlinesearchlimit(MyInput->GetInt("newton",name,"linesearchlimit",1,1000,getlinesearchlimit()));
			amount_of_stored_vectors = MyInput->GetInt("newton",name,"m",1,64,0);
			if (amount_of_stored_vectors==0) amount_of_stored_vectors = MyInput->GetInt("newton",name,"number_of_stored_vectors",1,64,1);
		     	reverseDirectionRange = MyInput->GetInt("newton",name,"reverse_direction_range",2,INT_MAX,50);
		     	MyInput->DontCombineParam("newton",name,"method","min_accuracy_for_hessian");
			MyInput->DontCombineParam("newton",name,"method","transfer_hessian");
			MyInput->DontCombineParam("newton",name,"method","print_common_info");
			MyInput->DontCombineParam("newton",name,"method","print_verbose_info");
			MyInput->DontCombineParam("newton",name,"method","print_improvement_info");
			MyInput->DontCombineParam("newton",name,"method","print_exportable_info");
			MyInput->DontCombineParam("newton",name,"method","reset_limit_gradient_inversion");
			MyInput->DontCombineParam("newton",name,"method","max_fr_reverse_direction");
			MyInput->DontCombineParam("newton",name,"method","n_iterations_for_hessian");
			reset_CG = MyInput->GetInt("newton",name,"reset_CG",1,INT_MAX,10);
			break;
		case(5) :
			LM_BFGS = true; pseudohessian = false;  CG=false;  SD =true; TN = false;BRR = false; Pure_DIIS = false;
			setDIIS(MyInput->GetInt("newton",name,"DIIS",0,64,0));
			historyTest = MyInput->GetBoolean("newton",name,"history_check",false);
			setdeltamax(MyInput->GetReal("newton",name,"deltamax",1e-15,DBL_MAX,0.5));
			setdeltamin(MyInput->GetReal("newton",name,"deltamin",0,getdeltamax(),getdeltamax()/100000));
			//setlinetolerance(MyInput->GetReal("newton",name,"linetolerance",0,1,0.9));
			setlinesearchlimit(MyInput->GetInt("newton",name,"linesearchlimit",1,1000,getlinesearchlimit()));
			amount_of_stored_vectors = MyInput->GetInt("newton",name,"m",1,64,0);
			if (amount_of_stored_vectors==0) amount_of_stored_vectors = MyInput->GetInt("newton",name,"number_of_stored_vectors",1,64,1);
			reverseDirectionRange = MyInput->GetInt("newton",name,"reverse_direction_range",2,INT_MAX,50);
			MyInput->DontCombineParam("newton",name,"method","min_accuracy_for_hessian");
			MyInput->DontCombineParam("newton",name,"method","transfer_hessian");
			MyInput->DontCombineParam("newton",name,"method","print_common_info");
			MyInput->DontCombineParam("newton",name,"method","print_verbose_info");
			MyInput->DontCombineParam("newton",name,"method","print_improvement_info");
			MyInput->DontCombineParam("newton",name,"method","print_exportable_info");
			MyInput->DontCombineParam("newton",name,"method","reset_limit_gradient_inversion");
			MyInput->DontCombineParam("newton",name,"method","max_fr_reverse_direction");
			MyInput->DontCombineParam("newton",name,"method","n_iterations_for_hessian");
			break;
		case(6) :
			LM_BFGS = true; pseudohessian = false;  CG=false; SD=false; TN = true;BRR = false;Pure_DIIS = false;
			setDIIS(MyInput->GetInt("newton",name,"DIIS",0,64,0));
			historyTest = MyInput->GetBoolean("newton",name,"history_check",false);
			//setdeltamax(MyInput->GetReal("newton",name,"deltamax",1e-15,DBL_MAX,0.5));
			//setdeltamin(MyInput->GetReal("newton",name,"deltamin",0,getdeltamax(),0));
			//setlinetolerance(MyInput->GetReal("newton",name,"linetolerance",0,1,0.9));
			amount_of_stored_vectors = MyInput->GetInt("newton",name,"m",0,64,0);
			if (amount_of_stored_vectors==0) amount_of_stored_vectors = MyInput->GetInt("newton",name,"number_of_stored_vectors",1,64,64);
		       reverseDirectionRange = MyInput->GetInt("newton",name,"reverse_direction_range",2,INT_MAX,50);
		       MyInput->DontCombineParam("newton",name,"method","min_accuracy_for_hessian");
			MyInput->DontCombineParam("newton",name,"method","transfer_hessian");
			MyInput->DontCombineParam("newton",name,"method","print_common_info");
			MyInput->DontCombineParam("newton",name,"method","print_verbose_info");
			MyInput->DontCombineParam("newton",name,"method","print_improvement_info");
			MyInput->DontCombineParam("newton",name,"method","print_exportable_info");
			MyInput->DontCombineParam("newton",name,"method","reset_limit_gradient_inversion");
			MyInput->DontCombineParam("newton",name,"method","max_fr_reverse_direction");
			MyInput->DontCombineParam("newton",name,"method","n_iterations_for_hessian");
			//n_iterations_CG = MyInput->GetInt("newton",name,"n_iterations_CG",1,100,10);
			setlinesearchlimit(MyInput->GetInt("newton",name,"linesearchlimit",1,1000,getlinesearchlimit()));
			break;
		case(7) :
			BRR = true; pseudohessian = false;  CG=false; SD=false; TN = false; LM_BFGS=false;Pure_DIIS = false;
			setDIIS(MyInput->GetInt("newton",name,"DIIS",0,64,0));
			historyTest = MyInput->GetBoolean("newton",name,"history_check",false);
			setdeltamax(MyInput->GetReal("newton",name,"deltamax",1e-15,DBL_MAX,0.5));
			setdeltamin(MyInput->GetReal("newton",name,"deltamin",0,getdeltamax(),0));
			//setlinetolerance(MyInput->GetReal("newton",name,"linetolerance",0,1,0.9));
			amount_of_stored_vectors = MyInput->GetInt("newton",name,"m",0,64,0);
			if (amount_of_stored_vectors==0) amount_of_stored_vectors = MyInput->GetInt("newton",name,"number_of_stored_vectors",1,64,64);
		    reverseDirectionRange = MyInput->GetInt("newton",name,"reverse_direction_range",2,INT_MAX,50);
		    MyInput->DontCombineParam("newton",name,"method","min_accuracy_for_hessian");
			MyInput->DontCombineParam("newton",name,"method","transfer_hessian");
			MyInput->DontCombineParam("newton",name,"method","print_common_info");
			MyInput->DontCombineParam("newton",name,"method","print_verbose_info");
			MyInput->DontCombineParam("newton",name,"method","print_improvement_info");
			MyInput->DontCombineParam("newton",name,"method","print_exportable_info");
			MyInput->DontCombineParam("newton",name,"method","reset_limit_gradient_inversion");
			MyInput->DontCombineParam("newton",name,"method","max_fr_reverse_direction");
			MyInput->DontCombineParam("newton",name,"method","n_iterations_for_hessian");
			//n_iterations_CG = MyInput->GetInt("newton",name,"n_iterations_CG",1,100,10);
			setlinesearchlimit(MyInput->GetInt("newton",name,"linesearchlimit",1,1000,getlinesearchlimit()));
			break;
		case(8) :
			Pure_DIIS = true; BRR = false; pseudohessian = false;  CG=false; SD=false; TN = false; LM_BFGS=false;
			e_info = true;
			setDIIS(MyInput->GetInt("newton",name,"DIIS",0,64,0));
			historyTest = MyInput->GetBoolean("newton",name,"history_check",false);
			setdeltamax(MyInput->GetReal("newton",name,"deltamax",1e-15,DBL_MAX,0.1));
			//setdeltamin(MyInput->GetReal("newton",name,"deltamin",0,getdeltamax(),0));
			//setlinetolerance(MyInput->GetReal("newton",name,"linetolerance",0,1,0.9));
			amount_of_stored_vectors = MyInput->GetInt("newton",name,"m",0,64,0);
			if (amount_of_stored_vectors==0) amount_of_stored_vectors = MyInput->GetInt("newton",name,"number_of_stored_vectors",1,64,5);
		    reverseDirectionRange = MyInput->GetInt("newton",name,"reverse_direction_range",2,INT_MAX,50);
		    MyInput->DontCombineParam("newton",name,"method","min_accuracy_for_hessian");
			MyInput->DontCombineParam("newton",name,"method","transfer_hessian");
			MyInput->DontCombineParam("newton",name,"method","print_common_info");
			MyInput->DontCombineParam("newton",name,"method","print_verbose_info");
			MyInput->DontCombineParam("newton",name,"method","print_improvement_info");
			MyInput->DontCombineParam("newton",name,"method","print_exportable_info");
			MyInput->DontCombineParam("newton",name,"method","reset_limit_gradient_inversion");
			MyInput->DontCombineParam("newton",name,"method","max_fr_reverse_direction");
			MyInput->DontCombineParam("newton",name,"method","n_iterations_for_hessian");
			//n_iterations_CG = MyInput->GetInt("newton",name,"n_iterations_CG",1,100,10);
			setlinesearchlimit(MyInput->GetInt("newton",name,"linesearchlimit",1,1000,getlinesearchlimit()));
			break;
		case(9) :
			use_vector_iterations = true;
			MyInput->DontCombineParam("newton",name,"method","history_check");
			MyInput->DontCombineParam("newton",name,"method","deltamax");
			MyInput->DontCombineParam("newton",name,"method","deltamin");
			MyInput->DontCombineParam("newton",name,"method","linetolerance");
			MyInput->DontCombineParam("newton",name,"method","linesearchlimit");
			MyInput->DontCombineParam("newton",name,"method","min_accuracy_for_hessian");
			MyInput->DontCombineParam("newton",name,"method","reverse_direction_range");
			MyInput->DontCombineParam("newton",name,"method","max_fr_reverse_direction");
			MyInput->DontCombineParam("newton",name,"method","transfer_hessian");
			MyInput->DontCombineParam("newton",name,"method","n_iterations_for_hessian");
			amount_of_stored_vectors = MyInput->GetInt("newton",name,"m",0,64,0);
			if (amount_of_stored_vectors==0) amount_of_stored_vectors = MyInput->GetInt("newton",name,"number_of_stored_vectors",1,64,12);
		    	print_common_info = MyInput->GetBoolean("newton",name,"print_common_info",false);
		    	print_verbose_info = MyInput->GetBoolean("newton",name,"print_verbose_info",false);
			print_improvement_info = MyInput->GetBoolean("newton",name,"print_improvement_info",false);
		    	print_exportable_info = MyInput->GetBoolean("newton",name,"print_exportable_info",false);
		    	reset_limit_gradient_inversion = MyInput->GetReal("newton",name,"reset_limit_gradient_inversion",-1.0,1.0,-1.0);
			break;
		default:
			pseudohessian = true;
			break;
		}

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
		smallAlpha = MyInput->GetReal("newton",name,"small_alpha",0,1,1e-5);
		maxNumSmallAlpha = MyInput->GetInt("newton",name,"max_n_small_alpha",1,getiterationlimit(),50);
		if (MyInput->ValueSet("newton",name,"external_potential_file")) {
			externalPotentialFile = MyInput->GetText("newton",name,"external_potential_file");
			extPotSet = true;
		} else {
			extPotSet = false;
		}
		resetHessianCriterion = MyInput->GetReal("newton",name,"reset_hessian_criterion",0,1e10,1e3);
		if (MyInput->ValueSet("newton",name,"no_hessian_scaling")){
			Message(literal,"The flag 'no_hessian_scaling' is ignored, use "
				"'max_accuracy_for_hessian_scaling' to change default behaviour");
		}
		max_accuracy_for_hessian_scaling = MyInput->GetReal("newton",name,"max_accuracy_for_hessian_scaling",0,1,0.1);
		ignore_newton_direction = MyInput->GetBoolean("newton",name,"ignore_newton_direction",false);
	}
	else {
		name = "no name";
		setiterationlimit(1000);
		setDIIS(0);
		e_info = !silentForced;
		setdeltamax(0.5);
		initGuessSym = firstToLast;
		settolerance(1e-7);
		historyTest = false;
		pseudohessian = true;
		transferHessian = false;
		writeInitialGuess = false;
		smallAlpha = 1e-5;
		maxNumSmallAlpha = 50;
		reverseDirectionRange = 50;
		maxFrReverseDirection = 0.4;
		minAccuracyForHessian = 0;
		numIterationsForHessian = getiterationlimit()+100;
		extPotSet = false;
		resetHessianCriterion = 1e5;
		ignore_newton_direction = false;
		use_vector_iterations = false;
		LM_BFGS = false;
		CG = false;
		SD = false;
		TN = false;
		BRR = false;
		Pure_DIIS = false;

		print_iteration_info = 1;
		print_verbose_info = false;
		print_common_info = false;
		print_exportable_info = false;
		print_improvement_info = false;
		Scheutjens = true;
		reset_limit_gradient_inversion = -1;
		//amount_of_stored_vectors = 12;
		setlinetolerance(0.9);
	}
}
void
SF_Solve::CheckSolution() {
	bool printedMessage;

	errorOccurred = MolQ->CheckMolecules(gettolerance());

	if (charged) {
		Vector q = SegQ->GetCharge();
		Lat->SubtractBoundaries(q);
		double totalCharge = 0;
		for (int z=1; z<=M; z++) {
			totalCharge += q[z];
		}
		if (fabs(totalCharge) > gettolerance()*Lat->GetTotalNumLatticeSites()) {
			errorOccurred = true;
			Message(literal, "No valid solution: system not neutral.!");
			Sysout().Outtext("Total charge in system: ");
			Sysout().Outreal(totalCharge,9,0);
			Sysout().Outimage();
			Sysout().Outtext("Total number of lattice sites in system: ");
			Sysout().Outreal(Lat->GetTotalNumLatticeSites(),9,0);
			Sysout().Outimage();
		}
	}

	if (ReactionQ->IterateAlphas()) {
		if (ReactionQ->gettolerance() < ReactionQ->getaccuracy()) {
			errorOccurred = true;
			Message(literal, "No valid solution: alphabulk values not valid");
		}
	}
	printedMessage = false;
	Vector phiTotal = SegQ->GetPhiTotal();
	for (int z=1; z<=M; z++) {
		if (abs(phiTotal[z]-1.0) > gettolerance()*10 && Lat->WithinBoundaries(z)) {
			errorOccurred = true;
			if (!printedMessage) {
				Message(literal, "No valid solution: sum phi not unity.");
				cout.precision(int(-log10(gettolerance())+1));
				cout << "z " << z << " sumphi "<<phiTotal[z] << endl;
				printedMessage = true;
			}
		}
	}

	if (getaccuracy() != fabs(getaccuracy())) { //true if accuracy == NaN
		WriteOverflowError();
		errorOccurred = true;
	}
	if (getaccuracy() > gettolerance()) {
		errorOccurred = true;
	}
	if (errorOccurred) {cout << "Problem not solved."<< endl; } else {cout << "Problem solved." << endl;};
}
