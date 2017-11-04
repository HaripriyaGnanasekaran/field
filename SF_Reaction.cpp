#include "SF_Reaction.h"

SF_Reaction::SF_Reaction(Text name_,
						 SF_MoleculeList* MolQ_,
						 SF_SegmentList* SegQ_,
						 Input* MyInput_) {
	name = name_;
	SegQ = SegQ_;
	MolQ = MolQ_;
	MyInput = MyInput_;
	Array<Text> param(1,3);
	param[1] = "equation";
	param[2] = "K";
	param[3] = "pK";
	MyInput->CheckParameterNames("reaction",name,param);
	equation = MyInput->GetText("reaction",name,"equation");
	equation.Setpos(1);
	if (!equation.More()) {
		Message(fatal,MyInput,"Syntax error in reaction : "
			+ name + " : equation : " + equation + "\nMissing parameter");
	}
	// remove spaces from composition
	char test;
	Text t = Blanks(equation.Length());
	while (equation.More()) {
		test = equation.Getchar();
		if (!isspace(test)) t.Putchar(test);
	}
	equation = Copy(t.Frontstrip().Strip());
	equation.Setpos(1);
	test = equation.Getchar();
	numLeft = 0;
	// check syntax and count number of different entries on left and right
	while (equation.More() && test != '=') {
		if (isdigit(test)) {
			--equation;
			equation.Getint();
			if (!equation.More()) {
				Message(fatal,MyInput,"Syntax error in reaction : "
					+ name + " : equation : " + equation);
			}
			test = equation.Getchar();
			if (test != '(') {
				Message(fatal,MyInput,"Syntax error in reaction : "
					+ name + " : equation : " + equation);
			}
			while (equation.More() && test != ')') {
				test = equation.Getchar();
			}
			if (test != ')') {
				Message(fatal,MyInput,"Syntax error in reaction : "
					+ name + " : equation : " + equation);
			}
			numLeft++;
		} else {
			Message(fatal,MyInput,"Syntax error in reaction : "
				+ name + " : equation : " + equation);
		}
		if (equation.More()) {
			test = equation.Getchar();
		} else {
			Message(fatal,MyInput,"Syntax error in reaction : "
				+ name + " : equation : " + equation);
		}
	}
	if (test != '=') {
		Message(fatal,MyInput,"Syntax error in reaction : "
			+ name + " : equation : " + equation + "\nMissing '=' sign");
	}
	if (!equation.More()) {
		Message(fatal,MyInput,"Syntax error in reaction : "
			+ name + " : equation : " + equation + "\nPlease define a right hand side");
	}
	test = equation.Getchar();
	numRight = 0;
	while (equation.More()) {
		if (isdigit(test)) {
			--equation;
			equation.Getint();
			if (!equation.More()) {
				Message(fatal,MyInput,"Syntax error in reaction : "
					+ name + " : equation : " + equation);
			}
			test = equation.Getchar();
			if (test != '(') {
				Message(fatal,MyInput,"Syntax error in reaction : "
					+ name + " : equation : " + equation);
			}
			while (equation.More() && test != ')') {
				test = equation.Getchar();
			}
			if (test != ')') {
				Message(fatal,MyInput,"Syntax error in reaction : "
					+ name + " : equation : " + equation);
			}
			numRight++;
		} else {
			Message(fatal,MyInput,"Syntax error in reaction : "
				+ name + " : equation : " + equation);
		}
		if (equation.More()) {
			test = equation.Getchar();
		}
	}
	// fill left and right arrays
	LeftState.Dim(1,numLeft);
	LeftStoch.Dim(1,numLeft);
	LeftComplex.Dim(1,numLeft);
	RightState.Dim(1,numRight);
	RightStoch.Dim(1,numRight);
	RightComplex.Dim(1,numRight);
	numLeft = 0;
	Text stateName;
	equation.Setpos(1);
	test = equation.Getchar();
	while (equation.More() && test != '=') {
		--equation;
		numLeft++;
		LeftStoch[numLeft] = equation.Getint();
		test = equation.Getchar();
		test = equation.Getchar();
		stateName = Blanks(100);
		while (equation.More() && test != ')') {
			stateName.Putchar(test);
			test = equation.Getchar();
		}
		stateName = Copy(stateName.Strip());
		if (!SegQ->StateDefined(stateName)) {
			Message(fatal,MyInput,"Syntax error in reaction : "
				+ name + " : equation : " + equation + "\nCannot use '"
				+ stateName + "' in equation, " + STATE + "/" + SEGMENT
				+ " not defined, or deleted if not used in a molecule.");
		}
		LeftState[numLeft] = SegQ->GetState(stateName);
		test = equation.Getchar();
	}
	test = equation.Getchar();
	numRight = 0;
	while (equation.More()) {
		--equation;
		numRight++;
		RightStoch[numRight] = equation.Getint();
		test = equation.Getchar();
		test = equation.Getchar();
		stateName = Blanks(100);
		while (equation.More() && test != ')') {
			stateName.Putchar(test);
			test = equation.Getchar();
		}
		stateName = Copy(stateName.Strip());
		if (!SegQ->StateDefined(stateName)) {
			Message(fatal,MyInput,"Syntax error in reaction : "
				+ name + " : equation : " + equation + "\nCannot use '"
				+ stateName + "' in equation, " + STATE + "/" + SEGMENT
				+ " not defined, or deleted if not used in a molecule.");
		}
		RightState[numRight] = SegQ->GetState(stateName);
		if (equation.More()) {
			test = equation.Getchar();
		}
	}
	// eliminate double usage of one state at the same side.
	// missing code
	// check for segments with freedom frozen, only allowed at both sides
	SF_State* State;
	int i;
	for (i=1; i<=numLeft; i++) {
		State = LeftState[i];
		if (State->GetFreedom() == frozen) {
			Boolean error = true;
			for (int j=1; j<=numRight; j++) {
				if (SegQ->GetBaseSegment(State) == SegQ->GetBaseSegment(RightState[j])) {
					if (RightStoch[i] != LeftStoch[j]) {
						Message(fatal,MyInput,"Syntax error in reaction : "
						+ name + " : equation : " + equation + "\nCannot use '"
						+ State->GetName()
						+ "' with different stocheometry on both sides of the equation, "
						+ STATE + "/" + SEGMENT + " has freedom set to 'frozen'.");
					}
					error = false;
				}
			}
			if (error) {
				Message(fatal,MyInput,"Syntax error in reaction : "
				+ name + " : equation : " + equation + "\nCannot use '"
				+ State->GetName() + "' on only one side of the equation, "
				+ STATE + "/" + SEGMENT + " has freedom set to 'frozen'.");
			}
		}
	}
	for (i=1; i<=numRight; i++) {
		State = RightState[i];
		if (State->GetFreedom() == frozen) {
			Boolean error = true;
			for (int j=1; j<=numLeft; j++) {
				if (SegQ->GetBaseSegment(State) == SegQ->GetBaseSegment(LeftState[j])) {
					error = false;
				}
			}
			if (error) {
				Message(fatal,MyInput,"Syntax error in reaction : "
				+ name + " : equation : " + equation + "\nCannot use '"
				+ State->GetName() + "' on only one side of the equation, "
				+ STATE + "/" + SEGMENT + " has freedom set to 'frozen'.");
			}
		}
	}
	// check for chainsegments, only allowed at both sides
	SF_Segment* Segment;
	SF_Segment* Segment2;
	SF_Molecule* Molecule;
	SF_MolStructure* Chain;
	for (i=1; i<=numLeft; i++) {
		Segment = SegQ->GetBaseSegment(LeftState[i]);
		for (int j=1; j<=MolQ->GetNumMolecules(); j++) {
			Molecule = MolQ->GetMolecule(j);
			Chain = Molecule->GetMolStructure();
			for (int k=1; k<=Chain->GetNumDiffSegments(); k++) {
				Segment2 = SegQ->GetBaseSegment(Chain->GetDiffSegment(k));
				if (Chain->GetAvLength() > 1 && Segment == Segment2) {
					Boolean error = true;
					for (int l=1; l<=numRight; l++) {
						Segment2 = SegQ->GetBaseSegment(RightState[l]);
						if (Segment == Segment2) {
							if (LeftStoch[i] != RightStoch[l]) {
								Message(fatal,MyInput,"Syntax error in reaction : "
								+ name + " : equation : " + equation + "\nCannot use (states of) '"
								+ SEGMENT + " : " + Segment->GetName()
								+ "' with different stocheometry on both sides of the equation, "
								+ SEGMENT
								+ " belongs to a molecule with more than one segment in the composition.");
							}
							error = false;
						}
					}
					if (error) {
						Message(fatal,MyInput,"Syntax error in reaction : "
						+ name + " : equation : " + equation + "\nCannot use (states of)'"
						+ SEGMENT + " : " + Segment->GetName() + "' on only one side of the equation, "
						+ SEGMENT
						+ " belongs to a molecule with more than one segment in the composition.");
					}
				}
			}
		}
	}
	for (i=1; i<=numRight; i++) {
		Segment = SegQ->GetBaseSegment(RightState[i]);
		for (int j=1; j<=MolQ->GetNumMolecules(); j++) {
			Molecule = MolQ->GetMolecule(j);
			Chain = Molecule->GetMolStructure();
			for (int k=1; k<=Chain->GetNumDiffSegments(); k++) {
				Segment2 = SegQ->GetBaseSegment(Chain->GetDiffSegment(k));
				if (Chain->GetAvLength() > 1 && Segment == Segment2) {
					Boolean error = true;
					for (int l=1; l<=numLeft; l++) {
						Segment2 = SegQ->GetBaseSegment(LeftState[l]);
						if (Segment == Segment2) {
							error = false;
						}
					}
					if (error) {
						Message(fatal,MyInput,"Syntax error in reaction : "
						+ name + " : equation : " + equation + "\nCannot use (states of)'"
						+ SEGMENT + " : " + Segment->GetName() + "' on only one side of the equation, "
						+ SEGMENT
						+ " belongs to a molecule with more than one segment in the composition.");
					}
				}
			}
		}
	}
	// check for electroneutral reaction
	double charge = 0;
	for (i=1; i<=numLeft; i++) {
		State = LeftState[i];
		charge += LeftStoch[i]*State->GetValence();
	}
	for (i=1; i<=numRight; i++) {
		State = RightState[i];
		charge -= RightStoch[i]*State->GetValence();
	}
	if (fabs(charge) > 1e-10) {
		Message(fatal,MyInput,"Syntax error in reaction : "
		+ name + " : equation : " + equation
		+ "\nThe reaction is not electroneutral.");
	}
	// see if internal free energies are given for the states, if not require a K or pK value for the reaction
	internalFreeEnergiesGiven = false;
	for (i=1; i<=numLeft; i++) {
		State = LeftState[i];
		if (State->InternalFreeEnergyGiven()) {
			internalFreeEnergiesGiven = true;
		} else if (internalFreeEnergiesGiven) {
			Message(fatal, "Sorry don't know how to proceed, "
				"please also set internal_energy or internal_degeneration "
				"for state " + State->GetName());
		}
		if (CreatedOrDestroyed(State)) {
			LeftComplex[i] = true;
		} else {
			LeftComplex[i] = false;
		}
	}
	for (i=1; i<=numRight; i++) {
		State = RightState[i];
		if (State->InternalFreeEnergyGiven()) {
			internalFreeEnergiesGiven = true;
		} else if (internalFreeEnergiesGiven) {
			Message(fatal, "Sorry don't know how to proceed, "
				"please also set internal_energy or internal_degeneration "
				"for state " + State->GetName());
		}
		if (CreatedOrDestroyed(State)) {
			RightComplex[i] = true;
		} else {
			RightComplex[i] = false;
		}
	}
	// also catch conflict in the first state
	State = LeftState[1];
	if (State->InternalFreeEnergyGiven()) {
		internalFreeEnergiesGiven = true;
	} else if (internalFreeEnergiesGiven) {
		Message(fatal, "Sorry don't know how to proceed, "
			"please also set internal_energy or internal_degeneration "
			"for state " + State->GetName());
	}
	if (internalFreeEnergiesGiven) {
		if (MyInput->ValueSet("reaction",name,"K") || MyInput->ValueSet("reaction",name,"pK")) {
			Message(fatal, "Sorry don't know how to proceed, "
			"both internal_energy or internal_degeneration for the states "
			"are set together with a reaction constant, this I can't handle");
		}
	} else if (MyInput->ValueSet("reaction",name,"K")) {
		if (MyInput->ValueSet("reaction",name,"pK")) {
			Message(fatal,MyInput,"Cannot set both K and pK for reaction : " + name);
		}
		pK = -log10(MyInput->GetReal("reaction",name,"K",-DBL_MAX,DBL_MAX));
	} else if (MyInput->ValueSet("reaction",name,"pK")) {
		pK = MyInput->GetReal("reaction",name,"pK",-DBL_MAX,DBL_MAX);
	} else {
		Message(fatal,MyInput,"Please set the reaction constant parameter K or pK for reaction : " + name);
	}
	AssignStateRatios();
}
SF_Reaction::~SF_Reaction() {
}
Text
SF_Reaction::GetName() {
	return name;
}
void
SF_Reaction::GetOutput(Output* Out) {
	Out->PutText("reaction",name,"equation",equation);
	Out->PutReal("reaction",name,"pK",pK);
	Out->PutReal("reaction",name,"pK-eff bulk",pKeff);
	Out->PutReal("reaction",name,"residue",ComputeResidue());
}
double
SF_Reaction::ComputeResidue() {
	CalculatePKeff();
	SF_State* State;
/*					Sysout().Outtext("reaction");
					Sysout().Outchar('\t');
					Sysout().Outtext(name);
					Sysout().Outchar('\t');
					Sysout().Outtext(" :");
					Sysout().Outimage();	*/
	double valueLeft = 0;
	int i;
	int currNumStateLeft=0;
	for (i=1; i<=numLeft; i++) {
		State = LeftState[i];
		if (LeftComplex[i]) {
			Text molName = (State->GetMolState(1))->GetMolName();
			SF_Molecule* Mol = MolQ->GetMolecule(molName);
			double ln_f = MolQ->GetChemicalPotential(Mol) - log(Mol->GetPhiBulk());
			valueLeft += LeftStoch[i]*(ln_f + log(State->GetPhiBulk()))/log(10.);
//					Sysout().Outtext("value left");
//					Sysout().Outchar('\t');
//					Sysout().Outreal(ln_f,9,0);
//					Sysout().Outimage();
//					Sysout().Outreal(LeftStoch[i],9,0);
//					Sysout().Outimage();
//					Sysout().Outreal((ln_f + log(State->GetPhiBulk())),9,0);
//					Sysout().Outimage();
//					Sysout().Outreal(ln_f,9,0);
//					Sysout().Outimage();
		} else {
			for (int j=1; j<=LeftStoch[i]; j++) {
				valueLeft += log10((Ratios[++currNumStateLeft])->GetAlphaBulkRatio());
//					Sysout().Outtext("value left");
//					Sysout().Outchar('\t');
//					Sysout().Outint(j,0);
//					Sysout().Outchar('\t');
//					Sysout().Outreal(log10((Ratios[currNumStateLeft])->GetAlphaBulkRatio()),9,0);
//					Sysout().Outimage();
//					Sysout().Outtext(((Ratios[currNumStateLeft])->GetState1())->GetName());
//					Sysout().Outchar('\t');
//					Sysout().Outreal(((Ratios[currNumStateLeft])->GetState1())->GetAlphaBulk(),9,0);
//					Sysout().Outimage();
//					Sysout().Outtext(((Ratios[currNumStateLeft])->GetState2())->GetName());
//					Sysout().Outchar('\t');
//					Sysout().Outreal(((Ratios[currNumStateLeft])->GetState2())->GetAlphaBulk(),9,0);
//					Sysout().Outimage();
			}
		}
	}
	double valueRight = 0;
	for (i=1; i<=numRight; i++) {
		State = RightState[i];
		if (RightComplex[i]) {
			Text molName = (State->GetMolState(1))->GetMolName();
			valueRight += RightStoch[i]*MolQ->GetChemicalPotential(MolQ->GetMolecule(molName))/log(10.0);
		}
	}
//					Sysout().Outtext("endvalues");
//					Sysout().Outchar('\t');
//					Sysout().Outreal(valueLeft,9,0);
//					Sysout().Outchar('\t');
//					Sysout().Outreal(valueRight,9,0);
//					Sysout().Outchar('\t');
//					Sysout().Outreal((valueLeft-valueRight)/pK,9,0);
//					Sysout().Outchar('\t');
//					Sysout().Outreal(-1+(valueLeft-valueRight)/pK,9,0);
//					Sysout().Outimage();
	return -1+(valueLeft-valueRight)/pKeff;
}
int
SF_Reaction::GetNumDiffStates() {
	return numLeft + numRight;
}
SF_State*
SF_Reaction::GetDiffState(int number) {
	if (number > numLeft) return RightState[number-numLeft];
	else return LeftState[number];
}
void
SF_Reaction::AssignStateRatios() {
	int i;
	SF_State* State1;
	SF_State* State2;
	numRatios = 0;
	for (i=1; i<=numLeft; i++) {
		State1 = LeftState[i];
		if (!CreatedOrDestroyed(State1)) {
			numRatios+=LeftStoch[i];
		}
	}
	Ratios.Dim(1,numRatios);
	int currNumLeftState = 0;
	for (i=1; i<=numLeft; i++) {
		State1 = LeftState[i];
		if (!CreatedOrDestroyed(State1)) {
			for (int j=1; j<=LeftStoch[i]; j++) {
				currNumLeftState++;
				int currNumRightState = 0;
				for (int k=1; k<=numRight; k++) {
					State2 = RightState[k];
					if (!CreatedOrDestroyed(State2)) {
						for (int m=1; m<=RightStoch[k]; m++) {
							currNumRightState++;
							if (currNumRightState == currNumLeftState) {
								Ratios[currNumLeftState] = new SF_StateRatio(State1,State2);
							}
						}
					}
				}
			}
		}
	}
}
int
SF_Reaction::GetNumStateRatios() {
	return numRatios;
}
SF_StateRatio*
SF_Reaction::GetStateRatio(int number) {
	return Ratios[number];
}
void
SF_Reaction::SwapDirection() {
	Array<int> StochTemp = LeftStoch;
	LeftStoch = RightStoch;
	RightStoch = StochTemp;
	Array<SF_State*> StateTemp;
	LeftState = RightState;
	RightState = StateTemp;
	pK = -pK;
	// swap equation code missing
}
void
SF_Reaction::ReplaceRatio(SF_StateRatio* Ratio1, SF_StateRatio* Ratio2) {
	for (int i=1; i<=numRatios; i++) {
		if (Ratios[i] == Ratio1) {
			Ratios[i] = Ratio2;
		}
	}
}
Boolean
SF_Reaction::IterateAlphas() {
	Boolean value = false;
	SF_State* State;
	int i,j,k;
	SF_Segment* Segment=NULL;
	SF_Molecule* Mol;
	SF_MolStructure* Chain;
	for (i=1; i<=numLeft; i++) {
		State = LeftState[i];
		if (CreatedOrDestroyed(State)) value = true;
		Segment = SegQ->GetBaseSegment(State);
		for (j=1; j<=MolQ->GetNumMolecules(); j++) {
			Mol = MolQ->GetMolecule(j);
			Chain = Mol->GetMolStructure();
			for (k=1; k<=Chain->GetNumDiffSegments(); k++) {
				if (Segment == SegQ->GetBaseSegment(Chain->GetDiffSegment(k))) {
					if (Mol->GetFreedom() != fixedPhiBulk
					&&  Mol->GetFreedom() != secondGeneration
					&&  Mol->GetFreedom() != thirdGeneration) {
						value = true;
					}
				}
			}
		}
	}
	for (i=1; i<=numRight; i++) {
		State = RightState[i];
		if (CreatedOrDestroyed(State)) value = true;
		for (j=1; j<=MolQ->GetNumMolecules(); j++) {
			Mol = MolQ->GetMolecule(j);
			Chain = Mol->GetMolStructure();
			for (k=1; k<=Chain->GetNumDiffSegments(); k++) {
				if (Segment == SegQ->GetBaseSegment(Chain->GetDiffSegment(k))) {
					if (Mol->GetFreedom() != fixedPhiBulk
					&&  Mol->GetFreedom() != secondGeneration
					&&  Mol->GetFreedom() != thirdGeneration) {
						value = true;
					}
				}
			}
		}
	}
	return value;
}
void
SF_Reaction::ComputeInternalFreeEnergies() {
 	if (internalFreeEnergiesGiven) {
 		Message(debug,"wtf, SF_Reaction::ComputeInternalFreeEnergies()");
 	}
 	int count = 0;
	double G = 0;
	SF_State* State=NULL;
	int stoch = 0;
	int i;
	for (i=1; i<=numRight; i++) {
		if (!(RightState[i])->InternalFreeEnergySet()) {
			count++;
			State = RightState[i];
			stoch = RightStoch[i];
		} else {
			G += RightStoch[i]*(RightState[i])->GetInternalFreeEnergy();
		}
	}
	for (i=1; i<=numLeft; i++) {
		if (!(LeftState[i])->InternalFreeEnergySet()) {
			count++;
			State = LeftState[i];
			stoch = -LeftStoch[i];
		} else {
			G -= LeftStoch[i]*(LeftState[i])->GetInternalFreeEnergy();
		}
	}
	if (count != 1) {
		return;
	}
	State->SetInternalFreeEnergy((log(10.0)*pK-G)/stoch);
}
void
SF_Reaction::UpdateReactionConstant(double roomTemp,double temp) {
 	if (!internalFreeEnergiesGiven) {
		return;
 	}
	SF_State* State;
	int i;
	double valueLeft = 0;
	for (i=1; i<=numLeft; i++) {
		State = LeftState[i];
		State->UpdateInternalFreeEnergy(roomTemp,temp);
		valueLeft += LeftStoch[i]*State->GetInternalFreeEnergy();
	}
	double valueRight = 0;
	for (i=1; i<=numRight; i++) {
		State = RightState[i];
		State->UpdateInternalFreeEnergy(roomTemp,temp);
		valueRight += RightStoch[i]*State->GetInternalFreeEnergy();
	}
	pK = (valueRight - valueLeft)/log(10.0);
}
Boolean
SF_Reaction::CreatedOrDestroyed(SF_State* State) {
	int count = 0;
	int i;
	for (i=1; i<=numRight; i++) {
		if (SegQ->GetBaseSegment(State) == SegQ->GetBaseSegment(RightState[i])) {
			count += RightStoch[i];
		}
	}
	for (i=1; i<=numLeft; i++) {
		if (SegQ->GetBaseSegment(State) == SegQ->GetBaseSegment(LeftState[i])) {
			count -= LeftStoch[i];
		}
	}
	if (count != 0) return true;
	else return false;
}
void
SF_Reaction::CalculatePKeff() {
	SF_State* State;
	double valueLeft = 0;
	int i;
	for (i=1; i<=numLeft; i++) {
		State = LeftState[i];
		valueLeft += LeftStoch[i]*SegQ->ChemIntBulk(State);
	}
	double valueRight = 0;
	for (i=1; i<=numRight; i++) {
		State = RightState[i];
		valueRight += RightStoch[i]*SegQ->ChemIntBulk(State);
	}
	pKeff = pK + (valueRight - valueLeft)/log(10.0);
}

