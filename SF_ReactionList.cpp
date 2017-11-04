#include "SF_ReactionList.h"

SF_ReactionList::SF_ReactionList(SF_MoleculeList* MolQ_,
								 SF_SegmentList* SegQ_,
								 Input* MyInput_) {
	MolQ = MolQ_;
	SegQ = SegQ_;
	MyInput = MyInput_;
	phiBulksFixed = false;
	numReactions = MyInput->GetNumNames("reaction");
	ReactionQ.Dim(1,numReactions);
	Array<Text> names = MyInput->GetNames("reaction");
	int i;
	for (i=1; i<=numReactions; i++) {
		ReactionQ[i] = new SF_Reaction(names[i],MolQ,SegQ,MyInput);
	}
	pseudohessian = true;
	forcePseudohessian = false;
	settolerance(1e-20);
	oldAccuracy = 1;
//	e_info = true;
//	x_info = true;
//	g_info = true;
//	d_info = true;
//	g_info = true;
	setdeltamax(0.1);
	setdeltamin(1e-7);
	setiterationlimit(10000);
	smallAlpha = 1e-3;
	maxNumSmallAlpha = 200;
	SetNeutralizerNeeded();
	internalFreeEnergiesGiven = false;
	for (i=1; i<=numDiffStates; i++) {
		if ( (StateQ[i])->InternalFreeEnergyGiven() ) {
			internalFreeEnergiesGiven = true;
		} else if (internalFreeEnergiesGiven) {
			Message(fatal, "Sorry don't know how to proceed, "
				"please also set internal_energy or internal_degeneration "
				"for state " + (StateQ[i])->GetName());
		}
	}
	x.Dim(1,numItVar);
	residues.Dim(1,numItVar);
	residualMap.Dim(1,numItVar);
	MapItVar();
}
SF_ReactionList::~SF_ReactionList() {
	int i;
	for (i=1; i<=numReactions; i++) {
		delete ReactionQ[i];
	}
	for (i=1; i<=numSegmentRatios; i++) {
		delete SegmentRatioQ[i];
	}
	for (i=1; i<=numStateRatios; i++) {
		delete StateRatioQ[i];
	}
}
void
SF_ReactionList::GetOutput(Output* Out) {
	for (int i=1; i<=numReactions; i++) {
		(ReactionQ[i])->GetOutput(Out);
	}
}
int
SF_ReactionList::GetNumReactions() {
	return numReactions;
}
SF_Reaction*
SF_ReactionList::GetReaction(int number) {
	return ReactionQ[number];
}
Boolean
SF_ReactionList::NeutralizerNeeded() {
	return neutralizerNeeded;
}
void
SF_ReactionList::SetNeutralizerNeeded() {
	SF_Reaction* Reaction;
	int i;
	numDiffStates = 0;
	for (i=1; i<=numReactions; i++) {
		Reaction = ReactionQ[i];
		numDiffStates += Reaction->GetNumDiffStates();
	}
	StateQ.Dim(1,numDiffStates);
	int count = 0;
	for (i=1; i<=numReactions; i++) {
		Reaction = ReactionQ[i];
		int limit = Reaction->GetNumDiffStates();
		for (int j=1; j<=limit; j++) {
			StateQ[++count] = Reaction->GetDiffState(j);
		}
	}
	for (i=1; i<=numDiffStates; i++) {
		for (int j=i+1; j<=numDiffStates; j++) {
			if (StateQ[i] == StateQ[j]) {
				Array<SF_State*> StateQTemp(1,numDiffStates-1);
				for (int k=1; k<=numDiffStates; k++) {
					if (k<j) StateQTemp[k] = StateQ[k];
					if (k>j) StateQTemp[k-1] = StateQ[k];
				}
				numDiffStates--;
				j--;
				StateQ = StateQTemp;
			}
		}
	}
	int numStatesWithFixedPhiBulk = 0;
	int numStatesWithFixedAlphaBulk = 0;
	for (i=1; i<=numDiffStates; i++) {
		if ((StateQ[i])->PhiBulkFixed()) {
			numStatesWithFixedPhiBulk++;
		}
		if ((StateQ[i])->AlphaBulkFixed()) {
			numStatesWithFixedAlphaBulk++;
//			Sysout().Outtext((StateQ[i])->GetName());
//			Sysout().Outimage();
		}
	}
	SF_Segment* Segment;
	int numStatesInSegments = 0;
	int numSegmentsWithDifferentStates = 0;
	for (i=1; i<=SegQ->GetNumSegments(); i++) {
		Segment = SegQ->GetSegment(i);
		if (Segment->GetNumStates() > 1) {
			numSegmentsWithDifferentStates++;
			numStatesInSegments += Segment->GetNumStates();
			for (int j=1; j<=Segment->GetNumStates(); j++) {
				SF_State* State = Segment->GetState(j);
				Boolean found = false;
				for (int k=1; k<=numDiffStates; k++) {
					if (StateQ[k] == State) {
						found = true;
					}
				}
				if (!found) {Message(fatal,"A " + STATE + " named '" + State->GetName()+ "' was defined but never used in a reaction");}
			}

		}
	}
	int numEquations = numReactions + numStatesWithFixedPhiBulk + numStatesWithFixedAlphaBulk + 2;
	int numUnknowns = numStatesInSegments - numSegmentsWithDifferentStates + 1;
	if (!SegQ->Charged()) {
		numEquations--;
	}
	if (numStatesWithFixedPhiBulk + numStatesWithFixedAlphaBulk > 1) {
		Message(fatal,"At most one component can have a fixed phibulk or alphabulk");
	}
	if (numStatesWithFixedPhiBulk + numStatesWithFixedAlphaBulk > 0 && !SegQ->Charged()) {
		Message(fatal,"A fixed phibulk or alphabulk and no charge? Cannot proceed, a neutralizer is needed according to the phase rule");
	}
	numItVar = numEquations;
	if (numEquations > numUnknowns) {
		if (SegQ->Charged()) {
			neutralizerNeeded = true;
		} else {
			neutralizerNeeded = false;
		}
		if (numEquations > numUnknowns + 1) {
			Message(fatal,"Unable to solve this set of equations, more equations than unknowns");
		}
	} else {
		neutralizerNeeded = false;
		if (numEquations < numUnknowns) {
			Message(fatal,"Unable to solve this set of equations, more unknowns than equations");
		}
		if (MolQ->GetNeutralizer() != NULL) {
			Message(fatal,"neutralizer not needed or specify the bulk volumefraction of one of the states");
		}
	}
}
void
SF_ReactionList::ComputeAlphaBulk() {
	restartCount = 0;
	SetItVar();
	MapResiduals();
	int j;
	haltonFPE = false;
	for (j=1; j<=3; j++) {
		iterate(&x[1],numItVar);
	}
	if (getaccuracy() != fabs(getaccuracy()) && getminimum() < 1e-1) {
		forcePseudohessian = true;
		iterate(&x[1],numItVar);
	}
	forcePseudohessian = false;
	while (!(getaccuracy() >= 0) || getaccuracy() > gettolerance()) {
		restartCount++;
		if (restartCount > 20) {
			ErrorOutput();
		}
		if (getiterations() < getiterationlimit()) {
			for (int j=1; j<=numItVar; j++) {
				x[j] = 0;
			}
			CopyItVar();
		}
		Sysout().Outtext("Don't worry, I'll try again");
		Sysout().Outimage();
		restart = false;
		MapResiduals();
		for (j=1; j<=3; j++) {
			iterate(&x[1],numItVar);
		}
		if (!(getaccuracy() >= 0) || getaccuracy() > gettolerance()) {
			for (int i=1; i<=SegQ->GetNumSegments(); i++) {
				SF_Segment* Seg = SegQ->GetSegment(i);
				int numStates = Seg->GetNumStates();
				if (numStates > 1) {
					for (int z=1; z<=numStates; z++) {
						double dummy = 1.0/numStates + 1e-9;
						if (z == 2) dummy -= 5e-9;
						(Seg->GetState(z))->SetAlphaBulk(dummy);
					}
				}
			}
			Sysout().Outtext("Don't worry, I'll try again.");
			Sysout().Outimage();
			SetItVar();
			MapResiduals();
			restart = false;
			for (j=1; j<=3; j++) {
				iterate(&x[1],numItVar);
			}
			if (getaccuracy() != fabs(getaccuracy()) && getminimum() < 1e-1) {
				forcePseudohessian = true;
				iterate(&x[1],numItVar);
			}
		}
	}
	if (!internalFreeEnergiesGiven) {
		ComputeInternalFreeEnergies();
	}
}
void
SF_ReactionList::InnerIterate() {
	iterate(&x[1],numItVar);
	if (getaccuracy() != fabs(getaccuracy()) && getminimum() < 1e-1) {
		forcePseudohessian = true;
		iterate(&x[1],numItVar);
	}
	forcePseudohessian = false;
	if (getaccuracy() > gettolerance()) {
		for (int i=1; i<=SegQ->GetNumStates(); i++) {
			Sysout().Outint(i,0);
			Sysout().Outchar('\t');
			Sysout().Outtext(SegQ->GetState(i)->GetName());
			Sysout().Outchar('\t');
			Sysout().Outreal(SegQ->GetState(i)->GetPhiBulk(),9,0);
			Sysout().Outchar('\t');
			Sysout().Outreal(SegQ->GetState(i)->GetAlphaBulk(),9,0);
			Sysout().Outimage();
		}
	}
}
void
SF_ReactionList::IterateWithFixedPhiBulks() {
	int i,j;
	int panicCount = 0;
	Boolean panic = false;
	Boolean done = false;
	phiBulksFixed = true;
	SetItVar();
	numItVar--;
	if (neutralizerNeeded) {
//		int count = 0;
//		for (i=1; i<=numStateRatios; i++) {
//			SF_State* State = (StateRatioQ[i])->GetState1();
//			if (State->PhiBulkFixed() && !phiBulksFixed) {
//				count++;
//			}
//			if (State->AlphaBulkFixed() && !phiBulksFixed) {
//				count++;
//			}
//		}
//		Sysout().Outtext("count : ");
//		Sysout().Outint(count,0);
//		Sysout().Outimage();
//		if (count != 1) {
//			Message(warning,"Thermodynamics wrong, fix up SF_ReactionList::IterateWithFixedPhiBulks()");
//		}
		for (i=1; i<=numItVar; i++) {
			if (residualMap[i] == 1) {
				NeutralizingRatio = StateRatioQ[i];
//				Sysout().Outtext("neutralizing ratio : ");
//				Sysout().Outint(i,0);
//				Sysout().Outimage();
			}
		}
		numItVar--;
	}
	Matrix Response(0,numItVar,0,numItVar);
	FillResponseMatrix(Response);
//	cout << "iterate2" << endl;
//	for (i=1; i<=numItVar; i++) {
//		for (j=1; j<=numItVar; j++) {
//			Sysout().Outint(i,0);
//			Sysout().Outchar('\t');
//			Sysout().Outint(j,0);
//			Sysout().Outchar('\t');
//			Sysout().Outreal(Response[i][j],9,0);
//			Sysout().Outimage();
//		}
//	}
//	for (i=1; i<=numStateRatios; i++) {
//		Sysout().Outtext((StateRatioQ[i])->GetState1()->GetName());
//		Sysout().Outchar('\t');
//		Sysout().Outtext((StateRatioQ[i])->GetState2()->GetName());
//		Sysout().Outimage();
//	}
//	for (i=1; i<=MolQ->GetNumMolecules(); i++) {
//		Sysout().Outint(i,0);
//		Sysout().Outchar('\t');
//		Sysout().Outtext(MolQ->GetMolecule(i)->GetName());
//		Sysout().Outchar('\t');
//		Sysout().Outreal(MolQ->GetMolecule(i)->GetPhiBulk(),9,0);
//		Sysout().Outimage();
//	}
//	for (i=1; i<=SegQ->GetNumStates(); i++) {
//		Sysout().Outint(i,0);
//		Sysout().Outchar('\t');
//		Sysout().Outtext(SegQ->GetState(i)->GetName());
//		Sysout().Outchar('\t');
//		Sysout().Outreal(SegQ->GetState(i)->GetPhiBulk(),9,0);
//		Sysout().Outchar('\t');
//		Sysout().Outreal(SegQ->GetState(i)->GetAlphaBulk(),9,0);
//		Sysout().Outimage();
//	}
	Vector dontTouch(1,numItVar);
	for (i=1; i<=numItVar; i++) {
		Boolean found = false;
		Boolean foundAgain = false;
		for (j=1; j<=numItVar; j++) {
			if (Response[i][j] > 0) {
				if (found) {
					foundAgain = true;
				}
				found = true;
				residualMap[i] = j;
			}
		}
		if (!foundAgain) {
			dontTouch[i] = 1;
		}
	}
	while (!done) {
		for (i=1; i<=numItVar; i++) {
			if (dontTouch[i] < 1) {
				residualMap[i] = 0;
			}
		}
		for (i=1; i<=numItVar; i++) {
			if (dontTouch[i] > 0) continue;
			Boolean abort = false;
			Boolean ok = false;
			while (!abort) {
				abort = false;
				while (!ok) {
					residualMap[i] = int(rand()*double(numItVar)/(RAND_MAX + 1.0))+1;
					ok = true;
					for (j=1; j<=numItVar; j++) {
						if (residualMap[j] == residualMap[i] && (dontTouch[j] > 0 || j<=i-1)) {
							ok = false;
						}
					}
				}
				if (!neutralizerNeeded) SetNeutralizingRatio();
				if (Response[i][residualMap[i]] <= 0 && !panic) {
					abort = true;
				} else if (Response[i][residualMap[i]] > 0) {
				} else if (!panic) { // GetResponse(i,residualMap[i]) == NaN
					abort = true;
				}
				if (abort) ok = false;
                if (ok) abort = true;
			}
			if (abort && !ok) {
				i=0;
				panicCount++;
				if (panicCount > 1000) panic = true;
			}
		}
		done = true;
		for (i=1; i<=numItVar; i++) {
			if (Response[i][residualMap[i]] <= 0 && !panic) done = false;
			else if (Response[i][residualMap[i]] > 0) ;
			else if (!panic) done = false; // GetResponse(i,residualMap[i]) == NaN
		}
	}
//	for (i=1; i<=numItVar; i++) {
//		Sysout().Outint(residualMap[i],0);
//		Sysout().Outimage();
//	}
	iterate(&x[1],numItVar);
	if (neutralizerNeeded) {
		numItVar++;
	}
	numItVar++;
	phiBulksFixed = false;
}
void
SF_ReactionList::residuals(double *const ff, double *const) {
//	fedisableexcept(FE_OVERFLOW|FE_INVALID);
	Vector f(ff,1,numItVar);
	CopyItVar();
	ComputeResidues();
	for (int j=1; j<=numItVar; j++) {
		f[j] = GetResidue(j);
	}
//	feenableexcept(FE_OVERFLOW|FE_INVALID);
}
void
SF_ReactionList::inneriteration(double *const gg, double *const, double) {
	Vector g(gg,1,numItVar);
	if (getiterations() == 0) {
		oldAccuracy = 1;
	}
	if (getiterations() == 3) {
		settolerance(1e-20); // accuracy should be at its best: DBL_EPSILON
	}
	if (getalpha() < smallAlpha) smallAlphaCount++;
	else smallAlphaCount = 0;
	if (smallAlphaCount == maxNumSmallAlpha) {
		resethessian();
		if (e_info) cout << "hessian reset" << endl;
//		residuals(&g[1],&x[1]);
	}
	if (getaccuracy() < gettolerance()) {
		pseudohessian = true;
	} else if (getaccuracy() < 0.1 && getminimum() < 0.1) {
		pseudohessian = false;
	} else {
		pseudohessian = true;
	}
	if (forcePseudohessian) {
		pseudohessian = true;
	}
	if (getaccuracy() < 1e-13) { // minimum accuracy
		if (getaccuracy() >= oldAccuracy) {
			// no progress
			settolerance(gettolerance()*10);
		}
	}
	oldAccuracy = getaccuracy();
}
double
SF_ReactionList::GetResidue(int number) {
	return residues[residualMap[number]];
}
double
SF_ReactionList::ComputeAlphaResidue(SF_State* ThisState) {
	SF_Segment* Segment = SegQ->GetBaseSegment(ThisState);
	double value=0;
	double sumAlpha=0;
	int i;
	SF_State* State;
	for (i=1; i<=Segment->GetNumStates(); i++) {
		State = Segment->GetState(i);

		sumAlpha += State->GetAlphaBulk();
	}
	for (i=1; i<=Segment->GetNumStates(); i++) {
		State = Segment->GetState(i);
		if (State == ThisState) {
			value += State->GetAlphaBulk();
		} else {
			value += State->GetAlphaBulk()/sumAlpha;
		}
	}
	return -1 + value;
}
void
SF_ReactionList::ComputeComplexPartOfChemPots() {
	Matrix derivatives(1,SegQ->GetNumStates(),1,MolQ->GetNumMolecules());
	Message(literal,"Please wait, determining chemical potentials");
	Vector finished(1,SegQ->GetNumStates());
	Matrix errors(1,SegQ->GetNumStates(),1,MolQ->GetNumMolecules());
	int i;
	for (i=1; i<=MolQ->GetNumMolecules(); i++) {
		for (int n=1; n<=SegQ->GetNumStates(); n++) {
			if (SegQ->GetBaseSegment(SegQ->GetState(n))->GetNumStates() > 1) {
				errors[n][i] = 1e30;
			} else {
				errors[n][i] = 0;
			}
		}
	}
	for (i=1; i<=MolQ->GetNumMolecules(); i++) {
		double hh = 0.1;
		cout.precision(5);
		double relativeErr = 1e30;
		double oldRelativeErr = 1e30;
		Vector ans(1,SegQ->GetNumStates());
		Vector ansCandidate(1,SegQ->GetNumStates());
		Vector errCandidate(1,SegQ->GetNumStates());
		Vector err(1,SegQ->GetNumStates());
		do { // ComputeDerivatives is very sensitive to the initial step size hh, try many step sizes
			oldRelativeErr = relativeErr;
			relativeErr = 0;
			cout << "." << flush;
			ComputeDerivatives(ans,err,i,hh);
			int n;
			for (n=1; n<=SegQ->GetNumStates(); n++) {
				if ( ( (errors[n][i] > err[n]/fabs(ans[n])) || (errCandidate[n] > err[n]/fabs(ans[n])) ) && ans[n] != 0 && err[n] > 0) { // save best guess so far, err == 0 is not ok
					ansCandidate[n] = ans[n];
					errCandidate[n] = err[n]/fabs(ans[n]);
				}
				if (ans[n] == 0 && err[n] == 0 && hh == 0.1) { // here err == 0 is ok, even a large step gives ans == 0
					derivatives[n][i] = 0;
					errors[n][i] = 0;
				}
			}
			n = 0;
			// We need to capture the pathetic case: derivatives of states belonging to the same segment type should add up to zero.
			for (int j=1; j<=SegQ->GetNumSegments(); j++) {
				SF_Segment *Seg = SegQ->GetSegment(j);
				double total = 0;
				double max = 0;
				for (int k=1; k<=Seg->GetNumStates(); k++) {
					n++;
					total += ansCandidate[n];
					if (max < fabs(ansCandidate[n])) {
						max = fabs(ansCandidate[n]);
					}
				}
//				cout << Seg->GetName().MainC() << '\t' << total  << '\t' << max << endl;
				if (max == 0) ;
				else if (fabs(total/max) < 1e-8) {
					n -= Seg->GetNumStates();
					for (int k=1; k<=Seg->GetNumStates(); k++) {
						n++;
						derivatives[n][i] = ansCandidate[n];
						errors[n][i] = errCandidate[n];
					}
				} else if (fabs(total/max) < 1e-4 && errors[n-1][i] == 1e30 && max < 1e-5) { // what a disaster
					n -= Seg->GetNumStates();
					for (int k=1; k<=Seg->GetNumStates(); k++) {
						n++;
						derivatives[n][i] = ansCandidate[n];
						errors[n][i] = errCandidate[n];
					}
				}
			}
			for (n=1; n<=SegQ->GetNumStates(); n++) {
				if (errors[n][i] > relativeErr) {
					relativeErr = errors[n][i];
				}
//				cout << MolQ->GetMolecule(i)->GetName().MainC() << '\t' << SegQ->GetState(n)->GetName().MainC() << '\t' << ans[n]  << '\t' << err[n] << '\t' << derivatives[n][i] << '\t' << errors[n][i] << endl;
			}
			hh /= 10;
		} while ( ( (relativeErr > 1e-10) // acceptable accuracy not obtained
			|| oldRelativeErr > relativeErr ) // try to obtain higher accuracy
			&& (hh > 1e-30) // bail out if step gets too small
			);
	}
	cout << endl;
	double sumBulk = 0;
	double sumSystem = 0;
//	for (i=1; i<=MolQ->GetNumMolecules(); i++) {
//		SF_Molecule* Mol = MolQ->GetMolecule(i);
//		for (int n=1; n<=SegQ->GetNumStates(); n++) {
//			SF_State* State = SegQ->GetState(n);
//			cout << State->GetName().MainC() << '\t' << State->GetAlphaBulk() << '\t' << Mol->GetName().MainC() << '\t' << derivatives[n][i] << endl;
//		}
//	}
	for (i=1; i<=MolQ->GetNumMolecules(); i++) {
		SF_Molecule* Mol1 = MolQ->GetMolecule(i);
		double value = 0;
		double N1 = Mol1->GetChainLength();
		int stateNumber = 0;
		for (int j=1; j<=MolQ->GetNumMolecules(); j++) {
			SF_Molecule* Mol2 = MolQ->GetMolecule(j);
			double phiBulk2 = Mol2->GetPhiBulk();
			for (int k=1; k<=Mol2->GetMolStructure()->GetNumDiffSegments(); k++) {
				SF_MolSegment* MolSeg = (Mol2->GetMolStructure())->GetDiffSegment(k);
				for (int m=1; m<=MolSeg->GetNumStates(); m++) {
					SF_MolState* MolState = MolSeg->GetState(m);
					SF_State* BaseState = SegQ->GetState(MolState->GetName());
					int n;
					for (n=1; n<=SegQ->GetNumStates(); n++) {
						if (SegQ->GetState(n) == BaseState) {
							stateNumber = n;
						}
					}
					double phiRefK = MolSeg->GetPhiRef();
					double alphaBulkM = MolState->GetAlphaBulk();
					double dummy = 0;
					dummy += log(alphaBulkM) + 1;
					dummy += MolState->GetInternalFreeEnergy();
					dummy -= MolSeg->GetState(1)->GetInternalFreeEnergy();
					dummy += SegQ->ChemIntBulk(MolState);
//					if (*Copy(MolState->GetMolName()) == *Copy(Mol1->GetName()) ) {
// 		 				 dummy -= 2*SegQ->ChemIntRef(MolState);
// 					}
					value += phiBulk2*phiRefK*N1*derivatives[stateNumber][i]*dummy;
//					cout << derivatives[stateNumber][i]*dummy << endl;
				}
			}
		}
		sumBulk += value*Mol1->GetPhiBulk()/Mol1->GetChainLength();
		sumSystem += value*Mol1->ComputeTheta()/Mol1->GetChainLength();
		Sysout().Outtext(Mol1->GetName()+ " : ");
		Sysout().Outreal(value,9,0);
		Sysout().Outimage();
	}
	cout << "sum bulk: " << sumBulk << endl;
	cout << "sum sys  : " << sumSystem << endl;
}
void
SF_ReactionList::UpdateReactionConstants(double roomTemp, double temp) {
	if (!internalFreeEnergiesGiven) return;
	for (int i=1; i<=numReactions; i++) {
		(ReactionQ[i])->UpdateReactionConstant(roomTemp,temp);
	}
}
void
SF_ReactionList::ComputeResidues() {
	int i,j=1;
	SF_State* State;
	int numSegments = SegQ->GetNumSegments();
	for (i=1; i<=numSegments; i++) {
		if ((SegQ->GetSegment(i))->GetNumStates() > 1) {
			for (int k=1; k<=(SegQ->GetSegment(i))->GetNumStates(); k++) {
				State = (SegQ->GetSegment(i))->GetState(k);
				if (State->PhiBulkFixed() && !phiBulksFixed) {
//					cout << "phibulk: " << State->GetPhiBulk() << endl;
					residues[j++] = -1+State->GetPhiBulk() / State->GetPhiBulkFixed();
//					cout << "  " << residues[j-1] << endl;
				}
				if (State->AlphaBulkFixed() && !phiBulksFixed) {
					residues[j++] = -1+State->GetAlphaBulk() / State->GetAlphaBulkFixed();
				}
			}
		}
	}
	for (i=1; i<=numReactions; i++) {
//	cout << "ComputeResidues" << j << endl;
		residues[j++] = (ReactionQ[i])->ComputeResidue();
	}
	double bulkCharge = 0;
	int numStates = SegQ->GetNumStates();

	for (i=1; i<=numStates; i++) {

		State = SegQ->GetState(i);
		bulkCharge += State->GetValence()*State->GetPhiBulk();
	}
	double compCharge = 1;
	if (SegQ->Charged()) {
		if (neutralizerNeeded && !phiBulksFixed) {
			compCharge = (MolQ->GetNeutralizer())->GetBulkCharge();
			residues[j++] = bulkCharge/compCharge;
		} else {
			State = NeutralizingRatio->GetState1();
//	Sysout().Outtext(State->GetName());
//	Sysout().Outimage();
			compCharge = State->GetValence();
//	cout << compCharge << endl;
			State = NeutralizingRatio->GetState2();
//	Sysout().Outtext(State->GetName());
//	Sysout().Outimage();
			if (compCharge != 0 && State->GetValence() != 0) {
				compCharge *= -State->GetValence();
			} else if (compCharge == 0) {
				if (State->GetValence() != 0) {
					compCharge = -State->GetValence();
				} else {
					Message(debug,"no charge in NeutralizingRatio");
				}
			}
//	cout << "ComputeResidues" << j << endl;
			residues[j++] = bulkCharge/compCharge;
//	cout << "ComputeResidues" << j << endl;
		}
	}
	if (!phiBulksFixed) {
		residues[j] = 0;
		for (i=1; i<=numStates; i++) {
			State = SegQ->GetState(i);
			residues[j] += State->GetPhiBulk();
		}
		residues[j] = 1 - 1/residues[j];
	}
}
double
SF_ReactionList::GetItVar(int number) {
	int i;
	int count=0;
	for (i=1; i<=numStateRatios; i++) {
		(StateRatioQ[i])->ComputeAlphaBulkRatio();
	}
	for (i=1; i<=numStateRatios; i++) {
		count++;
		if (count == number) {
			return log10((StateRatioQ[i])->GetAlphaBulkRatio());
		}
	}
	if (neutralizerNeeded && SegQ->Charged()) {
		count++;
		if (count == number) {
			double dummy = (MolQ->GetNeutralizer())->GetPhiBulk();
			return log(dummy/(1-dummy));
		}
	}
	count++;
	if (count == number) {
		double dummy = (MolQ->GetSolvent())->GetPhiBulk();
		return log(dummy/(1-dummy));
	}
	Message(fatal,"program error in SF_ReactionList::GetItVar(int number)");
	return 0; // never get here
}
void
SF_ReactionList::SetItVar() {
	int i,j=1;
	double dummy;
	for (i=1; i<=numStateRatios; i++) {
		(StateRatioQ[i])->ComputeAlphaBulkRatio();
		x[j++] = log10((StateRatioQ[i])->GetAlphaBulkRatio());
	}
	if (neutralizerNeeded && !phiBulksFixed && SegQ->Charged()) {
		dummy = (MolQ->GetNeutralizer())->GetPhiBulk();
		if (dummy < 1e-300) {
			dummy = 1e-300;
		} else if (dummy > 1) {
			dummy = 0.9999999;
		}
		x[j++] = log(dummy/(1-dummy));
	}
	if (!phiBulksFixed) {
		dummy = (MolQ->GetSolvent())->GetPhiBulk();
		if (dummy < 1e-300) {
			dummy = 1e-300;
		} else if (dummy >= 1) {
			dummy = 0.9999999;
		}
		x[j++] = log(dummy/(1-dummy));
	}
}
void
SF_ReactionList::CopyItVar() {
	int i,j=1;
	for (i=1; i<=numStateRatios; i++) {
		(StateRatioQ[i])->SetAlphaBulkRatio(pow(10,x[j++]));
	}
	for (i=1; i<=numSegmentRatios; i++) {
		(SegmentRatioQ[i])->UpdateAlphaBulk();
	}
	if (neutralizerNeeded && !phiBulksFixed && SegQ->Charged()) {
//		cout << " jaaa" << endl;
		(MolQ->GetNeutralizer())->SetPhiBulk(exp(x[j])/(1+exp(x[j]))+1e-300);
//		cout << " neutral: " << (MolQ->GetNeutralizer())->GetPhiBulk() << endl;
		j++;
	}
	if (!phiBulksFixed) {
//		cout << " oook " << x[j] << endl;
		(MolQ->GetSolvent())->SetPhiBulk(exp(x[j])/(1+exp(x[j]))+1e-300);
//		cout << " solvent : " << (MolQ->GetSolvent())->GetPhiBulk() << endl;
		j++;
	}
	SegQ->UpdatePhiBulk();
}
void
SF_ReactionList::MapResiduals() {
	int i,j;
	int panicCount = 0;
	Boolean panic = false;
	Boolean done = false;
//	for (i=1; i<=numStateRatios; i++) {
//		Sysout().Outtext((StateRatioQ[i])->GetState1()->GetName());
//		Sysout().Outchar('\t');
//		Sysout().Outtext((StateRatioQ[i])->GetState2()->GetName());
//		Sysout().Outimage();
//	}
//	for (i=1; i<=MolQ->GetNumMolecules(); i++) {
//		Sysout().Outint(i,0);
//		Sysout().Outchar('\t');
//		Sysout().Outtext(MolQ->GetMolecule(i)->GetName());
//		Sysout().Outchar('\t');
//		Sysout().Outreal(MolQ->GetMolecule(i)->GetPhiBulk(),9,0);
//		Sysout().Outimage();
//	}
//	for (i=1; i<=SegQ->GetNumStates(); i++) {
//		Sysout().Outint(i,0);
//		Sysout().Outchar('\t');
//		Sysout().Outtext(SegQ->GetState(i)->GetName());
//		Sysout().Outchar('\t');
//		Sysout().Outreal(SegQ->GetState(i)->GetPhiBulk(),9,0);
//		Sysout().Outchar('\t');
//		Sysout().Outreal(SegQ->GetState(i)->GetAlphaBulk(),9,0);
//		Sysout().Outimage();
//	}
	numItVar--;
	if (neutralizerNeeded) {
		numItVar--;
	} else {
		residualMap[1] = numItVar;
		SetNeutralizingRatio();
	}
	Matrix Response(0,numItVar,0,numItVar);
	FillResponseMatrix(Response);
//	for (i=1; i<=numItVar; i++) {
//		for (j=1; j<=numItVar; j++) {
//			Sysout().Outint(i,0);
//			Sysout().Outchar('\t');
//			Sysout().Outint(j,0);
//			Sysout().Outchar('\t');
//			Sysout().Outreal(Response[i][j],9,0);
//			Sysout().Outimage();
//		}
//	}
	Vector dontTouch(1,numItVar);
	for (i=1; i<=numItVar; i++) {
		Boolean found = false;
		Boolean foundAgain = false;
		for (j=1; j<=numItVar; j++) {
			if (Response[i][j] > 0) {
				if (found) {
					foundAgain = true;
				}
				found = true;
				residualMap[i] = j;
			}
		}
		if (!foundAgain) {
			dontTouch[i] = 1;
		}
	}
	while (!done) {
		for (i=1; i<=numItVar; i++) {
			if (dontTouch[i] < 1) {
				residualMap[i] = 0;
			}
		}
		for (i=1; i<=numItVar; i++) {
			if (dontTouch[i] > 0) continue;
			Boolean abort = false;
			Boolean ok = false;
			while (!abort) {
				abort = false;
				while (!ok) {
					residualMap[i] = int(rand()*double(numItVar)/(RAND_MAX + 1.0))+1;
					ok = true;
					for (j=1; j<=numItVar; j++) {
						if (residualMap[j] == residualMap[i] && (dontTouch[j] > 0 || j<=i-1)) {
							ok = false;
						}
					}
				}
				if (!neutralizerNeeded) SetNeutralizingRatio();
				if (Response[i][residualMap[i]] <= 0 && !panic) {
					abort = true;
				} else if (Response[i][residualMap[i]] > 0) {
				} else if (!panic) { // GetResponse(i,residualMap[i]) == NaN
					abort = true;
				}
				if (abort) ok = false;
                if (ok) abort = true;
			}
			if (abort && !ok) {
				i=0;
				panicCount++;
				if (panicCount > 1000) panic = true;
			}
		}
		done = true;
		for (i=1; i<=numItVar; i++) {
			if (Response[i][residualMap[i]] <= 0 && !panic) done = false;
			else if (Response[i][residualMap[i]] > 0) ;
			else if (!panic) done = false; // GetResponse(i,residualMap[i]) == NaN
		}
	}
	if (neutralizerNeeded) {
		numItVar++;
		residualMap[numItVar] = numItVar;
	}
	numItVar++;
	residualMap[numItVar] = numItVar;
//	for (i=1; i<=numItVar; i++) {
//		Sysout().Outint(residualMap[i],0);
//		Sysout().Outimage();
//	}
}
void
SF_ReactionList::MapItVar(void) {
	int i;
	numStateRatios = 0;
	for (i=1; i<=numReactions; i++) {
		numStateRatios += (ReactionQ[i])->GetNumStateRatios();
	}
	StateRatioQ.Dim(1,numStateRatios);
	int count = 0;
	for (i=1; i<=numReactions; i++) {
		for (int j=1; j<=(ReactionQ[i])->GetNumStateRatios(); j++) {
			StateRatioQ[++count] = (ReactionQ[i])->GetStateRatio(j);
		}
	}
	SF_StateRatio* Ratio1;
	SF_StateRatio* Ratio2;
	for (i=1; i<=numStateRatios; i++) {
		Ratio1 = StateRatioQ[i];
		for (int j=i+1; j<=numStateRatios; j++) {
			Ratio2 = StateRatioQ[j];
			if (*Ratio1 == *Ratio2) {
				Array<SF_StateRatio*> RatioTemp(1,numStateRatios-1);
				int k;
				for (k=1; k<=numStateRatios; k++) {
					if (k<j) RatioTemp[k] = StateRatioQ[k];
					if (k>j) RatioTemp[k-1] = StateRatioQ[k];
				}
				numStateRatios--;
				j--;
				StateRatioQ = RatioTemp;
				for (k=1; k<=numReactions; k++) {
					(ReactionQ[k])->ReplaceRatio(Ratio2,Ratio1);
				}
				delete Ratio2;
			} else if (Ratio1->SameStates(*Ratio2)) {
				Message(fatal,"Sorry, reactions defined in peculiar order, swap the direction of a reaction and try again");
			}
		}
	}
	SF_Segment* Segment;
	numSegmentRatios = 0;
	for (i=1; i<=SegQ->GetNumSegments(); i++) {
		Segment = SegQ->GetSegment(i);
		if (Segment->GetNumStates() > 1) {
			numSegmentRatios++;
		}
	}
	SegmentRatioQ.Dim(1,numSegmentRatios);
	count = 0;
	int numRatios;
	for (i=1; i<=SegQ->GetNumSegments(); i++) {
		Segment = SegQ->GetSegment(i);
		if (Segment->GetNumStates() > 1) {
			numRatios = 0;
			int j;
			for (j=1; j<=numStateRatios; j++) {
				Ratio1 = StateRatioQ[j];
				if (SegQ->GetBaseSegment(Ratio1->GetState1()) == Segment) {
					numRatios++;
				}
			}
			Array<SF_StateRatio*> RatioQTemp(1,numRatios);
			numRatios = 0;
			for (j=1; j<=numStateRatios; j++) {
				Ratio1 = StateRatioQ[j];
				if (SegQ->GetBaseSegment(Ratio1->GetState1()) == Segment) {
					RatioQTemp[++numRatios] = Ratio1;
				}
			}
			SegmentRatioQ[++count] = new SF_StateRatioList(RatioQTemp,SegQ);
		}
	}
}
void
SF_ReactionList::SetNeutralizingRatio(void) {
	if (neutralizerNeeded || !SegQ->Charged()) return;
	int i;
	for (i=1; i<=numItVar; i++) {
		if (residualMap[i] == numItVar) {
			NeutralizingRatio = StateRatioQ[i];
//			Sysout().Outint(i,0);
//			Sysout().Outimage();
		}
	}
}
Boolean
SF_ReactionList::IterateAlphas() {
	Boolean value = false;
	for (int i=1; i<=numReactions; i++) {
		if ((ReactionQ[i])->IterateAlphas()) value = true;
	}
	return value;
}
void
SF_ReactionList::FillResponseMatrix(Matrix response) {
	int i,j;
	for (i=1; i<=numItVar; i++) {
		Boolean done = false;
		double step = 1e-5;
		while (!done) {
			done = true;
			double xValue = x[i];
			step *= 10;
			CopyItVar();
			ComputeResidues();
			for (j=1; j<=numItVar; j++) {
				response[i][j] = residues[j];
			}
			x[i] -= step*(xValue + step);
			CopyItVar();
			ComputeResidues();
			for (j=1; j<=numItVar; j++) {
				response[i][j] -= residues[j];
				response[i][j] /= step*(xValue + step);
			}
			x[i] += step*(xValue + step);
			CopyItVar();
			for (j=1; j<=numItVar; j++) {
				if (response[i][j] >= 0 || response[i][j] < 0) {
					;
				} else {
					done = false;
//					cout << "quitting" << endl;
//					exit(1);
				}
			}
		}
	}
}
void
SF_ReactionList::ComputeInternalFreeEnergies() {
//	for (int x=1; x<=SegQ->GetNumStates(); x++) {
//		SegQ->GetState(x)->SetInternalFreeEnergy(0);
//	}
//	return;
	if (internalFreeEnergiesGiven) return;
	SF_Segment* Water = FindWater();
	if (Water == NULL) {
		Message(implementation,"!Don't know how to proceed yet, thermodynamics will be wrong");
		cout << "It's possible to work around this problem by redefining your reactions by using internal energies" << endl;
		exit(1);
	}
	int i;
	SF_State* OH = NULL;
	for (i=1; i<=Water->GetNumStates(); i++) {
		SF_State* State = Water->GetState(i);
		if (State->GetValence() > 0) {
			State->SetInternalFreeEnergy(0);
		} else if (State->GetValence() < 0) {
			OH = State;
		} else {
			State->SetInternalFreeEnergy(0);
		}
	}
	for (i=1; i<=SegQ->GetNumSegments(); i++) {
		SF_Segment* Seg = SegQ->GetSegment(i);
		if (Seg->GetNumStates() > 1 && Seg != Water) {
			Seg->GetState(1)->SetInternalFreeEnergy(0);
		}
	}
	for (i=1; i<=numReactions; i++) {
		SF_Reaction* Reaction = ReactionQ[i];
		Reaction->ComputeInternalFreeEnergies();
	}
	double intFreeEnergyMinus = OH->GetInternalFreeEnergy();
	for (i=1; i<=SegQ->GetNumSegments(); i++) {
		SF_Segment* Seg = SegQ->GetSegment(i);
		if (Seg->GetNumStates() == 1) {
			SF_State* State = Seg->GetState(1);
			if (State->GetValence() < 0) {
				State->SetInternalFreeEnergy(intFreeEnergyMinus*State->GetValence());
			}
		}
	}
}
SF_Segment*
SF_ReactionList::FindWater() {
	SF_Segment* Seg = NULL;
	for (int i=1; i<=SegQ->GetNumSegments(); i++) {
		Seg = SegQ->GetSegment(i);
		if (Seg->GetNumStates() > 2) {
			return Seg;
		}
	}
	return NULL;
}
void
SF_ReactionList::ErrorOutput() {
	Text message;
	Message(literal,"Unable to find the alphas for the bulk phase\n"
		"There are two possibilities:\n"
		" 1. Your inputfile is wrong, you are trying the impossible. \n"
		" 2. The routine to find the alphas is wrong. Please report your"
		" inputfile to the author: vanmale@fenk.wau.nl \n");
	Message(literal,"Please try to compute the alphas analytically for the following");
	Array<Text> names = MyInput->GetNames("reaction");
	int i;
	for (i=1; i<=numReactions; i++) {
		message = "reaction : " + names[i] + " : equation : " +
			MyInput->GetText("reaction",names[i],"equation") +"\n";
		if (MyInput->ValueSet("reaction",names[i],"pK")) {
			message = message + "reaction : " + names[i] + " : pK : " +
				MyInput->GetText("reaction",names[i],"pK") +"\n";
		}
		if (MyInput->ValueSet("reaction",names[i],"K")) {
			message = message + "reaction : " + names[i] + " : K : " +
				MyInput->GetText("reaction",names[i],"K") +"\n";
		}
		Message(literal,message);
	}
	Text dummy;
	for (i=1; i<=numDiffStates; i++) {
		SF_State* State = StateQ[i];
		message = Copy("");
		if (State->AlphaBulkFixed()) {
			message = message + STATE + " : " + State->GetName() + " : alphabulk : ";
			dummy = Blanks(20);
			dummy.Putreal(State->GetAlphaBulkFixed(),9);
			message = message + Copy(dummy.Frontstrip()) + "\n";
		}
		if (State->PhiBulkFixed()) {
			message = message + STATE + " : " + State->GetName() + " : phibulk : ";
			dummy = Blanks(20);
			dummy.Putreal(State->GetPhiBulkFixed(),9);
			message = message + Copy(dummy.Frontstrip()) + "\n";
		}
		message = message + STATE + " : " + State->GetName() + " : valence : ";
		dummy = Blanks(20);
		dummy.Putreal(State->GetValence(),9);
		message = message + Copy(dummy.Frontstrip()) + "\n";
		Message(literal,message);
	}
	if (MolQ->GetNeutralizer() != NULL) {
		message = "neutralizer : " + (MolQ->GetNeutralizer())->GetName() + "\n";
		dummy = Blanks(20);
		dummy.Putreal((MolQ->GetNeutralizer())->GetBulkCharge(),9);
		message = message + "charge of neutralizer in bulk phase : " + Copy(dummy.Frontstrip()) + "\n";
		Message(literal,message);
	}
	Message(fatal,"Unable to find the alphas for the bulk phase, see message above");
}
void
SF_ReactionList::ComputeDerivatives(Vector ans, Vector err, int i, double hh) {
	const double CON = 1.4;
	const double CON2 = CON*CON;
	const double BIG = 1.0e30;
	const int NTAB = 25;
	const double SAFE = 1.0;
	int j;
	SF_Molecule* Mol = MolQ->GetMolecule(i);
	int numStates = SegQ->GetNumStates();
	Vector phiBulkOld(1,MolQ->GetNumMolecules());
	for (i=1; i<=MolQ->GetNumMolecules(); i++) {
		SF_Molecule* Mol = MolQ->GetMolecule(i);
		phiBulkOld[i] = Mol->GetPhiBulk();
	}
	Vector alphaBulkOld(1,numStates);
	Vector moreThanOneState(1,numStates);
	for (i=1; i<=numStates; i++) {
		SF_State* State = SegQ->GetState(i);
		alphaBulkOld[i] = State->GetAlphaBulk();
		if (SegQ->GetBaseSegment(State)->GetNumStates() > 1) {
			moreThanOneState[i] = 1;
		}
	}
	double errt,fac;
	Cube a(1,numStates,1,NTAB,1,NTAB);
	Vector finished(1,numStates);
	SetPhiBulk(Mol,phiBulkOld,alphaBulkOld,hh);
	setiterationlimit(1000);
	IterateWithFixedPhiBulks();
	if ((getaccuracy() != fabs(getaccuracy())) || (getaccuracy() > gettolerance())) {
		SetPhiBulk(Mol,phiBulkOld,alphaBulkOld,0);
		IterateWithFixedPhiBulks();
		for (i=1; i<=numStates; i++) {
			ans[i] = 0;
			err[i] = BIG;
		}
		return;
	}
	for (j=1; j<=numStates; j++) {
		if (moreThanOneState[j] > 0) {
			SF_State* State = SegQ->GetState(j);
			a[j][1][1] = State->GetAlphaBulk();
		}
	}
	SetPhiBulk(Mol,phiBulkOld,alphaBulkOld,-hh);
	IterateWithFixedPhiBulks();
	if ((getaccuracy() != fabs(getaccuracy())) || (getaccuracy() > gettolerance())) {
		SetPhiBulk(Mol,phiBulkOld,alphaBulkOld,0);
		IterateWithFixedPhiBulks();
		for (i=1; i<=numStates; i++) {
			ans[i] = 0;
			err[i] = BIG;
		}
		return;
	}
	setiterationlimit(10000);
	for (j=1; j<=numStates; j++) {
		if (moreThanOneState[j] > 0) {
			SF_State* State = SegQ->GetState(j);
			a[j][1][1] -= State->GetAlphaBulk();
			a[j][1][1] /= 2*hh;
			err[j] = BIG;
			finished[j] = 1;
		}
	}
	Boolean done = false;
	for (int k=2; (k<=NTAB) && !done; k++) {
		hh /= CON;
		SetPhiBulk(Mol,phiBulkOld,alphaBulkOld,hh);
		IterateWithFixedPhiBulks();
		done = true;
		for (j=1; j<=numStates; j++) {
			if (moreThanOneState[j] > 0 && finished[j] > 0) {
				done = false;
			}
		}
		while (getaccuracy() != fabs(getaccuracy()) || getaccuracy() > gettolerance()) {
			SetPhiBulk(Mol,phiBulkOld,alphaBulkOld,hh);
			IterateWithFixedPhiBulks();
		}
		for (j=1; j<=numStates; j++) {
			if (moreThanOneState[j] > 0) {
				SF_State* State = SegQ->GetState(j);
				a[j][1][k] = State->GetAlphaBulk();
			}
		}
		SetPhiBulk(Mol,phiBulkOld,alphaBulkOld,-hh);
		IterateWithFixedPhiBulks();
		while (getaccuracy() != fabs(getaccuracy()) || getaccuracy() > gettolerance()) {
			SetPhiBulk(Mol,phiBulkOld,alphaBulkOld,-hh);
			IterateWithFixedPhiBulks();
		}
		for (j=1; j<=numStates; j++) {
			if (moreThanOneState[j] > 0 && finished[j] > 0) {
				SF_State* State = SegQ->GetState(j);
				a[j][1][k] -= State->GetAlphaBulk();
				a[j][1][k] /= 2*hh;
				fac = CON2;
				for (int l=2; l<=k; l++) {
					a[j][l][k]=(a[j][l-1][k]*fac-a[j][l-1][k-1])/(fac-1.0);
					fac=CON2*fac;
					double z = fabs(a[j][l][k]-a[j][l-1][k]);
					double y = fabs(a[j][l][k]-a[j][l-1][k-1]);
					if (z>y) errt = z;
					else errt = y;
					if (errt <= err[j]) {
						err[j] = errt;
						ans[j] = a[j][l][k];
					}
				}
				if (fabs(a[j][k][k]-a[j][k-1][k-1]) > SAFE*err[j]) {
					finished[j] = 0;
				}
			}
		}
	}
	SetPhiBulk(Mol,phiBulkOld,alphaBulkOld,0);
	IterateWithFixedPhiBulks();
	while (getaccuracy() != fabs(getaccuracy()) || getaccuracy() > gettolerance()) {
		SetPhiBulk(Mol,phiBulkOld,alphaBulkOld,0);
		IterateWithFixedPhiBulks();
	}
}
void
SF_ReactionList::SetPhiBulk(SF_Molecule* MolRef, const Vector phiBulkOld, const Vector alphaBulkOld,  double hh) {
	int i;
	for (i=1; i<=SegQ->GetNumStates(); i++) {
		SF_State* State = SegQ->GetState(i);
		State->SetAlphaBulk(alphaBulkOld[i]);
	}
	for (i=1; i<=MolQ->GetNumMolecules(); i++) {
		SF_Molecule* Mol = MolQ->GetMolecule(i);
		if (Mol == MolRef) {
//			Sysout().Outreal(phiBulkOld[i]+(1-phiBulkOld[i])*hh,9,0);
//			Sysout().Outimage();
			Mol->SetPhiBulk(phiBulkOld[i]+(1-phiBulkOld[i])*hh);
		} else {
			if (!Mol->GetMolStructure()->SomeSegmentsPinned() && !Mol->GetMolStructure()->SomeSegmentsGrafted()) {
//				Sysout().Outreal(phiBulkOld[i]-phiBulkOld[i]*hh,9,0);
//				Sysout().Outimage();
				Mol->SetPhiBulk(phiBulkOld[i]-phiBulkOld[i]*hh);
			}
		}
	}
}
