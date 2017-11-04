#include "SF_SegmentList.h"
#include <iostream>

SF_SegmentList::SF_SegmentList(Lattice* Lat_,Input* MyInput_) {
	MyInput = MyInput_;
	Lat = Lat_;
	M = Lat->GetTotalNumLayers();
    int i,j,k,m;
    int stateCount;
    int stateCount2;
    int numstates;
    int numstates2;
    int count;

	SF_State* State;
	SF_Segment* Segment;
    Array<SegmentParam*> SegParam;
    SegmentParam* SegPar;
    SegmentParam* SegPar2;
    SegmentParam* StatePar=NULL;
    SegmentParam* StatePar2=NULL;
    Array<Text> segmentNames;
    Array<Text> stateNames;
	int numAllowedChiParam;
	int numAllowedPfParam;
	Array<Text> allowedChiParam;
	Array<Text> allowedPfParam;

	numSegments = MyInput->GetNumNames(SEGMENT);
    if (numSegments == 0) {
		Message(fatal,MyInput, "No segments (" + SEGMENT + " : ...) are defined");
    }
    SegParam.Dim(1,numSegments);
    segmentNames.Dim(1,numSegments);
    segmentNames = MyInput->GetNames(SEGMENT);
    // determine the number of states in segments
    numStates = 0;
	for (i=1; i<=numSegments; i++) {
		if (GetNumStatesDefinedInSegment(segmentNames[i]) != 0) {
			numStates += GetNumStatesDefinedInSegment(segmentNames[i]);
		} else {
			numStates++;
		}
    }
    // get the names of the states in segments
	stateNames.Dim(1,numStates);
    stateCount=0;
    Text stateId;
    Text number;
    for (i=1; i<=numSegments; i++) {
		if (GetNumStatesDefinedInSegment(segmentNames[i]) != 0) {
            count = GetNumStatesDefinedInSegment(segmentNames[i]);
            for (j=1; j<=count; j++) {
			    number = Blanks(1000);
			    number.Putint(j);
			    number = Copy(number.Frontstrip());
			    stateId = "state" + number;
			    stateNames[stateCount + j] = MyInput->GetText(SEGMENT,segmentNames[i],stateId);
			}
			if (count == 1) {
				if (*Copy(stateNames[stateCount + 1]) != *Copy(segmentNames[i])) {
					Message(fatal,MyInput, "Error in inputfile.\nSegment '" + segmentNames[i] +
						"' has one state defined,\nwith a different name: '" + stateNames[stateCount + 1] +
						"'.\nPlease use the same name, or don't define the state.");
				}
			}
			if (count > 1) {
				for (j=1; j<=count; j++) {
					Text stateName = Copy(stateNames[stateCount+j]);
					for (k=1; k<=numSegments; k++) {
						if (*stateName == *Copy(segmentNames[k])) {
							Message(fatal,MyInput, "Name clash for '"
							+ stateName + "': defined to be a segment and a state\n"
							"When defining segments with different states,\n"
							"a state name cannot be equal to a segment name.\n"
							"This is done to avoid name clashes when assigning chi parameters");
						}
					}
				}
				for (j=1; j<=count; j++) {
					Text stateName = Copy(stateNames[stateCount+j]);
					for (k=j+1; k<=count; k++) {
						Text stateName2 = Copy(stateNames[stateCount+k]);
						if (*stateName == *stateName2) {
							Message(fatal,MyInput, "In segment '" + segmentNames[i]
							+ "' at least two states ('" + stateName + "') are equal.\n");
						}
					}
				}
			}
            stateCount += count;
        }
    }
    // create segmentParam.
    for (i=1; i<=numSegments; i++) {
    	SegParam[i] = new SegmentParam(segmentNames[i], MyInput, SEGMENT);
    }
    // check if all defined states are linked to a segment
	int numDefStates = MyInput->GetNumNames(STATE);
   	Array<Text> defStates;
   	defStates.Dim(1,numDefStates);
   	defStates = MyInput->GetNames(STATE);
   	Boolean used;
    for (k=1; k<=numDefStates; k++) {
    	used = false;
    	for (i=1; i<=numSegments; i++) {
    		SegPar = SegParam[i];
    		count = SegPar->NumStates();
    		for (j=1; j<=count; j++) {
    			Text stateName = SegPar->StateName(j);
    			if (*Copy(defStates[k]) == *Copy(stateName)) {
    				used = true;
    			}
    		}
    	}
    	if (!used) {
			Message(fatal,MyInput, "A state (" + STATE + " : " + defStates[k]
			+ " : ...)\nis defined but it is never used in a segment.");
    	}
    }
	// initialize chi matrix
	chi.Dim(1,numStates,1,numStates);
	//initialize Pf matrix
	Pf.Dim(1,numSegments*numSegments);
	Pfx.Dim(1,numSegments*numSegments);
	Pfy.Dim(1,numSegments*numSegments);

	// bizarre code to read chi matrix and do various syntax checks
	stateCount = 0;
	for (i=1; i<=numSegments; i++) {
		SegPar = SegParam[i];
		numstates = SegPar->NumStates();
		if (numstates == 0) numstates = 1;
		for (j=1; j<=numstates; j++) {
			stateCount++;
			stateCount2 = 0;
			for (k=1; k<=numSegments; k++) {
				SegPar2 = SegParam[k];
				numstates2 = SegPar2->NumStates();
				if (numstates2 == 0) numstates2 = 1;
				for (m=1; m<=numstates2; m++) {
					stateCount2++;
					Boolean ChiSeg1Seg2Set = false;
					Boolean ChiSeg2Seg1Set = false;
					Boolean ChiSeg1State2Set = false;
					Boolean ChiState2Seg1Set = false;
					Boolean ChiState1Seg2Set = false;
					Boolean ChiSeg2State1Set = false;
					Boolean ChiState1State2Set = false;
					Boolean ChiState2State1Set = false;
					ChiSeg1Seg2Set = MyInput->ValueSet(SEGMENT,SegPar->Name(),"chi - " + SegPar2->Name());
					ChiSeg2Seg1Set = MyInput->ValueSet(SEGMENT,SegPar2->Name(),"chi - " + SegPar->Name());
					if (SegPar2->NumStates() != 0) {
						StatePar2 = SegPar2->StatePTR(m);
						ChiSeg1State2Set = MyInput->ValueSet(SEGMENT,SegPar->Name(),"chi - " + StatePar2->Name());
						ChiState2Seg1Set = MyInput->ValueSet(STATE,StatePar2->Name(),"chi - " + SegPar->Name());
					}
					if (SegPar->NumStates() != 0) {
						StatePar = SegPar->StatePTR(j);
						ChiState1Seg2Set = MyInput->ValueSet(STATE,StatePar->Name(),"chi - " + SegPar2->Name());
						ChiSeg2State1Set = MyInput->ValueSet(SEGMENT,SegPar2->Name(),"chi - " + StatePar->Name());
						if (SegPar2->NumStates() != 0) {
							StatePar2 = SegPar2->StatePTR(m);
							ChiState1State2Set = MyInput->ValueSet(STATE,StatePar->Name(),"chi - " + StatePar2->Name());
							ChiState2State1Set = MyInput->ValueSet(STATE,StatePar2->Name(),"chi - " + StatePar->Name());
						}
					}
					Text mess;
					if ((ChiSeg1Seg2Set || ChiSeg1State2Set || ChiState1Seg2Set || ChiState1State2Set) &&
						(ChiSeg2Seg1Set || ChiSeg2State1Set || ChiState2Seg1Set || ChiState2State1Set)) {
						mess = "Asymmetric chi values, unable to set both\n'";
						if (ChiSeg1Seg2Set || ChiSeg1State2Set) {
							mess = mess + SEGMENT + " : " + SegPar->Name() + " : chi - ";
						} else {
							mess = mess + STATE + " : " + StatePar->Name() + " : chi - ";
						}
						if (ChiSeg1Seg2Set || ChiState1Seg2Set) {
							mess = mess + SegPar2->Name();
						} else {
							mess = mess + StatePar2->Name();
						}
						mess = mess + "' and\n'";
						if (ChiSeg2Seg1Set || ChiSeg2State1Set) {
							mess = mess + SEGMENT + " : " + SegPar2->Name() + " : chi - ";
						} else {
							mess = mess + STATE + " : " + StatePar2->Name() + " : chi - ";
						}
						if (ChiSeg2Seg1Set || ChiState2Seg1Set) {
							mess = mess + SegPar->Name();
						} else {
							mess = mess + StatePar->Name();
						}
						mess = mess + "'";
						Message(fatal,MyInput,mess);
					}
					if (ChiSeg1Seg2Set && ChiSeg1State2Set) {
						mess = "Unable to set both\n'";
						mess = mess + SEGMENT + " : " + SegPar->Name() + " : chi - " + SegPar2->Name() + "'and\n'";
						mess = mess + SEGMENT + " : " + SegPar->Name() + " : chi - " + StatePar2->Name() + "'";
						Message(fatal,MyInput,mess);
					}
					if (ChiSeg1Seg2Set && ChiState1Seg2Set) {
						mess = "Unable to set both\n'";
						mess = mess + SEGMENT + " : " + SegPar->Name() + " : chi - " + SegPar2->Name() + "'and\n'";
						mess = mess + STATE + " : " + StatePar->Name() + " : chi - " + SegPar2->Name() + "'";
						Message(fatal,MyInput,mess);
					}
					if (ChiSeg1Seg2Set && ChiState1State2Set) {
						mess = "Unable to set both\n'";
						mess = mess + SEGMENT + " : " + SegPar->Name() + " : chi - " + SegPar2->Name() + "'and\n'";
						mess = mess + STATE + " : " + StatePar->Name() + " : chi - " + StatePar2->Name() + "'";
						Message(fatal,MyInput,mess);
					}
					if (ChiState1Seg2Set && ChiSeg1State2Set) {
						mess = "Unable to set both\n'";
						mess = mess + STATE + " : " + StatePar->Name() + " : chi - " + SegPar2->Name() + "'and\n'";
						mess = mess + SEGMENT + " : " + SegPar->Name() + " : chi - " + StatePar2->Name() + "'";
						Message(fatal,MyInput,mess);
					}
					if (ChiState1Seg2Set && ChiState1State2Set) {
						mess = "Unable to set both\n'";
						mess = mess + STATE + " : " + StatePar->Name() + " : chi - " + SegPar2->Name() + "'and\n'";
						mess = mess + STATE + " : " + StatePar->Name() + " : chi - " + StatePar2->Name() + "'";
						Message(fatal,MyInput,mess);
					}
					if (ChiSeg1State2Set && ChiState1State2Set) {
						mess = "Unable to set both\n'";
						mess = mess + SEGMENT + " : " + SegPar->Name() + " : chi - " + StatePar2->Name() + "'and\n'";
						mess = mess + STATE + " : " + StatePar->Name() + " : chi - " + StatePar2->Name() + "'";
						Message(fatal,MyInput,mess);
					}
					if (ChiSeg1Seg2Set) chi[stateCount][stateCount2] = MyInput->GetReal(SEGMENT,SegPar->Name(),"chi - " + SegPar2->Name(),-DBL_MAX,DBL_MAX);
					if (ChiSeg1State2Set) chi[stateCount][stateCount2] = MyInput->GetReal(SEGMENT,SegPar->Name(),"chi - " + StatePar2->Name(),-DBL_MAX,DBL_MAX);
					if (ChiState1Seg2Set) chi[stateCount][stateCount2] = MyInput->GetReal(STATE,StatePar->Name(),"chi - " + SegPar2->Name(),-DBL_MAX,DBL_MAX);
					if (ChiState1State2Set) chi[stateCount][stateCount2] = MyInput->GetReal(STATE,StatePar->Name(),"chi - " + StatePar2->Name(),-DBL_MAX,DBL_MAX);
					if (chi[stateCount][stateCount2] != 0) {
						chi[stateCount2][stateCount] = chi[stateCount][stateCount2];
					}
				}
			}
		}
	}
	int segCount=0,segCount2;
	Pf_count=0;
	for (i=1; i<=numSegments; i++) {
		segCount++;
		SegPar = SegParam[i];
		segCount2=0;
		for (k=1; k<=numSegments; k++) {
			segCount2++;
			SegPar2 = SegParam[k];
			Boolean PfSeg1Seg2Set = false;
			Boolean PfSeg2Seg1Set = false;

			PfSeg1Seg2Set = MyInput->ValueSet(SEGMENT,SegPar->Name(),"Pf - " + SegPar2->Name());
			PfSeg2Seg1Set = MyInput->ValueSet(SEGMENT,SegPar2->Name(),"Pf - " + SegPar->Name());

			if (PfSeg1Seg2Set && PfSeg2Seg1Set) {
				Pf_count++;
				if (MyInput->GetReal(SEGMENT,SegPar->Name(),"Pf - " + SegPar2->Name(),0,1)!=
					MyInput->GetReal(SEGMENT,SegPar2->Name(),"Pf - " + SegPar->Name(),0,1))
					Message(fatal,"Pf values should be symmetric. Found asymmetric values instead") ;
				Pf[Pf_count] = MyInput->GetReal(SEGMENT,SegPar->Name(),"Pf - " + SegPar2->Name(),0,1);
				Pfx[Pf_count]=SegPar->Name();
				Pfy[Pf_count]=SegPar2->Name();
			}
			if (PfSeg2Seg1Set && !PfSeg1Seg2Set) {
				Pf_count++;
				Pf[Pf_count] = MyInput->GetReal(SEGMENT,SegPar2->Name(),"Pf - " + SegPar->Name(),0,1);
				Pfx[Pf_count]=SegPar2->Name();
				Pfy[Pf_count]=SegPar->Name();
			}
			if (PfSeg1Seg2Set && !PfSeg2Seg1Set) {
				Pf_count++;
				Pf[Pf_count] = MyInput->GetReal(SEGMENT,SegPar->Name(),"Pf - " + SegPar2->Name(),0,1);
				Pfx[Pf_count]=SegPar->Name();
				Pfy[Pf_count]=SegPar2->Name();
			}
		}
	}

	// check if diagonal elements of chi matrix are zero.
	stateCount = 0;
	for (i=1; i<=numSegments; i++) {
		SegPar = SegParam[i];
		numstates = SegPar->NumStates();
		if (numstates == 0) numstates = 1;
		for (j=0; j<numstates; j++) {
			stateCount++;
			if (chi[stateCount][stateCount] != 0) {
				Message(fatal,MyInput,"Programming error in class SegmentParam, diagonal element of chi matrix is not zero.");
    		}
		}
	}

	// get possible input chi parameters for validation input parameters at creation of class SF_Segment
	// get possible input Pf parameters for validation input parameters at creation of class SF_Segment
	numAllowedChiParam = 0;
	numAllowedPfParam = 0;
	for (i=1; i<=numSegments; i++) {
		SegPar = SegParam[i];
		numAllowedChiParam++;
		numAllowedPfParam++;
		numstates = SegPar->NumStates();
		numAllowedChiParam += numstates;
	}
	int chiParamCount = 0;
	int PfParamCount = 0;
	allowedChiParam.Dim(1,numAllowedChiParam);
	allowedPfParam.Dim(1,numAllowedPfParam);
	for (i=1; i<=numSegments; i++) {
		SegPar = SegParam[i];
		allowedChiParam[++chiParamCount] = SegPar->Name();
		allowedPfParam[++PfParamCount] = SegPar->Name();
		numstates = SegPar->NumStates();
		for (j=1; j<=numstates; j++) {
			allowedChiParam[++chiParamCount] = SegPar->StateName(j);
		}
	}
	// create segments and dump in SegmentQ
	SegmentQ.Dim(1,numSegments);
	for (i=1; i<=numSegments; i++) {
		// hier keuze surface/normal maken
		// + keuze chemisch verschillende states of niet, neen, test later
		SegmentQ[i] = new SF_Segment(segmentNames[i],allowedChiParam,allowedPfParam,Lat,MyInput);
	}
	// read states and dump in StateQ
	StateQ.Dim(1,numStates);
	stateCount = 0;
	for (i=1; i<=numSegments; i++) {
		Segment = SegmentQ[i];
		if (GetNumStatesDefinedInSegment(segmentNames[i]) != 0) {
            count = GetNumStatesDefinedInSegment(segmentNames[i]);
			for (j=1; j<=count; j++) {
				StateQ[stateCount+j] = Segment->GetState(j);
			}
			stateCount += count;
		} else {
			StateQ[++stateCount] = Segment->GetState(1);
		}
	}
	// check if surfaces are defined in Lattice.
	int numFrozen = 0;
	for (i=1; i<=numSegments; i++) {
		Segment = SegmentQ[i];
		if (Segment->GetFreedom() == frozen) {
			numFrozen++;
		}
	}
	Array<LatticeRange*> FrozenRanges;
	FrozenRanges.Dim(1,numFrozen);
	count = 0;
	for (i=1; i<=numSegments; i++) {
		Segment = SegmentQ[i];
		if (Segment->GetFreedom() == frozen) {
			FrozenRanges[++count] = Segment->GetLatRange();
		}
	}
	Lat->CheckBoundaries(FrozenRanges);
	for (i=1; i<=numSegments; i++) {
		delete SegParam[i];
	}
	charged = false;
	for (i=1; i<=numStates; i++) {
		State = StateQ[i];
		if (State->GetValence() != 0) charged = true;
	}
	if (charged) {
		psi.Dim(1,M);
		negValence = GetNegValence();
		posValence = GetPosValence();
		ComputeNumPos();
		ComputeNumNeg();
		if (Lat->GetSiteDistance() > 0.001) {
			Message(warning,"'lat : " + Lat->GetName() + " : bondlength' "
			"or 'distance' is rather high for a charged system, it is "
			"usually in the order of 1e-9. Your system may not become neutral.");
		}
	}
}
SF_SegmentList::~SF_SegmentList() {
	for (int i=1; i<=numSegments; i++) {
		delete SegmentQ[i];
	}
}
void
SF_SegmentList::GetOutput(Output* Out) const {
	int i,j;
	SF_State* State1;
	SF_State* State2;
	for (i=1; i<=numStates; i++) {
		State1 = StateQ[i];
		for (j=1; j<=numStates; j++) {
			State2 = StateQ[j];
			if (i!=j) {
				Out->PutReal("chi list",
							 State1->GetName(),
							 "chi - " + State2->GetName(),
							 chi[i][j]);
			}
		}
	}
	for (i=1; i<=numSegments; i++) {
		(SegmentQ[i])->GetOutput(Out);
	}
	for (i=1; i<=GetPf_count(); i++){{
			if (Pf[i] !=0) Out->PutReal("Pf list",Pfx[i], "Pf - " + Pfy[i],Pf[i]);
		}
	}
}
Boolean
SF_SegmentList::SegmentDefined(const Text name) const {
	for (int i=1; i<=numSegments; i++) {
		if (*Copy((SegmentQ[i])->GetName()) == *name) {
			return true;
		}
	}
	return false;
}
int
SF_SegmentList::GetNumSegments() const {
	return numSegments;
}
SF_Segment*
SF_SegmentList::GetSegment(const int number) const {
	return SegmentQ[number];
}
SF_Segment*
SF_SegmentList::GetSegment(const Text name) const {
	for (int i=1; i<=numSegments; i++) {
		if (*Copy((SegmentQ[i])->GetName()) == *name) {
			return SegmentQ[i];
		}
	}
	Message(fatal,"Programming error, SF_SegmentList::GetSegment(Text name) segment does not exist");
	return NULL; // never get here
}
SF_Segment*
SF_SegmentList::GetBaseSegment(const SF_State* State) const {
	int numSegStates;
	SF_Segment* Segment;
	for (int i=1; i<=numSegments; i++) {
		Segment = SegmentQ[i];
		numSegStates = Segment->GetNumStates();
		for (int j=1; j<=numSegStates; j++) {
			if (State == Segment->GetState(j)) {
				return Segment;
			}
		}
	}
	Message(fatal,"Programming error, SF_SegmentList::GetBaseSegment(SF_State* State) state does not exist");
	return NULL; // never get here
}
SF_Segment*
SF_SegmentList::GetBaseSegment(const SF_MolSegment* MolSeg) const {
	int numMolSegments;
	SF_Segment* Segment;
	for (int i=1; i<=numSegments; i++) {
		Segment = SegmentQ[i];
		numMolSegments = Segment->GetNumMolSegments();
		for (int j=1; j<=numMolSegments; j++) {
			if (MolSeg == Segment->GetMolSegment(j)) {
				return Segment;
			}
		}
	}
	Message(fatal,"Programming error, SF_SegmentList::GetBaseSegment(SF_MolSegment* MolSeg) state does not exist");
	return NULL; // never get here
}
int
SF_SegmentList::GetNumStates() const {
	return numStates;
}
SF_State*
SF_SegmentList::GetState(const int number) const {
	return StateQ[number];
}
SF_State*
SF_SegmentList::GetState(const Text name) const {
	for (int i=1; i<=numStates; i++) {
		if (*Copy((StateQ[i])->GetName()) == *name) {
			return StateQ[i];
		}
	}
	Message(fatal,"Programming error, SF_SegmentList::GetState(Text name) state does not exist");
	return NULL; // never get here
}
Boolean
SF_SegmentList::StateDefined(const Text name) const {
	for (int i=1; i<=numStates; i++) {
		if (*Copy((StateQ[i])->GetName()) == *name) {
			return true;
		}
	}
	return false;
}
SF_MolSegment*
SF_SegmentList::NewMolSegment(Text segName,Text molName) {
	SF_Segment* Segment;
	for (int i=1; i<=numSegments; i++) {
		Segment = SegmentQ[i];
		if (*Copy(Segment->GetName()) == *segName) {
			return Segment->NewMolSegment(molName);
		}
	}
	Message(fatal,MyInput,"Cannot use segment '" + segName + "' for 'mol : " +
	molName + " : composition', \nsegment is not defined.");
	return NULL; // never get here
}

void
SF_SegmentList::PutSides(const SF_State* State) const {
	int chiParamNum = GetChiParamNum(State);
	SF_State* State2;

	for (int i=1; i<=numStates; i++) {
		double chi12 = chi[chiParamNum][i];
		if (chi12 != 0) {
			State2 = StateQ[i];
			if (!State2->RhoSet()) {
				if (Lat->GetNumGradients()==3) {
					State2->PutRHO();
				} else {
					Vector phi=State2->GetPhi();
					Lat->SetBoundaries(phi);
					for (int z=1; z<=M; z++) {
						State2->PutRHO(Lat->SideFraction(phi,z),z);
					}
				}
			}
		}
	}
}

double
SF_SegmentList::ChemInt(const SF_State* State, const int z) const {
	double chemInt = 0;
	int chiParamNum = GetChiParamNum(State);
	SF_State* State2;
	for (int i=1; i<=numStates; i++) {
		double chi12 = chi[chiParamNum][i];
		if (chi12 != 0) {
			State2 = StateQ[i];
			if (!State2->RhoSet()) {
				if (Lat->GetNumGradients()==3) {
					State2->PutRHO();
				} else {
					Vector phi=State2->GetPhi();
					for (int z=1; z<=M; z++) {
						State2->PutRHO(Lat->SideFraction(phi,z),z);
					}
				}
			}
			chemInt += chi12*State2->GetRHO(z);
			//chemInt += chi12 * Lat->SideFraction(State2->GetPhi(),z);
		}
	}
	return chemInt;
}

double
SF_SegmentList::ChemInt(const SF_MolState* State, const int z) const {
	double chemInt = 0;
	int chiParamNum = GetChiParamNum(State);
	SF_State* State2;
	for (int i=1; i<=numStates; i++) {
		double chi12 = chi[chiParamNum][i];
		if (chi12 != 0) {
			State2 = StateQ[i];
			if (!State2->RhoSet()) {
				if (Lat->GetNumGradients()==3) {
					State2->PutRHO();
				} else {
					Vector phi=State2->GetPhi();
					for (int z=1; z<=M; z++) {
						State2->PutRHO(Lat->SideFraction(phi,z),z);
					}
				}
			}
			chemInt += chi12*State2->GetRHO(z);
			//chemInt += chi12 * Lat->SideFraction(State2->GetPhi(),z);
		}
	}
	return chemInt;
}
double
SF_SegmentList::ChemInt(const SF_Segment* Segment, const int z) const {
	double chemInt = 0;
	int numStatesSeg = Segment->GetNumStates();
	SF_State* State;
	for (int i=1; i<=numStatesSeg; i++) {
		State = Segment->GetState(i);
		chemInt += ChemInt(State,z);
	}
	return chemInt;
}

double
SF_SegmentList::ChemInt(const SF_MolSegment* Segment, const int z) const {
	double chemInt = 0;
	int numStatesSeg = Segment->GetNumStates();
	SF_MolState* State;
	for (int i=1; i<=numStatesSeg; i++) {
		State = Segment->GetState(i);
		chemInt += ChemInt(State,z);
	}
	return chemInt;
}
double
SF_SegmentList::ChemIntBulk(const SF_State* State) const {
	if (State->GetFreedom() == frozen) return 0;
	double chemInt = 0;
	int chiParamNum = GetChiParamNum(State);
	SF_State* State2;
	for (int i=1; i<=numStates; i++) {
		double chi12 = chi[chiParamNum][i];
		if (chi12 != 0) {
			State2 = StateQ[i];
			chemInt += chi12 * State2->GetPhiBulk();
		}
	}
	return chemInt;
}
double
SF_SegmentList::ChemIntBulk(const SF_MolState* State) const {
	double chemInt = 0;
	int chiParamNum = GetChiParamNum(State);
	SF_State* State2;
	for (int i=1; i<=numStates; i++) {
		double chi12 = chi[chiParamNum][i];
		if (chi12 != 0) {
			State2 = StateQ[i];
			chemInt += chi12 * State2->GetPhiBulk();
		}
	}
	return chemInt;
}
double
SF_SegmentList::ChemIntBulk(const SF_Segment* Segment) const {
	if (Segment->GetFreedom() == frozen) return 0;
	double chemInt = 0;
	int numStatesSeg = Segment->GetNumStates();
	SF_State* State;
	for (int i=1; i<=numStatesSeg; i++) {
		State = Segment->GetState(i);
		chemInt += ChemIntBulk(State);
	}
	return chemInt;
}
double
SF_SegmentList::ChemIntBulk(const SF_MolSegment* Segment) const {
	double chemInt = 0;
	int numStatesSeg = Segment->GetNumStates();
	SF_MolState* State;
	for (int i=1; i<=numStatesSeg; i++) {
		State = Segment->GetState(i);
		chemInt += ChemIntBulk(State);
	}
	return chemInt;
}
double
SF_SegmentList::ChemIntRef(const SF_MolState* MolState) const {
	double chemInt = 0;
	int chiParamNum = GetChiParamNum(MolState);
	SF_State* State2;
	SF_MolState* MolState2;
	int numMolStates;
	for (int i=1; i<=numStates; i++) {
		if (chi[chiParamNum][i] != 0) {

			State2 = StateQ[i];
			numMolStates = State2->GetNumMolStates();
			for (int j=1; j<=numMolStates; j++) {
				MolState2 = State2->GetMolState(j);
				if (*Copy(MolState->GetMolName()) == *Copy(MolState2->GetMolName())) {
					chemInt += chi[chiParamNum][i] * MolState2->GetPhiRef();
				}
			}
		}
	}
	return chemInt;
}
double
SF_SegmentList::ChemIntRef(const SF_MolSegment* MolSegment) const {
	return ChemIntRef(MolSegment->GetState(1));
}
Boolean
SF_SegmentList::ChemIntStatesEqual(const SF_State* State1, const SF_State* State2) const {
	int i;
	int stateNum1 = GetChiParamNum(State1);
	int stateNum2 = GetChiParamNum(State2);
	Boolean equal = true;
	for (i=1; i<=numStates; i++) {
		if (chi[i][stateNum1] != chi[i][stateNum2]) {
			equal = false;
		}
	}
	return equal;
}
void
SF_SegmentList::UpdateSWF() {
	SF_Segment* Segment;
	for (int i=1; i<=numSegments; i++) {
		Segment = SegmentQ[i];
		if (Segment->GetFreedom() != frozen) {
			Segment->UpdateSWF();
		} else {
			Segment->ResetRhoSet();
		}

	}
}
void
SF_SegmentList::UpdatePhiBulk() {
	SF_Segment* Segment;
	for (int i=1; i<=numSegments; i++) {
		Segment = SegmentQ[i];
		if (Segment->GetFreedom() != frozen) {
			Segment->UpdatePhiBulk();
		}
	}
}
void
SF_SegmentList::UpdatePhiStates() {
	SF_Segment* Segment;
	for (int i=1; i<=numSegments; i++) {
		Segment = SegmentQ[i];
		Segment->UpdatePhiStates();
	}
}
Vector
SF_SegmentList::GetPhiTotal() const {
	Vector phi;
	Vector phiTotal(1,M);
	SF_State* State;
	for (int i=1; i<=numStates; i++) {
		State = StateQ[i];
		phi = State->GetPhi();
		for (int z=1; z<=M; z++) {
			phiTotal[z] += phi[z];
		}
	}
	return phiTotal;
}
Boolean
SF_SegmentList::Charged() const {
	return charged;
}
Boolean
SF_SegmentList::ReactionsPresent() const {
	SF_Segment* Segment;
	for (int i=1; i<= numSegments; i++) {
		Segment = SegmentQ[i];
		if (Segment->GetNumStates() > 1) return true;
	}
	return false;
}
Vector
SF_SegmentList::GetCharge() const {
	Vector phi;
	Vector phiTotal(1,M);
	Vector charge(1,M);
	if (!charged) {
		return charge;
	}
	SF_State* State;
	double valence;
	for (int i=1; i<=numStates; i++) {
		State = StateQ[i];
		valence = State->GetValence();
		phi = State->GetPhi();
		for (int z=1; z<=M; z++) {
			charge[z] += phi[z]*valence;
			phiTotal[z] += phi[z];
		}
	}
	for (int z=1; z<=M; z++) {
		charge[z] /= phiTotal[z];
	}
	Lat->MultiplyWithLatticeSites(charge);
	return charge;
}
Vector
SF_SegmentList::GetAverageEpsilon() { //const?
	Vector epsilon(1,M);
	SF_State* State;
	UpdatePhiStates(); //kan weg?
	Vector phi;
	double eps;
	for (int i=1; i<=numStates; i++) {
		State = StateQ[i];
		phi = State->GetPhi();
		eps = State->GetEpsilon();
		for (int z=1; z<=M; z++) {
			epsilon[z] += phi[z]*eps;
		}
	}
	return epsilon;
}
Vector
SF_SegmentList::GetElectricPotential() const {
	if (!charged) {
		Vector A(1,M);
		return A;
	} else {
		return psi;
	}
}
Vector
SF_SegmentList::ComputeElectricPotential() {
	if (!charged) {
		Vector A(1,M);
		return A;
	} else {
		Vector psi0(1,M);
		Vector psiold(1,M);
		Vector posEnergy = GetPosEnergy();
		Vector negEnergy = GetNegEnergy();
		for (int z=1; z<=M; z++) {
			psi0[z] = (BOLTZMANN*TEMPERATURE/ELEM_CHARGE)*(posEnergy[z]/numPos - negEnergy[z]/numNeg)
				/(posValence/numPos-negValence/numNeg);
		}
		Array<Text> solveNames = MyInput->GetNames("newton");
		Text solveName = solveNames[1];
		int numchargeiter = MyInput -> GetInt("newton",solveName,"numchargeiter",1,INT_MAX,1);
		Vector epsAv = GetAverageEpsilon();
		Vector charge = GetCharge();
		double preFactor = (ELEM_CHARGE/EPS0)*Lat->GetSiteDistance()/Lat->GetSiteSurface();
		Lat->ElectricPotential(psi,psi0,epsAv,charge,preFactor);
		preFactor = (ELEM_CHARGE/EPS0)/(Lat->GetSiteDistance());
		for (int i=2; i <= numchargeiter;i++){
			for (int j=1; j<=M; j++) psiold[j] = psi[j];
			Lat->ElectricPotential(psi,psiold,epsAv,charge,preFactor);
		}
	}
	return psi;
}
Vector
SF_SegmentList::ComputeElectricPotential(Vector psi0) {


	if (!charged) {
		Vector A(1,M);
		return A;
	} else {
		Vector psiold(1,M);
		Vector epsAv = GetAverageEpsilon();
		Vector charge = GetCharge();
		Array<Text> solveNames = MyInput->GetNames("newton");
		Text solveName = solveNames[1];
		int numchargeiter = MyInput -> GetInt("newton",solveName,"numchargeiter",1,INT_MAX,1);
		double preFactor = (ELEM_CHARGE/EPS0)/(Lat->GetSiteDistance());
		Lat->ElectricPotential(psi,psi0,epsAv,charge,preFactor);
		for (int i=2; i <= numchargeiter;i++){
			for (int j=1; j<=M; j++) psiold[j] = psi[j];
			Lat->ElectricPotential(psi,psiold,epsAv,charge,preFactor);
		}
	}
	return psi;
}
Matrix
SF_SegmentList::GetChiMatrix() const {
	return chi;
}

Vector
SF_SegmentList::GetPfVector() const {
	return Pf;
}

Array <Text>
SF_SegmentList::GetPfx() const{
	return Pfx;
}

Array<Text>
SF_SegmentList::GetPfy() const{
	return Pfy;
}

int
SF_SegmentList::GetPf_count() const {
	return Pf_count;
}

void
SF_SegmentList::DeleteUnusedSegments() {
	SF_Segment* Segment;
	SF_Segment* SegmentToBeDeleted;
	int i,j,k,l,m;
	for (i=1; i<=numSegments; i++) {
		Segment =  SegmentQ[i];
		if (Segment->GetFreedom() != frozen
		&& Segment->GetNumMolSegments() == 0) {
			Message(warning,SEGMENT + " : " + Segment->GetName()
				+ " is not used in a molecule nor is it used as a 'frozen' segment in the system. It will be deleted");
			SegmentToBeDeleted = Segment;
			Array<SF_Segment*> SegmentQTemp(1,numSegments-1);
			for (j=1; j<=numSegments; j++) {
				if (j<i) SegmentQTemp[j] = SegmentQ[j];
				if (j>i) SegmentQTemp[j-1] = SegmentQ[j];
			}
			numSegments--;
			i--;
			SegmentQ = SegmentQTemp;
			SF_State* State;
			SF_State* SegState;
			int numSegStates = SegmentToBeDeleted->GetNumStates();
			for (j=1; j<=numSegStates; j++) {
				SegState = SegmentToBeDeleted->GetState(j);
				for (k=1; k<=numStates; k++) {
					State = StateQ[k];
					Array<SF_State*> StateQTemp(1,numStates-1);
					Matrix chiTemp(1,numStates-1,1,numStates-1);
					if (State == SegState) {
						for (l=1; l<= numStates; l++) {
							if (l<k) StateQTemp[l] = StateQ[l];
							if (l>k) StateQTemp[l-1] = StateQ[l];
						}
						for (l=1; l<=numStates; l++) {
							for (m=1; m<=numStates; m++) {
								if (l<k && m<k) chiTemp[l][m] = chi[l][m];
								if (l<k && m>k) chiTemp[l][m-1] = chi[l][m];
								if (l>k && m<k) chiTemp[l-1][m] = chi[l][m];
								if (l>k && m>k) chiTemp[l-1][m-1] = chi[l][m];
							}
						}
						numStates--;
						k--;
						StateQ = StateQTemp;
						chi = chiTemp;
					}
				}
			}
			delete SegmentToBeDeleted;
		}
	}
	if (charged) {
		negValence = GetNegValence();
		posValence = GetPosValence();
		ComputeNumPos();
		ComputeNumNeg();
	}
}
Vector
SF_SegmentList::GetPosEnergy() const {
	Vector posPotential(1,M);
	SF_State* State;
	LatticeRange* Range;
	Vector SWF;
	for (int i=1; i<=numStates; i++) {
		State = StateQ[i];
		if (State->GetFreedom() != frozen && State->GetValence() > 0) {
			SWF = State->GetSWF();
			Range = State->GetLatRange();
			for (int z=1; z<=M; z++) {
				if (Range->InRange(z) && SWF[z] > 0) {
					posPotential[z] += -ChemInt(State,z) -log(SWF[z]);
				}
			}
		}
	}
	return posPotential;
}
Vector
SF_SegmentList::GetNegEnergy() const {
	Vector negPotential(1,M);
	SF_State* State;
	LatticeRange* Range;
	Vector SWF;
	for (int i=1; i<=numStates; i++) {
		State = StateQ[i];
		if (State->GetFreedom() != frozen && State->GetValence() < 0) {
			SWF = State->GetSWF();
			Range = State->GetLatRange();
			for (int z=1; z<=M; z++) {
				if (Range->InRange(z) && SWF[z] > 0) {
					negPotential[z] += -ChemInt(State,z) -log(SWF[z]);
				}
			}
		}
	}
	return negPotential;
}
double
SF_SegmentList::GetPosValence() const {
	double totValence = 0;
	SF_State* State;
	for (int i=1; i<=numStates; i++) {
		State = StateQ[i];
		if (State->GetFreedom() != frozen && State->GetValence() > 0) {
			totValence += State->GetValence();
		}
	}
	return totValence;
}
double
SF_SegmentList::GetNegValence() const {
	double totValence = 0;
	SF_State* State;
	for (int i=1; i<=numStates; i++) {
		State = StateQ[i];
		if (State->GetFreedom() != frozen && State->GetValence() < 0) {
			totValence += State->GetValence();
		}
	}
	return totValence;
}
void
SF_SegmentList::ComputeNumPos() {
	numPos = 0;
	SF_State* State;
	for (int i=1; i<=numStates; i++) {
		State = StateQ[i];
		if (State->GetFreedom() != frozen && State->GetValence() > 0) {
			numPos++;
		}
	}
}
void
SF_SegmentList::ComputeNumNeg() {
	numNeg = 0;
	SF_State* State;
	for (int i=1; i<=numStates; i++) {
		State = StateQ[i];
		if (State->GetFreedom() != frozen && State->GetValence() < 0) {
			numNeg++;
		}
	}
}
int
SF_SegmentList::GetNumStatesDefinedInSegment(const Text name) const {
	Text number;
	Text stateId;
	int i=1;
	stateId = "state1";
	while (MyInput->ValueSet(SEGMENT,name,stateId)) {
		number = Blanks(100);
		number.Putint(++i);
		number = Copy(number.Frontstrip());
		stateId = "state" + number;
	}
	return i-1;
}
int
SF_SegmentList::GetChiParamNum(const SF_State* State) const {
	for (int i=1; i<=numStates; i++) {
		if (StateQ[i] == State) return i;
	}
	Message(fatal,"programming error in call to "
		"SF_SegmentList::GetChiParamNum(SF_State* State)");
	return -1; // never get here
}

int
SF_SegmentList::GetChiParamNum(const SF_MolState* MolState) const {
	SF_State* State;
	for (int i=1; i<=numStates; i++) {
		State = StateQ[i];
		for (int j=1; j<=State->GetNumMolStates(); j++) {
			if (State->GetMolState(j) == MolState) return i;
		}
	}
	Message(fatal,"programming error in call to "
		"SF_SegmentList::GetChiParamNum(SF_MolState* MolState)");
	return -1; // never get here
}



