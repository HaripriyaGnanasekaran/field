#include "SF_Segment.h"
#include <iostream>

SF_Segment::SF_Segment(Text name_,
					   Array<Text> chiParam,
					   Array<Text> PfParam,
                       Lattice* Lat_,
					   Input* MyInput_) {
	name = name_;
	Lat = Lat_;
    MyInput = MyInput_;
	int numChiParam = chiParam.Upperbound() - chiParam.Lowerbound() +1;
	int numPfParam = PfParam.Upperbound() - PfParam.Lowerbound()+1;
    Text stateId;
    Text number;
	stateId = "state1";
	int i=1;
	while (MyInput->ValueSet(SEGMENT,name,stateId)) {
		number = Blanks(100);
		number.Putint(++i);
		number = Copy(number.Frontstrip());
		stateId = "state" + number;
	}
	numStates = i-1;
	SpotSize =1; //default value;
	SpotType =2; //default value
	int numParam = numChiParam+numPfParam+numStates;
	if (numStates == 0) numParam++;
	Array <Text> param;
	param.Dim(1,numParam+13);
	for (i=1; i<=numChiParam; i++) {
		param[i] = "chi - " + chiParam[i];

	}
	int j=1;
	i--;
	for (j=1; j<=numPfParam; j++) {
		param[i+j] = "Pf - " + PfParam[j];
	}
	for (i=numChiParam+numPfParam+1; i<=numChiParam + numPfParam + numStates; i++) {
		number = Blanks(100);
		number.Putint(i-numChiParam-numPfParam);
		number = Copy(number.Frontstrip());
		param[i] = "state" + number;
	}
	param[numParam+1] = "pinned_range";
	param[numParam+2] = "grafted_range";
	param[numParam+3] = "frozen_range";
	param[numParam+4] = "valence";
	param[numParam+5] = "epsilon";
	param[numParam+6] = "freedom";
	param[numParam+7] = "MC";
	param[numParam+8] = "LD";
	param[numParam+9] = "spot_size";
	param[numParam+10] = "spot_type";
	param[numParam+11] = "move";
	param[numParam+12] = "spot_mix";
	param[numParam+13] = "spot_frac";
	//param[numParam+12] = "alpha_Mayer_Saupe";
	MyInput->CheckParameterNames(SEGMENT,name,param);
	Array<Text> freedoms(1,4);
	freedoms[1] = "free";
	freedoms[2] = "pinned";
	freedoms[3] = "grafted";
	freedoms[4] = "frozen";
    int freedomChoice = MyInput->GetChoice(SEGMENT,name,"freedom",freedoms,1);
	Array<Text> latname;
	latname = MyInput->GetNames("lat");
	int lattyp = 1; // default standard lattice
	if (MyInput->ValueSet("lat",latname[1],"latticetype")) {
		Array<Text> lattype(1,4);
		lattype[1] = "standard";
		lattype[2] = "stencils";
		lattype[3] = "FCC";
		lattype[4] = "HEX";
		lattyp = MyInput->GetChoice("lat", latname[1], "latticetype", lattype, 1);
	}

	if (freedomChoice != 1 && lattyp == 4){
		Message(warning,"Pinned, frozen and grafted ranges will be distorted when a hexagonal lattice is used. Choose standard, stencils or FCC instead to keep the correct shapes");
	}
    if (freedomChoice == 1) {
    	freedom = loose;
		MyInput->DontCombineParam(SEGMENT,name,"freedom","grafted_range");
		MyInput->DontCombineParam(SEGMENT,name,"freedom","frozen_range");
		MyInput->DontCombineParam(SEGMENT,name,"freedom","pinned_range");
    		LatRange = Lat->NewLatticeRange("free");
    } else if (freedomChoice == 2) {
    	freedom = pinned;
		MyInput->DontCombineParam(SEGMENT,name,"freedom","grafted_range");
		MyInput->DontCombineParam(SEGMENT,name,"freedom","frozen_range");
		if (!(MyInput->ValueSet(SEGMENT,name,"spot_frac"))){
			MyInput->AlwaysCombineParam(SEGMENT,name,"freedom","pinned_range");
		}

		Array<Text> spottype(1,2);
		spottype[1] = "sphere";
		spottype[2] = "node";

		SpotType = MyInput->GetChoice(SEGMENT,name,"spot_type", spottype, 2);

		if (SpotType==2) {
			SpotSize = MyInput->GetInt(SEGMENT,name,"spot_size",1,27,1);
			if (SpotSize == 1 || SpotSize == 2 ||SpotSize == 4 || SpotSize == 8 || SpotSize == 27) {
			} else {Message(fatal,"Only spotsizes 1, 2, 4, 8 and 27 supported" );}
		}
		if (SpotType==1) {
			Message(warning,"Spherical particles in development. Proceed when you are knowing what you are doing." );
			SpotSize = MyInput->GetInt(SEGMENT,name,"spot_size",1,50,5);
		}
		Text range;
		if (MyInput->ValueSet(SEGMENT,name,"pinned_range")) {
    		 range = MyInput->GetText(SEGMENT,name,"pinned_range");
		}
		else {
			Text name2 = MyInput->GetText(SEGMENT,name,"spot_mix");
			Message(literal,SEGMENT + " spot_mix " + name2 + ": SF_segment before range : " );
			
			range = MyInput->GetText(SEGMENT,name2,"pinned_range");
			Message(literal,SEGMENT + " spot_mix " + name2 + ": SF_segment after range : "+ range );
		}
		//Message(literal,SEGMENT + " spot_mix " + name + ": end changed code block 1: " );
		LatRange = Lat->NewLatticeRange(range);
		//Message(literal,SEGMENT + " spot_mix " + name + ": end changed code block 2: " );
		Text GNP;
		GNP= Blanks(100);
		GNP.Putreal(LatRange ->GetVolPos(),8);
		//Message(literal,SEGMENT + " spot_mix " + name + ": end changed code block 3: " );
		GNP= Copy(GNP.Strip().Frontstrip());
		if (LatRange ->GetNumPos()>1) Message(literal,SEGMENT + " : " + name + ": number of pinned sites : " + GNP);
		


	} else if (freedomChoice == 3) {
		freedom = grafted;
		MyInput->DontCombineParam(SEGMENT,name,"freedom","frozen_range");
		MyInput->DontCombineParam(SEGMENT,name,"freedom","pinned_range");
 		MyInput->AlwaysCombineParam(SEGMENT,name,"freedom","grafted_range");
    		Text range = MyInput->GetText(SEGMENT,name,"grafted_range");
		LatRange = Lat->NewLatticeRange(range);
	} else if (freedomChoice == 4) {
		freedom = frozen;
		MyInput->DontCombineParam(SEGMENT,name,"freedom","grafted_range");
		MyInput->DontCombineParam(SEGMENT,name,"freedom","pinned_range");
 		MyInput->AlwaysCombineParam(SEGMENT,name,"freedom","frozen_range");
    		Text range = MyInput->GetText(SEGMENT,name,"frozen_range");

 		Array<Text> spottype(1,2);
		spottype[1] = "sphere";
		spottype[2] = "node";
		SpotType = MyInput->GetChoice(SEGMENT,name,"spot_type", spottype, 2);

		if (SpotType==2) {
			SpotSize = MyInput->GetInt(SEGMENT,name,"spot_size",1,27,1);
			if (SpotSize == 1 || SpotSize == 2 ||SpotSize == 4 || SpotSize == 8 || SpotSize == 27) {
			} else {Message(fatal,"Only spotsizes 1, 2, 4, 8 and 27 supported" );}
		}
		if (SpotType==1) {
			Message(warning,"Spherical particles in development. Proceed when you are knowing what you are doing." );
			SpotSize = MyInput->GetInt(SEGMENT,name,"spot_size",1,50,5);
		}

		LatRange = Lat->NewLatticeRange(range);
		Text GNP;
		GNP= Blanks(100);
		GNP.Putreal(LatRange ->GetVolPos(),8);
		GNP= Copy(GNP.Strip().Frontstrip());
		if (LatRange ->GetNumPos()>0) Message(literal,SEGMENT + " : " + name + ": number of frozen sites : " + GNP);
	}
	double epsilon = 80;
	if (MyInput->ValueSet(SEGMENT,name,"epsilon")) {
		epsilon = MyInput->GetReal(SEGMENT,name,"epsilon",-DBL_MAX,DBL_MAX,80);
		if (epsilon < 1) {
			Message(warning,SEGMENT + " : " + name + " : epsilon is smaller than 1, are you sure!?");
		}
	}
	MC = MyInput->GetBoolean(SEGMENT,name,"MC",false);
	LD = MyInput->GetBoolean(SEGMENT,name,"LD",false);
	//Mayer_Saupe = MyInput->ValueSet(SEGMENT,name,"alpha_Mayer_Saupe");
	if (MC || LD) {
		if (freedom != pinned && freedom != frozen) {
			Message(fatal,MyInput,"Segment must be pinnned or frozen to use MC or LD");

		} else {
			Array<Text> moves(1,15);
			moveX=moveY=moveZ=false;
			moves[1] = "x";
			moves[2] = "y";
			moves[3] = "z";
			moves[4] = "xy";
			moves[5] = "yx";
			moves[6] = "xz";
			moves[7] = "zx";
			moves[8] = "yz";
			moves[9] = "zy";
			moves[10] = "xyz";
			moves[11] = "xzy";
			moves[12] = "zxy";
			moves[13] = "zyx";
			moves[14] = "yxz";
			moves[15] = "yzx";

   			int moveChoice = MyInput->GetChoice(SEGMENT,name,"move",moves,10);
			if (moveChoice==1) {moveX=true;}
			if (moveChoice==2) {moveY=true;}
			if (moveChoice==3) {moveZ=true;}
			if (moveChoice==4 || moveChoice==5) {moveX=moveY=true;}
			if (moveChoice==6 || moveChoice==7) {moveX=moveZ=true;}
			if (moveChoice==8 || moveChoice==9) {moveZ=moveY=true;}
			if (moveChoice>9) {moveX=moveY=moveZ=true;}
		}
	}
	valence = 0;
	if (MyInput->ValueSet(SEGMENT,name,"valence")) {
		valence = MyInput->GetReal(SEGMENT,name,"valence",-DBL_MAX,DBL_MAX,0);
	}
	phiBulk = 0;
    if (numStates == 0) {
    	numStates = 1;
    	StateQ.Dim(1,1);
    	StateQ[1] = new SF_State(name,valence,epsilon,0,0,false,1,0,freedom,LatRange,Lat);
    } else {
    	StateQ.Dim(1,numStates);
    	for (i=1; i<=numStates; i++) {
    		Text stateName = MyInput->GetText(SEGMENT,name,param[numChiParam+numPfParam+i]);
    		Array<Text> statePar(1,numChiParam+6);
    		for (int j=1; j<= numChiParam; j++) {
    			statePar[j] = "chi - " + chiParam[j];
    		}
    		statePar[numChiParam+1] = "valence";
    		statePar[numChiParam+2] = "epsilon";
    		statePar[numChiParam+3] = "phibulk";
    		statePar[numChiParam+4] = "alphabulk";
    		statePar[numChiParam+5] = "internal_degeneration";
    		statePar[numChiParam+6] = "internal_energy";
    		MyInput->CheckParameterNames(STATE,stateName,statePar);
    		double thisValence = valence;
    		double thisEpsilon = epsilon;
   			double thisPhiBulk = 0;
   			double thisAlphaBulk = 0;
   			Boolean thisInternalFreeEnergyGiven = false;
   			double thisInternalDegeneration = 0;
   			double thisInternalEnergy = 0;
			if (MyInput->ValueSet(STATE,stateName,"valence")) {
				thisValence = MyInput->GetReal(STATE,stateName,"valence",-DBL_MAX,DBL_MAX);
			}
			if (MyInput->ValueSet(STATE,stateName,"epsilon")) {
    			thisEpsilon = MyInput->GetReal(STATE,stateName,"epsilon",-DBL_MAX,DBL_MAX);
				if (thisEpsilon < 1) {
					Message(warning,SEGMENT + " : " + name + " : epsilon is smaller than 1, are you sure!?");
				}
    		}
    		if (MyInput->ValueSet(STATE,stateName,"phibulk")) {
    			thisPhiBulk = MyInput->GetReal(STATE,stateName,"phibulk",0,1);
    		}
    		if (MyInput->ValueSet(STATE,stateName,"alphabulk")) {
    			thisAlphaBulk = MyInput->GetReal(STATE,stateName,"alphabulk",0,1);
    		}
    		if (MyInput->ValueSet(STATE,stateName,"internal_degeneration")) {
    			thisInternalFreeEnergyGiven = true;
    			thisInternalDegeneration = MyInput->GetReal(STATE,stateName,"internal_degeneration",1,INT_MAX);
    			thisInternalEnergy = 0;
    		}
    		if (MyInput->ValueSet(STATE,stateName,"internal_energy")) {
    			thisInternalFreeEnergyGiven = true;
    			thisInternalEnergy = MyInput->GetReal(STATE,stateName,"internal_energy",-DBL_MAX,DBL_MAX);
    			if (!MyInput->ValueSet(STATE,stateName,"internal_degeneration")) {
    				thisInternalDegeneration = 1;
    			}
    		}
    		StateQ[i] = new SF_State(stateName,thisValence,thisEpsilon,thisPhiBulk,
    									thisAlphaBulk,thisInternalFreeEnergyGiven,thisInternalDegeneration,
    									thisInternalEnergy,freedom,LatRange,Lat);
    	}
    }
}
SF_Segment::~SF_Segment() {
	delete LatRange;
	for (int i=1; i<=numStates; i++) {
		delete StateQ[i];
	}
}

bool
SF_Segment::GetMC() const {
	return MC;
}

double
SF_Segment::GetValence() const {
	return valence;
}


int
SF_Segment::GetSpotSize() const {
	return SpotSize;
}

int
SF_Segment::GetSpotType() const {
	return SpotType;
}

bool
SF_Segment::GetLD() const {
	return LD;
}

//double
//SF_Segment::GetMayerSaupe() const {
//	if (Mayer_Saupe) return MyInput->GetReal(SEGMENT,name,"alpha_Mayer_Saupe",-100.0,100.0);
//	else return 0;
//}

bool
SF_Segment::MoveAllowed(int direction) const {
	switch (direction) {
		case 1: return moveX;
			break;
		case 2: return moveY;
			break;
		case 3: return moveZ;
			break;
		default:
			return false;
	}
}

Text
SF_Segment::GetName() const {
	return name;
}
void
SF_Segment::GetOutput(Output* Out) const {
	SF_State* State;
	if (numStates > 1) {
		Text stateId;
		Text number;
		for (int i=1; i<=numStates; i++) {
			number = Blanks(100);
			number.Putint(i);
			number = Copy(number.Frontstrip());
			stateId = "state" + number;
			Out->PutText(SEGMENT,name,stateId,(StateQ[i])->GetName());
		}
	}
	if (freedom == loose) {
		Out->PutText(SEGMENT,name,"freedom","free");
	} else if (freedom == pinned) {
		Out->PutText(SEGMENT,name,"freedom","pinned");
		Out->PutText(SEGMENT,name,"pinned_range",LatRange->GetOutput());
	} else if (freedom == grafted) {
		Out->PutText(SEGMENT,name,"freedom","grafted");
		Out->PutText(SEGMENT,name,"grafted_range",LatRange->GetOutput());
	} else if (freedom == frozen) {
		Out->PutText(SEGMENT,name,"freedom","frozen");
		Out->PutText(SEGMENT,name,"frozen_range",LatRange->GetOutput());
	}
	Out->PutReal(SEGMENT,name,"phibulk",phiBulk);
	int M = Lat->GetTotalNumLayers();
	double theta = 0;
	Vector phi = GetPhi();
	Lat->SubtractBoundaries(phi);
	Lat->MultiplyWithLatticeSites(phi);
	int z;
	for (z=1; z<=M; z++) {
		theta += phi[z];
	}
	Out->PutReal(SEGMENT,name,"theta",theta);
	double alpha_state;
	if (numStates > 1) {
		for (int i=1; i<=numStates; i++) {
			alpha_state = StateQ[i]->GetTheta()/theta;
			Out->PutReal(SEGMENT,name,"alpha-"+StateQ[i]->GetName(),alpha_state);
			}
	}

//	This is the wrong way to compute excess FL
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
	Out->PutReal(SEGMENT,name,"theta excess",theta);
	if (numStates == 1) {
		State = StateQ[1];
		Out->PutReal(SEGMENT,name,"valence",State->GetValence());
		Out->PutReal(SEGMENT,name,"epsilon",State->GetEpsilon());
		Out->PutReal(SEGMENT,name,"int.free energy",State->GetInternalFreeEnergy());
	}
	Out->PutProfile(SEGMENT,name,"phi",GetPhi());
	if (freedom != frozen) {
		Out->PutProfile(SEGMENT,name,"G",GetSWF());
	}
	if (numStates > 1) {
		for (int i=1; i<=numStates; i++) {
			(StateQ[i])->GetOutput(Out);
		}
	}
}
Vector
SF_Segment::GetPhi() const {
	int M = Lat->GetTotalNumLayers();
	int z,i,j;
	SF_MolSegment* MolSeg;
	Vector molPhi;
	Vector phi(1,M);
	for (i=1; i<=MolSegmentQ.Cardinal(); i++) {
		MolSeg = (SF_MolSegment*) MolSegmentQ[i];
		Array<DensityPart> DensityPartQ;
		DensityPartQ = MolSeg->GetDensityPartQ();
		int numDensPart = DensityPartQ.Upperbound() - DensityPartQ.Lowerbound() + 1;
		for (j=1; j<=numDensPart; j++) {
			molPhi = MolSeg->GetPhi(DensityPartQ[j]);
			for (z=1; z<=M; z++) {
				phi[z] += molPhi[z];
			}
		}
	}
	if (freedom == frozen) {
		for (z=1; z<=M; z++) {
			//if (LatRange->InRange(z)) phi[z] = 1;
			phi[z]=LatRange->GetRangeValue(z);
		}
		Lat->SetBoundaries(phi);
	}
	return phi;
}
double
SF_Segment::GetPhiBulk() const {
	return phiBulk;
}
void
SF_Segment::UpdatePhiBulk() {
	SF_MolSegment *MolSeg;
	phiBulk = 0;
	int i;
	for (i=1; i<=MolSegmentQ.Cardinal(); i++) {
		MolSeg = (SF_MolSegment*) MolSegmentQ[i];
		MolSeg->UpdatePhiBulk();
		phiBulk += MolSeg->GetPhiBulk();
	}
	for (i=1; i<=numStates; i++) {
		(StateQ[i])->UpdatePhiBulk();
	}
}
SegmentFreedom
SF_Segment::GetFreedom() const {
	return freedom;
}
LatticeRange*
SF_Segment::GetLatRange() const {
	return LatRange;
}
int
SF_Segment::GetNumStates(void) const {
	return numStates;
}
SF_State*
SF_Segment::GetState(int number) const {
	return StateQ[number];
}
Vector
SF_Segment::GetSWF() const {
	return swf;
}

void
SF_Segment::ClearAllPos() {
	LatRange->ClearAllPos();
	int i;
	SF_MolSegment* MolSegment;
	for (i=1; i<=MolSegmentQ.Cardinal(); i++) {
		MolSegment = (SF_MolSegment *) MolSegmentQ[i];
		MolSegment->ClearAllPos();
	}
	if (numStates == 1) {
		(StateQ[1])->ClearAllPos();
	} else {
		SF_State* State;
		for (i=1; i<=numStates; i++) {
			State = StateQ[i];
			State->ClearAllPos();
		}
	}
}

void
SF_Segment::DelPos(int r) {
	if (LatRange->InRange(r)) LatRange->SetPos(r,SpotType,SpotSize,0.0);
	int i;
	SF_MolSegment* MolSegment;
	for (i=1; i<=MolSegmentQ.Cardinal(); i++) {
		MolSegment = (SF_MolSegment *) MolSegmentQ[i];
		MolSegment->DelPos(r,SpotType,SpotSize);
	}
	if (numStates == 1) {
		(StateQ[1])->DelPos(r,SpotType,SpotSize);
	} else {
		SF_State* State;
		for (i=1; i<=numStates; i++) {
			State = StateQ[i];
			State->DelPos(r,SpotType,SpotSize);
		}
	}
}

void
SF_Segment::UpdatePos(double x, double y, double z, double* submask) {
	LatRange->SetPos(x,y,z,submask);
	int i;
	SF_MolSegment* MolSegment;
	for (i=1; i<=MolSegmentQ.Cardinal(); i++) {
		MolSegment = (SF_MolSegment *) MolSegmentQ[i];
		MolSegment->UpdatePos(x,y,z,submask);
	}
	if (numStates == 1) {
		(StateQ[1])->UpdatePos(x,y,z,submask);
	} else {
		SF_State* State;
		for (i=1; i<=numStates; i++) {
			State = StateQ[i];
			State->UpdatePos(x,y,z,submask);
		}
	}
}


void
SF_Segment::UpdatePos(int r) {
	double waarde = 1.0;
	UpdatePos(r,waarde);
}
void
SF_Segment::UpdatePos(int r,double waarde) {
	LatRange->SetPos(r,SpotType,SpotSize,waarde);

	int i;
	SF_MolSegment* MolSegment;
	for (i=1; i<=MolSegmentQ.Cardinal(); i++) {
		MolSegment = (SF_MolSegment *) MolSegmentQ[i];
		MolSegment->UpdatePos(r,SpotType,SpotSize);
	}
	if (numStates == 1) {
		(StateQ[1])->UpdatePos(r,SpotType,SpotSize);
	} else {
		SF_State* State;
		for (i=1; i<=numStates; i++) {
			State = StateQ[i];
			State->UpdatePos(r,SpotType,SpotSize);
		}
	}
}



void
SF_Segment::ResetRhoSet(){
	SF_State* State;
	for (int i=1; i<=numStates; i++) {
		State = StateQ[i];
		State->ResetRhoSet();
	}
}

void
SF_Segment::UpdateSWF() {
	int i;
	SF_MolSegment* MolSegment;
	for (i=1; i<=MolSegmentQ.Cardinal(); i++) {
		MolSegment = (SF_MolSegment *) MolSegmentQ[i];
		MolSegment->UpdateSWF();
	}
	if (numStates == 1) {
		(StateQ[1])->UpdateSWF();
		swf = (StateQ[1])->GetSWF();
	} else {
		int M = Lat->GetTotalNumLayers();
		SF_State* State;
		int z;
		Vector swfState;
		double alphaBulk;
		swf.Dim(1,M);
		for (i=1; i<=numStates; i++) {
			State = StateQ[i];
			State->UpdateSWF();
			swfState = State->GetSWF();
			Lat->SetBoundaries(swfState);
			alphaBulk = State->GetAlphaBulk();
			for (z=1; z<=M; z++) {
				swf[z] += alphaBulk*swfState[z];
			}
		}
	}
}
void
SF_Segment::UpdatePhiStates() {
	SF_MolSegment* MolSegment;
	for (int i=1; i<=MolSegmentQ.Cardinal(); i++) {
		MolSegment = (SF_MolSegment *) MolSegmentQ[i];
		MolSegment->UpdatePhiStates();
	}
}
int
SF_Segment::GetNumMolSegments() const {
	return MolSegmentQ.Cardinal();
}
SF_MolSegment*
SF_Segment::GetMolSegment(int number) const {
	return (SF_MolSegment*) MolSegmentQ[number];
}
SF_MolSegment*
SF_Segment::GetMolSegment(Text molName) const {
	SF_MolSegment* MolSeg;
	for (int i=1; i<=MolSegmentQ.Cardinal(); i++) {
		MolSeg = (SF_MolSegment*) MolSegmentQ[i];
		if (*Copy(MolSeg->GetMolName()) == *molName) {
			return MolSeg;
		}
	}
	return NULL;
}
SF_MolSegment*
SF_Segment::NewMolSegment(Text molName) {
	SF_MolSegment* MolSegment;
	int numMolSeg = MolSegmentQ.Cardinal();
	if (freedom == frozen) {
		Message(fatal,MyInput,"Cannot use segment '" + name +
			"' for 'mol : " + molName +
			" : composition', \nsegment has restriction 'frozen'.");
	}
	int i;
	for (i=1; i<=numMolSeg; i++) {
		MolSegment = GetMolSegment(i);
		if (*Copy(MolSegment->GetMolName()) == *molName) return MolSegment;
	}
	MolSegment = new SF_MolSegment(name,molName,Lat,MyInput);
	MolSegment->Into(MolSegmentQ);
	SF_State* State;
	SF_MolState* MolState;
	for (i=1; i<=numStates; i++) {
		State = StateQ[i];
		MolState = MolSegment->GetState(i);
		State->AddMolState(MolState);
	}
    return MolSegment;
}



