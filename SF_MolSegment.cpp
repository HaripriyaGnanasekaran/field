#include "SF_MolSegment.h"

SF_MolSegment::SF_MolSegment(Text name_,
							 Text molName_,
							 Lattice* Lat_,
							 Input* MyInput) {
	name = name_;
	molName = molName_;
	Lat = Lat_;
	MyInput->SetDefaultWarningsOff();
	numStates = 0;
    Text stateId;
    Text number;
	stateId = "state1";
	//Mayer_Saupe = MyInput->ValueSet(SEGMENT,name,"alpha_Mayer_Saupe");
	int i=1;
	while (MyInput->ValueSet(SEGMENT,name,stateId)) {
		number = Blanks(100);
		number.Putint(++i);
		number = Copy(number.Frontstrip());
		stateId = "state" + number;
	}
	numStates = i-1;
	Array<Text> freedoms(1,4);
	freedoms[1] = "free";
	freedoms[2] = "pinned";
	freedoms[3] = "grafted";
	freedoms[4] = "frozen";
    int freedomChoice = MyInput->GetChoice(SEGMENT,name,"freedom",freedoms,1);
    if (freedomChoice == 1) {
    	freedom = loose;
    	LatRange = Lat->NewLatticeRange("free");
    } else if (freedomChoice == 2) {
		MC = MyInput->GetBoolean(SEGMENT,name,"MC",false);
		LD = MyInput->GetBoolean(SEGMENT,name,"LD",false);
		Text range;
    	freedom = pinned;
    	if (MyInput->ValueSet(SEGMENT,name,"pinned_range")) {
    		range = MyInput->GetText(SEGMENT,name,"pinned_range");
		}
    	else {
			Message(literal,SEGMENT + " spot_mix " + name + ": SF_molsegment before newcode : " );
			Text name2 = MyInput->GetText(SEGMENT,name,"spot_mix");
			Message(literal,SEGMENT + " spot_mix " + name2 + ": SF_molsegment before range : " );
			range = MyInput->GetText(SEGMENT,name2,"pinned_range");
			Message(literal,SEGMENT + " spot_mix " + name2 + ": SF_molsegment after range : " );
		}
		LatRange = Lat->NewLatticeRange(range);


	} else if (freedomChoice == 3) {
    	freedom = grafted;
    	Text range = MyInput->GetText(SEGMENT,name,"grafted_range");
		LatRange = Lat->NewLatticeRange(range);
		Boolean found = false;
		int M = Lat->GetTotalNumLayers();
		for (int z=1; z<=M; z++) {
			if (LatRange->InRange(z)) {
				if (found) {
					Message(fatal,"Grafted segments are only supported for a range of 1 layer. "
					"You probably want to give them freedom : pinned");
				}
				found = true;
			}
		}
	} else if (freedomChoice == 4) {
		Message(fatal,MyInput,"Cannot create MolSegment '" +
			name + "' with freedom 'frozen'");
	}
	double valence = MyInput->GetReal(SEGMENT,name,"valence",-DBL_MAX,DBL_MAX,0);
    double epsilon = MyInput->GetReal(SEGMENT,name,"epsilon",-DBL_MAX,DBL_MAX,80);
    if (numStates == 0) {
    	numStates = 1;
    	StateQ.Dim(1,1);
    	StateQ[1] = new SF_MolState(name,molName,valence,epsilon,1.0,freedom,LatRange,Lat);
    } else {
    	int i;
    	StateQ.Dim(1,numStates);
    	for (i=1; i<=numStates; i++) {
    		number = Blanks(100);
			number.Putint(i);
			number = Copy(number.Frontstrip());
    		Text stateParam = "state" + number;
    		Text stateName = MyInput->GetText(SEGMENT,name,stateParam);
    		double stateValence = valence;
    		double stateEpsilon = epsilon;
 			if (MyInput->ValueSet(STATE,stateName,"valence")) {
				stateValence = MyInput->GetReal(STATE,stateName,"valence",-DBL_MAX,DBL_MAX,0);
 			}
    		if (MyInput->ValueSet(STATE,stateName,"epsilon")) {
    			stateEpsilon = MyInput->GetReal(STATE,stateName,"epsilon",-DBL_MAX,DBL_MAX,80);
    		}
    		StateQ[i] = new SF_MolState(stateName,molName,stateValence,
    				stateEpsilon,1.0,freedom,LatRange,Lat);
    	}
   }
	MyInput->SetDefaultWarningsOn();
	numPhi=0;
	phiBulk = 0;
	phiRef = 0;
}
SF_MolSegment::~SF_MolSegment() {
	Out();
	int i;
	for (i=1; i<=numStates; i++) {
		delete StateQ[i];
	}
	delete LatRange;
}
Text
SF_MolSegment::GetName() const {
	return name;
}
Text
SF_MolSegment::GetMolName() const {
	return molName;
}
void
SF_MolSegment::GetOutput(Output* Out) const {
	int i;
	int M = Lat->GetTotalNumLayers();
	//int gradients=Lat->GetNumGradients();
	for (i=1; i<=numPhi; i++) {
		Boolean outputAvZ = false;
		LatticeRange* LayerZero=NULL;
		if (Lat->GetNumGradients() == 1) {
			outputAvZ = true;
			LayerZero = new LatticeRange1D(1,1,M);
		}

		switch (DensityPartQ[i]) {
			case total:
				if (outputAvZ) {
					double avZ = Lat->MomentUnweighted(PhiQ[i],1,LayerZero,-0.5);
					avZ /= Lat->MomentUnweighted(PhiQ[i],0,LayerZero,-0.5);
					Out->PutReal("mol",molName,"av z phi-"+name,avZ);
					Vector phiExc = GetPhi(total);
					int z;
					for (z=1; z<=M; z++) {
						phiExc[z] -= phiBulk;
					}
					LatticeRange* LayerZero = new LatticeRange1D(1,1,M);
					Lat->MultiplyWithLatticeSites(phiExc);
					double R = Lat->MomentUnweighted(phiExc,1,LayerZero,-0.5)/Lat->MomentUnweighted(phiExc,0,LayerZero,-0.5);
					double R2 = Lat->MomentUnweighted(phiExc,2,LayerZero,-0.5)/Lat->MomentUnweighted(phiExc,0,LayerZero,-0.5);
					Out->PutReal("mol",molName,"1st moment exc. n-"+name,R);
					Out->PutReal("mol",molName,"2nd moment exc. n-"+name,pow(R2,0.5));

					Lat->DivideByLatticeSites(phiExc);
					R = Lat->MomentUnweighted(phiExc,1,LayerZero,-0.5)/Lat->MomentUnweighted(phiExc,0,LayerZero,-0.5);
					R2 = Lat->MomentUnweighted(phiExc,2,LayerZero,-0.5)/Lat->MomentUnweighted(phiExc,0,LayerZero,-0.5);
					delete LayerZero;
					//Out->PutReal("mol",molName,"1st moment exc. n-"+name,R);
					for (z=1; z<=M; z++) {
						phiExc[z] += phiBulk;
					}
					Out->PutReal("mol",molName,"1st moment exc. phi-"+name,R);
					Out->PutReal("mol",molName,"2nd moment exc. phi-"+name,pow(R2,0.5));
				}

				Out->PutProfile("mol",molName,"phi-"+name,PhiQ[i]);

				
				//if (Lat->GetMatrixApproach()==secondOrder) {
				//	if (gradients==1) {
				//		Out->PutProfile("mol",molName,"phi_"+name+"_x",GetPhi(x_dir));
				//		Out->PutProfile("mol",molName,"phi_"+name+"_yz",GetPhi(yz_dir));
				//	} else {
				//		Out->PutProfile("mol",molName,"phi_"+name+"_x",GetPhi(x_dir));
				//		Out->PutProfile("mol",molName,"phi_"+name+"_y",GetPhi(y_dir));
				//		Out->PutProfile("mol",molName,"phi_"+name+"_z",GetPhi(z_dir));
				//	}
				//}
				break;
			case unconstrained:
				if (outputAvZ) {
					double avZ = Lat->MomentUnweighted(PhiQ[i],1,LayerZero,-0.5);
					avZ /= Lat->MomentUnweighted(PhiQ[i],0,LayerZero,-0.5);
					Out->PutReal("mol",molName,"av z phi-"+name+" unconstrained",avZ);
				}
				Out->PutProfile("mol",molName,"phi-"+name+" unconstrained",PhiQ[i]);
				break;
			case constrained:
				if (outputAvZ) {
					double avZ = Lat->MomentUnweighted(PhiQ[i],1,LayerZero,-0.5);
					avZ /= Lat->MomentUnweighted(PhiQ[i],0,LayerZero,-0.5);
					Out->PutReal("mol",molName,"av z phi-"+name+" constrained",avZ);
				}
				Out->PutProfile("mol",molName,"phi-"+name+" constrained",PhiQ[i]);
				break;
			case renorm:
				if (outputAvZ) {
					double avZ = Lat->MomentUnweighted(PhiQ[i],1,LayerZero,-0.5);
					avZ /= Lat->MomentUnweighted(PhiQ[i],0,LayerZero,-0.5);
					Out->PutReal("mol",molName,"av z phi-"+name+" renorm",avZ);
				}
				Out->PutProfile("mol",molName,"phi-"+name+" renorm",PhiQ[i]);
				break;
			default:
				Message(implementation,"Unknown densitypart "
					"requested in SF_Molecule::GetOutput");
				break;
		}
		if (Lat->GetNumGradients() == 1) {
			delete LayerZero;
		}
	}
	if (numStates>1) {
		for (i=1; i<=numStates; i++) {
			Out->PutReal("mol",molName,"alphaAv-"+(StateQ[i])->GetName(),GetAlphaAv(i));
			Out->PutProfile("mol",molName,"alpha-"+(StateQ[i])->GetName(),GetAlpha(i));
		}
	}
}
double
SF_MolSegment::GetPhiBulk() const {
	return phiBulk;
}
void
SF_MolSegment::SetPhiBulk(double phiBulk_) {
	phiBulk = phiBulk_;
	int i;
	SF_MolState* State;
	for (i=1; i<=numStates; i++) {
		State = StateQ[i];
		State->SetPhiBulk(phiBulk);
	}
}
Vector
SF_MolSegment::GetPhiBulkBoundaries() const {
	return phiBulkBoundaries;
}
void
SF_MolSegment::SetPhiBulkBoundaries(Vector phiBulkBoundaries_) {
	phiBulkBoundaries = phiBulkBoundaries_;
	int i;
	SF_MolState* State;
	for (i=1; i<=numStates; i++) {
		State = StateQ[i];
		State->SetPhiBulkBoundaries(phiBulkBoundaries);
	}
}
double
SF_MolSegment::GetPhiRef() const {
	return phiRef;
}
void
SF_MolSegment::SetPhiRef(double phiRef_) {
	phiRef = phiRef_;
	(StateQ[1])->SetPhiRef(phiRef);
}
SegmentFreedom
SF_MolSegment::GetFreedom() const {
	return freedom;
}

Boolean
SF_MolSegment::GetMC() const {
	return MC;
}

//Boolean
//SF_MolSegment::MayerSaupeSet() const {
//	//cout<<"MSSet in SF_MolSegment"<<endl;
//	return Mayer_Saupe;
//}

Boolean
SF_MolSegment::GetLD() const {
	return LD;
}

double
SF_MolSegment::GetPhiGrafted() const {
	if (freedom != grafted) {
		Message(fatal,"Programming error in call to SF_MolSegment::GetPhiGrafted()");
	}
	int M = Lat->GetTotalNumLayers();
	double value = 0;
	Vector phi = GetPhi(total);
	for (int z=1; z<=M; z++) {
		if (LatRange->InRange(z)) {
			return phi[z];
		}
	}
	// never get here
	Message(fatal,"Programming error in use of SF_MolSegment::GetPhiGrafted()");
	return value;
}
LatticeRange*
SF_MolSegment::GetLatRange() const {
	return LatRange;
}
int
SF_MolSegment::GetNumStates() const {
	return numStates;
}
SF_MolState*
SF_MolSegment::GetState(int number) const {
	return StateQ[number];
}
Vector
SF_MolSegment::GetAlpha(int numState) const {
	int M = Lat->GetTotalNumLayers();
	Vector alpha(1,M);
	SF_MolState* State = StateQ[numState];
	Vector swfState = State->GetSWF();
	for (int z=1; z<=M; z++) {
		if (swf[z]>0) {
			alpha[z] = State->GetAlphaBulk()*swfState[z]/swf[z];
		} else {
			alpha[z] = 0;
		}
	}
	return alpha;
}
double
SF_MolSegment::GetAlphaAv(int numState) const {
	int M = Lat->GetTotalNumLayers();
	int layerCount = 0;
	double alphaAv = 0;
	double alphabulk;
	SF_MolState* State = StateQ[numState];
	Vector swfState = State->GetSWF();
	alphabulk = State->GetAlphaBulk();
	for (int z=1; z<=M; z++) {
		// error: on boundaries: swfState can be > 0 while swf == 0
		// a better measure for AlphaAv is thetaState/theta
		if (swf[z]>0) {
			alphaAv += alphabulk*swfState[z]/swf[z];
			layerCount++;
		}
	}
	return alphaAv/layerCount;
}
Vector
SF_MolSegment::GetSWF() const {
	return swf;
}



void
SF_MolSegment::ClearAllPos() {
	LatRange->ClearAllPos();

	if (numStates == 1) {
		(StateQ[1])->ClearAllPos();
	} else {
		SF_MolState* State;
		for (int i=1; i<=numStates; i++) {
			State = StateQ[i];
			State->ClearAllPos();
		}
	}
}


void
SF_MolSegment::DelPos(int r,int SpotType,int SpotSize) {
	if (LatRange->InRange(r)) LatRange->SetPos(r,SpotType,SpotSize,0.0);

	if (numStates == 1) {
		(StateQ[1])->DelPos(r,SpotType,SpotSize);
	} else {
		SF_MolState* State;
		for (int i=1; i<=numStates; i++) {
			State = StateQ[i];
			State->DelPos(r,SpotType,SpotSize);
		}
	}
}

void
SF_MolSegment::UpdatePos(double x, double y, double z, double *submask) {
	LatRange->SetPos(x,y,z,submask);
	if (numStates == 1) {
		(StateQ[1])->UpdatePos(x,y,z,submask);
	} else {
		SF_MolState* State;
		for (int i=1; i<=numStates; i++) {
			State = StateQ[i];
			State->UpdatePos(x,y,z,submask);
		}
	}
}

void
SF_MolSegment::UpdatePos(int r,int SpotType, int SpotSize) {
	LatRange->SetPos(r,SpotType,SpotSize,1.0);
	if (numStates == 1) {
		(StateQ[1])->UpdatePos(r,SpotType,SpotSize);
	} else {
		SF_MolState* State;
		for (int i=1; i<=numStates; i++) {
			State = StateQ[i];
			State->UpdatePos(r,SpotType,SpotSize);
		}
	}
}

void
SF_MolSegment::UpdateSWF() {
	int i;
	if (numStates == 1) {
		(StateQ[1])->UpdateSWF();
		swf = (StateQ[1])->GetSWF();
	} else {
		int M = Lat->GetTotalNumLayers();
		SF_MolState* State;
		int z;
		Vector swfState;
		double alphaBulk;
		swf.Dim(1,M);
		for (i=1; i<=numStates; i++) {
			State = StateQ[i];
			State->UpdateSWF();
			swfState = State->GetSWF();
			alphaBulk = State->GetAlphaBulk();
			for (z=1; z<=M; z++) {
				swf[z] += alphaBulk*swfState[z];
			}
		}
	}
}
void
SF_MolSegment::UpdatePhiBulk() {
	for (int i=1; i<=numStates; i++) {
		(StateQ[i])->SetPhiBulk(phiBulk);
	}
}
void
SF_MolSegment::UpdatePhiStates() {
	SF_MolState* State;
	int M = Lat->GetTotalNumLayers();
	for (int i=1; i<=numStates; i++) {
		State = StateQ[i];
		Vector swfState = State->GetSWF();
		for (int j=1; j<=numPhi; j++) {
			Vector phiState = State->GetPhi(DensityPartQ[j]);
			Vector phiSeg = GetPhi(DensityPartQ[j]);
			for (int z=1; z<=M; z++) {
				if (swf[z]>0) {
					phiState[z] = State->GetAlphaBulk()*swfState[z]*phiSeg[z]/swf[z];
				} else {
					phiState[z] = 0;
				}
			}
		}
	}
}
Vector
SF_MolSegment::GetPhi(DensityPart densityPart) const {
	int i;
	for (i=1; i<=numPhi; i++) {
		if (DensityPartQ[i] == densityPart) {
			//Lat->SubtractBoundaries(PhiQ[i]);
			//Lat->RestoreBoundaries(PhiQ[i]);
			return PhiQ[i];
		}
	}
	if (densityPart == total) {
		int M = Lat->GetTotalNumLayers();
		Vector phiTotal(1,M);
		Vector phi;
		for (i=1; i<=numPhi; i++) {
			phi = PhiQ[i];
			for (int z=1; z<=M; z++) {
				phiTotal[z] += phi[z];
			}
		}
		return phiTotal;
	}
	Message(debug,"Programming error in call to "
		"SF_MolSegment::GetPhi(DensityPart densityPart), "
		"Create phi before using it");
	Vector a;
	return a;
}
// allocate the Phi[z] vectors, descends to SF_MolState
void
SF_MolSegment::CreatePhi(DensityPart densityPart) {
	int i;
	SF_MolState* State;
	for (i=1; i<=numStates; i++) {
		State = StateQ[i];
		State->CreatePhi(densityPart);
	}
	if (numPhi == 0) {
		PhiQ.Dim(1,1);
		DensityPartQ.Dim(1,1);
		Vector phi;
		if (numStates > 1) phi.Dim(1,Lat->GetTotalNumLayers());
		else {
			State = StateQ[1];
			phi = State->GetPhi(densityPart);
		}
		PhiQ[1] = phi;
		DensityPartQ[1] = densityPart;
		numPhi = 1;
	} else {
		Array<Vector> PhiQtemp(1,numPhi+1);
		Array<DensityPart> DensityPartQtemp(1,numPhi+1);
		for (i=1; i<=numPhi; i++) {
			PhiQtemp[i] = PhiQ[i];
			DensityPartQtemp[i] = DensityPartQ[i];
		}
		Vector phi;
		if (numStates > 1) phi.Dim(1,Lat->GetTotalNumLayers());
		else {
			State = StateQ[1];
			phi = State->GetPhi(densityPart);
		}
		numPhi++;
		PhiQtemp[numPhi] = phi;
		DensityPartQtemp[i] = densityPart;
		PhiQ = PhiQtemp;
		DensityPartQ = DensityPartQtemp;
	}
}
void
SF_MolSegment::DeletePhi(DensityPart densityPart) {
	int i;
	int number = 0;
	for (i=1; i<=numPhi; i++) {
		if (DensityPartQ[i] == densityPart) number = i;
	}
	if (number == 0) Message(fatal,"Programming error in call to "
		"SF_MolSegment::DeletePhi(DensityPart densityPart), "
		"Create phi before deleting it");
	Array<Vector> PhiQtemp(1,numPhi-1);
	Array<DensityPart> DensityPartQtemp(1,numPhi-1);
	for (i=1; i<=numPhi; i++) {
		if (i < number) {
			PhiQtemp[i] = PhiQ[i];
			DensityPartQtemp[i] = DensityPartQ[i];
		} else if (i > number) {
			PhiQtemp[i-1] = PhiQ[i];
			DensityPartQtemp[i-1] = DensityPartQ[i];
		}
	}
	numPhi--;
	PhiQ = PhiQtemp;
	DensityPartQ = DensityPartQtemp;
	SF_MolState* State;
	for (i=1; i<=numStates; i++) {
		State = StateQ[i];
		State->DeletePhi(densityPart);
	}
}
Array<DensityPart>
SF_MolSegment::GetDensityPartQ() const {
	return DensityPartQ;
}





