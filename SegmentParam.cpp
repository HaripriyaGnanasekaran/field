#include "SegmentParam.h"

SegmentParam::SegmentParam() {
}
SegmentParam::SegmentParam(Text name_, 
						   Input* MyInput, 
						   Text ID) {
	name = name_;
    Text stateName;
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
	if (numStates != 0) States.Dim(1,numStates);
	for (i=1; i<=numStates; i++) {
		number = Blanks(1000);
		number.Putint(i);
		number = Copy(number.Frontstrip());
		stateId = "state" + number;
		stateName = MyInput->GetText(ID,name,stateId);
		States[i] = new SegmentParam(stateName,MyInput,STATE);
	}
	Array<Text> paramTotal = MyInput->GetParameters(ID,name);
	int numParamTotal = paramTotal.Upperbound();
	numChiParam = 0;
	for (i=1; i<=numParamTotal; i++) {
		if (ChiParameter(paramTotal[i])) numChiParam++;
	}
	int j=0;
	if (numChiParam != 0) {
		chiParam.Dim(1,numChiParam);
		for (i=1; i<=numParamTotal; i++) {
			if (ChiParameter(paramTotal[i])) chiParam[++j] = StripChi(paramTotal[i]);
		}
	}

	numParam = numParamTotal-numChiParam;
	j=0;
	if (numParam != 0) {
		param.Dim(1,numParam);
		for (i=1; i<=numParamTotal; i++) {
			if (!ChiParameter(paramTotal[i])) param[++j] = Copy(paramTotal[i]);
		}
	}
	for (i=1; i<=numChiParam; i++) {
		Text test = chiParam[i];
		for (j=1; j<=NumStates(); j++) {
			if (*Copy(test) == *Copy(StateName(j))) {
				Message(fatal,MyInput, "Parameter '" +
					SEGMENT + " : " + name + " : chi - " + StateName(j) +
					"' cannot be set.\n'" + StateName(j) + 
					"' is a state of segment '" + name + "'.");
			}
		}
		if (*test == *name) {
			Message(fatal,MyInput, "Parameter '" +
				ID + " : " + name + " : chi - " + test + "' cannot be set.\n");
		}
	}

	SegmentParam* State;
	for (i=1; i<=NumStates(); i++) {
		State = States[i];
		for (j=1; j<=State->NumChiParam(); j++) {
			if (*name == *Copy(State->ChiParam(j))) {
				Message(fatal,MyInput, "Parameter '" +
					STATE + " : " + State->Name() + " : chi - " + name +
					"' cannot be set.\n'" + State->Name() +
					"' is a state of segment '" + name + "'.");
			}
			if (*Copy(State->Name()) == *Copy(State->ChiParam(j))) {
				Message(fatal,MyInput, "Parameter '" +  
					STATE + " : " + State->Name() + " : chi - " + 
					State->ChiParam(j) + "' cannot be set.\n");
			}
			if (ChiParamSet(State->ChiParam(j))) {
				Message(fatal,MyInput, "Parameter 'chi - " + State->ChiParam(j) + 
					"' is defined in segment '" + SEGMENT + " : "  + name + 
					"'\nand one of its states, define either in segment or in states."); 
			}
		}
		for (j=1; j<=State->NumParam(); j++) {
			if (ParamSet(State->Param(j))) {
				Message(fatal,MyInput, "Parameter '" + State->Param(j) + 
					"' is defined in segment '" + SEGMENT + " : "  + name + 
					"'\nand one of its states, define either in segment or in states."); 
			}
		}
			
	}

}
SegmentParam::~SegmentParam() {
	int i;
	for (i=1; i<=numStates; i++) {
		delete States[i];
	}
}
Text
SegmentParam::Name() const {
	return name;
}
int
SegmentParam::NumParam() const {
	return numParam;
}
Text
SegmentParam::Param(int number) const {
	return param[number];
}
int
SegmentParam::NumStates() const {
	return numStates;
}
Text
SegmentParam::StateName(int number) const {
	if (numStates == 0) return name;
	SegmentParam* State;
	State = States[number];
	return State->Name();
}
SegmentParam*
SegmentParam::StatePTR(int number) const {
	SegmentParam* State;
	State = States[number];
	return State;
}
int
SegmentParam::NumChiParam() const {
	return numChiParam;
}

Text
SegmentParam::ChiParam(int number) const {
	return chiParam[number];
}

Boolean
SegmentParam::ParamSet(Text test) const {
	int i;
	for (i=1; i<=numParam; i++) {
		if (*Copy(param[i]) == *test) return true;
	}
	return false;
}
Boolean
SegmentParam::ChiParamSet(Text test) const {
	int i;
	for (i=1; i<=numChiParam; i++) {
		if (*Copy(chiParam[i]) == *test) return true;
	}
	return false;
}

Boolean
SegmentParam::ChiParameter(Text test) const {
	if (test.Length() < 6) return false;
	Text t;
	t = test.Sub(1,6);
	if (*t == *Copy("chi - ")) return true;
	else return false;
}
Text
SegmentParam::StripChi(Text test) const {
	int length = test.Length();
	test = test.Sub(7,length - 6);
	return test;
}

