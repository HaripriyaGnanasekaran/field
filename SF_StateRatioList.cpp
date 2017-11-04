#include <iostream>
#include "SF_StateRatioList.h"

SF_StateRatioList::SF_StateRatioList(Array<SF_StateRatio*> StateRatioQ_,
									 SF_SegmentList* SegQ_) {
	StateRatioQ = StateRatioQ_;
	numStateRatios = StateRatioQ.Upperbound();
	SegQ = SegQ_;
	SF_State* State = (StateRatioQ[1])->GetState1();
	Segment = SegQ->GetBaseSegment(State);
	for (int i=1; i<=numStateRatios; i++) {
		State = (StateRatioQ[1])->GetState1();
		if (SegQ->GetBaseSegment(State) != Segment) {
			Message(fatal,"Error defining stateratiolist");
		}
		State = (StateRatioQ[1])->GetState2();
		if (SegQ->GetBaseSegment(State) != Segment) {
			Message(fatal,"Error defining stateratiolist");
		}
	}
}
SF_StateRatioList::~SF_StateRatioList() {
}
void
SF_StateRatioList::UpdateAlphaBulk() {
	double alphaBulk;
	SF_State* State;
	int i;
	int numStates = Segment->GetNumStates();
	for (i=1; i<=numStates; i++) {
		State = Segment->GetState(i);
		State->SetAlphaBulk(0);
	}
	SF_StateRatio* Ratio;
	Ratio = StateRatioQ[1];
	State = Ratio->GetState1();
	State->SetAlphaBulk(1);
	State = Ratio->GetState2();
	State->SetAlphaBulk(1/Ratio->GetAlphaBulkRatio());
	SF_State* State1;
	SF_State* State2;
	for (i=2; i<=numStateRatios; i++) {
		Ratio = StateRatioQ[i];
		State1 = Ratio->GetState1();
		State2 = Ratio->GetState2();
		if (State1->GetAlphaBulk() > 0) {
			if (State2->GetAlphaBulk() > 0) {
				cout << State2->GetName().MainC() <<  '\t' << State2->GetAlphaBulk() << endl;
				Message(fatal,"error in SF_StateRatioList::UpdateAlphaBulk()");
			} else {
//				fedisableexcept(FE_OVERFLOW|FE_INVALID);
				alphaBulk = State1->GetAlphaBulk()/
Ratio->GetAlphaBulkRatio();
				State2->SetAlphaBulk(alphaBulk);
//				feenableexcept(FE_OVERFLOW|FE_INVALID);
			}
		} else {
			if (State2->GetAlphaBulk() <= 0) {
				cout << State2->GetName().MainC() <<  '\t' << State2->GetAlphaBulk() << endl;
				Message(fatal,"error in SF_StateRatioList::UpdateAlphaBulk()");
			} else {
				//yansen nov 2003
				State1->SetAlphaBulk(Ratio->GetAlphaBulkRatio()*State2->GetAlphaBulk());
			}
		}
	}
	double sumAlpha = 0;
	for (i=1; i<=Segment->GetNumStates(); i++) {
		sumAlpha += (Segment->GetState(i))->GetAlphaBulk();
	}
	for (i=1; i<=Segment->GetNumStates(); i++) {
		State = Segment->GetState(i);
		State->SetAlphaBulk(State->GetAlphaBulk()/sumAlpha);
	}
}
