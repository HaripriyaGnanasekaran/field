#include "SF_StateRatio.h"

SF_StateRatio::SF_StateRatio(SF_State* State1_, SF_State* State2_) {
	State1 = State1_;
	State2 = State2_;
	if (State1 == State2) {
		Sysout().Outtext("Error defining state ratios, equal states on both sides of equation");
		Sysout().Outimage();
		exit(-1);
	}
	ComputeAlphaBulkRatio();
}
SF_StateRatio::~SF_StateRatio() {
}
SF_State*
SF_StateRatio::GetState1() const {
	return State1;
}
SF_State*
SF_StateRatio::GetState2() const {
	return State2;
}
double 
SF_StateRatio::GetAlphaBulkRatio(void) const {
	return value;
}
void 
SF_StateRatio::SetAlphaBulkRatio(double value_) {
	value = value_;
}
Boolean
SF_StateRatio::operator==(const SF_StateRatio& That) const {
	if (That.GetState1() == State1 && That.GetState2() == State2) {
		return true;
	}
	else return false;
}
Boolean
SF_StateRatio::SameStates(const SF_StateRatio& That) const {
	if (That.GetState1() == State1) {
		if (That.GetState2() == State2) {
			return true;
		} else {
			return false;
		}
	} else if (That.GetState1() == State2) {
		if (That.GetState2() == State1) {
			return true;
		} else {
			return false;
		}
	} else {
		return false;
	}
}
void
SF_StateRatio::ComputeAlphaBulkRatio() {
	value = State1->GetAlphaBulk()/State2->GetAlphaBulk();
}
