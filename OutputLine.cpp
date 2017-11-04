#include "OutputLine.h"

OutputLine::OutputLine(const Text Element1_,
					   const Text Element2_,
					   const Text Element3_,
					   const Text Value_) {
	Element1 = Element1_;
	Element2 = Element2_;
	Element3 = Element3_;
	Value = Value_;
	profile = false;
	vector = false;
}
OutputLine::OutputLine(const Text Element1_,
					   const Text Element2_,
					   const Text Element3_,
					   const Text Value_,
					   const Vector profileValue_) {
	Element1 = Element1_;
	Element2 = Element2_;
	Element3 = Element3_;
	Value = Value_;
	profile = true;
	vector = false;
	profileValue = profileValue_;
}
OutputLine::OutputLine(const Text Element1_,
					   const Text Element2_,
					   const Text Element3_,
					   const Text Value_,
					   const Vector profileValue_,
					   const int max_,
					   const int min_) {
	Element1 = Element1_;
	Element2 = Element2_;
	Element3 = Element3_;
	Value = Value_;
	profile = false;
	vector = true;
	profileValue = profileValue_;
	max = max_;
	min = min_;
}
OutputLine::~OutputLine() {
	Out();
}
Text
OutputLine::GetElement1() const {
	return Element1;
}
Text
OutputLine::GetElement2() const {
	return Element2;
}
Text
OutputLine::GetElement3() const {
	return Element3;
}
Text
OutputLine::GetValue() const {
	return Value;
}
bool
OutputLine::IsProfile() const {
	return profile;
}
Vector
OutputLine::GetProfile() const {
	return profileValue;
}
bool
OutputLine::IsVector() const {
	return vector;
}
Vector
OutputLine::GetVector() const {
	return profileValue;
}
int
OutputLine::GetMax() const {
	return max;
}
int
OutputLine::GetMin() const {
	return min;
}



