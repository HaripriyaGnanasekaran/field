#include "Lat1DSphere1stOS.h"

Lat1DSphere1stOS::Lat1DSphere1stOS(Input* MyInput_, Text name_)
	: Lat1DFlat(MyInput_,name_), Lat1DCylinder(MyInput_,name_), Lat1DSphere(MyInput_,name_) {
}
Lat1DSphere1stOS::~Lat1DSphere1stOS() {
}
