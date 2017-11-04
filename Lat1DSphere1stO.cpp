#include "Lat1DSphere1stO.h"

Lat1DSphere1stO::Lat1DSphere1stO(Input* MyInput_, Text name_)
	: Lat1DFlat(MyInput_,name_), Lat1DCylinder(MyInput_,name_), Lat1DSphere(MyInput_,name_) {	
}
Lat1DSphere1stO::~Lat1DSphere1stO() {
}
