#include "Lat1DFlat.h"

Lat1DFlat::Lat1DFlat() {}
Lat1DFlat::Lat1DFlat(Input* MyInput_, Text name_) {
	MyInput = MyInput_;
	name = name_;
	if (*Copy(MyInput->GetText("lat",name,"geometry")) == *Copy("flat")) {
		Array<Text> param(1,13);
		param[1] = "dimensions";
		param[2] = "gradients";
		param[3] = "geometry";
		param[4] = "n_layers";
		param[5] = "lowerbound";
		param[6] = "upperbound";
		param[7] = "lambda";
		param[8] = "latticetype";
		param[9] = "bondlength";
		param[10] = "distance";
		param[11] = "layer_adjustment";
		param[12] = "n_layers_extended";
		param[13] = "lambda_chi";
		MyInput->CheckParameterNames("lat",name,param);
	}
	numLayers = MyInput->GetInt("lat",name,"n_layers",1,INT_MAX);
	M = numLayers+2;
	Array<Text> bounds(1,5);
	bounds[1] = "surface";
	bounds[2] = "mirror1";
	bounds[3] = "mirror2";
	bounds[4] = "periodic";
	bounds[5] = "bulk";
// default mirror1
	bound1 = MyInput->GetChoice("lat", name, "lowerbound", bounds, 2);
// default mirror1
	bound2 = MyInput->GetChoice("lat", name, "upperbound", bounds, 2);
	int numBulkBounds = 0;
	if (bound1 == 5) {
		bound1 = 2;
		bulkBound1 = true;
		numBulkBounds++;
	} else {
		bulkBound1 = false;
	}
	if (bound2 == 5) {
		bound2 = 2; //set real value later
		bulkBound2 = true;
		numBulkBounds++;
	} else {
		bulkBound2 = false;
	}
	if (numBulkBounds > 0) {
		BulkBoundaries.Dim(1,numBulkBounds);
	} else {
		BulkBoundaries.Dim(0,0);
	}
	numBulkBounds = 0;
	if (bulkBound1) {
		BulkBoundaries[++numBulkBounds] = Copy("lowerbound");
	}
	if (bulkBound2) {
		BulkBoundaries[++numBulkBounds] = Copy("upperbound");
	}
	if (bound1 == 4 && bound2 != 4) {
		if (MyInput->ValueSet("lat", name, "upperbound")) {
			Message(fatal,MyInput,"In 'lat : " + name
			+ "' lowerbound is periodic, while upperbound "
			"is not, this is not possible.");
		} else {
			bound2 = 4;
		}
	}
	if (bound2 == 4 && bound1 != 4) {
		if (MyInput->ValueSet("lat", name, "lowerbound")) {
			Message(fatal,MyInput,"In 'lat : " + name
			+ "' upperbound is periodic, while lowerbound "
			"is not, this is not possible.");
		} else {
			bound1 = 4;
		}
	}
	if (bound1 == 4) {
		bound1 = M-1;
	}
	if (bound2 == 4) {
		bound2 = 2;
	} else {
		bound2 = M - bound2 + 1; // set real value bound2
	}
	if (MyInput->ValueSet("lat",name,"n_layers_extended")) {
		int numLayersExtended = MyInput->GetInt("lat",name,
					"n_layers_extended",numLayers,INT_MAX);
		numExtraLatticeSites = numLayersExtended - numLayers;
		if (bound2 == M - 2) { // mirror2
			numExtraLatticeSites += 0.5;
		}
		Message(literal,"!!!WARNING!!! 'n_layers_extended' is experimental\n"
		"Check the thermodynamics carefully. Multistates probably wrong.");
	} else {
		numExtraLatticeSites = 0;
	}
	dimensions = MyInput->GetReal("lat", name, "dimensions", 1, DBL_MAX,3);
	Array<Text> latTypes(1,2);
	latTypes[1] = "standard";
	latTypes[2] = "stencils";
	int type;
	type = MyInput->GetChoice("lat", name, "latticetype",latTypes,1);
	latType = Copy(latTypes[type]);
	if (type == 1) {
		l1 = MyInput->GetReal("lat",name,"lambda",0,0.5);
	}
	if (type == 2) {
		l1 = 0.126041;
		if (MyInput->ValueSet("lat", name, "lambda")) {
			Message(fatal,"'lat : " + name + " : lambda' can only be set"
				" when 'latticetype' is set to standard");
		}
	}
	if (MyInput->ValueSet("lat",name,"lambda_chi")) {
		l1_chi = MyInput->GetReal("lat",name,"lambda_chi",0,0.5);
	} else {
		l1_chi = l1;
	}
	if (MyInput->ValueSet("lat",name,"bondlength")) {
		if (MyInput->ValueSet("lat",name,"distance")) {
			Message(fatal, "cannot set both 'lat : "
			+ name + " : bondlength' and 'distance'");
		}
		bondLength = MyInput->GetReal("lat",name,"bondlength",0,DBL_MAX,3e-10);
		siteDistance = bondLength/sqrt(2*dimensions*l1);
	} else if (MyInput->ValueSet("lat",name,"distance")) {
		Message(warning,"The use of 'lat : " + name
			+ " : distance' is deprecated, please use bondlength instead.");
		siteDistance = MyInput->GetReal("lat",name,"distance",0,DBL_MAX,3e-10);
		bondLength = siteDistance*sqrt(2*dimensions*l1);
	} else {
		bondLength = MyInput->GetReal("lat",name,"bondlength",0,DBL_MAX,3e-10);
		siteDistance = bondLength/sqrt(2*dimensions*l1);
	}
	layerAdjustment = MyInput->GetReal("lat",name,"layer_adjustment",0,2,1);
	if (MyInput->ValueSet("lat",name,"layer_adjustment")) {
		adjustOuterLayer = true;
		if (bound1 != 2 || bound2 != numLayers+1) {
			Message(fatal,"don't know how to adjust layers with "
			"these boundary conditions");
		}
	} else {
		adjustOuterLayer = false;
	}
}
Lat1DFlat::~Lat1DFlat() {
}
Text
Lat1DFlat::GetName() const {
	return name;
}
void
Lat1DFlat::GetOutput(Output* Out) const {
	Out->PutInt("lat",name,"gradients",1);
	if (dimensions - int(dimensions) > 1e-7) {
		Out->PutReal("lat",name,"dimensions",dimensions);
	} else {
		Out->PutInt("lat",name,"dimensions",int(dimensions));
	}
	Out->PutText("lat",name,"geometry","flat");
	Out->PutInt("lat",name,"n_layers",numLayers);
	Out->PutReal("lat",name,"n_extra_lattice_sites",numExtraLatticeSites);
	if (GetTotalNumLatticeSites() - int(GetTotalNumLatticeSites()) > 1e-7) {
		Out->PutReal("lat",name,"total number of lattice sites",
			     GetTotalNumLatticeSites());
	} else {
		Out->PutInt("lat",name,"total number of lattice sites",
			    int(GetTotalNumLatticeSites()));
	}
	Array<Text> bounds(1,5);
	bounds[1] = "surface";
	bounds[2] = "mirror1";
	bounds[3] = "mirror2";
	bounds[4] = "periodic";
	bounds[5] = "bulk";
	MyInput->SetDefaultWarningsOff();
	int bound1out = MyInput->GetChoice("lat", name, "lowerbound", bounds, 2);
	int bound2out = MyInput->GetChoice("lat", name, "upperbound", bounds, 2);
	MyInput->SetDefaultWarningsOn();
	Out->PutText("lat",name,"lowerbound",bounds[bound1out]);
	Out->PutText("lat",name,"upperbound",bounds[bound2out]);
	Out->PutReal("lat",name,"lambda",l1);
	Out->PutReal("lat",name,"lambda_chi",l1);
	Out->PutText("lat",name,"latticetype",latType);
	Out->PutReal("lat",name,"distance (m)",siteDistance);
	Text dim;
	if (dimensions - int(dimensions) > 1e-7) {
		dim = Blanks(100);
		dim.Putreal(dimensions-1,9);
		dim = Copy(dim.Strip().Frontstrip());
	} else {
		dim = Blanks(100);
		dim.Putint(int(dimensions-1));
		dim = Copy(dim.Strip().Frontstrip());
	}
	Out->PutReal("lat",name,"site surface (m^" + dim + ")",GetSiteSurface());
	Out->PutReal("lat",name,"volume (m)",GetVolume());
	Out->PutReal("lat",name,"layer_adjustment",layerAdjustment);
}
int
Lat1DFlat::GetTotalNumLayers() const {
	return numLayers + 2;
}

int
Lat1DFlat::Get_N_comp_ranges() const{
	return 1;
}


int
Lat1DFlat::Get_BL(int bn) const {
	switch (bn) {
	case 1:
		return bound1;
		break;
	case 2:
		return bound2;
		break;
	default:
		Message(fatal,"Asking a boundary layer outsize range in Lat1DFlat.");
		return 0;
		break;
	}
}



int
Lat1DFlat::GetR0(int co_or) const{
if 	(co_or ==1) return 1; else return 0;
}

int
Lat1DFlat::GetRm(int co_or) const{
	if 	(co_or ==1) return numLayers; else return 0;
}

int
Lat1DFlat::GetNumGradients() const {
	return 1;
}
int
Lat1DFlat::GetNumLayers(const int numGradient) const {
	if (numGradient == 1) return numLayers + 2;
	Message(fatal,"Argument out of range in call to "
			"'Lat1DFlat::GetNumLayers(int D)'");
	return -1; // never get here
}
double
Lat1DFlat::GetNumLatticeSites(const int) const {
	return 1;
}
double
Lat1DFlat::GetLambda(const int z1, const int z2) const {
	if (z1 == z2) return 1-2*l1;
	if (z1 - z2 == 1) return l1;
	if (z1 - z2 == -1) return l1;
	return 0;
}
void
Lat1DFlat::MultiplyWithLatticeSites(Vector A) const {
	A[2] *= layerAdjustment;
	A[M-1] *= 2-layerAdjustment;
}
void
Lat1DFlat::DivideByLatticeSites(Vector A) const {
	A[2] /= layerAdjustment;
	A[M-1] /= 2-layerAdjustment;
}
double
Lat1DFlat::GetVolume() const {
	return GetTotalNumLatticeSites()*siteDistance;
}
double
Lat1DFlat::GetTotalNumLatticeSites() const {
	double value = M;
	if (bound1 != 1) value--;
	if (bound1 == 3) value -= 0.5;
	if (bound2 != M) value--;
	if (bound2 == M-2) value -= 0.5;
	return value;
}
double
Lat1DFlat::GetNumExtraLatticeSites() const {
	return numExtraLatticeSites;
}
double
Lat1DFlat::GetSiteDistance() const {
	return siteDistance;
}
double
Lat1DFlat::GetSiteSurface() const {
	return pow(siteDistance,dimensions-1);
}
LatticeRange*
Lat1DFlat::NewLatticeRange(Text range) const {
	int i;
	int minX,maxX;
	char test;
	//mmm if the range is a file return LatticeRangeFile object
	Text t = Copy(range);
	if (t.Getchar()=='<') {
		t = t.Rest();
		return new LatticeRangeFile(t.Frontstrip(),GetTotalNumLayers());
	}
	//mmm
	if (*range == *Copy("free")) {
		return new LatticeRange1D(1, M, M);
	}
	Text low = Copy(range.Scanto(';'));
	Text high;
	if (low.Length() == range.Length()) high = Copy(low);
	else high = Copy(range.Scanto(';'));
	if (range.More()) Message(fatal,MyInput,"Error reading range '" +
		range + "', expecting one or two (separated by ';') coordinates");
	int length = low.Length();
	if (length == 0) {
		Message(fatal,MyInput,"Error reading range '"
		+ range + "', no value found for first coordinate");
	}
	if (*Copy(low) == *Copy("lowerbound")) minX = 1;
	else if (*Copy(low) == *Copy("upperbound")) minX = M;
	else if (*Copy(low) == *Copy("lastlayer")) minX = M-1;
	else {
		for (i=1; i<= length; i++) {
			test = low.Getchar();
			if (!isdigit(test)) {
				Message(fatal,MyInput,"Error reading range '" + range
				+ "', first coordinate should be a positive integer, "
				"'lowerbound', 'upperbound' or 'lastlayer'");
			}
		}
		low.Setpos(1);
		minX = low.Getint();
		if (minX < 0) {
			Message(fatal,MyInput,"Error reading range '" + range +
			"', coordinate cannot be negative");
		}
		if (minX > numLayers+1) {
			Message(fatal,MyInput,"Error reading range '" + range +
			"', coordinate exceeds n_layers+1");
		}
		minX++;
	}
	length = high.Length();
	if (length == 0) {
		Message(fatal,MyInput,"Error reading range '" + range
		+ "', no value found for last coordinate");
	}
	if (*Copy(high) == *Copy("upperbound")) maxX = M;
	else if (*Copy(high) == *Copy("lowerbound")) maxX = 1;
	else if (*Copy(high) == *Copy("lastlayer")) maxX = M-1;
	else {
		for (i=1; i<= length; i++) {
			test = high.Getchar();
			if (!isdigit(test)) {
				Message(fatal,MyInput,"Error reading range '" + range
				+ "', second coordinate should be a positive integer, "
				"'lowerbound', 'upperbound' or 'lastlayer'");
			}
		}
		high.Setpos(1);
		maxX = high.Getint();
		if (maxX < 0) {
			Message(fatal,MyInput,"Error reading range '" + range +
			"', coordinate cannot be negative");
		}
		if (maxX > numLayers+1) {
			Message(fatal,MyInput,"Error reading range '" + range +
			"', coordinate exceeds n_layers+1");
		}
		maxX++;
	}
	if (minX > maxX) {
		Message(hint,MyInput,"In range '" + range
		+ "' the first coordinate is higher than the last.\n"
		"These values will be swapped");
		int dummy = minX;
		minX = maxX;
		maxX = dummy;
	}
	return new LatticeRange1D(minX, maxX, M);
}
Boolean
Lat1DFlat::WithinBoundaries(const int z) const {
	if (z == 1) return false;
	if (z == M) return false;
	return true;
}
void
Lat1DFlat::SetBoundaries(Vector in) const {
	in[1] = in[bound1];
	in[M] = in[bound2];
}
void
Lat1DFlat::SetBulkBoundaries(Vector in, const Vector bulkValues) const {
	if (bulkBound1) {
		in[1] = bulkValues[1];
	}
	if (bulkBound2) {
		in[M] = bulkValues[2];
	}
}
void
Lat1DFlat::SubtractBoundaries(Vector in) const {
	if (bound1 != 1) in[1] = 0;
	if (bound1 == 3) in[2] /= 2;
	if (bound2 != M) in[M] = 0;
	if (bound2 == M-2) in[M-1] /= 2;
}
void
Lat1DFlat::RestoreBoundaries(Vector in) const {
	in[1] = in[bound1];
	if (bound1 == 3) in[2] *= 2;
	in[M] = in[bound2];
	if (bound2 == M-2) in[M-1] *= 2;
}
Array<Text>
Lat1DFlat::GetNamesBulkBoundaries() const {
	return BulkBoundaries;
}
Boolean
Lat1DFlat::BulkBoundary(const int z) const {
	if (z == 1 && bulkBound1 == true) {
		return true;
	}
	if (z == M && bulkBound2 == true) {
		return true;
	}
	return false;
}
int
Lat1DFlat::NumBulkBoundary(const int z) const {
	if (z == 1 && bulkBound1 == true) {
		return 1;
	}
	if (z == M && bulkBound2 == true) {
		return 2;
	}
	return 0;
}
void
Lat1DFlat::CheckBoundaries(Array<LatticeRange*> LatRanges) const {
	int first = LatRanges.Lowerbound();
	int last = LatRanges.Upperbound();
	int i,j;
	LatticeRange* LatRange;
	LatticeRange* LatRange2;
	for (i=first; i<=last; i++) {
		LatRange = LatRanges[i];
		if (LatRange->InRange(1)) {
			if (bound1 != 1) {
				Message(fatal,MyInput,Copy("'lat : ") + name
				+ " : lowerbound' is not set to 'surface' "
				"while a frozen segment is defined there");
			}
		}
		if (LatRange->InRange(M)) {
			if (bound2 != numLayers + 2) {
				Message(fatal,MyInput,Copy("'lat : ") + name
				+ " : upperbound' is not set to 'surface' "
				"while a frozen segment is defined there");
			}
		}
		for (j=i+1; j<=last; j++) {
		LatRange2 = LatRanges[j];
			for (int z=1; z<=M; z++) {
				if (LatRange->InRange(z) && LatRange2->InRange(z)) {
					Message(fatal,MyInput,"Two frozen segments "
					"are overlapping on the lattice");
				}
			}
		}
	}
	Boolean frozenFound = false;
	if (bound1 == 1) {
		for (i=first; i<=last; i++) {
			LatRange = LatRanges[i];
			if (LatRange->InRange(1)) frozenFound = true;
		}
		if (!frozenFound) {
			Message(fatal,MyInput,Copy("'lat : ") + name
			+ " : lowerbound' is set to 'surface' while "
			"no frozen segment is defined there");
		}
	}
	frozenFound = false;
	if (bound2 == M) {
		for (i=first; i<=last; i++) {
			LatRange = LatRanges[i];
			if (LatRange->InRange(M)) frozenFound = true;
		}
		if (!frozenFound) {
			Message(fatal,MyInput,Copy("'lat : ") + name
			+ " : upperbound' is set to 'surface' while "
			"no frozen segment is defined there");
		}
	}
}

void
Lat1DFlat:: Sides(Vector P, Vector Q) const{
	Message(fatal,"Sides not implemented");
}
double
Lat1DFlat::SideFraction(const Vector in, int z) const {
	if (z == 1) {
		if (bound1 == 1) return in[z] + l1_chi*(- in[z] - in[z] + in[z+1]);
		else if (bulkBound1) return in[z];
		else z = bound1;
	}
	if (z == M) {
		if (bound2 == M) return in[z] + l1_chi*(- in[z] - in[z] + in[z-1]);
		else if (bulkBound2) return in[z];
		else z = bound2;
	}
	return in[z] + l1_chi*(in[z-1] - in[z] - in[z] + in[z+1]);

}
void
Lat1DFlat::ElectricPotential(Vector &psi,Vector &psi0,const Vector &eps,const Vector &charge,const double preFactor) const {
	int z;
	double epsmin,epsplus;
	double pf = 2*preFactor;

	if (bound2 == M) {
		epsmin = eps[M-1]+eps[M];
		psi[M] = psi0[M-1] + pf*charge[M]/epsmin;
	} else {
		psi[M] = psi0[bound2];
	}
	epsplus = eps[1]+eps[2];
	if (bound1 == 1) {
		psi[1] = psi0[2] + pf*charge[1]/epsplus;
	} else {
		psi[1] = psi0[bound1];
	}
	for (z=2; z<M; z++) {
		epsmin = epsplus;
		epsplus = eps[z]+eps[z+1];
		psi[z] = (epsmin*psi0[z-1] + pf*charge[z] + epsplus*psi0[z+1])
				/(epsmin+epsplus);
	}
}
void
Lat1DFlat::ElectricFieldSquared(Vector &ESquared,
		const Vector psi,
		const double preFactor) const {
	int z;
	double Emin, Eplus,pf;
	pf = preFactor/2;
	Eplus = psi[1]-psi[2];
	Eplus *= Eplus;
	if (bound1 == 1) {
		ESquared[1] = pf*Eplus;
	}
	for (z=2; z<M; z++) {
		Emin = Eplus;
		Eplus = psi[z]-psi[z+1];
		Eplus *= Eplus;
		ESquared[z] = pf*(Emin + Eplus);
	}
	if (bound2 == M) {
		ESquared[M] = pf*Eplus;
	}
}
Vector
Lat1DFlat::Div1Grad2(const Vector in1, const Vector in2) const {
	Vector value(1,M);
	for (int z=1; z<=M; z++) {
		if (z > 1 && in1[z-1] != 0 && in1[z] != 0) {
			value[z] += (in1[z-1]+in1[z])*(in2[z-1]-in2[z]);
		}
		if (z < M && in1[z+1] != 0 && in1[z] != 0) {
			value[z] += (in1[z+1]+in1[z])*(in2[z+1]-in2[z]);
		}
		value[z] /= 2*siteDistance*siteDistance;
	}
	return value;
}
Vector
Lat1DFlat::FourierTransform(const Vector /*in*/) const {
	Message(implementation,"Lat1DFlat::FourierTransform");
	Vector value(1,M);
	return value;
}
double
Lat1DFlat::Moment(const Vector in,
				  const int moment,
				  const LatticeRange* Range,
				  const double offset) const {
	double value = 0;
	Vector copy(1,M);
	int z;
	for (z=1; z<=M; z++) {
		copy[z] = in[z];
	}
	SubtractBoundaries(copy);
	double Z;
	int zRange=0;
	Boolean found = false;
	Boolean foundagain = false;
	for (z=1; z<=M; z++) {
		if (Range->InRange(z)) {
			if (found) {
				foundagain = true;
			}
			zRange = z;
			found = true;
		}
	}
	if (foundagain) {
		Message(warning,"A moment is calculated wrong, because of a "
		"range that is more than one lattice layer");
	}
	if (!found) zRange = 0;
	for (z=1; z<=M; z++) {
		Z = z-zRange+offset+0.5;
		value += copy[z]*(pow(Z,moment+1)/(moment+1)-pow((Z-1),moment+1)/(moment+1));
	}
	return value;
}
double
Lat1DFlat::MomentUnweighted(const Vector in,
							const int moment,
							const LatticeRange* Range,
							const double offset) const {
	return Moment(in,moment,Range,offset);
}
Vector
Lat1DFlat::RenormPhi(Vector a,
						 Vector,
						 const Vector,
						 const double) const {
	Message(fatal,"Cannot have a molecule with freedom 'thirdGeneration' "
		"on a lattice with one gradient");
	return a; //never get here
}
double
Lat1DFlat::GetLayerAdjustment(void) const {
	return layerAdjustment;
}
void
Lat1DFlat::SetLayerAdjustment(double value) {
	adjustOuterLayer = true;
	layerAdjustment = value;
}

void Lat1DFlat::SetBx1(double *P) const {
	P[0]=P[bound1-1];
}

void Lat1DFlat::SetBx1(double *P,double *Q) const {
	if (bound1>2) {
		SetBx1(P); SetBx1(Q);
	}
	else {
		P[0]=Q[bound1-1];
		Q[0]=P[bound1-1];
	}
}

void Lat1DFlat::SetBxm(double *P) const {
	int mmx=M-1;
	P[mmx]=P[bound2-1];
}

void Lat1DFlat::SetBxm(double *P,double *Q) const {
	int mmx=M-1;
	if (bound2<numLayers-1) {
		SetBxm(P); SetBxm(Q);
	}
	else {
		P[mmx]=Q[bound2-1];
		Q[mmx]=P[bound2-1];
	}
}

void Lat1DFlat::removeboundaries(Vector P) const {double *pP=&P[1]; removeboundaries(pP); }
void Lat1DFlat::removeboundaries(double *P) const {
	int mmx=M-1;
	P[0]=P[mmx]=0;
}