#include "Lat2DFlat.h"

// lattice keyword
const Text Lat2DFlat::lat="lat";

// allowed keywords
const Text Lat2DFlat::dimensionsS="dimensions";
const Text Lat2DFlat::gradients="gradients";
const Text Lat2DFlat::geometry="geometry";
const Text Lat2DFlat::n_layers_x="n_layers_x";
const Text Lat2DFlat::n_layers_y="n_layers_y";
const Text Lat2DFlat::lowerbound_x="lowerbound_x";
const Text Lat2DFlat::upperbound_x="upperbound_x";
const Text Lat2DFlat::lowerbound_y="lowerbound_y";
const Text Lat2DFlat::upperbound_y="upperbound_y";
const Text Lat2DFlat::distance="distance";
const Text Lat2DFlat::bondlength="bondlength";
const Text Lat2DFlat::lambda="lambda";
const Text Lat2DFlat::lambda0="lambda0";
const Text Lat2DFlat::lambda1="lambda1";
const Text Lat2DFlat::lambda2="lambda2";
const Text Lat2DFlat::latticetype="latticetype";

// choice options
const Text Lat2DFlat::standard="standard";
const Text Lat2DFlat::stencils="stencils";
const Text Lat2DFlat::hexagonal="hexagonal";

// default values
const double Lat2DFlat::l1stand=1.0/6;
const double Lat2DFlat::l0sten=4.34637428282995E-01;
const double Lat2DFlat::l1sten=1.16014619191836E-01;
const double Lat2DFlat::l2sten=2.53260237374153E-02;

Lat2DFlat::Lat2DFlat(Input* MyInput_, Text name_) {

	MyInput = MyInput_;
	name = name_;
	if (MyInput->GetText(lat,name,geometry) == "flat") {
		checkParameters();
	}
	numLayersX = MyInput->GetInt(lat,name,n_layers_x,1,INT_MAX);
	numLayersY = MyInput->GetInt(lat,name,n_layers_y,1,INT_MAX);
	MX = numLayersX + 2;
	MY = numLayersY + 2;
	M = MX*MY;
	Array<Text> bounds(1,5);
	bounds[1] = "surface";
	bounds[2] = "mirror1";
	bounds[3] = "mirror2";
	bounds[4] = "periodic";
	bounds[5] = "bulk";
	boundX1 = MyInput->GetChoice(lat, name, lowerbound_x, bounds, 2); // default mirror1
	boundX2 = MyInput->GetChoice(lat, name, upperbound_x, bounds, 2); // default mirror1
	boundY1 = MyInput->GetChoice(lat, name, lowerbound_y, bounds, 2); // default mirror1
	boundY2 = MyInput->GetChoice(lat, name, upperbound_y, bounds, 2); // default mirror1
	int numBulkBounds = 0;
	if (boundX1 == 5) {
		boundX1 = 2;
		bulkBoundX1 = true;
		numBulkBounds++;
	} else {
		bulkBoundX1 = false;
	}
	if (boundX2 == 5) {
		boundX2 = 2; //set real value later
		bulkBoundX2 = true;
		numBulkBounds++;
	} else {
		bulkBoundX2 = false;
	}
	if (boundY1 == 5) {
		boundY1 = 2;
		bulkBoundY1 = true;
		numBulkBounds++;
	} else {
		bulkBoundY1 = false;
	}
	if (boundY2 == 5) {
		boundY2 = 2; //set real value later
		bulkBoundY2 = true;
		numBulkBounds++;
	} else {
		bulkBoundY2 = false;
	}
	if ((bulkBoundX1 || bulkBoundX2) && (bulkBoundY1 || bulkBoundY2)) {
		Message(fatal, MyInput,"In 'lat : " + name
			+ "' a 'bulk' boundary is defined in two perpendicular directions. This is not allowed.");
	}
	if (numBulkBounds > 0) {
		BulkBoundaries.Dim(1,numBulkBounds);
	} else {
		BulkBoundaries.Dim(0,0);
	}
	numBulkBounds = 0;
	if (bulkBoundX1) {
		BulkBoundaries[++numBulkBounds] = Copy(lowerbound_x);
	}
	if (bulkBoundX2) {
		BulkBoundaries[++numBulkBounds] = Copy(upperbound_x);
	}
	if (bulkBoundY1) {
		BulkBoundaries[++numBulkBounds] = Copy(lowerbound_y);
	}
	if (bulkBoundY2) {
		BulkBoundaries[++numBulkBounds] = Copy(upperbound_y);
	}
	if (boundX1 == 4 && boundX2 != 4) {
		if (MyInput->ValueSet(lat, name, upperbound_x)) {
			Message(fatal,MyInput,"In 'lat : " + name
				+ "' lowerbound_x is periodic, while upperbound_x is not, "
				"this is not possible.");
		} else {
			boundX2 = 4;
		}
	}
	if (boundX2 == 4 && boundX1 != 4) {
		if (MyInput->ValueSet(lat, name, lowerbound_x)) {
			Message(fatal,MyInput,"In 'lat : " + name + "' upperbound_x is "
				"periodic, while lowerbound_x is not, this is not possible.");
		} else {
			boundX1 = 4;
		}
	}
	if (boundY1 == 4 && boundY2 != 4) {
		if (MyInput->ValueSet(lat, name, upperbound_y)) {
			Message(fatal,MyInput,"In 'lat : " + name + "' lowerbound_y is "
				"periodic, while upperbound_y is not, this is not possible.");
		} else {
			boundY2 = 4;
		}
	}
	if (boundY2 == 4 && boundY1 != 4) {
		if (MyInput->ValueSet(lat, name, lowerbound_y)) {
			Message(fatal,MyInput,"In 'lat : " + name + "' upperbound_y is "
				"periodic, while lowerbound_y is not, this is not possible.");
		} else {
			boundY1 = 4;
		}
	}
	if (boundX1 == 4) boundX1 = numLayersX+1;
	if (boundX2 == 4) boundX2 = 2;
	else boundX2 = numLayersX + 3 - boundX2;
	if (boundY1 == 4) boundY1 = numLayersY+1;
	if (boundY2 == 4) boundY2 = 2;
	else boundY2 = numLayersY + 3 - boundY2;
	dimensions = MyInput->GetReal(lat, name, dimensionsS, 2, DBL_MAX,3);

	int lattyp = 1; // default standard lattice
	if (MyInput->ValueSet(lat,name,latticetype)) {
		Array<Text> lattype(1,2);
		lattype[1] = standard;
		lattype[2] = stencils;
		lattyp = MyInput->GetChoice(lat, name, latticetype, lattype, 1);
	}
	MyInput->DontCombineParam(lat,name,lambda,lambda0);
	MyInput->DontCombineParam(lat,name,lambda,lambda1);
	MyInput->DontCombineParam(lat,name,lambda,lambda2);
	MyInput->AlwaysCombineParam(lat,name,lambda0,lambda1);
	MyInput->AlwaysCombineParam(lat,name,lambda0,lambda2);
	if (lattyp == 1) {
		latType = standard;
		l1 = MyInput->GetReal(lat,name,lambda,0,0.5,l1stand);
		l1x=l1;
		l1y=l1;
		l0 = 1 - 4*l1;
		l2 = 0;
		if (MyInput->ValueSet(lat, name,lambda0)) {
			Message(fatal, MyInput, "'lambda0' cannot be defined with "
				"latticetype 'standard'.");
		}
	} else {
		latType = stencils;
		if (MyInput->ValueSet(lat,name,"HII")) {
			if (MyInput->GetBoolean(lat,name,"HII",false)) {lattyp =3;}
		}
		if (lattyp ==3) {
			l2=1/9.0;
			l1y=2.0*numLayersX*sqrt(3.0)/(3.0*(numLayersX*sqrt(3.0)+ numLayersY)) -2*l2;
			l1x=2.0* numLayersY/(3.0*(numLayersX*sqrt(3.0)+ numLayersY)) -2*l2;
			l1=(l1x+l1y)/2.0;
			l0=1-2*l1x-2*l1y-4*l2;

		} else {
			l0 = MyInput->GetReal(lat,name,lambda0,0,0.5,l0sten);
			l1 = MyInput->GetReal(lat,name,lambda1,0,0.5,l1sten);
			l1x=l1;
			l1y=l1;
			l2 = MyInput->GetReal(lat,name,lambda2,0,0.5,l2sten);
		}
		// check whether total of lambdas is ~1
		// i1 will be smallest computer number > 1
		double i1 = 4*l1 + 4*l2 + l0 - 1;
		i1 = (i1 < 0) ? -i1 : i1;
		if (i1 > DBL_EPSILON) {
			// error: total of lambdas in not ~1
			Message(fatal, MyInput, "Sum of lambda's is not 1.");
		}
		// make sure the sum of the lambdas is exactly 1.
		l0 = 1 - 4*l1 - 4*l2;
	}

	if (MyInput->ValueSet(lat,name,bondlength)) {
		if (MyInput->ValueSet(lat,name,distance)) {
			Message(fatal, "cannot set both 'lat : " + name
				+ " : bondlength' and 'distance'");
		}
		bondLength = MyInput->GetReal(lat, name, bondlength, 0, DBL_MAX, 3e-10);
		siteDistance = bondLength/sqrt(2*dimensions*l1);
	} else if (MyInput->ValueSet(lat, name, distance)) {
		Message(warning,"The use of 'lat : " + name
			+ " : distance' is depreciated, please use bondlength instead.");
		siteDistance = MyInput->GetReal(lat, name, distance, 0, DBL_MAX, 3e-10);
		bondLength = siteDistance*sqrt(2*dimensions*l1);
	} else {
		bondLength = MyInput->GetReal(lat, name, bondlength, 0, DBL_MAX, 3e-10);
		siteDistance = bondLength/sqrt(2*dimensions*l1);
	}
}

Lat2DFlat::~Lat2DFlat() {
}

// protected function
void Lat2DFlat::checkParameters() const {
	Array<Text> param(1,16);
	param[1] = dimensionsS;
	param[2] = gradients;
	param[3] = geometry;
	param[4] = n_layers_x;
	param[5] = n_layers_y;
	param[6] = lowerbound_x;
	param[7] = upperbound_x;
	param[8] = lowerbound_y;
	param[9] = upperbound_y;
	param[10] = distance;
	param[11] = bondlength;
	param[12] = lambda;
	param[13] = lambda0;
	param[14] = lambda1;
	param[15] = lambda2;
	param[16] = latticetype;
	MyInput->CheckParameterNames(lat,name,param);
}

Text
Lat2DFlat::GetName() const {
	return name;
}

void
Lat2DFlat::GetOutput(Output* Out) const {
	Out->PutInt(lat,name,"gradients",2);
	if (dimensions - int(dimensions) > 1e-7) {
		Out->PutReal(lat,name,"dimensions",dimensions);
	} else {
		Out->PutInt(lat,name,"dimensions",int(dimensions));
	}
	Out->PutText(lat,name,"geometry","flat");
	Out->PutInt(lat,name,"n_layers_x",numLayersX);
	Out->PutInt(lat,name,"n_layers_y",numLayersY);
	if (GetTotalNumLatticeSites() - int(GetTotalNumLatticeSites()) > 1e-7) {
		Out->PutReal(lat,name,"total number of lattice sites",GetTotalNumLatticeSites());
	} else {
		Out->PutInt(lat,name,"total number of lattice sites",int(GetTotalNumLatticeSites()));
	}
	Array<Text> bounds(1,5);
	bounds[1] = "surface";
	bounds[2] = "mirror1";
	bounds[3] = "mirror2";
	bounds[4] = "periodic";
	bounds[5] = "bulk";
	MyInput->SetDefaultWarningsOff();
	int boundX1out = MyInput->GetChoice(lat, name, "lowerbound_x", bounds, 2); // default mirror1
	int boundX2out = MyInput->GetChoice(lat, name, "upperbound_x", bounds, 2); // default mirror1
	int boundY1out = MyInput->GetChoice(lat, name, "lowerbound_y", bounds, 2); // default mirror1
	int boundY2out = MyInput->GetChoice(lat, name, "upperbound_y", bounds, 2); // default mirror1
	MyInput->SetDefaultWarningsOn();
	Out->PutText(lat,name,"lowerbound_x",bounds[boundX1out]);
	Out->PutText(lat,name,"upperbound_x",bounds[boundX2out]);
	Out->PutText(lat,name,"lowerbound_y",bounds[boundY1out]);
	Out->PutText(lat,name,"upperbound_y",bounds[boundY2out]);

	if (MyInput->ValueSet(lat,name,lambda)) {
		Out->PutReal(lat,name,lambda,l1);
	} else {
		Out->PutReal(lat, name, lambda0, l0);
		Out->PutReal(lat, name, lambda1, l1);
		Out->PutReal(lat, name, lambda2, l2);
		Out->PutReal(lat, name, "lambda1x", l1x);
		Out->PutReal(lat, name, "lambda1y", l1y);
	}
	Out->PutText(lat,name,latticetype,latType);

	Out->PutReal(lat,name,"distance (m)",siteDistance);
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
	Out->PutReal(lat,name,"site surface (m^" + dim + ")",GetSiteSurface());
	Out->PutReal(lat,name,"volume (m^2)",GetVolume());
}

int
Lat2DFlat::GetTotalNumLayers() const {
	return M;
}

int
Lat2DFlat::Get_N_comp_ranges() const{
	return 1;
}


int
Lat2DFlat::Get_BL(int bn) const {
	switch (bn) {
	case 1:
		return boundX1;
		break;
	case 2:
		return boundX2;
		break;
	case 3:
		return boundY1;
		break;
	case 4:
		return boundY2;
		break;
	default:
		Message(fatal,"Asking a boundary layer outsize range in Lat2DFlat.");
		return 0;
		break;
	}
}


int
Lat2DFlat::GetR0(int co_or) const{
if 	(co_or ==1||co_or==2) return 1; else return 0;
}

int
Lat2DFlat::GetRm(int co_or) const{
	if 	(co_or ==1) return numLayersX; else
		if (co_or ==2) return numLayersY; else return 0;
}

int
Lat2DFlat::GetNumGradients() const {
	return 2;
}
int
Lat2DFlat::GetNumLayers(const int numGradient) const {
	if (numGradient == 1) return MX;
	if (numGradient == 2) return MY;
	Message(fatal,"Argument out of range in call to 'Lat2DFlat::GetNumLayers(int D)'");
	return -1; // never get here
}
double
Lat2DFlat::GetNumLatticeSites(const int) const {
	return 1;
}
double
Lat2DFlat::GetLambda(const int z1, const int z2) const {
// Watch out: this function doesn't take boundaries into account.
	int diff = abs(z1 - z2);
	if (diff == 0) {
		return l0;
	} else if (diff == 1 || diff == MY) {
		return l1;
	} else if (diff == MY+1 || diff == MY-1) {
		return l2;
	} else return 0;
}
void
Lat2DFlat::MultiplyWithLatticeSites(Vector) const {
}
void
Lat2DFlat::DivideByLatticeSites(Vector) const {
}
double
Lat2DFlat::GetVolume() const {
	return GetTotalNumLatticeSites()*siteDistance*siteDistance;
}
double
Lat2DFlat::GetTotalNumLatticeSites() const {
	double valueX = MX;
	if (boundX1 != 1) valueX--;
	if (boundX1 == 3) valueX -= 0.5;
	if (boundX2 != MX) valueX--;
	if (boundX2 == MX-2) valueX -= 0.5;
	double valueY = MY;
	if (boundY1 != 1) valueY--;
	if (boundY1 == 3) valueY -= 0.5;
	if (boundY2 != MY) valueY--;
	if (boundY2 == MY-2) valueY -= 0.5;
	return valueX*valueY;
}
double
Lat2DFlat::GetNumExtraLatticeSites() const {
	return 0;
}
double
Lat2DFlat::GetSiteDistance() const {
	return siteDistance;
}
double
Lat2DFlat::GetSiteSurface() const {
	return pow(siteDistance,dimensions-1);
}
LatticeRange*
Lat2DFlat::NewLatticeRange(Text range) const {
	int i;
	Boolean error = false;
	int minX,maxX;
	int minY,maxY;
	Text xValue, yValue;
	char test;
	//mmm if the range is a file return a LatticeRangeFile object
	Text t = Copy(range);
	if (t.Getchar()=='<') {
		t = t.Rest();
		return new LatticeRangeFile(t.Frontstrip(),GetTotalNumLayers());
	}
	//mmm

	t.Setpos(1);
	if (t.Getchar()=='=') {
		t = t.Rest();
		return new LatticeRangeFile(t.Frontstrip(),numLayersX,numLayersY);
	}

	if (*range == *Copy("free")) return new LatticeRange2D(1,MX,MX,1,MY,MY);
	Text low = Copy(range.Scanto(';'));
	Text high;
	if (low.Length() == range.Length()) high = Copy(low);
	else high = Copy(range.Scanto(';'));
	if (range.More()) {
		Message(fatal,MyInput,"Error reading range '"
		+ range + "', expecting two coordinates (separated by ';')");
	}
	if (low.Length() == 0) {
		Message(fatal,MyInput,"Error reading range '" + range
		+ "', no value found for first coordinate");
	}
	xValue = Copy(low.Scanto(','));
	xValue = Copy(xValue.Strip());
	xValue = Copy(xValue.Frontstrip());
	if (!low.More()) error = true;
	yValue = Copy(low.Scanto(','));
	yValue = Copy(yValue.Strip());
	yValue = Copy(yValue.Frontstrip());
	if (low.More()) error = true;
	if (xValue.Length() == 0) error = true;
	if (yValue.Length() == 0) error = true;
	if (error) {
		Message(fatal,MyInput,"Error reading the first coordinate (" + low
			+ ") of range '" + range + "', expecting two positive integers,"
			"'lowerbound', 'upperbound' or 'lastlayer' separated by a comma");
	}
	if (*xValue == *Copy("lowerbound")) {
		minX = 0;
	} else if (*xValue == *Copy("upperbound")) {
		minX = MX-1;
	} else if (*xValue == *Copy("lastlayer")) {
		minX = MX-2;
	} else {
		int length = xValue.Length();
		for (i=1; i<= length; i++) {
			test = xValue.Getchar();
			if (!isdigit(test)) {
				Message(fatal,MyInput,"Error reading range '" + range
				+ "', expecting an integer, 'lowerbound', 'upperbound' or 'lastlayer' for first x coordinate");
			}
		}
		xValue.Setpos(1);
		minX = xValue.Getint();
		if (minX < 0) {
			Message(fatal,MyInput,"Error reading range '" + range
			+ "', value cannot be negative");
		}
		if (minX > numLayersX+1) {
			Message(fatal,MyInput,"Error reading range '" + range
			+ "', value for first x coordinate exceeds n_layers_x+1");
		}
	}
	if (*yValue == *Copy("lowerbound")) {
		minY = 0;
	} else if (*yValue == *Copy("upperbound")) {
		minY = MY-1;
	} else if (*yValue == *Copy("lastlayer")) {
		minY = MY-2;
	} else {
		int length = yValue.Length();
		for (i=1; i<= length; i++) {
			test = yValue.Getchar();
			if (!isdigit(test)) {
				Message(fatal,MyInput,"Error reading range '" + range
				+ "', expecting an integer, 'lowerbound', 'upperbound' or 'lastlayer' for first y coordinate");
			}
		}
		yValue.Setpos(1);
		minY = yValue.Getint();
		if (minY < 0) {
			Message(fatal,MyInput,"Error reading range '" + range
			+ "', value cannot be negative");
		}
		if (minY > numLayersY+1) {
			Message(fatal,MyInput,"Error reading range '" + range
			+ "', value for first y coordinate exceeds n_layers_y+1");
		}
	}
	minX++;
	minY++;
	if (high.Length() == 0) {
		Message(fatal,MyInput,"Error reading range '" + range
		+ "', no value found for first coordinate");
	}
	xValue = Copy(high.Scanto(','));
	xValue = Copy(xValue.Strip());
	xValue = Copy(xValue.Frontstrip());
	if (!high.More()) error = true;
	yValue = Copy(high.Scanto(','));
	yValue = Copy(yValue.Strip());
	yValue = Copy(yValue.Frontstrip());
	if (high.More()) error = true;
	if (xValue.Length() == 0) error = true;
	if (yValue.Length() == 0) error = true;
	if (error) {
		Message(fatal,MyInput,"Error reading the second coordinate of range '" + range
		+ "', expecting two positive integers, 'lowerbound', 'upperbound' or 'lastlayer' separated by a comma");
	}
	if (*xValue == *Copy("upperbound")) {
		maxX = MX-1;
	} else if (*xValue == *Copy("lowerbound")) {
		maxX = 0;
	} else if (*xValue == *Copy("lastlayer")) {
		maxX = MX-2;
	} else {
		int length = xValue.Length();
		for (i=1; i<= length; i++) {
			test = xValue.Getchar();
			if (!isdigit(test)) {
				Message(fatal,MyInput,"Error reading range '" + range
				+ "', expecting an integer, 'lowerbound', 'upperbound' or 'lastlayer' for second x coordinate");
			}
		}
		xValue.Setpos(1);
		maxX = xValue.Getint();
		if (maxX < 0) {
			Message(fatal,MyInput,"Error reading range '" + range
			+ "', value cannot be negative");
		}
		if (maxX > numLayersX+1) {
			Message(fatal,MyInput,"Error reading range '" + range
			+ "', value for second x coordinate exceeds n_layers_x+1");
		}
	}
	if (*yValue == *Copy("upperbound")) {
		maxY = MY-1;
	} else if (*yValue == *Copy("lowerbound")) {
		maxY = 0;
	} else if (*yValue == *Copy("lastlayer")) {
		maxY = MY-2;
	} else {
		int length = yValue.Length();
		for (i=1; i<= length; i++) {
			test = yValue.Getchar();
			if (!isdigit(test)) {
				Message(fatal,MyInput,"Error reading range '" + range
				+ "', expecting an integer, 'lowerbound', 'upperbound' or 'lastlayer' for second y coordinate");
			}
		}
		yValue.Setpos(1);
		maxY = yValue.Getint();
		if (maxY < 0) {
			Message(fatal,MyInput,"Error reading range '" + range
			+ "', value cannot be negative");
		}
		if (maxY > numLayersY+1) {
			Message(fatal,MyInput,"Error reading range '" + range
			+ "', value for second y coordinate exceeds n_layers_y+1");
		}
	}
	maxX++;
	maxY++;
	if (minX > maxX) {
		Message(hint,MyInput,"In range '" + range
		+ "' the first x coordinate is higher than the last.\nThese values will be swapped");
		int dummy = minX;
		minX = maxX;
		maxX = dummy;
	}
	if (minY > maxY) {
		Message(hint,MyInput,"In range '" + range
		+ "' the first y coordinate is higher than the last.\nThese values will be swapped");
		int dummy = minY;
		minY = maxY;
		maxY = dummy;
	}
	return new LatticeRange2D(minX, maxX, MX, minY, maxY, MY);
}
Boolean
Lat2DFlat::WithinBoundaries(const int z) const {
	if (z <= MY) return false;
	if (z >= M - MY) return false;
	if (z%MY == 1) return false;
	if (z%MY == 0) return false;
	return true;
}
void
Lat2DFlat::SetBoundaries(Vector in) const {
	int z;
	// set boundaries without corners
	// set lower boundary
	if (boundX1 != 1) {
		for (z=2; z<MY; z++) {
			in[z] = in[(boundX1-1)*MY+z];
		}
	}
	// set upper boundary
	if (boundX2 != MX) {
		for (z=2; z<MY; z++) {
			in[M-MY+z] = in[(boundX2-1)*MY+z];
		}
	}
	// set left boundary
	if (boundY1 != 1) {
		for (z=MY; z<M-MY; z+=MY) {
			in[z+1] = in[z+boundY1];
		}
	}
	// set right boundary
	if (boundY2 != MY) {
		for (z=MY+MY; z<M; z+=MY) {
			in[z] = in[z-MY+boundY2];
		}
	}

	// Set corners.
	// If one boundary is periodic or mirror, it's easy.
	// Else the corner is defined already.
	if (boundY1 != 1 && !bulkBoundY1) {
		in[1] = in[boundY1];
		in[M-MY+1] = in[M-MY+boundY1];
	}
	if (boundY2 != MY && !bulkBoundY2) {
		in[MY] = in[boundY2];
		in[M] = in[M-MY+boundY2];
	}
	if (boundX1 != 1 && !bulkBoundX1) {
		in[1] = in[(boundX1-1)*MY+1];
		in[MY] = in[boundX1*MY];
	}
	if (boundX2 != MX && !bulkBoundX2) {
		in[M-MY+1] = in[(boundX2-1)*MY+1];
		in[M] = in[boundX2*MY];
	}
}
void
Lat2DFlat::SetBulkBoundaries(Vector in, const Vector bulkValues) const {
	int z;
	if (bulkBoundX1) {
		for (z=2; z<MY; z++) {
			in[z] = bulkValues[1];
		}
		if (boundY1 != 1) {
			in[1] = bulkValues[1];
		}
		if (boundY2 != MY) {
			in[MY] = bulkValues[1];
		}
	}
	if (bulkBoundX2) {
		for (z=M-MY+2; z<M; z++) {
			in[z] = bulkValues[2];
		}
		if (boundY1 != 1) {
			in[M-MY+1] = bulkValues[2];
		}
		if (boundY2 != MY) {
			in[M] = bulkValues[2];
		}
	}
	if (bulkBoundY1) {
		for (z=MY+1; z<M-MY; z+=MY) {
			in[z] = bulkValues[3];
		}
		if (boundX1 != 1) {
			in[1] = bulkValues[3];
		}
		if (boundX2 != MX) {
			in[M-MY] = bulkValues[3];
		}
	}
	if (bulkBoundY2) {
		for (z=MY+MY; z<M; z+=MY) {
			in[z] = bulkValues[4];
		}
		if (boundX1 != 1) {
			in[MY] = bulkValues[4];
		}
		if (boundX2 != MX) {
			in[M] = bulkValues[4];
		}
	}
}
void
Lat2DFlat::SubtractBoundaries(Vector in) const {
	int z=0;
	for (int x=1; x<=MX; x++) {
		for (int y=1; y<=MY; y++) {
			z++;
			if (x == 1 && boundX1 != 1) in[z] = 0;
			else if (x == 2 && boundX1 == 3) in[z] /= 2;
			else if (x == MX-1 && boundX2 == MX-2) in[z] /= 2;
			else if (x == MX && boundX2 != MX) in[z] = 0;
			if (y == 1 && boundY1 != 1) in[z] = 0;
			else if (y == 2 && boundY1 == 3) in[z] /= 2;
			else if (y == MY-1 && boundY2 == MY-2) in[z] /= 2;
			else if (y == MY && boundY2 != MY) in[z] = 0;
		}
	}
}
void
Lat2DFlat::RestoreBoundaries(Vector in) const {
	int z;
	// Restore mirror2 values.
	if (boundY1 == 3) {
		for (z=MY+2; z<M-MY; z+=MY) {
			in[z] *= 2;
		}
	}
	if (boundY2 == MY-2) {
		for (z=MY+MY-1; z<M-MY; z+=MY) {
			in[z] *= 2;
		}
	}
	if (boundX1 == 3) {
		for (z=MY+2; z<MY+MY; z++) {
			in[z] *= 2;
		}
	}
	if (boundX2 == MX-2) {
		for (z=M-MY-MY+2; z<M-MY; z++) {
			in[z] *= 2;
		}
	}
	// Restore other boundaries.
	SetBoundaries(in);
}

Array<Text>
Lat2DFlat::GetNamesBulkBoundaries() const {
	return BulkBoundaries;

}

Boolean
Lat2DFlat::BulkBoundary(const int z) const {
	// corner can never be a bulkBoundary except in case of bulkerror
	if (z == 1 || z == MY || z == M-MY+1 || z == M) {
		return false;
	}
	if (z < MY && bulkBoundX1 == true) return true;
	if (z > M-MY && bulkBoundX2 == true) return true;
	if (z%MY == 1 && bulkBoundY1 == true) return true;
	if (z%MY == 0 && bulkBoundY2 == true) return true;
	return false;
}
int
Lat2DFlat::NumBulkBoundary(const int z) const {
	if (z == 1 || z == MY || z == M-MY+1 || z == M) {
		return false;
	}
	if (z <= MY && bulkBoundX1 == true) return 1;
	if (z >= M - MY && bulkBoundX2 == true) return 2;
	if (z%MY == 1 && bulkBoundY1 == true) return 3;
	if (z%MY == 0 && bulkBoundY2 == true) return 4;
	return false;
}

/*
Checks whether the surface exists on every place it is defined:
Do the input latticeranges correspond to the set boundaries?
*/
void
Lat2DFlat::CheckBoundaries(Array<LatticeRange*> LatRanges) const {
	// !function does not yet check for frozen segments defined outside boundaries!
	int first = LatRanges.Lowerbound();
	int last = LatRanges.Upperbound();
	int i,j,z;
	LatticeRange* LatRange;
	int x,y;
	Boolean frozenFound = false;
	if (boundX1 == 1) {
		for (i=first; i<=last; i++) {
			LatRange = LatRanges[i];
			for (y=2; y<=MY-1; y++) {
				z = y;
				if (LatRange->InRange(z)) frozenFound = true;
			}
		}
		if (!frozenFound) {
			Message(fatal,MyInput,"'lat : " + name
			+ " : lowerbound_x' is set to 'surface' while no frozen segment is defined there");
		}
	}
	frozenFound = false;
	if (boundX2 == MX) {
		for (i=first; i<=last; i++) {
			LatRange = LatRanges[i];
			for (y=2; y<=MY-1; y++) {
				z = y + M - MX;
				if (LatRange->InRange(z)) frozenFound = true;
			}
		}
		if (!frozenFound) {
			Message(fatal,MyInput,"'lat : " + name
			+ " : upperbound_x' is set to 'surface' while no frozen segment is defined there");
		}
	}
	frozenFound = false;
	if (boundY1 == 1) {
		for (i=first; i<=last; i++) {
			LatRange = LatRanges[i];
			for (x=2; x<=MX-1; x++) {
				z = (x-1)*MY + 1;
				if (LatRange->InRange(z)) frozenFound = true;
			}
		}
		if (!frozenFound) {
			Message(fatal,MyInput,"'lat : " + name
			+ " : lowerbound_y' is set to 'surface' while no frozen segment is defined there");
		}
	}
	frozenFound = false;
	if (boundY2 == MY) {
		for (i=first; i<=last; i++) {
			LatRange = LatRanges[i];
			for (x=2; x<=MX-1; x++) {
				z = x*MY;
				if (LatRange->InRange(z)) frozenFound = true;
			}
		}
		if (!frozenFound) {
			Message(fatal,MyInput,"'lat : " + name
			+ " : upperbound_y' is set to 'surface' while no frozen segment is defined there");
		}
	}
	Text number;
	if (boundX1 == 1) {
		for (y=2; y<=MY-1; y++) {
			frozenFound = false;
			for (i=first; i<=last; i++) {
				LatRange = LatRanges[i];
				z = y;
				if (LatRange->InRange(z)) frozenFound = true;
			}
			if (!frozenFound) {
				number = Blanks(100);
				number.Putint(y-1);
				number = Copy(number.Frontstrip());
				Message(fatal,MyInput,"'lat : " + name
				+ " : lowerbound_x' is set to 'surface' while no frozen segment is defined there at y-coordinate "
				+ number);
			}
		}
	}
	if (boundX2 == MX) {
		for (y=2; y<=MY-1; y++) {
			frozenFound = false;
			for (i=first; i<=last; i++) {
				LatRange = LatRanges[i];
				z = y + M - MX;
				if (LatRange->InRange(z)) frozenFound = true;
			}
			if (!frozenFound) {
				number = Blanks(100);
				number.Putint(y-1);
				number = Copy(number.Frontstrip());
				Message(fatal,MyInput,"'lat : " + name
				+ " : upperbound_x' is set to 'surface' while no frozen segment is defined there at y-coordinate "
				+ number);
			}
		}
	}

	if (boundY1 == 1) {
		for (x=2; x<=MX-1; x++) {
			frozenFound = false;
			for (i=first; i<=last; i++) {
				LatRange = LatRanges[i];
				z = (x-1)*MY + 1;
				if (LatRange->InRange(z)) frozenFound = true;
			}
			if (!frozenFound) {
				number = Blanks(100);
				number.Putint(x-1);
				number = Copy(number.Frontstrip());
				Message(fatal,MyInput,"'lat : " + name
				+ " : lowerbound_y' is set to 'surface' while no frozen segment is defined there at x-coordinate "
				+ number);
			}
		}
	}
	if (boundY2 == MY) {
		for (x=2; x<=MX-1; x++) {
			frozenFound = false;
			for (i=first; i<=last; i++) {
				LatRange = LatRanges[i];
				z = x*MY;
				if (LatRange->InRange(z)) frozenFound = true;
			}
			if (!frozenFound) {
				number = Blanks(100);
				number.Putint(x-1);
				number = Copy(number.Frontstrip());
				Message(fatal,MyInput,"'lat : " + name
				+ " : upperbound_y' is set to 'surface' while no frozen segment is defined there at x-coordinate "
				+ number);
			}
		}
	}
	if (boundX1 == 1 && boundY1 == 1) {
		frozenFound = false;
		for (i=first; i<=last; i++) {
			if (LatRanges[i]->InRange(1)) {
				frozenFound = true;
			}
		}
		if (!frozenFound) {
			Message(fatal, MyInput, "'lat : " + name
			+ " : lowerbound_x and lowerbound_y are set while no frozen segment is defined at (1,1).");
		}
	}
	if (boundX2 == MX && boundY1 == 1) {
		frozenFound = false;
		for (i=first; i<=last; i++) {
			if (LatRanges[i]->InRange(M-MY+1)) {
				frozenFound = true;
			}
		}
		if (!frozenFound) {
			Message(fatal, MyInput, "'lat : " + name
			+ " : upperbound_x and lowerbound_y are set while no frozen segment is defined at (MX,1).");
		}
	}
	if (boundX1 == 1 && boundY2 == MY) {
		frozenFound = false;
		for (i=first; i<=last; i++) {
			if (LatRanges[i]->InRange(MY)) {
				frozenFound = true;
			}
		}
		if (!frozenFound) {
			Message(fatal, MyInput, "'lat : " + name
			+ " : lowerbound_x and upperbound_y are set while no frozen segment is defined at (1,MY).");
		}
	}
	if (boundX2 == MX && boundY2 == MY) {
		frozenFound = false;
		for (i=first; i<=last; i++) {
			if (LatRanges[i]->InRange(M)) {
				frozenFound = true;
			}
		}
		if (!frozenFound) {
			Message(fatal, MyInput, "'lat : " + name
			+ " : upperbound_x and upperbound_y are set while no frozen segment is defined at (MX,MY).");
		}
	}
	// Check for overlap in LatRanges.
	for (i=first; i<last; i++) {
		for (j=i+1; j<=last; j++) {
			for (z=1; z<=M; z++) {
				if (LatRanges[i]->InRange(z) && LatRanges[j]->InRange(z)) {
					Message(fatal, MyInput, "'lat : " + name
						+ "' contains two overlapping lattice ranges:\n"
						+ LatRanges[i]->GetOutput() + "\n"
						+ LatRanges[i]->GetOutput());
				}
			}
		}
	}
}

double
Lat2DFlat::SideFraction(const Vector in, int z) const {
	// The input vector contains the data for which the sidefraction needs to be calculated.
	// The input data ranges from 1 to MX*MY

	// It is not yet possible to have a bulkboundary in a x and y position.
	// Positions for where to put
	// the code for this are given at the flag `bulkerror'.
	assert(!((bulkBoundX1 || bulkBoundX2) && (bulkBoundY1 || bulkBoundY2)));

	// If z is in a corner...
	// Either return a value or give z a value not in a corner.
	if (z==1) {
		if (!bulkBoundX1 && boundX1 != 1) { // If xbound is mirror or periodic...
			z = 1+(boundX1-1)*MY;
		} else if (!bulkBoundY1 && boundY1 != 1) { // Same for yborder.
			z = boundY1;
		} else if (boundX1 == 1 && boundY1 == 1) { // If pos. 1 belongs to system...
			return l0*in[1] + l1y*in[2]+l1x*in[MY+1] + l2*in[MY+2];
		} else if (boundX1 == 1) { // bulkBoundY1 == true
			return (l0+l1+l1)*in[1] + (l1+l2+l2)*in[MY+1];
		} else if (boundY1 == 1) { // bulkBoundX1 == true
			return (l0+l1+l1)*in[1] + (l1+l2+l2)*in[2];
		} else { // bulkBoundX1 == true && bulkBoundY1 == true
			// bulkerror
		}
	} else if (z==MY) {
		if (!bulkBoundX1 && boundX1 != 1) {
			z = boundX1*MY;
		} else if (!bulkBoundY2 && boundY2 != MY) {
			z = boundY2;
		} else if (boundX1 == 1 && boundY2 == MY) {
			return l0*in[z] + l1y*in[z-1]+l1x*in[z+MY] + l2*in[z+MY-1];
		} else if (boundX1 == 1) {
			return (l0+l1+l1)*in[MY] + (l1+l2+l2)*in[MY+MY];
		} else if (boundY2 == MY) {
			return (l0+l1+l1)*in[MY] + (l1+l2+l2)*in[MY-1];
		} else {
			// bulkerror
		}
	} else if (z==M+1-MY) {
		if (!bulkBoundX2 && boundX2 != MX) {
			z = 1+(boundX2-1)*MY;
		} else if (!bulkBoundY1 && boundY1 != 1) {
			z = M-MY+boundY1;
		} else if (boundX2 == MX && boundY1 == 1) {
			return l0*in[z] + l1x*in[z-MY]+l1y*in[z+1] + l2*in[z-MY+1];
		} else if (boundX2 == MX) {
			return (l0+l1+l1)*in[z] + (l1+l2+l2)*in[z-MY];
		} else if (boundY1 == 1) {
			return (l0+l1+l1)*in[z] + (l1+l2+l2)*in[z+1];
		} else {
			// bulkerror
		}
	} else if (z==M) {
		if (!bulkBoundX2 && boundX2 != MX) {
			z = boundX2*MY;
		} else if (!bulkBoundY2 && boundY2 != MY) {
			z = M-MY+boundY2;
		} else if (boundX2 == MX && boundY2 == MY) {
			return l0*in[z] + l1y*in[z-1]+l1x*in[z-MY] + l2*in[z-MY-1];
		} else if (boundX2 == MX) {
			return (l0+l1+l1)*in[z] + (l1+l2+l2)*in[z-MY];
		} else if (boundY2 == MY) {
			return (l0+l1+l1)*in[z] + (l1+l2+l2)*in[z-1];
		} else {
			// bulkerror
		}
	}

	// If z is in boundary...
	// Either return a value or give z a value not in a boundary.
	if (z%MY == 1) { // If in left boundary...
		if (bulkBoundY1) { // This goes wrong in case of bulkerror.
			return (l0+l1+l1)*in[z] + (l1+l2+l2)*(in[z-MY]+in[z+MY]);
		} else if (boundY1 == 1) {
			return l0*in[z] + l1x*(in[z-MY]+in[z+MY])+l1y*in[z+1]
				+ l2*(in[z-MY+1]+in[z+MY+1]);
		} else { // If boundary is periodical or mirror...
			z += boundY1-1;
		}
	} else if (z%MY == 0) { // If in right boundary...
		if (bulkBoundY2) {
			return (l0+l1+l1)*in[z] + (l1+l2+l2)*(in[z-MY]+in[z+MY]);
		} else if (boundY2 == MY) {
			return l0*in[z] + l1x*(in[z-MY]+in[z+MY])+l1y*in[z-1]
				+ l2*(in[z-MY-1]+in[z+MY-1]);
		} else {
			z += boundY2-MY;
		}
	} else if (z <= MY) { // If in lower boundary...
		if (bulkBoundX1) {
			return (l0+l1+l1)*in[z] + (l1+l2+l2)*(in[z-1]+in[z+1]);
		} else if (boundX1 == 1) {
			return l0*in[z] + l1y*(in[z-1]+in[z+1]) +l1x*in[z+MY]
				+ l2*(in[z+MY-1]+in[z+MY+1]);
		} else {
			z += (boundX1-1)*MY;
		}
	} else if (z > M-MY) { // If in upper boundary...
		if (bulkBoundX2) {
			return (l0+l1+l1)*in[z] + (l1+l2+l2)*(in[z-1]+in[z+1]);
		} else if (boundX2 == MX) {
			return l0*in[z] + l1y*(in[z-1]+in[z+1])+l1x*in[z-MY]
				+ l2*(in[z-MY-1]+in[z-MY+1]);
		} else {
			z += boundX2*MY - M;
		}
	}

	// Now z is not in a boundary, we can calculate the sidefraction.
	return l0*in[z] + l1x*(in[z-MY] + in[z+MY])+ l1y*(in[z-1] + in[z+1])
		+ l2*(in[z-MY-1] + in[z-MY+1] + in[z+MY-1] + in[z+MY+1]);
}

void
Lat2DFlat:: Sides(Vector P, Vector Q) const{
	Message(fatal,"Sides not implemented");
}

void
Lat2DFlat::ElectricPotential(Vector &psi,
		Vector &psi0,
		const Vector &eps,
		const Vector &charge,
		const double preFactor) const {
	int x,y,z;
	double epsxmin, epsxplus, epsymin, epsyplus;
	double pf = 2*preFactor;

	// bulkboundaries are not (yet) implemented!
	// apply boundary conditions
	SetBoundaries(psi0);

	// lower x with its corners
	if (boundX1 == 1) {
		epsyplus = eps[1]+eps[2];
		if (boundY1 == 1) {
			epsxplus = eps[1]+eps[MY+1];
			psi[1] = (pf*charge[1] + epsxplus*psi0[MY+1]
				+ epsyplus*psi0[2])/(epsxplus + epsyplus);
		}
		for (z=MY-1; z>1; z--) {
			epsymin = epsyplus;
			epsxplus=epsyplus=eps[z];
			epsxplus += eps[z+MY];
			epsyplus += eps[z+1];
			psi[z] = (pf*charge[z] + epsxplus*psi0[z+MY]
				+ epsymin*psi0[z-1] + epsyplus*psi0[z+1])
				/(epsxplus + epsymin + epsyplus);
		}
		if (boundY2 == MY) {
			epsxplus = eps[MY]+eps[MY+MY];
			psi[MY] = (pf*charge[z] + epsxplus*psi0[MY+MY]
				+ epsyplus*psi0[MY-1])/(epsxplus + epsyplus);
		}
	}
	// set phi in middle of lattice and at the y edges
	z = MY+1;
	for (x=2; x<MX; x++) {
		epsyplus = eps[z]+eps[z+1];
		if (boundY1 == 1) {
			epsxmin = epsxplus = eps[z];
			epsxmin += eps[z-MY];
			epsxplus += eps[z+MY];
			psi[z] = (pf*charge[z] + epsxmin*psi0[z-MY] + epsxplus*psi0[z+MY]
				+ epsyplus*psi0[z+1])/(epsxmin + epsxplus + epsyplus);
		}
		z++;
		for (y=2; y<MY; y++) {
			epsymin = epsyplus;
			epsxmin=epsxplus=epsyplus=eps[z];
			epsxmin += eps[z-MY];
			epsxplus += eps[z+MY];
			epsyplus += eps[z+1];
			psi[z] = (pf*charge[z] + epsxmin*psi0[z-MY] + epsxplus*psi0[z+MY]
				+ epsymin*psi0[z-1] + epsyplus*psi0[z+1])
				/(epsxmin + epsxplus + epsymin + epsyplus);
			z++;
		}
		if (boundY2 == MY) {
			epsxmin=epsxplus=eps[z];
			epsxmin += eps[z-MY];
			epsxplus += eps[z+MY];
			psi[z] = (pf*charge[z] + epsxmin*psi0[z-MY] + epsxplus*psi0[z+MY]
				+ epsyplus*psi0[z-1])/(epsxmin + epsxplus + epsyplus);
		}
		z++;
	}
	// upper x with its corners
	if (boundX2 == MX) {
		epsyplus = eps[M-MY+1]+eps[M-MY+2];
		if (boundY1 == 1) {
			epsxmin = eps[M-MY+1]+eps[M-MY-MY+1];
			psi[M-MY+1] = (pf*charge[M-MY+1] + epsxmin*psi0[M-MY-MY+1]
				+ epsyplus*psi0[M-MY+2])/(epsxmin + epsyplus);
		}
		for (z=M-MY+2; z<M; z++) {
			epsymin = epsyplus;
			epsxmin=epsyplus=eps[z];
			epsxmin += eps[z-MY];
			epsyplus += eps[z+1];
			psi[z] = (pf*charge[z] + epsxmin*psi0[z-MY]
				+ epsymin*psi0[z-1] + epsyplus*psi0[z+1])
				/(epsxmin + epsymin + epsyplus);
		}
		if (boundY2 == MY) {
			epsxmin = eps[M-MY]+eps[M];
			psi[M] = (pf*charge[M] + epsxmin*psi0[M-MY]
				+ epsyplus*psi0[M-1])/(epsxmin + epsyplus);
		}
	}
}
void
Lat2DFlat::ElectricFieldSquared(Vector &ESquared,
		const Vector psi,
		const double preFactor) const {
	int x, y, z;
	double pf, Exmin, Explus, Eymin, Eyplus;

	pf = preFactor/2;
	if (boundX1 == 1) {
		Eymin = psi[MY-1]-psi[MY];
		Eymin *= Eymin;
		for (z=MY-1; z>1; z--) {
			Eyplus = Eymin;
			Eymin = psi[z]-psi[z-1];
			Eymin *= Eymin;
			Explus = psi[z]-psi[z+MY];
			Explus *= Explus;
			ESquared[z] = pf*(Explus + Eymin + Eyplus);
		}
	}
	z = MY+1;
	for (x=2; x<MX; x++) {
		Eyplus = psi[z]-psi[z+1];
		Eyplus *= Eyplus;
		if (boundY1 == 1) {
			Exmin = psi[z]-psi[z-MY];
			Exmin *= Exmin;
			Explus = psi[z]-psi[z+MY];
			Explus *= Explus;
			ESquared[z] = pf*(Exmin + Explus + Eyplus);
		}
		z++;
		for (y=2; y<MY; y++) {
			Eymin = Eyplus;
			Eyplus = psi[z]-psi[z+1];
			Eyplus *= Eyplus;
			Exmin = psi[z]-psi[z-MY];
			Exmin *= Exmin;
			Explus = psi[z]-psi[z+MY];
			Explus *= Explus;
			ESquared[z] = pf*(Exmin + Explus + Eymin + Eyplus);
			z++;
		}
		if (boundY2 == MY) {
			Eymin = Eyplus;
			Exmin = psi[z]-psi[z-MY];
			Exmin *= Exmin;
			Explus = psi[z]-psi[z+MY];
			Explus *= Explus;
			ESquared[z] = pf*(Exmin + Explus + Eymin);
		}
		z++;
	}
	if (boundX2 == MX) {
		Eyplus = psi[M-MY+2]-psi[M-MY+1];
		Eyplus *= Eyplus;
		for (z=M-MY+2; z<M; z++) {
			Eymin = Eyplus;
			Eyplus = psi[z]-psi[z+1];
			Eyplus *= Eyplus;
			Exmin = psi[z]-psi[z-MY];
			Exmin *= Exmin;
			ESquared[z] = pf*(Exmin + Eymin + Eyplus);
		}
	}
}
Vector
Lat2DFlat::Div1Grad2(const Vector in1, const Vector in2) const {
	Vector value(1,M);
	for (int z=1; z<=M; z++) {
		if (z%MY != 1 && in1[z] != 0 && in1[z-1] != 0) {
			value[z] += (in1[z-1]+in1[z])*(in2[z-1]-in2[z]);
		}
		if (z%MY != 0 && in1[z] != 0 && in1[z+1] != 0) {
			value[z] += (in1[z+1]+in1[z])*(in2[z+1]-in2[z]);
		}
		if (z > MY && in1[z] != 0 && in1[z-MY] != 0) {
			value[z] += (in1[z-MY]+in1[z])*(in2[z-MY]-in2[z]);
		}
		if (z < M - MY && in1[z] != 0 && in1[z+MY] != 0) {
			value[z] += (in1[z+MY]+in1[z])*(in2[z+MY]-in2[z]);
		}
		value[z] /= 4*siteDistance*siteDistance;
	}
	return value;
}
Vector
Lat2DFlat::FourierTransform(const Vector /*in*/) const {
	Vector value(1,M);
	Message(implementation,"Lat2DFlat::FourierTransform");
	return value;
}
double
Lat2DFlat::Moment(const Vector,
				  const int,
				  const LatticeRange*,
				  const double) const {
	double value = 0;
	// als meer dan 1 coordinate in Latrange : error
	Message(implementation,"Lat2DFlat::Moment");
	return value;
}
double
Lat2DFlat::MomentUnweighted(const Vector,
							const int,
							const LatticeRange*,
							const double) const {
	double value = 0;
	// als meer dan 1 coordinate in Latrange : error
	Message(implementation,"Lat2DFlat::MomentUnweighted");
	return value;
}
Vector
Lat2DFlat::RenormPhi(Vector a,
					 Vector,
					 const Vector,
					 const  double) const {
	Message(fatal,"Cannot have a molecule with freedom 'thirdGeneration' "
		"on a flat lattice with two gradients");
	return a; //never get here
}
double
Lat2DFlat::GetLayerAdjustment(void) const {
	return 1;
}
void
Lat2DFlat::SetLayerAdjustment(double) {
	Message(fatal,"don't know how to eliminate lattice artifact in 2 gradients");
}

void Lat2DFlat::SetBy1(double *P)const {
	int jx=MY;
	for (int x=0; x<MX; x++) P[jx*x]=P[jx*x+boundY1-1];
}

void Lat2DFlat::SetBy1(double *P,double *Q)const {
	int jx=MY;

	if (boundY1>2) {
		SetBy1(P); SetBy1(Q);
	}
	else {
		for (int x=0; x<MX; x++) {
			P[jx*x]=Q[jx*x+boundY1-1];
			Q[jx*x]=P[jx*x+boundY1-1];
		}
	}
}

void Lat2DFlat::SetBym(double *P) const {
	int jx=MY;
	int mmy = MY-1;
	for (int x=0; x<MX; x++) P[jx*x+mmy]=P[jx*x+boundY2-1];
}

void Lat2DFlat::SetBym(double *P,double *Q) const {
	int mmy = MY-1;
	int jx=MY;
	if (boundY2<numLayersY-1) {
		SetBym(P); SetBym(Q);
	}
	else {
		for (int x=0; x<MX; x++) {
			P[jx*x+mmy]=Q[jx*x+boundY2-1];
			Q[jx*x+mmy]=P[jx*x+boundY2-1];
		}
	}
}

void Lat2DFlat::SetBx1(double *P) const {
	int jx=MY;
	for (int y=0; y<MY; y++) P[y]=P[(boundX1-1)*jx+y];
}

void Lat2DFlat::SetBx1(double *P,double *Q) const {
	int jx=MY;
	if (boundX1>2) {
		SetBx1(P); SetBx1(Q);
	}
	else {
		for (int y=0; y<MY; y++) {
			P[y]=Q[(boundX1-1)*jx+y];
			Q[y]=P[(boundX1-1)*jx+y];
		}
	}
}

void Lat2DFlat::SetBxm(double *P) const {
	int mmx=MX-1;
	int jx=MY;
	for (int y=0; y<MY; y++) P[mmx*jx+y]=P[(boundX2-1)*jx+y];
}

void Lat2DFlat::SetBxm(double *P,double *Q) const {
	int mmx=MX-1;
	int jx=MY;
	if (boundX2<numLayersX-1) {
		SetBxm(P); SetBxm(Q);
	}
	else {
		for (int y=0; y<MY; y++) {
			P[mmx*jx+y]=Q[(boundX2-1)*jx+y];
			Q[mmx*jx+y]=P[(boundX2-1)*jx+y];
		}
	}
}

void Lat2DFlat::removeboundaries(Vector P) const {double *pP=&P[1]; removeboundaries(pP); }
void Lat2DFlat::removeboundaries(double *P) const {
	int jx=MY; int mmx=MX-1; int mmy=MY-1;
	for (int x=0; x<=mmx; x++) P[x*jx]=P[x*jx+mmy]=0;
	for (int y=0; y<=mmy; y++) for (int k=0; k<5; k++) P[y]=P[mmx*jx+y]=0;
}