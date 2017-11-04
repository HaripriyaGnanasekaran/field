#include "Lat3D.h"
#include "tools.h"
#include <iostream>
// lattice keyword
const Text Lat3D::lat="lat";

// allowed keywords
const Text Lat3D::dimensionsS="dimensions";
const Text Lat3D::gradients="gradients";
const Text Lat3D::geometry="geometry";
const Text Lat3D::n_layers_x="n_layers_x";
const Text Lat3D::n_layers_y="n_layers_y";
const Text Lat3D::n_layers_z="n_layers_z";
const Text Lat3D::lowerbound_x="lowerbound_x";
const Text Lat3D::upperbound_x="upperbound_x";
const Text Lat3D::lowerbound_y="lowerbound_y";
const Text Lat3D::upperbound_y="upperbound_y";
const Text Lat3D::lowerbound_z="lowerbound_z";
const Text Lat3D::upperbound_z="upperbound_z";
const Text Lat3D::distance="distance";
const Text Lat3D::bondlength="bondlength";
const Text Lat3D::lambda="lambda";
const Text Lat3D::lambda0="lambda0";
const Text Lat3D::lambda1="lambda1";
const Text Lat3D::lambda2="lambda2";
const Text Lat3D::lambda3="lambda3";
const Text Lat3D::latticetype="latticetype";
const Text Lat3D::stiff_range="stiff_range";
const Text Lat3D::number_of_computation_ranges="number_of_computation_ranges";

// choice options
const Text Lat3D::standard="standard";
const Text Lat3D::stencils="stencils";
const Text Lat3D::FCC="FCC";
const Text Lat3D::HEX="HEX";

// default values
const double Lat3D::l1stand=1.0/6;
const double Lat3D::l0sten=1.0/27;
const double Lat3D::l1sten=1.0/27;
const double Lat3D::l2sten=1.0/27;
const double Lat3D::l3sten=1.0/27;

Lat3D::Lat3D(Input* MyInput_, Text name_) {
	MyInput = MyInput_;
	name = name_;
	N_comp_ranges=0;
	if (*Copy(MyInput->GetText(lat,name,geometry)) == *Copy("flat")) {
		if (MyInput->ValueSet("lat",name,"number_of_computation_ranges")) {
			N_comp_ranges = MyInput ->GetInt("lat",name,"number_of_computation_ranges",2,99);
		};
		checkParameters();
	}

	numLayersX = MyInput->GetInt(lat,name,n_layers_x,1,INT_MAX);
	numLayersY = MyInput->GetInt(lat,name,n_layers_y,1,INT_MAX);
	numLayersZ = MyInput->GetInt(lat,name,n_layers_z,1,INT_MAX);

	MX = numLayersX + 2;
	mmx = MX-1;
	MY = numLayersY + 2;
	mmy = MY-1;
	MZ = numLayersZ + 2;
	mmz = MZ-1;
	M = MX*MY*MZ; // by setting M the arrays --when there is just one comp range--- will be dimensioned from 1 to M. I will use them from 0 to M-1
	jx = MY*MZ;
	jy = MZ;
	jz = 1;
	Array<Text> bounds(1,4);
	bounds[1] = "mirror1";
	bounds[2] = "mirror2";
	bounds[3] = "periodic";
	bounds[4] = "bulk";
	boundX1 = MyInput->GetChoice(lat, name, lowerbound_x, bounds, 3); // default periodic
	boundX2 = MyInput->GetChoice(lat, name, upperbound_x, bounds, 3); // default periodic
	boundY1 = MyInput->GetChoice(lat, name, lowerbound_y, bounds, 3); // default periodic
	boundY2 = MyInput->GetChoice(lat, name, upperbound_y, bounds, 3); // default periodic
	boundZ1 = MyInput->GetChoice(lat, name, lowerbound_z, bounds, 3); // default periodic
	boundZ2 = MyInput->GetChoice(lat, name, upperbound_z, bounds, 3); // default periodic
	int numBulkBounds = 0;
	if (boundX1 == 4) {
		boundX1 = 1;
		bulkBoundX1 = true;
		numBulkBounds++;
	} else {
		bulkBoundX1 = false;
	}
	if (boundX2 == 4) {
		boundX2 = 1; //set real value later
		bulkBoundX2 = true;
		numBulkBounds++;
	} else {
		bulkBoundX2 = false;
	}
	if (boundY1 == 4) {
		boundY1 = 1;
		bulkBoundY1 = true;
		numBulkBounds++;
	} else {
		bulkBoundY1 = false;
	}
	if (boundY2 == 4) {
		boundY2 = 1; //set real value later
		bulkBoundY1 = true;
		numBulkBounds++;
	} else {
		bulkBoundY2 = false;
	}
	if (boundZ1 == 4) {
		boundZ1 = 1;
		bulkBoundZ1 = true;
		numBulkBounds++;
	} else {
		bulkBoundZ1 = false;
	}
	if (boundZ2 == 4) {
		boundZ2 = 1; //set real value later
		bulkBoundZ1 = true;
		numBulkBounds++;
	} else {
		bulkBoundZ2 = false;
	}

	if ((bulkBoundX1 || bulkBoundX2) && (bulkBoundY1 || bulkBoundY2)&& (bulkBoundZ1 || bulkBoundZ2)) {
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
	if (bulkBoundZ1) {
		BulkBoundaries[++numBulkBounds] = Copy(lowerbound_z);
	}
	if (bulkBoundZ2) {
		BulkBoundaries[++numBulkBounds] = Copy(upperbound_z);
	}
	if (boundX1 == 3 && boundX2 != 3) {
		if (MyInput->ValueSet(lat, name, upperbound_x)) {
			Message(fatal,MyInput,"In 'lat : " + name
				+ "' lowerbound_x is periodic, while upperbound_x is not, "
				"this is not possible.");
		} else {
			boundX2 = 3;
		}
	}
	if (boundX2 == 3 && boundX1 != 3) {
		if (MyInput->ValueSet(lat, name, lowerbound_x)) {
			Message(fatal,MyInput,"In 'lat : " + name + "' upperbound_x is "
				"periodic, while lowerbound_x is not, this is not allowed.");
		} else {
			boundX1 = 3;
		}
	}
	if (boundY1 == 3 && boundY2 != 3) {
		if (MyInput->ValueSet(lat, name, upperbound_y)) {
			Message(fatal,MyInput,"In 'lat : " + name + "' lowerbound_y is "
				"periodic, while upperbound_y is not, this is not allowed.");
		} else {
			boundY2 = 3;
		}
	}
	if (boundY2 == 3 && boundY1 != 3) {
		if (MyInput->ValueSet(lat, name, lowerbound_y)) {
			Message(fatal,MyInput,"In 'lat : " + name + "' upperbound_y is "
				"periodic, while lowerbound_y is not, this is not allowed.");
		} else {
			boundY1 = 3;
		}
	}
	if (boundZ1 == 3 && boundZ2 != 3) {
		if (MyInput->ValueSet(lat, name, upperbound_z)) {
			Message(fatal,MyInput,"In 'lat : " + name + "' lowerbound_z is "
				"periodic, while upperbound_z is not, this is not allowed.");
		} else {
			boundZ2 = 3;
		}
	}
	if (boundZ2 == 3 && boundZ1 != 3) {
		if (MyInput->ValueSet(lat, name, lowerbound_z)) {
			Message(fatal,MyInput,"In 'lat : " + name + "' upperbound_z is "
				"periodic, while lowerbound_z is not, this is not allowed.");
		} else {
			boundZ1 = 3;
		}
	}


	if (boundX1 == 3) { boundX1 = numLayersX; boundX2 = 1;}
	else boundX2 = numLayersX + 1 - boundX2;
	if (boundY1 == 3) { boundY1 = numLayersY; boundY2 = 1;}
	else boundY2 = numLayersY + 1 - boundY2;
	if (boundZ1 == 3) { boundZ1 = numLayersZ; boundZ2 = 1;}
	else boundZ2 = numLayersZ + 1 - boundZ2;

	dimensions = MyInput->GetReal(lat, name, dimensionsS, 2, DBL_MAX,3);

	int lattyp = 1; // default standard lattice
	if (MyInput->ValueSet(lat,name,latticetype)) {
		Array<Text> lattype(1,4);
		lattype[1] = standard;
		lattype[2] = stencils;
		lattype[3] = FCC;
		lattype[4] = HEX;
		lattyp = MyInput->GetChoice(lat, name, latticetype, lattype, 1);
	}
	MyInput->DontCombineParam(lat,name,lambda,lambda0);
	MyInput->DontCombineParam(lat,name,lambda,lambda1);
	MyInput->DontCombineParam(lat,name,lambda,lambda2);
	MyInput->DontCombineParam(lat,name,lambda,lambda3);
	MyInput->AlwaysCombineParam(lat,name,lambda0,lambda1);
	MyInput->AlwaysCombineParam(lat,name,lambda0,lambda2);
	MyInput->AlwaysCombineParam(lat,name,lambda0,lambda3);
	MyInput->AlwaysCombineParam(lat,name,lambda1,lambda2);
	MyInput->AlwaysCombineParam(lat,name,lambda1,lambda3);
	MyInput->AlwaysCombineParam(lat,name,lambda2,lambda3);
	if (lattyp == 1) {
		latType = "BCC";
		l1 = 1.0/6.0;
		l0 = 0;
		l2 = 0;
		l3 = 0;
	} else if (lattyp == 3 ) {
		latType = FCC;
		l1 = 1.0/27.0;
		l0 = 1.0/27.0;
		l2 = 1.0/27.0;
		l3 = 1.0/27.0;
	} else if (lattyp == 4) {
		latType = HEX;
	} else {
		latType = stencils;
		l0 = MyInput->GetReal(lat,name,lambda0,0,0.5,l0sten);
		l1 = MyInput->GetReal(lat,name,lambda1,0,0.5,l1sten);
		l2 = MyInput->GetReal(lat,name,lambda2,0,0.5,l2sten);
		l3 = MyInput->GetReal(lat,name,lambda3,0,0.5,l3sten);
		// check whether total of lambdas is ~1
		// i1 will be smallest computer number > 1
		double i1 = l0 + 6.0*l1 + 12.0*l2 + 8.0*l3 - 1.0;
		i1 = (i1 < 0) ? -i1 : i1;
		if (i1 > DBL_EPSILON) {
			// error: total of lambdas in not ~1
			Message(fatal, MyInput, "Sum of lambda's is not 1.");
		}
		// make sure the sum of the lambdas is exactly 1.
		l3 = (1.0 - l0- 6.0*l1 - 12.0*l2)/8.0;
	}
	bondLength = MyInput->GetReal(lat, name, bondlength, 0, DBL_MAX, 3e-10);
	siteDistance = bondLength;

	if (MyInput->ValueSet("lat",name,"stiff_range")) {
		Text range = MyInput->GetText("lat",name,"stiff_range");
		stiff = NewLatticeRange(range);
	}

	if (N_comp_ranges==0) N_comp_ranges++;
	CR_Q.Dim(0,N_comp_ranges);
	for (int i=1; i<=N_comp_ranges; i++) {
		Text t,range;

		Mmax=0;
		if (i<10) {t= Blanks(1); t.Putint(i);}
		if (i>9 && i<100) {t= Blanks(2); t.Putint(i);}
		if (MyInput -> ValueSet("lat",name,"computation_range_"+t)) {
			range=MyInput->GetText("lat",name,"computation_range_"+t);
			CR_Q[i] = dynamic_cast<LatticeRange3D *>(NewLatticeRange(range));
			if (CR_Q[i]==0) Message(fatal, "Not a 3D lattice range");
			if (CR_Q[i]->Get_M()>Mmax) Mmax=CR_Q[i]->Get_M();

			int x0=CR_Q[i]->Get_x0();
			int y0=CR_Q[i]->Get_y0();
			int z0=CR_Q[i]->Get_z0();
			//int xm=CR_Q[i]->Get_xm();
			int ym=CR_Q[i]->Get_ym();
			int zm=CR_Q[i]->Get_zm();

			if (i>1) {
				int x1,y1,z1;
				if (boundX1!=1 || boundX2 !=numLayersX) {
					Message(fatal,"can not combine boundary condition in x-direction with nultiple computation ranges");
				}
				x1=CR_Q[i-1]->Get_xm();
				if (x0-x1 != 1){
					Message(fatal,"computation ranges must be adjacent in the x-direction. computation_range_" + t);
				}
				y1=CR_Q[i-1]->Get_ym();
				z1=CR_Q[i-1]->Get_zm();
				if (y0>y1 || z0>z1){
					Message(fatal,"computation range does not touch previous one in computation_range_" + t);
				}
				y1=CR_Q[i-1]->Get_y0();
				z1=CR_Q[i-1]->Get_z0();
				if (y1>ym || z1>zm) {
					Message(fatal,"computation range does not touch previous one in computation_range_" + t);
				}
			}
		} else {
			if (N_comp_ranges==1) {
				Text t1 = Blanks(10);
				Text t2 = Blanks(10);
				Text t3 = Blanks(10);
				Text t4 = Blanks(10);
				Text t5 = Blanks(10);
				Text t6 = Blanks(10);
				t1.Putint(1);
				t1 = Copy(t1.Frontstrip());
				t2.Putint(1);
				t2 = Copy(t2.Frontstrip());
				t3.Putint(1);
				t3 = Copy(t3.Frontstrip());

				t4.Putint(numLayersX);
				t4 = Copy(t4.Frontstrip());
				t5.Putint(numLayersY);
				t5 = Copy(t5.Frontstrip());
				t6.Putint(numLayersZ);
				t6 = Copy(t6.Frontstrip());

				Text t = t1 + "," + t2 + "," + t3 + ";" + t4 + "," + t5 + "," + t6;
				CR_Q[1]=dynamic_cast<LatticeRange3D *>(NewLatticeRange(t));
				if (CR_Q[1]==0) Message(fatal, "Not a 3D lattice range");
				Mmax=CR_Q[1]->Get_M();
			} else
				Message(fatal,"can't find latice range for computations cor computation_range_1");
		}
	}
}

Lat3D::~Lat3D() {
}

// protected function
void Lat3D::checkParameters() const {
	int num_params = 21+N_comp_ranges;

	Array<Text> param(1,num_params);
	param[1] = dimensionsS;
	param[2] = gradients;
	param[3] = geometry;
	param[4] = n_layers_x;
	param[5] = n_layers_y;
	param[6] = n_layers_z;
	param[7] = lowerbound_x;
	param[8] = upperbound_x;
	param[9] = lowerbound_y;
	param[10] = upperbound_y;
	param[11] = lowerbound_z;
	param[12] = upperbound_z;
	param[13] = bondlength;
	param[14] = lambda;
	param[15] = lambda0;
	param[16] = lambda1;
	param[17] = lambda2;
	param[18] = lambda3;
	param[19] = latticetype;
	param[20] = stiff_range;
	param[21] = number_of_computation_ranges;
	for (int i=1; i<=N_comp_ranges; i++)
	{ Text t;
		if (N_comp_ranges <10) t= Blanks(1); else t=Blanks(2);
		t.Putint(i);
		param[21+i] = "computation_range_" + t;
	}
	MyInput->CheckParameterNames(lat,name,param);
}

Text
Lat3D::GetName() const {
	return name;
}

void
Lat3D::GetOutput(Output* Out) const {
	Out->PutInt(lat,name,"gradients",3);
	if (dimensions - int(dimensions) > 1e-7) {
		Out->PutReal(lat,name,"dimensions",dimensions);
	} else {
		Out->PutInt(lat,name,"dimensions",int(dimensions));
	}
	Out->PutText(lat,name,"geometry","flat");
	Out->PutInt(lat,name,"n_layers_x",numLayersX);
	Out->PutInt(lat,name,"n_layers_y",numLayersY);
	Out->PutInt(lat,name,"n_layers_z",numLayersZ);
	if (GetTotalNumLatticeSites() - int(GetTotalNumLatticeSites()) > 1e-7) {
		Out->PutReal(lat,name,"total number of lattice sites",GetTotalNumLatticeSites());
	} else {
		Out->PutInt(lat,name,"total number of lattice sites",int(GetTotalNumLatticeSites()));
	}
	Array<Text> bounds(1,4);
	bounds[1] = "mirror1";
	bounds[2] = "mirror2";
	bounds[3] = "periodic";
	bounds[4] = "bulk";
	MyInput->SetDefaultWarningsOff();
	int boundX1out = MyInput->GetChoice(lat, name, "lowerbound_x", bounds, 1); // default mirror1
	int boundX2out = MyInput->GetChoice(lat, name, "upperbound_x", bounds, 1); // default mirror1
	int boundY1out = MyInput->GetChoice(lat, name, "lowerbound_y", bounds, 1); // default mirror1
	int boundY2out = MyInput->GetChoice(lat, name, "upperbound_y", bounds, 1); // default mirror1
	int boundZ1out = MyInput->GetChoice(lat, name, "lowerbound_z", bounds, 1); // default mirror1
	int boundZ2out = MyInput->GetChoice(lat, name, "upperbound_z", bounds, 1); // default mirror1
	MyInput->SetDefaultWarningsOn();
	Out->PutText(lat,name,"lowerbound_x",bounds[boundX1out]);
	Out->PutText(lat,name,"upperbound_x",bounds[boundX2out]);
	Out->PutText(lat,name,"lowerbound_y",bounds[boundY1out]);
	Out->PutText(lat,name,"upperbound_y",bounds[boundY2out]);
	Out->PutText(lat,name,"lowerbound_z",bounds[boundZ1out]);
	Out->PutText(lat,name,"upperbound_z",bounds[boundZ2out]);

	if (MyInput->ValueSet(lat,name,lambda)) {
		Out->PutReal(lat,name,lambda,l1);
	} else {
		Out->PutReal(lat, name, lambda0, l0);
		Out->PutReal(lat, name, lambda1, l1);
		Out->PutReal(lat, name, lambda2, l2);
		Out->PutReal(lat, name, lambda3, l3);
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
	Out->PutReal(lat,name,"volume (m^3)",GetVolume());

}

int
Lat3D::GetTotalNumLayers() const {
	int Mtot=0;
	for (int i=1; i<=N_comp_ranges; i++) Mtot+=CR_Q[i]->Get_M();
	return Mtot;
}

int
Lat3D::Get_BL(int bn) const {
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
	case 5:
		return boundZ1;
		break;
	case 6:
		return boundZ2;
		break;
	default:
		Message(fatal,"Asking a boundary layer outsize range in Lat3D.");
		return 0;
		break;
	}
}

int
Lat3D::Get_N_comp_ranges() const {
	return N_comp_ranges;
}

int
Lat3D::GetR0(int co_or) const{
	co_or--;
	int comp_range = (co_or)/3 +1;
	if (co_or%3+1 ==1) return CR_Q[comp_range]->Get_x0();
	if (co_or%3+1 ==2) return CR_Q[comp_range]->Get_y0();
	if (co_or%3+1 ==3) return CR_Q[comp_range]->Get_z0();
	return 0;
}

int
Lat3D::GetRm(int co_or) const{
	co_or--;
	int comp_range = (co_or)/3 +1;
	if (co_or%3+1 ==1) return CR_Q[comp_range]->Get_xm();
	if (co_or%3+1 ==2) return CR_Q[comp_range]->Get_ym();
	if (co_or%3+1 ==3) return CR_Q[comp_range]->Get_zm();
	return 0;
}

int
Lat3D::GetNumGradients() const {
	return 3;
}
int
Lat3D::GetNumLayers(const int numGradient) const {
	//This will return the lattice dimensions not the computational ranges...
	if (numGradient == 1) return MX;
	if (numGradient == 2) return MY;
	if (numGradient == 3) return MZ;
	Message(fatal,"Argument out of range in call to 'Lat3D::GetNumLayers(int D)'");
	return -1; // never get here
}
double
Lat3D::GetNumLatticeSites(const int) const {
	return 1;
}
double
Lat3D::GetLambda(const int z1, const int z2) const {
	if (N_comp_ranges>1)
		Message(fatal,"GetLambda needed while number of computation ranges larger than unity. Ask Frans to correct this mistake");
	int diff = abs(z1 - z2);
	if (diff == 0) {return 0.0;
	} else {if (diff == jz || diff == jy || diff == jx) {
		return l1;
	} else {if (diff == jx+jy || diff == jx-jy || diff == jx+jz || diff == jx-jz || diff == jy+jz || diff == jy-jz ) {
		return l2;
	} else {if (diff == jx+jy+jz || diff == jx+jy-jz || diff == jx-jy+jz || diff == jx-jy-jz) {
		return l3;
	} else {return 0.0;
	}}}}}
void
Lat3D::MultiplyWithLatticeSites(Vector) const {
}
void
Lat3D::DivideByLatticeSites(Vector) const {
}
double
Lat3D::GetVolume() const {
	return GetTotalNumLatticeSites()*siteDistance*siteDistance*siteDistance;
}
double
Lat3D::GetTotalNumLatticeSites() const {
	if (N_comp_ranges == 1){
		double valueX = MX;

		if (boundX1 == 2) {valueX -= 1.5;} else { valueX -=1.0;}
		if (boundX2 == numLayersX-1) {valueX -= 1.5;} else { valueX -=1.0;}
		double valueY = MY;
		if (boundY1 == 2) {valueY -= 1.5;} else { valueY -=1.0;};
		if (boundY2 == numLayersY-1) {valueY -= 1.5;} else { valueY -=1.0;}
		double valueZ = MZ;
		if (boundZ1 == 2) {valueZ -= 1.5;} else { valueZ -=1.0;};
		if (boundZ2 == numLayersZ-1) {valueZ -= 1.5;} else { valueZ -=1.0;}
		return valueX*valueY*valueZ;
	} else {
		double Total=0;
		double valueX,valueY,valueZ;
		for (int i=1; i<=N_comp_ranges; i++){
			valueX=CR_Q[i]->GetNumLayers_x(); //boundary condition inherited from the lattice boundary conditions.
			valueY=CR_Q[i]->GetNumLayers_y()+2;
			valueZ=CR_Q[i]->GetNumLayers_z()+2;

			if (boundY1 == 2) {valueY -= 1.5;} else { valueY -=1.0;};
			if (boundY2 == numLayersY-1) {valueY -= 1.5;} else { valueY -=1.0;}
			if (boundZ1 == 2) {valueZ -= 1.5;} else { valueZ -=1.0;};
			if (boundZ2 == numLayersZ-1) {valueZ -= 1.5;} else { valueZ -=1.0;}
			Total += valueX*valueY*valueZ;
		}
		return Total;
	}
}

double
Lat3D::GetNumExtraLatticeSites() const {
	return 0;
}

double
Lat3D::GetSiteDistance() const {
	return siteDistance;
}
double
Lat3D::GetSiteSurface() const {
	return pow(siteDistance,dimensions-1);
}
LatticeRange*
Lat3D::NewLatticeRange(Text range) const {
	int i;
	Boolean error = false;
	int minX,maxX;
	int minY,maxY;
	int minZ,maxZ;
	Text xValue, yValue, zValue;
	char test;
	//mmm if the range is a file return a LatticeRangeFile object
	Text t = Copy(range);

	if (t.Getchar()=='<') {
		//Message(literal, " lat3cpp in if0 ");
		t = t.Rest();
		//Message(literal, " lat3cpp in if1 ");
		return new LatticeRangeFile(t.Frontstrip(),GetTotalNumLayers());
	}
	//mmm

	t.Setpos(1);

	if (t.Getchar()=='=') {
		t = t.Rest();
		return new LatticeRangeFile(t.Frontstrip(),numLayersX,numLayersY,numLayersZ);
	}

	if (*range == *Copy("free")) return new LatticeRange3D(0,mmx,mmx,0,mmy,mmy,0,mmz,mmz);
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
	if (!low.More()) error = true;
	zValue = Copy(low.Scanto(','));
	zValue = Copy(zValue.Strip());
	zValue = Copy(zValue.Frontstrip());
	if (low.More()) error = true;
	if (xValue.Length() == 0) error = true;
	if (yValue.Length() == 0) error = true;
	if (zValue.Length() == 0) error = true;
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

	if (*zValue == *Copy("lowerbound")) {
		minZ = 0;
	} else if (*zValue == *Copy("upperbound")) {
		minZ = MZ-1;
	} else if (*zValue == *Copy("lastlayer")) {
		minZ = MZ-2;
	} else {
		int length = zValue.Length();
		for (i=1; i<= length; i++) {
			test = zValue.Getchar();
			if (!isdigit(test)) {
				Message(fatal,MyInput,"Error reading range '" + range
				+ "', expecting an integer, 'lowerbound', 'upperbound' or 'lastlayer' for first z coordinate");
			}
		}
		zValue.Setpos(1);
		minZ = zValue.Getint();
		if (minZ < 0) {
			Message(fatal,MyInput,"Error reading range '" + range
			+ "', value cannot be negative");
		}
		if (minZ > numLayersZ+1) {
			Message(fatal,MyInput,"Error reading range '" + range
			+ "', value for first z coordinate exceeds n_layers_z+1");
		}
	}

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
	if (!high.More()) error = true;
	zValue = Copy(high.Scanto(','));
	zValue = Copy(zValue.Strip());
	zValue = Copy(zValue.Frontstrip());
	
	if (high.More()) error = true;
	if (xValue.Length() == 0) error = true;
	if (yValue.Length() == 0) error = true;
	if (zValue.Length() == 0) error = true;
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
	if (*zValue == *Copy("upperbound")) {
		maxZ = MZ-1;
	} else if (*zValue == *Copy("lowerbound")) {
		maxZ = 0;
	} else if (*zValue == *Copy("lastlayer")) {
		maxZ = MZ-2;
	} else {
		int length = zValue.Length();
		for (i=1; i<= length; i++) {
			test = zValue.Getchar();
			if (!isdigit(test)) {
				Message(fatal,MyInput,"Error reading range '" + range
				+ "', expecting an integer, 'lowerbound', 'upperbound' or 'lastlayer' for second z coordinate");
			}
		}
		zValue.Setpos(1);
		maxZ = zValue.Getint();
		if (maxZ < 0) {
			Message(fatal,MyInput,"Error reading range '" + range
			+ "', value cannot be negative");
		}
		if (maxZ > numLayersZ+1) {
			Message(fatal,MyInput,"Error reading range '" + range
			+ "', value for second z coordinate exceeds n_layers_z+1");
		}
	}
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
	if (minZ > maxZ) {
		Message(hint,MyInput,"In range '" + range
		+ "' the first z coordinate is higher than the last.\nThese values will be swapped");
		int dummy = minZ;
		minZ = maxZ;
		maxZ = dummy;
	}
	return new LatticeRange3D(minX, maxX, mmx, minY, maxY, mmy, minZ, maxZ, mmz);
}

Boolean
Lat3D::WithinBoundaries(const int position) const {
	int x,y,z,pos;
	pos=position-1;
	if (N_comp_ranges == 1){
		z = ((pos % jx) % jy)/jz;
		y = ((pos % jx)-z*jz)/jy;
		x = (pos-y*jy-z*jz)/jx;
		if (x < 1 || x > numLayersX || y < 1 || y > numLayersY || z < 1 || z > numLayersZ) {return false;} else {return true;}
	} else {
		int jx,jy;
		for (int i=1; i<=N_comp_ranges; i++){
			if (pos < CR_Q[i]->Get_M()) {
				jx=CR_Q[i]->Get_jx();
				jy=CR_Q[i]->Get_jy();
				z = ((pos % jx) % jy)/jz;
				y = ((pos % jx)-z*jz)/jy;
				x = (pos-y*jy-z*jz)/jx;
				if (x < 1 || x > CR_Q[i]->GetNumLayers_x() ||
					y < 1 || y > CR_Q[i]->GetNumLayers_y() ||
					z < 1 || z > CR_Q[i]->GetNumLayers_z()) {return false;} else {return true;}
			} else {
				pos-=CR_Q[i]->Get_M();
			}
		}
	}
	return false;
}

void Lat3D::SetBz1(Vector P) const {double *pP=&P[1]; SetBz1(pP);}
void Lat3D::SetBz1(double *P) const {
	int i,x,y,px=0,py;
	if (N_comp_ranges==1) {

		px=-jx;
		for (x=0; x<MX; x++)  {
			px += jx; py = -jy; for (y=0; y<MY; y++) {
				py += jy; P[px+py+0]=P[px+py+boundZ1];
			}
		}
	}
	else {
		//int Previous=0;
		int Current=0,Next=0;
		int NlayersX,NlayersY;
		//int NlayersZ;
		int Jx,Jy,BZ1;
		for (i=1; i<=N_comp_ranges; i++){
			//Previous=Current;
			Current=Next; Next +=CR_Q[i]->Get_M();
			NlayersX=CR_Q[i]->GetNumLayers_x()+1;
			NlayersY=CR_Q[i]->GetNumLayers_y()+1;
			//NlayersZ=CR_Q[i]->GetNumLayers_z()+1;
			Jx=CR_Q[i]->Get_jx();
			Jy=CR_Q[i]->Get_jy();

			BZ1=boundZ1;
			px=-Jx;
			for (x=0; x<=NlayersX; x++)  {
				px += Jx; py=-Jy; for (y=0; y<=NlayersY; y++) {
					py+=Jy;
					P[Current+px+py+0]=P[Current+px+py+BZ1];
				}
			}
		}
	}
}

void Lat3D::SetBz1(double *P,double *Q) const {
	int i,x,y,px=0,py;
	if (boundZ1>2) {
		SetBz1(P); SetBz1(Q);
	}
	else if (N_comp_ranges==1) {
		px=0;
		for (x=0; x<MX; x++)  {
			px += jx; py = -jy; for (y=0; y<MY; y++) {
				py += jy; P[px+py+0]=Q[px+py+boundZ1];
				Q[px+py+0]=P[px+py+boundZ1];
			}
		}
	}
	else {
		//int Previous=0;
		int Current=0,Next=0;
		//int NlayersX;
		int NlayersY;
		//int NlayersZ;
		int Jx,Jy,BZ1;
		for (i=1; i<=N_comp_ranges; i++) {
			//Previous=Current;
			Current=Next; Next +=CR_Q[i]->Get_M();
			//NlayersX=CR_Q[i]->GetNumLayers_x()+1;
			NlayersY=CR_Q[i]->GetNumLayers_y()+1;
			//NlayersZ=CR_Q[i]->GetNumLayers_z()+1;
			Jx=CR_Q[i]->Get_jx();
			Jy=CR_Q[i]->Get_jy();

			BZ1=boundZ1;
			px=-Jx;
			for (x=0; x<MX; x++)  {
				px += Jx; py=-Jy; for (y=0; y<=NlayersY; y++) {
					py+=Jy;
					P[Current+px+py+0]=Q[Current+px+py+BZ1];
					Q[Current+px+py+0]=P[Current+px+py+BZ1];
				}
			}
		}
	}
}

void Lat3D::SetBzm(Vector P) const {double *pP=&P[1]; SetBzm(pP);}
void Lat3D::SetBzm(double *P) const {
	int i,x,y,px=0,py;
	if (N_comp_ranges==1) {

		px=-jx;
		for (x=0; x<MX; x++){
			px += jx; py = -jy;  for (y=0; y<MY; y++){
				py += jy; P[px+py+mmz]=P[px+py+boundZ2];
			}
		}
	}
	else {
		//int Previous=0;
		int Current=0,Next=0;
		int NlayersX;
		int NlayersY,NlayersZ,Jx,Jy,BZ2;
		for (i=1; i<=N_comp_ranges; i++) {
			//Previous=Current;
			Current=Next; Next +=CR_Q[i]->Get_M();
			NlayersX=CR_Q[i]->GetNumLayers_x()+1;
			NlayersY=CR_Q[i]->GetNumLayers_y()+1;
			NlayersZ=CR_Q[i]->GetNumLayers_z()+1;
			Jx=CR_Q[i]->Get_jx();
			Jy=CR_Q[i]->Get_jy();

			BZ2=boundZ2;
			px=-Jx;
			for (x=0; x<=NlayersX; x++)  {
				px += Jx; py=-Jy; for (y=0; y<=NlayersY; y++) {
					py+=Jy; P[Current+px+py+(NlayersZ)]=P[Current+px+py+BZ2];
				}
			}
		}
	}
}

void Lat3D::SetBzm(double *P,double *Q) const {
	int i,x,y,px=0,py;
	if (boundZ2<numLayersZ-2) {
		SetBzm(P); SetBzm(Q);
	}
	else if (N_comp_ranges==1) {
		px=-jx;
		for (x=0; x<MX; x++){
			px += jx; py = -jy;  for (y=0; y<MY; y++){
				py += jy;
				P[px+py+mmz]=Q[px+py+boundZ2];
				Q[px+py+mmz]=P[px+py+boundZ2];
			}
		}
	}
	else {
		//int Previous=0;
		int Current=0,Next=0;
		int NlayersX,NlayersY,NlayersZ,Jx,Jy,BZ2;
		for (i=1; i<=N_comp_ranges; i++) {
			//Previous=Current;
			Current=Next; Next +=CR_Q[i]->Get_M();
			NlayersX=CR_Q[i]->GetNumLayers_x()+1;
			NlayersY=CR_Q[i]->GetNumLayers_y()+1;
			NlayersZ=CR_Q[i]->GetNumLayers_z()+1;
			Jx=CR_Q[i]->Get_jx();
			Jy=CR_Q[i]->Get_jy();

			BZ2=boundZ2;
			px=-Jx;
			for (x=0; x<=NlayersX; x++)  {
				px += Jx; py=-Jy; for (y=0; y<=NlayersY; y++) {
					py+=Jy;
					P[Current+px+py+(NlayersZ)]=Q[Current+px+py+BZ2];
					Q[Current+px+py+(NlayersZ)]=P[Current+px+py+BZ2];
				}
			}
		}
	}
}

void Lat3D::SetBy1(Vector P)const {SetBy1(&P[1]);}
void Lat3D::SetBy1(double *P)const {
	int i,x,z,px=0;
	if (N_comp_ranges==1) {
		px=-jx;
		for (x=0; x<MX; x++) 	{
			px += jx; for (z=0; z<MZ; z++)	P[px+0+z]=P[px+boundY1*jy+z];
		}
	}
	else {
		//int Previous=0;
		int Current=0,Next=0;
		int NlayersX;
		//int NlayersY;
		int NlayersZ;
		int Jx,Jy,BY1;
		for (i=1; i<=N_comp_ranges; i++) {
			//Previous=Current;
			Current=Next; Next +=CR_Q[i]->Get_M();
			NlayersX=CR_Q[i]->GetNumLayers_x()+1;
			//NlayersY=CR_Q[i]->GetNumLayers_y()+1;
			NlayersZ=CR_Q[i]->GetNumLayers_z()+1;
			Jx=CR_Q[i]->Get_jx();
			Jy=CR_Q[i]->Get_jy();
			//if (boundY1<3) BY1=boundY1; else BY1=NlayersY-2;
			BY1=boundY1;
			px=-Jx;
			for (x=0; x<=NlayersX; x++)  {
				px += Jx; for (z=0; z<=NlayersZ; z++) {
					P[Current+px+0+z]=P[Current+px+BY1*Jy+z];
				}
			}
		}
	}
}

void Lat3D::SetBy1(double *P,double *Q)const {
	int i,x,z,px=0;
	if (boundY1>2) {
		SetBy1(P); SetBy1(Q);
	}
	else if (N_comp_ranges==1) {
		px=-jx;
		for (x=0; x<MX; x++) 	{
			px += jx; for (z=0; z<MZ; z++) {
				P[px+0+z]=Q[px+boundY1*jy+z];
				Q[px+0+z]=P[px+boundY1*jy+z];
			}
		}
	}
	else {
		//int Previous=0;
		int Current=0,Next=0;
		int NlayersX;
		//int NlayersY;
		int NlayersZ;
		int Jx,Jy,BY1;
		for (i=1; i<=N_comp_ranges; i++) {
			//Previous=Current;
			Current=Next; Next +=CR_Q[i]->Get_M();
			NlayersX=CR_Q[i]->GetNumLayers_x()+1;
			//NlayersY=CR_Q[i]->GetNumLayers_y()+1;
			NlayersZ=CR_Q[i]->GetNumLayers_z()+1;
			Jx=CR_Q[i]->Get_jx();
			Jy=CR_Q[i]->Get_jy();
			//if (boundY1<3) BY1=boundY1; else BY1=NlayersY-2;
			BY1=boundY1;
			px=-Jx;
			for (x=0; x<=NlayersX; x++)  {
				px += Jx; for (z=0; z<=NlayersZ; z++) {
					P[Current+px+0+z]=Q[Current+px+BY1*Jy+z];
					Q[Current+px+0+z]=P[Current+px+BY1*Jy+z];
				}
			}
		}
	}
}

void Lat3D::SetBym(Vector P) const {SetBym(&P[1]);}
void Lat3D::SetBym(double *P) const {
	int i,x,z,px=0;
	//Message(debug,"SetBYm");
	if (N_comp_ranges==1) {
		px=-jx;
		for (x=0; x<MX; x++) {
			px += jx; for (z=0; z<MZ; z++)
				P[px+mmy*jy+z]=P[px+boundY2*jy+z];
		}
	}
	else {
		//int Previous=0;
		int Current=0,Next=0;
		int NlayersX,NlayersY,NlayersZ,Jx,Jy,BY2;
		for (i=1; i<=N_comp_ranges; i++) {
			//Previous=Current;
			Current=Next; Next +=CR_Q[i]->Get_M();
			NlayersX=CR_Q[i]->GetNumLayers_x()+1;
			NlayersY=CR_Q[i]->GetNumLayers_y()+1;
			NlayersZ=CR_Q[i]->GetNumLayers_z()+1;
			Jx=CR_Q[i]->Get_jx();
			Jy=CR_Q[i]->Get_jy();

			BY2=boundY2;
			px=-Jx;
			for (x=0; x<=NlayersX; x++)  {
				px += Jx; for (z=0; z<=NlayersZ; z++) {
					P[Current+px+(NlayersY)*Jy+z]=P[Current+px+BY2*Jy+z];
				}
			}
		}
	}
}
void Lat3D::SetBym(double *P,double *Q) const {
	int x,z,px=0;
	//Message(debug,"SetBYm PQ");
	if (boundY2<numLayersY-1) {
		SetBym(P); SetBym(Q);
	}
	else if (N_comp_ranges==1) {
		px=-jx;
		for (x=0; x<MX; x++) {
			px += jx; for (z=0; z<MZ; z++)	{
				P[px+mmy*jy+z]=Q[px+boundY2*jy+z];
				Q[px+mmy*jy+z]=P[px+boundY2*jy+z];
			}
		}
	}
	else {
		//int Previous=0;
		int Current=0,Next=0;
		int i,NlayersX,NlayersY,NlayersZ,Jx,Jy,BY2;
		for (i=1; i<=N_comp_ranges; i++) {
			//Previous=Current;
			Current=Next; Next +=CR_Q[i]->Get_M();
			NlayersX=CR_Q[i]->GetNumLayers_x()+1;
			NlayersY=CR_Q[i]->GetNumLayers_y()+1;
			NlayersZ=CR_Q[i]->GetNumLayers_z()+1;
			Jx=CR_Q[i]->Get_jx();
			Jy=CR_Q[i]->Get_jy();

			BY2=boundY2;
			px=-Jx;
			for (x=0; x<=NlayersX; x++)  {
				px += Jx; for (z=0; z<=NlayersZ; z++) {
					P[Current+px+(NlayersY)*Jy+z]=Q[Current+px+BY2*Jy+z];
					Q[Current+px+(NlayersY)*Jy+z]=P[Current+px+BY2*Jy+z];
				}
			}
		}
	}
}

void Lat3D::SetBx1(Vector P) const {SetBx1(&P[1]);}
void Lat3D::SetBx1(double *P) const {
	int i,y,z,py=0;
	//Message(debug,"SetBX1");
	if (N_comp_ranges==1) {
		py=-jy;
		for (y=0; y<MY; y++) {
			py +=jy; for (z=0; z<MZ; z++) P[0+py+z]=P[boundX1*jx+py+z];
		}
	}
	else {
		int Previous=0;
		int Current=0,Next=0;
		int JX,JY,X0,Y0,Z0,x0,y0,z0;
		//int xM,yM,zM;
		//int XM;
		int YM,ZM,BX1,Jx,Jy;
		//int NlayersX;
		int NlayersY;
		int NlayersZ;
		for (i=1; i<=N_comp_ranges; i++) {
			Previous=Current;
			Current=Next; Next +=CR_Q[i]->Get_M();
			//NlayersX=CR_Q[i]->GetNumLayers_x()+1;
			NlayersY=CR_Q[i]->GetNumLayers_y()+1;
			NlayersZ=CR_Q[i]->GetNumLayers_z()+1;
			if (i>1) {
				Jx=CR_Q[i]->Get_jx();
				Jy=CR_Q[i]->Get_jy();
				x0=CR_Q[i]->Get_x0();
				y0=CR_Q[i]->Get_y0();
				z0=CR_Q[i]->Get_z0();
				//xM=CR_Q[i]->Get_xm();
				//yM=CR_Q[i]->Get_ym();
				//zM=CR_Q[i]->Get_zm();
				X0=CR_Q[i-1]->Get_x0();
				Y0=CR_Q[i-1]->Get_y0();
				Z0=CR_Q[i-1]->Get_z0();
				//XM=CR_Q[i-1]->Get_xm();
				YM=CR_Q[i-1]->Get_ym();
				ZM=CR_Q[i-1]->Get_zm();
				JX=CR_Q[i-1]->Get_jx();
				JY=CR_Q[i-1]->Get_jy();
				py=-Jy;
				for (y=0; y<=NlayersY; y++) {
					py+=Jy; for (z=0; z<=NlayersZ; z++){
						if (y0+y-1 >= Y0 && y0+y-1<=YM && z0+z-1 >= Z0 && z0+z-1<=ZM) {
							P[Current+0+py+z]=P[Previous+(x0-X0)*JX+(y0+y-Y0)*JY+(z0+z-Z0)];
						}
						else {

							BX1=boundX1;
							P[Current+0+py+z]=P[Current+BX1*Jx+py+z];
						}
					}
				}
			}
			else {
				Jy=CR_Q[i]->Get_jy();
				Jx=CR_Q[i]->Get_jx();

				BX1=boundX1;
				py=-Jy;
				for (y=0; y<=NlayersY; y++) {
					py +=Jy; for (z=0; z<=NlayersZ; z++) {
						P[Previous+0+py+z]=P[Previous+BX1*Jx+py+z];
					}
				}
			}
		}
	}
}
void Lat3D::SetBx1(double *P,double *Q) const {
	int i,y,z,py=0;
	//Message(debug,"SetBX1 PQ");
	if (boundX1>2) {
		SetBx1(P); SetBx1(Q);
	}
	else if (N_comp_ranges==1) {
		py=-jy;
		for (y=0; y<MY; y++) {
			py +=jy; for (z=0; z<MZ; z++) {
				P[0+py+z]=Q[boundX1*jx+py+z];
				Q[0+py+z]=P[boundX1*jx+py+z];
			}
		}
	}
	else {
		int Previous=0,Current=0,Next=0;
		int JX,JY,X0,Y0,Z0,x0,y0,z0;
		//int xM,yM,zM;
		//int XM;
		int YM,ZM,BX1,Jx,Jy;
		//int NlayersX;
		int NlayersY;
		int NlayersZ;
		for (i=1; i<=N_comp_ranges; i++) {
			Previous=Current; Current=Next; Next +=CR_Q[i]->Get_M();
			//NlayersX=CR_Q[i]->GetNumLayers_x()+1;
			NlayersY=CR_Q[i]->GetNumLayers_y()+1;
			NlayersZ=CR_Q[i]->GetNumLayers_z()+1;
			if (i>1) {
				Jx=CR_Q[i]->Get_jx();
				Jy=CR_Q[i]->Get_jy();
				x0=CR_Q[i]->Get_x0();
				y0=CR_Q[i]->Get_y0();
				z0=CR_Q[i]->Get_z0();
				//xM=CR_Q[i]->Get_xm();
				//yM=CR_Q[i]->Get_ym();
				//zM=CR_Q[i]->Get_zm();
				X0=CR_Q[i-1]->Get_x0();
				Y0=CR_Q[i-1]->Get_y0();
				Z0=CR_Q[i-1]->Get_z0();
				//XM=CR_Q[i-1]->Get_xm();
				YM=CR_Q[i-1]->Get_ym();
				ZM=CR_Q[i-1]->Get_zm();
				JX=CR_Q[i-1]->Get_jx();
				JY=CR_Q[i-1]->Get_jy();
				py=-Jy;
				for (y=0; y<=NlayersY; y++) {
					py+=Jy; for (z=0; z<=NlayersZ; z++){
						if (y0+y-1 >= Y0 && y0+y-1<=YM && z0+z-1 >= Z0 && z0+z-1<=ZM) {
							P[Current+0+py+z]=P[Previous+(x0-X0)*JX+(y0+y-Y0)*JY+(z0+z-Z0)];
							Q[Current+0+py+z]=Q[Previous+(x0-X0)*JX+(y0+y-Y0)*JY+(z0+z-Z0)];
						}
						else {

							BX1=boundX1;
							P[Current+0+py+z]=Q[Current+BX1*Jx+py+z];
							Q[Current+0+py+z]=P[Current+BX1*Jx+py+z];
						}
					}
				}
			}
			else {
				Jy=CR_Q[i]->Get_jy();
				Jx=CR_Q[i]->Get_jx();

				BX1=boundX1;
				py=-Jy;
				for (y=0; y<=NlayersY; y++) {
					py +=Jy; for (z=0; z<=NlayersZ; z++) {
						P[Previous+0+py+z]=Q[Previous+BX1*Jx+py+z];
						Q[Previous+0+py+z]=P[Previous+BX1*Jx+py+z];
					}
				}
			}
		}
	}
}

void Lat3D::SetBxm(Vector P) const {SetBxm(&P[1]);}
void Lat3D::SetBxm(double *P) const {
	int i,y,z,py=0;
	//Message(debug,"SetBXm");
	if (N_comp_ranges==1) {
		py=-jy;
		for (y=0; y<MY; y++) {
			py +=jy; for (z=0; z<MZ; z++)
			P[mmx*jx+py+z]=P[boundX2*jx+py+z];
		}
	}
	else {
		//int Previous=0;
		int Current=0,Next=0;
		int JX=0,JY;
		//int X0;
		int Y0,Z0;
		//int x0;
		int y0,z0;
		//int xM,yM,zM;
		//int XM;
		int YM,ZM,BX2;
		int Jx=0,Jy=0,NlayersX,NlayersY,NlayersZ;
		for (i=1; i<=N_comp_ranges; i++) {
			//Previous=Current;
			Current=Next;  Next +=CR_Q[i]->Get_M();
			NlayersX=CR_Q[i]->GetNumLayers_x()+1;
			NlayersY=CR_Q[i]->GetNumLayers_y()+1;
			NlayersZ=CR_Q[i]->GetNumLayers_z()+1;
			if (i<N_comp_ranges) {
				Jx=CR_Q[i]->Get_jx();
				Jy=CR_Q[i]->Get_jy();
				//x0=CR_Q[i]->Get_x0();
				y0=CR_Q[i]->Get_y0();
				z0=CR_Q[i]->Get_z0();
				//xM=CR_Q[i]->Get_xm();
				//yM=CR_Q[i]->Get_ym();
				//zM=CR_Q[i]->Get_zm();
				//X0=CR_Q[i+1]->Get_x0();
				Y0=CR_Q[i+1]->Get_y0();
				Z0=CR_Q[i+1]->Get_z0();
				//XM=CR_Q[i+1]->Get_xm();
				YM=CR_Q[i+1]->Get_ym();
				ZM=CR_Q[i+1]->Get_zm();
				JX=CR_Q[i+1]->Get_jx();
				JY=CR_Q[i+1]->Get_jy();
				py=-Jy;
				for (y=0; y<=NlayersY; y++) {
					py+=Jy; for (z=0; z<=NlayersZ; z++){
						if (y0+y-1 >= Y0 && y0+y-1<=YM && z0+z-1 >= Z0 && z0+z-1<=ZM )  {
							P[Current+(NlayersX)*Jx+py+z]=P[Next+JX+(y0+y-Y0)*JY+(z0+z-Z0)];
						}
						else {

							BX2=boundX2;
							P[Current+(NlayersX)*Jx+py+z]=P[Current+BX2*Jx+py+z];
						}
					}
				}
			}
			else{

				BX2=boundX2;
				py=-Jy;
				for (y=0; y<=NlayersY; y++) {
					py +=Jy; for (z=0; z<=NlayersZ; z++) {
						P[Current+(NlayersX)*Jx+py+z]=P[Current+BX2*Jx+py+z];
					}
				}
			}
		}
	}
}
void Lat3D::SetBxm(double *P,double *Q) const {
	int i,y,z,py=0;
	//Message(debug,"SetBXm PQ");
	if (boundX2<numLayersX-1) {
		SetBxm(P); SetBxm(Q);
	}
	else if (N_comp_ranges==1) {
		py=-jy;
		for (y=0; y<MY; y++) {
			py +=jy; for (z=0; z<MZ; z++) {
				P[mmx*jx+py+z]=Q[boundX2*jx+py+z];
				Q[mmx*jx+py+z]=P[boundX2*jx+py+z];
			}
		}
	}
	else {
		//int Previous=0;
		int Current=0,Next=0;
		int JX=0,JY;
		//int X0;
		int Y0,Z0;
		//int x0;
		int y0,z0;
		//int xM,yM,zM;
		//int XM;
		int YM,ZM,BX2;
		int Jx=0,Jy=0,NlayersX,NlayersY,NlayersZ;
		for (i=1; i<=N_comp_ranges; i++) {
			//Previous=Current;
			Current=Next; Next +=CR_Q[i]->Get_M();
			NlayersX=CR_Q[i]->GetNumLayers_x()+1;
			NlayersY=CR_Q[i]->GetNumLayers_y()+1;
			NlayersZ=CR_Q[i]->GetNumLayers_z()+1;
			if (i<N_comp_ranges) {
				Jx=CR_Q[i]->Get_jx();
				Jy=CR_Q[i]->Get_jy();
				//x0=CR_Q[i]->Get_x0();
				y0=CR_Q[i]->Get_y0();
				z0=CR_Q[i]->Get_z0();
				//xM=CR_Q[i]->Get_xm();
				//yM=CR_Q[i]->Get_ym();
				//zM=CR_Q[i]->Get_zm();
				//X0=CR_Q[i+1]->Get_x0();
				Y0=CR_Q[i+1]->Get_y0();
				Z0=CR_Q[i+1]->Get_z0();
				//XM=CR_Q[i+1]->Get_xm();
				YM=CR_Q[i+1]->Get_ym();
				ZM=CR_Q[i+1]->Get_zm();
				JX=CR_Q[i+1]->Get_jx();
				JY=CR_Q[i+1]->Get_jy();
				py=-Jy;
				for (y=0; y<=NlayersY; y++) {
					py+=Jy; for (z=0; z<=NlayersZ; z++){
						if (y0+y-1 >= Y0 && y0+y-1<=YM && z0+z-1 >= Z0 && z0+z-1<=ZM )  {
							P[Current+(NlayersX)*Jx+py+z]=P[Next+JX+(y0+y-Y0)*JY+(z0+z-Z0)];
							Q[Current+(NlayersX)*Jx+py+z]=Q[Next+JX+(y0+y-Y0)*JY+(z0+z-Z0)];
						}
						else {

							BX2=boundX2;
							P[Current+(NlayersX)*Jx+py+z]=Q[Current+BX2*Jx+py+z];
							Q[Current+(NlayersX)*Jx+py+z]=P[Current+BX2*Jx+py+z];
						}
					}
				}
			}
			else{

				BX2=boundX2;
				py=-Jy;
				for (y=0; y<=NlayersY; y++) {
					py +=Jy; for (z=0; z<=NlayersZ; z++) {
						P[Current+(NlayersX)*Jx+py+z]=Q[Current+BX2*Jx+py+z];
						Q[Current+(NlayersX)*Jx+py+z]=P[Current+BX2*Jx+py+z];
					}
				}
			}
		}
	}
}

void Lat3D::SetBoundaries(Vector P) const {double *pP=&P[1]; SetBoundaries(pP);}
void Lat3D::SetBoundaries(double *P) const {
	int i,x,y,z,px=0,py=0;
	//Message(debug,"setboundaries");
	SetBz1(P); SetBzm(P); SetBy1(P); SetBym(P); SetBx1(P); SetBxm(P);
	//int Previous=0;
	int Current=0,Next=0,NlayersX,NlayersY,NlayersZ,BX1,BX2,BY1,BY2,BZ1,BZ2,mmx,mmy,mmz,Jx,Jy;
	for (i=1; i<=N_comp_ranges; i++) {
		//Previous=Current;
		Current=Next; Next +=CR_Q[i]->Get_M();
		NlayersX=CR_Q[i]->GetNumLayers_x(); mmx=NlayersX+1;
		NlayersY=CR_Q[i]->GetNumLayers_y(); mmy=NlayersY+1;
		NlayersZ=CR_Q[i]->GetNumLayers_z(); mmz=NlayersZ+1;
		Jx=CR_Q[i]->Get_jx();
		Jy=CR_Q[i]->Get_jy();

		BX1=boundX1;
		BY1=boundY1;
		BZ1=boundZ1;
		BX2=boundX2;
		BY2=boundY2;
		BZ2=boundZ2;

		for (z=1; z<mmz; z++) {
			P[Current+z]              =P[Current+BX1*Jx+BY1*Jy+z];
			P[Current+mmx*Jx+z]       =P[Current+BX2*Jx+BY1*Jy+z];
			P[Current+mmy*Jy+z]       =P[Current+BX1*Jx+BY2*Jy+z];
			P[Current+mmx*Jx+mmy*Jy+z]=P[Current+BX2*Jx+BY2*Jy+z];
		}
		py=0;
		for (y=1; y<mmy; y++) {py+=Jy;
			P[Current+py]            =P[Current+BX1*Jx+py+BZ1];
			P[Current+py+mmz]        =P[Current+BX1*Jx+py+BZ2];
			P[Current+mmx*Jx+py]     =P[Current+BX2*Jx+py+BZ1];
			P[Current+mmx*Jx+py+mmz] =P[Current+BX2*Jx+py+BZ2];
		}
		px=0;
		for (x=1; x<mmx; x++) {px+=Jx;
			P[Current+px]            =P[Current+px+BY1*Jy+BZ1];
			P[Current+px+mmy*Jy]     =P[Current+px+BY2*Jy+BZ1];
			P[Current+px+mmz]        =P[Current+px+BY1*Jy+BZ2];
			P[Current+px+mmy*Jy+mmz] =P[Current+px+BY2*Jy+BZ2];
		}

		P[Current+0] =                P[Current+BX1*Jx+BY1*Jy+BZ1];
		P[Current+mmz] =              P[Current+BX1*Jx+BY1*Jy+BZ2];
		P[Current+mmy*Jy] =           P[Current+BX1*Jx+BY2*Jy+BZ1];
		P[Current+mmx*Jx] =		      P[Current+BX2*Jx+BY1*Jy+BZ1];
		P[Current+mmy*Jy+mmz] =       P[Current+BX1*Jx+BY2*Jy+BZ2];
		P[Current+mmx*Jx+mmz] =       P[Current+BX2*Jx+BY1*Jy+BZ2];
		P[Current+mmx*Jx+mmy*Jy] =    P[Current+BX2*Jx+BY2*Jy+BZ1];
		P[Current+mmx*Jx+mmy*Jy+mmz]= P[Current+BX2*Jx+BY2*Jy+BZ2];
	}
}

void
Lat3D::SetBulkBoundaries(Vector in, const Vector bulkValues) const {
    //Message(fatal,"Bulkboundaries not implemented in 3d");
}

void Lat3D::removeboundaries(Vector P) const {double *pP=&P[1]; removeboundaries(pP); }
void Lat3D::removeboundaries(double *P) const {
	int i,x,y,z;
	//Message(debug,"removeboundaries");
	//int Previous=0;
	int Current=0,Next=0,NlayersX,NlayersY,NlayersZ,mmx,mmy,mmz,jx,jy;
		for (i=1; i<=N_comp_ranges; i++) {
			//Previous=Current;
			Current=Next; Next +=CR_Q[i]->Get_M();
			NlayersX=CR_Q[i]->GetNumLayers_x(); mmx=NlayersX+1;
			NlayersY=CR_Q[i]->GetNumLayers_y(); mmy=NlayersY+1;
			NlayersZ=CR_Q[i]->GetNumLayers_z(); mmz=NlayersZ+1;
			jx=CR_Q[i]->Get_jx();
			jy=CR_Q[i]->Get_jy();

			for (x=0; x<=mmx; x++) for (y=0; y<=mmy; y++) P[Current+x*jx+y*jy  ]=P[Current+x*jx+y*jy+mmz]=0;
			for (x=0; x<=mmx; x++) for (z=0; z<=mmz; z++) P[Current+x*jx+     z]=P[Current+x*jx+mmy*jy+z]=0;
			for (y=0; y<=mmy; y++) for (z=0; z<=mmz; z++) P[Current+     y*jy+z]=P[Current+mmx*jx+y*jy+z]=0;
		}
}

void Lat3D::SubtractBoundaries(Vector P) const {double *pP=&P[1]; SubtractBoundaries(pP);}
void Lat3D::SubtractBoundaries(double *P) const {
	int i,x,y,z;
	//int Previous=0;
	int Current=0,Next=0,NlayersX,NlayersY,NlayersZ,mmx,mmy,mmz,jx,jy,BX1,BX2,BY1,BY2,BZ1,BZ2;
	//Message(debug,"SubtractBoundaries");
	for (i=1; i<=N_comp_ranges; i++) {
		//Previous=Current;
		Current=Next; Next +=CR_Q[i]->Get_M();
		NlayersX=CR_Q[i]->GetNumLayers_x(); mmx=NlayersX+1;
		NlayersY=CR_Q[i]->GetNumLayers_y(); mmy=NlayersY+1;
		NlayersZ=CR_Q[i]->GetNumLayers_z(); mmz=NlayersZ+1;
		jx=CR_Q[i]->Get_jx();
		jy=CR_Q[i]->Get_jy();

		BX1=boundX1;
		BY1=boundY1;
		BZ1=boundZ1;
		BX2=boundX2;
		BY2=boundY2;
		BZ2=boundZ2;


		if (BZ1 ==2 || BZ2 == NlayersZ-1){
			for (x=1; x<mmx; x++) for (y=1; y<mmy; y++) {
				if (BZ1 ==2)        {P[Current+x*jx+y*jy+1]/=2;}
				if (BZ2 ==NlayersZ-1) {P[Current+x*jx+y*jy+NlayersZ]/=2; }
			}
		}
		if (BY1 ==2 || BY2 == NlayersY-1){
			for (x=1; x<mmx; x++) for (z=1; z<mmz; z++) {
				if (BY1 ==2)        {P[Current+x*jx+jy+z]/=2;}
				if (BY2 ==NlayersY-1) {P[Current+x*jx+NlayersY*jy+z]/=2; }
			}
		}
		if (BX1 ==2 || BX2 == NlayersX-1){
			for (y=1; y<mmy; y++) for (z=1; z<mmz; z++) {
				if (BX1 ==2)        {P[Current+jx+y*jy+z]/=2;}
				if (BX2 ==NlayersX-1) {P[Current+NlayersX*jx+y*jy+z]/=2; }
			}
		}
	}
	removeboundaries(P);
}


//void
//Lat2DFlat::SubtractBoundaries(Vector in) const {
//	int z=0;
//	for (int x=1; x<=MX; x++) {
//		for (int y=1; y<=MY; y++) {
//			z++;
//			if (x == 1 && boundX1 != 1) in[z] = 0;
//			else if (x == 2 && boundX1 == 3) in[z] /= 2;
//			else if (x == MX-1 && boundX2 == MX-2) in[z] /= 2;
//			else if (x == MX && boundX2 != MX) in[z] = 0;
//			if (y == 1 && boundY1 != 1) in[z] = 0;
//			else if (y == 2 && boundY1 == 3) in[z] /= 2;
//			else if (y == MY-1 && boundY2 == MY-2) in[z] /= 2;
//			else if (y == MY && boundY2 != MY) in[z] = 0;
//		}
//	}
//}


void Lat3D::RestoreBoundaries(Vector P) const {double *pP=&P[1];RestoreBoundaries(pP); }
void Lat3D::RestoreBoundaries(double *P) const {
	int i,x,y,z;
	//int Previous=0;
	int Current=0,Next=0,NlayersX,NlayersY,NlayersZ,mmx,mmy,mmz,jx,jy,BX1,BX2,BY1,BY2,BZ1,BZ2;
	//Message(debug,"RestoreBoundaries");
	for (i=1; i<=N_comp_ranges; i++) {
		//Previous=Current;
		Current=Next; Next +=CR_Q[i]->Get_M();
		NlayersX=CR_Q[i]->GetNumLayers_x(); mmx=NlayersX+1;
		NlayersY=CR_Q[i]->GetNumLayers_y(); mmy=NlayersY+1;
		NlayersZ=CR_Q[i]->GetNumLayers_z(); mmz=NlayersZ+1;
		jx=CR_Q[i]->Get_jx();
		jy=CR_Q[i]->Get_jy();

		BX1=boundX1;
		BY1=boundY1;
		BZ1=boundZ1;
		BX2=boundX2;
		BY2=boundY2;
		BZ2=boundZ2;

		if (BZ1 ==2 || BZ2 == NlayersZ-1){
			for (x=1; x<mmx; x++) for (y=1; y<mmy; y++) {
				if (BZ1 ==2)        {P[Current+x*jx+y*jy+1]*=2;}
				if (BZ2 ==NlayersZ-1) {P[Current+x*jx+y*jy+NlayersZ]*=2; }
			}
		}
		if (BY1 ==2 || BY2 == NlayersY-1){
			for (x=1; x<mmx; x++) for (z=1; z<mmz; z++) {
				if (BY1 ==2)        {P[Current+x*jx+jy+z]*=2;}
				if (BY2 ==NlayersY-1) {P[Current+x*jx+NlayersY*jy+z]*=2; }
			}
		}
		if (BX1 ==2 || BX2 == NlayersX-1){
			for (y=1; y<mmy; y++) for (z=1; z<mmz; z++) {
				if (BX1 ==2)        {P[Current+jx+y*jy+z]*=2;}
				if (BX2 ==NlayersX-1) {P[Current+NlayersX*jx+y*jy+z]*=2; }
			}
		}
	}
	SetBoundaries(P);
}

Array<Text>
Lat3D::GetNamesBulkBoundaries() const {
	return BulkBoundaries;
}

Boolean
Lat3D::BulkBoundary(const int z) const {
	//Message(fatal,"Bulkboundaries not implemented in 3d");
	return false;
}

int
Lat3D::NumBulkBoundary(const int z) const {
	//"Bulkboundaries not implemented in 3d"
	return 0;
}

/*
Checks whether the surface exists on every place it is defined:
Do the input latticeranges correspond to the set boundaries?
*/
void
Lat3D::CheckBoundaries(Array<LatticeRange*> LatRanges) const {
	// not necessary.
}

void Lat3D::Sides(Vector P, Vector Q) const {
	double *pP=&P[1], *pQ=&Q[1]; Sides(pP,pQ);
}

void Lat3D::Sides(double *P, double *Q) const { //this procedure is not tested yet!!
	int PrevM=0,M,jx,jy;
	double lambda=1.0/12.0;
	SetBoundaries(Q);
	if (latType == HEX)
	for (int i=1; i<=N_comp_ranges; i++){
		if (i>1) PrevM+=CR_Q[i-1]->Get_M();
		M=CR_Q[i]->Get_M();
		jx=CR_Q[i]->Get_jx();
		jy=CR_Q[i]->Get_jy();
		P+=PrevM;
		Q+=PrevM;
		set   (P,Q+jx,lambda,M-jx);
		addset(P+jx,Q,lambda,M-jx);
		addset(P+jy,Q+jx,lambda,M-jx-jy);
	       addset(P+jy,Q,lambda,M-jy);
		addset(P,Q+jy,lambda,M-jy);
		addset(P+jx,Q+jy,lambda,M-jx-jy);
		addset(P+jy,Q+1,lambda,M-jy-1);
		addset(P,Q+1,lambda,M-1);
		addset(P+jx,Q+1,lambda,M-jx-1);
		addset(P+1,Q+jx,lambda,M-jx-1);
		addset(P+1,Q,lambda,M-1);
		addset(P+1,Q+jy,lambda,M-jy-1);
		P-=PrevM;
		Q-=PrevM;
	}  else {

		for (int i=1; i<=N_comp_ranges; i++){
			if (i>1) PrevM+=CR_Q[i-1]->Get_M();
			M=CR_Q[i]->Get_M();
			jx=CR_Q[i]->Get_jx();
			jy=CR_Q[i]->Get_jy();
			P+=PrevM;
			Q+=PrevM;
			set(P,Q,l0,M);
			addset(P,Q+jx,l1,M-jx); addset(P+jx,Q,l1,M-jx);
	      	 	addset(P,Q+jy,l1,M-jy); addset(P+jy,Q,l1,M-jy);
			addset(P,Q+1,l1,M-1); addset(P+1,Q,l1,M-1);
			addset(P,Q+jx+jy,l2,M-jx-jy); addset(P,Q+jx-jy,l2,M-jx+jy);
			addset(P+jx-jy,Q,l2,M-jx+jy); addset(P+jx+jy,Q,l2,M-jx-jy);
			addset(P,Q+jx+1,l2,M-jx-1); addset(P,Q+jx-1,l2,M-jx+1);
			addset(P+jx-1,Q,l2,M-jx+1); addset(P+jx+1,Q,l2,M-jx-1);
			addset(P,Q+jy+1,l2,M-jy-1); addset(P,Q+jy-1,l2,M-jy+1);
			addset(P+jy-1,Q,l2,M-jy+1); addset(P+jy+1,Q,l2,M-jy-1);
			addset(P,Q+jx+jy+1,l3,M-jx-jy-1);addset(P,Q+jx+jy-1,l3,M-jx-jy+1);
			addset(P,Q+jx-jy+1,l3,M-jx+jy-1);addset(P,Q+jx-jy-1,l3,M-jx+jy+1);
			addset(P+jx-jy-1,Q,l3,M-jx+jy+1);addset(P+jx-jy+1,Q,l3,M-jx+jy-1);
			addset(P+jx+jy-1,Q,l3,M-jx-jy+1);addset(P+jx+jy+1,Q,l3,M-jx-jy-1);

			P-=PrevM;
			Q-=PrevM;
		}
		removeboundaries(P);
	}
}


double
Lat3D::SideFraction(Vector P, int position) const {
	// The input vector contains the data for which the sidefraction needs to be calculated.
	double result=0;
	double *pP=&P[1];
	int x,y,z,pos;

	if (latType == HEX) {Message(fatal, "side fraction in HEX not implemented. Ask Johan"); }

	pos=position-1;
	if (N_comp_ranges==1) {
		z = ((pos % jx) % jy);
		y = ((pos % jx)-z)/jy;
		x = (pos-y*jy-z)/jx;

		if (x > 0 && x <= numLayersX && y > 0 && y <= numLayersY && z > 0 && z <= numLayersZ) {
			if (x==1)			{pP[y*jy+z*jz]=    pP[boundX1*jx+y*jy+z];}
			if (x==numLayersX)	{pP[mmx*jx+y*jy+z]=pP[boundX2*jx+y*jy+z];}
			if (y==1)			{pP[x*jx+z*jz]=    pP[x*jx+boundY1*jy+z];}
			if (y==numLayersY)	{pP[x*jx+mmy*jy+z]=pP[x*jx+boundY2*jy+z];}
			if (z==1)			{pP[x*jx+y*jy]=    pP[x*jx+y*jy+boundZ1];}
			if (z==numLayersZ)	{pP[x*jx+y*jy+mmz]=pP[x*jx+y*jy+boundZ2];}
			if (l2 > 0) {
				if (x==1 && y==1) {                  pP[z]              =pP[boundX1*jx+boundY1*jy+z];}
				if (x==1 && z==1) {                  pP[y*jy]           =pP[boundX1*jx+y*jy+boundZ1];}
				if (y==1 && z==1) {                  pP[x*jx]           =pP[x*jx+boundY1*jy+boundZ1];}
				if (x==1 && y==numLayersY) {         pP[mmy*jy+z]       =pP[boundX1*jx+boundY2*jy+z];}
				if (x==1 && z==numLayersZ) {         pP[y*jy+mmz]       =pP[boundX1*jx+y*jy+boundZ2];}
				if (y==1 && z==numLayersZ) {         pP[x*jx+mmz]       =pP[x*jx+boundY1*jy+boundZ2];}
				if (y==1 && x==numLayersX) {         pP[mmx*jx+z]       =pP[boundX2*jx+boundY1*jy+z];}
				if (z==1 && y==numLayersY) {         pP[x*jx+mmy*jy]    =pP[x*jx+boundY2*jy+boundZ1];}
				if (z==1 && x==numLayersX) {         pP[mmx*jx+y*jy]    =pP[boundX2*jx+y*jy+boundZ1];}
				if (x==numLayersX && y==numLayersY) {pP[mmx*jx+mmy*jy+z]=pP[boundX2*jx+boundY2*jy+z];}
				if (x==numLayersX && z==numLayersZ) {pP[mmx*jx+y*jy+mmx]=pP[boundX2*jx+y*jy+boundZ2];}
				if (y==numLayersY && z==numLayersZ) {pP[x*jx+mmy*jy+mmz]=pP[x*jx+boundY2*jy+boundZ2];}
			}
			if (l3 > 0) {
				if (x==1 && y==1 && z==1) {                           pP[0] =               pP[boundX1*jx+boundY1*jy+boundZ1];}
				if (x==1 && y==1 && z==numLayersZ) {                  pP[mmz] =             pP[boundX1*jx+boundY1*jy+boundZ2];}
				if (x==1 && y==numLayersY && z==1) {                  pP[mmy*jy] =          pP[boundX1*jx+boundY2*jy+boundZ1];}
				if (x==numLayersX && y==1 && z==1) {                  pP[mmx*jx] =          pP[boundX2*jx+boundY1*jy+boundZ1];}
				if (x==1 && y==numLayersY && z==numLayersZ) {         pP[mmy*jy+mmz] =      pP[boundX1*jx+boundY2*jy+boundZ2];}
				if (x==numLayersX && y==1 && z==numLayersZ) {         pP[mmx*jx+mmz] =      pP[boundX2*jx+boundY1*jy+boundZ2];}
				if (x==numLayersX && y==numLayersY && z==1) {         pP[mmx*jx+mmy*jy] =   pP[boundX2*jx+boundY2*jy+boundZ1];}
				if (x==numLayersX && y==numLayersY && z==numLayersZ) {pP[mmx*jx+mmy*jy+mmz]=pP[boundX2*jx+boundY2*jy+boundZ2];}
			}
			result= l1*(pP[pos-jx]+pP[pos+jx]+pP[pos-jy]+pP[pos+jy]+pP[pos-1]+pP[pos+1]);
			if (l0 > 0) {result += l0*pP[pos]; }
			if (l2 > 0){
				result+= l2*(pP[pos-jx-jy]+pP[pos-jx+jy]+pP[pos-jy-1]+pP[pos-jy+1]+pP[pos-jx-1]+pP[pos-jx+1]+
						 pP[pos+jx-jy]+pP[pos+jx+jy]+pP[pos+jy-1]+pP[pos+jy+1]+pP[pos+jx-1]+pP[pos+jx+1]);
			}
			if (l3 > 0){
				result+= l3*(pP[pos-jx-jy-1]+pP[pos-jx-jy+1]+pP[pos-jx+jy-1]+pP[pos-jx+jy+1])+
				         pP[pos+jx-jy-1]+pP[pos+jx-jy+1]+pP[pos+jx-jy-1]+pP[pos+jx-jy+1];
			}
		}
	}
	else {
		int i,jx,jy,NlayersX,NlayersY,NlayersZ,mmx,mmy,mmz,BX1,BX2,BY1,BY2,BZ1,BZ2;
		//int Previous=0;
		int Current=0,Next=0;
		for (i=1; i<=N_comp_ranges; i++){
			//Previous=Current;
			Current=Next; Next +=CR_Q[i]->Get_M();
			if (pos>=0 && pos < CR_Q[i]->Get_M()) {

				jx=CR_Q[i]->Get_jx();
				jy=CR_Q[i]->Get_jy();

				z = ((pos % jx) % jy);
				y = ((pos % jx)-z)/jy;
				x = (pos-y*jy-z)/jx;


				NlayersX=CR_Q[i]->GetNumLayers_x(); mmx=NlayersX+1;
				NlayersY=CR_Q[i]->GetNumLayers_y(); mmy=NlayersY+1;
				NlayersZ=CR_Q[i]->GetNumLayers_z(); mmz=NlayersZ+1;

				BX1=boundX1;
				BY1=boundY1;
				BZ1=boundZ1;
				BX2=boundX2;
				BY2=boundY2;
				BZ2=boundZ2;
				if (x > 0 && x <= NlayersX && y > 0 && y <= NlayersY && z > 0 && z <= NlayersZ) {
					if (x==1)			SetBx1(pP);
					if (x==NlayersX)	SetBxm(pP);
					if (y==1)			{pP[Current+x*jx+z]=       pP[Current+x*jx+BY1*jy+z];}
					if (y==NlayersY)	{pP[Current+x*jx+mmy*jy+z]=pP[Current+x*jx+BY2*jy+z];}
					if (z==1)			{pP[Current+x*jx+y*jy]=    pP[Current+x*jx+y*jy+BZ1];}
					if (z==NlayersZ)	{pP[Current+x*jx+y*jy+mmz]=pP[Current+x*jx+y*jy+BZ2];}
					if (l2 > 0) {
						if (x==1 && y==1) {              pP[Current+z]              =pP[Current+BX1*jx+BY1*jy+z];}
						if (x==1 && z==1) {              pP[Current+y*jy]           =pP[Current+BX1*jx+y*jy+BZ1];}
						if (y==1 && z==1) {              pP[Current+x*jx]           =pP[Current+x*jx+BY1*jy+BZ1];}
						if (x==1 && y==NlayersY) {       pP[Current+mmy*jy+z]       =pP[Current+BX1*jx+BY2*jy+z];}
						if (x==1 && z==NlayersZ) {       pP[Current+y*jy+mmz]       =pP[Current+BX1*jx+y*jy+BZ2];}
						if (y==1 && z==NlayersZ) {       pP[Current+x*jx+mmz]       =pP[Current+x*jx+BY1*jy+BZ2];}
						if (y==1 && x==NlayersX) {       pP[Current+mmx*jx+z]       =pP[Current+BX2*jx+BY1*jy+z];}
						if (z==1 && y==NlayersY) {       pP[Current+x*jx+mmy*jy]    =pP[Current+x*jx+BY2*jy+BZ1];}
						if (z==1 && x==NlayersX) {       pP[Current+mmx*jx+y*jy]    =pP[Current+BX2*jx+y*jy+BZ1];}
						if (x==NlayersX && y==NlayersY) {pP[Current+mmx*jx+mmy*jy+z]=pP[Current+BX2*jx+BY2*jy+z];}
						if (x==NlayersX && z==NlayersZ) {pP[Current+mmx*jx+y*jy+mmx]=pP[Current+BX2*jx+y*jy+BZ2];}
						if (y==NlayersY && z==NlayersZ) {pP[Current+x*jx+mmy*jy+mmz]=pP[Current+x*jx+BY2*jy+BZ2];}
					}
					if (l3 > 0) {
						if (x==1 && y==1 && z==1) {                     pP[Current+0] =               pP[Current+BX1*jx+BY1*jy+BZ1];}
						if (x==1 && y==1 && z==NlayersZ) {              pP[Current+mmz] =             pP[Current+BX1*jx+BY1*jy+BZ2];}
						if (x==1 && y==NlayersY && z==1) {              pP[Current+mmy*jy] =          pP[Current+BX1*jx+BY2*jy+BZ1];}
						if (x==NlayersX && y==1 && z==1) {              pP[Current+mmx*jx] =          pP[Current+BX2*jx+BY1*jy+BZ1];}
						if (x==1 && y==NlayersY && z==NlayersZ) {       pP[Current+mmy*jy+mmz] =      pP[Current+BX1*jx+BY2*jy+BZ2];}
						if (x==NlayersX && y==1 && z==NlayersZ) {       pP[Current+mmx*jx+mmz] =      pP[Current+BX2*jx+BY1*jy+BZ2];}
						if (x==NlayersX && y==NlayersY && z==1) {       pP[Current+mmx*jx+mmy*jy] =   pP[Current+BX2*jx+BY2*jy+BZ1];}
						if (x==NlayersX && y==NlayersY && z==NlayersZ) {pP[Current+mmx*jx+mmy*jy+mmz]=pP[Current+BX2*jx+BY2*jy+BZ2];}
					}
					result= l1*(pP[Current+pos-jx]+pP[Current+pos+jx]+pP[Current+pos-jy]+
								    pP[Current+pos+jy]+pP[Current+pos-1]+pP[Current+pos+1]);
					if (l0 > 0) {result += l0*pP[Current+pos];}
					if (l2 > 0){
						result+= l2*(pP[Current+pos-jx-jy]+pP[Current+pos-jx+jy]+pP[Current+pos-jy-1]+
									     pP[Current+pos-jy+1]+pP[Current+pos-jx-1]+pP[Current+pos-jx+1]+
									     pP[Current+pos+jx-jy]+pP[Current+pos+jx+jy]+pP[Current+pos+jy-1]+
									     pP[Current+pos+jy+1]+pP[Current+pos+jx-1]+pP[Current+pos+jx+1]);
					}
					if (l3 > 0){
						result+= l3*(pP[Current+pos-jx-jy-1]+pP[Current+pos-jx-jy+1]+pP[Current+pos-jx+jy-1]+pP[Current+pos-jx+jy+1])+
							     pP[Current+pos+jx-jy-1]+pP[Current+pos+jx-jy+1]+pP[Current+pos+jx-jy-1]+pP[Current+pos+jx-jy+1];
					}
					pos = -1; //to prevent that this procedure is entered again.

				}
			}
			else {
				pos-=CR_Q[i]->Get_M();
			}
		}
	}
    return result;
}

void
Lat3D::ElectricPotential(Vector &psi,Vector &psi0, const Vector &eps, const Vector &charge, const double preFactor) const {
	double *p_psi = &psi[1];
	double *p_psi0 = &psi0[1];
	const double *p_eps = &eps[1];
	const double *p_charge = &charge[1];
	ElectricPotential(p_psi,p_psi0,p_eps,p_charge,preFactor);
}
void
Lat3D::ElectricPotential(double *psi, double *psi0, const double *eps0, const double *charge, const double preFactor) const {
	SetBoundaries(psi0);
	double *eps = const_cast<double*>(eps0);
	SetBoundaries(eps);
	Vector eps_tmp,eps_tot,Eps_;
	eps_tmp.Dim(1,Mmax); eps_tot.Dim(1,Mmax);
	double *epstmp = &eps_tmp[1]; double *epstot = &eps_tot[1];

	int i,M,jx,jy,prevM=0;

	for (i=1; i<=N_comp_ranges; i++){
		if (i>1) prevM+=CR_Q[i-1]->Get_M();
		M=CR_Q[i]->Get_M();
		jx=CR_Q[i]->Get_jx();
		jy=CR_Q[i]->Get_jy();
		eps+=prevM;
		psi+=prevM;
		psi0+=prevM;
		charge+=prevM;

		add(epstmp+jx,eps,   eps+jx,M-jx);   cp(epstot,epstmp,M);    times(psi+jx,epstmp+jx,psi0,   M-jx);
		add(epstmp,   eps+jx,eps,   M-jx); add2(epstot,epstmp,M); addTimes(psi,   epstmp,   psi0+jx,M-jx);
		add(epstmp+jy,eps,   eps+jy,M-jy); add2(epstot,epstmp,M); addTimes(psi+jy,epstmp+jy,psi0,   M-jy);
		add(epstmp,   eps+jy,eps,   M-jy); add2(epstot,epstmp,M); addTimes(psi,   epstmp,   psi0+jy,M-jy);
		add(epstmp+1, eps,   eps+jz,M-1);  add2(epstot,epstmp,M); addTimes(psi+1, epstmp+1, psi0,   M-1);
		add(epstmp,   eps+1, eps,   M-1);  add2(epstot,epstmp,M); addTimes(psi,   epstmp,   psi0+1, M-1);
		addset(psi,charge,2*preFactor,M); div(psi,epstot,M);
		eps-=prevM;
		psi-=prevM;
		psi0-=prevM;
		charge-=prevM;
	}
	SetBoundaries(psi);
}
void
Lat3D::ElectricFieldSquared(Vector &E2,const Vector psi0,const double preFactor) const {
	double *p_E2=&E2[1];
	const double *p_psi0 = &psi0[1];
	ElectricFieldSquared(p_E2,p_psi0,preFactor);
}

void
Lat3D::ElectricFieldSquared(double *E2,const double *psi, const double preFactor) const {
	//Vector psi_0(1,GetTotalNumLayers());
	double *psi0 = const_cast<double*>(psi);
	//double *psi0=&psi_0[1]; cp(psi0,psi,GetTotalNumLayers());
    SetBoundaries(psi0);
    int i,M,jx,jy,prevM=0;
    for (i=1; i<=N_comp_ranges; i++){
    	if (i>1) prevM=CR_Q[i-1]->Get_M();
    	M=CR_Q[i]->Get_M();
    	jx=CR_Q[i]->Get_jx();
    	jy=CR_Q[i]->Get_jy();
    	E2+=prevM;
    	psi0+=prevM;

        min2(E2+jx,psi0+jx,psi0,M-jx);
        addmin2(E2,psi0,psi0+jx,M-jx);
        addmin2(E2+jy,psi0+jy,psi0,M-jy);
        addmin2(E2,psi0,psi0+jy,M-jy);
        addmin2(E2+1,psi0+1,psi0,M-1);
        addmin2(E2,psi0,psi0+1,M-1);
        norm(E2,preFactor/2,M);
    	E2-=prevM;
    	psi0-=prevM;
    }
    SetBoundaries(E2);
}
Vector
Lat3D::Div1Grad2(const Vector in1, const Vector in2) const {
	Vector value(1,M);
      Message(implementation,"Lat3D::Div1Grad2");
	return value;
}

Vector
Lat3D::FourierTransform(const Vector /*in*/) const {
	Vector value(1,M);
	Message(implementation,"Lat3D::FourierTransform");
	return value;
}
double
Lat3D::Moment(const Vector,
				  const int,
				  const LatticeRange*,
				  const double) const {
	double value = 0;
	// als meer dan 1 coordinate in Latrange : error
	Message(implementation,"Lat3D::Moment");
	return value;
}
double
Lat3D::MomentUnweighted(const Vector,
							const int,
							const LatticeRange*,
							const double) const {
	double value = 0;
	// als meer dan 1 coordinate in Latrange : error
	Message(implementation,"Lat3D::MomentUnweighted");
	return value;
}
Vector
Lat3D::RenormPhi(Vector a,
					 Vector,
					 const Vector,
					 const  double) const {
	Message(fatal,"Cannot have a molecule with freedom 'thirdGeneration' "
		"on a flat lattice with three gradients");
	return a; //never get here
}
double
Lat3D::GetLayerAdjustment(void) const {
	return 1;
}
void
Lat3D::SetLayerAdjustment(double) {
	Message(fatal,"don't know how to eliminate lattice artifact in 3 gradients");
}
