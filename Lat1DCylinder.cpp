#include "Lat1DCylinder.h"
Lat1DCylinder::Lat1DCylinder() {
}
Lat1DCylinder::Lat1DCylinder(Input* MyInput_, Text name_)
 	: Lat1DFlat(MyInput_,name_) {
	Array<Text> param(1,15);
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
	param[11] = "n_sites_first_layer";
	param[12] = "offset_first_layer";
	param[13] = "layer_adjustment";
	param[14] = "n_layers_extended";
	param[15] = "lambda_chi";
	MyInput->CheckParameterNames("lat",name,param);
	MyInput->DontCombineParam("lat",name,"n_sites_first_layer",
				  "offset_first_layer");
	MyInput->DontCombineParam("lat",name,"n_sites_first_layer",
				  "layer_adjustment");
//	MyInput->DontCombineParam("lat",name,"offset_first_layer",
//				  "layer_adjustment");
	if (adjustOuterLayer && (MyInput->ValueSet("lat", name,
				 "n_sites_first_layer"))) {
		Message(implementation, "bool adjustOuterLayer may not be true when "
			"n_sites_first_layer is set.");
	}
	lambda_1.Dim(1,M);
	lambda1.Dim(1,M);

	if (l1 != l1_chi) {
		lambda_1_chi.Dim(1,M);
		lambda1_chi.Dim(1,M);
	}
	L.Dim(1,M);
	if (*Copy(MyInput->GetText("lat",name,"geometry")) == *Copy("cylindrical")) {
		// layerAdjustment must be >= 2*l1 else l0[2] < 0
		// and layerAdjustment must be <= 2-2*l1 else l0[M-1] < 0
		layerAdjustment = MyInput->GetReal("lat",name,"layer_adjustment",0,2,1);
		if (MyInput->ValueSet("lat", name, "n_sites_first_layer")) {
			double v = MyInput->GetReal("lat",name,"n_sites_first_layer",PI,DBL_MAX,PI);
			layerAdjustment = sqrt(v/PI);
		}
		offsetLayer1 = MyInput->GetReal("lat",name,"offset_first_layer",0,DBL_MAX,0);
		if (offsetLayer1 <= 0 && bound1 == 1) {
			Message(fatal, "A calculation with a surface at the lowerbound in a curved "
				"lattice cannot make sense without setting offset_first_layer");
		}
		calculateLambdas();
	}
}
Lat1DCylinder::~Lat1DCylinder() {
}
void
Lat1DCylinder::GetOutput(Output* Out) const {
	Out->PutInt("lat",name,"gradients",1);
	if (dimensions - int(dimensions) > 1e-7) {
		Out->PutReal("lat",name,"dimensions",dimensions);
	} else {
		Out->PutInt("lat",name,"dimensions",int(dimensions));
	}
	Out->PutText("lat",name,"geometry","cylindrical");
	Out->PutInt("lat",name,"n_layers",numLayers);
	Out->PutReal("lat",name,"total number of lattice sites",
		     GetTotalNumLatticeSites());
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
	Out->PutReal("lat",name,"lambda_chi",l1_chi);
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
	Out->PutReal("lat",name,"n_sites_first_layer",L[2]);
	MyInput->SetDefaultWarningsOff();
	Out->PutText("lat",name,"offset_first_layer",MyInput->GetText("lat",name,"offset_first_layer","0"));
	MyInput->SetDefaultWarningsOff();
	Out->PutReal("lat",name,"layer_adjustment",layerAdjustment);
	Out->PutProfile("lat",name,"lattice sites",L);
	Out->PutProfile("lat",name,"lambda_1",lambda_1);
	Out->PutProfile("lat",name,"lambda1",lambda1);
	if (l1 != l1_chi) {
		Out->PutProfile("lat",name,"lambda_1_chi",lambda_1_chi);
		Out->PutProfile("lat",name,"lambda1_chi",lambda1_chi);
	}

}

double
Lat1DCylinder::GetNumLatticeSites(const int z) const {
	return L[z];
}
double
Lat1DCylinder::GetLambda(const int z1, const int z2) const {
	if (z1 == z2) return 1-lambda1[z1]-lambda_1[z1];
	if (z1 - z2 == 1) return lambda_1[z1];
	if (z1 - z2 == -1) return lambda1[z1];
	return 0;
}
void
Lat1DCylinder::MultiplyWithLatticeSites(Vector A) const {
	int z;
	for (z=1; z<=M; z++) {
		A[z] *= L[z];
	}
}
void
Lat1DCylinder::DivideByLatticeSites(Vector A) const {
	int z;
	for (z=1; z<=M; z++) {
		if (A[z] != 0) {
			A[z] /= L[z];
		}
	}
}
double
Lat1DCylinder::GetVolume() const {
	return GetTotalNumLatticeSites()*pow(siteDistance,dimensions-1);
}
double
Lat1DCylinder::GetTotalNumLatticeSites() const {
	double value = 0;
	int z;
	for (z=1; z<=M; z++) {
		value += L[z];
	}
	if (bound1 != 1) value -= L[1];
	if (bound1 == 3) value -= 0.5*L[2];
	if (bound2 != M) value -= L[M];
	if (bound2 == M-2) value -= 0.5*L[M-1];
	return value;
}
double
Lat1DCylinder::SideFraction(const Vector in, int z) const {
	if (z == 1) {
		if (bound1 == 1) return in[z] + lambda1_chi[z]*(- in[z] - in[z] + in[z+1]);
		else z = bound1;
	}
	if (z == M) {
		if (bound2 == M) return in[z] + lambda_1_chi[z]*(- in[z] - in[z] + in[z-1]);
		else z = bound2;
	}
	return in[z] + lambda1_chi[z]*(in[z+1] - in[z])
		+ lambda_1_chi[z]*(in[z-1] - in[z]);
}

void
Lat1DCylinder::Sides(Vector P, Vector Q) const {
	Message(fatal,"Sides not implemented. ");
}
void
Lat1DCylinder::ElectricPotential(Vector &psi, Vector &psi0,
		const Vector &eps, const Vector &charge,
		const double preFactor) const {
	//the charge is for the entire layer, not 1 lattice site
	int z;
	double Z, pf, epsmin, epsplus;

	SetBoundaries(psi0);
	pf = preFactor/PI;
	Z = offsetLayer1+layerAdjustment-1; // radius to outer edge of layer 1
	if (offsetLayer1 <= 0) {
		Z++;
		epsplus = Z*(eps[2]+eps[3]);
		psi[2] = psi0[3] + pf*charge[2]/epsplus;
		z = 3;
	} else {
		epsplus = Z*(eps[1]+eps[2]);
		if (bound1 == 1) {
			psi[1] = psi0[2] + pf*charge[1]/epsplus;
		}
		z = 2;
	}
	for (; z<M; z++) {
		Z++;
		epsmin = epsplus;
		epsplus = Z*(eps[z]+eps[z+1]);
		psi[z] = (pf*charge[z] + psi0[z-1]*epsmin + psi0[z+1]*epsplus)
			/ (epsmin + epsplus);
	}
	if (bound2 == M) {
		psi[M] = psi0[M-1] + pf*charge[M]/epsplus;
	}
}
void
Lat1DCylinder::ElectricFieldSquared(Vector &ESquared,
		const Vector psi,
		const double preFactor) const {
	int z;
	double Emin, Eplus, pf, Z;

	pf = preFactor/2;
	Z = offsetLayer1+layerAdjustment-1.5; // radius to center of layer 1
	if (offsetLayer1 <= 0) {
		Z++;
		Eplus = psi[2]-psi[3];
		Eplus *= (Z+0.5)*Eplus;
		ESquared[2] = pf*Eplus/Z;
		z = 3;
	} else {
		Eplus = psi[1]-psi[2];
		Eplus *= (Z+0.5)*Eplus;
		if (bound1 == 1) {
			ESquared[1] = pf*Eplus/Z;
		}
		z = 2;
	}
	for (; z<M; z++) {
		Z++;
		Emin = Eplus;
		Eplus = psi[z]-psi[z+1];
		Eplus *= (Z+0.5)*Eplus;
		ESquared[z] = pf*(Emin + Eplus)/Z;
	}
	if (bound2 == M) {
		ESquared[M] = pf*Eplus/Z;
	}
}
Vector
Lat1DCylinder::Div1Grad2(const Vector in1, const Vector in2) const {
	Vector value(1,M);
	for (int z=1; z<=M; z++) {
		double Z = z + offsetLayer1 - 2 + layerAdjustment;
		if (z > 1 && in1[z-1] != 0 && in1[z] != 0) {
			value[z] += ((Z-1)*in1[z-1]+Z*in1[z])*(in2[z-1]-in2[z]);
		}
		if (z < M && in1[z+1] != 0 && in1[z] != 0) {
			value[z] += ((Z+1)*in1[z+1]+Z*in1[z])*(in2[z+1]-in2[z]);
		}
		value[z] /= 2*siteDistance*siteDistance*Z;
	}
	return value;
}
Vector
Lat1DCylinder::FourierTransform(const Vector /*in*/) const {
	Message(implementation,"Lat1DCylinder::FourierTransform");
	Vector value(1,M);
	return value;
}
double
Lat1DCylinder::Moment(const Vector in,
					  const int m,
					  const LatticeRange* Range,
					  const double offset) const {
	double value = 0;
	Vector copy(1,M);
	int z;
	for (z=1; z<=M; z++) {
		copy[z] = in[z];
	}
	SubtractBoundaries(copy);
	if (m==0) {
		MultiplyWithLatticeSites(copy);
		for (z=1; z<=M; z++) {
			value += copy[z];
		}
		return value;
	}
	double Z;
	double Z2;
	int zRange=0;
	bool found = false;
	for (z=1; z<=M; z++) {
		if (Range->InRange(z)) {
			if (found) {
				Message(warning,"A moment is calculated wrong, because of a"
				"range that is more than one lattice layer");
			}
			zRange = z;
			found = true;
		}
	}
	double weight = 0;
	for (z=1; z<=M; z++) {
		Z = z-zRange+offset+0.5;
		Z2 = z-2+offsetLayer1+layerAdjustment;
		if (z == 1 && offsetLayer1 < 1) {
			weight = (pow(Z,m+1)*Z2-pow((Z-offsetLayer1),m+1)*(Z2-offsetLayer1))/(m+1.0);
			weight -= (pow(Z,m+2)-pow((Z-offsetLayer1),m+2))/((m+1)*(m+2.0));
		} else if (z == 2) {
			weight = (pow(Z,m+1)*Z2-pow((Z-layerAdjustment),m+1)*(Z2-layerAdjustment))/(m+1.0);
			weight -= (pow(Z,m+2)-pow((Z-layerAdjustment),m+2))/((m+1)*(m+2.0));
		} else {
			weight = (pow(Z,m+1)*Z2-pow((Z-1),m+1)*(Z2-1))/(m+1.0);
			weight -= (pow(Z,m+2)-pow((Z-1),m+2))/((m+1)*(m+2.0));
		}
		value += copy[z]*weight*2*PI;
	}
	return value;
}
double
Lat1DCylinder::MomentUnweighted(const Vector in,
								const int m,
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
	double weight = 0;
	for (z=1; z<=M; z++) {
		if (layerAdjustment != 1 && zRange == 1) {
			Z = z-2+offsetLayer1+layerAdjustment;
		} else {
			Z = z-zRange+offset+0.5;
		}
		if (z == 1 && offsetLayer1 < 1) {
			weight = (pow(Z,m+1)-pow(Z-offsetLayer1,m+1))/(m+1.0);
		} else if (z == 2) {
			weight = (pow(Z,m+1)-pow(Z-layerAdjustment,m+1))/(m+1.0);
		} else {
			weight = (pow(Z,m+1)-pow(Z-1,m+1))/(m+1.0);
		}
		value += copy[z]*weight;
	}
	return value;
}
void
Lat1DCylinder::SetLayerAdjustment(double value) {
	adjustOuterLayer = true;
	layerAdjustment = value;
	calculateLambdas();
}
void
Lat1DCylinder::calculateLambdas() {
	Vector V(1,M);
	Vector S(1,M);
	S[1] = offsetLayer1*PI*2;
	V[1] = offsetLayer1*offsetLayer1*PI;
	L[1] = V[1];
	int z;
	double x;
	for (z=2; z<=M; z++) {
		x = z + offsetLayer1 - 2 + layerAdjustment;
		if (z == M-1 && adjustOuterLayer) {
			x += -layerAdjustment + 1;
		}
		double D = dimensions - 1;
		V[z] = (PI* pow(2,(D-1))/D) * pow(x,D);
		L[z] = V[z] - V[z-1];
		S[z] = pow(2,(D-1))*PI*pow(x,(D-1));
		lambda_1[z] = l1*S[z-1]/L[z];
		lambda1[z] = l1*S[z]/L[z];
		if (l1 != l1_chi) {
			lambda_1_chi[z] = l1_chi*S[z-1]/L[z];
			lambda1_chi[z] = l1_chi*S[z]/L[z];
		}
	}
	if (L[1] > PI) L[1] = V[1] - (offsetLayer1-1)*(offsetLayer1-1)*PI;
	if (L[1] > 0) lambda1[1] = L[2]*lambda_1[2]/L[1];
	if (L[1] > 0 && l1 != l1_chi) lambda1_chi[1] = L[2]*lambda_1_chi[2]/L[1];
	if (MyInput->ValueSet("lat",name,"n_layers_extended")) {
		int numLayersExtended = MyInput->GetInt("lat",name,"n_layers_extended",numLayers,INT_MAX);
		double numExtraLayers = numLayersExtended - numLayers;
		if (bound2 == M - 2) { // mirror2
			numExtraLayers += 0.5;
		}
		x = M + offsetLayer1 - 3 + layerAdjustment;
		double x2 = x + numExtraLayers;
		double D = dimensions - 1;
		numExtraLatticeSites = (PI* pow(2,(D-1))/D) * pow(x2,D) - (PI* pow(2,(D-1))/D) * pow(x,D);;
	} else {
		numExtraLatticeSites = 0;
	}
	if (l1 == l1_chi) {
		lambda_1_chi = lambda_1;
		lambda1_chi = lambda1;
	}
}
