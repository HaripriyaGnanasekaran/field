#include "Lat1DSphere.h"

Lat1DSphere::Lat1DSphere(Input* MyInput_, Text name_)
	: Lat1DCylinder(MyInput_,name_) {

	if (!MyInput->ValueSet("lat",name,"n_sites_first_layer")) {
		layerAdjustment = MyInput->GetReal("lat",name,"layer_adjustment",0,2,1);
	} else {
		double vol = MyInput->GetReal("lat",name,"n_sites_first_layer",1,DBL_MAX,1);
		layerAdjustment = pow(3*vol/(4*PI),1.0/3.0);
	}
	offsetLayer1 = MyInput->GetReal("lat",name,"offset_first_layer",0,DBL_MAX,0);
	if (offsetLayer1 <= 0 && bound1 == 1) {
		Message(fatal, "A calculation with a surface at the lowerbound in a curved "
			"lattice cannot make sense without setting offset_first_layer");
	}
	calculateLambdas();
}
Lat1DSphere::~Lat1DSphere() {
}
void
Lat1DSphere::GetOutput(Output* Out) const {
	Out->PutInt("lat",name,"gradients",1);
	if (dimensions - int(dimensions) > 1e-7) {
		Out->PutReal("lat",name,"dimensions",dimensions);
	} else {
		Out->PutInt("lat",name,"dimensions",int(dimensions));
	}
	Out->PutText("lat",name,"geometry","spherical");
	Out->PutInt("lat",name,"n_layers",numLayers);
	Out->PutReal("lat",name,"total number of lattice sites",GetTotalNumLatticeSites());
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
	if (l1 != l1_chi) {
		Out->PutReal("lat",name,"lambda_chi",l1_chi);
	}
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
void
Lat1DSphere::ElectricPotential(Vector &psi, Vector &psi0,
		const Vector &eps, const Vector &charge,
		const double preFactor) const {
	//the charge is for the entire layer, not 1 lattice site
	int z;
	double Z, pf, epsmin, epsplus;

	pf = preFactor/2/PI;
	Z = offsetLayer1+layerAdjustment-1; //outer radius of layer 1
	if (offsetLayer1 <= 0) {
		Z++;
		epsplus = Z*Z*(eps[2]+eps[3]);
		psi[2] = psi0[3] + pf*charge[2]/epsplus;
		z = 3;
	} else {
		epsplus = Z*Z*(eps[1]+eps[2]);
		if (bound1 == 1) {
			psi[1] = psi0[2] + pf*charge[1]/epsplus;
		}
		z = 2;
	}
	for (; z<M; z++) {
		Z++;
		epsmin = epsplus;
		epsplus = Z*Z*(eps[z]+eps[z+1]);
		psi[z] = (pf*charge[z] + psi0[z-1]*epsmin + psi0[z+1]*epsplus)
			/ (epsmin + epsplus);
	}
	if (bound2 == M) {
		psi[M] = psi0[M-1] + pf*charge[M]/epsplus;
	}
}
void
Lat1DSphere::ElectricFieldSquared(Vector &ESquared,
		const Vector psi,
		const double preFactor) const {
	int z;
	double Z,pf,Eplus,Emin;

	pf = preFactor*PI*2;
	Z = offsetLayer1+layerAdjustment-1; //outer radius of layer z
	Eplus = Z*(psi[1]-psi[2]);
	Eplus *= Eplus;
	if (bound1 == 1) {
		ESquared[1] = pf*Eplus/L[1];
	}
	for (z=2; z<M; z++) {
		Z++;
		Emin = Eplus;
		Eplus = Z*(psi[z]-psi[z+1]);
		Eplus *= Eplus;
		ESquared[z] = pf*(Emin + Eplus)/L[z];
	}
	if (bound2 == M) {
		ESquared[M] = pf*Eplus/L[M];
	}
}
Vector
Lat1DSphere::Div1Grad2(const Vector in1, const Vector in2) const {
	Vector value(1,M);
	double Z;
	for (int z=1; z<=M; z++) {
		Z = z -1 + offsetLayer1;
		if (z > 1 && in1[z-1] != 0 && in1[z] != 0) {
			value[z] += ((Z-1)*(Z-1)*in1[z-1]+Z*Z*in1[z])*(in2[z-1]-in2[z]);
		}
		if (z < M && in1[z+1] != 0 && in1[z] != 0) {
			value[z] += ((Z+1)*(Z+1)*in1[z+1]+Z*Z*in1[z])*(in2[z+1]-in2[z]);
		}
		value[z] /= 2*siteDistance*siteDistance*Z*Z;
	}
	return value;
}
Vector
Lat1DSphere::FourierTransform(const Vector /*in*/) const {
	Message(implementation,"Lat1DSphere::FourierTransform");
	Vector value(1,M);
	return value;
}
double
Lat1DSphere::Moment(const Vector in,
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
	if (dimensions != 3) {
		Message(warning, "assuming moment is taken from center "
				"to outside for dimensions != 3");
		double Z;
		double D = dimensions;
		double weight = 0;
		for (z=1; z<=M; z++) {
			Z = z-2+offsetLayer1+layerAdjustment;
			if (z == 1 && offsetLayer1 < 1) {
				weight = (pow(Z,D+m)-pow((Z-offsetLayer1),D+m))/(m+D);
			} else if (z == 2) {
				weight = (pow(Z,D+m)-pow((Z-layerAdjustment),D+m))/(m+D);
			} else {
				weight = (pow(Z,D+m)-pow((Z-1),D+m))/(m+D);
			}
			value += copy[z]*weight*pow(2,D-1)*PI;
		}
		return value;
	}
	double Z;
	double Z2;
	int zRange = 0;
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
		Z = z-zRange+offset+0.5;
		Z2 = z-2+offsetLayer1+layerAdjustment;
		if (layerAdjustment != 1 && zRange == 1) {
			Z = Z2;
		}
		if (z == 1 && offsetLayer1 < 1) {
			weight = (pow(Z,m+1)*Z2*Z2-pow((Z-offsetLayer1),m+1)
					  *(Z2-offsetLayer1)*(Z2-offsetLayer1))/(m+1.0);
			weight -= (pow(Z,m+2)*Z2-pow((Z-offsetLayer1),m+2)
					   *(Z2-offsetLayer1))*2/((m+1)*(m+2.0));
			weight += (pow(Z,m+3)-pow((Z-offsetLayer1),m+3))*2/((m+1)*(m+2)*(m+3.0));
		} else if (z == 2) {
			weight = (pow(Z,m+1)*Z2*Z2-pow((Z-layerAdjustment),m+1)
					  *(Z2-layerAdjustment)*(Z2-layerAdjustment))/(m+1.0);
			weight -= (pow(Z,m+2)*Z2-pow((Z-layerAdjustment),m+2)
					   *(Z2-layerAdjustment))*2/((m+1)*(m+2.0));
			weight += (pow(Z,m+3)-pow((Z-layerAdjustment),m+3))*2/((m+1)*(m+2)*(m+3.0));
		} else {
			weight = (pow(Z,m+1)*Z2*Z2-pow((Z-1),m+1)*(Z2-1)*(Z2-1))/(m+1.0);
			weight -= (pow(Z,m+2)*Z2-pow((Z-1),m+2)*(Z2-1))*2/((m+1)*(m+2.0));
			weight += (pow(Z,m+3)-pow((Z-1),m+3))*2/((m+1)*(m+2)*(m+3.0));
		}			
		value += copy[z]*weight*4*PI;
	}
	return value;
}
void
Lat1DSphere::SetLayerAdjustment(double value) {
	layerAdjustment = value;
	adjustOuterLayer = true;
	calculateLambdas();
}

void
Lat1DSphere::calculateLambdas() {
	Vector V(1,M);
	Vector S(1,M);
	S[1] = offsetLayer1*offsetLayer1*PI*4;
	V[1] = (4.0/3.0)*pow(offsetLayer1,3)*PI;
	L[1] = V[1];
	int z;
	double x;
	for (z=2; z<=M; z++) {
		x = z + offsetLayer1 - 2 + layerAdjustment;
		if (z == M-1 && adjustOuterLayer) {
			x += -layerAdjustment + 1;
		}
		double D = dimensions;
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
	if (L[1] > 4.0*PI/3.0) L[1] = V[1] - pow((offsetLayer1-1),3)*PI*4.0/3.0;
	if (L[1] > 0) lambda1[1] = L[2]*lambda_1[2]/L[1];
	if (L[1] > 0 && l1 != l1_chi) lambda1_chi[1] = L[2]*lambda_1_chi[2]/L[1];
	if (MyInput->ValueSet("lat",name,"n_layers_extended")) {
		int numLayersExtended = MyInput->GetInt("lat",name,
								"n_layers_extended",numLayers,INT_MAX);
		double numExtraLayers = numLayersExtended - numLayers;
		if (bound2 == M - 2) { // mirror2
			numExtraLayers += 0.5;
		}
		x = M + offsetLayer1 - 3 + layerAdjustment;
		double x2 = x + numExtraLayers;
		double D = dimensions;
		numExtraLatticeSites = (PI* pow(2,(D-1))/D) * pow(x2,D) 
							   - (PI* pow(2,(D-1))/D) * pow(x,D);;
	} else {
		numExtraLatticeSites = 0;
	}
	if (l1 == l1_chi) {
		lambda_1_chi = lambda_1;
		lambda1_chi = lambda1;
	}
}
