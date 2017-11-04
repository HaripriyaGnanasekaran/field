#include "Lat2DCylinder.h"

// allowed keywords
const Text Lat2DCylinder::n_sites_first_layer_x="n_sites_first_layer_x";
const Text Lat2DCylinder::offset_first_layer_x="offset_first_layer_x";

Lat2DCylinder::Lat2DCylinder(Input* MyInput_, Text name_)
		: Lat2DFlat(MyInput_, name_)  {
	int z;
	double x, D;
	lambda1_1.Dim(1,MX);
	lambda11.Dim(1,MX);
	L.Dim(1,MX);
	Vector V(1,MX);
	Vector S(1,MX);

	checkParameters();
	MyInput->DontCombineParam(lat,name,offset_first_layer_x,
		n_sites_first_layer_x);
	if (dimensions != 3 && (boundX1 == 1 || boundX2 == MX)
		&& (bulkBoundY1 || bulkBoundY2)) { // Deze combinatie kan niet werken!
		Message(fatal, MyInput, "Please try this in three spatial dimensions!");
	}
	offsetLayer1 = MyInput->GetReal(lat,name,offset_first_layer_x,0,DBL_MAX,0);
	if (offsetLayer1 <= 0 && boundX1 == 1) {
		Message(fatal, Text("A calculation with a surface at the lowerbound ")
			+ "in the cylindrical lattice cannot make sense without setting "
			+ offset_first_layer_x + ".");
	}
	layerAdjustment = 1;
	if (MyInput->ValueSet(lat, name, n_sites_first_layer_x)) {
		double v = MyInput->GetReal(lat,name,n_sites_first_layer_x,1,DBL_MAX,PI);
		layerAdjustment = sqrt(v/PI);
	}
	S[1] = offsetLayer1*PI*2;
	V[1] = offsetLayer1*offsetLayer1*PI;
	L[1] = V[1];
	for (z=2; z<=MX; z++) {
		x = z + offsetLayer1 - 2 + layerAdjustment;
		D = dimensions - 1;
		V[z] = (PI* pow(2,(D-1))/D) * pow(x,D);
		L[z] = V[z] - V[z-1];
		S[z] = pow(2,(D-1))*PI*pow(x,(D-1));
		lambda1_1[z] = l1*S[z-1]/L[z];
		lambda11[z] = l1*S[z]/L[z];
	}
	// lambda1_1[1] = 0; These probabilities are used to calculate lambda0(z).
	lambda11[MX] = 0;
	if (L[1] > PI) {
		L[1] = V[1] - (offsetLayer1-1)*(offsetLayer1-1)*PI;
	}
	if (L[1] > 0) {
		lambda11[1] = L[2]*lambda1_1[2]/L[1];
	}
}
Lat2DCylinder::~Lat2DCylinder() {
}
// static protected function
void
Lat2DCylinder::checkParameters() const {
	Array<Text> param(1,18);
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
	param[17] = n_sites_first_layer_x;
	param[18] = offset_first_layer_x;
	MyInput->CheckParameterNames(lat,name,param);
}

void
Lat2DCylinder::GetOutput(Output* Out) const {
	Out->PutInt(lat,name,gradients,2);
	if (dimensions - int(dimensions) > 1e-7) {
		Out->PutReal(lat,name,"dimensions",dimensions);
	} else {
		Out->PutInt(lat,name,"dimensions",int(dimensions));
	}
	Out->PutText(lat,name,geometry,"flat");
	Out->PutInt(lat,name,n_layers_x,numLayersX);
	Out->PutInt(lat,name,n_layers_y,numLayersY);
	Out->PutInt(lat,name,"total number of lattice sites",
		(int)GetTotalNumLatticeSites());
	Array<Text> bounds(1,5);
	bounds[1] = "surface";
	bounds[2] = "mirror1";
	bounds[3] = "mirror2";
	bounds[4] = "periodic";
	bounds[5] = "bulk";
	MyInput->SetDefaultWarningsOff();
	// default mirror1
	int boundX1out = MyInput->GetChoice(lat, name, lowerbound_x, bounds, 2);
	int boundX2out = MyInput->GetChoice(lat, name, upperbound_x, bounds, 2);
	int boundY1out = MyInput->GetChoice(lat, name, lowerbound_y, bounds, 2);
	int boundY2out = MyInput->GetChoice(lat, name, upperbound_y, bounds, 2);
	MyInput->SetDefaultWarningsOn();
	Out->PutText(lat,name,lowerbound_x,bounds[boundX1out]);
	Out->PutText(lat,name,upperbound_x,bounds[boundX2out]);
	Out->PutText(lat,name,lowerbound_y,bounds[boundY1out]);
	Out->PutText(lat,name,upperbound_y,bounds[boundY2out]);
	outputLambdas(Out);
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
	Out->PutReal(lat,name,"volume (m)",GetVolume());
	Out->PutReal(lat,name,n_sites_first_layer_x,L[2]);
	MyInput->SetDefaultWarningsOff();
	Out->PutText(lat,name,offset_first_layer_x,MyInput->GetText(lat,name,
		offset_first_layer_x,"0"));
	MyInput->SetDefaultWarningsOn();
	Out->PutVector(lat,name,"lattice sites",L,MX);
	Out->PutVector(lat,name,"lambda1_1",lambda1_1,MX);
	Out->PutVector(lat,name,"lambda11",lambda11,MX);
}
void
Lat2DCylinder::outputLambdas(Output *Out) const {
	Out->PutReal(lat,name,lambda,l1);
}
double
Lat2DCylinder::GetNumLatticeSites(const int z) const {
	return L[(z-1)/MY+1];
}
double
Lat2DCylinder::GetLambda(const int z1, const int z2) const {
	int x = (z1-1)/MY+1;
	if (z1 == z2) return 1-lambda11[x]-lambda1_1[x]-2*l1;
	if (z1 - z2 == MY) return lambda1_1[x];
	if (z1 - z2 == -MY) return lambda11[x];
	if (z1 - z2 == 1) return l1;
	if (z1 - z2 == -1) return l1;
	return 0;
}
void
Lat2DCylinder::MultiplyWithLatticeSites(Vector A) const {
	int x,y,z=1;
	for (x=1; x<=MX; x++) {
		for (y=1; y<=MY; y++) {
			A[z] = A[z]*L[x];
			z++;
		}
	}
}
void
Lat2DCylinder::DivideByLatticeSites(Vector A) const {
	int x,y,z=1;
	for (x=1; x<=MX; x++) {
		for (y=1; y<=MY; y++) {
			if (L[x] != 0) {
				A[z] = A[z]/L[x];
			}
			z++;
		}
	}
}
double
Lat2DCylinder::GetVolume() const {
	return GetTotalNumLatticeSites()*pow(siteDistance,dimensions-1);
}
double
Lat2DCylinder::GetTotalNumLatticeSites() const {
	double valueX = 0;
	for(int x=1; x<=MX; x++) {
		valueX += L[x];
	}
	if (boundX1 != 1) valueX -= L[1];
	if (boundX1 == 3) valueX -= 0.5*L[2];
	if (boundX2 != MX) valueX -= L[MX];
	if (boundX2 == MX-2) valueX -= 0.5*L[MX-1];
	double valueY = MY;
	if (boundY1 != 1) valueY--;
	if (boundY1 == 3) valueY -= 0.5;
	if (boundY2 != MY) valueY--;
	if (boundY2 == MY-2) valueY -= 0.5;
	return valueX*valueY;
}

void
Lat2DCylinder::Sides(Vector P,Vector Q) const {
	Message(fatal,"Sides not implemented. ");
}

double
Lat2DCylinder::SideFraction(const Vector in, int z) const {

/* The input vector contains the data for which the sidefraction needs to be
 calculated. The input data ranges from 1 to MX*MY. */

/*It is not yet possible to have a bulkboundary in a x and y position. Positions
for where to put the code for this are given at the flag `bulkerror'. */
	assert(!((bulkBoundX1 || bulkBoundX2) && (bulkBoundY1 || bulkBoundY2)));
	int x = (z-1)/MY+1;
	double inz = in[z];
	// If z is in a corner...
	// Either return a value or give z a value not in a corner.
	if (z==1) {
		if (!bulkBoundX1 && boundX1 != 1) { // If xbound is mirror or periodic..
			z = 1+(boundX1-1)*MY;
		} else if (!bulkBoundY1 && boundY1 != 1) { // Same for yborder.
			z = boundY1;
		} else if (boundX1 == 1 && boundY1 == 1) { //If pos. 1 belongs to system
			return getLambda0(1)*inz + l1*in[2] + lambda11[1]*in[MY+1];
		} else if (boundX1 == 1) { // bulkBoundY1 == true
			return inz*(getLambda0(1) + l1 + l1) + in[MY+1]*lambda11[1];
		} else if (boundY1 == 1) { // bulkBoundX1 == true
			return (lambda11[1] + getLambda0(1) + lambda1_1[1])*inz + l1*in[2];
		} else { // bulkBoundX1 == true && bulkBoundY1 == true
			// bulkerror
		}
	} else if (z==MY) {
		if (!bulkBoundX1 && boundX1 != 1) {
			z = boundX1*MY;
		} else if (!bulkBoundY2 && boundY2 != MY) {
			z = boundY2;
		} else if (boundX1 == 1 && boundY2 == MY) {
			return getLambda0(z)*inz + l1*in[z-1] + lambda11[1]*in[z+MY];
		} else if (boundX1 == 1) {
			return (getLambda0(1) + l1 + l1)*inz + lambda11[1]*in[MY+MY];
		} else if (boundY2 == MY) {
			return (lambda11[1] + getLambda0(1) + lambda1_1[1])*inz +l1*in[z-1];
		} else {
			// bulkerror
		}
	} else if (z==M+1-MY) {
		if (!bulkBoundX2 && boundX2 != MX) {
			z = 1+(boundX2-1)*MY;
		} else if (!bulkBoundY1 && boundY1 != 1) {
			z = M-MY+boundY1;
		} else if (boundX2 == MX && boundY1 == 1) {
			return getLambda0(MX)*inz + l1*in[z+1] + lambda1_1[MX]*in[z-MY];
		} else if (boundX2 == MX) {
			return (getLambda0(MX) + l1 + l1)*inz + lambda1_1[MX]*in[z-MY];
		} else if (boundY1 == 1) {
			return (lambda11[MX] + getLambda0(MX) + lambda1_1[MX])*inz
				+ l1*in[z+1];
		} else {
			// bulkerror
		}
	} else if (z==M) {
		if (!bulkBoundX2 && boundX2 != MX) {
			z = boundX2*MY;
		} else if (!bulkBoundY2 && boundY2 != MY) {
			z = M-MY+boundY2;
		} else if (boundX2 == MX && boundY2 == MY) {
			return getLambda0(MX)*inz + l1*in[z-1] + lambda1_1[MX]*in[z-MY];
		} else if (boundX2 == MX) {
			return (getLambda0(MX)+l1+l1)*inz + lambda1_1[MX]*in[z-MY];
		} else if (boundY2 == MY) {
			return (lambda11[MX] + getLambda0(MX) + lambda1_1[MX])*inz
				+ l1*in[z-1];
		} else {
			// bulkerror
		}
	}
	// If z is in boundary...
	// Either return a value or give z a value not in a boundary.
	if (z%MY == 1) { // If in left boundary...
		if (bulkBoundY1) { // This goes wrong in case of bulkerror.
			return inz + lambda11[x]*(in[z+MY]-inz)
				+ lambda1_1[x]*(in[z-MY]-inz);
		} else if (boundY1 == 1) {
			return getLambda0(x)*inz + lambda11[x]*in[z+MY]
				+ lambda1_1[x]*in[z-MY] + l1*in[z+1];
		} else { // If boundary is periodical or mirror...
			z += boundY1-1;
		}
	} else if (z%MY == 0) { // If in right boundary...
		if (bulkBoundY2) {
			return in[z] + lambda11[x]*(in[z+MY]-inz)
				+ lambda1_1[x]*(in[z-MY]-inz);
		} else if (boundY2 == MY) {
			return getLambda0(x)*inz + lambda11[x]*in[z+MY]
				+ lambda1_1[x]*in[z-MY] + l1*in[z-1];
		} else {
			z += boundY2-MY;
		}
	} else if (z <= MY) { // If in lower boundary...
		if (bulkBoundX1) {
			return inz + l1*(in[z-1]+in[z+1]-2*inz);
		} else if (boundX1 == 1) {
			return getLambda0(1)*inz + l1*(in[z-1]+in[z+1])
				+ lambda11[1]*in[z+MY];
		} else {
			z += (boundX1-1)*MY;
		}
	} else if (z > M-MY) { // If in upper boundary...
		if (bulkBoundX2) {
			return inz + l1*(in[z-1]+in[z+1]-2*inz);
		} else if (boundX2 == MX) {
			return getLambda0(MX)*inz + l1*(in[z-1]+in[z+1])
				+ lambda1_1[MX]*in[z-MY];
		} else {
			z += boundX2*MY - M;
		}
	}

	// Now z is not in a boundary, we can calculate the sidefraction.
	return inz + lambda11[x]*(in[z+MY]-inz) + lambda1_1[x]*(in[z-MY]-inz)
		+ l1*(in[z-1] + in[z+1] - inz - inz);
}
void
Lat2DCylinder::ElectricPotential(Vector &psi,
		Vector &psi0,
		const Vector &eps,
		const Vector &charge,
		const double preFactor) const {
	int x,y,z;
	double Z, pf, xminpf, xpluspf;
	double epsxmin, epsxplus, epsymin, epsyplus;

	// bulkboundarie support is not yet implemented!
	// the charge is for the entire layer, not 1 lattice site
	SetBoundaries(psi0);

	// lower x and its corners
	Z = L[3]/2/PI-2; // radius to center of layer 1
	if (offsetLayer1 <= 0) {
		Z++;
		pf = preFactor/PI/Z;
		xpluspf = 1+0.5/Z;
		z = MY+MY;
		epsymin = eps[z]+eps[z-1];
		if (boundY2 == MY) {
			epsxplus = xpluspf*(eps[z] + eps[z+MY]);
			psi[z] = (pf*charge[z] + epsymin*psi0[z-1]
				+ epsxplus*psi0[z+MY])/(epsymin+epsxplus);
		}
		for (z=2*MY-1; z>MY+1; z--) {
			epsyplus = epsymin;
			epsymin = eps[z]+eps[z-1];
			epsxplus = xpluspf*(eps[z] + eps[z+MY]);
			psi[z] = (pf*charge[z] + epsxplus*psi0[z+MY]
				+ epsymin*psi0[z-1] + epsyplus*psi0[z+1])
				/(epsxplus + epsymin + epsyplus);
		}
		if (boundY1 == 1) {
			epsxplus = xpluspf*(eps[z] + eps[z+MY]);
			psi[z] = (pf*charge[z] + epsxplus*psi0[z+MY]
				+ epsymin*psi0[z+1])/(epsxplus+epsymin);
		}
		x = 3;
		z = MY+MY+1;
	} else {
		if (boundX1 == 1) {
			pf = preFactor/PI/Z;
			xminpf = 1-0.5/Z;
			xpluspf = 1+0.5/Z;
			epsymin = eps[MY]+eps[MY-1];
			if (boundY2 == MY) {
				epsxplus = xpluspf*(eps[MY] + eps[MY+MY]);
				psi[MY] = (pf*charge[MY] + epsymin*psi0[MY-1]
					+ epsxplus*psi0[MY+MY])/(epsymin+epsxplus);
			}
			for (z=MY-1; z>1; z--) {
				epsyplus = epsymin;
				epsymin = eps[z]+eps[z-1];
				epsxplus = xpluspf*(eps[z] + eps[z+MY]);
				psi[z] = (pf*charge[z] + epsxplus*psi0[z+MY]
					+ epsymin*psi0[z-1] + epsyplus*psi0[z+1])
					/(epsxplus + epsymin + epsyplus);
			}
			if (boundY1 == 1) {
				epsxplus = xpluspf*(eps[1] + eps[MY+1]);
				psi[1] = (pf*charge[1] + epsxplus*psi0[MY+1]
					+ epsymin*psi0[2])/(epsxplus+epsymin);
			}
		}
		x = 2;
		z = MY+1;
	}
	// middle of lattice and lower and upper y
	for (; x<MX; x++) {
		Z++;
		pf = preFactor/PI/Z;
		xminpf = 1-0.5/Z;
		xpluspf = 1+0.5/Z;
		epsyplus = eps[z]+eps[z+1];
		if (boundY1==1) {
			epsxmin = xminpf*(eps[z-MY] + eps[z]);
			epsxplus = xpluspf*(eps[z] + eps[z+MY]);
			psi[z] = (pf*charge[z] + epsxmin*psi0[z-MY] + epsxplus*psi0[z+MY]
				+ epsyplus*psi0[z+1])
				/(epsxmin + epsxplus + epsyplus);
		}
		z++;
		for (y=2; y<MY; y++) {
			epsymin = epsyplus;
			epsyplus = eps[z]+eps[z+1];
			epsxmin = xminpf*(eps[z-MY] + eps[z]);
			epsxplus = xpluspf*(eps[z] + eps[z+MY]);
			psi[z] = (pf*charge[z] + epsxmin*psi0[z-MY] + epsxplus*psi0[z+MY]
				+ epsymin*psi0[z-1] + epsyplus*psi0[z+1])
				/(epsxmin + epsxplus + epsymin + epsyplus);
			z++;
		}
		if (boundY2==MY) {
			epsxmin = xminpf*(eps[z-MY] + eps[z]);
			epsxplus = xpluspf*(eps[z] + eps[z+MY]);
			psi[z] = (pf*charge[z] + epsxmin*psi0[z-MY] + epsxplus*psi0[z+MY]
				+ epsyplus*psi0[z-1]) / (epsxmin + epsxplus + epsyplus);
		}
		z++;
	}
	// upper x and its corners
	if (boundX2 == MX) {
		Z++;
		pf = preFactor/PI/Z;
		xminpf = 1-0.5/Z;
		epsyplus = eps[M-MY+1]+eps[M-MY+2];
		if (boundY1 == 1) {
			epsxmin = xminpf*(eps[M-MY-MY+1]+eps[M-MY+1]);
			psi[M-MY+1] = (pf*charge[M-MY+1] + epsxmin*psi0[M-MY-MY+1]
				+epsyplus*psi0[M-MY+2])/(epsxmin + epsyplus);
		}
		for (z=M-MY+2; z<M; z++) {
			epsymin = epsyplus;
			epsyplus = eps[z]+eps[z+1];
			epsxmin = xminpf*(eps[z-MY] + eps[z]);
			psi[z] = (pf*charge[z] + epsxmin*psi0[z-MY]
				+ epsymin*psi0[z-1] + epsyplus*psi0[z+1])
				/(epsxmin + epsymin + epsyplus);
		}
		if (boundY2 == MY) {
			epsxmin = xminpf*(eps[M]+eps[M-MY]);
			psi[M] = (pf*charge[M] + epsxmin*psi0[M-MY]
				+ epsyplus*psi0[M-1])/(epsxmin + epsyplus);
		}
	}
}
void
Lat2DCylinder::ElectricFieldSquared(Vector &ESquared,
		const Vector psi,
		const double preFactor) const {
	int x, y, z;
	double pf, Z, Exmin, Explus, Eymin, Eyplus, xminpf, xpluspf;

	pf = preFactor/2;
	Z = L[3]/2/PI-2; // radius to center of layer 1
	if (offsetLayer1 <= 0) {
		Z++;
		xminpf = 1-0.5/Z;
		xpluspf = 1+0.5/Z;
		Eymin = psi[MY+MY]-psi[MY+MY-1];
		Eymin *= Eymin;
		for (z=2*MY-1; z>MY+1; z--) {
			Eyplus = Eymin;
			Eymin = psi[z]-psi[z-1];
			Eymin *= Eymin;
			Explus = psi[z]-psi[z+MY];
			Explus *= xpluspf*Explus;
			ESquared[z] = pf*(Explus + Eymin + Eyplus);
		}
		x = 3;
		z = MY+MY+1;
	} else {
		if (boundX1 == 1) {
			xminpf = 1-0.5/Z;
			xpluspf = 1+0.5/Z;
			Eymin = psi[MY]-psi[MY-1];
			Eymin *= Eymin;
			for (z=MY-1; z>1; z--) {
				Eyplus = Eymin;
				Eymin = psi[z]-psi[z-1];
				Eymin *= Eymin;
				Explus = psi[z]-psi[z+MY];
				Explus *= xpluspf*Explus;
				ESquared[z] = pf*(Explus + Eymin + Eyplus);
			}
		}
		x = 2;
		z = MY+1;
	}
	for (; x<MX; x++) {
		Z++;
		xminpf = 1-0.5/Z;
		xpluspf = 1+0.5/Z;
		Eyplus = psi[z]-psi[z+1];
		Eyplus *= Eyplus;
		if (boundY1 == 1) {
			Exmin = psi[z]-psi[z-MY];
			Exmin *= xminpf*Exmin;
			Explus = psi[z]-psi[z+MY];
			Explus *= xpluspf*Explus;
			ESquared[z] = pf*(Exmin + Explus + Eyplus);
		}
		z++;
		for (y=2; y<MY; y++) {
			Eymin = Eyplus;
			Eyplus = psi[z]-psi[z+1];
			Eyplus *= Eyplus;
			Exmin = psi[z]-psi[z-MY];
			Exmin *= xminpf*Exmin;
			Explus = psi[z]-psi[z+MY];
			Explus *= xpluspf*Explus;
			ESquared[z] = pf*(Exmin + Explus + Eymin + Eyplus);
			z++;
		}
		if (boundY2 == MY) {
			Exmin = psi[z]-psi[z-MY];
			Exmin *= xminpf*Exmin;
			Explus = psi[z]-psi[z+MY];
			Explus *= xpluspf*Explus;
			ESquared[z] = pf*(Exmin + Explus + Eyplus);
		}
		z++;
	}
	if (boundX2 == MX) {
		Z++;
		xminpf = 1-0.5/Z;
		xpluspf = 1+0.5/Z;
		Eyplus = psi[M-MY+2]-psi[M-MY+1];
		Eyplus *= Eyplus;
		for (z=M-MY+2; z<M; z++) {
			Eymin = Eyplus;
			Eyplus = psi[z]-psi[z+1];
			Eyplus *= Eyplus;
			Exmin = psi[z]-psi[z-MY];
			Exmin *= xminpf*Exmin;
			ESquared[z] = pf*(Exmin + Eymin + Eyplus);
		}
	}
}
Vector
Lat2DCylinder::Div1Grad2(const Vector in1, const Vector in2) const {
	Vector value(1,MX*MY);
	Message(implementation,"Lat2DCylinder::Div1Grad2");
	for (int z=1; z<=MX*MY; z++) {
		double Z = ((z -1)/MY+1) + offsetLayer1;
		if (z%MY != 1 && in1[z] != 0 && in1[z-1] != 0) {
			value[z] += (in1[z-1]+in1[z])*(in2[z-1]-in2[z]);
		}
		if (z%MY != 0 && in1[z] != 0 && in1[z+1] != 0) {
			value[z] += (in1[z+1]+in1[z])*(in2[z+1]-in2[z]);
		}
		if (z > MY && in1[z] != 0 && in1[z-MY] != 0) {
		value[z] += ((Z-1)*in1[z-MY]+Z*in1[z])*(in2[z-MY]-in2[z]);
		}
		if (z < (MX - 1) * MY && in1[z] != 0 && in1[z+MY] != 0) {
			value[z] += ((Z+1)*in1[z+MY]+Z*in1[z])*(in2[z+MY]-in2[z]);
		}
		value[z] /= 4*siteDistance*siteDistance*Z;
	}
	return value;
}
Vector
Lat2DCylinder::FourierTransform(const Vector in) const {
	Message(implementation,"Lat2DCylinder::FourierTransform(Vector in)");
	return in;
}
double
Lat2DCylinder::Moment(const Vector in, const int moment,
	const LatticeRange* Range, const double g) const {

	// if Range is NULL, get moment in a specific gradient (see below)
	if (Range == NULL) {
		return Moment(in, moment, (int)g);
	}
	
	// function return unnormalized moment

	/* since the moment in a space with D>1 is a vector, we have a problem
	in 2D classes with this function, since it returns a double. Fortunately in
	this cylindrical lattice all momenta are along the z-axis. So we can
	just return the length of the vector. */

	/* the offset supplied with in this function, thus, is the offset along the
	z-axis. */

	/* only the 0th, 1st and 2nd moment have been implemented sofar. */

	if (moment < 0 || moment > 2) {
		Message(implementation,
			Text("Lat2DCylinder::Moment: only moment 0, 1, 2"));
		return sqrt(-1.0);
	}

	if (dimensions != 3) {
		Message(fatal, Text("Sorry, no implemented for D !=3."));
	}

	/* in this lattice there are two possibilities:
	1 - either all layers in the r-direction are equally wide, but the first
		layer doesn't start at the center of the lattice
 	2 - the first layer has a differt width in the r-direction
	if offsetLayer1 == 0 we have situation 2 */

	/* I (jos) don't understand the reason for the latticerange parameter, so
	I'm ignoring it for now. */

	double value = 0;

	Vector copy(1,M);
	int z,x,y;
	for (z=1; z<=M; z++) {
		copy[z] = in[z];
	}
	SubtractBoundaries(copy);
	if (moment == 0) {
		MultiplyWithLatticeSites(copy);
		for (z=1; z<=M; z++) {
			value += copy[z];
		}
		return value;
	}
	if (moment == 1) {
		MultiplyWithLatticeSites(copy);
		// first project `copy' on the z-axis
		Vector alongz(1,MY);
		for (x=MX; x>0; x--) {
			for (y=MY; y>0; y--) {
				alongz[y] += copy[(x-1)*MY+y];
			}
		}
		// determine the moment along the z-axis
		for (y=MY; y>0; y--) {
			value += alongz[y]*(y*y-(y-1)*(y-1));
		}
		return value/2;
	}
	double zoff = Moment(in, 1, Range, -0.5)/Moment(in, 0, Range, -0.5);
	double roff;
	// if you're here: moment == 2
	if (offsetLayer1 == 0) {
		roff = layerAdjustment - 2;
		for (x=MX-1; x>2; x--) { // all layers except x==1
			for (y=MY-1; y>1; y--) {
				value += copy[(x-1)*MY+y] *
				( -zoff*zoff/2 + 4/3.0*roff + y/2.0 - 3*roff*roff/2 + 4/3.0*x
				+ zoff*y - 2*zoff*x*y - zoff/2 + 3*roff*x*x - 3/2.0*x*x
				- 3*roff*x + 3*roff*roff*x + zoff*zoff*roff + zoff*zoff*x
				+ roff*roff*roff + x*x*x - 2*y*zoff*roff - y*y/2.0 - 5/12.0
				+ zoff*x + zoff*roff - x*y + y*y*roff + x*y*y - roff*y);
			}
		}
		for (y=MY-1; y>1; y--) {
			value += copy[MY+y] * layerAdjustment * layerAdjustment *
			( layerAdjustment*layerAdjustment/4 + y*y/2.0 - y/2.0 + 1/6.0 - zoff*y + zoff/2
			+ zoff*zoff/2);
		}
		return value*2*PI;
	}
	// if you're here: moment == 2 && offsetLayer1 != 0
	roff = offsetLayer1 - 1;
	for (x=MX-1; x>1; x--) { // all layers except x==1
		for (y=MY-1; y>1; y--) {
			value += copy[(x-1)*MY+y] *
			( -zoff*zoff/2 + 4/3.0*roff + y/2.0 - 3*roff*roff/2 + 4/3.0*x
			+ zoff*y - 2*zoff*x*y - zoff/2 + 3*roff*x*x - 3/2.0*x*x - 3*roff*x
			+ 3*roff*roff*x + zoff*zoff*roff + zoff*zoff*x + roff*roff*roff
			+ x*x*x - 2*y*zoff*roff - y*y/2.0 - 5/12.0 + zoff*x + zoff*roff
			- x*y + y*y*roff + x*y*y - roff*y);
		}
	}
	return value*2*PI;
}
/**
 * This function returns the moment of a vector in a specified gradient.
 * If the gradient = 0, then no gradient is specified and the normal moment is
 * returned.
 * If the gradient = 1, the moment of the gradient along the cylinder axis is
 * returned.
 * If the gradient = 2, the gradient perpendicular to the cylinder axis is
 * returned.
 * The ran and offs arguments are ignored.
 * The edges of the system are not taken into account.
**/
double
Lat2DCylinder::Moment(const Vector in, const int m,
	const int gradient) const {
	int x,y,z;
	double off,value;
	Vector copy, // copy of incoming vector
		proj; // projection of in on specified gradient

	if (gradient != 1 && gradient != 2) {
		Message(fatal, Text("Gradient outside range in "
			"Lat2DCylinder::Moment."));
	}

	copy.Dim(1,M);
	for (z=1; z<=M; z++) {
		copy[z] = in[z];
	}
	SubtractBoundaries(copy);

	if (gradient==2 || m==0) { // gradient in z-direction
		MultiplyWithLatticeSites(copy);
		proj.Dim(1,MY);
		value=0;
		off=0;
		for (y=MY; y>0; y--) {
			proj[y]=0;
			for (x=MX; x>0; x--) {
				proj[y] += copy[(x-1)*MY+y];
			}
			off += proj[y]*(pow(y-1.,2)-pow(y-2.,2));
			value += proj[y];
		}
		if (m==0) {
			return value;
		}
		off /= 2;
		if (m==1) {
			return off;
		}
		off /= value;
		off += 1;
		value = 0;
		for (y=MY; y>0; y--) {
			value += proj[y]*(pow(y-off,m+1)-pow(y-off-1,m+1));
		}
		value /= (m+1);
	} else { // gradient in r-direction
		proj.Dim(1,MX);
		for (x=MX; x>0; x--) {
			proj[x]=0;
			for (y=MY; y>0; y--) {
				proj[x] += copy[(x-1)*MY+y];
			}
		}
		if (offsetLayer1 == 0) { //layer 1 may have a different width
			off = layerAdjustment-2;
			value = proj[2]*pow(off+2,m+2);
		} else { // layer 1 starts at layerAdjustment-2, omit layer 1
			off = offsetLayer1-1;
			value = proj[1]*(pow(off+1,m+2)-pow(off,m+2))
				+ proj[2]*(pow(off+2,m+2)-pow(off+1,m+2));
		}
		for (x=MX; x>2; x--) {
			value += proj[x]*(pow(off+x,m+2)-pow(off+x-1,m+2));
		}
		value *= 2*PI/(m+2);
	}
	return value;
}
Vector
Lat2DCylinder::RenormPhi(Vector phiRenorm,
						 Vector phiBound,
						 const Vector phiRestricted,
						 const double thetaRenorm) const {
	int x,y,z;
	double theta = 0;
	for (z=1; z<=MX*MY; z++) {
		theta += phiRestricted[z];
	}
	SubtractBoundaries(phiRestricted);
	Vector Renorm(1,MY);
	z=0;
	for (y=1; y<=MY; y++) {
		Renorm[y] = 0;
		for (x=1; x<=MX; x++) {
			z++;
			Renorm[y] += phiRestricted[z]*L[x];
		}
		Renorm[y] *= thetaRenorm/phiBound[y+MY];
	}
	RestoreBoundaries(phiRestricted);
	// bogus function !!
	SubtractBoundaries(phiRenorm);

	return Renorm;
}
double
Lat2DCylinder::getLambda0(int x) const {
	return 1 - 2 * l1 - lambda11[x] - lambda1_1[x];
}
