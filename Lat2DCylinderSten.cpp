#include "Lat2DCylinderSten.h"

Lat2DCylinderSten::Lat2DCylinderSten(Input *MyInput_, Text name_)
		: Lat2DCylinder(MyInput_, name_) {
	int z;
	double x, D;

	Vector V(1,MX);
	Vector S(1,MX);
	lambda21.Dim(1,MX);
	lambda2_1.Dim(1,MX);

//	layerAdjustment = sqrt(volumeFirstLayer/PI);
	S[1] = offsetLayer1*PI*2;
	V[1] = offsetLayer1*offsetLayer1*PI;
	L[1] = V[1];
	// TODO: define lambda for z=1
	for (z=2; z<=MX; z++) {
		x = z + offsetLayer1 - 2 + layerAdjustment;
		D = dimensions - 1;
		V[z] = (PI* pow(2,(D-1))/D) * pow(x,D);
		L[z] = V[z] - V[z-1];
		S[z] = pow(2,(D-1))*PI*pow(x,(D-1));
		lambda2_1[z] = l2*S[z-1]/L[z];
		lambda21[z] = l2*S[z]/L[z];
	}
	if (L[1] > PI) {
		L[1] = V[1] - (offsetLayer1-1)*(offsetLayer1-1)*PI;
	}
	if (L[1] > 0) {
		lambda21[1] = L[2]*lambda2_1[2]/L[1];
	}
}

Lat2DCylinderSten::~Lat2DCylinderSten() {
}
void
Lat2DCylinderSten::GetOutput(Output* Out) const {
	Lat2DCylinder::GetOutput(Out);
	Out->PutVector(lat,name,"lambda2_1",lambda2_1,MX);
	Out->PutVector(lat,name,"lambda21",lambda21,MX);
}
void
Lat2DCylinderSten::outputLambdas(Output *Out) const {
	if (MyInput->ValueSet(lat,name,lambda)) {
		Out->PutReal(lat,name,lambda,l1);
	} else {
		Out->PutReal(lat, name, lambda0, l0);
		Out->PutReal(lat, name, lambda1, l1);
		Out->PutReal(lat, name, lambda2, l2);
	}
}
double
Lat2DCylinderSten::GetLambda(const int z1, const int z2) const {
	int x = (z1-1)/MY+1;
	int d = z1 - z2;
	if (!d) {
		return getLambda0(x);
	} else if (d == 1 || d == -1) {
		return l1;
	} else if (d == MY) {
		return lambda11[x];
	} else if (d == -MY) {
		return lambda1_1[x];
	} else if (d == MY+1 || d == MY-1) {
		return lambda21[x];
	} else if (d == -MY+1 || d == -MY-1) {
		return lambda2_1[x];
	}
	return 0;
}

void
Lat2DCylinderSten::Sides(Vector P,Vector Q) const {
	Message(fatal,"Sides not implemented. ");
}

double
Lat2DCylinderSten::SideFraction(const Vector in, int z) const {
	// The input vector contains the data for which the sidefraction needs to be calculated.
	// The input data ranges from 1 to MX*MY

	// It is not yet possible to have a bulkboundary in a x and y position.
	// Positions for where to put
	// the code for this are given at the flag `bulkerror'.
	assert(!((bulkBoundX1 || bulkBoundX2) && (bulkBoundY1 || bulkBoundY2)));
	int x = (z-1)/MY+1;
	double inz = in[z];
	// If z is in a corner...
	// Either return a value or give z a value not in a corner.
	if (z==1) {
		if (!bulkBoundX1 && boundX1 != 1) { // If xbound is mirror or periodic...
			z = 1+(boundX1-1)*MY;
		} else if (!bulkBoundY1 && boundY1 != 1) { // Same for yborder.
			z = boundY1;
		} else if (boundX1 == 1 && boundY1 == 1) { // If pos. 1 belongs to system...
			return getLambda0(1)*inz + l1*in[2] + lambda11[1]*in[MY+1]
				+ lambda21[MY+2];
		} else if (boundX1 == 1) { // bulkBoundY1 == true
			return inz*(getLambda0(1) + l1 + l1)
				+ in[MY+1]*(2*lambda21[1] + lambda11[1]);
		} else if (boundY1 == 1) { // bulkBoundX1 == true
			return (lambda11[1] + getLambda0(1) + lambda1_1[1])*inz
				+ (lambda21[1] + l1 + lambda2_1[1])*in[2];
		} else { // bulkBoundX1 == true && bulkBoundY1 == true
			// bulkerror
		}
	} else if (z==MY) {
		if (!bulkBoundX1 && boundX1 != 1) {
			z = boundX1*MY;
		} else if (!bulkBoundY2 && boundY2 != MY) {
			z = boundY2;
		} else if (boundX1 == 1 && boundY2 == MY) {
			return getLambda0(z)*inz + l1*in[z-1] + lambda11[1]*in[z+MY]
				+ lambda2_1[1]*in[z+MY-1];
		} else if (boundX1 == 1) {
			return (getLambda0(1) + l1 + l1)*inz
				+ (2*lambda21[1]+lambda11[1])*in[MY+MY];
		} else if (boundY2 == MY) {
			return (lambda11[1] + getLambda0(1) + lambda1_1[1])*inz
				+ (lambda21[1] + l1 + lambda2_1[1])*in[z-1];
		} else {
			// bulkerror
		}
	} else if (z==M+1-MY) {
		if (!bulkBoundX2 && boundX2 != MX) {
			z = 1+(boundX2-1)*MY;
		} else if (!bulkBoundY1 && boundY1 != 1) {
			z = M-MY+boundY1;
		} else if (boundX2 == MX && boundY1 == 1) {
			return getLambda0(MX)*inz + l1*in[z+1]
				+ lambda1_1[MX]*in[z-MY] + lambda2_1[MX]*in[z-MY+1];
		} else if (boundX2 == MX) {
			return (getLambda0(MX) + l1 + l1)*inz
				+ (lambda1_1[MX] + 2*lambda2_1[MX])*in[z-MY];
		} else if (boundY1 == 1) {
			return (lambda11[MX] + getLambda0(MX) + lambda1_1[MX])*inz
				+ (lambda21[MX] + l1 + lambda2_1[MX])*in[z+1];
		} else {
			// bulkerror
		}
	} else if (z==M) {
		if (!bulkBoundX2 && boundX2 != MX) {
			z = boundX2*MY;
		} else if (!bulkBoundY2 && boundY2 != MY) {
			z = M-MY+boundY2;
		} else if (boundX2 == MX && boundY2 == MY) {
			return getLambda0(MX)*inz + l1*in[z-1]
				+ lambda1_1[MX]*in[z-MY] + lambda2_1[MX]*in[z-MY-1];
		} else if (boundX2 == MX) {
			return (getLambda0(MX)+l1+l1)*inz
				+ (lambda1_1[MX] + 2*lambda2_1[MX])*in[z-MY];
		} else if (boundY2 == MY) {
			return (lambda11[MX] + getLambda0(MX) + lambda1_1[MX])*inz
				+ (lambda21[MX] + l1 + lambda2_1[MX])*in[z-1];
		} else {
			// bulkerror
		}
	}
	// If z is in boundary...
	// Either return a value or give z a value not in a boundary.
	if (z%MY == 1) { // If in left boundary...
		if (bulkBoundY1) { // This goes wrong in case of bulkerror.
			return inz + (lambda11[x]+2*lambda21[x])*(in[z+MY]-inz)
				+ (lambda1_1[x]+2*lambda2_1[x])*(in[z-MY]-inz);
		} else if (boundY1 == 1) {
			return getLambda0(x)*inz + lambda11[x]*in[z+MY]
				+ lambda1_1[x]*in[z-MY] + l1*in[z+1] + lambda21[x]*in[z+MY+1]
				+ lambda2_1[x]*in[z-MY+1];
		} else { // If boundary is periodical or mirror...
			z += boundY1-1;
		}
	} else if (z%MY == 0) { // If in right boundary...
		if (bulkBoundY2) {
			return in[z] + (lambda11[x]+2*lambda21[x])*(in[z+MY]-inz)
				+ (lambda1_1[x]+2*lambda2_1[x])*(in[z-MY]-inz);
		} else if (boundY2 == MY) {
			return getLambda0(x)*inz + lambda11[x]*in[z+MY]
				+ lambda1_1[x]*in[z-MY] + l1*in[z-1] + lambda21[x]*in[z+MY-1]
				+ lambda2_1[x]*in[z-MY-1];
		} else {
			z += boundY2-MY;
		}
	} else if (z <= MY) { // If in lower boundary...
		if (bulkBoundX1) {
			return inz
				+ (lambda21[1] + l1 + lambda2_1[1])*(in[z-1]+in[z+1]-2*inz);
		} else if (boundX1 == 1) {
			return getLambda0(1)*inz + l1*(in[z-1]+in[z+1])
				+ lambda11[1]*in[z+MY] + lambda21[1]*(in[z+MY-1]+in[z+MY+1]);
		} else {
			z += (boundX1-1)*MY;
		}
	} else if (z > M-MY) { // If in upper boundary...
		if (bulkBoundX2) {
			return inz
				+ (lambda21[MX] + l1 + lambda2_1[MX])*(in[z-1]+in[z+1]-2*inz);
		} else if (boundX2 == MX) {
			return getLambda0(MX)*inz + l1*(in[z-1]+in[z+1])
				+ lambda1_1[MX]*in[z-MY]
				+ lambda2_1[MX]*(in[z-MY-1]+in[z-MY+1]);
		} else {
			z += boundX2*MY - M;
		}
	}

	// Now z is not in a boundary, we can calculate the sidefraction.
	return inz + lambda11[x]*(in[z+MY]-inz) + lambda1_1[x]*(in[z-MY]-inz)
		+ l1*(in[z-1] + in[z+1] - inz - inz)
		+ lambda21[x]*(in[z+MY+1] + in[z+MY-1] - inz - inz)
		+ lambda2_1[x]*(in[z-MY+1] + in[z-MY-1] - inz - inz);
}
double
Lat2DCylinderSten::getLambda0(int x) const {
	return 1 - 2 * l1 - lambda11[x] - lambda1_1[x] - 2 * lambda21[x]
		- 2 * lambda2_1[x];
}
