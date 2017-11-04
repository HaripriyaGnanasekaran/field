#include "LatticeRange1D.h"

LatticeRange1D::LatticeRange1D() {
}
LatticeRange1D::~LatticeRange1D() {
}
LatticeRange1D::LatticeRange1D(int lowX_,
							   int highX_,
							   int maxX_) {
	lowX = lowX_;
	highX = highX_;
	maxX = maxX_;
	if (lowX > highX) {
		Message(debug,"error in call to LatticeRange1D::LatticeRange1D"
			"lowX is greater than highX");
	}
	if (lowX < 1) {
		Message(debug,"error in call to LatticeRange1D::LatticeRange1D"
			"lowX is smaller than 1");
	}
	if (highX > maxX) {
		Message(debug,"error in call to LatticeRange1D::LatticeRange1D"
			"highX is greater than maxX");
	}
}
int
LatticeRange1D::Dimensions() const {
	return 1;
}
void
LatticeRange1D::MapVector(Vector out, Vector in) const {
	int z;
	for (z=1; z<lowX; z++) {
		out[z] = 0;
	}
	for (z=lowX; z<=highX; z++) {
		out[z] = in[z];
	}
	for (z=highX + 1; z<=maxX; z++) {
		out[z] = 0;
	}
}
bool
LatticeRange1D::InRange(int x) const {
	if (x < lowX) return false;
	if (x > highX) return false;
	else return true;
}

double
LatticeRange1D::GetRangeValue(int x) const {
	if (InRange(x)) return 1.0;
	return 0.0;
}

int
LatticeRange1D::GetNumLayers() const {
	return highX-lowX+1;
}

int
LatticeRange1D::GetNumLayers_x() const {
	return highX-lowX+1;
}


int
LatticeRange1D::GetNumPos() {
	int aantal=0;
	//Message(fatal,"GetNumPos not used in 1D");
	return aantal;
}

double
LatticeRange1D::GetVolPos() {
	return 0.0;
}

int
LatticeRange1D::GetPos(int pos_number) {
    int z=0;
    Message(fatal,"GetPos not used in 1D");
	return z;
}

bool
LatticeRange1D::SetPosLocal(int pos, double waarde) {
    Message(fatal,"SetPosLocal not used in 1D");
	return false;
}

bool
LatticeRange1D::ClearAllPos()  {
	  Message(fatal,"ClearAllPos not used in 1D");
		return false;
}

bool
LatticeRange1D::SetPos(int pos, int SpotType, int SpotSize, double waarde) {
    Message(fatal,"SetPos not used in 1D");
	return false;
}

bool
LatticeRange1D::SetPos(int pos,  int SpotSize, double waarde) {
    Message(fatal,"SetPos not used in 1D");
	return false;
}

bool
LatticeRange1D::SetPos(double x, double y, double z, double *submask) {
    Message(fatal,"SetPos not used in 1D");
	return false;
}


bool
LatticeRange1D::ChangePos(int pos_new, int pos_old){
	Message(fatal,"ChangePos not used in 1D");
	return false;
}

int
LatticeRange1D::GetNumLayers_y() const {
Message(fatal,"GetNumLayers y not used in 1D");
	return -1;
}

int
LatticeRange1D::GetNumLayers_z() const {
	Message(fatal,"GetNumLayers z not used in 1D");
		return -1;
}

int
LatticeRange1D::Get_jx() const {
	Message(fatal,"GetJx not used in 1D");
		return -1;
}

int
LatticeRange1D::Get_jy() const {
	Message(fatal,"GetJy not used in 1D");
		return -1;;
}

int
LatticeRange1D::Get_jz() const {
	Message(fatal,"GetJz not used in 1D");
		return -1;
}

int
LatticeRange1D::Get_x0() const {
	return lowX;
}

int
LatticeRange1D::Get_y0() const {
	Message(fatal,"Get y0 not used in 1D");
		return -1;
}

int
LatticeRange1D::Get_z0() const {
	Message(fatal,"Get z0 not used in 1D");
		return -1;
}

int
LatticeRange1D::Get_xm() const {
	return highX;
}

int
LatticeRange1D::Get_ym() const {
	Message(fatal,"Get ym not used in 1D");
		return -1;
}

int
LatticeRange1D::Get_zm() const {
	Message(fatal,"Get zm not used in 1D");
		return -1;
}

int
LatticeRange1D::Get_M() const {
	return (highX-lowX+3);
}

Text
LatticeRange1D::GetOutput() const {
	Text t1 = Blanks(1000);
	Text t2 = Blanks(1000);
	t1.Putint(lowX-1);
	t1 = Copy(t1.Frontstrip());
	t2.Putint(highX-1);
	t2 = Copy(t2.Frontstrip());
	Text t = t1 + ";" + t2;
	if (*Copy(t1) == *Copy(t2)) {
		return t1;
	}
	return t;
}
