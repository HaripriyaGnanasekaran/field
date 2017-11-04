#include "LatticeRange2D.h"

LatticeRange2D::LatticeRange2D() {
}
LatticeRange2D::~LatticeRange2D(void) {
}
LatticeRange2D::LatticeRange2D(int lowX_,int highX_,int maxX_,int lowY_,int highY_, int maxY_) {
	lowX = lowX_;
	highX = highX_;
	maxX = maxX_;
	lowY = lowY_;
	highY = highY_;
	maxY = maxY_;
	if (lowX > highX) {
		Message(debug,"error in call to LatticeRange2D::LatticeRange2D"
			"lowX is greater than highX");
	}
	if (lowX < 1) {
		Message(debug,"error in call to LatticeRange2D::LatticeRange2D"
			"lowX is smaller than 1");
	}
	if (highX > maxX) {
		Message(debug,"error in call to LatticeRange2D::LatticeRange2D"
			"highX is greater than maxX");
	}
	if (lowY > highY) {
		Message(debug,"error in call to LatticeRange2D::LatticeRange2D"
			"lowY is greater than highY");
	}
	if (lowY < 1) {
		Message(debug,"error in call to LatticeRange2D::LatticeRange2D"
			"lowY is smaller than 1");
	}
	if (highY > maxY) {
		Message(debug,"error in call to LatticeRange2D::LatticeRange2D"
			"highY is greater than maxY");
	}
}
int
LatticeRange2D::Dimensions() const {
	return 2;
}
void
LatticeRange2D::MapVector(Vector out, Vector in) const {
	int x,y,z;

	for (x=1; x<lowX; x++) {
		z = (x-1)*maxY+1;
		for (y=1; y<=maxY; y++) {	out[z++] = 0;	}
	}
	for (x=lowX; x<=highX; x++) {
		z = (x-1)*maxY+1;
		for (y=1; y<lowY; y++) {
			out[z++] = 0;
		}
		for (y=lowY; y<=highY; y++) {
			out[z] = in[z];
			z++;
		}
		for (y=highY + 1; y<=maxY; y++) {
			out[z++] = 0;
		}
	}
	for (x=highX + 1; x<=maxX; x++) {
		z = (x-1)*maxY+1;
		for (y=1; y<=maxY; y++) {
			out[z++] = 0;
		}
	}
}

bool
LatticeRange2D::InRange(int z) const {
	int x = (z-1)/maxY + 1;
	if (x < lowX) return false;
	if (x > highX) return false;
	int y = z%maxY;
	if (y == 0) y = maxY;
	if (y < lowY) return false;
	if (y > highY) return false;
	return true;
}

double
LatticeRange2D::GetRangeValue(int z) const {
	if (InRange(z)) return 1.0;
	return 0.0;
}

int
LatticeRange2D::GetNumLayers() const {
	return (highX-lowX+1)*(highY-lowY+1);
}

int
LatticeRange2D::GetNumLayers_x() const {
	return highX-lowX+1;
}

int
LatticeRange2D::GetNumLayers_y() const {
	return highY-lowY+1;
}

int
LatticeRange2D::GetNumLayers_z() const {
	Message(fatal,"GetNumLayers z not used in 2D");
		return -1;
}

int
LatticeRange2D::GetNumPos()  {
	int aantal=0;
	//Message(fatal,"GetNumPos not used in 2D");
	return aantal;
}

double
LatticeRange2D::GetVolPos() {
	return 0.0;
}

int
LatticeRange2D::GetPos(int pos_number)  {
    int z=0;
    Message(fatal,"GetPos not used in 2D");
	return z;
}

bool
LatticeRange2D::SetPosLocal(int pos, double waarde)  { ;
    Message(fatal,"SetPosLocal not used in 2D");
	return false;
}

bool
LatticeRange2D::ClearAllPos()  {
	  Message(fatal,"ClearAllPos not used in 2D");
		return false;
}

bool
LatticeRange2D::SetPos(int pos, int SpotType, int SpotSize, double waarde)  { ;
    Message(fatal,"SetPos not used in 2D");
	return false;
}

bool
LatticeRange2D::SetPos(int pos, int SpotSize, double waarde)  { ;
    Message(fatal,"SetPos not used in 2D");
	return false;
}


bool
LatticeRange2D::SetPos(double x, double y, double z, double *submask)  { ;
    Message(fatal,"SetPos not used in 2D");
	return false;
}

bool
LatticeRange2D::ChangePos(int pos_new, int pos_old){
	Message(fatal,"ChangePos not used in 2D");
	return false;
}
int
LatticeRange2D::Get_jx() const {
	Message(fatal,"GetJx not used in 2D");
		return -1;
}

int
LatticeRange2D::Get_jy() const {
	Message(fatal,"GetJy not used in 2D");
		return -1;
}

int
LatticeRange2D::Get_jz() const {
	Message(fatal,"Get jz not used in 1D");
		return -1;
}

int
LatticeRange2D::Get_x0() const {
	return lowX;
}

int
LatticeRange2D::Get_y0() const {
		return lowY;
}

int
LatticeRange2D::Get_z0() const {
	Message(fatal,"Get z0 not used in 2D");
		return -1;
}

int
LatticeRange2D::Get_xm() const {
	return highX;
}

int
LatticeRange2D::Get_ym() const {
		return highY;
}

int
LatticeRange2D::Get_zm() const {
	Message(fatal,"Get zm not used in 1D");
		return -1;
}

int
LatticeRange2D::Get_M() const {
	return (highX-lowX+3)*(highY-lowY+3);
}

Text
LatticeRange2D::GetOutput() const {
	Text t1 = Blanks(1000);
	Text t2 = Blanks(1000);
	Text t3 = Blanks(1000);
	Text t4 = Blanks(1000);
	t1.Putint(lowX-1);
	t1 = Copy(t1.Frontstrip());
	t2.Putint(lowY-1);
	t2 = Copy(t2.Frontstrip());
	t3.Putint(highX-1);
	t3 = Copy(t3.Frontstrip());
	t4.Putint(highY-1);
	t4 = Copy(t4.Frontstrip());
	Text t = t1 + "," + t2 + ";" + t3 + "," + t4;
	return t;
}
