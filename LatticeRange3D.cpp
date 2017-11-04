#include "LatticeRange3D.h"
#include <iostream>
#include <ostream>
#include <sstream>

LatticeRange3D::LatticeRange3D() {
}
LatticeRange3D::~LatticeRange3D(void) {
}
LatticeRange3D::LatticeRange3D(int lowX_, int highX_, int maxX_, int lowY_,  int highY_,  int maxY_, int lowZ_,  int highZ_,  int maxZ_) {
	lowX = lowX_;
	highX = highX_;
	maxX = maxX_;
	lowY = lowY_;
	highY = highY_;
	maxY = maxY_;
	lowZ = lowZ_;
	highZ = highZ_;
	maxZ = maxZ_;
	shifty = (maxZ+1);
    shiftx = (maxY+1)*shifty;
	if (lowX > highX) {
		Message(debug,"error in call to LatticeRange3D::LatticeRange3D"
			"lowX is greater than highX");
	}
	if (lowX < 0) {
		Message(debug,"error in call to LatticeRange3D::LatticeRange3D"
			"lowX is smaller than 1");
	}
	if (highX > maxX) {
		Message(debug,"error in call to LatticeRange3D::LatticeRange3D"
			"highX is greater than maxX");
	}
	if (lowY > highY) {
		Message(debug,"error in call to LatticeRange3D::LatticeRange3D"
			"lowY is greater than highY");
	}
	if (lowY < 0) {
		Message(debug,"error in call to LatticeRange3D::LatticeRange3D"
			"lowY is smaller than 1");
	}
	if (highY > maxY) {
		Message(debug,"error in call to LatticeRange3D::LatticeRange3D"
			"highY is greater than maxY");
	}
	if (lowZ > highZ) {
		Message(debug,"error in call to LatticeRange3D::LatticeRange3D"
			"lowZ is greater than highZ");
	}
	if (lowZ < 0) {
		Message(debug,"error in call to LatticeRange3D::LatticeRange3D"
			"lowZ is smaller than 1");
	}
	if (highZ > maxZ) {
		Message(debug,"error in call to LatticeRange3D::LatticeRange3D"
			"highZ is greater than maxZ");
	}
}
int
LatticeRange3D::Dimensions() const {
	return 3;
}


void LatticeRange3D::MapVector(Vector out, Vector in) const {double *p_out = &out[1], *p_in = &in[1]; MapVector(p_out,p_in);}
void LatticeRange3D::MapVector(double *out, double *in) const {
	int x,y,z,i=0;

	for (x=0; x<=maxX; x++) for (y=0; y<=maxY; y++) for(z=0; z<=maxZ; z++){
		if ((x < lowX || x > highX) || (y<lowY || y>highY) || (z<lowZ || z>highZ)) {out[i]=0;} else {out[i]=in[i];} i++;
	}
}

bool
LatticeRange3D::InRange(int position) const {
	int x;
    	position--;
	x = (position) / shiftx;
	if ( (x<lowX) || (x>highX) || ( (position%shifty)<lowZ ) || ( (position%shifty)>highZ ) || ( (position%shiftx/shifty)<lowY ) || ( (position%shiftx/shifty)>highY ) ) {
		return false;
	}
	return true;
}

double
LatticeRange3D::GetRangeValue(int position) const {
	if (InRange(position))  return 1.0;
	return 0.0;
}

Boolean
LatticeRange3D::InRangeOld(int position) const {
	int x,y,z;
    position--; // position must be one removed since in SFbox all arrays start with unity.
	z = ((position % shiftx) % shifty);
	y = ((position % shiftx)/*-z */)/shifty;
	x = (position-y*shifty /*-z */)/shiftx;
	if (x < lowX || x > highX || y<lowY || y>highY || z<lowZ || z>highZ) {return false;} else {return true;}
}

int
LatticeRange3D::GetNumLayers() const {
	return (highX-lowX+1)*(highY-lowY+1)*(highZ-lowZ+1);
}

int
LatticeRange3D::GetNumLayers_x() const {
	return highX-lowX+1;
}

int
LatticeRange3D::GetNumLayers_y() const {
	return highY-lowY+1;
}

int
LatticeRange3D::GetNumLayers_z() const {
	return highZ-lowZ+1;
}


int
LatticeRange3D::GetNumPos()  {
	int aantal=0;
	//Message(fatal,"GetNumPos not used in 3D without mask file");
	return aantal;
}

double
LatticeRange3D::GetVolPos() {
	return 0.0;
}

int
LatticeRange3D::GetPos(int pos_number)  {
    int z=0;
    Message(fatal,"GetPos not used in 3D without mask file");
	return z;
}

bool
LatticeRange3D::SetPosLocal(int pos, double waarde)  {
    Message(fatal,"SetPosLocal not used in 3D without mask file");
	return false;
}

bool
LatticeRange3D::ClearAllPos()  {
	  Message(fatal,"ClearAllPos not used without mask file");
		return false;
}


bool
LatticeRange3D::SetPos(int pos, int SpotType, int SpotSize, double waarde)  {
    Message(fatal,"SetPos not used in 3D without mask file");
	return false;
}

bool
LatticeRange3D::SetPos(int pos, int SpotSize, double waarde)  {
    Message(fatal,"SetPos not used in 3D without mask file");
	return false;
}

bool
LatticeRange3D::SetPos(double x, double y, double z, double *submask)  {
    Message(fatal,"SetPos not used in 3D without mask file");
	return false;
}

bool
LatticeRange3D::ChangePos(int pos_new, int pos_old){
	Message(fatal,"ChangePos not used in 3D without mask file");
	return false;
}

int
LatticeRange3D::Get_jx() const {
	return (highY-lowY+3)*(highZ-lowZ+3);
}

int
LatticeRange3D::Get_jy() const {
	return highZ-lowZ+3;
}

int
LatticeRange3D::Get_jz() const {
	return 1;
}

int
LatticeRange3D::Get_x0() const {
	return lowX;
}

int
LatticeRange3D::Get_y0() const {
	return lowY;
}

int
LatticeRange3D::Get_z0() const {
	return lowZ;
}

int
LatticeRange3D::Get_xm() const {
	return highX;
}

int
LatticeRange3D::Get_ym() const {
	return highY;
}

int
LatticeRange3D::Get_zm() const {
	return highZ;
}

int
LatticeRange3D::Get_M() const {
	return (highX-lowX+3)*(highY-lowY+3)*(highZ-lowZ+3);
}


Text
LatticeRange3D::GetOutput() const {
	Text t1 = Blanks(1000);
	Text t2 = Blanks(1000);
	Text t3 = Blanks(1000);
	Text t4 = Blanks(1000);
	Text t5 = Blanks(1000);
	Text t6 = Blanks(1000);
	t1.Putint(lowX);
	t1 = Copy(t1.Frontstrip());
	t2.Putint(lowY);
	t2 = Copy(t2.Frontstrip());
	t3.Putint(lowZ);
	t3 = Copy(t3.Frontstrip());

	t4.Putint(highX);
	t4 = Copy(t4.Frontstrip());
	t5.Putint(highY);
	t5 = Copy(t5.Frontstrip());
	t6.Putint(highZ);
	t6 = Copy(t6.Frontstrip());

	Text t = t1 + "," + t2 + "," + t3 + ";" + t4 + "," + t5 + "," + t6;
	return t;
}
