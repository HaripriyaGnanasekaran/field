#ifndef LAT2DCYL1STOSTENSxH
#define LAT2DCYL1STOSTENSxH

#include "Lat2DCylinderSten.h"
#include "Lat2DCyl1stOS.h"
#include "Input.h"

class Lat2DCyl1stOStenS : virtual public Lat2DCylinderSten, virtual public Lat2DCyl1stOS {
public:
///
	Lat2DCyl1stOStenS(Input*, Text);
///
	~Lat2DCyl1stOStenS(); 
	void GetLatticeInfo(int*) const; 

///
	void PropagateG(const Vector, const Vector, Vector) const;
	void PropagateG(Matrix, const Vector, const int) const;
	void PropagateG(Vector, const Vector) const;
	void PropagateG(Matrix, const Vector, const int,const double) const;
	void PropagateG(Vector, const Vector, const double) const;
	void PropagateG(const Vector, const Vector, Vector,const double) const;

///
	void Propagate2G(Matrix, Matrix, const Vector, const int, const LatticeRange*) const;
///
	void Propagate2G(Vector, Vector, const Vector, const LatticeRange*) const;	
///
	void Propagate2G(Matrix, Matrix, const Vector, const int, const LatticeRange*, const double) const;
///
	void Propagate2G(Vector, Vector, const Vector, const LatticeRange*, const double) const;
private:
	double logl1;
	Vector logl0;
	Vector logl11;
	Vector logl1_1;
	Vector logl21;
	Vector logl2_1;
};

#endif
