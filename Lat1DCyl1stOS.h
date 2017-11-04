#ifndef LAT1DCYL1STOSxH
#define LAT1DCYL1STOSxH

#include "Lat1DCyl1stO.h"
#include "Lat1DFlat1stOS.h"

///
class Lat1DCyl1stOS : virtual public Lat1DCylinder, virtual public Lat1DFlat1stOS {
  public:
///
	Lat1DCyl1stOS(Input*, Text);
///
	Lat1DCyl1stOS() {};
///
	virtual ~Lat1DCyl1stOS();
	void GetLatticeInfo(int*) const; 

/// following functions are affected by an overflow protection
	void PropagateG(Matrix, const Vector, const int) const;
	void PropagateG(const Vector, const Vector, Vector) const;
	void PropagateG(Vector, const Vector) const;


///
	void Propagate2G(Matrix, Matrix, const Vector, const int, const LatticeRange*) const;
///
	void Propagate2G(Vector, Vector, const Vector, const LatticeRange*) const;
///
	double ComputeLnGN(const Vector) const;
///
	void NormPhiRestr(Vector, const Vector, double) const;
};

#endif
