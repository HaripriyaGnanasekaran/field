#ifndef LAT1DCYL1STOxH
#define LAT1DCYL1STOxH

#include "Lat1DCylinder.h"
#include "Lat1DFlat1stO.h"
#include "Input.h"

///
class Lat1DCyl1stO : virtual public Lat1DCylinder, virtual public Lat1DFlat1stO {
  public:
///
	Lat1DCyl1stO() {};
///
	Lat1DCyl1stO(Input*, Text);
///
	virtual ~Lat1DCyl1stO();
	void GetLatticeInfo(int*) const; 


	/// following functions are affected by an overflow protection
	void PropagateG(Matrix ,const Vector ,const int) const;
	void PropagateG(const Vector, const Vector, Vector) const;
	void PropagateG(Vector ,const Vector) const;
	void PropagateG(Matrix ,const Vector ,const int, const double) const;
	void PropagateG(const Vector, const Vector, Vector, const double) const;
	void PropagateG(Vector ,const Vector, const double) const;
///
	void Propagate2G(Matrix ,Matrix ,const Vector, const int, const LatticeRange*) const;
///
	void Propagate2G(Vector ,Vector ,const Vector ,const LatticeRange*) const;
///
	void Propagate2G(Matrix ,Matrix ,const Vector, const int, const LatticeRange*, const double) const;
///
	void Propagate2G(Vector ,Vector ,const Vector ,const LatticeRange*, const double) const;
///
	double ComputeLnGN(const Vector) const;
///
	void NormPhiRestr(Vector, const Vector, double) const;
};

#endif
