#ifndef LAT2DCYL1STOxH
#define LAT2DCYL1STOxH

#include "Lat2DCylinder.h"
#include "Lat2DFlat1stO.h"
#include "Input.h"

///
class Lat2DCyl1stO : virtual public Lat2DCylinder, virtual public Lat2DFlat1stO {
public:
///
	Lat2DCyl1stO(Input*, Text);
///
	~Lat2DCyl1stO();
	
// following functions are affected by an overflow protection
///
	void GetLatticeInfo(int*) const; 


	void PropagateG(const Vector, const Vector, Vector) const;
	void PropagateG(Matrix, const Vector, const int) const;
	void PropagateG(Vector, const Vector) const;
	void PropagateG(const Vector, const Vector, Vector, const double) const;
	void PropagateG(Matrix, const Vector, const int, const double) const;
	void PropagateG(Vector, const Vector, const double) const;

	void Propagate2G(Matrix, Matrix, const Vector, const int, const LatticeRange*) const;
/** Overloads Lat1stO::\Ref{Lat1stO::Propagate2G}.*/
	void Propagate2G(Vector, Vector, const Vector, const LatticeRange*) const;
/** Overloads Lat1stO::\Ref{Lat1stO::Propagate2G}.*/
	void Propagate2G(Matrix, Matrix, const Vector, const int, const LatticeRange*, const double) const;
/** Overloads Lat1stO::\Ref{Lat1stO::Propagate2G}.*/
	void Propagate2G(Vector, Vector, const Vector, const LatticeRange*, const double) const;	
	
	
	
	/** Overloads Lat1stO::\Ref{Lat1stO::ComputeLnGN}.*/
	double ComputeLnGN(const Vector) const;
protected:
///
	Lat2DCyl1stO() {};
private:
///
	Lat2DCyl1stO(const Lat2DCyl1stO &);
///
	Lat2DCyl1stO &operator=(const Lat2DCyl1stO &);
};

#endif
