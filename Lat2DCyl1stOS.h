#ifndef LAT2DCYL1STOSxH
#define LAT2DCYL1STOSxH

#include "Lat2DCylinder.h"
#include "Lat2DFlat1stOS.h"
#include "Input.h"

/** The overflow protected version of the class Lat2D1stO.
 *
 * The general idea behind overflow protection is to use calculate with log(G)
 * instead of with G.  Please read the documentation for overflow protection
 * (try Jan's thesis).
 *
 * Approach
 *
 * During propagation, we need to select a \f$\lambda\f$ that is always
 * positive, non-zero, because we need the log of this \f$\lambda\f$. A natural
 * choice would be \f$\lambda_{0,0}\f$, but this lambda can be zero for the
 * first layer. Also, the inward pointing \f$\lambda_{0,-1} is not eligable,
 * since it's zero for the first layer.  During the calculation we need the
 * ratio between all lambda's and the chosen lambda's near to it. Since
 * \f$\lambda_{1,0}\f$ and \f$\lambda_{-1,0}\f$ are independant of the distance
 * from the axis, it's most convenient to use one of these. We choose
 * \f$\lambda_{1,0}\f$.
 *
 *
 */
class Lat2DCyl1stOS : virtual public Lat2DCylinder, virtual public Lat2DFlat1stOS {
public:
///
	Lat2DCyl1stOS(Input*, Text);
///
	~Lat2DCyl1stOS(); 
///
	void PropagateG(const Vector, const Vector, Vector) const;
	void PropagateG(Matrix, const Vector, const int) const;
	void PropagateG(Vector, const Vector) const;
	void PropagateG(const Vector, const Vector, Vector, const double) const;
	void PropagateG(Matrix, const Vector, const int, const double) const;
	void PropagateG(Vector, const Vector, const double) const;

	void GetLatticeInfo(int*) const; 

///
	void Propagate2G(Matrix, Matrix, const Vector, const int, const LatticeRange*) const;
///
	void Propagate2G(Vector, Vector, const Vector, const LatticeRange*) const;
///
	void Propagate2G(Matrix, Matrix, const Vector, const int, const LatticeRange*, const double) const;
///
	void Propagate2G(Vector, Vector, const Vector, const LatticeRange*, const double) const;
///	
	double ComputeLnGN(const Vector) const;
protected:
///
	Lat2DCyl1stOS() {};
private:
	double logl1;
	Vector logl11;
	Vector logl0;
	Vector logl1_1;
};

#endif
