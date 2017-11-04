#ifndef LAT2DCYL1STOSTENxH
#define LAT2DCYL1STOSTENxH

#include "Lat2DCyl1stO.h"
#include "Lat2DCylinderSten.h"

/** This class enhances Lat2DCyl1stO by adding the ability to use stencil
	coefficients.*/
class Lat2DCyl1stOSten : virtual public Lat2DCylinderSten, virtual public Lat2DCyl1stO {
public:
/** The public constructor.
	@param input The inputfile parsed for the current calculation by the Input
		class.
	@param name_ The name given to the lattice in the input file. */
	Lat2DCyl1stOSten(Input* input, Text name_);
/// The destructor.
	~Lat2DCyl1stOSten();
	
///
	void GetLatticeInfo(int*) const; 


	void PropagateG(const Vector, const Vector, Vector) const;
	void PropagateG(Matrix, const Vector, const int) const;
	void PropagateG(Vector, const Vector) const;
	void PropagateG(const Vector, const Vector, Vector,const double) const;
	void PropagateG(Matrix, const Vector, const int, const double) const;
	void PropagateG(Vector, const Vector, const double) const;

	void Propagate2G(Matrix, Matrix, const Vector, const int, const LatticeRange*) const;
	void Propagate2G(Vector, Vector, const Vector, const LatticeRange*) const;


/** Overloads Lat1stO::\Ref{Lat1stO::Propagate2G}.*/
	void Propagate2G(Matrix, Matrix, const Vector, const int, const LatticeRange*, const double) const;
/** Overloads Lat1stO::\Ref{Lat1stO::Propagate2G}.*/
	void Propagate2G(Vector, Vector, const Vector, const LatticeRange*, const double) const;	
	
protected:
/// The default constructor can only be called by another Lattice class.
	Lat2DCyl1stOSten(void) {};
private:
/** To avoid unwanted lattice creation the copy constructor has been made
	private. */
	Lat2DCyl1stOSten(const Lat2DCyl1stOSten &);
/** To avoid unwanted lattice corruption the assignment operator has been
	made private. */
	Lat2DCyl1stOSten &operator=(const Lat2DCyl1stOSten &);
};

#endif
