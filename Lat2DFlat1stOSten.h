#ifndef Lat2DFlat1stOStenxH
#define Lat2DFlat1stOStenxH

#include "Lat2DFlat1stO.h"

/** This class enhances Lat2DFlat1stO by adding the ability to use stencil
	coefficients.*/
class Lat2DFlat1stOSten : public Lat2DFlat1stO {
  public:
/** The public constructor.
	@param input The inputfile parsed for the current calculation by the Input
		class.
	@param name_ The name given to the lattice in the input file. */
	Lat2DFlat1stOSten(Input* input, Text name_);
/// The destructor.
	~Lat2DFlat1stOSten();
	void GetLatticeInfo(int*) const; 

// following functions are affected by an overflow protection
///
	void PropagateG(const Vector, const Vector, Vector) const;
	void PropagateG(Matrix, const Vector, const int) const;
	void PropagateG(Vector, const Vector) const;
	void PropagateG(const Vector, const Vector, Vector,const double) const;
	void PropagateG(Matrix, const Vector, const int, const double) const;
	void PropagateG(Vector, const Vector, const double) const;


/** Overloads Lat1stO::\Ref{Lat1stO::Propagate2G}.*/
	void Propagate2G(Matrix, Matrix, const Vector, const int, const LatticeRange*) const;
/** Overloads Lat1stO::\Ref{Lat1stO::Propagate2G}.*/
	void Propagate2G(Vector, Vector, const Vector, const LatticeRange*) const;	
/** Overloads Lat1stO::\Ref{Lat1stO::Propagate2G}.*/
	void Propagate2G(Matrix, Matrix, const Vector, const int, const LatticeRange*,const double) const;
/** Overloads Lat1stO::\Ref{Lat1stO::Propagate2G}.*/
	void Propagate2G(Vector, Vector, const Vector, const LatticeRange*,const double) const;
	
protected:
/// The default constructor can only be called by another Lattice class.
	Lat2DFlat1stOSten() {};
/** Overloads Lat2DFlat1stO::\Ref{Lat2DFlat1stO::UpdateBoundaries}.*/
	void UpdateBoundaries(Vector) const;
/** Overloads Lat2DFlat1stO::\Ref{Lat2DFlat1stO::UpdateBoundaries}.*/
	void UpdateBoundaries(Matrix,const int) const;
private:
/** To avoid unwanted lattice creation the copy constructor has been made
	private. */
	Lat2DFlat1stOSten(const Lat2DFlat &);
/**	To avoid unwanted lattice corruption the assignment operator has been
	made private. */
	Lat2DFlat1stOSten &operator=(const Lat2DFlat &);
};

#endif
