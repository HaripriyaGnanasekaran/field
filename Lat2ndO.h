#ifndef LAT2NDOxH
#define LAT2NDOxH

#include "Lattice.h"

///
class Lat2ndO : public virtual Lattice {
  public:
///
	virtual ~Lat2ndO() {};
///
	MatrixApproach GetMatrixApproach() const {return secondOrder;}
///
	virtual void MakeSafe(Vector) const = 0;

///
	virtual void RestoreFromSafe(Vector) const = 0;
	virtual void GetLatticeInfo(int* ) const = 0;
///
	virtual void PropagateF(Matrix, const Vector, const int, const double) const = 0; //intermediate and store
	virtual void PropagateF(Vector, const Vector, const double) const = 0; //intermediate
	virtual void PropagateF(Matrix, const Vector, const int, const double, const bool) const = 0; //for inhomogeneous stiffness;intermediate and store
	virtual void PropagateF(Vector, const Vector, const double, const bool) const = 0; //for inhomogeneous stiffness;intermediate
	//virtual void PropagateG(const Vector Gi, const Vector G, Vector Gout) const = 0;

	virtual void PropagateF(Matrix, const Vector, const int) const = 0; //end and store
	virtual void PropagateF(Vector, const Vector) const = 0; //end;

	virtual void PropagateB(Vector, const Vector, const double) const = 0; //intermediate segments
	virtual void PropagateB(Vector, const Vector, const double, const bool) const = 0; //for inhomogeneous stiffness; intermediate segments
	virtual void PropagateB(Vector, const Vector) const = 0; //start
///
	virtual void Init2G(Vector, Vector, const Vector, const LatticeRange*) const = 0;
///

	virtual void PropagateG(const Vector Gi, const Vector G, Vector Gout) const = 0;
	virtual void PropagateG(const Vector Gi, const Vector G, Vector Gout,const double) const = 0;

	virtual void Propagate2G(Matrix, Matrix, const Vector, const int, const LatticeRange*) const = 0;
	virtual void Propagate2G(Matrix, Matrix, const Vector, const int, const LatticeRange*,const double) const = 0;
///
	virtual void Propagate2G(Vector, Vector, const Vector, const LatticeRange*) const = 0;
	virtual void Propagate2G(Vector, Vector, const Vector, const LatticeRange*,const double) const = 0;
///
	virtual void ConnectG(const Vector, const Matrix, const int, Vector) const = 0;
	virtual void ConnectG(const Vector, const Matrix, const int, Vector, double*) const = 0;
///
	virtual void ConnectG(const Vector, const Vector, Vector) const = 0;
	virtual void ConnectG(const Vector, const Vector, Vector, double*) const = 0;
///
	virtual void Connect2G(const Vector, const Matrix, const int, const Vector, const Matrix, const int, Vector) const = 0;
///
	virtual void Connect2G(const Vector, const Vector, const Vector, const Vector, Vector) const = 0;
///
	virtual void CorrectDoubleCountG(Vector, const Vector) const = 0;
	virtual void CorrectDoubleCountG(Vector, double*, const Vector) const = 0;
/** The natural logarithm of the sum of Gi over the lattice is calculated.
	\f[ \sum_{z=1}^M \left( Gi(z) L(z) \right) \f]
	@param Gi A Vector with indexes ranging from 1 to the amount of layers in the lattice, L.*/
	virtual double ComputeLnGN(const Vector Gi) const = 0;
///
	virtual void NormPhiFree(Vector, const double) const = 0;
	virtual void NormPhiFree(Vector, double*, const double) const = 0;
///
	virtual void NormPhiRestr(Vector, const Vector, double) const = 0;
	virtual void NormPhiRestr(Vector, double*, const Vector, double) const = 0;
};

#endif
