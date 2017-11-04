#ifndef LAT1D2NDOxH
#define LAT1D2NDOxH

#include "Lat1DFlat.h"
#include "Lat2ndO.h"
#include "Input.h"

///
class Lat1D2ndO : virtual public Lat1DFlat, virtual public Lat2ndO {
public:
///
	Lat1D2ndO() {};
///
	Lat1D2ndO(Input*, Text);
///
	~Lat1D2ndO();

	/// following functions are affected by an overflow protection
	Boolean OverflowProtection() const;
///
	void MakeSafe(Vector) const;

///
	void RestoreFromSafe(Vector) const;
	void GetLatticeInfo(int*) const;
///
	void PropagateF(Matrix,  Vector, const int, const double) const;
	void PropagateF(Vector,  Vector, const double) const;
	void PropagateF(Matrix,  Vector, const int, const double, const bool ) const;
	void PropagateF(Vector,  Vector, const double, const bool) const;
	void PropagateF(Matrix,  Vector, const int) const;
	void PropagateF(Vector,  Vector) const;


	void PropagateB(Vector,  Vector, const double) const;
	void PropagateB(Vector,  Vector, const double, const bool) const;
	void PropagateB(Vector,  Vector) const;
///
	void ConnectG(Vector, Matrix, int, Vector) const;
	void ConnectG(Vector, Matrix, int, Vector, double*) const;
	void ConnectG(Vector, Vector, Vector) const;
	void ConnectG(Vector, Vector, Vector, double*) const;
	void Init2G(Vector, Vector, const Vector, const LatticeRange*) const;
///
	void Propagate2G(Matrix, Matrix, const Vector, const int, const LatticeRange*) const;
	void Propagate2G(Vector, Vector, const Vector, const LatticeRange*) const;
	void Propagate2G(Matrix, Matrix, const Vector, const int, const LatticeRange*, const double) const;
	void Propagate2G(Vector, Vector, const Vector, const LatticeRange*, const double) const;


	void PropagateG(Matrix, const Vector, const int) const;
	void PropagateG(Vector, const Vector) const;
	void PropagateG(const Vector Gi, const Vector G, Vector Gout) const;
	void PropagateG(Matrix, const Vector, const int,const double) const;
	void PropagateG(Vector, const Vector,const double) const;
	void PropagateG(const Vector Gi, const Vector G, Vector Gout, const double) const;
///
	void Connect2G(const Vector, const Matrix, const int, const Vector, const Matrix, const int, Vector) const;
	void Connect2G(const Vector, const Vector, const Vector, const Vector, Vector) const;
///
	void CorrectDoubleCountG(Vector, const Vector) const;
	void CorrectDoubleCountG(Vector, double*, const Vector) const;
///
	double ComputeLnGN(const Vector) const;
///
	void NormPhiFree(Vector, const double) const;
	void NormPhiFree(Vector, double*, const double) const;
///
	void NormPhiRestr(Vector, const Vector, double) const;
	void NormPhiRestr(Vector, double*, const Vector, double) const;
protected:
	static Vector Gx0;

///
	void UpdateBoundaries(Vector) const;
///
	void UpdateBoundaries(Matrix,const int) const;

};

#endif
