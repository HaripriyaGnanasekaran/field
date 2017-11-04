#ifndef LAT3D1STOxH
#define LAT3D1STOxH

#include "Lat3D.h"
#include "Lat1stO.h"
#include "Input.h"

///
class Lat3D1stO : virtual public Lat3D, virtual public Lat1stO {
public:
///
	Lat3D1stO() {};
///
	Lat3D1stO(Input*, Text);
///
	~Lat3D1stO();
	
	/// following functions are affected by an overflow protection
	Boolean OverflowProtection() const;
///
	void MakeSafe(Vector) const;
	void GetLatticeInfo(int*) const; 

	void Init2G(Vector, Vector, const Vector, const LatticeRange*) const;
///
	void RestoreFromSafe(Vector) const;
///
	void PropagateG(Matrix, const Vector, const int) const;
	void PropagateG(Vector, const Vector) const;
	void PropagateG(const Vector Gi, const Vector G, Vector Gout) const;
	void PropagateG(Matrix, const Vector, const int,const double) const;
	void PropagateG(Vector, const Vector,const double) const;
	void PropagateG(const Vector Gi, const Vector G, Vector Gout, const double) const;

///
	void Propagate2G(Matrix, Matrix, const Vector, const int, const LatticeRange*) const;
///
	void Propagate2G(Vector, Vector, const Vector, const LatticeRange*) const;
///
	void Propagate2G(Matrix, Matrix, const Vector, const int, const LatticeRange*, const double) const;
///
	void Propagate2G(Vector, Vector, const Vector, const LatticeRange*, const double) const;	
///
	void ConnectG(const Vector, const Matrix, const int, Vector) const;
///
	void ConnectG(const Vector, const Vector, Vector) const;
///
	void Connect2G(const Vector, const Matrix, const int, const Vector, const Matrix, const int, Vector) const;
///
	void Connect2G(const Vector, const Vector, const Vector, const Vector, Vector) const;
///
	void CorrectDoubleCountG(Vector, const Vector) const;
///
	double ComputeLnGN(const Vector) const;
///
	void NormPhiFree(Vector, const double) const;
///
	void NormPhiRestr(Vector, const Vector, double) const;
protected:
	static Vector Gx;
	static Vector Gt;
	static Vector Gp;
	static Vector G2; 

///
	void UpdateBoundaries(Vector) const;
///
	void UpdateBoundaries(Matrix,const int) const;
	
};

#endif
