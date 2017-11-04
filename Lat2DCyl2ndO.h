#ifndef LAT2DCYL2NDOxH
#define LAT2DCYL2NDOxH

#include "Lat2DCylinder.h"
#include "Lat2ndO.h"
#include "Input.h"

///
class Lat2DCyl2ndO : virtual public Lat2DCylinder, virtual public Lat2ndO {
public:
///
	Lat2DCyl2ndO() {};
///
	Lat2DCyl2ndO(Input*, Text);
///
	~Lat2DCyl2ndO();

	/// following functions are affected by an overflow protection
	MatrixApproach GetMatrixApproach() const {return secondOrder;}
	void GetLatticeInfo(int*) const;
	Boolean OverflowProtection() const;
	void MakeSafe(Vector) const;
	void RestoreFromSafe(Vector) const;
	void LReflect(double *Pout, double *Pin,int pos) const;
	void UReflect(double *Pout, double *Pin,int pos) const;
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
	void Init2G(Vector, Vector, const Vector, const LatticeRange*) const;
///
	void Propagate2G(Matrix, Matrix, const Vector, const int, const LatticeRange*) const;
///
	void Propagate2G(Vector, Vector, const Vector, const LatticeRange*) const;
///
	void PropagateG(Matrix, const Vector, const int) const;
///
	void PropagateG(Vector, const Vector) const;
///
	void PropagateG(const Vector Gi, const Vector G, Vector Gout) const;

	void PropagateG(Matrix, const Vector, const int,const double) const;
///
	void PropagateG(Vector, const Vector,const double) const;
///
	void PropagateG(const Vector Gi, const Vector G, Vector Gout, const double) const;
///
	void Propagate2G(Matrix, Matrix, const Vector, const int, const LatticeRange*, const double) const;
///
	void Propagate2G(Vector, Vector, const Vector, const LatticeRange*, const double) const;


	void ConnectG(const Vector, const Matrix, const int, Vector) const;
	void ConnectG(const Vector, const Matrix, const int, Vector, double*) const;
	//void ConnectG(const Vector, const Matrix, const int, Vector, Vector) const;
///
	void ConnectG(const Vector, const Vector, Vector) const;
	void ConnectG(const Vector, const Vector, Vector, double*) const;
	//void ConnectG(const Vector, const Vector, Vector, Vector) const;
///
	void Connect2G(const Vector, const Matrix, const int, const Vector, const Matrix, const int, Vector) const;
///
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
	double *H,*l1,*l_1,*l11,*l_11;

protected:
	static Vector Gx0,Hv;

///
	void UpdateBoundaries(Vector) const ;
///
	void UpdateBoundaries(Matrix,const int) const;

};

#endif
