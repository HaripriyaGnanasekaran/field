#ifndef LAT1DSPHERE2NDOxH
#define LAT1DSPHERE2NDOxH

#include "Lat1DSphere.h"
#include "Lat1DCyl2ndO.h"
#include "Input.h"
#include "Lat2ndO.h"

///
class Lat1DSphere2ndO : virtual public Lat1DCyl2ndO, public Lat1DSphere, virtual public Lat2ndO  {
  public:
///
	Lat1DSphere2ndO(Input*, Text);
///
	virtual ~Lat1DSphere2ndO();

    MatrixApproach GetMatrixApproach() const {return secondOrder;}
	Boolean OverflowProtection(void) const;
	void MakeSafe(Vector) const;
	void RestoreFromSafe(Vector) const;
	void GetLatticeInfo(int*) const;
	void LReflect(double *Pout, double *Pin,int pos) const;
	void UReflect(double *Pout, double *Pin,int pos) const;

	void PropagateF(Matrix,  Vector, const int, const double) const ;
	void PropagateF(Vector,  Vector, const double) const ;
	void PropagateF(Matrix,  Vector, const int, const double, const bool ) const;
	void PropagateF(Vector,  Vector, const double, const bool) const ;
	void PropagateF(Matrix,  Vector, const int) const;
	void PropagateF(Vector,  Vector) const;

	void PropagateB(Vector,  Vector, const double) const;
	void PropagateB(Vector,  Vector, const double, const bool) const;
	void PropagateB(Vector,  Vector) const;

	void ConnectG(Vector GiA, Matrix GiB, int s, Vector out) const;
	void ConnectG(Vector GiA, Matrix GiB, int s, Vector out, double*) const;
	void ConnectG(Vector GiA, Vector GiB, Vector out) const;
	void ConnectG(Vector GiA, Vector GiB, Vector out, double*) const;


	void PropagateG(Matrix ,const Vector ,const int) const;
	void PropagateG(const Vector, const Vector, Vector) const;
	void PropagateG(Vector ,const Vector) const;
	void PropagateG(Matrix ,const Vector ,const int, const double) const;
	void PropagateG(const Vector, const Vector, Vector, const double) const;
	void PropagateG(Vector ,const Vector, const double) const;

	void Propagate2G(Matrix ,Matrix ,const Vector, const int, const LatticeRange*) const;
	void Propagate2G(Vector ,Vector ,const Vector ,const LatticeRange*) const;
	void Propagate2G(Matrix ,Matrix ,const Vector, const int, const LatticeRange*, const double) const;
	void Propagate2G(Vector ,Vector ,const Vector ,const LatticeRange*, const double) const;

	void Connect2G(const Vector, const Matrix, const int, const Vector, const Matrix, const int, Vector) const;
	void Connect2G(const Vector, const Vector, const Vector, const Vector, Vector) const;

	void Init2G(Vector, Vector, const Vector, const LatticeRange*) const;
///
	double ComputeLnGN(const Vector) const;

	void CorrectDoubleCountG(Vector in, const Vector G) const;
	void CorrectDoubleCountG(Vector in, double*, const Vector G) const;
///
	void NormPhiRestr(Vector, const Vector, double) const;
	void NormPhiRestr(Vector, double*, const Vector, double) const;

	void NormPhiFree(Vector phi, const double C) const;
	void NormPhiFree(Vector phi, double*, const double C) const;
	int M;
	double *H;
	double *gx0;

protected:
	//static Vector Gx0,Hv;

	void UpdateBoundaries(Vector) const;
	void UpdateBoundaries(Matrix,const int) const;
};

#endif
