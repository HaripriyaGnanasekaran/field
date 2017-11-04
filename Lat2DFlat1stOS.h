#ifndef LAT2DFLAT1STOSxH
#define LAT2DFLAT1STOSxH

#include "Lat2DFlat.h"
#include "Lat1stO.h"
#include "Input.h"
#include "misc.h"

///
class Lat2DFlat1stOS : virtual public Lat2DFlat, virtual public Lat1stO {
  public:
///
	Lat2DFlat1stOS() {};
///
	Lat2DFlat1stOS(Input*, Text);
///
	~Lat2DFlat1stOS();
	
	/// following functions are affected by an overflow protection
	Boolean OverflowProtection() const;
///
	void MakeSafe(Vector) const;
	void GetLatticeInfo(int*) const; 

///	
	void Init2G(Vector, Vector, const Vector, const LatticeRange*) const;
///
	void RestoreFromSafe(Vector) const;
///
	void PropagateG(const Vector, const Vector, Vector) const;
	void PropagateG(Matrix, const Vector, const int) const;
	void PropagateG(Vector, const Vector) const;
	void PropagateG(const Vector, const Vector, Vector, const double) const;
	void PropagateG(Matrix, const Vector, const int, const double) const;
	void PropagateG(Vector, const Vector, const double) const;

///
	void Propagate2G(Matrix, Matrix, const Vector, const int, const LatticeRange*) const;
///
	void Propagate2G(Vector, Vector, const Vector, const LatticeRange*) const;
///
	void Propagate2G(Matrix, Matrix, const Vector, const int, const LatticeRange*, const double) const;
///
	void Propagate2G(Vector, Vector, const Vector, const LatticeRange*, const double) const;
	
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
	static Vector Gi1Prime;
	static Vector Gi2Prime;
///
	void UpdateBoundaries(Vector) const;
///
	void UpdateBoundaries(Matrix,const int) const;
///
	double getNewGi(const double G, const double b,
		const double *n, short nn) const;
///
	void addLntoLn(const double g, const double a,
		double &value, const bool add) const;
///
	double Log_L1L0, Log_L0, Log_L1;
};

#endif
