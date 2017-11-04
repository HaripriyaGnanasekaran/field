#ifndef LAT1DFLAT1STOSxH
#define LAT1DFLAT1STOSxH

#include "Lat1DFlat1stO.h"
#include "misc.h"

///
class Lat1DFlat1stOS :  virtual public Lat1DFlat, virtual public Lat1stO {
  public:
///
	Lat1DFlat1stOS() {};
///
	Lat1DFlat1stOS(Input*, Text);
///
	virtual ~Lat1DFlat1stOS();

	/// following functions are affected by an overflow protection
	Boolean OverflowProtection() const;
///
	void MakeSafe(Vector) const;
///
	void RestoreFromSafe(Vector) const;
	void GetLatticeInfo(int*) const; 

///	
	void Init2G(Vector, Vector, const Vector, const LatticeRange*) const;
///
	void PropagateG(const Vector, const Vector, Vector) const;
	void PropagateG(Matrix, const Vector, const int) const;
	void PropagateG(Vector, const Vector) const;
	void PropagateG(const Vector, const Vector, Vector,const double f) const;
	void PropagateG(Matrix, const Vector, const int, const double f) const;
	void PropagateG(Vector, const Vector, const double f) const;
///
	void Propagate2G(Matrix, Matrix, const Vector, const int, const LatticeRange*) const;
///
	void Propagate2G(Vector, Vector, const Vector, const LatticeRange*) const;
///
	void Propagate2G(Matrix, Matrix, const Vector, const int, const LatticeRange*,const double f) const;
///
	void Propagate2G(Vector, Vector, const Vector, const LatticeRange*,const double f) const;	
		
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
///
	double LogL1L0, LogL0, LogL1;
};

#endif
