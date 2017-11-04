#ifndef LAT1DFLATxH
#define LAT1DFLATxH

#include "Lattice.h"
#include "LatticeRange1D.h"
#include "LatticeRangeFile.h"
#include "Input.h"

///
class Lat1DFlat : public virtual Lattice {
  public:
///
	Lat1DFlat();
///
	Lat1DFlat(Input*, Text);
///
	virtual ~Lat1DFlat();

///
	Text GetName() const;
///
	void GetOutput(Output*) const;
///
	int GetTotalNumLayers(void) const;
	int Get_N_comp_ranges(void) const;
	int Get_BL(const int) const;
	int GetR0(const int) const;
	int GetRm(const int) const;
///
	int GetNumGradients(void) const;
///
	int GetNumLayers(const int) const; ///takes int numGradient
	double GetNumLatticeSites(const int) const; ///takes int numLayer
/// takes two coordinates as arguments, returns the lambda.
/// only use this in chain propagators
	double GetLambda(const int, const int) const;
///
	void MultiplyWithLatticeSites(Vector) const;
///
	void DivideByLatticeSites(Vector) const;
///
	double GetVolume(void) const;
///
	double GetTotalNumLatticeSites(void) const;
///	Extra lattice sites are an explicit bulk phase, can be used for
/// interaction curves, not included in totalNumLatticeSites
	double GetNumExtraLatticeSites(void) const;
///
	double GetSiteDistance(void) const;
///
	double GetSiteSurface(void) const;
///
	LatticeRange* NewLatticeRange(Text) const;
///
	Boolean WithinBoundaries(const int) const;
///
	void SetBoundaries(Vector) const;
///
	void SetBulkBoundaries(Vector, const Vector) const;
///
	void SubtractBoundaries(Vector) const;
///
	void RestoreBoundaries(Vector) const;
///
	Array<Text> GetNamesBulkBoundaries(void) const;
///
	Boolean BulkBoundary(const int) const;
///
	int NumBulkBoundary(const int) const;
///
	void CheckBoundaries(Array<LatticeRange*>) const;
///
	double SideFraction(const Vector, int) const;

	void Sides(Vector,Vector) const;
///
	void ElectricPotential(Vector &, Vector &, const Vector &, const Vector &, const double) const;
///
	void ElectricFieldSquared(Vector &, const Vector, const double) const;
///
	Vector Div1Grad2(const Vector, const Vector) const;
///
	Vector FourierTransform(const Vector) const;
///
	double Moment(const Vector, const int, const LatticeRange*, const double = 0) const;
///
	double MomentUnweighted(const Vector, const int, const LatticeRange*, const double = 0) const;
///
	Vector RenormPhi(Vector, Vector, const Vector, const double) const;
///
	double GetLayerAdjustment(void) const;
///
	void SetLayerAdjustment(double);
  protected:
///
	Text name;
///
	int numLayers;
///	Extra lattice sites are an explicit bulk phase, can be used for
/// interaction curves, not included in totalNumLatticeSites
	double numExtraLatticeSites;
///
	int M;
///
	int bound1, bound2;
///
	Boolean bulkBound1, bulkBound2;
///
	Array<Text> BulkBoundaries;
///
	Text latType;
	double l1; /// lambda1 for the chain propagator
	double l1_chi; /// lambda1 for the Flory-Huggins interactions
///
	double dimensions;
///
	double siteDistance;
///
	double bondLength;
	Boolean adjustOuterLayer; /// needed to correct for lattice artifact
	double layerAdjustment; /// radius of first layer in system (z=2)
	Input* MyInput;

	void SetBx1(double *P) const;

	void SetBx1(double *P,double *Q) const;

	void SetBxm(double *P) const;

	void SetBxm(double *P,double *Q) const;
	void removeboundaries(Vector P) const;
	void removeboundaries(double *P) const;

};

#endif
