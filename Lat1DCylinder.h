#ifndef LAT1DCYLINDERxH
#define LAT1DCYLINDERxH

#include "Lat1DFlat.h"

///
class Lat1DCylinder : virtual public Lat1DFlat {
  public:
///
	Lat1DCylinder();
///
	Lat1DCylinder(Input*,Text);
///
	virtual ~Lat1DCylinder();

///
	void GetOutput(Output*) const;
///
	double GetNumLatticeSites(const int) const; ///takes int numLayer
	double GetLambda(const int, const int) const;
///
	void MultiplyWithLatticeSites(Vector) const;
///
	void DivideByLatticeSites(Vector) const;
///
	double GetVolume(void) const;
///
	double GetTotalNumLatticeSites(void) const;
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
	void SetLayerAdjustment(double);
  protected:
///
	Vector lambda1;
///
	Vector lambda_1;
///
	Vector lambda1_chi;
///
	Vector lambda_1_chi;
///
	Vector L;
///
	double offsetLayer1;
///
	void calculateLambdas();
};

#endif
