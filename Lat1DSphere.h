#ifndef LAT1DSPHERExH
#define LAT1DSPHERExH

#include "Lat1DCylinder.h"

///
class Lat1DSphere : virtual public Lat1DCylinder {
  protected:
///
	void calculateLambdas();
  public:
///
	Lat1DSphere(Input*,Text);
///
	virtual ~Lat1DSphere();

///
	void GetOutput(Output*) const;
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
	void SetLayerAdjustment(double);
};

#endif
