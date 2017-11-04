#ifndef LAT2DCYLINDERxH
#define LAT2DCYLINDERxH

#include "Lat2DFlat.h"

///
class Lat2DCylinder : virtual public Lat2DFlat {
public:
///
	Lat2DCylinder(Input*,Text);
///
	~Lat2DCylinder();

///
	void GetOutput(Output*) const;
/// takes int numLayer
	double GetNumLatticeSites(const int) const;
///
	double GetLambda(const int,const int) const;
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
	Vector Div1Grad2(const Vector,const Vector) const;
///
	Vector FourierTransform(const Vector) const;
///
	double Moment(const Vector, const int, const LatticeRange*, const double = 0) const;
///
	double Moment(const Vector in, const int m, const int gradient) const;
///
	Vector RenormPhi(Vector, Vector, const Vector, const double) const;
protected:
///
	static const Text n_sites_first_layer_x, offset_first_layer_x;
///
	Vector lambda11;
///
	Vector lambda1_1;
///
	Vector L;
///
	double offsetLayer1;
/// the width of the first layer in the x (perp. to axis) direction
	double layerAdjustment;
///
	Lat2DCylinder() {};
///
	void checkParameters() const;
///
	double getLambda0(int) const;
///
	void outputLambdas(Output *) const;
private:
///
	Lat2DCylinder(const Lat2DCylinder &);
///
	Lat2DCylinder &operator=(const Lat2DCylinder &);
///
	void calculateLambdas();
};

void addtoLnGN(const double g, const double a, double &value, const bool add);
#endif
