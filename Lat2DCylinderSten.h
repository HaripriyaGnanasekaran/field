#ifndef LAT2DCYLINDERSTENxH
#define LAT2DCYLINDERSTENxH

#include "Lat2DCylinder.h"

/** This class enhances Lat2DCylinderSten by adding the ability to use stencil
	coefficients.*/
class Lat2DCylinderSten : virtual public Lat2DCylinder {
public:
/** The public constructor.
	@param input The inputfile parsed for the current calculation by the Input
		class.
	@param name_ The name given to the lattice in the input file. */
	Lat2DCylinderSten(Input* input, Text name_);
/// The destructor.
	~Lat2DCylinderSten();
/** Overloads Lattice::GetOutput.
	The same output is given as by Lat2DCylinder::GetOutput
	with the value for #lambda2 and the vectors
	Lat2DCylinderSten::lambda21 and Lat2DCylinderSten::lambda2_1.*/
	void GetOutput(Output* Out) const;
/** Overloads Lattice::\Ref{Lattice::GetLambda}.
	The interaction with the next nearest neighbour is incorparated.*/
	double GetLambda(const int, const int) const;
/** Overloads Lattice::\Ref{Lattice::SideFraction}.
	The interaction with the next nearest neighbour is incorparated.*/
	double SideFraction(const Vector, int) const;
	void Sides(Vector,Vector) const;
protected:
/// Contains the value \f$\lambda_2\f$ for interaction with the surrounding layer.
	Vector lambda21;
/// Contains the value \f$\lambda_2\f$ for interaction with the surrounded layer.
	Vector lambda2_1;

/// The default constructor can only be called by another Lattice class.
	Lat2DCylinderSten(void) {};
/** Internal function to calculate \f$\lambda_0\f$ as a function of x. 
	@param x The x-coordinate of the layer for which \f$\lambda_0\f$ will be
	calculated.*/
	double getLambda0(int x) const;
/** Outputs the value for lambdas for the stencil coefficients. */
	void outputLambdas(Output *) const;
private:
/** To avoid unwanted lattice creation the copy constructor has been made
	private.*/
	Lat2DCylinderSten(const Lat2DCylinderSten &);
/** To avoid unwanted lattice corruption the assignment operator has been made
	private.*/
	Lat2DCylinderSten &operator=(const Lat2DCylinderSten &);
};

#endif
