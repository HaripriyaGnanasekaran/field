#ifndef LAT2DFLATxH
#define LAT2DFLATxH

#include <assert.h>
#include "Lattice.h"
#include "LatticeRange2D.h"
#include "LatticeRangeFile.h"
#include "Input.h"

///
class Lat2DFlat : public virtual Lattice {
public:
/** The public constructor.
	@param input The inputfile parsed for the current calculation by the Input
		class.
	@param name_ The name given to the lattice in the input file. */
	Lat2DFlat(Input* input, Text name_);
/// The destructor.
	~Lat2DFlat();

/** Overloads Lattice::\Ref{Lattice::GetName}.*/
	Text GetName() const;

/** Overloads Lattice::\Ref{Lattice::GetOutput}.*/
	void GetOutput(Output*) const;
/** Overloads Lattice::\Ref{Lattice::GetTotalNumLayers}.*/
	int GetTotalNumLayers(void) const;
	int Get_N_comp_ranges(void) const;
	int Get_BL(const int) const;
	int GetR0(const int) const;
	int GetRm(const int) const;
/** Overloads Lattice::\Ref{Lattice::GetNumGradients}.*/
	int GetNumGradients(void) const;
/** Overloads Lattice::\Ref{Lattice::GetNumLayers}.*/
	int GetNumLayers(const int) const; //takes int numGradient
/** Overloads Lattice::\Ref{Lattice::GetNumLatticeSites}.*/
	double GetNumLatticeSites(const int) const;//takes int numLayer
/** Overloads Lattice::\Ref{Lattice::GetLambda}.*/
	double GetLambda(const int, const int) const;
/** Overloads Lattice::\Ref{Lattice::MultiplyWithLatticeSites}.*/
	void MultiplyWithLatticeSites(Vector) const;
/** Overloads Lattice::\Ref{Lattice::DivideByLatticeSites}.*/
	void DivideByLatticeSites(Vector) const;
/** Overloads Lattice::\Ref{Lattice::GetVolume}.*/
	double GetVolume(void) const;
/** Overloads Lattice::\Ref{Lattice::GetTotalNumLatticeSites}.*/
	double GetTotalNumLatticeSites() const;
/** Overloads Lattice::\Ref{Lattice::GetNumExtraLatticeSites}.*/
	double GetNumExtraLatticeSites() const;
/** Overloads Lattice::\Ref{Lattice::GetSiteDistance}.*/
	double GetSiteDistance(void) const;
/** Overloads Lattice::\Ref{Lattice::GetSiteSurface}.*/
	double GetSiteSurface(void) const;
/** Overloads Lattice::\Ref{Lattice::NewLatticeRange}.*/
	LatticeRange* NewLatticeRange(Text) const;
/** Overloads Lattice::\Ref{Lattice::WithinBoundaries}.*/
	Boolean WithinBoundaries(const int) const;
/** Overloads Lattice::\Ref{Lattice::SetBoundaries}.*/
	void SetBoundaries(Vector) const;
/** Overloads Lattice::\Ref{Lattice::SetBulkBoundaries}.*/
	void SetBulkBoundaries(Vector, const Vector) const;
/** Overloads Lattice::\Ref{Lattice::SubtractBoundaries}.*/
	void SubtractBoundaries(Vector) const;
/** Overloads Lattice::\Ref{Lattice::RestoreBoundaries}.*/
	void RestoreBoundaries(Vector) const;
/** Overloads Lattice::\Ref{Lattice::GetNamesBulkBoundaries}.*/
	Array<Text> GetNamesBulkBoundaries(void) const;
/** Overloads Lattice::\Ref{Lattice::BulkBoundary}.*/
	Boolean BulkBoundary(const int) const;
/** Overloads Lattice::\Ref{Lattice::NumBulkBoundary}.*/
	int NumBulkBoundary(const int) const;
/** Overloads Lattice::\Ref{Lattice::CheckBoundaries}.*/
	void CheckBoundaries(Array<LatticeRange*>) const;
/** Overloads Lattice::\Ref{Lattice::SideFraction}.*/
	double SideFraction(const Vector, int) const;
	void  Sides(Vector,Vector) const;
/** Overloads Lattice::\Ref{Lattice::ElectricPotential}.*/
	void ElectricPotential(Vector &, Vector &, const Vector &, const Vector &, const double) const;
/** Overloads Lattice::\Ref{Lattice::ElectricFieldSquared}.*/
	void ElectricFieldSquared(Vector &, const Vector, const double) const;
/** Overloads Lattice::\Ref{Lattice::Div1Grad2}.*/
	Vector Div1Grad2(const Vector, const Vector) const;
/** Overloads Lattice::\Ref{Lattice::FourierTransform}.*/
	Vector FourierTransform(const Vector) const;
/** Overloads Lattice::\Ref{Lattice::Moment}.*/
	double Moment(const Vector, const int, const LatticeRange*, const double = 0) const;
/** Overloads Lattice::\Ref{Lattice::MomentUnweighted}.*/
	double MomentUnweighted(const Vector, const int, const LatticeRange*, const double = 0) const;
/** Overloads Lattice::\Ref{Lattice::RenormPhi}.*/
	Vector RenormPhi(Vector, Vector, const Vector, const  double) const;
/** Overloads Lattice::\Ref{Lattice::GetLayerAdjustment} and is empty.*/
	double GetLayerAdjustment(void) const;
/** Overloads Lattice::\Ref{Lattice::SetLayerAdjustment} and is empty.*/
	void SetLayerAdjustment(double);

protected:
/** This constant Text contains the string ``lat''. */
	static const Text lat;
/** This constant Text contains the string ``dimensions''. */
	static const Text dimensionsS;
/** This constant Text contains the string ``gradients''. */
	static const Text gradients;
/** This constant Text contains the string ``geometry''. */
	static const Text geometry;
/** This constant Text contains the string ``n\_layers\_x''. */
	static const Text n_layers_x;
/** This constant Text contains the string ``n\_layers\_y''. */
	static const Text n_layers_y;
/** This constant Text contains the string ``lowerbound\_x''. */
	static const Text lowerbound_x;
/** This constant Text contains the string ``upperbound\_x''. */
	static const Text upperbound_x;
/** This constant Text contains the string ``lowerbound\_y''. */
	static const Text lowerbound_y;
/** This constant Text contains the string ``upperbound\_y''. */
	static const Text upperbound_y;
/** This constant Text contains the string ``distance''. */
	static const Text distance;
/** This constant Text contains the string ``bondlength''. */
	static const Text bondlength;
/** This constant Text contains the string ``lambda''. */
	static const Text lambda;
/** This constant Text contains the string ``lambda0''. */
	static const Text lambda0;
/** This constant Text contains the string ``lambda1''. */
	static const Text lambda1;
/** This constant Text contains the string ``lambda2''. */
	static const Text lambda2;
/** This constant Text contains the string ``latticetype''. */
	static const Text latticetype;

/** This constant Text contains the string ``standard''. */
	static const Text standard;
/** This constant Text contains the string ``stencils''. */
	static const Text stencils;
	static const Text hexagonal;

/// The default value for \f$\lambda_1 = 1/6\f$.
	static const double l1stand;
/** The default value for \f$\lambda_0 = 0,434637428282995\f$ when using
	stencil coefficients.*/
	static const double l0sten;
/** The default value for \f$\lambda_1 = 0,116014619191836\f$ when using
	stencil coefficients.*/
	static const double l1sten;
/** The default value for \f$\lambda_2 = 0,0253260237374153\f$ when using
	stencil coefficients.*/
	static const double l2sten;

/// This Text contains the name of this Lattice as defined in the input file.
	Text name;
/** The number of layers in the x-direction without the borders. In effect
#numLayersX = MX - 2#. */
	int numLayersX;
/** The number of layers in the x-direction without the borders. In effect
#numLayersY = MY - 2#. */
	int numLayersY;
/// The total number of layers in the x-direction (with borders).
	int MX;
/// The total number of layers in the y-direction (with borders).
	int MY;
/// The total number of layers in the lattice, including borders.
	int M;
/** The setting for the lower bound of the system in the x-direction.
	This integer can have different values. The value in the integer refers to
	the x coordinate at which the values for the border can be found:
	- \b 1	The lower border values are defined as \e surface or as \e bulk,
		depending on the value of Boolean #bulkBoundX1. If
		::bulkBoundX1 is false, frozen segments should be defined at the lower
		boundary. If ::bulkBoundX1 is set to true, then a bulk value can be
		defined at the border (used for calculations on dynamics).
	- \b 2	The \e mirror1 boundary is defined. The lower bound has
		the same values as the layers at intern x coordinate 2.
	- \b 3	The \e mirror2 boundary is defined. The lower bound has the
		same values as the layers at x coordinate 3 and the layers at
		x-coordinate 2 are counted as being half in the system.
	- \b MX-1	The \e periodic boundary is defined. The lower bound
		has the same values as the layers at x-coordinate MX-1. */
	int boundX1;
/** The setting for the upper bound of the system in the x-direction.
	The possible values for this parameter are 2, MX-2, MX-1 and MX, which
	have similar meanings to #boundX1.*/
	int boundX2;
/** The setting for the lower bound of the system in the y-direction.
	The possible values for this parameter are 1, 2, 3 and MY-1, which
	have similar meanings to #boundX1.*/
	int boundY1;
/** The setting for the upper bound of the system in the y-direction.
	The possible values for this parameter are 2, MY-2, MY-1 and MY, which
	have similar meanings to #boundX1.*/
	int boundY2;
/** Is true if when a bulk is defined at the lower x-boundary in calculations
	on dynamics. */
	Boolean bulkBoundX1;
/** Is true if when a bulk is defined at the upper x-boundary in calculations
	on dynamics. */
	Boolean bulkBoundX2;
/** Is true if when a bulk is defined at the lower y-boundary in calculations
	on dynamics. */
	Boolean bulkBoundY1;
/** Is true if when a bulk is defined at the upper y-boundary in calculations
	on dynamics. */
	Boolean bulkBoundY2;
/// A Text array containing descriptions of the set boundaries.
	Array<Text> BulkBoundaries;
/// A Text describing the type of Lattice (standard or stencils).
	Text latType;
/// The value for \f$\lambda_0\f$ set by the #lambda0# parameter.
	double l0;
/// The value for \f$\lambda_1\f$ set by the #lambda1# parameter.
	double l1;
	double l1x;
	double l1y;
/** The value for \f$\lambda_2\f$ set by the #lambda2# parameter. This value is
	zero if latType is #standerd#.*/
	double l2;
/// The number of dimensions in the calculation, default: 3.
	double dimensions;
/// The distance between adjecent layers.
	double siteDistance;
///
	double bondLength;
/// A pointer to the processed input file.
	Input* MyInput;
	void SetBy1(double *P)const;

	void SetBy1(double *P,double *Q)const;

	void SetBym(double *P) const;

	void SetBym(double *P,double *Q) const;

	void SetBx1(double *P) const;

	void SetBx1(double *P,double *Q) const;
	void SetBxm(double *P,double *Q) const;
	void SetBxm(double *P) const;
	void removeboundaries(Vector P) const;
	void removeboundaries(double *P) const;



/// The default constructor can only be called by another Lattice class.
	Lat2DFlat() {};
///
  	virtual void checkParameters() const;
private:
/** To avoid unwanted lattice creation the copy constructor has been made
	private. */
	Lat2DFlat(const Lat2DFlat &);
/**	To avoid unwanted lattice corruption the assignment operator has been
	made private. */
	Lat2DFlat &operator=(const Lat2DFlat &);



};

#endif
