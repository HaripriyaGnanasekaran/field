#ifndef LATTICExH
#define LATTICExH

#include "Input.h"
#include "Output.h"
#include "LatticeRange.h"
#include "Message.h"
#include "SF_Box.h"
///
enum MatrixApproach {firstOrder, 
					 scaf, 
					 secondOrder, 
					 ris, 
					 scafRis, 
					 secondOrderScaf};

/**
	The class Lattice is an abstract class. It is the base class for all
	lattice classes in the SFBox. A class is abstract if it contains one
	or more pure virtual functions. A virtual function is preceded by
	the keyword 'virtual'. A function is a pure virtual function if it's
	declaration is followed by '= 0'. This means the function is not
	defined for this class and therefor no instance of this class can be
	made. Derived classes can only be instantiated, if they override all
	pure virual functions of the base class.

	Virtual functions use dynamic binding to resolve to the correct
	function. If class A has a virtual function a() and a non-virtual
	function b, and if class B inherits class A and overloads both of A's
	functions, than a call to function a() via an object of class B
	referenced by a A pointer, will call function B::a(), whereas a call to
	the non-virtual function b(), via the same pointer, is a call to A::b().
*/
class Lattice {
public:
///
	virtual ~Lattice() {};
	virtual int Get_N_comp_ranges (void) const = 0;
	virtual int Get_BL(const int) const = 0; 
	virtual int GetR0 (const int) const = 0;
	virtual int GetRm (const int) const = 0;

///
	virtual Text GetName() const = 0;
///
	virtual void GetOutput(Output*) const = 0;
///
	virtual MatrixApproach GetMatrixApproach() const = 0;
///
	virtual Boolean OverflowProtection() const = 0;
///
	virtual int GetTotalNumLayers(void) const = 0;
///
	virtual int GetNumGradients(void) const = 0;
///
	virtual int GetNumLayers(const int) const = 0; //takes int numGradient
///
	virtual double GetNumLatticeSites(const int) const = 0; //takes int numLayer
///
	virtual double GetLambda(const int, const int) const = 0; //takes 2* int numLayer
///
	virtual void MultiplyWithLatticeSites(Vector) const = 0;
///
	virtual void DivideByLatticeSites(Vector) const = 0;
///
	virtual double GetVolume(void) const = 0;

///
	virtual double GetTotalNumLatticeSites(void) const = 0;
///
	virtual double GetNumExtraLatticeSites(void) const = 0;
///
	virtual double GetSiteDistance(void) const = 0;
///
	virtual double GetSiteSurface(void) const = 0;
///
	virtual LatticeRange* NewLatticeRange(Text) const = 0;
///
	virtual Boolean WithinBoundaries(const int) const = 0; //takes int numLayer
///
	virtual void SetBoundaries(Vector) const = 0;
///
	virtual void SetBulkBoundaries(Vector, const Vector) const = 0; //second Vector with bulk values
///
	virtual void SubtractBoundaries(Vector) const = 0;
///
	virtual void RestoreBoundaries(Vector) const = 0;
///
	virtual Array<Text> GetNamesBulkBoundaries(void) const = 0;
///
	virtual Boolean BulkBoundary(const int) const = 0; 
///
	virtual int NumBulkBoundary(const int) const = 0; 
///
	virtual void CheckBoundaries(Array<LatticeRange*>) const = 0;
///
	virtual double SideFraction(const Vector, int) const = 0; //takes int numLayer
	
	virtual void Sides(Vector,Vector) const = 0; 

///
	virtual void ElectricPotential(Vector &, Vector &, const Vector &, const Vector &, const double) const = 0;
///
	virtual void ElectricFieldSquared(Vector &, const Vector, const double) const = 0;
///
	virtual Vector Div1Grad2(const Vector, const Vector) const = 0;
///
	virtual Vector FourierTransform(const Vector) const = 0;
///
	virtual double Moment(const Vector, const int, const LatticeRange*, const double = 0) const = 0; 
///
	virtual double MomentUnweighted(const Vector, const int, const LatticeRange*, const double = 0) const = 0; 
///
	virtual Vector RenormPhi(Vector, Vector, const Vector, const double) const = 0;
///
	virtual double GetLayerAdjustment(void) const = 0;
///
	virtual void SetLayerAdjustment(double) = 0;
};

#endif
