#ifndef SF_STATExH
#define SF_STATExH

#include "SF_MolState.h"
#include "SF_Box.h"

///
class SF_State {
  public:
///
	SF_State();
///
	SF_State(Text,double,double,double,double,Boolean,double,double,SegmentFreedom,LatticeRange*,Lattice*);
///
	~SF_State();

///
	Text GetName(void) const;
///
	void GetOutput(Output*) const;
///
	double GetValence(void) const;
///
	double GetEpsilon(void) const;
///
	Vector GetAlpha(void) const;
///
	Boolean AlphaBulkFixed(void) const;
///
	double GetAlphaBulk(void) const;
///
	double GetAlphaBulkFixed(void) const;
///
	void SetAlphaBulk(double);
///
	Boolean PhiBulkFixed(void) const;
///
	double GetPhiBulk(void) const;
///
	double GetPhiBulkFixed(void) const;
///
	void UpdatePhiBulk(void);
///
	double GetTheta(void) const;
///
	Boolean InternalFreeEnergyGiven(void) const;
///
	void UpdateInternalFreeEnergy(double, double);
///
	void SetInternalFreeEnergy(double);
///
	Boolean InternalFreeEnergySet(void) const;
///
	double GetInternalFreeEnergy(void) const;
///
	SegmentFreedom GetFreedom(void) const;
///
	LatticeRange* GetLatRange(void) const;
///
	void SetSWF(Vector);

///
	Vector GetSWF(void) const;
	Vector GetRHO(void) const;
	double GetRHO(int) const;
	void PutRHO(void);
	void PutRHO(double,int);

	Boolean RhoSet(void);
	void ResetRhoSet();
///
	void DelPos(int,int,int);
	void UpdatePos(int,int,int);
	void UpdatePos(double,double,double,double*);
	void UpdateSWF(void);
	void ClearAllPos(void);


///
	void SetExternalPotential(Vector);
///
	Vector GetExternalPotential(void);
///
	Boolean ExternalPotential(void);
///
	Vector GetPhi(void) const;
///
	int GetNumMolStates(void) const;
///
	SF_MolState* GetMolState(int) const;
///
	void AddMolState(SF_MolState*);
  private:
///
	Text name;
///
	SegmentFreedom freedom;
///
	LatticeRange* LatRange;
///
	double valence;
///
	double epsilon;
///
	Boolean BoolAlphaBulkFixed;
///
	double alphaBulk;
///
	double alphaBulkFixed;
///
	Boolean BoolPhiBulkFixed;
///
	double phiBulk;
///
	double phiBulkFixed;
///
	Boolean equilibrium;
	Boolean Rho_up_to_date;
///
	Head MolStateQ;
///
	Vector swf;
	Vector rho;
///
	Vector swfFull;
///
	Boolean extPot;
///
	Vector externalPotential;
///
	double internalFreeEnergyGiven;
///
	double internalDegeneration;
///
	double internalEnergy;
///
	double internalFreeEnergy;
///
	Boolean internalFreeEnergySet;
///
	Lattice* Lat;
};

#endif
