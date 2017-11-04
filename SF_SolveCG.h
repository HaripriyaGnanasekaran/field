#ifndef SF_SolveCGxH
#define SF_SolveCGxH

#include "SF_Solve.h"

///
class SF_SolveCG : public SF_Solve {
  public:
///
	SF_SolveCG(Boolean,  SF_ReactionList*, SF_MoleculeList*, SF_SegmentList*,  Lattice*,  Input*);
///
	virtual ~SF_SolveCG();
	
///
	void GetOutput(Output*) const;	
	void GetResiduals(Vector) const;

protected:
///	
	void UpdateSWF(void);
	Vector u_prime;
	void UpdateItVar(Boolean);
	void residuals(double *const, double *const);
	void inneriteration(double *const, double *const, double);
	void ProcessInput(void);
 

};

#endif
