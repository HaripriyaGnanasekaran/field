#ifndef SF_SOLVECOPELxH
#define SF_SOLVECOPELxH

#include "SF_Solve.h"

///
class SF_SolveCopel : public SF_Solve {
  public:
///
	SF_SolveCopel(Boolean,
				  SF_ReactionList*,
				  SF_MoleculeList*,
				  SF_SegmentList*,
				  Lattice*,
				  Input*);
///
	virtual ~SF_SolveCopel();
	
///
	void GetOutput(Output*) const;
///
	void GetResiduals(Vector) const;
  protected:
///
	void UpdateSWF(void);
///
	void UpdateItVar(Boolean);
///
	void residuals(double *const, double *const);
};

#endif
