#ifndef SF_SOLVE_JOHANxH
#define SF_SOLVE_JOHANxH

#include "SF_MoleculeList.h"
#include "SF_ReactionList.h"
#include "LatticeRangeFile.h"


#include "SF_Solve.h"

class SF_Solve_Johan : public SF_Solve {
  public:
///
	SF_Solve_Johan(Boolean,
				  SF_ReactionList*,
				  SF_MoleculeList*,
				  SF_SegmentList*,
				  Lattice*,
				  Input*);

///
	virtual ~SF_Solve_Johan();


///

///
	virtual void GetResiduals(Vector) const;
///

};


#endif
