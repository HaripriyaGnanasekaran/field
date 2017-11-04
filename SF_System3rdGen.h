#ifndef SF_SYSTEM3RDGENxH
#define SF_SYSTEM3RDGENxH

#include "SF_System.h"
#include "SF_Solve3rdGen.h"

///
class SF_System3rdGen : public SF_System {
  public:
///
	SF_System3rdGen(Input*, Text, Boolean);
///
	virtual ~SF_System3rdGen(void);
///
	void Go(Output*,Lattice*)  const ;
///
	void GetOutput(Output*,Lattice*) const ;
///
	SF_Solve* GetSolve(void) const;
  protected:
///
  	Lat1stO* LatRN;
///
  	SF_SegmentList* SegQRN;
///
  	SF_MoleculeList* MolQRN;
///
  	SF_ReactionList* ReactionQRN;
///
	SF_Solve3rdGen* Solve3rdGen;
///
	SF_Solve3rdGen* NewSolve(void) const;
};

#endif

	
	
	
	
	
	
	
	
	
