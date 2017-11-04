#ifndef SF_SOLVE3RDGENxH
#define SF_SOLVE3RDGENxH

#include "SF_Solve.h"

///
class SF_Solve3rdGen : public SF_Solve {
  public:
///
	SF_Solve3rdGen(Boolean, 
				   SF_ReactionList*,
				   SF_MoleculeList*,
				   SF_SegmentList*,
				   Lattice*,
				   SF_ReactionList*,
				   SF_MoleculeList*,
				   SF_SegmentList*,
				   Lattice*,
				   Input*);
///
	virtual ~SF_Solve3rdGen();
	
///
	void Iterate(void);
///
	void GetResiduals(Vector) const;
  protected:
///
	SF_ReactionList* ReactionQRN;
///
	SF_MoleculeList* MolQRN;
///
	SF_SegmentList* SegQRN;
///
	Lattice* LatRN;
///
	Vector psi0RN;
///
	Array<Vector> SWFRN;
///
	int MRN;
///
	int MolNum3rdGen;
///
	Array<SF_State*> StateQRN;
///
	void residuals(double *const, double *const);
///
	void inneriteration(double *const, double *const, double);
///
	void ComputePhi(void);
///
	void UpdateSWF(void);
///
	void UpdateItVar(Boolean);
///
	Boolean SegItVar(const Lattice*, const int) const;
///
	Boolean PotItVar(const Lattice*, const int) const;
};

#endif
