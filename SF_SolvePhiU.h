#ifndef SF_SOLVEPhiUxH
#define SF_SOLVEPhiUxH

#include "SF_Solve.h"

///
class SF_SolvePhiU : public SF_Solve {
  public:
///
	SF_SolvePhiU(Boolean,  SF_ReactionList*, SF_MoleculeList*, SF_SegmentList*,  Lattice*,  Input*);
///
	virtual ~SF_SolvePhiU();
	Vector phiTotal_;
	Vector phi_tmp;
///
	void GetOutput(Output*) ;
///
	void GetResiduals(Vector) ;


  protected:
///
	void UpdateSWF(void);
//	Vector u_prime;

	Array<Vector> PHI;

///
	void UpdateItVar(Boolean);
///
	void residuals(double *const, double *const);

	void inneriteration(double *const, double *const, double);
	//void ProcessInput(void);
};

#endif
