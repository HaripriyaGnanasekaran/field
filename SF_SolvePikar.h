#ifndef SF_SOLVEPIKARxH
#define SF_SOLVEPIKARxH

#include "SF_Solve.h"

///
class SF_SolvePikar : public SF_Solve {
  public:
///
	SF_SolvePikar(Boolean,  SF_ReactionList*, SF_MoleculeList*, SF_SegmentList*,  Lattice*,  Input*);
///
	virtual ~SF_SolvePikar();
	
///
	void GetOutput(Output*) const;	
	void GetResiduals(Vector) const;

	int iterationlimit_picard,its;
	double tolerance_picard,accuracy_picard,deltamax_picard;


protected:
///	
	double *gp;
	void UpdateSWF(void);
	Vector u_prime;
	double NORM(double*, int ) ;
	void iterate(double* ,int) ;
	void UpdateItVar(Boolean);
	void residuals(double *const, double *const);
	void inneriteration(double *const, double *const, double);
	void ProcessInput(void);
	void itinfo(void);  

};

#endif
