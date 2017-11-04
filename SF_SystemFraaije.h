#ifndef SF_SYSTEMFRAAIJExH
#define SF_SYSTEMFRAAIJExH

#include "SF_System.h"
#include "SF_SolveFraaije.h"

///
class SF_SystemFraaije : public SF_System {
  public:
///
	SF_SystemFraaije() {};
///
	SF_SystemFraaije(Input*, Text, Boolean);
///
	virtual ~SF_SystemFraaije(void);

///
	void Go(Output*,Lattice*);
///
	void SetInitialGuess(SF_System*,Lattice*);
///
	void GetOutput(Output*,Lattice*) const;
///
	SF_Solve* GetSolve(void) const;
  protected:
///
	SF_SolveFraaije* SolveFraaije;	
///
	SF_SegmentList* SegNewQ;
  ///
  	SF_MoleculeList* MolNewQ;
  ///
  	SF_ReactionList* ReactionNewQ;
///
	SF_SolveFraaije* NewSolve(void) const;
///
	void UpdateChi(void);
///
	double timeNow;
///
	double timeStep;
///
	double timeLimit;
///
	double timeStepIncr;
///
	double error;
///
	double maxError;
///
	double minError;
///
	double stopCriterion;
///
	double stop;
///
	int outputTimerLimit;
};

#endif

	
	
	
	
	
	
	
	
	
