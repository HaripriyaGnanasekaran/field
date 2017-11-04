#ifndef SF_REACTIONxH
#define SF_REACTIONxH

#include "SF_SegmentList.h"
#include "SF_MoleculeList.h"
#include "SF_StateRatio.h"

///
class SF_Reaction {
  public:
///
	SF_Reaction(Text,SF_MoleculeList*,SF_SegmentList*,Input*);
///
	~SF_Reaction(void);
///
	Text GetName(void);
///
	void GetOutput(Output*);
///
	double ComputeResidue();
///
	int GetNumDiffStates();
	
///
	SF_State* GetDiffState(int);
///
	void AssignStateRatios(void);
///
	int GetNumStateRatios(void);
///
	SF_StateRatio* GetStateRatio(int);
///
	void SwapDirection(void);
///
	int GetNumDiffSegments(void);
///
	void ReplaceRatio(SF_StateRatio*,SF_StateRatio*);
///
	Boolean IterateAlphas(void);
///
	void ComputeInternalFreeEnergies();
///
	void UpdateReactionConstant(double,double);
  private:
///
	Text name;
///
	SF_MoleculeList* MolQ;
///
	SF_SegmentList* SegQ;
///
	Input* MyInput;
///
	Text equation;
///
	double pK;
///
	double pKeff;
///
	int numLeft;
///
	int numRight;
///
	Array<SF_State*> LeftState;
///
	Array<SF_State*> RightState;
///
	Array<Boolean> LeftComplex;
///
	Array<Boolean> RightComplex;
///
	Array<int> LeftStoch;
///
	Array<int> RightStoch;
///
	int numRatios;
///
	Array<SF_StateRatio*> Ratios;
/// 
	Boolean internalFreeEnergiesGiven;
///
	Boolean CreatedOrDestroyed(SF_State*);
///
	void CalculatePKeff(void);
};


#endif
