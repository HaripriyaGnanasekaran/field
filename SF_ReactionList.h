#ifndef SF_REACTIONLISTxH
#define SF_REACTIONLISTxH

#include "sfnewton.h"
#include "SF_Reaction.h"
#include "SF_StateRatioList.h"

///
class SF_ReactionList : public SFNewton {
  public:
///
    SF_ReactionList(SF_MoleculeList*, SF_SegmentList*, Input*);
///
    virtual ~SF_ReactionList(void);
///
    void GetOutput(Output*);
///
    int GetNumReactions(void);
///
    SF_Reaction* GetReaction(int);
///
    Boolean NeutralizerNeeded(void);
///
    void ComputeAlphaBulk();
///
    void InnerIterate(void);
///
    void IterateWithFixedPhiBulks(void);
///
    void residuals(double *const, double *const);
///
    void inneriteration(double *const, double *const, double);
///
    void ComputeResidues(void);
///
    double GetResidue(int);
///
    double GetItVar(int);
///
    void SetNeutralizingRatio(void);
///
    Boolean IterateAlphas(void);
///
    double ComputeAlphaResidue(SF_State*);
///
    void ComputeComplexPartOfChemPots(void);
///
    void UpdateReactionConstants(double,double);
  private:
///
    SF_MoleculeList* MolQ;
///
    SF_SegmentList* SegQ;
///
    Input* MyInput;
///
    int numReactions;
///
    Array<SF_Reaction*> ReactionQ;
///
    int numDiffStates;
///
    Array<SF_State*> StateQ;
///
    int numSegments;
///
    int numStateRatios;
///
    Array<SF_StateRatio*> StateRatioQ;
///
    int numSegmentRatios;
///
    Array<SF_StateRatioList*> SegmentRatioQ;
///
    Boolean neutralizerNeeded;
///
    SF_StateRatio* NeutralizingRatio;
///
    int numItVar;
///
    Vector x;
///
    Vector residues;
/// 
	Boolean forcePseudohessian;
///
    Array<int> residualMap;
///
    double smallAlpha;
///
    int smallAlphaCount;
///
    int maxNumSmallAlpha;
///
	double oldAccuracy;
///
    Boolean restart;
///
    int restartCount;
///
    Boolean phiBulksFixed;
///
    Boolean internalFreeEnergiesGiven;
///
    void CopyItVar(void);
///
    void SetItVar(void);
///
    void MapResiduals(void);
///
    void MapItVar(void);
///
    void SolveConflict(int);
///
    int GetMaxResponse(int);
///
    void FillResponseMatrix(Matrix);
///
    void SetNeutralizerNeeded(void);
///
    void ComputeInternalFreeEnergies(void);
///
    SF_Segment* FindWater(void);
///
    void ErrorOutput(void);
/// 
    void ComputeDerivatives(Vector,Vector,int,double);
/// 
    void SetPhiBulk(SF_Molecule*, const Vector, const Vector, double);
};

#endif
