#ifndef SF_SOLVExH
#define SF_SOLVExH

#include "SF_MoleculeList.h"
#include "SF_ReactionList.h"
#include "LatticeRangeFile.h"

///
enum InitialGuessSymmetry {firstToLast,boundsToMiddle};

///
class SF_Solve : public SFNewton {
  public:
///
	SF_Solve() {};
///
	SF_Solve(Boolean,  SF_ReactionList*,  SF_MoleculeList*,  SF_SegmentList*,  Lattice*,  Input*);
///
	virtual ~SF_Solve();

///
	Text GetName() const;
///
	virtual void GetOutput(Output*) const;
///
	void Iterate(void);
	void ReComputePhi(bool);
///
	Boolean ErrorOccurred(void);
///
	void SetErrorOccurred(Boolean = true);
///
	virtual void SetInitialGuess(SF_Solve*);

	virtual void MoveInitialGuess(int* ,int* ,int);
///
	virtual Boolean NewSolveNeeded(SF_Solve*);
///
	virtual void UpdateSolve(SF_Solve*);
///
	virtual void WriteInitialGuessToFile(void) const;
Array<bool> ItVarArray;
///
	int GetNumItVar(void) const;
///
	Vector GetItVar(void) const;
///
	virtual void GetResiduals(Vector) const;
///
	double GetTolerance(void) const;

///
	void Jump(int); // All iteration variables are shifted over a given number of layers
/////
//	Boolean GetOverflowProtectionNeeded(void) const;
/////
//	Boolean Restart(void) const;
  protected:
///
	Input* MyInput;
///
	SF_ReactionList* ReactionQ;
///
	SF_MoleculeList* MolQ;
///
	SF_SegmentList* SegQ;
///
	Lattice* Lat;
///
	Text name;
///
	Boolean errorOccurred;
	Boolean reset_hess;
	Boolean reset_pseudohessian;


///
	Text initialGuessInputFile;
///
	Text initialGuessOutputFile;
///
	Boolean extPotSet;
///
	Text externalPotentialFile;
///
	Vector psi0;
///
	Boolean potentialItVar;
///
	Boolean charged;
///
	Vector ESquared;
///
	Boolean diffEpsilon;
///
	Boolean multiState;
///
	InitialGuessSymmetry initGuessSym;
///
	Boolean writeInitialGuess;
///
	Boolean compute;
///
	double phiBulkSolvent;
///
	double phiBulkNeutralizer;
///
	Boolean neutralizerPresent;
///
	Boolean iterateAlphas;
///
	int oldNumIterations;
///
	int oldNumFunctions;
///
	int smallAlphaCount;
///
	double smallAlpha;
///
	int maxNumSmallAlpha;
///
	int reverseDirectionRange;
///
	double maxFrReverseDirection;
///
	Vector reverseDirection;
///
	double minAccuracyForHessian;
///
	int numIterationsForHessian;
	int Newtonmethod;
///
	int numIterationsSinceHessian;
///
	Boolean transferHessian;
///
	Boolean numIterationsUpdated;
///
	double minAccuracySoFar;
///
	double resetHessianCriterion;
///
	int numItVar;
///
	Vector x;
	Vector U_offset;
///
	int M;
///
	Boolean historyTest,Offset;
///
	int numDiffItVar;
///
	Array<int> NumStatesForItVar;
///
	Array<Vector> SWF;
///
	int numLiquidStates;
///
	Array<SF_State*> StateQ;
///
//	LatticeRangeFile ItVarRange;
	//Array<bool> ItVarArray;
/////
//	Boolean overflowProtectionNeeded;
/////
//	Boolean restart;
///
	void residuals(double *const, double *const);
	void residuals(Vector);
///
	void inneriteration(double *const, double *const, double);
	void SetUOffset();
///
	virtual void ComputePhi(void);
///
	virtual void UpdateSWF(void);
///
	virtual void UpdateItVar(Boolean);
///
	virtual void SetInitialGuessFromFile(void);
///
	virtual void SetInitialGuessFromPrevious(const Lattice*, const SF_SegmentList*, const SF_MoleculeList*);
///
	virtual void ComputeInitialGuess(void);
	void NoInitialGuess(void);
///
	virtual void ComputeAlphaBulk(void);
///
	virtual void SetExternalPotentials(void);
///
	virtual void SetFrozenSWFZero(void);
///
	void CopyInitialGuess(const Vector, Vector, const int) const;
///
	void CopyInitialGuess(const Vector, Vector, const int, const int) const;
///
	void CopyInitialGuess(const Vector, Vector, const int, const int, const int) const;
///
	void Copy1Dto1DFirstToLast(const Vector, Vector, const int) const;
///
	void Copy1Dto2DFirstToLast(const Vector, Vector, const int) const;
	void Copy1Dto3DFirstToLast(const Vector, Vector, const int) const;

///
	void Copy2Dto2DFirstToLast(const Vector, Vector, const int, const int) const;
///
	void Copy3Dto3DFirstToLast(const Vector, Vector, const int, const int, const int) const;
///
	void Copy2Dto1DFirstToLast(const Vector, Vector, const int, const int) const;
///
	void Copy1Dto1DBoundsToMiddle(const Vector, Vector, const int) const;
///
	void Copy1Dto2DBoundsToMiddle(const Vector, Vector, const int) const;
///
	void Copy2Dto2DBoundsToMiddle(const Vector, Vector, const int, const int) const;
///
	void Copy2Dto1DBoundsToMiddle(const Vector, Vector, const int, const int) const;
///
	void WriteOverflowError(void) const;
///
	Boolean SegItVar(const int) const;
///
	Boolean PotItVar(const int) const;
///
	virtual void ProcessInput();
///
	void CheckSolution(void);

};

inline
Boolean
SF_Solve::SegItVar(const int z) const {
	return ItVarArray[z];
}

#endif
