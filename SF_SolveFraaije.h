#ifndef SF_SOLVEFRAAIJExH
#define SF_SOLVEFRAAIJExH

#include "SF_Solve.h"

///
class SF_SolveFraaije : public SF_Solve {
  public:
///
	SF_SolveFraaije(Boolean, 
			 SF_ReactionList*, 
			 SF_MoleculeList*, 
			 SF_SegmentList*, 
			 SF_ReactionList*, 
			 SF_MoleculeList*, 
			 SF_SegmentList*, 
			 Lattice*, 
			 Input*);
///
	virtual ~SF_SolveFraaije();
	
///
	Text GetName();
///
	void GetOutput(Output*) const;
///
	virtual void Iterate(Boolean,double);
///	void WriteInitialGuessToFile(void);
	virtual void GetResiduals(Vector) const ;
///
	void UpdateDensities();
///
	virtual double CalcStop() const;
///
	virtual double CalcError() const;
  protected:
///
	int numBulkBounds;
///
	Boolean firstIterate;
///
	SF_ReactionList* ReactionNewQ;
///
	SF_MoleculeList* MolNewQ;
///
	SF_SegmentList* SegNewQ;
///
	Array<Vector> SWFNew;
///
	void UpdateItVar(Boolean);
///
	void UpdateSWF();
///
	Array<SF_MolState*> StateQ;
///
	Array<SF_MolState*> StateNewQ;
///
	Array<Vector> dPhi;
///
	Array<Vector> dPhiNew;
///	void SetInitialGuessFromFile(void);
	void ComputePhi(void);
///
	void CalculateDPhiNew(void);
///
	Vector CalcFluxState(const SF_MolState*, const int);
///
	Vector CalcFluxTwoStates(const Vector, const Vector, const SF_MolState*);
///
	Vector CalcUPrime(const SF_MolState*);
///
	void SetInitialGuessFromPrevious(const Lattice*, const SF_SegmentList*, const SF_MoleculeList*);
///
	double timeStep;
///
	void SetFrozenSWFZero(void);
///
	Boolean BulkItVar(const int) const;
///
	void residuals(double *const, double *const);
  private:
};

#endif
