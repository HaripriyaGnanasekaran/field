#ifndef SF_MOLSEGMENTxH
#define SF_MOLSEGMENTxH

#include "SF_MolState.h"
#include "LatticeRange1D.h"
///
class SF_SegmentList;

class SF_MolSegment: public Link{
  public:
///
	SF_MolSegment(Text,Text,Lattice*,Input*);
///
	virtual ~SF_MolSegment(void);

///
	Text GetName(void) const;
///
	Text GetMolName(void) const;
 ///
 	void GetOutput(Output*) const;
///
	double GetPhiBulk(void) const;
///
	void SetPhiBulk(double);
///
	Vector GetPhiBulkBoundaries(void) const;
///
	void SetPhiBulkBoundaries(Vector);
///
	double GetPhiRef(void) const;
///
	void SetPhiRef(double);
///
	SegmentFreedom GetFreedom(void) const;
	Boolean GetMC(void) const;
	Boolean GetLD(void) const;
///
	double GetPhiGrafted(void) const;
///
	LatticeRange* GetLatRange(void) const;
///
	int GetNumStates(void) const;
///
	SF_MolState* GetState(int) const;
///
	Vector GetAlpha(int) const;
///
	double GetAlphaAv(int) const;
///
	Vector GetSWF(void) const;
///
	void DelPos(int,int,int);
	void UpdatePos(int,int,int);
	void UpdatePos(double,double,double,double*);
	void UpdateSWF(void);
	void ClearAllPos(void);
///
	void UpdatePhiBulk(void);
///
	void UpdatePhiStates(void);
///
	Vector GetPhi(DensityPart) const;
///
	void CreatePhi(DensityPart);
///
	void DeletePhi(DensityPart);
	//Boolean MayerSaupeSet() const;
///
	Array<DensityPart> GetDensityPartQ() const;
  private:
///
	Text name;
///
	Text molName;
///
	SegmentFreedom freedom;
	bool MC,LD;
	//Boolean Mayer_Saupe;

///
	LatticeRange* LatRange;
///
	double phiBulk;
///
	Vector phiBulkBoundaries;
///
	double phiRef;
///
	Vector pka;
///
	Vector swf;
///
	Vector swfFull;
///
	int numStates;
///
	Array<SF_MolState*> StateQ;
///
	int numPhi;
///
	Array<Vector> PhiQ;
///
	Array<DensityPart> DensityPartQ;
///
	Lattice* Lat;
};

#endif
