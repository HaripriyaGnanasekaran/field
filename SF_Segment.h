#ifndef SF_SEGMENTxH
#define SF_SEGMENTxH

#include "SF_MolSegment.h"
#include "SF_State.h"

///
class SF_Segment {
///
	friend class SF_SegmentList;
  public:
///
	SF_Segment(void);
///
	SF_Segment(Text,Array<Text>,Array<Text>,Lattice*,Input*);
///
	~SF_Segment(void);

///
	Text GetName(void) const;
///
	void GetOutput(Output*) const;
///
	Vector GetPhi(void) const;
///
	double GetPhiBulk(void) const;
///
	void UpdatePhiBulk(void);
///
	SegmentFreedom GetFreedom(void) const;
///
	LatticeRange* GetLatRange(void) const;
///
	int GetNumStates(void) const;
///
	SF_State* GetState(int) const;
///	void DeleteState(int);
	Vector GetSWF(void) const;
///
	void UpdateSWF(void);
	void ResetRhoSet(void);
	void UpdatePos(int);
	void UpdatePos(int,double);

	void UpdatePos(double,double,double,double*);
	void DelPos(int);
	int SpotSize;
	int SpotType;

	bool GetMC(void) const;
	double GetValence(void) const;
	bool GetLD(void) const;
	//double GetMayerSaupe() const;
	int GetSpotSize(void) const;
	int GetSpotType(void) const;
	bool MoveAllowed(int) const;
	void ClearAllPos();
///
	void UpdatePhiStates(void);
///
	int GetNumMolSegments(void) const;
///
	SF_MolSegment* GetMolSegment(int) const;
///
	SF_MolSegment* GetMolSegment(Text) const; ///Text molName
///
	SF_MolSegment* NewMolSegment(Text);  ///Text molName

  private:
	Text name;
///
	SegmentFreedom freedom;
	bool MC,LD,moveX,moveY,moveZ;
	//Boolean Mayer_Saupe;
///
	LatticeRange* LatRange;
///
	double phiBulk;
///
	Vector pka;
///
	Vector swf;
///
	int numStates;
///
	Array<SF_State*> StateQ;
///
	Head MolSegmentQ;
///
	Lattice* Lat;
///
	Input* MyInput;
	double valence;
};

#endif
