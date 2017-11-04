#ifndef SF_SEGMENTLISTxH
#define SF_SEGMENTLISTxH

#include "SegmentParam.h"
#include "SF_Segment.h"
#include "Input.h"

///
class SF_SegmentList {
  public:
///
	SF_SegmentList(Lattice*,Input*);
///
	~SF_SegmentList();
///
	void GetOutput(Output*) const;
///
	Boolean SegmentDefined(const Text) const;
///
	int GetNumSegments(void) const;

///
	SF_Segment* GetSegment(const int) const;
///
	SF_Segment* GetSegment(const Text) const;
///
	SF_Segment* GetBaseSegment(const SF_State*) const;
///
	SF_Segment* GetBaseSegment(const SF_MolSegment*) const;
///
	int GetNumStates(void) const;
///
	Boolean StateDefined(const Text) const;
///
	SF_State* GetState(const int) const;
///
	SF_State* GetState(const Text) const;
///
	SF_MolSegment* NewMolSegment(Text,Text); /// segname,molname
	double ChemInt(const SF_Segment*, const int) const; ///int layer
	double ChemInt(const SF_MolSegment*, const int) const; ///int layer
	double ChemInt(const SF_State*, const int) const; ///int layer
	double ChemInt(const SF_MolState*, const int) const; ///int layer
	double ChemIntBulk(const SF_Segment*) const;
///
	double ChemIntBulk(const SF_MolSegment*) const;
///
	double ChemIntBulk(const SF_State*) const;
///
	double ChemIntBulk(const SF_MolState*) const;
///
	double ChemIntRef(const SF_MolSegment*) const;
///
	double ChemIntRef(const SF_MolState*) const;
///
	void PutSides(const SF_State*) const;
///
	Boolean ChemIntStatesEqual(const SF_State*, const SF_State*) const;
///
	void UpdateSWF(void);
///
	void UpdatePhiBulk(void);
///
	void UpdatePhiStates(void);
///
	Vector GetPhiTotal(void) const;
///
	Boolean Charged(void) const;
///
	Boolean ReactionsPresent(void) const;
///
	Vector GetCharge(void) const;
///
	Vector GetAverageEpsilon(void); ///const?
	Vector GetElectricPotential(void) const;
///
	Vector ComputeElectricPotential();
///
	Vector ComputeElectricPotential(Vector);
///
	Matrix GetChiMatrix() const;
	Vector GetPfVector() const; 
	Array <Text> GetPfx() const;
	Array <Text> GetPfy() const; 
		int Pf_count;
	int GetPf_count() const; 
	int GetChiParamNum(const SF_State*) const;
	

///
	void DeleteUnusedSegments(void);
  private:
///
	Input* MyInput;
///
	Lattice* Lat;
///
	int M;

///
	int numSegments;
///
	Array<SF_Segment*> SegmentQ;
///
	int numStates;
///

///
	Matrix chi;
	Vector Pf;
	Array<Text> Pfx;
	Array<Text> Pfy; 
///
	Text* SegStateNames;
///
	Boolean charged;
///
	Vector psi;
///
	int numPos;
///
	int numNeg;
///
	double posValence;
///
	double negValence;
    ///
    Boolean CheckName(Text) const;
///
	int GetNumStatesDefinedInSegment(const Text) const;
///
	Vector GetPosEnergy(void) const;
///
	Vector GetNegEnergy(void) const;
///
	double GetPosValence(void) const;
///
	double GetNegValence(void) const;
///
	void ComputeNumPos(void);
///
	void ComputeNumNeg(void);
///
	Array<SF_State*> StateQ;
///
	int GetChiParamNum(const SF_MolState*) const; 

};

#endif
