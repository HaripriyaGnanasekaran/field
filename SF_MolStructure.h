#ifndef SF_MOLSTRUCTURExH
#define SF_MOLSTRUCTURExH

#include "SF_SegmentBlock.h"
#include "SF_SegmentList.h"
#include "SF_ParserBranched.h"
///
//! \brief SF_MolStructure is used to process molecular compositions.

/*!	A linear composition is defined by a string of segment names between brackets followed by numbers.
	Fragments enclosed in square brackets are considered as branches connected to the preceding segment.
	Consequent fragments in square brackets without a segment in between are used to represent multiple
	branches connected to the segment preceeding the first fragment. Any composition allowed for linear
	or branched chains may be used as a branch.

	Examples:
		mol : surfactant : composition : (C)12((O)1(C)2)8(O)1
		mol : ethanol : composition : (H)1(C)1[(H)1][(H)1](C)1[(H)1][(H)1](O)1(H)1
		mol : dendrimer : composition : @dend(3,(ux)1(u)2,3,(br)1(u)2,2,(br)1(u)2(ue)1,3)
*/
class SF_MolStructure {
  public:
///
	SF_MolStructure(Text,
					Text,
					SF_SegmentList*,
					Lattice*,
					Input*);
///
	~SF_MolStructure();
///
	double GetAvLength() const;
        bool found;	
///
	int GetMaxLength() const;

//! Returns a pointer to the segment number r in SegBlockQ[1], works for linear chains.
	SF_MolSegment* GetSegment(int) const;

//! Returns a pointer to the segment number r in the link number linknum, the segments are numbered starting from the parent node.
	SF_MolSegment* GetSegment(int linknum, int r) const;
///
	int GetNumDiffSegments() const;
///
	SF_MolSegment* GetDiffSegment(int) const;
///
	int GetNumSegmentBlocks(void) const;

	Text GetComposition(void) const;
///
	SF_SegmentBlock* GetSegmentBlock(int) const;
	SF_SegmentBlock* GetSegmentBlock(int i, int j) const;
	SF_SegmentBlock* GetSegmentBlockComb(int i) const;
///
	double GetAvNumSegments(const SF_MolSegment*) const;

//!	Warning: rewritten to use MolType field.
	MoleculeType GetMoleculeType() const;
///
	void SetPhiBulk(double);
///
	void SetPhiBulkBoundaries(Vector);
///
	double ChemIntRef() const;

//!	Note: returns false for non-linear chains
	Boolean Symmetric() const;
///
	Boolean SomeSegmentsPinned() const;
	Boolean SomeSegmentsMC() const;
	Boolean SomeSegmentsLD() const;
///
	Boolean SomeSegmentsGrafted() const;
///
	SF_MolSegment* GetGraftedSegment() const;

///
	Boolean AllInteractionsEqual() const;
///
	int GetNumGenerations() const;
	int GetNumArms() const;
///
	int GetNumRepeat(int) const;
	int GetSideType() const;


//!	Returns a pointer to SF_Molstructure::LinkNodeTree
	SF_LinkNodeTree* GetLinkNodeTree() const;


	//MoleculeType MolType;
  protected:

	SF_LinkNodeTree* LinkNodeTree;


  private:

//!	A helper class, reads branched compositions, builds LinkNodeTree and SegBlockQ.
	friend class SF_ParserBranched;

//!	Builds DiffSegmentQ, used to be a part of the constructor dealing with dendrimers
	void SetDiffSegmentQ();

//!	Calculates the molecule length, used to be a part of the constructor dealing with dendrimers
	//void SetMolLength();

//!	Sets PhiRef for the segments in the DiffSegmentQ, used to be a part of the constructor dealing with dendrimers
	void SetPhiRef();

///
	Text molName;
///
	double molLength;
	int side_type;
	//bool COMB;
///
	int numDiffSegments;
///
	Text composition;
///
	double chemIntRef;
///
	Array<SF_SegmentBlock*> SegBlockQ;
	Array<SF_SegmentBlock*> SegBlock_1Q;
	Array<SF_SegmentBlock*> SegBlock_2Q;
	Array<SF_SegmentBlock*> SegBlock_combQ;

///
	Array<int> repeatQ;
///
	int numGenerations;
	int numArms;
	int lin,den,a_den,bra;
///
	Lattice* Lat;
///
	SF_SegmentList* SegQ;
	Array<Text> AliasNames;
	Array<Text> AliasComps;
	int numAlias;
///
	MoleculeType MolType;
	void SetMolLength();
///
	Array<SF_MolSegment*> DiffSegmentQ;
///
	Boolean someSegmentsPinned;
	Boolean someSegmentsMC;
	Boolean someSegmentsLD;

///
	Boolean someSegmentsGrafted;
///
	SF_MolSegment* GraftedSegment;

};

#endif
