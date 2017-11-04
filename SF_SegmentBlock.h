#include "SF_Box.h"
#include "SF_SegmentList.h"
//include "PolydispFlory.h"
///
#ifndef SF_SEGMENT_BLOCKxH
#define SF_SEGMENT_BLOCKxH
class SF_SegmentBlock : public Link {
  public:
///
	SF_SegmentBlock(Text,Text,SF_SegmentList*,Input*,Text,Text);
///
	virtual ~SF_SegmentBlock();

///
	Text GetBlockName(void) const;
///
	Text GetBlockRepeat(void) const;
///
	double GetAvLength(void) const;
///
	int GetMinLength(void) const;
///
	int GetMaxLength(void) const;
///
	double GetAvRepeat(void) const;
///
	int GetMinRepeat(void) const;
///
	int GetMaxRepeat(void) const;
///
	SF_MolSegment* GetSegment(int) const;
///
	int GetNumBlocks(void) const;
///
	SF_SegmentBlock* GetBlock(int) const;
///
	int GetNumDiffSegments(void) const;
///
	SF_MolSegment* GetDiffSegment(int) const;
///
	int GetAvNumSegments(const SF_MolSegment*) const;
///
	Boolean Symmetric(void) const;
///
	Boolean PolydisperseChain(void) const;
///
	Boolean Polydisperse(void) const;
///
	double GetFraction(int) const;
///
	Array<Text> GetFractionNames(void) const;
  private:
///
	Text blockName;
///
	Text nameBlockBefore;
///
	Text nameBlockAfter;
///
	Input* MyInput;
///
	SF_SegmentList* SegQ;
///
	double avLength;
///
	int minLength;
///
	int maxLength;
///
	double avRepeat;
///
	int minRepeat;
///
	int maxRepeat;
///
	int numBlocks;
///
	int numDiffSegments;
///
	Text repeated;
///
	Text repeat;
///
	Boolean polydisperse;
///	Polydispersity* Polydisp;
	Array<SF_SegmentBlock*> BlockQ;
///
	Array<SF_MolSegment*> DiffSegmentQ;
///
	int CountNumBlocks(Text&,Text) const;
///
	void AssignOneBlock(Text,Text);
///
	void AssignMultiBlock(Text,Text);
///
	void AssignDiffSegmentQ(void);
///
	Text ReadRepeat(Text&,Text) const;
///	Polydispersity* NewPolydispersity(Text, Text);
	Boolean operator==(const SF_SegmentBlock&) const;
///
	Boolean operator!=(const SF_SegmentBlock&) const;
};
#endif
