#ifndef SEGMENTPARAMxH
#define SEGMENTPARAMxH

#include "SF_Box.h"

///
class SegmentParam {
  public:
///
	SegmentParam();
///
	SegmentParam(Text,Input*, Text);
///
	~SegmentParam();
	
///
	Text Name() const;
///
	int NumParam() const; 
///
	Text Param(int) const; 
///
	int NumStates() const;
    ///
    Text StateName(int) const;
    ///
    SegmentParam* StatePTR(int) const;
    ///
    int NumChiParam() const;
    ///
    Text ChiParam(int) const;
///
	Boolean ParamSet(Text) const;
    ///
    Boolean ChiParamSet(Text) const;
  private:
///
	Text name;
///
	int numStates;
///
	Array<SegmentParam*> States;
///
	int numParam;
///
	Array<Text> param;
///
	int numChiParam;
///
	Array<Text> chiParam;
///
    Boolean ChiParameter(Text) const; 
///
    Text StripChi(Text) const; 
};

#endif
