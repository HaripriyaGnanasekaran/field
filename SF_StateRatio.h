#ifndef SF_STATERATIOxH
#define SF_STATERATIOxH

#include "SF_State.h"

///
class SF_StateRatio {
  public:
///
	SF_StateRatio(SF_State*, SF_State*);
///
	~SF_StateRatio();
///
	SF_State* GetState1(void) const;
///
	SF_State* GetState2(void) const;
///
	double GetAlphaBulkRatio(void) const;
///
	void SetAlphaBulkRatio(double);
///
	void ComputeAlphaBulkRatio();
 ///
 	Boolean operator==(const SF_StateRatio&) const;
 ///
 	Boolean SameStates(const SF_StateRatio&) const;
 private:
///
	SF_State* State1;
///
	SF_State* State2;
///
	double value;
};

#endif
