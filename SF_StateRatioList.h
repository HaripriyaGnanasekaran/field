#ifndef SF_STATERATIOLISTxH
#define SF_STATERATIOLISTxH

#define _GNU_SOURCE 1
#ifdef USE_FENV
#include <fenv.h>
#endif

#include "SF_StateRatio.h"
#include "SF_SegmentList.h"

///
class SF_StateRatioList {
  public:
///
	SF_StateRatioList(Array<SF_StateRatio*>,SF_SegmentList*);
///
	~SF_StateRatioList(void);
///
	void UpdateAlphaBulk(void);
  private:
///
	int numStateRatios;
///
	Array<SF_StateRatio*> StateRatioQ;
///
	SF_SegmentList* SegQ;
///
	SF_Segment* Segment;
};

#endif
