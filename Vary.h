#ifndef VARYxH
#define VARYxH

#include <math.h>
#include <limits.h>

#include <float.h>
#include "Message.h"
#include "Input.h"

///
const int maxNumCalculations = 10000000;

///
class Vary {
  public:
///
	Vary(Input*);
///
	~Vary();
///
	Boolean NextCalculation(void);
///
	void SkipToLast(void);
///
	Text GetVariableName(void);
///
	Text GetCurrValue(void);
/// 
  private:
///
	Input* MyInput;
///
	Text className;
///
	Text classInstance;
///
	Text parameter;
///
	Text outputName;
///
	int numCalculations;
///
	int currNumCalculation;
///
	int type;
///
	int scale;
///
	double startValue;
///
	double endValue;
///
	double step;
///
	double currValue;
};

#endif
