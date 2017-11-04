#ifndef INPUTxH
#define INPUTxH

#include <fenk.h>
#include <stdlib.h>
#include <ctype.h>

///
static const int imageLength = 255;

///
class Input {
  public:
///
	Input(Text);
///
	~Input();

///
	Text GetFileName() const;

///
	void SetAllMessagesOn();
///
	void SetAllMessagesOff();

///
	void SetDefaultWarningsOn();
///
	void SetDefaultWarningsOff();

///
	Boolean NextCalculation();
///
	int GetCurrNumCalculation() const;

///
	void CheckFirstArguments(Array<Text>) const ;

///
	int GetNumNames(Text) const;
///
	int GetNumNames(Text,int,int) const;
///
	Array<Text> GetNames(Text) const ;

///
	Array<Text> GetParameters(Text,Text) const;
///
	void CheckParameterNames(Text, Text, Array<Text>) const;
///
	void DontCombineParam(Text, Text, Text, Text) const;
///
	void AlwaysCombineParam(Text, Text, Text, Text) const;

///
	Boolean ValueSet(Text,Text,Text) const;
	/// checks whether given parameter is set to given value, accepts more than
	/// one parameter in input (for template files)
	Boolean ParamSetToValue(Text,Text,Text,Text) const;
	/// accepts wild card '*' in input file (for template files)
	Boolean ValueSetWild(Text,Text,Text) const;
	/// only used in class Vary so far
	int LastNumCalcValueSet(Text,Text,Text) const;

	// no default set in following functions
	// if no value is found or if the value is outside the bounds,
	/// an error message is generated
	int GetInt(Text, Text, Text, int, int) const;
///
	double GetReal(Text, Text, Text, double, double) const;
///
	Text GetText(Text, Text, Text) const;
///
	int GetChoice(Text, Text, Text, Array<Text>) const;
///
	Boolean GetBoolean(Text, Text, Text) const;

	// following functions provide a default as last argument
	// if the value is outside the bounds, an error message is generated
	/// if no value is found the default is returned
	int GetInt(Text, Text, Text, int, int, int) const;
///
	double GetReal(Text, Text, Text, double, double, double) const;
///
	Text GetText(Text, Text, Text, Text) const;
///
	int GetChoice(Text, Text, Text, Array<Text>, int) const;
///
	Boolean GetBoolean(Text, Text, Text, Boolean) const;
///
	unsigned int TotalNumCalculations() const;

	/// last argument is the new value of the variable
	void SetVariable(Text, Text, Text, Text);
Boolean CheckInt(Text) const;
  private:
///
	Text fileName;
///
	Text (*input)[5];
///
	int lastLineNr;
///
	int numberOfCalculations;
///
	int lastLineSystem;
///
	int firstLineSystem;
///
	int largestMessNumToScreen;
///
	int currNumCalculation;
///
	void Message(int, Text, int) const;
///
	Boolean CheckReal(Text) const;
///

///
	int GetNumParameters(Text,Text) const;
///
	int GetLineNr(Text,Text,Text) const;
///
	Boolean WildCompare(const Text wild, const Text full) const;
	/// following members private to avoid implementing them
	Input();
///
	Input(const Input&);
///
	Input& operator=(const Input&);
};
#endif
