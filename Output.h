#ifndef OUTPUTxH
#define OUTPUTxH

#include "Input.h"
#include "OutputLine.h"
#include "Vary.h"
#include <fstream>
#ifdef _WIN32
typedef unsigned uint;
#endif
using namespace std;


///
	Text parseFileName(const Text &filename, const Input *);
	//static const int imageLength = 255;
///
class Output {
  public:
///
	Output(Input*,Vary*,Array<int>,int,Array<Text>);
///
	~Output();

///
	void PutText(const Text, const Text, const Text, const Text);
///
	void PutBoolean(const Text, const Text, const Text, const bool);
///
	void PutInt(const Text, const Text, const Text, const int);
///
	void PutReal(const Text, const Text, const Text, const double);
///
	void PutProfile(const Text, const Text, const Text, const Vector);
///
	void PutVector(const Text, const Text, const Text, const Vector, const int, const int = 1);
///
	void WriteOutput(void);
	void WriteMCOutput(int);
///
	void Clear(void);
///
	Array<Text> GetFileNames(void);
///
	void SetErrorOccurred(void);
	int counter;
  protected:
///
	void OutputAna(Text, bool);
///
	void OutputKal(Text, bool);
	void OutputAve(Text, bool);
///
	void OutputProfiles(Text);
	void OutputVectors(Text);
	void OutputVtkProfiles(Text);
///
	void OutputTemplate(Text, bool);
///
	bool ValueInTemplate(const OutputLine*);
///
	bool ValueForLayerInTemplate(const OutputLine*, int);
///
	Text getNumberedFilename(const Text) const;
///
	Input* MyInput;
///
	Vary* MyVary;
///
	Input* Template;
///
	bool templateSet;
///
	int accuracy,MCS,MC_output_interval;
///
	bool writeBounds;
///
	bool append;
	bool blank_line;
///
	bool write;
///
	bool outputProfiles;
///
	bool skipErrors;
///
	bool errorOccurred;
///

  private:
///
	Head OutputLineQ;
/// The number of gradients in the system.
	int gradients;
/// Array with the dimensions of each gradient (1, ...) in the system.
	Array<int> grads;
///
	int accuracyDefault;
///
	Array<Text> fileNames;
///
	Array<Text> fileNamesOld;
///
	Array<bool> firstWriteQ;
///
	bool FileExists(Text);
///
	void open(ofstream &, Text &filename, bool firstwrite=true);
	///
	Boolean CheckReal(Text) const;

};

#endif

