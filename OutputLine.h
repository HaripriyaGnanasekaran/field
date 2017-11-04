#ifndef OUTPUTLINExH
#define OUTPUTLINExH

#include "Message.h"

///
class OutputLine : public Link {
  public:
///
	OutputLine();
///
	OutputLine(const Text,const Text,const Text,const Text);
///
	OutputLine(const Text,const Text,const Text,const Text,const Vector);
///
	OutputLine(const Text,const Text,const Text,const Text,const Vector,const int, const int);
///
	virtual ~OutputLine(void);
///
	Text GetElement1(void) const;
///
	Text GetElement2(void) const;
///
	Text GetElement3(void) const;
///
	Text GetValue(void) const;
///
	bool IsProfile(void) const;
///
	Vector GetProfile(void) const;
///
	bool IsVector(void) const;
///
	Vector GetVector(void) const;
///
	int GetMax(void) const;
///
	int GetMin(void) const;
  private:
///
	Text Element1;
///
	Text Element2;
///
	Text Element3;
///
	Text Value;
///
	bool profile;
///
	bool vector;
///
	Vector profileValue;
///
	int max;
///
	int min;
};

#endif

