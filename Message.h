#ifndef MESSAGExH
#define MESSAGExH

#include "Input.h"

///
enum Error {fatal, warning, hint, debug,implementation, literal};
///
extern Boolean warnings;
///
extern Boolean silent;
///
extern Boolean silentForced;

///
void Message(Error,Input*,Text);
///
void Message(Error,Text);
///
Text GetVersion(void);
Text GetDate(void);

#endif
