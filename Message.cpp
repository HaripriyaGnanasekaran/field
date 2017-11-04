#include "Message.h"

Boolean warnings;
Boolean silent;
Boolean silentForced;

void
Message(Error error, Input* MyInput, Text mess) {
	if (error == fatal) Sysout().Outtext("Fatal error ");
	if (error == warning && !silent) Sysout().Outtext("Warning ");
	if (error == hint && !silent) Sysout().Outtext("For your information ");
	if (error == debug) Sysout().Outtext("Debug message ");
	if (error == implementation && !silent) Sysout().Outtext("Missing function:");
	if (!silent && error != fatal) {
		Sysout().Outimage();
	}
	if (error == fatal || (error == warning && !silent)) {
		Sysout().Outtext("reading input file '" + MyInput->GetFileName() + "', ");
		Sysout().Outtext("in calculation number ");
		Sysout().Outint(MyInput->GetCurrNumCalculation(),0);
		Sysout().Outimage();
	}
	if (!silent || error == fatal || error == debug) {
		Sysout().Outtext(mess);
		Sysout().Outimage();
		Sysout().Outimage();
	}
	if (error == fatal) exit(-1);
}

void
Message(Error error, Text mess) {
	if (error == fatal) Sysout().Outtext("Fatal error ");
	if (error == warning && !silent) Sysout().Outtext("Warning ");
	if (error == hint && !silent) Sysout().Outtext("For your information ");
	if (error == debug) Sysout().Outtext("Debug message ");
	if (error == implementation && !silent) Sysout().Outtext("!!!!Missing function: ");
	if (!silent || error == fatal || error == debug || error == literal) {
		Sysout().Outimage();
		Sysout().Outtext(mess);
		Sysout().Outimage();
		Sysout().Outimage();
	}
	if (error == fatal) exit(-1);
}

Text
GetVersion(void) {
	return "2.1.0.0";
}

Text
GetDate(void) {
	return "2015-7913";
}
