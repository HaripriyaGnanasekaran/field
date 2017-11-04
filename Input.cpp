#include "Input.h"

Input::Input() {
	input = new Text[1][5];
}
Input::Input(Text t) {
	firstLineSystem = 0;
	lastLineSystem = 0;
	currNumCalculation = 0;
	largestMessNumToScreen = 3;
	fileName = t;
	Text lineContents;
	int elementNr, lineNr;
	const char test = ':';
	Text empty = "";
	Boolean more; 
	Infile file(t);
	if (!file.Open(Blanks(imageLength)) ) {
		Message(-1,"file '" + fileName + "' can't be opened or the file may not exist in this directory",-1);
	}
	lastLineNr = 1;
	numberOfCalculations=1;
	more = file.Inrecord();
	Boolean noDataSinceStart=false;
	while (!file.Endfile()) {
		lineContents=notext;
		while (more) {
			lineContents = lineContents + file.Image;
			more = file.Inrecord();
		}
		lineContents =
			lineContents + file.Image.Sub(1,file.Image.Pos()-1).Strip();
		more = file.Inrecord();
		lastLineNr++;
		if (lineContents.Length() > 4) { 
			if  (*lineContents.Sub(1,5) == *Copy("start")) {
				noDataSinceStart=true;
			} else if (lineContents.Length() > 0 && noDataSinceStart == true) {
				noDataSinceStart=false;
				numberOfCalculations++;
			}
		}
	}
	if (lastLineNr == 1) Message(-1,"inputfile is empty",-1);
	file.Close();
	input = new Text[lastLineNr+1][5];

	file.Open(Blanks(imageLength));
	more = file.Inrecord();
	lineNr = 0;
	while (!file.Endfile()) {
		lineNr++;
		lineContents = notext;
		while ( more ) {
			lineContents = lineContents + file.Image;
			more = file.Inrecord();
		}
		lineContents =
			lineContents + file.Image.Sub(1,file.Image.Pos()-1).Strip();

		Text dummy = "";
		Boolean endLine = false;
		while(lineContents.More()) {
			char test2 = lineContents.Getchar();
//mmm temporary allow '#' comments
if (test2 == '#') { 
  endLine = true;
} else {
			if (test2 == '/' && lineContents.More()) {
				test2 = lineContents.Getchar();
				if (test2 == '/') {
					endLine = true;
				} else if (!endLine) {
					Text dummy2 = Blanks(1);
					dummy2.Putchar(test2);
					dummy = dummy + "/" + Copy(dummy2.Frontstrip().Strip());
				}
			} else if (!endLine) {
				Text dummy2 = Blanks(1);
				dummy2.Putchar(test2);
				dummy = dummy + Copy(dummy2);
			}
}
		}
		lineContents = Copy(dummy);
		elementNr = 0;
		if (!lineContents.More()) {
			input[lineNr][1] = Copy("");
		}
		while (lineContents.More()) {
			elementNr++;
			if (elementNr > 4) Message(-1, "line contains too many elements", lineNr);

			input[lineNr][elementNr] = Copy(lineContents.Scanto(test));
			input[lineNr][elementNr] = Copy(input[lineNr][elementNr].Strip());
			input[lineNr][elementNr] = Copy(input[lineNr][elementNr].Frontstrip());
		}
		if (elementNr == 1 && *input[lineNr][1] == *Copy("start")); // ok
		else if (elementNr < 4 && elementNr != 0) Message(-1, "line contains too few elements", lineNr);
		if (*input[lineNr][1] == *empty && *input[lineNr][2] == *empty && *input[lineNr][3] == *empty && *input[lineNr][4] == *empty) ; // ok, empty line
		else if (*input[lineNr][1] == *empty || *input[lineNr][2] == *empty || *input[lineNr][3] == *empty || *input[lineNr][4] == *empty)
			if (*input[lineNr][1] != *Copy("start")) Message(-1, "statement missing", lineNr); //not ok
		more = file.Inrecord();
	}
	// check comments
	const char comment = '/';
	char check;
	Text c;
	for (lineNr = 1; lineNr <= lastLineNr; lineNr++) {
		c = Copy(input[lineNr][1]);
		if (c.Length() > 1) {
			check = c.Getchar();
			if (check == comment) {
				check = c.Getchar();
				if (check == comment) {
					input[lineNr][1] = Copy("");
					input[lineNr][2] = Copy("");
					input[lineNr][3] = Copy("");
					input[lineNr][4] = Copy("");
				}
			}
		}
	}
}
Input::~Input() {
	delete[] input;
}
Text Input::GetFileName() const {
	return fileName;
	
}
void Input::SetAllMessagesOff() {
	largestMessNumToScreen = -1;
}
void Input::SetAllMessagesOn() {
	largestMessNumToScreen = 3;
}

void Input::SetDefaultWarningsOff() {
	if (largestMessNumToScreen > 2) {
		largestMessNumToScreen = 2;
	}
}
void Input::SetDefaultWarningsOn() {
	if (largestMessNumToScreen == 2) {
		largestMessNumToScreen = 3;
	}
}
void
Input::CheckFirstArguments(Array<Text> t1) const {
	int i,lineNr;
	Boolean error;
	Text mess;
	int l1 = t1.Upperbound() - t1.Lowerbound() + 1;
	mess = "First argument should be : ";
	if (l1 == 1) mess = mess + t1[t1.Lowerbound()];
	else {
		mess = mess +  t1[t1.Lowerbound()];
		for (i=2; i<l1; i++) mess = mess + ", " + t1[i];
		mess = mess + " or " + t1[t1.Upperbound()];
	}
	for (lineNr = 1; lineNr <= lastLineNr; lineNr++) {
		error = true;
		for (i = 1; i <= l1; i++) {
			if (*input[lineNr][1] == *Copy(t1[i])) error = false;
		}
		if (*input[lineNr][1] == *Copy("") ) error = false;
		if (*input[lineNr][1] == *Copy("start") ) error = false;
		if (error) {
			mess = mess + ", not " + input[lineNr][1] + ".\n";
			Message(-1, mess, lineNr);
		}
	}
}
Boolean
Input::NextCalculation() {
	int lineNr;
	if (lastLineSystem == lastLineNr) return false;
	firstLineSystem = lastLineSystem+1;
	for (lineNr = firstLineSystem; lineNr <= lastLineNr; lineNr++) {
		if (lineNr == lastLineNr) lastLineSystem = lineNr;
		if (*input[lineNr][1] == *Copy("start")) {
			lastLineSystem = lineNr;
			for (lineNr = lineNr+1; lineNr <= lastLineNr; lineNr++) {
				if (*input[lineNr][1] == *Copy("")) lastLineSystem++;
				else {
					currNumCalculation++;
					return true;
				}
			}
		}
	}
	currNumCalculation++;
	return true;
}
int
Input::GetCurrNumCalculation() const {
	return currNumCalculation;
}
int
Input::GetNumNames(Text comp, int minNr, int maxNr) const {
	int lineNr,count,i;
	Text mess, Number;
	Number = Blanks(9);
	Text *names = new Text[maxNr];
	Boolean duplicate = false;
	count = 0;
	for (lineNr = 1; lineNr <= lastLineSystem; lineNr++) {
		if (*input[lineNr][1] == *comp) {
			duplicate = false;
			for (i=1; i<=count; i++) {
				if (*names[i-1] == *input[lineNr][2]) duplicate = true;
			}
			if (!duplicate) {
				count++;
				if (count > maxNr) {
					mess = "Too many instances of class '" + comp
					+ "' defined \n" + "Maximum number is: ";
					Number.Putint(maxNr);
					Message(-1, mess + Number, lineNr);
				} else names[count-1] = input[lineNr][2];
			}
		}
	}
	if (count < minNr) {
		mess = "Too little instances of class '" + comp + "' defined.\nDefine at least ";
		Number.Putint(minNr);
		Number = Copy(Number.Strip());
		Number = Copy(Number.Frontstrip());
		Message(-1, mess + Number, 0);
	}
	delete[] names;
	return count;
}
int
Input::GetNumNames(Text comp) const {
	int lineNr,count,i;
	Text mess, Number;
	Number = Blanks(9);
	int maxNr = lastLineNr;
	Text *names = new Text[maxNr];
	Boolean duplicate = false;
	count = 0;
	for (lineNr = 1; lineNr <= lastLineSystem; lineNr++) {
		if (*input[lineNr][1] == *comp) {
			duplicate = false;
			for (i=1; i<=count; i++) {
				if (*names[i-1] == *input[lineNr][2]) duplicate = true;
			}
			if (!duplicate) {
				count++;
				names[count-1] = input[lineNr][2];
			}
		}
	}
	delete[] names;
	return count;
}
Array <Text>
Input::GetNames(Text comp) const {
	int lineNr,count,i;
	Text mess, Number;
	int number = GetNumNames(comp);
	Array<Text> names;
	names.Dim(1,number);
	Boolean duplicate = false;
	count = 0;

	for (lineNr = 1; lineNr <= lastLineSystem; lineNr++) {
		if (*input[lineNr][1] == *comp) {
			duplicate = false;
			for (i=1; i<=count; i++) {
				if (*Copy(names[i]) == *input[lineNr][2]) duplicate = true;
			}
			if (!duplicate) {
				count++;
				if (count > number) {
					mess = "Too many instances of class '" + comp +
						"' defined.\nMaximum number is: ";
					Number = Blanks(9);
					Number.Putint(number);
					Number = Copy(Number.Frontstrip());
					mess = mess + Number +
					"\nCould also be a programming error in call to Input::GetNames";
					Message(-1, mess, lineNr);
				} else names[count] = input[lineNr][2];
			}
		}
	}
	return names;
}
Array <Text>
Input::GetParameters(Text comp, Text name) const {
	int lineNr,count,i;
	Text mess, Number;
	Array<Text> names;
	int number = GetNumParameters(comp,name);
	names.Dim(1,number);
	Boolean duplicate = false;
	count = 0;
	for (lineNr = 1; lineNr <= lastLineSystem; lineNr++) {
		if (*input[lineNr][1] == *comp && *input[lineNr][2] == *name) {
			duplicate = false;
			for (i=1; i<=count; i++) {
				if (*Copy(names[i]) == *input[lineNr][3]) duplicate = true;
			}
			if (!duplicate) {
				count++;
				if (count > number) {
					mess = "Too many parameters in class '" + comp + " : " + name +
						"' defined.\nMaximum number is: ";
					Number = Blanks(9);
					Number.Putint(number);
					Number = Copy(Number.Frontstrip());
					mess = mess + Number +
						"\nCould also be a programming error in call to Input::GetParameters";
					Message(-1, mess, lineNr);
				} else names[count] = input[lineNr][3];
			}
		}
	}
	return names;
}
void
Input::CheckParameterNames(Text comp, Text name, Array<Text> choices) const {
	int lineNr, i;
	Text mess;
	Boolean error;
	int numChoices = choices.Upperbound();
	if (choices.Lowerbound() != 1) {
		Sysout().Outtext("Please use an array with lowerbound = 1 for call to Input::CheckParameterNames");
		Sysout().Outimage();
		exit(-1);
	}
	for (lineNr = firstLineSystem; lineNr <= lastLineSystem; lineNr++) {
		if (*input[lineNr][1] == *comp && *input[lineNr][2] == *name) {
			error = true;
			for (i=1; i<=numChoices; i++) {
				if (*input[lineNr][3] == *Copy(choices[i])) error = false;
			}
			if (error) {
				mess = "The parameter '" + input[lineNr][3] + "' is unknown in the current context: '";
				mess = mess + comp + " : " + name + " : " + "'.";
				if (numChoices == 1) mess = mess + "Only " + choices[1] + " can be set.";
				else {
					mess = mess + "Choose from ";
					mess = mess + choices[1];
					for (i=2; i<=numChoices-1; i++) mess = mess + ", " + choices[i];
					mess = mess + " or " + choices[numChoices];
				}
				Message(-1, mess, lineNr);
			}
		}
	}
}
void
Input::DontCombineParam(Text comp, Text name, Text param1, Text param2) const {
	int lineNr;
	int lineNrError=0;
	Text mess, value1;
	Boolean param1Set = false;
	Boolean param2Set = false;
	for (lineNr = firstLineSystem; lineNr <= lastLineSystem; lineNr++) {
		if (*input[lineNr][1] == *comp && *input[lineNr][2] == *name) {
			if (*input[lineNr][3] == *Copy(param1)) {
				value1 = input[lineNr][4];
				param1Set = true;
			}
			if (*input[lineNr][3] == *Copy(param2)) {
				param2Set = true;
				lineNrError = lineNr;
			}
		}
	}
	if (param1Set && param2Set) {
		mess = "Don't set " + comp + " : " + name + " : " + param2 +
			" when "  + comp + " : " + name + " : " + param1 +
			" is set to to " + value1;
		Message(-1,mess,lineNrError);
	}
}
void
Input::AlwaysCombineParam(Text comp, Text name, Text param1, Text param2) const {
	int lineNr;
	//int lineNrError;
	Text mess, value1;
	Boolean param1Set = false;
	Boolean param2Set = false;
	for (lineNr = firstLineSystem; lineNr <= lastLineSystem; lineNr++) {
		if (*input[lineNr][1] == *comp && *input[lineNr][2] == *name) {
			if (*input[lineNr][3] == *Copy(param1)) {
				value1 = input[lineNr][4];
				param1Set = true;
			}
			if (*input[lineNr][3] == *Copy(param2)) {
				param2Set = true;
				//lineNrError = lineNr;
			}
		}
	}
	if (param1Set && !param2Set) {
		mess = "Always set " + comp + " : " + name + " : " + param2 +
			" when "  + comp + " : " + name + " : " + param1 +
			" is set to to " + value1;
		Message(-1,mess,0);
	}
}
Boolean
Input::ValueSet(Text comp, Text name, Text param) const {
	int lineNr;
	Boolean duplicate = false, present = false;
	Text mess;
	for (lineNr = firstLineSystem; lineNr <= lastLineSystem; lineNr++) {
		if (*input[lineNr][1] == *comp && *input[lineNr][2] == *name && *input[lineNr][3] == *param) {
			if (duplicate) {
				mess = "No more than one value of " + comp + " : " + name + " : "
					+ param + " can be set in one calculation";
				Message(-1, mess, lineNr);
			}
			duplicate = true;
			present = true;
		}
	}
	if (!present) { //no value in current system, take last value
		for (lineNr = 1; lineNr <= lastLineSystem; lineNr++) {
			if (*input[lineNr][1] == *comp && *input[lineNr][2] == *name && *input[lineNr][3] == *param) {
				present = true;
			}
		}
	}
	return present;
}
Boolean
Input::ParamSetToValue(Text comp, Text name, Text param, Text value) const {
	int lineNr;
	Text mess;
	for (lineNr = firstLineSystem; lineNr <= lastLineSystem; lineNr++) {
		if (*input[lineNr][1] == *comp && *input[lineNr][2] == *name
				&& *input[lineNr][3] == *param && *input[lineNr][4] == *value) {
			return true;
		}
	}
	return false;
}
/**
  This function returns a boolean indicating whether a line with wildcards
  corresponds to a line in a template file.
**/
Boolean
Input::ValueSetWild(Text comp, Text name, Text param) const {
	int lineNr;
	Boolean present = false;
	Text mess;
	for (lineNr = firstLineSystem; lineNr <= lastLineSystem; lineNr++) {
		if (WildCompare(input[lineNr][1],comp) &&
				WildCompare(input[lineNr][2],name) &&
				WildCompare(input[lineNr][3],param)) {
			present = true;
		}
	}
	if (!present) { //no value in current system, take last value
		for (lineNr = 1; lineNr <= lastLineSystem; lineNr++) {
			if (WildCompare(input[lineNr][1],comp) &&
					WildCompare(input[lineNr][2],name) &&
					WildCompare(input[lineNr][3],param)) {
				present = true;
			}
		}
	}
	return present;
}
int
Input::LastNumCalcValueSet(Text comp, Text name, Text param) const {
	int lineNr;
	for (lineNr = firstLineSystem; lineNr <= lastLineSystem; lineNr++) {
		if (*input[lineNr][1] == *comp && *input[lineNr][2] == *name && *input[lineNr][3] == *param) {
			return currNumCalculation;
		}
	}
	//no value in current system, search last value
	int numCalc = 1;
	int numCalculation = 0;
	for (lineNr = 1; lineNr <= lastLineSystem; lineNr++) {
		if (*input[lineNr][1] == *Copy("start")) {
			numCalc++;
		}
		if (*input[lineNr][1] == *comp && *input[lineNr][2] == *name && *input[lineNr][3] == *param) {
			numCalculation = numCalc;
		}
	}
	return numCalculation;
}
int
Input::GetInt(Text comp, Text name, Text param, int lowBound, int upBound) const {
	int value;
	Text mess;
	int lineNr = GetLineNr(comp,name,param);
	if (lineNr != -1) {
		if (!CheckInt(input[lineNr][4])) {
			mess = "Integer value needed for: " + comp + " : " + name + " : " + param + "\n"
				+ input[lineNr][4] + " is not interpreted in this way";
			Message(-1,mess, lineNr);
		}
		input[lineNr][4].Setpos(1);
		value = input[lineNr][4].Getint();
		Text number;
		number = Blanks(50);
		if (value < lowBound) {
			mess = "Illegal value for: " + comp + " : " + name + " : " + param + "\n"
				+ input[lineNr][4] + " is too low, value should be at least ";
			number.Putint(lowBound);
			number = Copy(number.Frontstrip());
			Message(-1,mess + number,lineNr);
		}
		if (value > upBound) {
			mess = "Illegal value for: " + comp + " : " + name + " : " + param + "\n"
				+ input[lineNr][4] + " is too high, value should be less than ";
			number.Putint(upBound);
			number = Copy(number.Frontstrip());
			Message(-1,mess + number, lineNr);
		}
		return value;
	}
	if (*name == *Copy("")) name = "(no name given)";
	mess = "No value found for: " + comp + " : " + name + " : " + param;
	Message(-1,mess, 0);
	return 0; // never get here
}
double
Input::GetReal(Text comp, Text name, Text param, double lowBound, double upBound) const {
	Text mess;
	double value;
	int lineNr = GetLineNr(comp,name,param);
	if (lineNr != -1) {
		if (!CheckReal(input[lineNr][4])) {
			mess = "Real value needed for: " + comp + " : " + name + " : " + param + "\n"
				+ input[lineNr][4] + " is not interpreted in this way, valid examples: 12.54 and 1.45e-4";
			Message(-1,mess, lineNr);
		}
		input[lineNr][4].Setpos(1);
		value = input[lineNr][4].Getreal();
		Text number;
		number = Blanks(9);
		if (value < lowBound) {
			mess = "Illegal value for: " + comp + " : " + name + " : " + param + "\n"
				+ input[lineNr][4] + " is too low, value should be at least ";
			number.Putreal(lowBound,3);
			number = Copy(number.Frontstrip());
			Message(-1,mess + number,lineNr);
		}
		if (value > upBound) {
			mess = "Illegal value for: " + comp + " : " + name + " : " + param + "\n"
				+ input[lineNr][4] + " is too high, value should be less than ";
			number.Putreal(upBound,3);
			number = Copy(number.Frontstrip());
			Message(-1,mess + number, lineNr);
		}
		return value;
	}
	if (*name == *Copy("")) name = "(no name given)";
	mess = "No value found for: " + comp + " : " + name + " : " + param;
	Message(-1,mess, 0);
	return 0; // never get here
}
Text
Input::GetText(Text comp, Text name, Text param) const {
	int lineNr = GetLineNr(comp,name,param);
	
	if (lineNr != -1) {
		return (input[lineNr][4]);
	}
	if (*name == *Copy("")) name = "(no name given)";
	Text mess = "No value found for: " + comp + " : " + name + " : " + param;
	Message(-1,mess, 0);
	return ""; // never get here
}
Boolean
Input::GetBoolean(Text comp, Text name, Text param) const {
	Text mess;
	int lineNr = GetLineNr(comp,name,param);
	if (lineNr != -1) {
		Text value = Copy(input[lineNr][4]);
		if (*value == *Copy("true")) return true;
		else if (*value == *Copy("false")) return false;
		else {
			mess = "Wrong value for boolean: "  + comp + " : " + name + " : " + param;
			mess = mess + "\nExpecting: 'true' or 'false'.";
			Message(-1, mess, lineNr);
		}
	}
	if (*name == *Copy("")) name = "(no name given)";
	mess = "No value found for: " + comp + " : " + name + " : " + param;
	Message(-1,mess, 0);
	return false; // never get here
}
int
Input::GetChoice(Text comp, Text name, Text param, Array<Text> choice) const {
	int i;
	Text mess;
	int numChoices = choice.Upperbound();
	if (choice.Lowerbound() != 1) {
		Sysout().Outtext("please use an array with lowerbound = 1 for call to Input::GetChoice");
		Sysout().Outimage();
		exit(-1);
	}
	int lineNr = GetLineNr(comp,name,param);
	if (lineNr != -1) {
		for (i=1; i <= numChoices; i++) {
			if (*input[lineNr][4] == *choice[i]) {
				return i;
			}
		}
		mess = "Wrong value for " + comp + " : " + name + " : " + param + " do not use: " + input[lineNr][4] + ",\n";
		if (numChoices == 1) {
			mess = mess + "only " + choice[1] + " is allowed.";
		} else {
			mess = mess + "choose from ";
			mess = mess + choice[1];
			for (i=2; i<=numChoices-1; i++) mess = mess + ", " + choice[i];
			mess = mess + " or " + choice[numChoices];
		}
		Message(-1, mess, lineNr);
	}
	if (*name == *Copy("")) name = "(no name given)";
	mess = "No value set for " + comp + " : " + name + " : " + param;
	Message(-1, mess, 0);
	return 0; // never get here
}

int
Input::GetInt(Text comp,
			  Text name,
			  Text param,
			  int lowBound,
			  int upBound,
			  int def) const {
	int value;
	Text mess;
	int lineNr = GetLineNr(comp,name,param);
	if (lineNr != -1) {
		if (!CheckInt(input[lineNr][4])) {
			mess = "Integer value needed for: " + comp + " : " + name + " : " + param + "\n"
				+ input[lineNr][4] + " is not interpreted in this way";
			Message(-1,mess, lineNr);
		}
		input[lineNr][4].Setpos(1);
		value = input[lineNr][4].Getint();
		Text number;
		number = Blanks(50);
		if (value < lowBound) {
			mess = "Illegal value for: " + comp + " : " + name + " : " + param
			+ "\n" + input[lineNr][4] + " is too low, value should be at least ";
			number.Putint(lowBound);
			number = Copy(number.Frontstrip());
			Message(-1,mess + number,lineNr);
		}
		if (value > upBound) {
			mess = "Illegal value for: " + comp + " : " + name + " : " + param
			+ "\n" + input[lineNr][4] + " is too high, value should be less than ";
			number.Putint(upBound);
			number = Copy(number.Frontstrip());
			Message(-1,mess + number, lineNr);
		}
		return value;
	}
	if (*name == *Copy("")) name = "(no name given)";
	mess = comp + " : " + name + " : " + param + " : ";
	Text number;
	number = Blanks(50);
	number.Putint(def);
	number = Copy(number.Frontstrip());
	Message(3,mess + number, 0);
	return def;
}
double
Input::GetReal(Text comp,
			   Text name,
			   Text param,
			   double lowBound,
			   double upBound,
			   double def) const {
	if (def == -1) {
		Message(0, "checking",0);
	}
	Text mess;
	double value;
	Text number;
	number = Blanks(14);
	int lineNr = GetLineNr(comp,name,param);
	if (lineNr != -1) {
		if (!CheckReal(input[lineNr][4])) {
			mess = "Real value needed for: " + comp + " : " + name + " : " + param + "\n"
				+ input[lineNr][4] + " is not interpreted in this way, valid examples: "
				"12.54 and 1.45e-4";
			Message(-1,mess, lineNr);
		}
		input[lineNr][4].Setpos(1);
		value = input[lineNr][4].Getreal();
		if (value < lowBound) {
			mess = "Illegal value for: " + comp + " : " + name + " : " + param + "\n"
			+ input[lineNr][4] + " is too low, value should be at least ";
			number.Putreal(lowBound,8);
			Message(-1,mess + number,lineNr);
		}
		if (value > upBound) {
			mess = "Illegal value for: " + comp + " : " + name + " : " + param + "\n"
			+ input[lineNr][4] + " is too high, value should be less than ";
			number.Putreal(upBound,8);
			Message(-1,mess + number, lineNr);
		}
		return value;
	}
	if (def > upBound) {
		number = Blanks(14);
		number.Putreal(def,8);
		mess = "Bug: Default value for: " + comp + " : " + name + " : " + param + "\n"
		+ number + " is too high; value should be less than ";
		number = Blanks(14);
		number.Putreal(upBound,8);
		Message(-1, mess + number, lineNr);
	}
	if (def < lowBound) {
		number = Blanks(14);
		number.Putreal(def,8);
		mess = "Bug: Default value for: " + comp + " : " + name + " : " + param + "\n"
		+ number + " is too low, value should be at least ";
		number.Putreal(lowBound,8);
		Message(-1,mess + number,lineNr);
	}
	if (*name == *Copy("")) name = "(no name given)";
	mess = comp + " : " + name + " : " + param + " : ";
	number = Blanks(14);
	number.Putreal(def,8);
	Message(3,mess + number, 0);
	return def;
}
Text
Input::GetText(Text comp, Text name, Text param, Text def) const {
	Text mess;
	int lineNr = GetLineNr(comp,name,param);
	if (lineNr != -1) {
		return input[lineNr][4];
	}
	if (*name == *Copy("")) name = "(no name given)";
	mess = comp + " : " + name + " : " + param + " : " + def;
	Message(3,mess, 0);
	return def;
}
Boolean
Input::GetBoolean(Text comp, Text name, Text param, Boolean def) const {
	Text mess;
	int lineNr = GetLineNr(comp,name,param);
	if (lineNr != -1) {
		if (*input[lineNr][4] == *Copy("true")) {
			return true;
		} else if (*input[lineNr][4] == *Copy("false")) {
			return false;
		} else {
			mess = "Wrong value for boolean: "  + comp + " : " + name + " : " + param;
			mess = mess + "\nExpecting: 'true' or 'false'.";
			Message(-1, mess, lineNr);
		}
	}
	if (*name == *Copy("")) name = "(no name given)";
	mess = comp + " : " + name + " : " + param + " : ";
	if (def) mess = mess + "true";
	else mess = mess + "false";
	Message(3,mess, 0);
	return def;
}
int
Input::GetChoice(Text comp, Text name, Text param, Array<Text> choice, int def) const {
	Text mess;
	int numChoices = choice.Upperbound();
	if (choice.Lowerbound() != 1) {
		Sysout().Outtext("Programming error: please use an array with lowerbound = 1 for call to Input::GetChoice");
		Sysout().Outimage();
		exit(-1);
	}
	int lineNr = GetLineNr(comp,name,param);
	if (lineNr != -1) {
		int i;
		for (i=1; i <= numChoices; i++) {
			if (*input[lineNr][4] == *Copy(choice[i])) {
				return i;
			}
		}
		mess = "Wrong value for " + comp + " : " + name + " : " + param + " do not use: " + input[lineNr][4] + ",\n";
		if (numChoices == 1) mess = mess + "only " + choice[1] + " is allowed.";
		else {
			mess = mess + "choose from ";
			mess = mess + choice[1];
			for (i=2; i<=numChoices-1; i++) mess = mess + ", " + choice[i];
			mess = mess + " or " + choice[numChoices];
		}
		Message(-1, mess, lineNr);
	}
	if (*name == *Copy("")) name = "(no name given)";
	mess = comp + " : " + name + " : " + param + " : " + choice[def];
	Message(3, mess, 0);
	return def;
}
int
Input::GetNumParameters(Text comp, Text name) const {
	int lineNr,count,i;
	Text mess, Number;
	Number = Blanks(9);
	int maxNr = lastLineNr;
	Text *names = new Text[maxNr];
	Boolean duplicate = false;
	count = 0;
	for (lineNr = 1; lineNr <= lastLineSystem; lineNr++) {
		if (*input[lineNr][1] == *comp && *input[lineNr][2] == *name) {
			duplicate = false;
			for (i=1; i<=count; i++) {
				if (*names[i-1] == *input[lineNr][3]) duplicate = true;
			}
			if (!duplicate) {
				count++;
				names[count-1] = input[lineNr][3];
			}
		}
	}
	delete[] names;
	return count;
}
unsigned int
Input::TotalNumCalculations() const {
	return numberOfCalculations;
}
void
Input::SetVariable(Text comp, Text name, Text param, Text value) {
	int lineNr = GetLineNr(comp,name,param);
	if (lineNr == -1) {
		if (*name == *Copy("")) name = "(no name given)";
		Text mess = "Trying to change value of " + comp + " : " + name
			+ " : " + param + "\nThis parameter is not given in the inputfile";
		Message(-1,mess, 0);
	}
	input[lineNr][4] = value;
}
void
Input::Message(int messCode, Text message, int lineNr) const {
	if (messCode <= largestMessNumToScreen) {
		if (messCode != 3) Sysout().Outimage();
		if (messCode == -1) Sysout().Outtext("Fatal error");
		if (messCode == 0) Sysout().Outtext("Fatal error");
		if (messCode == 1) Sysout().Outtext("Non fatal error");
		if (messCode == 2) Sysout().Outtext("Warning:");
		if (messCode == 3) Sysout().Outtext("Warning: default used: ");
		if (messCode != 3) {
			Sysout().Outtext(" reading inputfile '");
			Sysout().Outtext(fileName);
			if (lineNr > 0) {
				Sysout().Outtext("' in line ");
				Sysout().Outint(lineNr,0);
				Sysout().Outimage();
			} else if (lineNr == 0) {
				Sysout().Outtext("' missing line");
				Sysout().Outimage();
			} else {
				Sysout().Outtext("'");
				Sysout().Outimage();
			}
		}
		Sysout().Outtext(message);
		Sysout().Outimage();
	}
	if (messCode == -1) {
		Sysout().Outtext("Aborting...");
		Sysout().Outimage();
		exit(messCode);
	}
}
Boolean
Input::CheckReal(Text t) const {
	Boolean error = false, expset = false, pointset = false, signset = false, lastpoint = false;
	int i;
	int length = t.Length();
	char test;
	for (i=1; i<=length; i++) {
		test = t.Getchar();
		if (test == '+' || test == '-') {
			if (i != 1) error = true; // sign should be the first character
			signset = true;
			lastpoint = false;
		} else if (test == '.') {
			if (expset) error = true; // no decimal point allowed in exponent
			if (pointset) error = true; // no two decimal points allowed in mantissa
			pointset = true;
			lastpoint = true;
		} else if (isdigit(test)) {
			lastpoint = false;
		} else if (test == 'e' || test == 'E') {
			if (expset) error = true;
			if (i == 1) error = true;
			if ((signset && !pointset) || (!signset && pointset))
				if (i <= 2) error = true;
			if (signset && pointset && i <= 3) error = true;
			expset = true;
			i++;
			if (i<= length) {
				test = t.Getchar();
				if (test == '+' || test == '-' || isdigit(test)) ; // ok
				else error = true;
			} else error = true;
			lastpoint = false;
		}
		else error = true; // illegal character
		if (error) break;
	}
	if (lastpoint) error = true;
	return !error;
}
Boolean
Input::CheckInt(Text t) const {
	Boolean error = false, expset = false, signset = false;
	int i;
	int length = t.Length();
	char test;
	for (i=1; i<=length; i++) {
		test = t.Getchar();
		if (test == '+' || test == '-') {
			if (i != 1) error = true; // sign should be the first character
			signset = true;
		} else if (isdigit(test)) {
			; // ok
		} else if (test == 'e' || test == 'E') {
			if (expset) error = true;
			if (i == 1) error = true;
			if (signset && i <= 2) error = true;
			expset = true;
			i++;
			if (i<= length) {
				test = t.Getchar();
				if (test == '+' || isdigit(test)) ; // ok
				else error = true;
			} else error = true;
		}
		else error = true; // illegal character
		if (error) break;
	}
	return !error;
}
int
Input::GetLineNr(Text comp, Text name, Text param) const {
	int lineNr;
	int lineNrWithValue=0;
	Text mess;
	Boolean duplicate = false, present = false;

	for (lineNr = firstLineSystem; lineNr <= lastLineSystem; lineNr++) {
		if (*input[lineNr][1] == *comp && *input[lineNr][2] == *name && *input[lineNr][3] == *param) {
			if (duplicate) {
				mess = "No more than one value of " + comp + " : " + name + " : " + param + " can be set in one calculation";
				Message(-1, mess, lineNr);
			}
			duplicate = true;
			present = true;
			lineNrWithValue = lineNr;
		}
	}
	if (!duplicate) { //no value in current system, take last value
		for (lineNr = 1; lineNr <= lastLineSystem; lineNr++) {
			if (*input[lineNr][1] == *comp && *input[lineNr][2] == *name && *input[lineNr][3] == *param) {
				present = true;
				lineNrWithValue = lineNr;
			}
		}
	}
	if (!present) {
		return -1;
	}
	return lineNrWithValue;
}
Boolean
Input::WildCompare(const Text wildOrig, const Text fullOrig) const {
	if (*Copy(wildOrig) == *Copy(fullOrig)) {
		return true;
	}
	Text wild = Copy(wildOrig);
	Text full = Copy(fullOrig);
	wild.Setpos(1);
	full.Setpos(1);
	char test = '*';
	if (!wild.More() && full.More()) {
		return false;
	}
	if (!wild.More() && !full.More()) {
		return true;
	}
	char t = wild.Getchar();
	if (t == test) {
		if (!wild.More()) {
			return true;
		}
		while (t == test && wild.More()) {
			t = wild.Getchar();
		}
		int pos = wild.Pos();
		while (full.More()) {
			wild.Setpos(pos);
			full.Scanto(t);
			while (full.More() && wild.More()) {
				char s = wild.Getchar();
				if (s == full.Getchar()) {
					continue;
				}
				if (s == test) {
					if (wild.More()) {
						Message(-1,"illegal use of "
 							"wild card in'" +
							wildOrig + "'" ,-1);
 			  	 	 	return false;
	 				} else {
						return true;
	   				}
				}
			}
		}
		if (wild.More()) {
			return false;
		} else {
			return true;
		}
	}
	if (!full.More()) {
		return false;
	}
 	while (full.More()) {
		if (t != full.Getchar()) {
			return false;
		}
		if (!wild.More()) {
			if (!full.More()) {
				return true;
			} else {
				return false;
			}
		}
		if (!full.More()) {
			return false;
		}
		t = wild.Getchar();
		if (t == test) {
			if (wild.More()) {
				Message(-1,"illegal use of "
 					"wild card in'" +
					wildOrig + "'" ,-1);
			} else {
				return true;
			}
		}
	}
	if (wild.More()) {
		return false;
	}
	return true;
}




