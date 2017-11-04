#include "Vary.h"

Vary::Vary(Input* MyInput_) { 
	MyInput = MyInput_;
	classInstance = "";
	outputName = "";
	int numClasses = MyInput->GetNumNames("var",0,1);
	if (numClasses == 0) {
		numCalculations = 1;
		currNumCalculation = 0;
		endValue = 0;
		step = 1;
		currValue = 0;
		scale = 0;
		type = 1;
	} else {
		Array<Text> parameters = MyInput->GetNames("var");
		className = parameters[1];
		parameters = MyInput->GetParameters("var",className);
		step = MyInput->GetReal("var",className,"step",-DBL_MAX,DBL_MAX);
		Array<Text> choice(1,2);
		choice[1] = "linear";
		choice[2] = "exponential";
	//	choice[3] = "cubic";
		scale = MyInput->GetChoice("var",className,"scale",choice,1); // default linear
		choice.Dim(1,2);
		choice[1] = "integer";
		choice[2] = "real";
		type = MyInput->GetChoice("var",className,"type",choice);
		outputName = MyInput->GetText("var",className,"output_name","");
		for (int i=parameters.Lowerbound(); i<=parameters.Upperbound(); i++) {
			if (*parameters[i] == *Copy("type") ||
				*parameters[i] == *Copy("scale") ||
				*parameters[i] == *Copy("end_value") ||
				*parameters[i] == *Copy("output_name") ||
				*parameters[i] == *Copy("step")) {
					continue;
			}
			if (*classInstance == *Copy("")) {
				classInstance = Copy(parameters[i]);
				parameter = MyInput->GetText("var",className,classInstance);
			} else {
				Message(fatal,"Error in 'var : '\nTwo possibilities:"
				"1. You made a typo in one of the parameters.\n"
				"   Choose from 'type', 'scale', 'end_value', 'output_name', or 'step'\n"
				"2. You try to vary more than one variable per input file");
			}
		}
		if (!MyInput->ValueSet(className,classInstance,parameter)) {
			Message(fatal,"Error in 'var : '\nYou can only vary a variable that is set in the inputfile");
		}
		if (MyInput->LastNumCalcValueSet("var",className,"end_value") <
				MyInput->LastNumCalcValueSet(className,classInstance,parameter) ) { // bail out
			numCalculations = 1;
			currNumCalculation = 0;
			endValue = 0;
			step = 1;
			currValue = 0;
			scale = 0;
			type = 1;
			return;
		}
		startValue = MyInput->GetReal(className,classInstance,parameter,-DBL_MAX,DBL_MAX);
		if (type == 1) { // integer
			endValue = MyInput->GetInt("var",className,"end_value",-INT_MAX,INT_MAX);
		} else { // real
			endValue = MyInput->GetReal("var",className,"end_value",-DBL_MAX,DBL_MAX);
		}
		currValue = startValue;
		if (scale == 1) { //linear
			numCalculations = int((endValue - startValue)/step) + 1;
			if (int((endValue - startValue)/step) < (endValue - startValue)/step) {
				numCalculations++;
			}
			if (numCalculations <= 1) {
				Message(fatal,"Error in 'var : '\nInvalid endvalue/step combination"
					". Maybe step should be of opposite sign.");
			}
			if (endValue - startValue < 0) {
				if (step > 0) step = -step;
			}
			if (endValue - startValue > 0) {
				if (step < 0) step = -step;
			}
		}
		if (scale == 2) { //exponential
			numCalculations = int(fabs((log10(fabs(endValue)) - log10(fabs(startValue)))*step)) + 1;
			if (int(fabs((log10(fabs(endValue)) - log10(fabs(startValue)))*step)) <
					fabs((log10(fabs(endValue)) - log10(fabs(startValue)))*step)) {
				numCalculations++;
			}
			if (fabs(log10(fabs(endValue)) - log10(fabs(startValue))) <= 0) {
				Message(fatal,"Error in 'var : '\nInvalid endvalue/step combination");
			}
			if (log(endValue) - log(startValue) < 0) {
				if (step > 0) step = -step;
			}
			if (log(endValue) - log(startValue) > 0) {
				if (step < 0) step = -step;
			}
		}
		if (numCalculations > maxNumCalculations) {
			Text number = Blanks(100);
			number.Putint(maxNumCalculations);
			number = Copy(number.Strip().Frontstrip());
			Message(fatal,"Maximum number of calculations (more than " + number + ") in 'var : ' exceeded");
		}
	}
	currNumCalculation = 0;
}
Vary::~Vary() {
	if (numCalculations > 1) {
		Text value = Blanks(100);
		if (type == 1) {
			value.Putint(int(startValue));
			value = Copy(value.Frontstrip());
			MyInput->SetVariable(className,classInstance,parameter,value);
		}
		if (type == 2) {
			value.Putreal(startValue,15);
			value = Copy(value.Frontstrip());
			MyInput->SetVariable(className,classInstance,parameter,value);
		}
	}
}
Boolean
Vary::NextCalculation() {
	if (currNumCalculation == 0) {
		currNumCalculation++;
		return true;
	}
	if (type == 1) {
		if (int(currValue) - int(endValue) >= 0 && step > 0) {
			return false;
		}
		if (int(currValue) - int(endValue) <= 0 && step < 0) {
			return false;
		}
	}
	if (scale == 1) { //linear
		if (type == 2) {
			if (fabs(currValue - endValue) < fabs(step/maxNumCalculations)) {
				return false;
			}
		}
		double oldValue = currValue;
		if (type == 1) {
			while (int(currValue) - int(oldValue) == 0) {
				currValue += step;
			}
		}
		if (type == 2) {
			currValue += step;
		}
		if (currValue > endValue && step > 0) {
			currValue = endValue;
		} else if (currValue < endValue && step < 0) {
			currValue = endValue;
		}
		Text value = Blanks(100);
		if (type == 1) {
			value.Putint(int(currValue));
			value = Copy(value.Frontstrip());
			MyInput->SetVariable(className,classInstance,parameter,value);
		}
		if (type == 2) {
			value.Putreal(currValue,15);
			value = Copy(value.Frontstrip());
			MyInput->SetVariable(className,classInstance,parameter,value);
		}
	}
	if (scale == 2) { //exponential
		if (type == 2) {
			if (fabs(currValue - endValue) < fabs(currValue/maxNumCalculations)) {
				return false;
			}
		}
		double oldValue = currValue;
		if (type == 1) {
			while (int(currValue) - int(oldValue) == 0) {
				currValue *= exp(log((double)10)/step);
			}
		}
		if (type == 2) {
			currValue *= exp(log((double)10)/step);
		}
		if (currValue > endValue && step > 0) {
			currValue = endValue;
		} else if (currValue < endValue && step < 0) {
			currValue = endValue;
		}
		Text value = Blanks(100);
		if (type == 1) {
			value.Putint(int(currValue));
			value = Copy(value.Frontstrip());
			MyInput->SetVariable(className,classInstance,parameter,value);
		}
		if (type == 2) {
			value.Putreal(currValue,15);
			value = Copy(value.Frontstrip());
			MyInput->SetVariable(className,classInstance,parameter,value);
		}
	}
	if (scale == 0) { // no variable varied
		currValue += step;
	}
	currNumCalculation++;
	return true;
}
void
Vary::SkipToLast() {
	if (scale == 1) { //linear
		currValue = endValue-step;
	}
	if (scale == 2) { //exponential
		currValue = endValue/exp(log((double)10)/step);
	}
}
Text
Vary::GetVariableName() {
	if (*outputName == *Copy("")) {
		return parameter;
	}
	return outputName;
}
Text
Vary::GetCurrValue() {
	// determine accuracy needed
	char *format, *num;
	unsigned int width=0;
	int expo=0;
	double st;
	Text number;
	if (scale == 0) {
		return Copy(""); 
	}
	if (scale == 1) { //linear
		if (type == 1) { // integers: no expoonent
			width=(int)log10((fabs(startValue)>fabs(endValue))?fabs(startValue):fabs(endValue))+1;
			width=(width>100)?100:width;
			expo=0;
			format = new char[6];
			num = new char[30];
			format[0]='%';
			format[1]='0';
			width=sprintf(format+2,"%u",width);
			format[2+width]='d';
			format[3+width]='\0';
			width=sprintf(num,format,(int)currValue);
			delete [] format;
			number=Copy(num);
			delete num;
		} else {		
			st=step;
			expo=(int)log10((fabs(startValue)>fabs(endValue))?fabs(startValue):fabs(endValue));
			width=expo+1;
			/*
			double sV=startValue;
			double eV=endValue;
			while (sV != 0 || eV != 0 || st != 0) {
				width--;
				sV-=(int)(sV/pow(10,width))*pow(10,width);
				eV-=(int)(eV/pow(10,width))*pow(10,width);
				st-=(int)(st/pow(10,width))*pow(10,width);
			}
			width=expo-width; //number of significant numbers
			*/
			width=5;
			st=currValue/pow((double)10,expo);
			format = new char[6];
			num = new char[30];
			format[0]='%';
			format[1]='.';
			width=sprintf(format+2,"%u",width);
			format[2+width]='f';
			format[3+width]='\0';
			width=sprintf(num,format,st);
			num[width]='e';
			width++;
			width+=sprintf(num+width,"%u",expo);
			number=Copy(num);
			delete [] format;
			delete [] num;
		}
	} else if (scale == 2) { // exponential
		number=Blanks(100);
		width = (int)log10(
			fabs(currValue/(currValue*exp(log((double)10)/step)	- currValue))) + 1;
		if (width == 1) width++;
		number.Putreal(currValue,width);
		number=number.Frontstrip();
	}
	return number;
}


