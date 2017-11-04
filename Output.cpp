#include <iostream>
#ifndef _WIN32
#include <libgen.h>
#else
#include <stdio.h>
#endif
#include "Output.h"

Output::Output(Input* MyInput_,
			   Vary* MyVary_,
			   Array<int> grads_,
			   int accuracy_,
			   Array<Text> fileNamesOld_) {
	MyInput = MyInput_;
	MyVary = MyVary_;
	grads = grads_;
	counter=0;
	gradients = grads[0];
	accuracy = accuracyDefault = accuracy_;
	fileNamesOld = fileNamesOld_;
	errorOccurred = false;
	MCS=0;

	int numFiles = MyInput->GetNumNames("output");
	fileNames = MyInput->GetNames("output");
	Array<Text> param(1,10);
	param[1] = "accuracy";
	param[2] = "write_bounds";
	param[3] = "write";
	param[4] = "append";
	param[5] = "type";
	param[6] = "template";
	param[7] = "write_profiles";
	param[8] = "skip_errors";
	param[9] = "MC_output_interval";
	param[10] = "write_blank_lines";
	Array<Text> types(1,7);
	types[1] = "ana";
	types[2] = "kal";
	types[3] = "vtk";
	types[4] = "ave";
	types[5] = "template";
	types[6] = "profiles";
	types[7] = "vectors";
	if (numFiles == 0) {
		Message(warning,MyInput,"no output defined");
	}
	int i;
	for (i=1; i<=numFiles; i++) {
		Text name = fileNames[i];
		MyInput->CheckParameterNames("output",name,param);
		accuracy = MyInput->GetInt("output",name,"accuracy",1,100,accuracyDefault);
		writeBounds = MyInput->GetBoolean("output",name,"write_bounds",true);
		write = MyInput->GetBoolean("output",name,"write",true);
		append = MyInput->GetBoolean("output",name,"append",false);
		outputProfiles = MyInput->GetBoolean("output",name,"write_profiles",true);
		blank_line = MyInput->GetBoolean("output",name,"write_blank_lines",true);
		skipErrors = MyInput->GetBoolean("output",name,"skip_errors",false);
		MyInput->GetChoice("output",name,"type",types);
		templateSet = MyInput->ValueSet("output",name,"template");
		if (templateSet) {
			// check existence and syntax of template file
			Template = new Input(MyInput->GetText("output",name,"template"));
			delete Template;
		}
	}
	if (fileNames.Upperbound() > fileNamesOld.Upperbound()) {
		firstWriteQ.Dim(1,fileNames.Upperbound());
		for (i=1; i<=fileNames.Upperbound(); i++) {
			firstWriteQ[i] = true;
			for (int j=1; j<=fileNamesOld.Upperbound(); j++) {
				if (*Copy(fileNames[i]) == *Copy(fileNamesOld[j])) {
					firstWriteQ[i] = false;
				}
			}
		}
	} else {
		firstWriteQ.Dim(1,fileNames.Upperbound());
		for (i=1; i<=fileNames.Upperbound(); i++) {
			firstWriteQ[i] = false;
		}
	}
	for (i=1; i<=numFiles; i++) {
		Text name = fileNames[i];
		append = MyInput->GetBoolean("output",name,"append",false);
		int choice = MyInput->GetChoice("output",name,"type",types);
		if (choice == 4) {
			name = getNumberedFilename(name);
		}
		if (!append && firstWriteQ[i]) {
			if (FileExists(name)) {
				Message(warning,"An output file named '" + parseFileName(name,MyInput)
					+ "' already "
					"exists, it will be overwritten after the calculation.");
			}
		}
	}
}
Output::~Output() {
	while (OutputLineQ.Cardinal() > 0) {
		delete (OutputLine*) OutputLineQ[1];
	}
}
void
Output::PutText(const Text elem1,
				const Text elem2,
				const Text elem3,
				const Text value) {
	(new OutputLine(elem1,elem2,elem3,value))->Into(OutputLineQ);
}
void
Output::PutBoolean(const Text elem1,
				   const Text elem2,
				   const Text elem3,
				   const bool value) {
	Text elem4;
	if (value) elem4 = "true";
	else elem4 = "false";
	(new OutputLine(elem1,elem2,elem3,elem4))->Into(OutputLineQ);
}
void
Output::PutInt(const Text elem1,
			   const Text elem2,
			   const Text elem3,
			   const int value) {
	Text elem4;
	elem4 = Blanks(100);
	elem4.Putint(value);
	elem4 = Copy(elem4.Strip().Frontstrip());
	(new OutputLine(elem1,elem2,elem3,elem4))->Into(OutputLineQ);
}
void
Output::PutReal(const Text elem1,
				const Text elem2,
				const Text elem3,
				const double value) {
	Text elem4;
	elem4 = Blanks(100);
	elem4.Putreal(value,accuracy);
	elem4 = Copy(elem4.Strip().Frontstrip());
	(new OutputLine(elem1,elem2,elem3,elem4))->Into(OutputLineQ);
}
void
Output::PutProfile(const Text elem1,
				   const Text elem2,
				   const Text elem3,
				   const Vector value) {
	(new OutputLine(elem1,elem2,elem3,"profile",value))->Into(OutputLineQ);
}
void
Output::PutVector(const Text elem1,
				  const Text elem2,
				  const Text elem3,
				  const Vector value,
				  const int max,
				  const int min) {
	(new OutputLine(elem1,elem2,elem3,"vector",value, max,min))->Into(OutputLineQ);
}

void
Output::WriteMCOutput(int mcs) {
	Array<Text> types(1,7);
	types[1] = "ana";
	types[2] = "kal";
	types[3] = "vtk";
	types[4] = "ave";
	types[5] = "template";
	types[6] = "profiles";
	types[7] = "vectors";
	int numFiles = MyInput->GetNumNames("output");
	MCS = mcs;
	int i;
	for (i=1; i<=numFiles; i++) {
		Text name = fileNames[i];
		accuracy = MyInput->GetInt("output",name,"accuracy",1,100,accuracyDefault);
		writeBounds = MyInput->GetBoolean("output",name,"write_bounds",true);
		write = MyInput->GetBoolean("output",name,"write",true);
		//append = MyInput->GetBoolean("output",name,"append",false);
		append = true;
		outputProfiles = MyInput->GetBoolean("output",name,"write_profiles",true);
		MC_output_interval = MyInput->GetInt("output",name,"MC_output_interval",1,10000,1000);
		skipErrors = MyInput->GetBoolean("output",name,"skip_errors",false);
		if (skipErrors && errorOccurred) {
			write = false;
		}
		int choice = MyInput->GetChoice("output",name,"type",types);
		templateSet = MyInput->ValueSet("output",name,"template");
		if (templateSet) {
			Template = new Input(MyInput->GetText("output",name,"template"));
			while (Template->NextCalculation()); // make Template ready for read
		}
		switch(choice) {
			case 1:
				append = MyInput->GetBoolean("output",name,"append",true);
				OutputAna(name,firstWriteQ[i]);
				break;
			case 2:
				append = MyInput->GetBoolean("output",name,"append",false);
				OutputKal(name,false);
				break;
			case 3:
				append = MyInput->GetBoolean("output",name,"append",false);
				OutputVtkProfiles(name);
				break;
			case 4:
				append = false;
				OutputAve(name,MCS==1);
				break;
			case 5:
				OutputTemplate(name,false);
				break;
			case 6:
				append = MyInput->GetBoolean("output",name,"append",false);
				OutputProfiles(name);
				break;
			case 7:
				append = MyInput->GetBoolean("output",name,"append",false);
				OutputVectors(name);
				break;
			default:
				//error
				break;
		};
		if (templateSet) {
			delete Template;
		}
	}
	fileNamesOld.Dim(1,numFiles);
	for (i=1; i<=numFiles; i++) {
		fileNamesOld[i] = Copy(fileNames[i]);
	}
}


void
Output::WriteOutput() {
	Array<Text> types(1,7);
	types[1] = "ana";
	types[2] = "kal";
	types[3] = "vtk";
	types[4] = "ave";
	types[5] = "template";
	types[6] = "profiles";
	types[7] = "vectors";
	int numFiles = MyInput->GetNumNames("output");
	int i;
	MCS = 0;
	for (i=1; i<=numFiles; i++) {
		Text name = fileNames[i];
		accuracy = MyInput->GetInt("output",name,"accuracy",1,100,accuracyDefault);
		writeBounds = MyInput->GetBoolean("output",name,"write_bounds",true);
		write = MyInput->GetBoolean("output",name,"write",true);
		blank_line = MyInput->GetBoolean("output",name,"write_blank_lines",true);
		append = MyInput->GetBoolean("output",name,"append",false);
		outputProfiles = MyInput->GetBoolean("output",name,"write_profiles",true);
		skipErrors = MyInput->GetBoolean("output",name,"skip_errors",false);
		if (skipErrors && errorOccurred) {
			write = false;
		}
		int choice = MyInput->GetChoice("output",name,"type",types);
		templateSet = MyInput->ValueSet("output",name,"template");
		if (templateSet) {
			Template = new Input(MyInput->GetText("output",name,"template"));
			while (Template->NextCalculation()); // make Template ready for read
		}
		switch(choice) {
			case 1:
				append = MyInput->GetBoolean("output",name,"append",true);
				OutputAna(name,firstWriteQ[i]);
				break;
			case 2:
				append = MyInput->GetBoolean("output",name,"append",false);
				OutputKal(name,firstWriteQ[i]);
				break;
			case 3:
				append = MyInput->GetBoolean("output",name,"append",false);
				OutputVtkProfiles(name);
				break;
			case 4:
				append = false;
				OutputAve(name,firstWriteQ[i]);
				break;
			case 5:
				OutputTemplate(name,firstWriteQ[i]);
				break;
			case 6:
				append = MyInput->GetBoolean("output",name,"append",false);
				OutputProfiles(name);
				break;
			case 7:
				append = MyInput->GetBoolean("output",name,"append",false);
				OutputVectors(name);
				break;
			default:
				//error
				break;
		};
		if (templateSet) {
			delete Template;
		}
	}
	fileNamesOld.Dim(1,numFiles);
	for (i=1; i<=numFiles; i++) {
		fileNamesOld[i] = Copy(fileNames[i]);
	}
}
void
Output::Clear() {
	while(OutputLineQ.Cardinal() > 0) {
		delete (OutputLine*) OutputLineQ[1];
	}
}
Array<Text>
Output::GetFileNames() {
	return fileNames;
}
void
Output::SetErrorOccurred() {
	errorOccurred = true;
}
void
Output::OutputAna(Text name, bool firstWrite) {
	ofstream Out;
	open(Out, name, firstWrite);
	if (!append && firstWrite) {
		Out << "system delimiter" << endl << endl;
	}
	OutputLine* Line;
	for (int i=1; i<=OutputLineQ.Cardinal(); i++) {
		bool write2 = true;
		Line = (OutputLine*) OutputLineQ[i];
		if (templateSet) {
			write2 = ValueInTemplate(Line);
		}
		if (write && write2) {
			if (!Line->IsProfile() || outputProfiles) {
				Out << Line->GetElement1().MainC() << " : ";
				Out << Line->GetElement2().MainC() << " : ";
				Out << Line->GetElement3().MainC() << " : ";
				Out << Line->GetValue().MainC() << endl;
			}
			if (Line->IsProfile() && outputProfiles) {
				Vector profile = Line->GetProfile();
				if (gradients == 1) {
					for (int z=1; z<=grads[1]; z++) {
						if (writeBounds || (z>1 && z<grads[1])) {
							Out << profile[z] << endl;
						}
					}
				} else if (gradients == 2) {
					int z=0;
					for (int x=1; x<=grads[1]; x++) {
						for (int y=1; y<=grads[2]; y++) {
							z++;
							if (writeBounds || (x>1 && x<grads[1] && y>1 && y<grads[2])) {
								Out << profile[z] << endl;
							}
						}
					}
				}
				else if (gradients == 3) {
					int N_comp_ranges=grads[4];
					if (N_comp_ranges==1){
						int i=0;
						for (int x=0; x<grads[1]; x++) {
							for (int y=0; y<grads[2]; y++) {
								for (int z=0; z<grads[3]; z++) {
									i++;
									if ((x>1 && x<grads[1]-1 && y>1 && y<grads[2]-1 && z>1 && z<grads[3]-1)) {
										Out << profile[i] << endl;
									}
								}
							}
						}
					}
					else {
						int i,x,y,z,jx,jy,x0,y0,z0,xm,ym,zm,Current=0,Next=0;
						for (i=1; i<=N_comp_ranges; i++){
							x0=grads[4+(i-1)*6+1]; xm=grads[4+(i-1)*6+4];
							y0=grads[4+(i-1)*6+2]; ym=grads[4+(i-1)*6+5];
							z0=grads[4+(i-1)*6+3]; zm=grads[4+(i-1)*6+6];

							jy=(zm-z0+3);
							jx=jy*(ym-y0+3);
							Current=Next; Next+=jx*(xm-x0+3);
							for (x=x0; x<=xm; x++) {
								for (y=y0; y<=ym; y++){
									for (z=z0; z<=zm; z++){
										//if (Current+x*jx+y*jy+z+1>350)
										//cout<< "("<< x << "," << y << "," << z << ") current" << Current << "jx" << jx << "jy"<< jy << endl;
										Out << "("<< x << "," << y << "," << z << ") "<<'\t'<< profile[Current+(x-x0+1)*jx+(y-y0+1)*jy+(z-z0+1)+1] << endl;
									}
								}
							}
						}
					}
				}
			} else if (Line->IsVector()) {
				Vector value = Line->GetVector();
				int min = Line->GetMin();
				int max = Line->GetMax();
				for (int z=min; z<=max; z++) {
					Out << value[z] << endl;
				}
			}
		}
	}
	if (write) {
		Out << "system delimiter" << endl << endl;
	}
#ifdef _WIN32
	Text temp = Copy("attrib -r \"" + name+"\"");
	system(temp.MainC());
#endif
}

void
Output::OutputAve(Text name, bool firstWrite) {
	double* Values=0;
	double* Val_err=0;
	double* Fluct=0;
	double* Fluct_err=0;
	double Value_old=0, Val_err_old=0, Fluct_old=0, Fluct_err_old=0;
	int* Numbers=0;
	int number=0;
	int lineNr=0;
	int elementNr=0;
	bool more;
	int num_of_variables=0;
	Text lineContents;
	Text TextOfValue;
	if (!write) return;
	if (!firstWrite) {
		name = parseFileName(name,MyInput);
		Infile file(name);
		file.Open(Blanks(imageLength));
		more = file.Inrecord();

		while (lineNr<6) {
			lineNr++;
			lineContents = notext;
			while ( more ) {
				lineContents = lineContents + file.Image;
				more = file.Inrecord();
			}
			lineContents =lineContents + file.Image.Sub(1,file.Image.Pos()-1).Strip();
			Text dummy = "";
			Boolean endLine = false;
			while(lineContents.More()) {
				char test2 = lineContents.Getchar();
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
			lineContents = Copy(dummy);
			elementNr = 0;
			if (lineNr==1) {
				while (lineContents.More()) {
					elementNr++;
					lineContents.Scanto('\t');
				}
				Numbers = new int[elementNr];
				Values= new double[elementNr];
				Val_err = new double[elementNr];
				Fluct = new double[elementNr];
				Fluct_err = new double [elementNr];
				num_of_variables=elementNr-1;

			} else {
				Text input;
				lineContents=Copy(lineContents.Frontstrip().Strip());
				lineContents.Setpos(1);
				input = Copy(lineContents.Scanto('\t')); //skip first field
				while (lineContents.More()) {
					elementNr++;
					input = Copy(lineContents.Scanto('\t'));
					input=Copy(input.Frontstrip().Strip());	input.Setpos(1);
					if (lineNr ==2) {
						if (CheckReal(input)) {
							input.Setpos(1);
							Values[elementNr] = input.Getreal();
						}
					}
					if (lineNr ==3) Val_err[elementNr] = input.Getreal();
					if (lineNr ==4) Fluct[elementNr] = input.Getreal();
					if (lineNr ==5) Fluct_err[elementNr] = input.Getreal();
					if (lineNr ==6) Numbers[elementNr]= input.Getint();
				}
			}
			more = file.Inrecord();
		}
		file.Close();
	}  else number=1;

	firstWrite=true;
	ofstream Out;
	open(Out, name, firstWrite);
	OutputLine* Line;
	int k=0;
	bool firstValueWritten = false;
	if ((firstWrite && !append) || !FileExists(name)) {
		for (int i=1; i<=OutputLineQ.Cardinal(); i++) {
			bool write2 = true;
			Line = (OutputLine*) OutputLineQ[i];
			if (templateSet && write2) {
				write2 = ValueInTemplate(Line);
			}
			if (write2) {
				k++;
				//if (firstValueWritten) {
					Out << '\t'; //start with tab.
				//}
				Out << Line->GetElement1().MainC() << " : ";
				Out << Line->GetElement2().MainC() << " : ";
				Out << Line->GetElement3().MainC() ;
				//firstValueWritten = true;
			}
		}
		Out << endl;
	}
	if (num_of_variables == 0) {
		num_of_variables = k;
		Numbers = new int[k+1];
		Values= new double[k+1];
		Val_err = new double[k+1];
		Fluct = new double[k+1];
		Fluct_err = new double [k+1];
	}
	k=0;
	firstValueWritten = false;
	for (int i=1; i<=OutputLineQ.Cardinal(); i++) {
		bool write2 = true;
		Line = (OutputLine*) OutputLineQ[i];
		if (templateSet && write2) {
			write2 = ValueInTemplate(Line);
		}
		if (write2) {

			if (!firstValueWritten) {
				Out << "Values" << '\t';
				firstValueWritten=true;
			} else Out << '\t';

			TextOfValue = Copy(Line->GetValue().MainC());
			double Value_n;
			if (CheckReal(TextOfValue)) {
				TextOfValue.Setpos(1);
				Value_n=TextOfValue.Getreal();
			} else {Value_n=123456789;}
			k++;

			if (number==1) {

				if (Value_n < 1000000) {
					Values[k]=Value_n; Numbers[k]=number;
				} else {
					Values[k]=0; Numbers[k]=0;
				}
				Val_err[k]=0;
				Fluct[k]=0;
				Fluct_err[k]=0;
			} else {
				double Fluct_err_n,Value2,Val_err_n;
				int Num_old;
				Value_old= Values[k];
				Val_err_old=Val_err[k];
				Fluct_old=Fluct[k];
				Fluct_err_old=Fluct_err[k];
				Num_old=Numbers[k];
				if (Value_n > 1000000) {
					//keep old values;
				} else {
					if (Num_old==0) {
						Values[k]=Value_n;
						Val_err[k]=0;
						Fluct[k]=0;
						Fluct_err[k]=0;
					} else {
						Values[k]=(Value_old*Num_old+ Value_n)/(Numbers[k]+1); 				 //cumulative moving average value;
						//if (Num_old>1) {
							Value2 = Fluct_old*Fluct_old+Value_old*Value_old;				 //Find old value of moving average of value*value
							Value2 = (Value2*(Num_old) + Value_n*Value_n)/(Numbers[k]+1);		 //cumulative moving average of value*value;
							if ((Value2-Values[k]*Values[k])>0) {
								Fluct[k]=pow(Value2-Values[k]*Values[k],0.5);			     //definition of fluctuations based on moving average (value*value) and (Value[k])
							}
							Val_err_n=pow(pow((Values[k]-Value_old),2),0.5);			      //estimate error in average from difference with previous result
							Val_err[k]=(Val_err_old*(Num_old)+Val_err_n)/(Numbers[k]+1);           //cumulative moving average of error in value;
							if ((Fluct[k]-Fluct_old)>1e-10) {
								Fluct_err_n=pow(pow((Fluct[k]-Fluct_old),2),0.5);		      //estimate error from difference of fluct[k] and the previous one.
								Fluct_err[k]=(Fluct_err_old*(Num_old)+Fluct_err_n)/(Numbers[k]+1); //cumulative moving average of error in fluct;
							}
						//}
					}
					Numbers[k]++;
				}
			}
			if (Value_n==123456789) {
				Out << Line->GetValue().MainC();
			} else {Out << Values[k];}
		}
	}
	Out << endl;
	for (lineNr=3; lineNr<7; lineNr++) {
		for (k=1; k<=num_of_variables; k++) {
			if (k==1) {
				//if (lineNr==2) Out << "Values";
				if (lineNr==3) Out << "Val_err";
				if (lineNr==4) Out << "Fluct";
				if (lineNr==5) Out << "Fluct_err";
				if (lineNr==6) Out << "Numbers";
			}
			Out << '\t';
			if (lineNr==2) Out << Values[k];
			if (lineNr==3) Out << Val_err[k];
			if (lineNr==4) Out << Fluct[k];
			if (lineNr==5) Out << Fluct_err[k];
			if (lineNr==6) Out << Numbers[k];
		}
		Out << endl;
	}

#ifdef _WIN32
	Text temp = Copy("attrib -r \"" + name+"\"");
	system(temp.MainC());
#endif
}

void
Output::OutputKal(Text name, bool firstWrite) {
	ofstream Out;
	open(Out, name, firstWrite);
	OutputLine* Line;

	bool firstValueWritten = false;
	if ((firstWrite && !append) || !FileExists(name)) {
		for (int i=1; i<=OutputLineQ.Cardinal(); i++) {
			bool write2 = true;
			Line = (OutputLine*) OutputLineQ[i];
			if (Line->IsProfile() || Line->IsVector()) {
				write2=false;
				if (templateSet && gradients == 1 && ValueInTemplate(Line)) {
					for (int z=1; z<=grads[1]; z++) {
						if (ValueForLayerInTemplate(Line,z-1)) {
							if (firstValueWritten) {
								Out << '\t';
							}
							Out << Line->GetElement1().MainC() << " : ";
							Out << Line->GetElement2().MainC() << " : ";
							Out << Line->GetElement3().MainC() << " : ";
							Out << z-1;
							firstValueWritten = true;
						}
					}
				}
			}
			Line = (OutputLine*) OutputLineQ[i];
			if (templateSet && write2) {
				write2 = ValueInTemplate(Line);
			}
			if (write2) {
				if (firstValueWritten) {
					Out << '\t';
				}
				Out << Line->GetElement1().MainC() << " : ";
				Out << Line->GetElement2().MainC() << " : ";
				Out << Line->GetElement3().MainC();
				firstValueWritten = true;
			}
		}
		Out << endl;
	}
	if (!write) return;
	firstValueWritten = false;
	for (int i=1; i<=OutputLineQ.Cardinal(); i++) {
		bool write2 = true;
		Line = (OutputLine*) OutputLineQ[i];
		if (Line->IsProfile() || Line->IsVector()) {
			if (templateSet && gradients == 1 && ValueInTemplate(Line)) {
				for (int z=1; z<=grads[1]; z++) {
					write2 = false;
					if (ValueForLayerInTemplate(Line,z-1)) {
						Vector profile = Line->GetProfile();
						if (firstValueWritten) {
							Out << '\t';
						}
						Out << profile[z];
						firstValueWritten = true;
					}
				}
			} else {
				write2 = false;
			}
		}
		if (templateSet && write2) {
			write2 = ValueInTemplate(Line);
		}
		if (write2) {
			if (firstValueWritten) {
				Out << '\t';
			}
			Out << Line->GetValue().MainC();
			firstValueWritten = true;
		}
	}

	if (write) {
		Out << endl;
	}
#ifdef _WIN32
	Text temp = Copy("attrib -r \"" + name+"\"");
	system(temp.MainC());
#endif

}

void
Output::OutputVtkProfiles(Text name) {
	if (!write) {
		return;
	}

	name=getNumberedFilename(name);
	ofstream Out;
	open(Out, name);
	OutputLine* Line;
	Out << "# vtk DataFile Version 3.0" << endl;
	Out << "vtk output" << endl;
	Out << "ASCII" << endl;
	Array<Text> latname;
	latname = MyInput->GetNames("lat");

	int lattyp = 1; // default standard lattice
	if (MyInput->ValueSet("lat",latname[1],"latticetype")) {
		Array<Text> lattype(1,4);
		lattype[1] = "standard";
		lattype[2] = "stencils";
		lattype[3] = "FCC";
		lattype[4] = "HEX";
		lattyp = MyInput->GetChoice("lat", latname[1], "latticetype", lattype, 1);
	}


	if (lattyp == 4) {Out << "DATASET STRUCTURED_GRID" << endl; }
	else {Out << "DATASET STRUCTURED_POINTS" << endl;}
	if (gradients == 1) {
		int n_layers = grads[1]-2;
		Out << "DIMENSIONS" << n_layers << endl;
		Out << "SPACING 1" << endl;
		Out << "ORIGIN 0" << endl;
		Out << "POINT_DATA " << n_layers << endl;
	} else if (gradients == 2) {
		int n_layers_x = grads[1]-2;
		int n_layers_y = grads[2]-2;
		Out << "DIMENSIONS " <<  n_layers_x << " " <<  n_layers_y << " 1" <<endl;
              Out << "SPACING 1 1 1" << endl;
		Out << "ORIGIN 0 0 0" << endl;
		Out << "POINT_DATA " << n_layers_x*n_layers_y << endl;
	}  else if (gradients == 3) {
		int n_layers_x = grads[1]-2;
		int n_layers_y = grads[2]-2;
		int n_layers_z = grads[3]-2;
		Out << "DIMENSIONS " <<  n_layers_x << " " <<  n_layers_y << " " << n_layers_z << endl;
		if (lattyp == 4) {
			Out << "POINTS " << n_layers_x*n_layers_y*n_layers_z << " double" << endl;
			double x,y,z;
			double yc1 = pow(0.75,0.5);
			double yc2 = pow(0.1875,0.5);
			for (int i=0; i<n_layers_z; i++) {
				for (int j=0; j<n_layers_y; j++){
					for (int k=0; k<n_layers_x; k++){
						x = k+0.5*j+0.5*i;
            					y = yc1*j+yc2*i;
            					z = 0.75*i;
						Out << x << " " << y << " " << z << endl;
					}
				}
			}
		}

		else {
			Out << "SPACING 1 1 1" << endl;
			Out << "ORIGIN 0 0 0" << endl;
		}
		Out << "POINT_DATA " << n_layers_x*n_layers_y*n_layers_z << endl;

	}
	Out << "SCALARS " ;

	for (int i=1; i<=OutputLineQ.Cardinal(); i++) {
		bool write2 = true;
		Line = (OutputLine*) OutputLineQ[i];
		if (!Line->IsProfile()) {
			continue;
		}
		if (templateSet) {
			write2 = ValueInTemplate(Line);
		}
		if (write2) {
			//Out << Line->GetElement1().MainC() << "_";
			//Out << Line->GetElement2().MainC() << "_";
			//Out << Line->GetElement3().MainC() << " double " << endl; //in deze regel kunnen spaties zitten en dan gaat 't mis.
			Out << "SFBox_profile double" << endl;
		}
	}
	Out << "LOOKUP_TABLE default "<< endl;

	if (gradients == 1) {
		for (int z=1; z<=grads[1]; z++) {
			if ((z>1 && z<grads[1]-1)) {
				for (int i=1; i<=OutputLineQ.Cardinal(); i++) {
					Line = (OutputLine*) OutputLineQ[i];
					if (!Line->IsProfile()) {
						continue;
					}
					bool write2 = true;
					Line = (OutputLine*) OutputLineQ[i];
					if (templateSet) {
						write2 = ValueInTemplate(Line);
					}
					if (write2) {
						Vector profile = Line->GetProfile();
						Out << profile[z];
					}
				}
				Out << endl;
			}
		}
	} else if (gradients == 2) {
		int i=1;
		int n_layers_x = grads[1]-2;
		int n_layers_y = grads[2]-2;
		int jx = (n_layers_y+2);
		bool write2 = false;
		Line = (OutputLine*) OutputLineQ[1];
		while (i<=OutputLineQ.Cardinal() && !write2) {
			Line = (OutputLine*) OutputLineQ[i];
			if (templateSet) write2 = ValueInTemplate(Line);
			i++;
		}
		if (!Line->IsProfile()) Message(warning,"Template file does not contain (known) profile in Vtk output");
		Vector profile = Line->GetProfile();

		for (int y=1; y<=n_layers_y; y++) {
			for (int x=1; x<=n_layers_x; x++) {
				Out << profile[jx*x+y] << endl;
			}
		}
	} else if (gradients == 3) {
		int n_layers_x = grads[1]-2;
		int n_layers_y = grads[2]-2;
		int n_layers_z = grads[3]-2;
		int jx = (n_layers_y+2)*(n_layers_z+2);
		int jy = n_layers_z+2;
		int jz = 1;
		int i =1;
		bool write2 = false;
		Line = (OutputLine*) OutputLineQ[1];
		while (i<=OutputLineQ.Cardinal() && !write2) {
			Line = (OutputLine*) OutputLineQ[i];
			if (templateSet) write2 = ValueInTemplate(Line);
			i++;
		}
		if (!Line->IsProfile()) Message(warning,"Template file does not contain (known) profile in Vtk output");
		Vector profile = Line->GetProfile();

		int N_comp_ranges=grads[4];
		if (N_comp_ranges==1) {
			for (int z=1; z<=n_layers_z; z++) {
				for (int y=1; y<=n_layers_y; y++) {
					for (int x=1; x<=n_layers_x; x++) {Out << profile[jx*x+jy*y+z*jz+1] << endl;}
				}
			}
		}
		else {
			int i,jx,jy,x0,y0,z0,xm,ym,zm,Current=0,Next=0;
			for (i=1; i<=N_comp_ranges; i++){
				x0=grads[4+(i-1)*6+1]; xm=grads[4+(i-1)*6+4];
				y0=grads[4+(i-1)*6+2]; ym=grads[4+(i-1)*6+5];
				z0=grads[4+(i-1)*6+3]; zm=grads[4+(i-1)*6+6];

				jy=(zm-z0+3);
				jx=jy*(ym-y0+3);
				Current=Next; Next+=jx*(xm-x0+3);

				for (int z=z0; z<=zm; z++) {
					for (int y=y0; y<=ym; y++) {
						for (int x=x0; x<=xm; x++) {
							for (int k=1; k<=OutputLineQ.Cardinal(); k++) {
								Out << profile[Current+(x-x0+1)*jx+(y-y0+1)*jy+(z-z0+1)+1] << endl;
									//possibly we need to put the elements one down...
									//writebounds true ignored;

							}
						}
					}
				}
			}
		}
	}
#ifdef _WIN32
	Text temp = Copy("attrib -r \"" + name+"\"");
	system(temp.MainC());
#endif
}


void
Output::OutputVectors(Text name) {
	bool error_made=false;
	int count = 0;
	int min=0;
	int max=0;

	if (!write) {
		return;
	}
	name=getNumberedFilename(name);
	ofstream Out;
	open(Out, name);
	OutputLine* Line;

	Out << "element" << '\t';
	for (int i=1; i<=OutputLineQ.Cardinal(); i++) {
		bool write2 = true;
		Line = (OutputLine*) OutputLineQ[i];
		if (!Line->IsVector()) {

			continue;
		}
		if (templateSet) {
			write2 = ValueInTemplate(Line);
		}
		if (write2) {
			Out << Line->GetElement1().MainC() << " : ";
			Out << Line->GetElement2().MainC() << " : ";
			Out << Line->GetElement3().MainC() << '\t';
			count++;
			int min_ = Line->GetMin(); if (min==0) min=min_;
			else {
				if (min != min_) {Message(warning,"Vectors should be equally long when outputted in one vector output file"); error_made=true;}
			}
			int max_ = Line->GetMax(); if (max==0) max=max_;
			else {
				if (max != max_) {Message(warning,"Vectors should be equally long when outputted in one vector output file"); error_made=true;}
			}
		}
	}
	Out << endl;

	if (!error_made) {
		Matrix X(min,max,1,count);
		Vector profile(min,max);
		int k=0;
		for (int i=1; i<=OutputLineQ.Cardinal(); i++) {
			Line = (OutputLine*) OutputLineQ[i];
			if (!Line->IsVector()) {
				continue;
			}
			bool write2 = true;
			Line = (OutputLine*) OutputLineQ[i];
			if (templateSet) {
				write2 = ValueInTemplate(Line);
			}
			if (write2) {
				k++;
				profile = Line->GetVector();
				for (int z=min; z<=max; z++) X[z][k] = profile[z];
			}
		}
		for (int z=min; z<=max; z++) {
			Out << z << '\t';
			for (int k=1; k<=count; k++) {
				Out << X[z][k] << '\t';
			}

			Out << endl;
		}
	}

#ifdef _WIN32
	Text temp = Copy("attrib -r \"" + name+"\"");
	system(temp.MainC());
#endif
}

void
Output::OutputProfiles(Text name) {
	int count = 0;
	int M=0;

	if (!write) {
		return;
	}
	name=getNumberedFilename(name);
	ofstream Out;
	open(Out, name);
	OutputLine* Line;
	if (gradients == 1) {
		Out << "layer" << '\t';
		M=grads[1];
	} else if (gradients == 2) {
		Out << "layer x" << '\t';
		Out << "layer y" << '\t';
		M=grads[1]*grads[2];
	}  else if (gradients == 3) {
		Out << "layer x" << '\t';
		Out << "layer y" << '\t';
		Out << "layer z" << '\t';
		M=grads[1]*grads[2]*grads[3];
	}

	for (int i=1; i<=OutputLineQ.Cardinal(); i++) {
		bool write2 = true;
		Line = (OutputLine*) OutputLineQ[i];
		if (!Line->IsProfile()) {
			continue;
		}
		if (templateSet) {
			write2 = ValueInTemplate(Line);
		}
		if (write2) {
			Out << Line->GetElement1().MainC() << " : ";
			Out << Line->GetElement2().MainC() << " : ";
			Out << Line->GetElement3().MainC() << '\t';
			count++;
		}
	}
	Out << endl;

	if (gradients == 1) {
		Matrix X(1,M,1,count);
		Vector profile(1,M);
		int k=0;
		for (int i=1; i<=OutputLineQ.Cardinal(); i++) {
			Line = (OutputLine*) OutputLineQ[i];
			if (!Line->IsProfile()) {
				continue;
			}
			bool write2 = true;
			Line = (OutputLine*) OutputLineQ[i];
			if (templateSet) {
				write2 = ValueInTemplate(Line);
			}
			if (write2) {
				k++;
				profile = Line->GetProfile();
				for (int z=1; z<=M; z++) X[z][k] = profile[z];
			}
		}
		for (int z=1; z<=grads[1]; z++) {
			if (writeBounds || (z>1 && z<grads[1])) {
				Out << z-1 << '\t';
				for (int k=1; k<=count; k++) {
					Out << X[z][k] << '\t';
				}
				//for (int i=1; i<=OutputLineQ.Cardinal(); i++) {
				//	Line = (OutputLine*) OutputLineQ[i];
				//	if (!Line->IsProfile()) {
				//		continue;
				//	}
				//	bool write2 = true;
				//	Line = (OutputLine*) OutputLineQ[i];
				//	if (templateSet) {
				//		write2 = ValueInTemplate(Line);
				//	}
				//	if (write2) {
				//		profile = Line->GetProfile();
				//		Out << profile[z] << '\t';
				//	}
				//}
				Out << endl;
			}
		}
	} else if (gradients == 2) {
		int z=0;
		Matrix X(1,M,1,count);
		Vector profile(1,M);
		int k=0;
		for (int i=1; i<=OutputLineQ.Cardinal(); i++) {
			Line = (OutputLine*) OutputLineQ[i];
			if (!Line->IsProfile()) {
				continue;
			}
			bool write2 = true;
			Line = (OutputLine*) OutputLineQ[i];
			if (templateSet) {
				write2 = ValueInTemplate(Line);
			}
			if (write2) {
				k++;
				profile = Line->GetProfile();
				for (int z=1; z<=M; z++) X[z][k] = profile[z];
			}
		}
		for (int x=1; x<=grads[1]; x++) {
			for (int y=1; y<=grads[2]; y++) {
				z++;
				if (writeBounds || (x>1 && x<grads[1] && y>1 && y<grads[2])) {
					Out << x-1 << '\t' << y-1 << '\t';
					for (int k=1; k<=count; k++) {
						Out << X[z][k] << '\t';
					}
					//for (int i=1; i<=OutputLineQ.Cardinal(); i++) {
					//	Line = (OutputLine*) OutputLineQ[i];
					//	if (!Line->IsProfile()) {
					//		continue;
					//	}
					//	bool write2 = true;
					//	Line = (OutputLine*) OutputLineQ[i];
					//	if (templateSet) {
					//		write2 = ValueInTemplate(Line);
					//	}
					//	if (write2) {
					//		Vector profile = Line->GetProfile();
					//		Out << profile[z] << '\t';
					//	}
					//}
					Out << endl;
				}
			}
			if (blank_line) Out << endl;
		}
		if (blank_line) Out << endl;
	} else if (gradients == 3) {
		int N_comp_ranges=grads[4];
		if (N_comp_ranges==1) {
			int pos=0;
			Matrix X(1,M,1,count);
			Vector profile(1,M);
			int k=0;
			for (int i=1; i<=OutputLineQ.Cardinal(); i++) {
				Line = (OutputLine*) OutputLineQ[i];
				if (!Line->IsProfile()) {
					continue;
				}
				bool write2 = true;
				Line = (OutputLine*) OutputLineQ[i];
				if (templateSet) {
					write2 = ValueInTemplate(Line);
				}
				if (write2) {
					k++;
					profile = Line->GetProfile();
					for (int z=1; z<=M; z++) X[z][k] = profile[z];
				}
			}
			for (int x=1; x<=grads[1]; x++) {
				for (int y=1; y<=grads[2]; y++) {
					for (int z=1; z<=grads[3]; z++) {
						pos++;
						if ( (x>1 && x<grads[1] && y>1 && y<grads[2] && z>1 && z<grads[3])) {
							Out << x-1 << '\t' << y-1 << '\t'<< z-1 << '\t';
							for (int k=1; k<=count; k++) {
								Out << X[pos][k] << '\t';
							}
							//for (int i=1; i<=OutputLineQ.Cardinal(); i++) {
							//	Line = (OutputLine*) OutputLineQ[i];
							//	if (!Line->IsProfile()) {
							//		continue;
							//	}
							//	bool write2 = true;
							//	Line = (OutputLine*) OutputLineQ[i];
							//	if (templateSet) {
							//		write2 = ValueInTemplate(Line);
							//	}
							//	if (write2) {
							//		Vector profile = Line->GetProfile();
							//		Out << profile[pos] << '\t';
							//	}
							//}
							Out << endl;
						}
					}
				}
			}
		}
		else {
			Message(warning,"Output option profiles in mulitple ranges not implemented");

			//int i,jx,jy,x0,y0,z0,xm,ym,zm,Current=0,Next=0;
			//for (i=1; i<=N_comp_ranges; i++){
			//	x0=grads[4+(i-1)*6+1]; xm=grads[4+(i-1)*6+4];
			//	y0=grads[4+(i-1)*6+2]; ym=grads[4+(i-1)*6+5];
			//	z0=grads[4+(i-1)*6+3]; zm=grads[4+(i-1)*6+6];
//
//				jy=(zm-z0+3);
//				jx=jy*(ym-y0+3);
//				Current=Next; Next+=jx*(xm-x0+3);
//
//				for (int x=x0; x<=xm; x++) {
//					for (int y=y0; y<=ym; y++) {
//						for (int z=z0; z<=zm; z++) {
//							Out << x << '\t' << y << '\t'<< z << '\t';
//							for (int k=1; k<=OutputLineQ.Cardinal(); k++) {
//								Line = (OutputLine*) OutputLineQ[k];
//								if (!Line->IsProfile()) {
//									continue;
//								}
//								bool write2 = true;
//								Line = (OutputLine*) OutputLineQ[k];
//								if (templateSet) {
//									write2 = ValueInTemplate(Line);
//								}
//								if (write2) {
//									Vector profile = Line->GetProfile();
//									Out << profile[Current+(x-x0+1)*jx+(y-y0+1)*jy+(z-z0+1)+1] << '\t';
//									//possibly we need to put the elements one down...
//									//writebounds true ignored;
//								}
//							}
//							Out << endl;
//						}
//					}
//				}
//			}
		}
	}
#ifdef _WIN32
	Text temp = Copy("attrib -r \"" + name+"\"");
	system(temp.MainC());
#endif
}
void
Output::OutputTemplate(Text name, bool firstWrite) {
	ofstream Out;
	open(Out, name);
	OutputLine* Line;
	for (int i=1; i<=OutputLineQ.Cardinal(); i++) {
		bool write2 = true;
		Line = (OutputLine*) OutputLineQ[i];
		if (templateSet) {
			write2 = ValueInTemplate(Line);
		}
		if (write && write2) {
			if (!Line->IsProfile() || outputProfiles) {
				Out << Line->GetElement1().MainC() << " : ";
				Out << Line->GetElement2().MainC() << " : ";
				Out << Line->GetElement3().MainC() << " : ";
				Out << Line->GetValue().MainC() << endl;
			}
		}
	}
	if (write) {
		Out << "start" << endl;
	}
#ifdef _WIN32
	Text temp = Copy("attrib -r \"" + name+"\"");
	system(temp.MainC());
#endif

}
bool
Output::ValueInTemplate(const OutputLine* Line) {
	return Template->ValueSetWild(Line->GetElement1(),
									   Line->GetElement2(),
									   Line->GetElement3());
}
bool
Output::ValueForLayerInTemplate(const OutputLine* Line, int z) {
	Text number = Blanks(100);
	number.Putint(z);
	number = Copy(number.Strip());
	number = Copy(number.Frontstrip());
	return   Template->ParamSetToValue(Line->GetElement1(),
									   Line->GetElement2(),
									   Line->GetElement3(),
									   number);
}
bool
Output::FileExists(Text name) {
#ifndef _WIN32
	ifstream In(parseFileName(name,MyInput).MainC());
	return !(!In);
#else
	FILE *fp=fopen(parseFileName(name,MyInput).MainC(),"r");
	if ( fp != NULL) fclose(fp);
	return fp != NULL;
#endif
}
Text
Output::getNumberedFilename(const Text name) const {
	int digits;
	Text newname, insert;
	char *namecpy, *extPtr;

	// make the part to be inserted
	if (MyInput->TotalNumCalculations() > 1) {
		digits = int(log10(1.*MyInput->TotalNumCalculations()))+1;
		insert = Blanks(100);
		insert.Putint(MyInput->GetCurrNumCalculation());
		insert = insert.Frontstrip();
		while (insert.Length() < digits) {
			insert = "0" + insert;
		}
		insert = "_" + insert;
	}
	if (MCS>0) {
	    digits = int(log10(MCS/MC_output_interval))+1;
	    insert = Blanks(100);
		insert.Putint(MCS/MC_output_interval);
		insert = insert.Frontstrip();
		while (insert.Length() < digits) {
			insert = "0" + insert;
		}
		insert = "_" + insert;
	}

	if (*(MyVary->GetVariableName()) != *Copy("")) {
		// find out how to format the extension
		// if it's an integer, it's straightforward
		insert = insert + "_" + MyVary->GetVariableName() + "="
			+ MyVary->GetCurrValue();
	}

	// insert the new part just before the file extention
	namecpy = new char[strlen(name.MainC())+1];
	strcpy(namecpy, name.MainC());
	extPtr=strchr(namecpy, '.');
	newname=Copy(name);
	if (extPtr != NULL) {
		while (strchr(extPtr+1, '.') != NULL) {
			extPtr=strchr(extPtr+1, '.');
		}
		newname=newname.Sub(1,extPtr-namecpy);
	}
	newname = newname + insert;
	if (extPtr != NULL) {
		newname = newname + extPtr;
	}
	delete [] namecpy;

/*	Text number = Blanks(100);

	number.Putint(MyInput->GetCurrNumCalculation());
	cout << "number" << MyInput->GetCurrNumCalculation() << endl;
	name = name + "_" + MyVary->GetVariableName() + MyVary->GetCurrValue()
				+ "_" + Copy(number.Strip().Frontstrip());*/
	return newname;
}

Boolean
Output::CheckReal(Text t) const {
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


/* This function opens the given output stream into a given filename.
if the boolean firstwrite is set, the file is created, else, it is appended to
*/
void
Output::open(ofstream & o, Text &name, bool firstWrite) {

	name = parseFileName(name,MyInput);

	if (append || !firstWrite) {
		o.open(name.MainC(),ios::app);
	} else {
		o.open(name.MainC());
	}
	o.precision(accuracy-1);
	o.setf(ios::showpoint);
	o.setf(ios::scientific);
}
/*
This function parses the filename. The parsed name is returned.
If the filename starts with 'filename' this part is replaced by the part of
the inputfile's name a dot, or in the absense of a dot the entire filename.
*/
Text
parseFileName(const Text &name, const Input *MyInput) {
	const char *in, *end, *orig;
	char *n;
	size_t p, size;
	Text newname;

	n = NULL;
	orig = name.MainC();
	if (!strncmp(orig, "filename", 8)) {
		orig += 8;
#ifndef _WIN32
		in = basename(MyInput->GetFileName().MainC());
#else
		//Get the basename from a DOS style name
		Text t = MyInput->GetFileName();
		for (int i=t.Length(); i>0; i--) {
			if ( *t.Sub(i,1) == "\\" ) {
				t = Copy(t.From(i+1));
				break;
			}
		}
		in = t.MainC();
#endif
		end = strrchr(in, '.');
		if (end == NULL) {
			end = in + strlen(in);
		}

		size = end-in + strlen(orig) + 1;
		n = new char[size];
		n[size-1] = '\0';
		for (p=0; p<(size_t)(end-in); p++) {
			n[p] = in[p];
		}
		for (p=0; p<strlen(orig); p++) {
			n[p+end-in] = orig[p];
		}
		newname = n;
		delete [] n;


		return newname;
	}
		return name;
}

