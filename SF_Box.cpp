#define _GNU_SOURCE 1
#include <signal.h>
#ifdef USE_FENV
#include <fenv.h>
#endif
#include <stdio.h>
#include "SF_Box.h"
#include "Vary.h"
#include "SF_System.h"
#include "SF_System3rdGen.h"
#include "SF_SystemFraaije.h"
#include "SF_SystemMC.h"
#include "SF_SystemLD.h"
#ifdef _OPENMP
#include <omp.h>
#endif
Boolean help = false;
int num_calculations=0;
Boolean check = false;

Boolean version = false;

void fpe_handler(int);

void fpeon() {
#if USE_FENV
	signal(SIGFPE, fpe_handler);
	feenableexcept(FE_DIVBYZERO | FE_INVALID
			| FE_UNDERFLOW
			 | FE_OVERFLOW
		);
	feclearexcept (FE_ALL_EXCEPT);
#endif
}
void
fpeoff() {
#if USE_FENV
	fedisableexcept(FE_ALL_EXCEPT);
#endif
}

void fpe_handler(int) {
#if USE_FENV
	int e;
	printf("caught FPE, exiting.\n");
	feclearexcept (FE_ALL_EXCEPT);
/*	e = fetestexcept(FE_ALL_EXCEPT);
	if (!e) {
		printf("no exception information set\n");
	}
	if (e & FE_DIVBYZERO) {
		printf("divide by zero\n");
	}
	if (e & FE_UNDERFLOW) {
		printf("underflow\n");
	}
	if (e & FE_OVERFLOW) {
		printf("overflow\n");
	}*/
//	exit(1);
#endif
}

Output*
NewOutput(SF_System* Sys,
		  Input* MyInput,
		  Vary* MyVary,
		  Array<Text> FileNames) {
	Lattice* MyLat= Sys->GetLattice();
	SF_Solve* Solve= Sys->GetSolve();
	int gradients = MyLat->GetNumGradients();
	int N_comp_ranges = MyLat->Get_N_comp_ranges();
	int N_var=0;
	if (N_comp_ranges>1) {N_var=6*N_comp_ranges;};
	Array<int> grad(0,4+N_var);
	for (int i=1; i<=gradients; i++) {
		grad[i] = MyLat->GetNumLayers(i);
	}
	grad[0]=gradients;
	grad[4]=N_comp_ranges;
	if (N_comp_ranges>1) {
		for (int i=1; i<=N_comp_ranges; i++){

			grad[4+(i-1)*6+1]=MyLat->GetR0((i-1)*3+1);
			grad[4+(i-1)*6+2]=MyLat->GetR0((i-1)*3+2);
			grad[4+(i-1)*6+3]=MyLat->GetR0((i-1)*3+3);
			grad[4+(i-1)*6+4]=MyLat->GetRm((i-1)*3+1);
			grad[4+(i-1)*6+5]=MyLat->GetRm((i-1)*3+2);
			grad[4+(i-1)*6+6]=MyLat->GetRm((i-1)*3+3);
		}
	}

	int accuracy = int(-log10(Solve->GetTolerance())+0.5);
	return new Output(MyInput,MyVary,grad,accuracy,FileNames);

}
SF_System*
NewSystem(Input* MyInput, Boolean compute) {
	int num = MyInput->GetNumNames("sys",0,1);
	if (num == 0) {
		return new SF_System(MyInput,"noname",compute);
	} else {
		Array<Text> sysNames = MyInput->GetNames("sys");
		Text sysName = sysNames[1];
		sysNames.Dim(1,6);
		sysNames[1] = "equilibrium";
		sysNames[2] = "third_generation";
		sysNames[3] = "dynamic";
		sysNames[4] = "steady-state";
		sysNames[5] = "MC_SCF";
		sysNames[6] = "LD_SCF";
		switch (MyInput->GetChoice("sys",sysName,"calculation_type",sysNames,1)) {
			case 1:
				return new SF_System(MyInput,sysName,compute);
				break;
			case 2:
				return new SF_System3rdGen(MyInput,sysName,compute);
				break;
			case 3:
				return new SF_SystemFraaije(MyInput,sysName,compute);
				break;
			case 4:
				Message(fatal,"steady-state is not implemented yet");
				break;
			case 5:
				return new SF_SystemMC(MyInput,sysName,compute);
				break;
			case 6:
				return new SF_SystemLD(MyInput,sysName,compute);
				break;
			default:
				Message(fatal,"Error in SF_Box::NewSystem, undefined sys type");
				break;
		}
	}
	return NULL; // never get here
}
void
Go(Input* MyInput, Boolean compute = true) {
	Output* Out;
	SF_System* Sys=NULL;
	SF_System* SysOld=NULL;
	Vary* variate;
	Array<Text> OutputFileNames(0,0);
	while (MyInput->NextCalculation()) {
		variate = new Vary(MyInput);
		Boolean lastChecked = false;
		Boolean varying = false;

		while (variate->NextCalculation()) {
			if (!compute) num_calculations++;
			if (MyInput->GetCurrNumCalculation()>1 && compute) {
				cout << endl;
				cout << "Next problem: "<< MyInput->GetCurrNumCalculation() << " out of " << num_calculations << endl;
			}
			if (!compute && !lastChecked) {
				variate->SkipToLast();
				lastChecked = true;
			}
			Sys = NewSystem(MyInput,compute);
 			if (MyInput->GetCurrNumCalculation() == 1 && !varying) {
 	  			SysOld = Sys;
 			}
			Sys->SetInitialGuess(SysOld,Sys->GetLattice());
			if (MyInput->GetCurrNumCalculation() > 1 || varying) {
	 			delete SysOld;
 			}
			Out = NewOutput(Sys,MyInput,variate,OutputFileNames);
			if (compute) {
				Sys->Go(Out,Sys->GetLattice());
			}
			OutputFileNames = Out->GetFileNames();
			delete Out;
			SysOld = Sys;
			varying = true;
		}
		delete variate;
	}
	delete Sys;
}
void
PrintError(Text arg0) {
	Sysout().Outtext("Usage: " + arg0 + " [options] + inputfile. Use '" +
		arg0 + " -h' for details.");
	Sysout().Outimage();
	exit(-1);
}
void
ReadFlags(Text flags, Text arg0) {
	flags = Copy(flags.Frontstrip().Strip());
	flags.Setpos(1);
	if (flags.Getchar() == '-') {
		if (!flags.More()) {
			PrintError(arg0);
		}
		for (int i=2; i<=flags.Length(); i++) {
			char flag = flags.Getchar();
			if (flag == 'h') {
				help = true;
			} else if (flag == 'v') {
				version = true;
			} else if (flag == 'c') {
				check = true;
			} else if (flag == 'w') {
				warnings = true;
			} else if (flag == 's') {
				silent = true;
			} else {
				PrintError(arg0);
			}
		}
	} else {
		PrintError(arg0);
	}
}
void
PrintHelp(Text arg0) {
	Sysout().Outtext("Usage:");
	Sysout().Outimage();
	Sysout().Outtext(arg0 + " INPUTFILE");
	Sysout().Outimage();
	Sysout().Outtext("    performs the calculations read from INPUTFILE");
	Sysout().Outimage();
	Sysout().Outtext(arg0 + " -v");
	Sysout().Outimage();
	Sysout().Outtext("    prints version");
	Sysout().Outimage();
	Sysout().Outtext(arg0 + " -h");
	Sysout().Outimage();
	Sysout().Outtext("    prints this message");
	Sysout().Outimage();
	Sysout().Outtext(arg0 + " -c INPUTFILE");
	Sysout().Outimage();
	Sysout().Outtext("    only checks the syntax of the calculations read from INPUTFILE");
	Sysout().Outimage();
	Sysout().Outtext(arg0 + " -w INPUTFILE");
	Sysout().Outimage();
	Sysout().Outtext("    print many warnings for default values");
	Sysout().Outimage();
	Sysout().Outtext(arg0 + " -s INPUTFILE");
	Sysout().Outimage();
	Sysout().Outtext("    silent mode: no output to screen");
	Sysout().Outimage();
	Sysout().Outimage();
	Sysout().Outtext("combinations of the '-' flags can also be made, example:");
	Sysout().Outimage();
	Sysout().Outtext(arg0 + " -cw INPUTFILE");
	Sysout().Outimage();
	Sysout().Outtext("    check only with many warnings");
	Sysout().Outimage();
	Sysout().Outimage();
	Sysout().Outtext("Information on how to construct inputfiles can be in helpfile ");
	Sysout().Outimage();
	Sysout().Outtext("contact: frans.leermakers@wur.nl");
	Sysout().Outimage();
}
void
PrintVersion(Text arg0) {
	Sysout().Outtext(arg0 + " version " + GetVersion() + " date: " + GetDate());
	Sysout().Outimage();
}
int
main (int argc, char *argv[]) {
//	fpeon();
	Text inputFileName;
	Text flags;
	Boolean filePresent = false;
	warnings = false;
	silent = false;
	if (argc == 2) {
		Text arg1 = Copy (argv[1]);
		arg1 = Copy(arg1.Frontstrip().Strip());
		if (arg1.Getchar() == '-') {
			ReadFlags(argv[1],argv[0]);
			if (check || warnings || silent) {
				PrintError(argv[0]);
			}
			if (warnings && silent) {
				PrintError(argv[0]);
			}
		} else {
			inputFileName = arg1;
			filePresent = true;
		}
	} else if (argc == 3) {
		ReadFlags(argv[1],argv[0]);
		if (help || version) {
			PrintError(argv[0]);
		}
		inputFileName = Copy (argv[2]);
		inputFileName = inputFileName.Strip();
		filePresent = true;
	} else {
		PrintError(argv[0]);
	}
	if (help) {
		PrintHelp(argv[0]);
		return 0;
	}
	if (version && !silent) {
		PrintVersion(argv[0]);
	}
	if (filePresent) {
		if (!version && !silent) {
			PrintVersion(argv[0]);
		}
		if (!silent) {
#ifdef _OPENMP
                     printf("Maximum number of threads: %d \n",omp_get_max_threads());
#endif
			Sysout().Outtext("Checking inputfile " + inputFileName + " for errors");
			Sysout().Outimage();
		}
		Input* MyInput;
		MyInput = new Input(inputFileName);
		Array<Text> firstArguments(1,10);
		firstArguments[1] = "lat";
		firstArguments[2] = "mon";
		firstArguments[3] = "state";
		firstArguments[4] = "mol";
		firstArguments[5] = "newton";
		firstArguments[6] = "output";
		firstArguments[7] = "sys";
		firstArguments[8] = "reaction";
		firstArguments[9] = "var";
		firstArguments[10] = "alias";
		MyInput->SetAllMessagesOn();
		if (!warnings) {
			MyInput->SetAllMessagesOff();
		}
		MyInput->CheckFirstArguments(firstArguments);
		Go(MyInput,false);

		if (!silent) {
			Sysout().Outtext("no errors found in " + inputFileName);
			Sysout().Outimage();
		}
		int num = MyInput->GetNumNames("sys",0,1);
		if (num >= 1) {
			Array<Text> sysNames = MyInput->GetNames("sys");
			Text sysName = sysNames[1];
			uint random_seed = MyInput->GetInt("sys",sysName,"random_seed",0,INT_MAX-2,2348725);
			srand(random_seed+2);	// srand(1) resets the random number iterator rather than defining a new randomseed therefore we add two so it still seems logical to the user.
		}
		delete MyInput;
	}

	if (filePresent && !check) {
		if (!silent) {
			Sysout().Outtext("Here we go...");
			Sysout().Outimage();
		}
		Input* MyInput;
		MyInput = new Input(inputFileName);
		warnings = false;
		silentForced = silent;
		silent = true;
		MyInput->SetAllMessagesOff();

		Go(MyInput,true);

		delete MyInput;
		if (!silentForced) {
			Sysout().Outtext(argv[0]);
			Sysout().Outtext(" has finished.");
			Sysout().Outimage();
		}
	}

	return 0;
}
