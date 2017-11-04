#include "SF_SystemFraaije.h"
#include "SF_SolveFraaije.h"

SF_SystemFraaije::SF_SystemFraaije(Input* MyInput_, 
								   Text name_, 
								   Boolean compute_) {
	MyInput = MyInput_;
	name = name_;
	compute = compute_;
	Array<Text> param(1,17);
	param[1] = "matrix_approach";
	param[2] = "overflow_protection";
	param[3] = "calculation_type";
	param[4] = "temperature";
	param[5] = "change_chi_with_temperature";
	param[6] = "warnings";
	param[7] = "extra_output";
	param[8] = "time";
	param[9] = "time_step";
	param[10] = "time_limit";
	param[11] = "time_step_incr";
	param[12] = "min_error";
	param[13] = "max_error";
	param[14] = "stop_criterion";
	param[15] = "output_timer_limit";
	param[16] = "iterate_lattice_artefact";
	param[17] = "random_seed";
	error = 1;
	stop = 1;
	MyInput->CheckParameterNames("sys",name,param);
	MyInput->SetAllMessagesOff();
	if (MyInput->ValueSet("sys",name,"warnings")) {
		if (MyInput->GetBoolean("sys",name,"warnings")) {
			MyInput->SetAllMessagesOn();
		}
	}
	Array<Text> approachParam(1,1);
	approachParam[1] = "first_order";
	int numApproach = MyInput->GetChoice("sys",name,"matrix_approach",approachParam,1);
	if (numApproach == 1) {
		approach = firstOrder;
	}
	overflowProtection = MyInput->GetBoolean("sys",name,"overflow_protection",false);	
	temperature = MyInput->GetReal("sys",name,"temperature",0,DBL_MAX,TEMPERATURE);
	changeChiWithTemperature = MyInput->GetBoolean("sys",name,"change_chi_with_temperature",false);
	timeNow = MyInput->GetReal("sys",name,"time",0,DBL_MAX,0);
	timeStep = MyInput->GetReal("sys",name,"time_step",0,DBL_MAX,0.01);
	timeLimit = MyInput->GetReal("sys",name,"time_limit",0,DBL_MAX,DBL_MAX);
	timeStepIncr = MyInput->GetReal("sys",name,"time_step_incr",0,DBL_MAX,1.02);
	minError = MyInput->GetReal("sys",name,"min_error",0,DBL_MAX,1e-4);
	maxError = MyInput->GetReal("sys",name,"max_error",minError,DBL_MAX,1e-3);
	stopCriterion = MyInput->GetReal("sys",name,"stop_criterion",0,1,1e-9);
	outputTimerLimit = MyInput->GetInt("sys",name,"output_timer_limit",0,10000,10);
	if (approach == firstOrder) {
		Text latName;
		MyInput->GetNumNames("lat",1,1);
		Array<Text> latNames = MyInput->GetNames("lat");
		latName = latNames[1];
		Lat1st = NewLat1stO(latName);
		SegQ = new SF_SegmentList(Lat1st,MyInput);
		MolQ = new SF_MolList1stO(SegQ,Lat1st,MyInput);
		ReactionQ = new SF_ReactionList(MolQ,SegQ,MyInput);
		SegNewQ = new SF_SegmentList(Lat1st,MyInput);
		MolNewQ = new SF_MolList1stO(SegNewQ,Lat1st,MyInput);
		ReactionNewQ = new SF_ReactionList(MolNewQ,SegNewQ,MyInput);
		SolveFraaije = NewSolve();
	}
	if (Lat1st->GetNumGradients() > 1) {
		Message(fatal,"dynamics not implemented for 2 gradients");
	}
	if (changeChiWithTemperature) {
		UpdateChi();
	}
	TEMPERATURE = temperature;
	iterateLatticeArtefact = MyInput->GetBoolean("sys",name,"iterate_lattice_artefact",false);

}
SF_SystemFraaije::~SF_SystemFraaije() {
	delete ReactionNewQ;
	delete MolNewQ;
	delete SegNewQ;
	delete SolveFraaije;
		Solve = new SF_Solve(false,ReactionQ,MolQ,SegQ,Lat1st,MyInput);
}
void
SF_SystemFraaije::Go(Output* Out,Lattice* Lat) {
	SolveFraaije->Iterate(true,timeStep);
	SolveFraaije->UpdateDensities();
	GetOutput(Out,Lat);
	Out->WriteOutput();
	stop = 1;
	error = 1;
	int outputTimer = 0;
	Boolean errorOutsideBounds = false;
	while (timeNow < timeLimit && stop > stopCriterion) {
		SolveFraaije->Iterate(false,timeStep);
		error = SolveFraaije->CalcError();
		stop = SolveFraaije->CalcStop();
		while ((error > maxError || error < minError) && timeStep <= timeNow && stop > stopCriterion) {
			if (error > maxError) {
				if (errorOutsideBounds) {
					timeStepIncr = 1 + (timeStepIncr-1)/2;
				}
				timeStep /= timeStepIncr;
				errorOutsideBounds = false;
			}
			if (error < minError) {
				if (errorOutsideBounds) {
					timeStepIncr = 1 + (timeStepIncr-1)*1.5;
				}
				errorOutsideBounds = true;
				timeStep *= timeStepIncr;
			}
			SolveFraaije->Iterate(false,timeStep);
			error = SolveFraaije->CalcError();
			stop = SolveFraaije->CalcStop();
			Sysout().Outreal(timeNow,4,0);
			Sysout().Outchar('\t');
			Sysout().Outreal(timeStep,4,0);
			Sysout().Outchar('\t');
			Sysout().Outreal(error,4,0);
			Sysout().Outchar('\t');
			Sysout().Outreal(stop,4,0);
			Sysout().Outimage();
		}
		SolveFraaije->UpdateDensities();
		timeNow += timeStep;
		Sysout().Outreal(timeNow,4,0);
		Sysout().Outimage();
		outputTimer++;
		if (outputTimer == outputTimerLimit || stop <= stopCriterion || timeNow >= timeLimit) {
			Out->Clear();
			GetOutput(Out,Lat);
			Out->WriteOutput();
			outputTimer = 0;
		}
	}
}
SF_Solve*
SF_SystemFraaije::GetSolve() const {
	return SolveFraaije;
}
void
SF_SystemFraaije::SetInitialGuess(SF_System* Old,Lattice* Lat) {
	SolveFraaije->SetInitialGuess(Old->GetSolve());
}
SF_SolveFraaije*
SF_SystemFraaije::NewSolve() const {
	return new SF_SolveFraaije(compute,
							   ReactionQ,
							   MolQ,
							   SegQ,
							   ReactionNewQ,
							   MolNewQ,
							   SegNewQ,
							   Lat1st,
							   MyInput);
}
void
SF_SystemFraaije::UpdateChi() {
	Matrix chi = SegQ->GetChiMatrix();
	Matrix chiNew = SegNewQ->GetChiMatrix();
	int max=chi.IMAX;
	for (int i=1; i<=max; i++) {
		for (int j=1; j<=max; j++) {
			chi[i][j] *= TEMPERATURE/temperature;
			chiNew[i][j] *= TEMPERATURE/temperature;
		}
	}
}
void
SF_SystemFraaije::GetOutput(Output* Out,Lattice* Lat) const {
	Out->PutText("sys",name,"inputfile",MyInput->GetFileName());
	Out->PutText("sys",name,"calculation_type","dynamic");
	if (approach == firstOrder) {
		Out->PutText("sys",name,"matrix_approach","first_order");
	}
	Out->PutBoolean("sys",name,"overflow_protection",overflowProtection);
	Out->PutBoolean("sys",name,"change_chi_with_temperature",changeChiWithTemperature);
	Out->PutReal("sys",name,"temperature",TEMPERATURE);
	Out->PutReal("sys",name,"free energy",MolQ->GetFreeEnergy());
	Out->PutReal("sys",name,"excess free energy",MolQ->GetExcessFreeEnergy());
	Out->PutReal("sys",name,"time",timeNow);
	Out->PutReal("sys",name,"time_step",timeStep);
	Out->PutReal("sys",name,"time_limit",timeLimit);
	Out->PutReal("sys",name,"time_step_incr",timeStepIncr);
	Out->PutReal("sys",name,"error",error);
	Out->PutReal("sys",name,"min_error",minError);
	Out->PutReal("sys",name,"max_error",maxError);
	Out->PutReal("sys",name,"stop",stop);
	Out->PutReal("sys",name,"stop_criterion",stopCriterion);
	Out->PutReal("sys",name,"output_timer_limit",outputTimerLimit);
//	Vector pressureProfile = MolQ->GetPressureProfile();
//	Lat->DivideByLatticeSites(pressureProfile);
//	Out->PutReal("sys",name,"first moment of pressure",Lat->Moment(,1,0));
//	Out->PutReal("sys",name,"second moment of pressure",Lat->Moment(MolQ->GetPressureProfile(),2,0));
	Out->PutProfile("sys",name,"free energy density",MolQ->GetFreeEnergyProfile());
	//  update needed
	//	Out->PutProfile("sys",name,"pressure",MolQ->GetPressureProfile());
	if (SegNewQ->Charged()) {
		Out->PutProfile("sys",name,"potential",SegNewQ->GetElectricPotential());
		Out->PutProfile("sys",name,"charge",SegNewQ->GetCharge());
		Out->PutProfile("sys",name,"epsilon",SegNewQ->GetAverageEpsilon());
	}
	if (MyInput->ValueSet("sys",name,"extra_output")) {
		GetExtraOutput(Out,GetLattice());
	}
	Lat1st->GetOutput(Out);
	MolNewQ->GetOutput(Out);
	SegNewQ->GetOutput(Out);
	ReactionNewQ->GetOutput(Out);
	SolveFraaije->GetOutput(Out);
	SolveFraaije->WriteInitialGuessToFile();
}
			 
