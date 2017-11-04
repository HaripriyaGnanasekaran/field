#include "SF_System3rdGen.h"

SF_System3rdGen::SF_System3rdGen(Input* MyInput_, 
					 			 Text name_,
					 			 Boolean compute_)
	: SF_System(MyInput_,name_,compute_) {
	if (approach == firstOrder) {
		Text latName;
		MyInput->GetNumNames("lat",1,1);
		Array<Text> latNames = MyInput->GetNames("lat");
		latName = latNames[1];
		Text n_layers_x = MyInput->GetText("lat",latName,"n_layers_x");
		MyInput->SetVariable("lat",latName,"n_layers_x","1");
		LatRN = NewLat1stO(latName);
		SegQRN = new SF_SegmentList(LatRN,MyInput);
		MolQRN = new SF_MolList1stO(SegQRN,LatRN,MyInput);
		ReactionQRN = new SF_ReactionList(MolQRN,SegQRN,MyInput);
		MyInput->SetVariable("lat",latName,"n_layers_x",n_layers_x);
		Solve3rdGen = NewSolve();
	}		
}
SF_System3rdGen::~SF_System3rdGen() {
	delete Solve3rdGen;	
	delete ReactionQRN;
	delete MolQRN;
	delete SegQRN;
	delete LatRN;
}
void
SF_System3rdGen::Go(Output* Out,Lattice* Lat) const {
	Solve3rdGen->Iterate();
	GetOutput(Out,Lat);
	Out->WriteOutput();
}
SF_Solve*
SF_System3rdGen::GetSolve() const {
	return Solve3rdGen;
}
void
SF_System3rdGen::GetOutput(Output* Out,Lattice* Lat) const {
	Out->PutText("sys",name,"inputfile",MyInput->GetFileName());
	Out->PutText("sys",name,"calculation_type","third_generation");
	if (approach == firstOrder) {
		Out->PutText("sys",name,"matrix_approach","first_order");
	}
	Out->PutBoolean("sys",name,"overflow_protection",overflowProtection);
	Out->PutBoolean("sys",name,"change_chi_with_temperature",changeChiWithTemperature);
	Out->PutReal("sys",name,"temperature",TEMPERATURE);
	Out->PutReal("sys",name,"free energy",MolQ->GetFreeEnergy());
	Out->PutReal("sys",name,"excess free energy",MolQ->GetExcessFreeEnergy());
	Out->PutProfile("sys",name,"free energy density",MolQ->GetFreeEnergyProfile());

	// update
	//	Out->PutProfile("sys",name,"pressure",MolQ->GetPressureProfile());
	if (SegQ->Charged()) {
		Out->PutProfile("sys",name,"potential",SegQ->GetElectricPotential());
		Out->PutProfile("sys",name,"charge",SegQ->GetCharge());
		Out->PutProfile("sys",name,"epsilon",SegQ->GetAverageEpsilon());
	}
	if (MyInput->ValueSet("sys",name,"extra_output")) {
		GetExtraOutput(Out,GetLattice());
	}
	Lat->GetOutput(Out);
	MolQ->GetOutput(Out);
	SegQ->GetOutput(Out);
	ReactionQ->GetOutput(Out);
	Solve->GetOutput(Out);
	Solve->WriteInitialGuessToFile();
}
SF_Solve3rdGen*
SF_System3rdGen::NewSolve() const {
	return new SF_Solve3rdGen(compute,
							  ReactionQ,
							  MolQ,
							  SegQ,
							  LatRN,
							  ReactionQRN,
							  MolQRN,
							  SegQRN,
							  LatRN,
							  MyInput);
}

			 
