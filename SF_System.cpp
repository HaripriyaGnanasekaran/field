#include "SF_System.h"
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <math.h>
#include "tools.h"
#ifdef CUDA
#include "Cuda_tools.h"
#endif

const double SF_System::ROOMTEMPERATURE = 298.15;

SF_System::SF_System(Input* MyInput_, Text name_, Boolean compute_) {
	Cuda_enabled=false;
	MyInput = MyInput_;
	name = name_;
	compute = compute_;

#ifdef CUDA
	Cuda_enabled = true;
#endif

	Array<Text> param(1,24);
	param[1] = "matrix_approach";
	param[2] = "overflow_protection";
	param[3] = "calculation_type";
	param[4] = "temperature";
	param[5] = "change_chi_with_temperature";
	param[6] = "extra_output";
	param[7] = "super_iteration";
	param[8] = "super_function";
	param[9] = "super_function_value";
	param[10] = "super_molecule";
	param[11] = "super_tolerance";
	param[12] = "iterate_lattice_artefact";
	param[13] = "calc_chem_pot_graft";
	param[14] = "lattice_artefact_tolerance";
	param[15] = "super_iteration_limit";
	param[16] = "core_monomer";
	param[17] = "probe_monomer";
	param[18] = "pinned_monomer";
	param[19] = "probe_chain";
	param[20] = "probe_surfactant";
	param[21] = "probe_surface";
	param[22] = "Gibbs_molecule";
	param[23] = "GPU";
	//param[24] = "MayerSaupe";
	param[24] = "random_seed";

	MyInput->CheckParameterNames("sys",name,param);

	Array<Text> approachParam(1,2);
	approachParam[1] = "first_order";
	approachParam[2] = "second_order";
	int numApproach = MyInput->GetChoice("sys",name,"matrix_approach",approachParam,1);
	if (numApproach == 1) {
		approach = firstOrder;
	}
	if (numApproach == 2) {
		approach = secondOrder;
	}
	//MayerSaupe = MyInput->GetBoolean("sys",name,"MayerSaupe",false);
	//if (MayerSaupe) {
	//	if (approach ==firstOrder) Message(fatal,"MayerSaupe can only be combined with second_order in matrix_approach");
	//}

	overflowProtection = MyInput->GetBoolean("sys",name,"overflow_protection",false);
	temperature = MyInput->GetReal("sys",name,"temperature",0.0001,DBL_MAX,ROOMTEMPERATURE);
	TEMPERATURE = temperature;

	changeChiWithTemperature = MyInput->GetBoolean("sys",name,"change_chi_with_temperature",false);
	if (approach == firstOrder) {
		Text latName;
		MyInput->GetNumNames("lat",1,1);
		Array<Text> latNames = MyInput->GetNames("lat");
		latName = latNames[1];
		Lat1st = NewLat1stO(latName);
		SegQ = new SF_SegmentList(Lat1st,MyInput);
		MolQ = new SF_MolList1stO(SegQ,Lat1st,MyInput);
		ReactionQ = new SF_ReactionList(MolQ,SegQ,MyInput);
		Solve = NewSolve(ReactionQ,MolQ,SegQ,Lat1st);
	}
	if (approach == secondOrder) {
		Text latName;
		MyInput->GetNumNames("lat",1,1);
		Array<Text> latNames = MyInput->GetNames("lat");
		latName = latNames[1];

		Lat2nd = NewLat2ndO(latName);
		SegQ = new SF_SegmentList(Lat2nd,MyInput);

		MolQ = new SF_MolList2ndO(SegQ,Lat2nd,MyInput);
		ReactionQ = new SF_ReactionList(MolQ,SegQ,MyInput);
		Solve = NewSolve(ReactionQ,MolQ,SegQ,Lat2nd);
	}
	if (changeChiWithTemperature) {
		UpdateChi();
	}
	ReactionQ->UpdateReactionConstants(ROOMTEMPERATURE, temperature);
	graftedMolecules = false;
	for (int i=1; i<=MolQ->GetNumMolecules(); i++) {
		SF_Molecule* Mol = MolQ->GetMolecule(i);
		if (Mol->GetMolStructure()->SomeSegmentsPinned() ||
			Mol->GetMolStructure()->SomeSegmentsGrafted() ||
			Mol->GetFreedom() == secondGeneration ||
			Mol->GetFreedom() == thirdGeneration) {
			graftedMolecules = true;
		}
	}
	if (MyInput->ValueSet("sys",name,"calc_chem_pot_graft")) {
		Message(literal,"Warning! \nDon't set 'sys : " + name + " : calc_chem_pot_graft,"
			" the chemical potential of grafted molecules is always calculated");
	}

	superIterate = MyInput->GetBoolean("sys",name,"super_iteration",false);
	if (superIterate) {
		MyInput->AlwaysCombineParam("sys",name,"super_iterate","super_function");
		MyInput->AlwaysCombineParam("sys",name,"super_iterate","super_function_value");
		MyInput->AlwaysCombineParam("sys",name,"super_iterate","super_molecule");
		MyInput->AlwaysCombineParam("sys",name,"super_iterate","super_iteration_limit");
		Text SuperMolName = MyInput->GetText("sys",name,"super_molecule");
		if (!MolQ->MoleculeDefined(SuperMolName)) {
			Message(fatal,"'sys : " + name + " : super_molecule : " + SuperMolName +
				"' is a molecule that is not defined");
		}
		SuperMol = MolQ->GetMolecule(SuperMolName);
		Array<Text> functionChoices(1,8);
		functionChoices[1] = "grand_potential";
		functionChoices[2] = "phibulk";
		functionChoices[3] = "FH-MU";
		functionChoices[4] = "theta";
		functionChoices[5] = "GPE";
		functionChoices[6] = "Equate_to_solvent";
		functionChoices[7] = "GibbsExc";
		functionChoices[8] = "P_Laplace";
		int choice = MyInput->GetChoice("sys",name,"super_function",functionChoices);
		if (choice == 1) superFunction = grandPotential;
		if (choice == 2) superFunction = phiBulk;
		if (choice == 3) superFunction = chemPot;
		if (choice == 4) superFunction = theta;
		if (choice == 5) superFunction = GPE;
		if (choice == 6) superFunction = ETS;
		if (choice == 7) {superFunction = GibbsExc; superFunctionValue = MyInput->GetReal("sys",name,"super_function_value",-DBL_MAX,DBL_MAX);}
		if (choice == 8) {superFunction = P_Laplace; superFunctionValue = MyInput->GetReal("sys",name,"super_function_value",-DBL_MAX,DBL_MAX);}
else
		superFunctionValue = MyInput->GetReal("sys",name,"super_function_value",-DBL_MAX,DBL_MAX);
		settolerance(MyInput->GetReal("sys",name,"super_tolerance",1e-13,1,1e-7));
		setiterationlimit(MyInput->GetInt("sys",name,"super_iteration_limit",0,INT_MAX,100));
	} else {
		MyInput->DontCombineParam("sys",name,"super_iterate","super_function");
		MyInput->DontCombineParam("sys",name,"super_iterate","super_function_value");
		MyInput->DontCombineParam("sys",name,"super_iterate","super_molecule");
		MyInput->DontCombineParam("sys",name,"super_iterate","super_tolerance");
		MyInput->DontCombineParam("sys",name,"super_iterate","super_iteration_limit");
		if (( MyInput->ValueSet("sys",name,"super_function")
			    || MyInput->ValueSet("sys",name,"super_function_value")
			    || MyInput->ValueSet("sys",name,"super_molecule")
				|| MyInput->ValueSet("sys",name,"super_tolerance")
				|| MyInput->ValueSet("sys",name,"super_iteration_limit"))
				&& !MyInput->ValueSet("sys",name,"super_iteration") ) {
				Message(fatal,"Set 'sys : " + name + " : super_iteration' to 'true'");
		}
	}
	iterateLatticeArtefact = MyInput->GetBoolean("sys",name,"iterate_lattice_artefact",false);
	if (iterateLatticeArtefact) {
		artefactTolerance = MyInput->GetReal("sys",name,"lattice_artefact_tolerance",0,1,sqrt(Solve->GetTolerance()));
		ArtefactIter = new NoArtefact(GetLattice(),MolQ,SegQ,Solve,artefactTolerance);
	} else {
		artefactTolerance = 0; // unused value
	}
}
SF_System::~SF_System() {
	if (iterateLatticeArtefact) {
		//delete ArtefactIter;
	}

	delete Solve;
	delete ReactionQ;
	delete MolQ;
	delete SegQ;
	delete GetLattice();
}



void
SF_System::Go(Output* Out, Lattice* Lat) {
	Vector oldPhiBulks;
	Array<MolFreedom> oldFreedoms;
	int numMol = MolQ->GetNumMolecules();
	//int numSeg = SegQ->GetNumSegments();
	//int numStates = SegQ->GetNumStates();

	if (iterateLatticeArtefact) {
		oldPhiBulks.Dim(1,numMol);
		oldFreedoms.Dim(1,numMol);
		if (Solve->e_info || Solve->s_info) {
			cout << "LATTICE ARTEFACT ITERATION started" << endl;
		}
		Solve->Iterate();
		for (int i=1; i<=numMol; i++) {
			SF_Molecule* Mol = MolQ->GetMolecule(i);
			oldFreedoms[i] = Mol->GetFreedom();
			oldPhiBulks[i] = Mol->GetPhiBulk();
			if (Mol->GetFreedom() != solvent && Mol->GetFreedom() != neutralizer) {
				Mol->SetFreedom(fixedTheta);
				Mol->SetTheta(Mol->ComputeTheta());
			}
		}
		Solve->samehessian = true;
		Solve->Iterate(); // this updates normalisations of volumefractions
		ArtefactIter->Go();
	} else {
		Solve->Iterate();

	}
	if (superIterate) {
		Vector x(1,1);
		setopenline("SUPER ITERATION: ");
		e_info = true;
		setdeltamax(0.5);
		if (superFunction == theta) {
			x[1] = log(SuperMol->GetPhiBulk()/(1-SuperMol->GetPhiBulk()));
		} else {
			x[1] = -log(SuperMol->GetTheta());
		}
		iterate(&x[1],1);
	}

	GetOutput(Out,Lat);

	Out->WriteOutput();
	if (iterateLatticeArtefact) {
		for (int i=1; i<=numMol; i++) {
			SF_Molecule* Mol = MolQ->GetMolecule(i);
			Mol->SetFreedom(oldFreedoms[i]);
			Mol->SetPhiBulk(oldPhiBulks[i]);
		}
	}
}
void
SF_System::SetInitialGuess(SF_System* Old, Lattice* Lat) {
	if (Solve->NewSolveNeeded(Old->GetSolve())) {
		Solve->SetInitialGuess(Old->GetSolve());
	} else {
		if (iterateLatticeArtefact) {
			delete ArtefactIter;
		}
		Old->Solve->UpdateSolve(Solve);
		// swap solves between the two systems
		SF_Solve* Temp = Old->Solve;
		Old->Solve = Solve;
		Solve = Temp;
		if (iterateLatticeArtefact) {
			ArtefactIter = new NoArtefact(Lat,MolQ,SegQ,Solve,artefactTolerance);
		}
	}
}
Lattice*
SF_System::GetLattice() const {
	if (approach == secondOrder) return Lat2nd;
	if (approach == firstOrder ) return Lat1st;
	return Lat1st; //do not reach this point;
}
SF_SegmentList*
SF_System::GetSegmentList() const {
	return SegQ;
}
SF_MoleculeList*
SF_System::GetMoleculeList() const {
	return MolQ;
}
SF_Solve*
SF_System::GetSolve() const {
	return Solve;
}

double
SF_System::GetGibbsExcess() {

	Lat1stO* Lat; //assuming first order here...
	Lat = Lat1st;
	if (Lat == NULL) Message(fatal,"Problem terminated because a lattice if 1st order was expected....");
	double V=Lat->GetTotalNumLatticeSites();
	int M = Lat->GetTotalNumLayers()-2;
	if (V < PI*M*M-1 || V > PI*M*M+1) Message(fatal,"Problem terminated because cylindrical geometry expected");

	Text probeName = MyInput->GetText("sys",name,"Gibbs_molecule");
	if (!MolQ->MoleculeDefined(probeName)) {Message(fatal,"sys : " + name + " : Gibbs_molecule : " + probeName + " is molecule that is not defined.");}
	SF_Molecule* Mol = MolQ->GetMolecule(probeName);

	Vector phi = Mol->GetPhi(total);
	double phi_low=1;
	double phi_high=0;;
	for (int z=3; z<M; z++) {
		if (phi[z] < phi_low) phi_low=phi[z];
		if (phi[z] > phi_high) phi_high = phi[z];
	}
	Lat->SubtractBoundaries(phi);

	double theta_exc = Mol->ComputeThetaExcess();
	double Dphi=phi_high-phi_low;
	double Gvolume=0;
	if (Dphi>0) Gvolume=theta_exc/Dphi; if (Gvolume < 0) Gvolume=-1.0*Gvolume;
	double RGibbs=0;
	//double offset = Lat->GetLayerAdjustment();

	if (Gvolume > 0) RGibbs=pow(Gvolume/PI,1.0/2.0);
	Lat->RestoreBoundaries(phi);
	double GExcess=0;
	Vector phisuper=SuperMol->GetPhi(total);
	double phi1=phisuper[1];
	double phim=phisuper[M];
	GExcess=SuperMol->ComputeThetaExcess()-PI*RGibbs*RGibbs*(phi1-phim);
	return -GExcess/(superFunctionValue*2*PI*RGibbs)+1.0;
}

void
SF_System::residuals(double *const f, double *const x) {
	samehessian=false;

	double ThetaSolvent;
	SF_Molecule* Mol;
	int numMol=0;
	if (superFunction == theta) {
		SuperMol->SetPhiBulk(exp(x[0])/(1+exp(x[0])));
	} else {
		SuperMol->SetTheta(exp(-x[0]));
	}
	Solve->samehessian = true;
	if (iterateLatticeArtefact) {
		ArtefactIter->Go();
	} else {
		Solve->Iterate();
	}
	switch(superFunction) {
		case grandPotential:
			f[0] = MolQ->GetGrandPotential() - superFunctionValue;
			break;
		case phiBulk:
			f[0] = SuperMol->GetPhiBulk()/superFunctionValue - 1;
			break;
		case chemPot:
			if (SuperMol->GetFreedom() == secondGeneration) {
				f[0] = MolQ->GetChemicalPotential(SuperMol,constrained)/superFunctionValue - 1;
			} else {
				f[0] = -1.0*(MolQ->GetChemicalPotential(SuperMol)/superFunctionValue - 1);
			}
			break;
		case theta:
			f[0] = SuperMol->ComputeTheta()/superFunctionValue - 1;
			break;
		case GPE:
			f[0] = MolQ->GetGrandPotentialExcess() - superFunctionValue;
			break;
		case ETS:
			numMol = MolQ->GetNumMolecules();
			ThetaSolvent =0;
			for (int i=1; i<=numMol; i++) {
				Mol = MolQ->GetMolecule(i);
				if (Mol->GetFreedom() == solvent) ThetaSolvent =Mol->ComputeTheta();
			}
			f[0] = ThetaSolvent/SuperMol->ComputeTheta() -1;
			break;
		case GibbsExc:
			f[0] = GetGibbsExcess() ;
			break;
		case P_Laplace:
			f[0] = -1000*(MolQ->GetGrandPotentialProfile())[2];
			break;
		default:
			Message(fatal,"error in SF_System::residuals()");
			break;
	}
}

Lat2ndO*
SF_System::NewLat2ndO(Text latName) const {
	Array<Text> geometries(1,3);
	geometries[1] = "flat";
	geometries[2] = "cylindrical";
	geometries[3] = "spherical";
	int geom = MyInput->GetChoice("lat",latName,"geometry",geometries,1); // default: flat (only choice)
	int grad = MyInput->GetInt("lat",latName,"gradients",1,3,1);

	Array<Text> lattype(1,1);
	lattype[1] = "standard";
	//lattype[2] = "stencils";
	//lattype[3] = "FCC";
	//lattype[4] = "HEX";
	//Boolean stencils = (MyInput->GetChoice("lat", latName, "latticetype", lattype, 1) == 2)
	//			|| MyInput->ValueSet("lat", latName, "lambda0");
	//Boolean FCC = MyInput->GetChoice("lat", latName, "latticetype", lattype, 1) == 3;
	//Boolean HEX = MyInput->GetChoice("lat", latName, "latticetype", lattype, 1) == 4;
	if (overflowProtection) {Message(fatal,"3D flat lattice with overflow protection not yet implemented.");}
	//if (FCC||HEX) {Message(fatal,"No 2nd order in FCC or HEX yet. Please use lattictype : standard");}}

	switch (geom) {
		case 1:
				if (grad==1) {
					return new Lat1D2ndO(MyInput,latName);
				}
				if (grad==2) {
					return new Lat2D2ndO(MyInput,latName);
				}
				if (grad==3) {
					return new Lat3D2ndO(MyInput,latName);
				}
			break;
		case 2:
				if (grad==1) {
					return new Lat1DCyl2ndO(MyInput,latName);
				}
				if (grad==2) {
					return new Lat2DCyl2ndO(MyInput,latName);
				}
				if (grad==3) Message(fatal,"in 3gradients no cylindrical geometry supported.");
			break;
		case 3:
				if (grad==1) {
					return new Lat1DSphere2ndO(MyInput,latName);
				} else {
					Message(fatal,"in 2 or 3 gradients no spherical geometry supported.");
				}
			break;
		default:
			Message(fatal,"Error in SF_System::LoadLat2ndO, undefined lattice");
			break;
	}
	return NULL; // program should never reach this line
}

Lat1stO*
SF_System::NewLat1stO(Text latName) const {
	Array<Text> geometries(1,3);
	geometries[1] = "flat";
	geometries[2] = "cylindrical";
	geometries[3] = "spherical";
	int grad = MyInput->GetInt("lat",latName,"gradients",1,3,1);
	int geom;
    if (grad<3) {geom = MyInput->GetChoice("lat",latName,"geometry",geometries,1);
    } else {geom = MyInput->GetChoice("lat",latName,"geometry",geometries,1);
		if (geom>1) Message(literal, " Geometry switched to flat: in 3G only flat geometry implemented. ");
		geom=1;
	}

	bool stencils=false;
	bool FCC=false, HEX=false;
	Array<Text> lattype(1,4);
	lattype[1] = "standard";
	lattype[2] = "stencils";
	lattype[3] = "FCC";
	lattype[4] = "HEX";
	stencils = (MyInput->GetChoice("lat", latName, "latticetype", lattype, 1) == 2)
				|| MyInput->ValueSet("lat", latName, "lambda0");

	FCC = MyInput->GetChoice("lat", latName, "latticetype", lattype, 1) == 3;
	HEX = MyInput->GetChoice("lat", latName, "latticetype", lattype, 1) == 4;
	switch (geom) {
		case 1:
			if (grad == 1 && !overflowProtection) {
				return new Lat1DFlat1stO(MyInput,latName);
			} else if (grad == 1 && overflowProtection) {
				return new Lat1DFlat1stOS(MyInput,latName);
			} else if (grad == 2 && !overflowProtection) {
				if (stencils) {
					return new Lat2DFlat1stOSten(MyInput,latName);
				} else {
					return new Lat2DFlat1stO(MyInput,latName);
				}
			} else if (grad == 2 && overflowProtection) {
				if (stencils) {
					Message(fatal,"2D flat lattice with stencils and "
						"overflow protection not yet implemented.");
				}
				return new Lat2DFlat1stOS(MyInput,latName);
			} else if (grad == 3) {
				if (overflowProtection) {
					Message(fatal,"3D flat lattice with overflow protection not yet implemented.");

				}

				if (FCC) { return new Lat3DFCC1stO(MyInput,latName);} else
				if (HEX) { return new Lat3DHEX1stO(MyInput,latName);} else
				{return new Lat3D1stO(MyInput,latName);}

			}
			break;
		case 2:
			if (grad == 1 && !overflowProtection) {
				return new Lat1DCyl1stO(MyInput,latName);
			} else if (grad == 1 && overflowProtection) {
				return new Lat1DCyl1stOS(MyInput,latName);
			} else if (grad == 2 && !overflowProtection) {
				if (stencils) {
					return new Lat2DCyl1stOSten(MyInput,latName);
				} else {
					return new Lat2DCyl1stO(MyInput,latName);
				}
			} else if (grad == 2 && overflowProtection) {
				if (stencils) {
					return new Lat2DCyl1stOStenS(MyInput,latName);
				}
				return new Lat2DCyl1stOS(MyInput,latName);
			}
			break;
		case 3:
			if (grad == 1 && !overflowProtection) return new Lat1DSphere1stO(MyInput,latName);
			if (grad == 1 && overflowProtection) return new Lat1DSphere1stOS(MyInput,latName);
			if (grad == 2 && !overflowProtection) Message(fatal,MyInput,"Can't make a 2D spherical lattice.");
			if (grad == 2 && overflowProtection) Message(fatal,MyInput,"Can't make a 2D spherical lattice.");
			if (grad == 3) Message(fatal,MyInput,"In 3D only flat lattice possible");
			break;
		default:
			Message(fatal,"Error in SF_System::LoadLat1stO, undefined lattice");
			break;
	}
	return NULL; // program should never reach this line
}
SF_Solve*
SF_System::NewSolve(SF_ReactionList* ReactionQTemp,SF_MoleculeList* MolQTemp,SF_SegmentList* SegQTemp,Lattice* LatTemp) const {
	int num = MyInput->GetNumNames("newton",0,1);
	if (num == 0) {
		return new SF_Solve(compute,ReactionQTemp,MolQTemp,SegQTemp,LatTemp,MyInput);
	} else {
		Array<Text> solveNames = MyInput->GetNames("newton");
		Text solveName = solveNames[1];
		solveNames.Dim(1,6);
		solveNames[1] = "standard";
		solveNames[2] = "copel";
		solveNames[3] = "phi-and-U";
		solveNames[4] = "Picard";
		solveNames[5] = "CG";
		solveNames[6] = "Johan";

		switch (MyInput->GetChoice("newton",solveName,"solver",solveNames,1)) {
			case 1:
				return new SF_Solve(compute,ReactionQTemp,MolQTemp,SegQTemp,LatTemp,MyInput);
				break;
			case 2:
				return new SF_SolveCopel(compute,ReactionQTemp,MolQTemp,SegQTemp,LatTemp,MyInput);
				break;
			case 3:
				Message(literal,"phi-and-U iteration is under construction; better not to proceed....");
				return new SF_SolvePhiU(compute,ReactionQTemp,MolQTemp,SegQTemp,LatTemp,MyInput);
				break;
			case 4:
				return new SF_SolvePikar(compute,ReactionQTemp,MolQTemp,SegQTemp,LatTemp,MyInput);
				break;
			case 5:
				return new SF_SolveCG(compute,ReactionQTemp,MolQTemp,SegQTemp,LatTemp,MyInput);
				break;
			case 6:
				return new SF_Solve_Johan(compute,ReactionQTemp,MolQTemp,SegQTemp,LatTemp,MyInput);
				break;
			default:
				Message(fatal,"Error in SF_System::NewSolve, undefined newton type");
				break;
		}
	}
	return NULL; // never get here
}
void
SF_System::UpdateChi() {
	Matrix chi = SegQ->GetChiMatrix();
	int max=chi.IMAX;
	for (int i=1; i<=max; i++) {
		for (int j=1; j<=max; j++) {
			chi[i][j] *= ROOMTEMPERATURE/temperature;
		}
	}
}
void
SF_System::AddChemPotGraft(Vector profile, Lattice* Lat) const {
	int M = Lat->GetTotalNumLayers();
	for (int i=1; i<=MolQ->GetNumMolecules(); i++) {
		SF_Molecule* Mol = MolQ->GetMolecule(i);
		if (Mol->GetFreedom() == secondGeneration) {
			Vector phi = Mol->GetPhi(constrained);
			Lat->SubtractBoundaries(phi);
			double N = Mol->GetChainLength();
			double mu = MolQ->GetChemicalPotential(Mol,constrained);
			for (int z=1; z<=M; z++) {
				profile[z] += phi[z]*mu/N;
			}
			Lat->RestoreBoundaries(phi);
		} else if (Mol->GetFreedom() == thirdGeneration) {
			Message(implementation,"SF_System::AddChemPotGraft for thirdgeneration");
		} else if (Mol->GetMolStructure()->SomeSegmentsPinned()
				|| Mol->GetMolStructure()->SomeSegmentsGrafted()) {
			Vector phi = Mol->GetPhi(total);
			Lat->SubtractBoundaries(phi);
			double N = Mol->GetChainLength();
			double mu = MolQ->GetChemicalPotential(Mol);


			for (int z=1; z<=M; z++) {
				profile[z] += phi[z]*mu/N;
			}

			Lat->RestoreBoundaries(phi);
		}
	}
}



bool
SF_System::AddChemPotGraftForHelmholtz(Vector profile, Lattice* Lat) const {
	bool apply;
	int NReactionStates=0;
	int imol=0;
	int jstate=0;
	apply=true;
	int M = Lat->GetTotalNumLayers();
	for (int i=1; i<=MolQ->GetNumMolecules(); i++) {
		double theta;

		SF_Molecule* Mol = MolQ->GetMolecule(i);
		if (Mol->GetFreedom() == secondGeneration || Mol->GetFreedom() == thirdGeneration) {
			apply=false;
		} else if (Mol->GetMolStructure()->SomeSegmentsPinned()
				|| Mol->GetMolStructure()->SomeSegmentsGrafted()) {

			int NReactions = ReactionQ->GetNumReactions();
			Vector phi = Mol->GetPhi(total);
			theta = Mol ->GetTheta();
			Lat->SubtractBoundaries(phi);
			double N = Mol->GetChainLength();
			double mu = MolQ->GetNumMolTimesChemPot(i)/theta*N;

			for (int z=1; z<=M; z++) {
				profile[z] += phi[z]*mu/N;
			}
			Lat->RestoreBoundaries(phi);
			SF_MolStructure* Chain;
			Chain = Mol->GetMolStructure();
			SF_MolStructure* Chain2;
			SF_MolSegment* MolSeg;
			SF_MolSegment* MolSeg2;
			SF_MolState* MolState;
			SF_MolState* MolState2;
			SF_Reaction* Reaction;
			SF_State* R1State;
			SF_State* R2State;
			SF_Segment* Segment1;
			SF_Segment* Segment2;
			SF_Molecule* Mol2;
			Vector alpha;
			bool left,right;

			for (int k=1; k<=Chain->GetNumDiffSegments(); k++) {
				MolSeg = Chain->GetDiffSegment(k);
				int NStates =MolSeg->GetNumStates();
				if (NStates>1){
					for (int l=1; l<=NStates; l++) {
						MolState = MolSeg->GetState(l);
						phi = Chain->GetDiffSegment(k)->GetPhi(total);
						alpha = Chain->GetDiffSegment(k)->GetAlpha(l);
						for (int m=1; m<=NReactions; m++){
							Reaction = ReactionQ -> GetReaction(m);
							NReactionStates=Reaction->GetNumDiffStates();
							for (int n=1; n<=NReactionStates; n++) {
								if (n<=NReactionStates/2) {left = true;} else {left = false;}
								R1State=Reaction->GetDiffState(n);
								Segment1 = SegQ->GetBaseSegment(R1State);
								if (R1State->GetName()==MolState->GetName()) {
									for (int o=1; o<=NReactionStates; o++) {
										if (o<=NReactionStates/2) {right = false;} else {right=true;}
										R2State=Reaction->GetDiffState(o);
										Segment2 = SegQ->GetBaseSegment(R2State);
										if (Segment1 != Segment2 && ((left && right) || (!left && !right))) {
											for (int p=1; p<=MolQ->GetNumMolecules(); p++) {
												Mol2=MolQ->GetMolecule(p);
												Chain2 = Mol2->GetMolStructure();
												if (Chain2->GetNumDiffSegments()==1){
													MolSeg2=Chain2->GetDiffSegment(1);
													if (MolSeg2->GetName()==Segment2->GetName()) {
														imol=p;
														int NStates2=MolSeg2->GetNumStates();
														for (int q=1; q<=NStates2; q++){
															MolState2=MolSeg2->GetState(q);
															if (MolState2->GetName()==R2State->GetName()){
																jstate=q;
															}
														}
													}
												}
											}
											if (imol>0&& jstate>0) {
												mu=MolQ->GetChemicalPotential(imol,total,jstate);
												for (int z=1; z<=M; z++) {
													profile[z] -= phi[z]*alpha[z]*mu;
												}
											}
											else {cout << "error in Helmholtz free energy po"<< endl; }
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
	return apply;
}

int
getmemuse() {
#ifndef WIN32
	char filename[25];
	char line[100];
	FILE *fd;

	sprintf(filename, "/proc/%i/status",getpid());
	if (!(fd = fopen(filename, "r"))) {
		return -1;
	}
	do {
		char* x=fgets(line,100,fd); x = x;
		if (!strncmp(line,"VmSize:",7)) {
			int exitvalue = atoll(line+7);
			fclose(fd);
			return exitvalue;
		}
	} while (!feof(fd));
	fclose(fd);
#endif
	return -1;
}

double
SF_System::GetFreeEnergyPo (Lattice* Lat){
	double Fpo_value=0;
	Vector profile;
	if (graftedMolecules == true) {
		if (SegQ->ReactionsPresent()) {
		// ugly, we really want a profile for free energy(po)
		// there's no easy way to obtain it when reactions are present
			Fpo_value = MolQ->GetGrandPotential2();
			for (int i=1; i<=MolQ->GetNumMolecules(); i++) {
				SF_Molecule* Mol = MolQ->GetMolecule(i);
				if (Mol->GetFreedom() == secondGeneration) {
					Fpo_value += MolQ->GetNumMolTimesChemPot(i,constrained);
				} else if (Mol->GetMolStructure()->SomeSegmentsPinned()
						|| Mol->GetMolStructure()->SomeSegmentsGrafted()) {
					Fpo_value += MolQ->GetNumMolTimesChemPot(i);
				}
			}
			return Fpo_value;
		} else { // ok, this is as things should be
			profile = MolQ->GetGrandPotentialProfile();
			AddChemPotGraft(profile,Lat);
			Lat->MultiplyWithLatticeSites(profile);
			int M = Lat->GetTotalNumLayers();
			for (int z=1; z<=M; z++) {
				Fpo_value += profile[z];
			}
			Lat->DivideByLatticeSites(profile);
			return Fpo_value;
		}
	} else return MolQ->GetFreeEnergy();
}

void
SF_System::GetOutput(Output* Out, Lattice* Lat) const {

	Out->PutText("sys",name,"inputfile",MyInput->GetFileName());
	Out->PutText("sys",name,"program version",GetVersion());
	Out->PutText("sys",name,"program date",GetDate());
	Out->PutText("sys",name,"calculation_type","equilibrium");
	Out->PutInt("sys",name,"calculation number", MyInput->GetCurrNumCalculation());
	if (approach == firstOrder) {
		Out->PutText("sys",name,"matrix_approach","first_order");
	}
	if (approach == secondOrder) {
		Out->PutText("sys",name,"matrix_approach","second_order");
	}
	//if (MayerSaupe) {
	//	Out->PutText("sys",name,"MayerSaupe","true");
	//}
	Out->PutBoolean("sys",name,"overflow_protection",overflowProtection);
	Out->PutBoolean("sys",name,"change_chi_with_temperature",changeChiWithTemperature);
	Out->PutReal("sys",name,"temperature",TEMPERATURE);
	Out->PutReal("sys",name,"grand potential",MolQ->GetGrandPotential());
	Out->PutReal("sys",name,"free energy",MolQ->GetFreeEnergy());
	Out->PutReal("sys",name,"grand potential2",MolQ->GetGrandPotential2());



	Vector profile;
	if (graftedMolecules == true) {
		if (SegQ->ReactionsPresent()) {
		// ugly, we really want a profile for free energy(po)
		// there's no easy way to obtain it when reactions are present
			double value = MolQ->GetGrandPotential();
			for (int i=1; i<=MolQ->GetNumMolecules(); i++) {
				SF_Molecule* Mol = MolQ->GetMolecule(i);
				if (Mol->GetFreedom() == secondGeneration) {
					value += MolQ->GetNumMolTimesChemPot(i,constrained);
				} else if (Mol->GetMolStructure()->SomeSegmentsPinned()
						|| Mol->GetMolStructure()->SomeSegmentsGrafted()) {
					value += MolQ->GetNumMolTimesChemPot(i);
				}
			}
			Out->PutReal("sys",name,"free energy(po)",value);
		} else { // ok, this is as things should be
			profile = MolQ->GetGrandPotentialProfile();
			AddChemPotGraft(profile,Lat);
			Lat->MultiplyWithLatticeSites(profile);
			int M = Lat->GetTotalNumLayers();
			double value = 0;
			for (int z=1; z<=M; z++) {
				value += profile[z];
			}
			Lat->DivideByLatticeSites(profile);
			Out->PutReal("sys",name,"free energy(po)",value);
		}
	}
	Vector profile2;
	bool correct=false;

	if (graftedMolecules == true) {
		profile2 = MolQ->GetGrandPotentialProfile();
		correct=AddChemPotGraftForHelmholtz(profile2,Lat);
		Lat->MultiplyWithLatticeSites(profile2);
		int M = Lat->GetTotalNumLayers();
		double value = 0;
		for (int z=1; z<=M; z++) {
			value += profile2[z];
		}
		Lat->DivideByLatticeSites(profile2);
		if (correct) Out->PutReal("sys",name,"Helmholtz energy(po)",value);
	}

	Out->PutReal("sys",name,"excess free energy",MolQ->GetExcessFreeEnergy());

	Out->PutReal("sys",name,"entropy",MolQ->GetEntropy());

	Out->PutReal("sys",name,"contact int.",MolQ->GetContactInteractions());
	Out->PutReal("sys",name,"max_mem_use",getmemuse()/1024.0);
	if (ReactionQ->GetNumReactions() > 0) {
		Out->PutReal("sys",name,"internal free energy",MolQ->GetInternalFreeEnergy());
	}
	if (SegQ->Charged()) {
		Out->PutReal("sys",name,"electric int",MolQ->GetElectricInteractions());
	}
	Out->PutProfile("sys",name,"free energy density",MolQ->GetFreeEnergyProfile());
	Out->PutProfile("sys",name,"grand potential density",MolQ->GetGrandPotentialProfile());

	if (graftedMolecules == true && !SegQ->ReactionsPresent()) {
		Lat->DivideByLatticeSites(profile);
		Out->PutProfile("sys",name,"free energy(po) density",profile);

	}
	if (graftedMolecules == true ) {
		Lat->DivideByLatticeSites(profile2);
		if (correct) Out->PutProfile("sys",name,"Helmholtz energy(po) density",profile2);
	}

	Out->PutProfile("sys",name,"excess free energy density",MolQ->GetExcessFreeEnergyProfile());
	Out->PutProfile("sys",name,"entropy density",MolQ->GetEntropyProfile());
	Out->PutProfile("sys",name,"contact int. density",MolQ->GetContactInteractionsProfile());
	if (ReactionQ->GetNumReactions() > 0) {
		Out->PutProfile("sys",name,"internal free energy density",MolQ->GetInternalFreeEnergyProfile());
	}
	if (SegQ->Charged()) {
		Out->PutProfile("sys",name,"electric int. density",MolQ->GetElectricInteractionsProfile());
	}
	if (SegQ->Charged()) {
		Vector psi = SegQ->GetElectricPotential();
		Lat->SetBoundaries(psi);
		Out->PutProfile("sys",name,"potential",psi);

		psi = SegQ->GetCharge();
		int M = Lat->GetTotalNumLayers();
		Vector charge_density(1,M);
		for (int z=1; z<=M; z++) charge_density[z]=psi[z];
		Lat->DivideByLatticeSites(charge_density);
		Out->PutProfile("sys",name,"charge density",charge_density);
		Out->PutProfile("sys",name,"epsilon",SegQ->GetAverageEpsilon());
	}
	if (MyInput->ValueSet("sys",name,"extra_output")) {
		GetExtraOutput(Out,Lat);
	}

	Lat->GetOutput(Out);
	MolQ->GetOutput(Out);
	SegQ->GetOutput(Out);
	ReactionQ->GetOutput(Out);
	Solve->GetOutput(Out);
	Solve->WriteInitialGuessToFile();
}
void
SF_System::GetExtraOutput(Output* Out, Lattice* Lat) const{
	Array<Text> choices(1,44);
	choices[1] = "pressure";
	choices[2] = "joanne";
	choices[3] = "second_gen";
	choices[4] = "article";
	choices[5] = "hairy";
	choices[6] = "MSc_thesis";
	choices[7] = "poster";
	choices[8] = "vesicle";
	choices[9] = "isolated";
	choices[10] = "spinodal";
	choices[11] = "emulsion";
	choices[12] = "two_phase";
	choices[13] = "two_phase_cyl";
	choices[14] = "two_phase_sph";
	choices[15] = "depletion";
	choices[16] = "pi_a";
	choices[17] = "bend";
	choices[18] = "multistate";
	choices[19] = "cecilia";
	choices[20] = "moduli";
	choices[21] = "pore";
	choices[22] = "Emilia";
	choices[23] = "Sergio";
	choices[24] = "DSM";
	choices[25] = "Egorov";
	choices[26] = "Brush3G";
	choices[27] = "Escape";
	choices[28] = "Moments2G";
	choices[29] = "Moments3G";
	choices[30] = "Probesize";
	choices[31] = "Soumi";
	choices[32] = "Loreal";
	choices[33] = "Frans";
	choices[34] = "GetPhiZ";
	choices[35] = "Dendron_overlap";
	choices[36] = "Helene";
	choices[37] = "Katya";
	choices[38] = "Sabine";
	choices[39] = "Johan";
	choices[40] = "Boris";
	choices[41] = "zeta_potential";
	choices[42] = "coexisting_bilayers";
	choices[43] = "interfacial_width";
	choices[44] = "line_tension";
	int choice = MyInput->GetChoice("sys",name,"extra_output",choices);
	switch(choice) {
		case 1 :
			GetOutputPressure(Out,Lat);
			break;
		case 2 :
			GetOutputJoanne(Out,Lat);
			break;
		case 3 :
			GetOutputSecondGen(Out,Lat);
			break;
		case 4 :
			Message(literal,"Sorry, 'extra_output : article' is not fully implemented");
			break;
		case 5 :
			GetOutputHairy(Out,Lat);
			break;
		case 6 :
			Message(literal,"Sorry, 'extra_output : MSc_thesis' is not fully implemented");
			break;
		case 7 :
			Message(literal,"Sorry, 'extra_output : poster' is not fully implemented");
			break;
		case 8 :
			GetOutputVesicle(Out,Lat);
			break;
		case 9 :
			GetOutputIsolated(Out,Lat);
			break;
		case 10 :
			GetOutputSpinodal(Out,Lat);
			break;
		case 11 :
			GetOutputEmulsion(Out,Lat);
			break;
		case 12 :
			GetOutputTwoPhase(Out,"flat",Lat);
			break;
		case 13 :
			GetOutputTwoPhase(Out,"cylinder",Lat);
			break;
		case 14 :
			GetOutputTwoPhase(Out,"sphere",Lat);
			break;
		case 15 :
			GetOutputDepletion(Out,Lat);
			break;
		case 16 :
			GetOutputPiA(Out,Lat);
			break;
		case 17 :
			GetOutputBend(Out,Lat);
			break;
		case 18 :
			GetOutputMultiState(Out,Lat);
			break;
		case 19 :
			GetCeciliaOutput(Out,Lat);
			break;
		case 20 :
			GetModuliOutput(Out,Lat);
			break;
		case 21 :
			GetPoreOutput(Out,Lat);
			break;
		case 22 :
			GetEmiliaOutput(Out,Lat);
			break;
		case 23 :
			GetSergioOutput(Out,Lat);
			break;
		case 24 :
			GetDSMOutput(Out,Lat); GetSergioOutput(Out,Lat);
			break;
		case 25 :
			GetEgorovOutput(Out,Lat);
			break;
		case 26 :
			GetBrush3GOutput(Out,Lat);
			break;
		case 27 :
			GetEscapeOutput(Out,Lat);
			break;
		case 28 :
			GetMoments2GOutput(Out,Lat);
			break;
		case 29 :
			GetMoments3GOutput(Out,Lat);
			break;
		case 30 :
			GetProbeSizeOutput(Out,Lat);
			break;
		case 31 :
			GetSoumiOutput(Out,Lat);
			break;
		case 32 :
			GetLorealOutput(Out,Lat);
			break;
		case 33 :
			GetFransOutput(Out,Lat);
			break;
		case 34 :
			GetPhiZ(Out,Lat);
			break;
		case 35 :
			GetDendronOverlap(Out,Lat);
			break;
		case 36 :
			GrandPotentialLeftRight(Out,Lat);
			break;
		case 37 :
			Katya(Out,Lat);
			break;
		case 38 :
			Sabine(Out,Lat);
			break;
		case 39 :
			Johan(Out,Lat);
			break;
		case 40 :
			Boris(Out,Lat);
			break;
		case 41 :
			Zeta(Out,Lat);
			break;
		case 42 :
			CoExistBilayers(Out,Lat);
			break;
		case 43 :
			InterfacialWidth(Out,Lat);
			break;
		case 44 :
			LineTension(Out,Lat);
			break;
		default :
			Message(fatal,"Error in SF_System::GetExtraOutput");
			break;
	}
}


void
SF_System:: LineTension(Output* Out,Lattice* Lat) const{
	Text dim;
	Text probeName = MyInput->GetText("sys",name,"probe_chain");
	if (!MolQ->MoleculeDefined(probeName)) {Message(fatal,"sys : " + name + " : probe_chain : " + probeName + " is not defined");}
	SF_Molecule* Mol = MolQ->GetMolecule(probeName);
	int n_layers_x, n_layers_y;
	Vector phi = Mol->GetPhi(total);
	Vector omega = MolQ->GetGrandPotentialProfile();
	if (Lat->GetNumGradients()!=2) Message(fatal,"The line tension module expected two gradient calculations. Calculations stopped");
	n_layers_x = Lat->GetNumLayers(1)-2;
	n_layers_y = Lat->GetNumLayers(2)-2;
	double gamma=0;
	int x=1;
	int i;
	int jx=n_layers_y+2;
	for (int y=1; y<=n_layers_y; y++) {
		gamma += omega[jx*x+y+1];
	}
	Out->PutReal("sys",name,"gamma(x=1)",gamma);
	Out->PutReal("sys",name,"gammaXA",gamma*n_layers_x);
	Vector pos(1,n_layers_y);
	Vector tension(1,n_layers_y);
	double low,high,theta,GP;
	for (int y=1; y<=n_layers_y; y++){
		low =1; high=0; theta=0; GP=0;
		for (int x=1; x<=n_layers_x; x++) {
			i=x*jx+y+1;
			if (phi[i] < low) low = phi[i];
			if (phi[i] > high) high = phi[i];
			theta +=phi[i];
			GP += omega[i];
		}
		theta -=n_layers_x*low;
		if ((high-low)>0) pos[y]=theta/(high-low);


		dim = Blanks(9);
		dim.Putint(y);
		dim = Copy(dim.Strip());
		dim = Copy(dim.Frontstrip());
		Out ->PutReal("sys",name,"pos_x["+ dim + "]",pos[y]);

		tension[y]=GP-omega[1*jx+y+1]*pos[y]-omega[n_layers_x*jx+y+1]*(n_layers_x-pos[y]);
	}
	double Tension=0;
	for (int y=1; y<=n_layers_y; y++) {
		Tension += tension[y];
	}
	Tension=Tension-n_layers_y*tension[n_layers_y/2];
	int yup = n_layers_y/2+1;
	int ymean = n_layers_y/2;
	int ydown = n_layers_y/2-1;
	double L = n_layers_y*pow(1+pow((pos[yup]-pos[ydown])/2,2),0.5);
	Out ->PutReal("sys",name,"dx/dy",(pos[yup]-pos[ydown])/2);
	Out ->PutReal("sys",name,"d2x/dy2",pos[yup]+pos[ydown]-2*pos[ymean]);
	Out ->PutReal("sys",name,"L",L);
	Out ->PutReal("sys",name,"tension",Tension/2);

	if ((pos[yup]-pos[ydown])!=0) {
		double dpos = (pos[yup]-pos[ydown]); if (dpos < 0) dpos *=-1;
		double angle = atan (2/dpos) * 180 / PI;
		Out ->PutReal("sys",name,"contact angle",angle);
	}	else Out ->PutReal("sys",name,"contact angle",90.0);

}


void
SF_System:: InterfacialWidth(Output* Out,Lattice* Lat) const{

	Text probeName = MyInput->GetText("sys",name,"Gibbs_molecule");
	if (!MolQ->MoleculeDefined(probeName)) {Message(fatal,"sys : " + name + " : Gibbs_molecule : " + probeName + " is molecule that is not defined");}
	SF_Molecule* Mol = MolQ->GetMolecule(probeName);

	Vector phi = Mol->GetPhi(total);
	int M = Lat->GetTotalNumLayers();
	double phi_low=1;
	double phi_high=0;
	double gradient_max =0;
	double gradient = 0;
	double width=0;
	for (int z=3; z<M; z++) {
		if (phi[z] < phi_low) phi_low=phi[z];
		if (phi[z] > phi_high) phi_high = phi[z];
		gradient = (phi[z]-phi[z-1])*(phi[z]-phi[z-1]);
		if (gradient > gradient_max) gradient_max=gradient;
	}
	if (gradient_max > 0) width =(phi_high-phi_low)/sqrt(gradient_max);
	Out->PutReal("sys",name,"width",width);
	Lat->SubtractBoundaries(phi);
	double theta_exc = Mol->ComputeThetaExcess();
	double Dphi=phi_high-phi_low;
	double Gvolume=0;
	if (Dphi>0) Gvolume=theta_exc/Dphi; if (Gvolume < 0) Gvolume=-1.0*Gvolume;
	Out->PutReal("sys",name,"VGibbs",Gvolume);
	Text latName;
	MyInput->GetNumNames("lat",1,1);
	Array<Text> latNames = MyInput->GetNames("lat");
	latName = latNames[1];
	Array<Text> geometries(1,3);
	geometries[1] = "flat";
	geometries[2] = "cylindrical";
	geometries[3] = "spherical";
	int geom = MyInput->GetChoice("lat",latName,"geometry",geometries,1);
	double RGibbs=0;
	double offset = Lat->GetLayerAdjustment();
	double offsetLayer1 = MyInput->GetReal("lat",latName,"offset_first_layer",0,DBL_MAX,0);
	Vector GPp = MolQ->GetGrandPotentialProfile();
	double GP=MolQ->GetGrandPotential();
	double RSOT=0;

	if (Gvolume > 0) {
		if (geom==1) {
			RGibbs=Gvolume;
		}
		if (geom==2) {
			Gvolume += PI*pow(offsetLayer1,2.0);
			GP +=GPp[2]*PI*pow(offsetLayer1,2.0);
			RGibbs=pow(Gvolume/PI,1.0/2.0);
			double X=GP/GPp[2]/PI; if (X<0) X=-1.0*X;
			RSOT=pow(X,0.5);
		}
		if (geom==3) {
			GP +=GPp[2]*PI*pow(offsetLayer1,3.0);
			Gvolume += 4.0/3.0*PI*pow(offsetLayer1,3.0);
			RGibbs=pow(Gvolume*3.0/(PI*4.0),1.0/3.0);
			double X=GP/GPp[2]*3.0/2.0/PI; if (X<0) X=-1.0*X;
			RSOT=pow(X,1.0/3.0);
		}
	}
	Out->PutReal("sys",name,"RGibbs",RGibbs);
	Out->PutReal("sys",name,"RSOT",RSOT);
	Out->PutReal("sys",name,"shift",offset-1.0);
	Lat->RestoreBoundaries(phi);
}

void
SF_System:: CoExistBilayers(Output* Out,Lattice* Lat) const{
	int M = Lat->GetTotalNumLayers();
	double theta_exc;
	double tension;
	Vector GP = MolQ->GetGrandPotentialProfile();
	Lat->SubtractBoundaries(GP);

	tension=0;
	for (int z=1; z<=M/2; z++) {
		tension +=GP[z];
	}
	Out->PutReal("sys",name,"tension_bilayer_z=1",tension);
	tension=0;
	for (int z=M/2+1; z<=M; z++) {
		tension +=GP[z];
	}
	Out->PutReal("sys",name,"tension_bilayer_z=M",tension);

	for (int i=1; i<=MolQ->GetNumMolecules(); i++) {
		SF_Molecule* Mol = MolQ->GetMolecule(i);
		Text molName = Mol->GetName();
		Vector phi = Mol->GetPhi(total);
		Lat->SubtractBoundaries(phi);
		double phibulk = Mol->GetPhiBulk();
		theta_exc=0;
		for (int z=1; z<=M/2; z++) {
			theta_exc += phi[z]-phibulk;
		}
		Out->PutReal("mol",molName,"Theta_excess_z=1",theta_exc);

		theta_exc=0;
		for (int z=M/2+1; z<=M; z++) {
			theta_exc += phi[z]-phibulk;
		}
		Out->PutReal("mol",molName,"Theta_excess_z=M",theta_exc);
	}
}

void
SF_System:: Zeta(Output* Out,Lattice* Lat) const{

	Text probeName = MyInput->GetText("sys",name,"probe_monomer");
	if (!SegQ->SegmentDefined(probeName)) {Message(fatal,"sys : " + name + " : probe_monomer : " + probeName + " is a monomer that is not defined");}
	SF_Segment* Mon = SegQ->GetSegment(probeName);

	Vector phi = Mon->GetPhi();
	Lat->MultiplyWithLatticeSites(phi);
	int M = Lat->GetTotalNumLayers();


	LatticeRange* LayerZero=NULL;
	if (Lat->GetNumGradients() == 1) {
		LayerZero = new LatticeRange1D(1,1,M);
	}
	double R = pow(Lat->MomentUnweighted(phi,2,LayerZero,0)/Lat->MomentUnweighted(phi,0,LayerZero,0),0);
	Vector psi = SegQ->GetElectricPotential();
	int r = R;
	int r1= R+1;
	Out->PutReal("sys",name,"plane_of_shear",R);
	Out->PutReal("sys",name,"zeta",psi[r]+(psi[r1]-psi[r])*(R-r));
}

void
SF_System:: Boris(Output* Out,Lattice* Lat) const{

	Text probeName = MyInput->GetText("sys",name,"probe_monomer");
	Text chainName = MyInput->GetText("sys",name,"probe_chain");
	if (!SegQ->SegmentDefined(probeName)) {Message(fatal,"sys : " + name + " : probe_monomer : " + probeName + " is a monomer that is not defined");}
	if (!MolQ->MoleculeDefined(chainName)) {Message(fatal,"sys : " + name + " : probe_chain : " + chainName + " is a molecule that is not defined");}
	SF_Segment* Mon = SegQ->GetSegment(probeName);
	SF_Molecule* Mol = MolQ->GetMolecule(chainName);

	Vector phi = Mon->GetPhi();
	Lat->MultiplyWithLatticeSites(phi);
	int n_layers = Lat->GetNumLayers(1)-2;
	Vector n(1,n_layers+2); for (int z=1; z<=n_layers+2; z++) n[z] = phi[z];
	Vector c(1,n_layers+2);
	Vector m(1,n_layers+2);
	c[2] = phi[2];
	for (int z=3; z<=n_layers+1; z++) c[z] = c[z-1]+phi[z];
	Out->PutProfile("sys",name,"cummulative_"+probeName,c);
	Out->PutProfile("sys",name,"n_"+probeName,n);
	Lat->DivideByLatticeSites(phi);
	Vector phi_ = Mol->GetPhi(total);
	for (int z=1; z<=n_layers+2; z++) m[z] = phi_[z];
	Lat->MultiplyWithLatticeSites(m);
	Out->PutProfile("sys",name,"n_"+chainName,m);
	Vector charge = SegQ->GetCharge();
	Out->PutProfile("sys",name,"charge",charge);
}

void
SF_System:: Johan(Output* Out,Lattice* Lat) const{
	Solve->ReComputePhi(true);
	if (!MyInput->ValueSet("sys",name,"probe_chain")) Message(fatal,"sys : probe_chain not defined");
	Text chainName = MyInput->GetText("sys",name,"probe_chain");
	if (!MolQ->MoleculeDefined(chainName)) {Message(fatal,"sys : " + name + " : probe_chain : " + chainName + " is a molecule that is not defined");}
	SF_Molecule* Mol = MolQ->GetMolecule(chainName);
	Vector phi = Mol->GetPhi(total);
	int n_layers_x = Lat->GetNumLayers(1)-2;
	int n_layers_y = Lat->GetNumLayers(2)-2;
	int n_layers_z = Lat->GetNumLayers(3)-2;
	int jx = (n_layers_y+2)*(n_layers_z+2);
	int jy = n_layers_z+2;
	double theta=0;
	for (int x=1; x<=n_layers_x; x++)
	for (int y=1; y<=n_layers_y; y++)
	for (int z=1; z<=n_layers_z; z++) theta+=phi[jx*x+jy*y+z];
	double Length = Mol->GetChainLength();
	cout << "number of loops " << theta/Length << endl;

}

void
SF_System:: Sabine(Output* Out,Lattice* Lat) const{

	if (Lat->GetNumGradients()==3) {
		Text dim;
		if (!MyInput->ValueSet("sys",name,"probe_monomer")) Message(fatal,"sys : probe_moniner not defined");
		Text probeName = MyInput->GetText("sys",name,"probe_monomer");
		if (!SegQ->SegmentDefined(probeName)) {Message(fatal,"sys : " + name + " : probe_monomer : " + probeName + " is a monomer that is not defined");}
		SF_Segment* Mon = SegQ->GetSegment(probeName);
		Vector phi = Mon->GetPhi();

		int n_layers_x = Lat->GetNumLayers(1)-2;
		int n_layers_y = Lat->GetNumLayers(2)-2;
		int n_layers_z = Lat->GetNumLayers(3)-2;
		int jx = (n_layers_y+2)*(n_layers_z+2);
		int jy = n_layers_z+2;
		Vector alpha(1,n_layers_x);
		Vector y_old(1,n_layers_x);
		Vector y_new(1,n_layers_x);

		double Phi,Phi_old;
		double Mean_y,d_y;
		double Mean_alpha,d_alpha;
		double yy;
		double aa;
		double max_y;
		for (int x=1; x<=n_layers_x; x++) {y_old[x]=0; y_new[x]=0; alpha[x]=0; }
		for (int z=1; z<=n_layers_z; z++) {
			Mean_y=0;
			Mean_alpha=0;
			for (int x=1; x<=n_layers_x; x++) {
				Phi=0; y_old[x]=y_new[x];
				for (int y=1; y<=n_layers_y; y++) {
					Phi_old = Phi;
					Phi = phi[jx*x+jy*y+z];
					if (Phi_old < 0.5 && Phi > 0.5) {
						y_new[x]=y-1+(0.5-Phi_old)/(Phi-Phi_old);
						Mean_y +=y_new[x]; y=n_layers_y;
						if (y_old[x] > 0) {
							if (y_new[x] == y_old[x]) alpha[x]=90.0; else
							alpha[x] = atan(1/(y_new[x]-y_old[x]))*180.0/PI;
							if (alpha[x] > 0) alpha[x]=180-alpha[x];
							if (alpha[x] < 0) alpha[x] =-alpha[x];
							Mean_alpha +=alpha[x];
						}
					}
				}
			}
			Mean_y = Mean_y/n_layers_x;
			Mean_alpha = Mean_alpha/n_layers_x;
			yy=0; aa=0; max_y=0;
			for (int x=1; x<=n_layers_x; x++){
				if (y_new[x] > max_y) max_y=y_new[x];
				yy += pow(y_new[x]-Mean_y,2);
				aa += pow(alpha[x]-Mean_alpha,2);
			}
			d_y=0; d_alpha=0;
			if (yy>0) d_y=sqrt(yy/n_layers_x);
			if (aa>0) d_alpha=sqrt(aa/n_layers_x);
			dim = Blanks(9);
			dim.Putint(z);
			dim = Copy(dim.Strip());
			dim = Copy(dim.Frontstrip());
			Out->PutReal("mon",probeName,"<alpha>[" + dim + "]",Mean_alpha);
			Out->PutReal("mon",probeName,"<y>[" + dim + "]",Mean_y);
			Out->PutReal("mon",probeName,"d_alpha[" + dim + "]",d_alpha);
			Out->PutReal("mon",probeName,"d_y[" + dim + "]",d_y);
			Out->PutReal("mon",probeName,"D_y(max)[" + dim + "]",max_y-Mean_y);
		}
	}
}

void
SF_System:: Katya(Output* Out,Lattice* Lat) const{
	int numSeg = SegQ->GetNumSegments();
	for (int i=1; i<=numSeg; i++) {
		SF_Segment* Seg = SegQ->GetSegment(i);
		Vector phi = Seg->GetPhi();

		if (Lat->GetNumGradients()==2) {
			int n_layers_x = Lat->GetNumLayers(1)-2;
			int n_layers_y = Lat->GetNumLayers(2)-2;
			int jx = n_layers_y+2;
			double teller=0, noemer=0;
			Lat->MultiplyWithLatticeSites(phi);
			for (int x=1; x<= n_layers_x; x++)
			for (int y=1; y<=n_layers_y; y++) {
				teller += x*phi[jx*x+y];
				noemer += phi[jx*x+y];
			}
			Lat->DivideByLatticeSites(phi);
			Out->PutReal("mon",Seg->GetName(),"<r>",teller/noemer);
		} else 	if (Lat->GetNumGradients()==3) {
			int n_layers_x = Lat->GetNumLayers(1)-2;
			int n_layers_y = Lat->GetNumLayers(2)-2;
			int n_layers_z = Lat->GetNumLayers(3)-2;
			int jx = (n_layers_y+2)*(n_layers_z+2);
			int jy = n_layers_z+2;

			double teller=0, noemer=0;
			for (int x=1; x<=n_layers_x; x++)
			for (int y=1; y<=n_layers_y; y++)
			for (int z=1; z<=n_layers_z; z++) {
				teller += pow((x-(n_layers_x+1)/2)*(x-(n_layers_x+1)/2)+(y-(n_layers_y+1)/2)*(y-(n_layers_y+1)/2),0.5)*phi[jx*x+jy*y+z];
				noemer += phi[jx*x+jy*y+z];
			}
			Out->PutReal("mon",Seg->GetName(),"<r>", teller/noemer);
		}
	}
}

void
SF_System:: GrandPotentialLeftRight(Output* Out,Lattice* Lat) const{
//	Vector w = MolQ->GetGrandPotentialProfile();
//	int n_layers = Lat->GetNumLayers(1)-2;
//	Vector Wleft(1,n_layers+2);
//	Vector Wright(1,n_layers+2);
//	for (int z=2; z<=n_layers+1; z++) {
//		Wleft[z] =Wleft[z-1]+w[z];
//	}
//	Out->PutProfile("sys",name,"GPLeft",Wleft);

	GetEmiliaOutput(Out,Lat);
	if (!MyInput->ValueSet("sys",name,"probe_chain")) Message(fatal,"sys : probe_chain not defined");
	Text chainName = MyInput->GetText("sys",name,"probe_chain");
	SF_Molecule* MolOil = MolQ->GetMolecule(chainName);
	Vector phi_oil = MolOil->GetPhi(total);
	Text GibbsMolName = MyInput->GetText("sys",name,"Gibbs_molecule"); //phi-water
	if (!MolQ->MoleculeDefined(GibbsMolName)) {Message(fatal,"'sys : " + name + " : Gibbs_molecule : " + GibbsMolName + " is a molecule that is not defined");}
	SF_Molecule* GibbsMol;
	GibbsMol = MolQ->GetMolecule(GibbsMolName);  //phi-oil
	Vector phi_water = GibbsMol->GetPhi(total);
	int n_layers = Lat->GetNumLayers(1)-2;

	int zlow=1;
	int zhigh=n_layers;
	double philow=phi_oil[1];
	double phiHigh=phi_water[n_layers];
	int z=1;
	while (phi_oil[z]>philow/2 && z <= n_layers) z++;
	zlow = z;
	z=n_layers;
	while (phi_water[z]>phiHigh/2 && z >= 1) z--;
	zhigh = z;
	int zhalf=zlow + (zhigh-zlow)/2;
	Out->PutInt("lat",Lat->GetName(),"z_CC",zhalf);

	for (int i=1; i<=MolQ->GetNumMolecules(); i++) {
		SF_Molecule* Mol = MolQ->GetMolecule(i);
		Vector phi = Mol->GetPhi(total);
		Out->PutReal("mol",Mol->GetName(),"Phi_CC",phi[zhalf+1]);
	}
}

void
SF_System::GetDendronOverlap(Output* Out,Lattice* Lat) const{
	if (Lat->GetNumGradients()==1) {
		int n_layers = Lat->GetNumLayers(1)-2;
		double d = n_layers/2.0+0.5;
		if (!MyInput->ValueSet("sys",name,"probe_chain")) Message(fatal,"sys : probe_chain not defined");
		Text chainName = MyInput->GetText("sys",name,"probe_chain");
		if (!MolQ->MoleculeDefined(chainName)) {Message(fatal,"sys : " + name + " : probe_chain : " + chainName + " is a molecule that is not defined");}
		SF_Molecule* MolPol = MolQ->GetMolecule(chainName);
		Vector phi = MolPol->GetPhi(total);
		double zPhi=0;
		double zzPhi=0;
		double Phi=0;
		for (int z=n_layers/2+1; z<=n_layers; z++) {
			zPhi += phi[z]*(z-d);
			zzPhi += phi[z]*(z-d)*(z-d);
			Phi +=phi[z];
		}
		//if (n_layers%2==1) Phi-=0.5*phi[n_layers/2+1];
		if (Phi > 0) Out->PutReal("mol",chainName,"L0",Phi); else Out->PutReal("mol",chainName,"L0",0);
		if (Phi > 0) Out->PutReal("mol",chainName,"L1",zPhi/Phi); else Out->PutReal("mol",chainName,"L1",0);
		if (Phi > 0) Out->PutReal("mol",chainName,"L2",sqrt(zzPhi/Phi)); else Out->PutReal("mol",chainName,"L2",0);

		zPhi =0;
		zzPhi =0;
		Phi =0;

		double THETA=0;
		for (int z=1; z<=n_layers; z++) {THETA +=phi[z];}
		int theta = THETA;
		if (THETA-1.0*theta > 0.5) THETA++;
		d=theta+0.5;

		for (int z=theta+1; z<=n_layers; z++) {
			zPhi += phi[z]*(z-d);
			zzPhi += phi[z]*(z-d)*(z-d);
			Phi +=phi[z];
		}
		if (Phi > 0) Out->PutReal("mol",chainName,"LM0",Phi); else Out->PutReal("mol",chainName,"L0",0);
		if (Phi > 0) Out->PutReal("mol",chainName,"LM1",zPhi/Phi); else Out->PutReal("mol",chainName,"L1",0);
		if (Phi > 0) Out->PutReal("mol",chainName,"LM2",sqrt(zzPhi/Phi)); else Out->PutReal("mol",chainName,"L2",0);

		int z_min_half=0;
		double phi_plus_half=0;
		double phi_min_half=0;
		double z_half=0;

		for (int z=1; z<=n_layers; z++) {
			if (phi[z]>0.5) {phi_plus_half = phi[z]; z_min_half=z; phi_min_half=phi[z+1];}
		}
		z_half = z_min_half+ (phi_plus_half-0.5)/(phi_plus_half-phi_min_half);
		zPhi =0;
		zzPhi =0;
		Phi =0;
		for (int z=z_half+1; z<=n_layers; z++) {
			zPhi += phi[z]*(z-z_half);
			zzPhi += phi[z]*(z-z_half)*(z-z_half);
			Phi +=phi[z];
		}
		if (Phi > 0) Out->PutReal("mol",chainName,"LH0",Phi); else Out->PutReal("mol",chainName,"L0",0);
		if (Phi > 0) Out->PutReal("mol",chainName,"LH1",zPhi/Phi); else Out->PutReal("mol",chainName,"L1",0);
		if (Phi > 0) Out->PutReal("mol",chainName,"LH2",sqrt(zzPhi/Phi)); else Out->PutReal("mol",chainName,"L2",0);

		if (MyInput->ValueSet("sys",name,"probe_monomer")) {
			Text monName = MyInput->GetText("sys",name,"probe_monomer");
			if (!SegQ->SegmentDefined(monName)) {Message(fatal,"sys : " + name + " : probe_monomer : " + monName + " is a monomer that is not defined");}
			SF_Segment* Mon = SegQ->GetSegment(monName);
			Vector phim = Mon->GetPhi();

			double MIN =1.0;
			double MAX =-1.0;
			bool maxfound = false;
			bool minfound = false;
			int zmin=0;

			for (int z=1; z<=n_layers; z++) {
				if (!maxfound) {
					if (phim[z] > MAX) {MAX = phim[z];} else {maxfound=true; }
				} else {
					if (!minfound) {if (phim[z] < MIN) {MIN = phim[z];} else {minfound=true; zmin =z;}}
				}
			}
			double theta1=0.0;
			double theta2=0.0;
			for (int z=1; z<=n_layers; z++) {
				if (z < zmin) theta1+=phim[z];
				theta2+=phim[z];
			}
			Out->PutReal("mon",monName,"f1",theta1/theta2);
		}

	}
}


void
SF_System::GetPhiZ(Output* Out,Lattice* Lat) const{

	if (Lat->GetNumGradients()==3) {
		Text dim;
		if (!MyInput->ValueSet("sys",name,"probe_chain")) Message(fatal,"sys : probe_chain not defined");
		Text chainName = MyInput->GetText("sys",name,"probe_chain");
		if (!MolQ->MoleculeDefined(chainName)) {Message(fatal,"sys : " + name + " : probe_chain : " + chainName + " is a molecule that is not defined");}
		SF_Molecule* MolPol = MolQ->GetMolecule(chainName);

		int n_layers_x = Lat->GetNumLayers(1)-2;
		int n_layers_y = Lat->GetNumLayers(2)-2;
		int n_layers_z = Lat->GetNumLayers(3)-2;
		int jx = (n_layers_y+2)*(n_layers_z+2);
		int jy = n_layers_z+2;
		//int jz = 1;
		double rho;
		Vector phi = MolPol->GetPhi(total);
		for (int z=1; z<=n_layers_z; z++) {
			rho=0;
			for (int x=1; x<=n_layers_x; x++) {
				for (int y=1; y<=n_layers_y; y++) {
					rho +=phi[jx*x+jy*y+z];
				}
			}
			dim = Blanks(9);
			dim.Putint(z);
			dim = Copy(dim.Strip());
			dim = Copy(dim.Frontstrip());
			Out->PutReal("mol",chainName,"rho["+ dim + "]",rho);
		}
	} else {if (Lat->GetNumGradients()==2) {
		Text dim;
		if (!MyInput->ValueSet("sys",name,"probe_chain")) Message(fatal,"sys : probe_chain not defined");
		Text chainName = MyInput->GetText("sys",name,"probe_chain");
		if (!MolQ->MoleculeDefined(chainName)) {Message(fatal,"sys : " + name + " : probe_chain : " + chainName + " is a molecule that is not defined");}
		SF_Molecule* MolPol = MolQ->GetMolecule(chainName);

		int n_layers_x = Lat->GetNumLayers(1)-2;
		int n_layers_y = Lat->GetNumLayers(2)-2;
		int jx = (n_layers_y+2);
		//int jy = 1;
		//int jz = 1;
		double rho;
		Vector phi = MolPol->GetPhi(total);Lat->MultiplyWithLatticeSites(phi);

		for (int y=1; y<=n_layers_y; y++) {
			rho=0;
			for (int x=1; x<=n_layers_x; x++) {
				rho +=phi[jx*x+y];
			}
			dim = Blanks(9);
			dim.Putint(y);
			dim = Copy(dim.Strip());
			dim = Copy(dim.Frontstrip());
			Out->PutReal("mol",chainName,"rho["+ dim + "]",rho);
		}
		Lat->DivideByLatticeSites(phi);
	}}
}

void
SF_System::GetFransOutput(Output* Out,Lattice* Lat) const{
	if (!MyInput->ValueSet("sys",name,"probe_monomer")) Message(fatal,"sys : probe_mon not defined");
	Text probeName = MyInput->GetText("sys",name,"probe_monomer");
	if (!SegQ->SegmentDefined(probeName)) {Message(fatal,"sys : " + name + " : probe_monomer : " + probeName + " is a monomer that is not defined");}

	if (!MyInput->ValueSet("sys",name,"pinned_monomer")) Message(fatal,"sys : pinned_mon not defined");
	Text pinnedName = MyInput->GetText("sys",name,"pinned_monomer");
	if (!SegQ->SegmentDefined(pinnedName)) {Message(fatal,"sys : " + name + " : pinned_monomer : " + probeName + " is a monomer that is not defined");}


	SF_Segment* Mon = SegQ->GetSegment(probeName);
	SF_Segment* MonPin = SegQ->GetSegment(pinnedName);


	if (Lat->GetNumGradients()==2) {
	Vector phi = MonPin->GetPhi();Lat->MultiplyWithLatticeSites(phi);
	Vector G = Mon->GetSWF(); Lat->SetBoundaries(G);
		int n_layers_y = Lat->GetNumLayers(2)-2;
		int jx = n_layers_y+2;
		double force_y=0;
		double a = 0, b=0, c=0;
		if (G[jx]>0) b = log(G[jx]); else b=0;
		if (G[jx+1]>0) c = log(G[jx+1]); else c=0;
		for (int y=1; y<=n_layers_y; y++) {
			a = b; b = c;
			if (G[jx+y+1]>0) c= log(G[jx+y+1]); else c=0;
			force_y += phi[jx+y]*(c-a)/2;
		}
		Out->PutReal("sys",name,"fy",force_y);
		Lat->DivideByLatticeSites(phi);
	}
	if (Lat->GetNumGradients()==3) {
	Vector phi = MonPin->GetPhi();
	Vector G = Mon->GetSWF(); Lat->SetBoundaries(G);

		int n_layers_x = Lat->GetNumLayers(1)-2;
		int n_layers_y = Lat->GetNumLayers(2)-2;
		int n_layers_z = Lat->GetNumLayers(3)-2;
		int jx = (n_layers_y+2)*(n_layers_z+2);
		int jy = n_layers_z+2;
		int jz = 1;

		double fx=0,fy=0,fz=0;
		for (int x=1; x<=n_layers_x; x++)
		for (int y=1; y<=n_layers_y; y++)
		for (int z=1; z<=n_layers_z; z++) {
			if (phi[jx*x+jy*y+jz*z]>0) {
				fx += phi[jx*x+jy*y+jz*z]*(log(G[jx*(x-1)+jy*y+jz*z]) - log(G[jx*(x+1)+jy*y+jz*z]))/2;
				fy += phi[jx*x+jy*y+jz*z]*(log(G[jx*x+jy*(y-1)+jz*z]) - log(G[jx*x+jy*(y+1)+jz*z]))/2;
				fz += phi[jx*x+jy*y+jz*z]*(log(G[jx*x+jy*y+jz*(z-1)]) - log(G[jx*x+jy*y+jz*(z+1)]))/2;
			}
		}

		Out->PutReal("sys",name,"fx",fx);
		Out->PutReal("sys",name,"fy",fy);
		Out->PutReal("sys",name,"fz",fz);
	}


}

void
SF_System::GetLorealOutput(Output* Out,Lattice* Lat) const{

	if (Lat->GetNumGradients()==2) {

		int n_layers_x = Lat->GetNumLayers(1)-2;
		int n_layers_y = Lat->GetNumLayers(2)-2;
		int jx = n_layers_y+2;
		if (Lat->GetVolume()!=n_layers_x*n_layers_y) {

			double FFV=0;
			double GP=MolQ->GetGrandPotential();
			Vector w = MolQ->GetGrandPotentialProfile();
			for (int y=1; y<=n_layers_y; y++) {
				FFV +=w[jx*n_layers_x+y];
			}
			Out->PutReal("sys",name,"grand potential Exc",GP-FFV*PI*n_layers_x*n_layers_x);

			for (int i=1; i<=MolQ->GetNumMolecules(); i++) {
				SF_Molecule* Mol = MolQ->GetMolecule(i);

				Vector phi = Mol->GetPhi(total);
				double Length = Mol->GetChainLength();
				double Area = PI*n_layers_x*n_layers_x;
				double Theta = Mol->ComputeThetaExcess();
				double Phibulk = Mol->GetPhiBulk();
				FFV=0;
				for (int y=1; y<=n_layers_y; y++) {
					FFV +=phi[jx*n_layers_x+y];
				}
				FFV -=n_layers_y*Phibulk;
				Out->PutReal("mol",Mol->GetName(),"n_Exc_Exc",(Theta-FFV*Area)/Length);
			}
		}
	}

	bool surface_name_set = false;
	bool mon_set = false;
	if (!MyInput->ValueSet("sys",name,"probe_chain")) Message(fatal,"sys : probe_chain not defined");
	Text chainName = MyInput->GetText("sys",name,"probe_chain");
	if (!MolQ->MoleculeDefined(chainName)) {Message(fatal,"sys : " + name + " : probe_chain : " + chainName + " is a molecule that is not defined");}

	if (!MyInput->ValueSet("sys",name,"probe_surfactant")) Message(fatal,"sys : probe_surfactant not defined");
	Text surfactantName = MyInput->GetText("sys",name,"probe_surfactant");
	if (!MolQ->MoleculeDefined(surfactantName)) {Message(fatal,"sys : " + name + " : probe_surfactant : " + surfactantName + " is a molecule that is not defined");}

	if (MyInput->ValueSet("sys",name,"probe_surface")) surface_name_set =true;

	Text surfaceName;
	if (surface_name_set) {
		surfaceName = MyInput->GetText("sys",name,"probe_surface");
		if (!MolQ->MoleculeDefined(surfaceName)) {
			if (!SegQ->SegmentDefined(surfaceName)) {Message(fatal,"sys : " + name + " : probe_surface : " + surfaceName + " is neither a molecule nor a segment in the system");}
		} else mon_set = true;
	}


	SF_Molecule* MolPol = MolQ->GetMolecule(chainName);
	SF_Molecule* MolSurf = MolQ->GetMolecule(surfactantName);


	double ThetaPol = MolPol->ComputeThetaExcess();
	double ThetaSurf = MolSurf->ComputeThetaExcess();


	double PolLength = MolPol->GetChainLength();
	double SurfLength = MolSurf->GetChainLength();


	double ChargePol = MolPol->GetNumberOfCharges();
	double ChargeSurf = MolSurf->GetNumberOfCharges();


	if (!surface_name_set) {
		Out->PutReal("sys",name,"ChargeRatio",ThetaPol/PolLength*ChargePol/(ThetaSurf/SurfLength*ChargeSurf));

	} else {
		if (!mon_set) {
			SF_Molecule* MolSurface = MolQ->GetMolecule(surfaceName);
			double ThetaSurface = MolSurface->ComputeThetaExcess();
			double SurfaceLength = MolSurface->GetChainLength();
			double ChargeSurface = MolSurface->GetNumberOfCharges();

			Out->PutReal("sys",name,"ChargeRatio_P/S",ThetaPol/PolLength*ChargePol/(ThetaSurf/SurfLength*ChargeSurf));
			Out->PutReal("sys",name,"ChargeRatio_P/Su",ThetaPol/PolLength*ChargePol/(ThetaSurface/SurfaceLength*ChargeSurface));
			Out->PutReal("sys",name,"ChargeRatio_S/Su",ThetaSurf/SurfLength*ChargeSurf/(ThetaSurface/SurfaceLength*ChargeSurface));


		} else {

			SF_Segment* MonH = SegQ->GetSegment(surfaceName);
			double ThetaSurface = 1;
			double SurfaceLength = 1;
			double ChargeSurface = MonH->GetValence();

			Out->PutReal("sys",name,"ChargeRatio_P/S",ThetaPol/PolLength*ChargePol/(ThetaSurf/SurfLength*ChargeSurf));
			Out->PutReal("sys",name,"ChargeRatio_P/H",ThetaPol/PolLength*ChargePol/(ThetaSurface/SurfaceLength*ChargeSurface));
			Out->PutReal("sys",name,"ChargeRatio_S/H",ThetaSurf/SurfLength*ChargeSurf/(ThetaSurface/SurfaceLength*ChargeSurface));
		}

	}

	Out->PutReal("sys",name,"excess free energy",MolQ->GetExcessFreeEnergy());
	Out->PutReal("sys",name,"F_soumi",MolQ->GetSoumiFreeEnergy());

}

void
SF_System::GetSoumiOutput(Output* Out,Lattice* Lat) const{
	if (Lat->GetNumGradients()==2) {
		if (!MyInput->ValueSet("sys",name,"probe_chain")) Message(fatal,"sys : probe_chain not defined");
		Text probeName = MyInput->GetText("sys",name,"probe_chain");
		if (!MolQ->MoleculeDefined(probeName)) {Message(fatal,"sys : " + name + " : probe_chain : " + probeName + " is a molecule that is not defined");}
		SF_Molecule* Mol = MolQ->GetMolecule(probeName);
		int n_layers_x = Lat->GetNumLayers(1)-2;
		int n_layers_y = Lat->GetNumLayers(2)-2;
		int jx = n_layers_y+2;

		Vector phi = Mol->GetPhi(total);
		Lat->MultiplyWithLatticeSites(phi);
		Text dim;
		double area=PI*(n_layers_x*n_layers_x);
		double rho;
		for (int y=1; y<=n_layers_y; y++) {
			rho=0;
			for (int x=1; x<=n_layers_x; x++) {
				rho +=phi[jx*x+y];
			}
			dim = Blanks(9);
			dim.Putint(y);
			dim = Copy(dim.Strip());
			dim = Copy(dim.Frontstrip());
			Out->PutReal("mol",probeName,"rho["+ dim + "]",rho/area);
		}
		Lat->DivideByLatticeSites(phi);
	}

}

void
SF_System::GetProbeSizeOutput(Output* Out,Lattice* Lat) const{
	if (!MyInput->ValueSet("sys",name,"probe_monomer")) Message(fatal,"sys : probe_monomer not defined");
	if (!MyInput->ValueSet("sys",name,"pinned_monomer")) Message(fatal,"sys : pinned_monomer not defined");
	Text probeName = MyInput->GetText("sys",name,"probe_monomer");
	Text pinnedName = MyInput->GetText("sys",name,"pinned_monomer");
	if (!SegQ->SegmentDefined(probeName)) {Message(fatal,"sys : " + name + " : probe_monomer : " + probeName + " is a monomer that is not defined");}
	if (!SegQ->SegmentDefined(pinnedName)) {Message(fatal,"sys : " + name + " : pinned_monomer : " + pinnedName + " is a monomer that is not defined");}
	SF_Segment* MonProbe = SegQ->GetSegment(probeName);
	SF_Segment* MonPinned = SegQ->GetSegment(pinnedName);

	if (Lat->GetNumGradients()==3) {

		int n_layers_x = Lat->GetNumLayers(1)-2;
		int n_layers_y = Lat->GetNumLayers(2)-2;
		int n_layers_z = Lat->GetNumLayers(3)-2;
		int jx = (n_layers_y+2)*(n_layers_z+2);
		int jy = n_layers_z+2;
		//int jz = 1;

		Vector phi_pinned = MonPinned->GetPhi();
		double xm=0,ym=0,zm=0;
		double theta=0;
		for (int x=1; x<=n_layers_x; x++) {
			for (int y=1; y<=n_layers_y; y++) {
				for (int z=1; z<=n_layers_z; z++) {
					theta +=phi_pinned[jx*x+jy*y+z+1];
					xm +=phi_pinned[jx*x+jy*y+z+1]*x;
					ym +=phi_pinned[jx*x+jy*y+z+1]*y;
					zm +=phi_pinned[jx*x+jy*y+z+1]*(z);
				}
			}
		}
		xm/=theta; ym/=theta; zm/=theta;

		Vector phi_probe = MonProbe->GetPhi();
		double M0=0, M1x=0, M2x=0, M1y=0, M2y=0, M1z=0, M2z=0;
		for (int x=1; x<=n_layers_x; x++) {
			for (int y=1; y<=n_layers_y; y++) {
				for (int z=1; z<=n_layers_z; z++) {
					M0 +=phi_probe[jx*x+jy*y+z+1];
					M1x +=phi_probe[jx*x+jy*y+z+1]*abs(x-xm);
					M2x +=phi_probe[jx*x+jy*y+z+1]*pow(x-xm,2);
					M1y +=phi_probe[jx*x+jy*y+z+1]*abs(y-ym);
					M2y +=phi_probe[jx*x+jy*y+z+1]*pow(y-ym,2);
					M1z +=phi_probe[jx*x+jy*y+z+1]*abs(z-zm);
					M2z +=phi_probe[jx*x+jy*y+z+1]*pow(z-zm,2);
				}
			}
		}
		Out->PutReal("mon",probeName,"mx",xm);
		Out->PutReal("mon",probeName,"my",ym);
		Out->PutReal("mon",probeName,"mz",zm);
		Out->PutReal("mon",probeName,"M0",M0);
		Out->PutReal("mon",probeName,"M1x",M1x);
		Out->PutReal("mon",probeName,"M2x",M2x);
		Out->PutReal("mon",probeName,"M1y",M1y);
		Out->PutReal("mon",probeName,"M2y",M2y);
		Out->PutReal("mon",probeName,"M1z",M1z);
		Out->PutReal("mon",probeName,"M2z",M2z);
		Out->PutReal("mon",probeName,"Rgx",M1x/M0);
		Out->PutReal("mon",probeName,"Rgy",M1y/M0);
		Out->PutReal("mon",probeName,"Rgz",M1z/M0);
		Out->PutReal("mon",probeName,"Fluc_x",M2x/M0-M1x*M1x/(M0*M0));
		Out->PutReal("mon",probeName,"Fluc_y",M2y/M0-M1y*M1y/(M0*M0));
		Out->PutReal("mon",probeName,"Fluc_z",M2z/M0-M1z*M1z/(M0*M0));
	} else {
		if (Lat->GetNumGradients()==2) {
			int n_layers_x = Lat->GetNumLayers(1)-2;
			int n_layers_y = Lat->GetNumLayers(2)-2;
			int jx = n_layers_y+2;
			Vector phi_pinned = MonPinned->GetPhi();
			Lat->MultiplyWithLatticeSites(phi_pinned);
			double xm=0,ym=0;
			double theta=0;
			for (int x=1; x<=n_layers_x; x++) {
				for (int y=1; y<=n_layers_y; y++) {
					theta +=phi_pinned[jx*x+y+1];
					xm +=phi_pinned[jx*x+y+1]*x;
					ym +=phi_pinned[jx*x+y+1]*y;
				}
			}
			ym/=theta;
			Lat->DivideByLatticeSites(phi_pinned);

			Vector phi_probe = MonProbe->GetPhi();
			Lat->MultiplyWithLatticeSites(phi_probe);
			double M0=0, M1x=0, M2x=0, M1y=0, M2y=0;
			for (int x=1; x<=n_layers_x; x++) {
				for (int y=1; y<=n_layers_y; y++) {
					M0 +=phi_probe[jx*x+y+1];
					M1x +=phi_probe[jx*x+y+1]*x;
					M2x +=phi_probe[jx*x+y+1]*pow(x,2.0);
					M1y +=phi_probe[jx*x+y+1]*y;
					M2y +=phi_probe[jx*x+y+1]*pow(y,2.0);
				}
			}
			Lat->DivideByLatticeSites(phi_probe);
			Out->PutReal("mon",probeName,"mx",xm);
			Out->PutReal("mon",probeName,"my",ym);
			Out->PutReal("mon",probeName,"M0",M0);
			Out->PutReal("mon",probeName,"M1x",M1x);
			Out->PutReal("mon",probeName,"M2x",M2x);
			Out->PutReal("mon",probeName,"M1y",M1y);
			Out->PutReal("mon",probeName,"M2y",M2y);
			Out->PutReal("mon",probeName,"Rgx",M1x/M0);
			Out->PutReal("mon",probeName,"Rgy",M1y/M0);
			Out->PutReal("mon",probeName,"Fluc_x",M2x/M0-M1x*M1x/(M0*M0));
			Out->PutReal("mon",probeName,"Fluc_y",M2y/M0-M1y*M1y/(M0*M0));
		}
	}
}

void
SF_System::GetMoments3GOutput(Output* Out,Lattice* Lat) const{
	// Here I assume that in z=1 the surface is placed; trains are in layer z=2;
	if (Lat->GetNumGradients()==3) {
		if (!MyInput->ValueSet("sys",name,"probe_monomer")) Message(fatal,"sys : probe_mon not defined");
		Text probeName = MyInput->GetText("sys",name,"probe_monomer");
		if (!SegQ->SegmentDefined(probeName)) {Message(fatal,"sys : " + name + " : probe_monomer : " + probeName + " is a monomer that is not defined");}
		SF_Segment* Mon = SegQ->GetSegment(probeName);
		int n_layers_x = Lat->GetNumLayers(1)-2;
		int n_layers_y = Lat->GetNumLayers(2)-2;
		int n_layers_z = Lat->GetNumLayers(3)-2;
		int jx = (n_layers_y+2)*(n_layers_z+2);
		int jy = n_layers_z+2;
		//int jz = 1;
		double M0=0, M1=0, M2=0;
		Vector phi = Mon->GetPhi();

		for (int x=1; x<=n_layers_x; x++) {
			for (int y=1; y<=n_layers_y; y++) {
				for (int z=1; z<=n_layers_z; z++) {
					M0 +=phi[jx*x+jy*y+z];
					M1 +=phi[jx*x+jy*y+z]*(z-1.5);
					M2 +=phi[jx*x+jy*y+z]*(z-1.5)*(z-1.5);
				}
			}
		}
		Out->PutReal("mon",probeName,"M0",M0);
		Out->PutReal("mon",probeName,"M1",M1/M0);
		Out->PutReal("mon",probeName,"M2",sqrt(M2/M0));
		Out->PutReal("mon",probeName,"Fluc",M2/M0-M1*M1/M0/M0);
	}
}

void
SF_System::GetMoments2GOutput(Output* Out,Lattice* Lat) const{
	if (Lat->GetNumGradients()==2) {
		if (!MyInput->ValueSet("sys",name,"probe_monomer")) Message(fatal,"sys : probe_mon not defined");
		Text probeName = MyInput->GetText("sys",name,"probe_monomer");
		if (!SegQ->SegmentDefined(probeName)) {Message(fatal,"sys : " + name + " : probe_monomer : " + probeName + " is a monomer that is not defined");}
		SF_Segment* Mon = SegQ->GetSegment(probeName);
		int n_layers_x = Lat->GetNumLayers(1)-2;
		int n_layers_y = Lat->GetNumLayers(2)-2;
		int jx = n_layers_y+2;
		double M0=0, M1=0, M2=0;
		Vector phi = Mon->GetPhi();
		Lat->MultiplyWithLatticeSites(phi);
		for (int x=1; x<=n_layers_x; x++) {
			for (int y=1; y<=n_layers_y; y++) {
				M0 +=phi[jx*x+y];
				M1 +=phi[jx*x+y]*(y-0.5);
				M2 +=phi[jx*x+y]*(y-0.5)*(y-0.5);
			}
		}
		Out->PutReal("mon",probeName,"M0",M0);
		Out->PutReal("mon",probeName,"M1",M1/M0);
		Out->PutReal("mon",probeName,"M2",sqrt(M2/M0));
		Out->PutReal("mon",probeName,"Fluc",M2/M0-M1*M1/M0/M0);

		Lat->DivideByLatticeSites(phi);
	}
}

void
SF_System::GetEscapeOutput(Output* Out,Lattice* Lat) const{
	if (Lat->GetNumGradients()==2) {
		if (!MyInput->ValueSet("sys",name,"probe_monomer")) Message(fatal,"sys : probe_mon not defined");
		Text probeName = MyInput->GetText("sys",name,"probe_monomer");
		if (!SegQ->SegmentDefined(probeName)) {Message(fatal,"sys : " + name + " : probe_monomer : " + probeName + " is a monomer that is not defined");}

		SF_Segment* Mon = SegQ->GetSegment(probeName);
		int n_layers_x = Lat->GetNumLayers(1)-2;
		int n_layers_y = Lat->GetNumLayers(2)-2;
		int jx = n_layers_y+2;


		double rmin1=0,rmin2=0,rmax1=0;
		double phimax1=0,phimax2=0,phimin1=0;
		bool min1=false,max1=false,min2=false;
		double F_1=0,F0=0,F1=0,Fmin1=0,Fmin2=0,Fmax1=0;
		double a,b,c,xx;


		Vector phi = Mon->GetPhi();
		Text dim;
		double theta;
		//double FM;
		for (int x=1; x<=n_layers_x; x++) {
			theta = 0; //FM=0;
			for (int y=1; y<=n_layers_y; y++) {theta += phi[jx*x+y];}
				F_1=F0; F0=F1; F1=-log(theta);
			    if ((theta>phimax1) && !min1 && x>3) {phimax1=theta; rmin1=x;}
			    else
			    {
					if (rmin1>0 && !min1) {
						min1=true;
						c=F0;
						a=0.5*(F_1+F1-2*F0);
						b=-0.5*(F_1-F1);
						xx=-b/(2*a);
						rmin1=rmin1+xx;
						Fmin1=a*xx*xx + b*xx + c;
						phimin1=theta;
					}
			    }
			if (min1) {
				if ((theta<phimin1) && !max1) {phimin1=theta; rmax1=x;}
				else {
					if (rmax1>0 && !max1) {
						max1=true;
						c=F0;
						a=0.5*(F_1+F1-2*F0);
						b=-0.5*(F_1-F1);
						xx=-b/(2*a);
						rmax1=rmax1+xx;
						Fmax1=a*xx*xx + b*xx + c;
						phimax2=theta;
					}
				}
			}
			if (max1 && min1) {
				if ((theta>phimax2) && !min2) {phimax2=theta; rmin2=x;}
				else {
					if (rmin2>0 && !min2) {
						min2=true;
						c=F0;
						a=0.5*(F_1+F1-2*F0);
						b=-0.5*(F_1-F1);
						xx=-b/(2*a);
						rmin2=rmin2+xx;
						Fmin2=a*xx*xx + b*xx + c;
					}
				}
			}

			dim = Blanks(9);
			dim.Putint(x);
			dim = Copy(dim.Strip());
			dim = Copy(dim.Frontstrip());

			Out ->PutReal("mon", probeName,"n[" + dim + "]" ,theta);
		}
		if (min1) {
			Out->PutReal("mon",probeName,"phimax1",phimax1);
			Out->PutReal("mon",probeName,"Fmin1",Fmin1);
			Out->PutReal("mon",probeName,"rmin1",rmin1);
		} else {
			Out->PutReal("mon",probeName,"phimax1",0);
			Out->PutReal("mon",probeName,"Fmin1",0);
			Out->PutReal("mon",probeName,"rmin1",0);
		}

		if (max1) {
			Out->PutReal("mon",probeName,"phimin1",phimin1);
			Out->PutReal("mon",probeName,"Fmax1",Fmax1);
			Out->PutReal("mon",probeName,"rmax1",rmax1);
		} else {
			Out->PutReal("mon",probeName,"phimin1",0);
			Out->PutReal("mon",probeName,"Fmax1",0);
			Out->PutReal("mon",probeName,"rmax1",0);
		}

		if (min2) {
			Out->PutReal("mon",probeName,"phimax2",phimax2);
			Out->PutReal("mon",probeName,"Fmin2",Fmin2);
			Out->PutReal("mon",probeName,"rmin2",rmin2);
		} else {
			Out->PutReal("mon",probeName,"phimax2",0);
			Out->PutReal("mon",probeName,"Fmin2",0);
			Out->PutReal("mon",probeName,"rmin2",0);
		}

	}
}

void
SF_System::GetBrush3GOutput(Output* Out,Lattice* Lat) const{
	if (Lat->GetNumGradients()==3) {
		if (!MyInput->ValueSet("sys",name,"probe_chain")) Message(fatal,"sys : probe_chain not defined");
		Text probeName = MyInput->GetText("sys",name,"probe_chain");
		if (!MolQ->MoleculeDefined(probeName)) {Message(fatal,"sys : " + name + " : probe_chain : " + probeName + " is a molecule that is not defined");}

		SF_Molecule* Mol = MolQ->GetMolecule(probeName);
		int n_layers_x = Lat->GetNumLayers(1)-2;
		int n_layers_y = Lat->GetNumLayers(2)-2;
		int n_layers_z = Lat->GetNumLayers(3)-2;
		int jx = (n_layers_y+2)*(n_layers_z+2);
		int jy = n_layers_z+2;
		int jz = 1;

		Vector phi = Mol->GetPhi(total);
		// Here  I assume that there is a surface at z = 1 and that segment in layer 2 is 0.5 distance from the surface.

		double H=0;
		double theta=0;

		for (int x=1; x<=n_layers_x; x++) {
			for (int y=1; y<=n_layers_y; y++) {
				for (int z=2; z<=n_layers_z; z++) {
					theta += phi[jx*x+jy*y+jz*z];
					H += (z-1.5)*phi[jx*x+jy*y+jz*z];
				}
			}
		}

		Out ->PutReal("mol", Mol->GetName(),"<H>",H/theta);
	}
}


void
SF_System::GetEgorovOutput(Output* Out,Lattice* Lat) const{
	if (Lat->GetNumGradients()==2) {
		if (!MyInput->ValueSet("sys",name,"probe_chain")) Message(fatal,"sys : probe_chain not defined");
		Text probeName = MyInput->GetText("sys",name,"probe_chain");
		if (!MolQ->MoleculeDefined(probeName)) {Message(fatal,"sys : " + name + " : probe_chain : " + probeName + " is a molecule that is not defined");}

		SF_Molecule* Mol = MolQ->GetMolecule(probeName);
		int n_layers_x = Lat->GetNumLayers(1)-2;
		int n_layers_y = Lat->GetNumLayers(2)-2;
		int jx = n_layers_y+2;

		Vector phi = Mol->GetPhi(total);
		Lat->MultiplyWithLatticeSites(phi);
		Text dim;
		double theta, FM;
		for (int x=1; x<=n_layers_x; x++) {
			theta = 0; FM=0;
			for (int y=1; y<=n_layers_y; y++) {theta += phi[jx*x+y]; FM += y*phi[jx*x+y];}
			dim = Blanks(9);
			dim.Putint(x);
			dim = Copy(dim.Strip());
			dim = Copy(dim.Frontstrip());

			Out ->PutReal("mol", Mol->GetName(),"H[" + dim + "]" ,FM/theta);
		}

		Lat->DivideByLatticeSites(phi);
	}
}


void
SF_System::GetPoreOutput(Output* Out,Lattice* Lat) const{
	int M = Lat->GetTotalNumLayers();
	SF_Molecule* Mol = MolQ->GetMolecule("lipid");
	if (Mol == NULL) {
		Message(literal, "Don't know how to generate vesicle output, "
			"please define a molecule named 'lipid'");
		return;
	}
	Text latName;
	MyInput->GetNumNames("lat",1,1);
	Array<Text> latNames = MyInput->GetNames("lat");
	latName = latNames[1];
	Array<Text> geometries(1,3);
	geometries[1] = "flat";
	geometries[2] = "cylindrical";
	geometries[3] = "spherical";
	int geom = MyInput->GetChoice("lat",latName,"geometry",geometries,1);
	if (geom !=2 ) {
		Message(literal, "Pores are studied in cylindrical geometry");
		return;
	}
	int dim = Lat->GetNumGradients();
	if (dim!=2) {
		Message(literal, "Pores are studied in two gradient system");
		return;
	}
	int n_layers_x = Lat->GetNumLayers(1)-2;
	int n_layers_y = Lat->GetNumLayers(2)-2;
	int jx = n_layers_y+2;

	Vector phi = Mol->GetPhi(total);
	Lat->SubtractBoundaries(phi);
	double phibulk = Mol->GetPhiBulk();
	int z;
	for (z=1; z<=M; z++) {
		phi[z] -= phibulk;
	}
	double theta_m=0;
	for (int y=1; y<=n_layers_y; y++) {theta_m += phi[jx*n_layers_x+y]; }
	Lat->MultiplyWithLatticeSites(phi);
	double theta=0;
	for (z=1; z<=M; z++) {
		theta += phi[z];
	}
	Lat->DivideByLatticeSites(phi);
	double R = sqrt(n_layers_x*n_layers_x-theta/theta_m/PI);
	Out ->PutReal("sys",name,"Pore Radius",R);
	Vector w = MolQ->GetGrandPotentialProfile();
	Lat->SubtractBoundaries(w);
	double tension=0;
	for (int y=1; y<=n_layers_y; y++) {tension += w[jx*n_layers_x+y]; }
	Lat->MultiplyWithLatticeSites(w);
	double GP=0;
	for (z=1; z<=M; z++) {
		GP+=w[z];
	}
	Lat->DivideByLatticeSites(w);
	double PoreG=GP-(n_layers_x*n_layers_x-R*R)*PI*tension;
	Out ->PutReal("sys",name,"Pore GrandPotential",2*PoreG);
	Out ->PutReal("sys",name,"Pore LineTension",2*PoreG/(2*PI*R));
	Out ->PutReal("sys",name,"Pore MembraneTension",tension*2);
}

void
SF_System::GetModuliOutput(Output* Out,Lattice* Lat) const{
	double kc=0,kappa=0,kmean=0,J=0;
	double large=10000000.0;
	int M = Lat->GetTotalNumLayers();
	SF_Molecule* Mol = MolQ->GetMolecule("lipid");
	if (Mol == NULL) {
		Message(literal, "Don't know how to generate vesicle output, "
			"please define a molecule named 'lipid'");
		return;
	}
	Vector phi = Mol->GetPhi(total);
	Vector w = MolQ->GetGrandPotentialProfile();
	Lat->SubtractBoundaries(w);
	Lat->MultiplyWithLatticeSites(w);
	double GP=0;
	Lat->SubtractBoundaries(phi);
	double phibulk = Mol->GetPhiBulk();
	int z;
	for (z=1; z<=M; z++) {
		phi[z] -= phibulk;
		GP+=w[z];
	}
	Lat->DivideByLatticeSites(w);
	LatticeRange* LayerZero = new LatticeRange1D(1,1,M);
	double R = Lat->MomentUnweighted(phi,1,LayerZero,-0.5)/Lat->MomentUnweighted(phi,0,LayerZero,-0.5);
	delete LayerZero;
	Text latName;
	MyInput->GetNumNames("lat",1,1);
	Array<Text> latNames = MyInput->GetNames("lat");
	latName = latNames[1];
	Array<Text> geometries(1,3);
	geometries[1] = "flat";
	geometries[2] = "cylindrical";
	geometries[3] = "spherical";

	int geom = MyInput->GetChoice("lat",latName,"geometry",geometries,1);
	switch (geom) {
		case 1:
			kc=large;
			kappa=large;
			kmean=large;
			J=0;
			break;
		case 2:
			kc=GP*R/PI;
			kappa=kc/2;
			kmean=large;
			J=1/R;
			break;
		case 3:
			kc=large;
			kappa=large;
			kmean=GP/(4.0*PI);
			J=2/R;
			break;
	}

	Out ->PutReal("sys",name,"k_c",kc);
	Out ->PutReal("sys",name,"kappa",kappa);
	Out ->PutReal("sys",name,"kmean",kmean);
	Out ->PutReal("sys",name,"J",J);
}


void
SF_System::GetDSMOutput(Output* Out,Lattice* Lat) const{

	int n_layers = Lat->GetNumLayers(1)-2;
	if (MyInput->ValueSet("sys",name,"core_monomer")) {
		Text CoreMonName = MyInput->GetText("sys",name,"core_monomer");
		if (!SegQ->SegmentDefined(CoreMonName)) {Message(fatal,"sys : " + name + " : core_monomer : " + CoreMonName + " is a monomer that is not defined");}
		SF_Segment* CoreMon;
		CoreMon = SegQ->GetSegment(CoreMonName);
		Vector phi = CoreMon->GetPhi();
		double z_value=0;
		bool target_found=false;
		int z=n_layers-1;
		while (!target_found && z > 0) {
			if (phi[z] > phi[1]/2) {
				target_found = true;
				z_value = z-1+ (phi[z]-phi[1]/2)/(phi[z]-phi[z+1]);
			}
			z=z-1;
		}
		Out->PutReal("sys",name,"core-size",z_value);
	}
	if (MyInput->ValueSet("sys",name,"probe_monomer")) {
		Text ProbeMonName = MyInput->GetText("sys",name,"probe_monomer");
		if (!SegQ->SegmentDefined(ProbeMonName)) {Message(fatal,"sys : " + name + " : probe_monomer : " + ProbeMonName + " is a monomer that is not defined");}
		SF_Segment* ProbeMon;
		ProbeMon = SegQ->GetSegment(ProbeMonName);
		Vector phi = ProbeMon->GetPhi();
		int zmax1=0,zmax2=0,zmax3=0,zmin1=0,zmin2=0;
		double phimax1=0,phimax2=0,phimax3=0,phimin1=0,phimin2=0;
		bool max1=false,min1=false,max2=false,max3=false,min2=false;

		for (int z=1; z<=n_layers; z++) {
			if ((phi[z+1]>phimax1) && !max1) {phimax1=phi[z+1]; zmax1=z;} else {if (zmax1>0) max1=true;}
			if ((phi[n_layers-z+2] > phimax3) && !max3) {phimax3 = phi[n_layers-z+2]; zmax3=n_layers-z+1;} else {if (zmax3>0) max3=true;}
		}
		if (zmax1<zmax3){
			phimin1=phimax1;
			for(int z=zmax1+1; z<zmax3; z++) {
				if (phi[z+1]<phimin1 && !min1) {phimin1=phi[z+1]; zmin1=z;} else {if (zmin1>0) min1=true;}
			}
			phimin2=phimax3;
			for(int z=zmax3-1; z> zmax1; z--) {
				if (phi[z+1]<phimin2 && !min2) {phimin2=phi[z+1]; zmin2=z;} else {if (zmin2>0) min2=true;}
			}
		}
		if (zmin1<zmin2){
			phimax2=phimin1;
			for(int z=zmin1+1; z< zmin2; z++) {
				if (phi[z+1]>phimax2 ) {phimax2=phi[z+1]; zmax2=z; max2=true;};
			}
		}


		if (max1) Out->PutReal("mon",ProbeMonName,"max1",phimax1); else Out->PutReal("mon",ProbeMonName,"max1",0);
		if (min1) Out->PutReal("mon",ProbeMonName,"min1",phimin1); else Out->PutReal("mon",ProbeMonName,"min1",0);
		if (max2) Out->PutReal("mon",ProbeMonName,"max2",phimax2); else Out->PutReal("mon",ProbeMonName,"max2",0);
		if (min2) Out->PutReal("mon",ProbeMonName,"min2",phimin2); else Out->PutReal("mon",ProbeMonName,"min2",0);
		if (max3) Out->PutReal("mon",ProbeMonName,"max3",phimax3); else Out->PutReal("mon",ProbeMonName,"max3",0);


		if (max1) Out->PutInt("mon",ProbeMonName,"z_max1",zmax1); else Out->PutInt("mon",ProbeMonName,"z_max1",0);
		if (min1) Out->PutInt("mon",ProbeMonName,"z_min1",zmin1); else Out->PutInt("mon",ProbeMonName,"z_min1",0);
		if (max2) Out->PutInt("mon",ProbeMonName,"z_max2",zmax2); else Out->PutInt("mon",ProbeMonName,"z_max2",0);
		if (min2) Out->PutInt("mon",ProbeMonName,"z_min2",zmin2); else Out->PutInt("mon",ProbeMonName,"z_min2",0);
		if (max3) Out->PutInt("mon",ProbeMonName,"z_max3",zmax3); else Out->PutInt("mon",ProbeMonName,"z_max3",0);

		if (max1) Out->PutReal("mon",ProbeMonName,"Fmin1",-log(phimax1)); else Out->PutReal("mon",ProbeMonName,"Fmin1",0);
		if (max2) Out->PutReal("mon",ProbeMonName,"Fmin2",-log(phimax2)); else Out->PutReal("mon",ProbeMonName,"Fmin2",0);
		if (max3) Out->PutReal("mon",ProbeMonName,"Fmin3",-log(phimax3)); else Out->PutReal("mon",ProbeMonName,"Fmin3",0);
		if (min1) Out->PutReal("mon",ProbeMonName,"Fmax1",-log(phimin1)); else Out->PutReal("mon",ProbeMonName,"Fmax1",0);
		if (min2) Out->PutReal("mon",ProbeMonName,"Fmax2",-log(phimin2)); else Out->PutReal("mon",ProbeMonName,"Fmax2",0);

		if (max1 && min1) Out->PutReal("mon",ProbeMonName,"DeltaF11",log(phimax1/phimin1)); else Out->PutReal("mon",ProbeMonName,"DeltaF11",0);
		if (max1 && min2) Out->PutReal("mon",ProbeMonName,"DeltaF12",log(phimax1/phimin2)); else Out->PutReal("mon",ProbeMonName,"DeltaF12",0);
		if (max2 && min1) Out->PutReal("mon",ProbeMonName,"DeltaF21",log(phimax2/phimin1)); else Out->PutReal("mon",ProbeMonName,"DeltaF21",0);
		if (max2 && min2) Out->PutReal("mon",ProbeMonName,"DeltaF22",log(phimax2/phimin2)); else Out->PutReal("mon",ProbeMonName,"DeltaF22",0);
		if (max3 && min1) Out->PutReal("mon",ProbeMonName,"DeltaF31",log(phimax3/phimin1)); else Out->PutReal("mon",ProbeMonName,"DeltaF31",0);
		if (max3 && min2) Out->PutReal("mon",ProbeMonName,"DeltaF32",log(phimax3/phimin2)); else Out->PutReal("mon",ProbeMonName,"DeltaF32",0);



		if (zmin1>0) {
			double phi1,phi2,phi3;
			Lat->MultiplyWithLatticeSites(phi);
			phi1=0; phi2=0; phi3=0;
			for (int z=1; z<zmin1; z++) {
				phi1=phi1+phi[z];
			}

			for (int z=zmin1; z<=n_layers; z++) {
				phi2=phi2+phi[z];
			}
			phi3=phi1+phi2;
			Out->PutReal("sys",name,"f1",phi1/phi3);
			Out->PutReal("sys",name,"f2",phi2/phi3);

			Lat->DivideByLatticeSites(phi);
		} else {
			Out->PutReal("sys",name,"f1",0);
			Out->PutReal("sys",name,"f2",0);
		}
	} else {
		Message(literal,"sys : " + name + " : probe_monomer : mon-name, is expected but not found. Most of DSM output skipped. ");
	}
}

void
SF_System::GetSergioOutput(Output* Out,Lattice* Lat) const{
	int numSeg = SegQ->GetNumSegments();
	int M=Lat->GetNumLayers(1);
	double phiBulk;
	for (int i=1; i<=numSeg; i++) {
		SF_Segment* Seg = SegQ->GetSegment(i);
		Vector phiExc = Seg->GetPhi();
		phiBulk = Seg->GetPhiBulk();
		int z;
		for (z=1; z<=M; z++) phiExc[z] -= phiBulk;
		LatticeRange* LayerZero = new LatticeRange1D(1,1,M);
		double R1 = Lat->MomentUnweighted(phiExc,1,LayerZero,-0.5)/Lat->MomentUnweighted(phiExc,0,LayerZero,-0.5);
		double R2 = Lat->MomentUnweighted(phiExc,2,LayerZero,-0.5)/Lat->MomentUnweighted(phiExc,0,LayerZero,-0.5);
		double R3 = Lat->MomentUnweighted(phiExc,3,LayerZero,-0.5)/Lat->MomentUnweighted(phiExc,0,LayerZero,-0.5);
		double R4 = Lat->MomentUnweighted(phiExc,4,LayerZero,-0.5)/Lat->MomentUnweighted(phiExc,0,LayerZero,-0.5);
		double R5 = Lat->MomentUnweighted(phiExc,5,LayerZero,-0.5)/Lat->MomentUnweighted(phiExc,0,LayerZero,-0.5);
		double R6 = Lat->MomentUnweighted(phiExc,6,LayerZero,-0.5)/Lat->MomentUnweighted(phiExc,0,LayerZero,-0.5);
		double R7 = Lat->MomentUnweighted(phiExc,7,LayerZero,-0.5)/Lat->MomentUnweighted(phiExc,0,LayerZero,-0.5);
		double R8 = Lat->MomentUnweighted(phiExc,8,LayerZero,-0.5)/Lat->MomentUnweighted(phiExc,0,LayerZero,-0.5);
		double R9 = Lat->MomentUnweighted(phiExc,9,LayerZero,-0.5)/Lat->MomentUnweighted(phiExc,0,LayerZero,-0.5);
		double R10 = Lat->MomentUnweighted(phiExc,10,LayerZero,-0.5)/Lat->MomentUnweighted(phiExc,0,LayerZero,-0.5);

		delete LayerZero;
		for (z=1; z<=M; z++) {phiExc[z] += phiBulk;	}
		Out->PutReal("mon",Seg->GetName(),"1st moment exc.",R1);
		Out->PutReal("mon",Seg->GetName(),"2nd moment exc.",pow(R2,0.5));
		double norm=1.0/3.0;
		Out->PutReal("mon",Seg->GetName(),"3nd moment exc.",pow(R3,norm));norm=0.25;
		Out->PutReal("mon",Seg->GetName(),"4nd moment exc.",pow(R4,norm));norm=0.2;
		Out->PutReal("mon",Seg->GetName(),"5nd moment exc.",pow(R5,norm));norm=1.0/6.0;
		Out->PutReal("mon",Seg->GetName(),"6nd moment exc.",pow(R6,norm));norm=1.0/7.0;
		Out->PutReal("mon",Seg->GetName(),"7nd moment exc.",pow(R7,norm));norm=0.125;
		Out->PutReal("mon",Seg->GetName(),"8nd moment exc.",pow(R8,norm));norm=1.0/9.0;
		Out->PutReal("mon",Seg->GetName(),"9nd moment exc.",pow(R9,norm));norm=0.1;
		Out->PutReal("mon",Seg->GetName(),"10nd moment exc.",pow(R10,norm));

	}

	int n_layers = Lat->GetNumLayers(1)-2;
	for (int i=1; i<=MolQ->GetNumMolecules(); i++) {
		SF_Molecule* Mol = MolQ->GetMolecule(i);
		if ( MyInput->ValueSet("mol",Mol->GetName(),"find_phi_value")) {
			double phi_value=MyInput->GetReal("mol",Mol->GetName(),"find_phi_value",0,1,0.01);
			Vector phi = Mol->GetPhi(total);
			double z_value=0;
			bool target_found=false;
			int z=n_layers-1;
			while (!target_found && z > 0) {
				if (phi[z] > phi_value) {
					target_found = true;
					z_value = z-1+ (phi[z]-phi_value)/(phi[z]-phi[z+1]);
				}
				z=z-1;
			}
			Out->PutReal("mol",Mol->GetName(),"z-value",z_value);
		}
	}
}


void
SF_System::GetEmiliaOutput(Output* Out,Lattice* Lat) const{

	int n_layers = Lat->GetNumLayers(1)-2;

       if (MyInput->ValueSet("sys",name,"Gibbs_molecule")) {
		double GibbsPlane;
		double T1,TM,Tt,TG;
		MyInput->AlwaysCombineParam("sys",name,"Gibbs_molecule","extra_output");
		Text GibbsMolName = MyInput->GetText("sys",name,"Gibbs_molecule");
		if (!MolQ->MoleculeDefined(GibbsMolName)) {Message(fatal,"'sys : " + name + " : Gibbs_molecule : " + GibbsMolName + " is a molecule that is not defined");}
		SF_Molecule* GibbsMol;
		GibbsMol = MolQ->GetMolecule(GibbsMolName);
		Vector phi = GibbsMol->GetPhi(total);
		Tt=0;
		for (int z=1; z<=n_layers; z++) Tt +=phi[z];
		T1 = phi[1];
		TM = phi[n_layers];
		if (T1==TM) {GibbsPlane = 0;} else {GibbsPlane = (Tt-TM*n_layers)/(T1-TM);}
		Out->PutReal("sys",name,"GibbsPlane",GibbsPlane);

		for (int i=1; i<=MolQ->GetNumMolecules(); i++) {

			SF_Molecule* Mol = MolQ->GetMolecule(i);
			Text molName = Mol->GetName();
			Vector phi = Mol->GetPhi(total);
			Tt=0;
			for (int z=1; z<=n_layers; z++) {Tt +=phi[z];}

			TG= Tt-GibbsPlane*phi[1]-(n_layers - GibbsPlane)*phi[n_layers];
			Out ->PutReal("mol",molName,"Theta_Gibbs",TG);
		}
	}
}


void
SF_System::GetCeciliaOutput(Output* Out,Lattice* Lat) const{

    	Vector GPProfile = MolQ->GetGrandPotentialProfile();
	int n_layers_x = Lat->GetNumLayers(1)-2;
	int n_layers_y = Lat->GetNumLayers(2)-2;
	int jx = n_layers_y+2;
	//int jy = 1;
       double tension=0;
	for (int y=1; y<=n_layers_y; y++) {tension +=GPProfile[jx+y];}
	Out->PutReal("sys",name,"tension(1)",tension);


       if (MyInput->ValueSet("sys",name,"Gibbs_molecule")) {
		double GibbsPlane;
		double T1,TM,Tt,TG;
		MyInput->AlwaysCombineParam("sys",name,"Gibbs_molecule"," ");
		Text GibbsMolName = MyInput->GetText("sys",name,"Gibbs_molecule");
		if (!MolQ->MoleculeDefined(GibbsMolName)) {Message(fatal,"'sys : " + name + " : Gibbs_molecule : " + GibbsMolName + " is a molecule that is not defined");}
		SF_Molecule* GibbsMol;
		GibbsMol = MolQ->GetMolecule(GibbsMolName);
		Vector phi = GibbsMol->GetPhi(total);
		Array<double> theta_x(n_layers_x+1);
		Tt=0;
		for (int x=1; x<=n_layers_x; x++) {
			theta_x[x] = 0;
			for (int y=1; y<=n_layers_y; y++) {theta_x[x] += phi[jx*x+y]; }
			Tt += theta_x[x];
		}
		T1=theta_x[1]; TM = theta_x[n_layers_x];

		if (T1==TM) {GibbsPlane = 0;} else {GibbsPlane = (Tt-TM*n_layers_x)/(T1-TM);}
		Out->PutReal("sys",name,"GibbsPlane",GibbsPlane);

		for (int i=1; i<=MolQ->GetNumMolecules(); i++) {

			SF_Molecule* Mol = MolQ->GetMolecule(i);
			Text molName = Mol->GetName();
			Vector phi = Mol->GetPhi(total);
			Tt=0;
			for (int x=1; x<=n_layers_x; x++) {
				theta_x[x] = 0;
				for (int y=1; y<=n_layers_y; y++) {theta_x[x] += phi[jx*x+y]; 	}
				Tt += theta_x[x];
			}

			TG= Tt-GibbsPlane*theta_x[1]-(n_layers_x - GibbsPlane)*theta_x[n_layers_x];
			Out ->PutReal("sys",name,"Theta_Gibbs" + molName,TG);
		}
	}
}

void
SF_System::GetOutputHairy(Output* Out, Lattice* Lat) const {
	LatticeRange* PinnedRange=NULL;
	for (int i=1; i<=MolQ->GetNumMolecules(); i++) {
		SF_Molecule* Mol = MolQ->GetMolecule(i);
		if (Mol->GetChainLength() > 1) {
			SF_MolStructure* Chain = Mol->GetMolStructure();
			int j;
			for (j=1; j<=Chain->GetNumDiffSegments(); j++) {
				SF_MolSegment* MolSeg = Chain->GetDiffSegment(j);
				if (MolSeg->GetFreedom() == pinned
					|| MolSeg->GetFreedom() == grafted) {
					PinnedRange = MolSeg->GetLatRange();
				}
			}
			for (j=1; j<=Chain->GetNumDiffSegments(); j++) {
				SF_MolSegment* MolSeg = Chain->GetDiffSegment(j);
				if (MolSeg->GetFreedom() == loose) {
					Vector phi = MolSeg->GetPhi(total);
					double segTheta = Lat->Moment(phi,0,PinnedRange,0);
					double distance = Lat->Moment(phi,1,PinnedRange,0)*Lat->GetSiteDistance()/segTheta;
					Out->PutReal("mol",Mol->GetName(),"av distance " + MolSeg->GetName(),distance);
					distance = Lat->Moment(phi,2,PinnedRange,0)*Lat->GetSiteDistance()/segTheta;
					Out->PutReal("mol",Mol->GetName(),"av dist^2 " + MolSeg->GetName(),distance);
				}
			}
		}
	}
}
void
SF_System::GetOutputPressure(Output* Out, Lattice* Lat) const {
	Vector pressureProfile = MolQ->GetGrandPotentialProfile();
	LatticeRange* LatRange = Lat->NewLatticeRange("lowerbound");
	Out->PutReal("sys",name,"P0",
		Lat->MomentUnweighted(pressureProfile,0,LatRange,-0.5));
	Out->PutReal("sys",name,"P1",
		Lat->MomentUnweighted(pressureProfile,1,LatRange,-0.5));
	Out->PutReal("sys",name,"P2",
		Lat->MomentUnweighted(pressureProfile,2,LatRange,-0.5));
	delete LatRange;
	LatRange = Lat->NewLatticeRange("upperbound");
	Out->PutReal("sys",name,"P1M",
		Lat->MomentUnweighted(pressureProfile,1,LatRange,+0.5));
	Out->PutReal("sys",name,"P2M",
		Lat->MomentUnweighted(pressureProfile,2,LatRange,+0.5));
	delete LatRange;
	AddChemPotGraft(pressureProfile,Lat);
	LatRange = Lat->NewLatticeRange("lowerbound");
	Out->PutReal("sys",name,"F0",
		Lat->MomentUnweighted(pressureProfile,0,LatRange,-0.5));
	Out->PutReal("sys",name,"F1",
		Lat->MomentUnweighted(pressureProfile,1,LatRange,-0.5));
	Out->PutReal("sys",name,"F2",
		Lat->MomentUnweighted(pressureProfile,2,LatRange,-0.5));
	delete LatRange;
	LatRange = Lat->NewLatticeRange("upperbound");
	Out->PutReal("sys",name,"F1M",
		Lat->MomentUnweighted(pressureProfile,1,LatRange,+0.5));
	Out->PutReal("sys",name,"F2M",
		Lat->MomentUnweighted(pressureProfile,2,LatRange,+0.5));
	delete LatRange;
}
void
SF_System::GetOutputJoanne(Output* Out, Lattice* Lat) const {
	int i,j,k,z,M;
	M = Lat->GetTotalNumLayers();
	for (i=1; i<=MolQ->GetNumMolecules(); i++) {
		SF_Molecule* Mol = MolQ->GetMolecule(i);
		if (Mol->GetChainLength() > 1) {
			Vector phi = Mol->GetPhi(total);
			Vector lnPhi(1,M);
			for (z=1; z<=M; z++) {
				lnPhi[z] = log(phi[z]);
			}
			Out->PutProfile("mol",Mol->GetName(),"ln(phi)",lnPhi);
			SF_MolStructure* Chain = Mol->GetMolStructure();
			LatticeRange* PinnedRange = NULL;
			for (j=1; j<=Chain->GetNumDiffSegments(); j++) {
				SF_MolSegment* MolSeg = Chain->GetDiffSegment(j);
				if (MolSeg->GetFreedom() == pinned
					|| MolSeg->GetFreedom() == grafted) {
					PinnedRange = MolSeg->GetLatRange();
				}
			}
			if (PinnedRange != NULL) {
				for (j=1; j<=Chain->GetNumDiffSegments(); j++) {
					SF_MolSegment* MolSeg = Chain->GetDiffSegment(j);
					if (MolSeg->GetFreedom() == loose) {
						phi = MolSeg->GetPhi(total);
						Vector lnPhi(1,M);
						for (z=1; z<=M; z++) {
							lnPhi[z] = log(phi[z]);
						}
						Out->PutProfile("mol",Mol->GetName(),"ln(phi-" + MolSeg->GetName() + ")",lnPhi);
						double Re = Lat->Moment(phi,1,PinnedRange,0);
						Re /= Lat->Moment(phi,0,PinnedRange,0);
						Out->PutReal("mol",Mol->GetName(),"Re-" + MolSeg->GetName(),Re);
						Out->PutReal("mol",Mol->GetName(),"ln(Re-" + MolSeg->GetName() + ")",log(Re));
					}
					double alphaAv = 0;
					for (k=1; k<=MolSeg->GetNumStates(); k++) {
						SF_MolState* State = MolSeg->GetState(k);
						phi = State->GetPhi(total);
						alphaAv += Lat->Moment(phi,0,PinnedRange,0)*State->GetValence();
					}
					phi = MolSeg->GetPhi(total);
					alphaAv /= Lat->Moment(phi,0,PinnedRange,0);
					Out->PutReal("mol",Mol->GetName(),"alpha_average-" + MolSeg->GetName(),alphaAv);
				}
			}
			Out->PutReal("mol",Mol->GetName(),"theta/N",Mol->GetTheta()/Mol->GetChainLength());
			Out->PutReal("mol",Mol->GetName(),"ln(theta/N)",log(Mol->GetTheta()/Mol->GetChainLength()));
		} else {
			Out->PutReal("mol",Mol->GetName(),"ln(phibulk)",log(Mol->GetPhiBulk()));
		}
	}
	GetOutputMultiState(Out, Lat);
}
void
SF_System::GetOutputSecondGen(Output* Out, Lattice* Lat) const {
	int i,j;
	double Rg,fluct;
	int M = Lat->GetTotalNumLayers();
	Vector phi;
	LatticeRange* PinnedRange = new LatticeRange1D(1,1,M);
	for (i=1; i<=MolQ->GetNumMolecules(); i++) {
		SF_Molecule* Mol = MolQ->GetMolecule(i);
		if (Mol->GetChainLength() > 1) {
			SF_MolStructure* Chain = Mol->GetMolStructure();
			for (j=1; j<=Chain->GetNumDiffSegments(); j++) {
				SF_MolSegment* MolSeg = Chain->GetDiffSegment(j);
				if (MolSeg->GetFreedom() == loose) {
					if (Mol->GetFreedom() == secondGeneration) {
						phi = MolSeg->GetPhi(constrained);
					} else {
						phi = MolSeg->GetPhi(total);
					}
					Rg = Lat->Moment(phi,2,PinnedRange,-0.5);
					double N = Lat->Moment(phi,0,PinnedRange,-0.5);
					Rg /= N;
					Out->PutReal("mol",Mol->GetName(),"Rg-" + MolSeg->GetName(),sqrt(Rg));
					fluct = Lat->Moment(phi,4,PinnedRange,-0.5);
					fluct /= N;
					fluct -= Rg*Rg;
					fluct /= Rg*Rg;
					Out->PutReal("mol",Mol->GetName(),"fluct-" + MolSeg->GetName(),fluct);
					if (Lat->GetNumGradients()==2) {
						Rg = Lat->Moment(phi,2,NULL,1);
						Rg /= N;
						Out->PutReal("mol", Mol->GetName(),"Rgx-"
							 + MolSeg->GetName(),sqrt(Rg));
						fluct = Lat->Moment(phi,4,NULL,1);
						fluct /= N;
						fluct -= Rg*Rg;
						fluct /= Rg*Rg;
						Out->PutReal("mol",Mol->GetName(),"fluctx-"
							 + MolSeg->GetName(), fluct);
						Rg = Lat->Moment(phi,2,NULL,2);
						Rg /= N;
						Out->PutReal("mol",Mol->GetName(),"Rgy-"
							+ MolSeg->GetName(), sqrt(Rg));
						fluct = Lat->Moment(phi,4,NULL,2);
						fluct /= N;
						fluct -= Rg*Rg;
						fluct /= Rg*Rg;
						Out->PutReal("mol",Mol->GetName(),"flucty-"
							+ MolSeg->GetName(), fluct);
						Rg = Lat->Moment(phi,1,NULL,2); //y coordinate of CoM
						Rg /= N;
						Out->PutReal("mol",Mol->GetName(),"CoM-"
							+ MolSeg->GetName(), Rg);
					}
				}
			}
			if (Mol->GetFreedom() == secondGeneration) {
				phi = Mol->GetPhi(constrained);
			} else {
				phi = Mol->GetPhi(total);
			}
			Rg = Lat->Moment(phi,2,PinnedRange,-0.5);
			double N = Lat->Moment(phi,0,PinnedRange,-0.5);
			Rg /= N;
			Out->PutReal("mol",Mol->GetName(),"Rg",sqrt(Rg));
			fluct = Lat->Moment(phi,4,PinnedRange,-0.5);
			fluct /= N;
			fluct -= Rg*Rg;
			fluct /= Rg*Rg;
			Out->PutReal("mol",Mol->GetName(),"fluct",fluct);
			if (Lat->GetNumGradients()==2) {
				Rg = Lat->Moment(phi,2,NULL,1);
				Rg /= N;
				Out->PutReal("mol", Mol->GetName(),"Rgx",
					sqrt(Rg));
				fluct = Lat->Moment(phi,4,NULL,1);
				fluct /= N;
				fluct -= Rg*Rg;
				fluct /= Rg*Rg;
				Out->PutReal("mol",Mol->GetName(),"fluctx",
					fluct);
				Rg = Lat->Moment(phi,2,NULL,2);
				Rg /= N;
				Out->PutReal("mol",Mol->GetName(),"Rgy",
					sqrt(Rg));
				fluct = Lat->Moment(phi,4,NULL,2);
				fluct /= N;
				fluct -= Rg*Rg;
				fluct /= Rg*Rg;
				Out->PutReal("mol",Mol->GetName(),"flucty",
					fluct);
				Rg = Lat->Moment(phi,1,NULL,2); //y coordinate of CoM
				Rg /= N;
				Out->PutReal("mol",Mol->GetName(),"CoM", Rg);
				// output charge CoM
				double CoM = Mol->GetChargeCoM(1);
				Out->PutReal("mol",Mol->GetName(),
					"rel. CoM pos. charges", CoM-Rg);
				CoM = Mol->GetChargeCoM(-1);
				Out->PutReal("mol",Mol->GetName(),
					"rel. CoM neg. charges", CoM-Rg);
			}
		}
	}
	delete PinnedRange;
}
void
SF_System::GetOutputVesicle (Output* Out, Lattice* Lat) const {
	int M = Lat->GetTotalNumLayers();
	SF_Molecule* Mol = MolQ->GetMolecule("lipid");
	if (Mol == NULL) {
		Message(literal, "Don't know how to generate vesicle output, "
			"please define a molecule named 'lipid'");
		return;
	}
	Vector phi = Mol->GetPhi(total);
	double phibulk = Mol->GetPhiBulk();
	int z;
	for (z=1; z<=M; z++) {
		phi[z] -= phibulk;
	}
	LatticeRange* LayerZero = new LatticeRange1D(1,1,M);
	double R = Lat->MomentUnweighted(phi,1,LayerZero,-0.5)/Lat->MomentUnweighted(phi,0,LayerZero,-0.5);
	Out->PutReal("sys",name,"R",R);
	delete LayerZero;
	LatticeRange* LatRange = new LatticeRange1D(int(R)+1,int(R)+1,M);
	double d = 	Lat->MomentUnweighted(phi,2,LatRange,int(R)-R-0.5)/Lat->MomentUnweighted(phi,0,LatRange,int(R)-R-0.5);
	d = sqrt(d);
	Out->PutReal("sys",name,"d",d);
	for (z=1; z<=M; z++) {
		phi[z] += phibulk;
	}
	Vector pressureProfile = MolQ->GetGrandPotentialProfile();
	Out->PutReal("sys",name,"P0",
		Lat->MomentUnweighted(pressureProfile,0,LatRange,int(R)-R-0.5));
	Out->PutReal("sys",name,"P1",
		Lat->MomentUnweighted(pressureProfile,1,LatRange,int(R)-R-0.5));
	Out->PutReal("sys",name,"P2",
		Lat->MomentUnweighted(pressureProfile,2,LatRange,int(R)-R-0.5));
	delete LatRange;
}
void
SF_System::GetOutputTwoPhase (Output* Out, Text geometry, Lattice* Lat) const {
	if (SegQ->ReactionsPresent()) {
		Message(literal, "Don't know how to generate two phase output "
			"for reactions.");
		return;
	}
	int M = Lat->GetTotalNumLayers();
	double deltaP = (MolQ->GetGrandPotentialProfile())[M-1]-(MolQ->GetGrandPotentialProfile())[2];
	Out->PutReal("sys",name,"deltaP",deltaP);
	double R=0;
	for (int i=1; i<=MolQ->GetNumMolecules(); i++) {
		SF_Molecule* Mol = MolQ->GetMolecule(i);
		Text molName = Mol->GetName();
		double mu1 = MolQ->GetChemicalPotentialLayer(Mol,2);
		double muInf = MolQ->GetChemicalPotential(Mol);
		double N = Mol->GetChainLength();
		R = FindGibbsDevPlane(Mol->GetPhi(total), geometry,Lat);
		Out->PutReal("mol",molName,"mu(z=1)",mu1);
		Out->PutReal("mol",molName,"mu(z=1)+N*deltaP-mu(inf)",mu1+N*deltaP-muInf);
		Out->PutReal("mol",molName,"RGibbs",R);
	}
	double gammaG = CalcGammaGibbs(R,geometry,Lat);
	Out->PutReal("sys",name,"gammaG",gammaG);
	double dgdR = CalcGammaGibbs(R+1.0e-7,geometry,Lat) - CalcGammaGibbs(R-1.0e-7 ,geometry,Lat);
	dgdR /= 2.0e-7;
	double J;
	Vector fProfile = MolQ->GetFreeEnergyProfile();
	if (*geometry == "flat") {
		J=0;
	}
	if (*geometry == "sphere") {
		J=2/R;
	}
	if (*geometry == "cylinder") {
		J=1/R;
	}
	Out->PutReal("sys",name,"gammaGJ+dg/dR",gammaG*J+dgdR);
	Out->PutReal("sys",name,"gammaGJ",gammaG*J);
	Out->PutReal("sys",name,"dg/dR",dgdR);
	double f = 0;
	for (int z=1; z<=M; z++) {
		f += fProfile[z];
	}
	Out->PutReal("sys",name,"f",f);
}
double
SF_System::FindGibbsDevPlane(Vector phi, Text geometry, Lattice* Lat) const {
	int M = Lat->GetTotalNumLayers();
	double eps, epsold;
	double tol = Solve->GetTolerance()*50;
	double Rstep = 0.1;
	double Rs = (M-1)/2.0;;
	double LRs=0;
	double offset = Lat->GetLayerAdjustment();
	int z;
	epsold = 0.0;
	Boolean reverse = false;
	if (phi[M] > phi[1]) {
		reverse = true;
	}
	do {
		eps = 0.0;
		for (z=1; z<M; z++) {
			if (*geometry == "flat") {
				LRs = Rs - z;
			} else if (*geometry == "sphere") {
				LRs = (4.0/3.0) * M_PI * (pow(Rs-2+offset,3) - pow(z-2+offset,3));
			} else if (*geometry == "cylinder") {
				LRs = M_PI * (pow(Rs-2+offset,2) - pow(z-2+offset,2));
			}
			if (z<=Rs) {
				eps += (phi[1] - phi[z])*Lat->GetNumLatticeSites(z);
			} else if (z>Rs) {
				eps += (phi[M] - phi[z])*Lat->GetNumLatticeSites(z);
			}
			if ((Rs>z) && (Rs<(z+1))) {
				z++;
				eps += LRs * (phi[1] - phi[z]);
				eps += (Lat->GetNumLatticeSites(z) - LRs) * (phi[M] - phi[z]);
			}

		}
		if (eps * epsold < 0) Rstep /= 2;
	    if (eps > tol && !reverse) Rs -= Rstep;
		if (eps < -tol && !reverse) Rs += Rstep;
	    if (eps > tol && reverse) Rs += Rstep;
		if (eps < -tol && reverse) Rs -= Rstep;
		epsold = eps;
	}
	while (fabs(eps) > tol);
	return Rs-2+offset;
}
double
SF_System::CalcGammaGibbs(double Rs, Text geometry, Lattice* Lat) const {
	int M = Lat->GetTotalNumLayers();
	Vector GrandPotential = MolQ->GetGrandPotentialProfile();
	double GammaGibbs = 0;
	double LRs=0;
	double offset = Lat->GetLayerAdjustment();
	Rs = Rs + 2 - offset;
	for (int z=1; z<=M; z++) {
		if (*geometry == "flat") {
			LRs = Rs - z;
		} else if (*geometry == "sphere") {
			LRs = (4.0/3.0) * M_PI * (pow(Rs-2+offset,3) - pow(z-2+offset,3));
		} else if (*geometry == "cylinder") {
			LRs = M_PI * (pow(Rs-2+offset,2) - pow(z-2+offset,2));
		}
		if (z<=Rs) {
			GammaGibbs += (GrandPotential[z] - GrandPotential[2])*Lat->GetNumLatticeSites(z);
		}
		if (z>Rs) {
			GammaGibbs += (GrandPotential[z] - GrandPotential[M-1])*Lat->GetNumLatticeSites(z);
		}
		if ((Rs>z) && (Rs<(z+1))) {
			z++;
			GammaGibbs += LRs * (GrandPotential[z] - GrandPotential[2]);
			GammaGibbs += (Lat->GetNumLatticeSites(z) - LRs) * (GrandPotential[z] - GrandPotential[M-1]);
		}
	}
	if (*geometry == "sphere") {
		GammaGibbs /= 4*M_PI*(Rs-2+offset)*(Rs-2+offset);
	}
	if (*geometry == "cylinder") {
		GammaGibbs /= 2*M_PI*(Rs-2+offset);
	}
	return GammaGibbs;
}
void
SF_System::GetOutputEmulsion (Output* Out, Lattice* Lat) const {
	if (SegQ->ReactionsPresent()) {
		Message(literal, "Don't know how to generate emulsion output "
			"for reactions.");
		return;
	}
	int M = Lat->GetTotalNumLayers();
	double deltaP = -(MolQ->GetGrandPotentialProfile())[2];
	Out->PutReal("sys",name,"deltaP",deltaP);
	SF_Molecule* Mol = MolQ->GetMolecule("lipid");
	if (Mol == NULL) {
		Message(literal, "Don't know how to generate emulsion output, "
			"please define a molecule named 'lipid'");
		return;
	}
	Vector phiExc = Mol->GetPhi(total);
	double phibulk = Mol->GetPhiBulk();
	int z;
	for (z=1; z<=M; z++) {
		phiExc[z] -= phibulk;
	}
	LatticeRange* LayerZero = new LatticeRange1D(1,1,M);
	double R = Lat->MomentUnweighted(phiExc,1,LayerZero,-0.5)/Lat->MomentUnweighted(phiExc,0,LayerZero,-0.5);
	delete LayerZero;
	for (z=1; z<=M; z++) {
		phiExc[z] += phibulk;
	}
	LatticeRange* LatRange = new LatticeRange1D(int(R)+1,int(R)+1,M);
	Vector grandPotentialProf = MolQ->GetGrandPotentialProfile();
	Out->PutReal("sys",name,"P0",
		Lat->MomentUnweighted(grandPotentialProf,0,LatRange,int(R)-R-0.5) + deltaP*R);
	Out->PutReal("sys",name,"P1",
		Lat->MomentUnweighted(grandPotentialProf,1,LatRange,int(R)-R-0.5) - 0.5*deltaP*R*R);
	Out->PutReal("sys",name,"P2",
		Lat->MomentUnweighted(grandPotentialProf,2,LatRange,int(R)-R-0.5) + 1.0*deltaP*R*R*R/3.0);
	delete LatRange;
	for (int i=1; i<=MolQ->GetNumMolecules(); i++) {
		SF_Molecule* Mol = MolQ->GetMolecule(i);
		Text molName = Mol->GetName();
		double mu1 = MolQ->GetChemicalPotentialLayer(Mol,2);
		double muInf = MolQ->GetChemicalPotential(Mol);
		double N = Mol->GetChainLength();
		Out->PutReal("mol",molName,"mu(z=1)",mu1);
		Out->PutReal("mol",molName,"mu(z=1)+N*deltaP-mu(inf)",mu1+N*deltaP-muInf);
	}
}
void
SF_System::GetOutputIsolated(Output* Out, Lattice* Lat) const {
	SF_Molecule *Mol = NULL;
	for (int i=1; i<=MolQ->GetNumMolecules(); i++) {
		if (MolQ->GetMolecule(i)->GetChainLength() > 1) {
			Mol = MolQ->GetMolecule(i);
			Mol->ContactNumberDistr(Out);
		}
	}
	if (Mol == NULL) {
		Message(literal, "Don't know how to generate isolated output, "
			"please define a polymer");
		return;

	}
}
void
SF_System::GetOutputSpinodal(Output* Out, Lattice* Lat) const {
	SF_Molecule *Mol = NULL;
	for (int i=1; i<=MolQ->GetNumMolecules(); i++) {
		if (MolQ->GetMolecule(i)->GetChainLength() > 1) {
			Mol = MolQ->GetMolecule(i);
			Mol->ContactNumberDistr(Out,true);
		}
	}
	if (Mol == NULL) {
		Message(literal, "Don't know how to generate spinodal output, "
			"please define a polymer");
		return;

	}
}
void
SF_System::GetOutputDepletion(Output* Out, Lattice* Lat) const {
	int M = Lat->GetTotalNumLayers();
	SF_Molecule *Mol = NULL;
	for (int i=1; i<=MolQ->GetNumMolecules(); i++) {
		Mol = MolQ->GetMolecule(i);
		if (Mol->GetChainLength() > 1) {
			Vector phiExc = Mol->GetPhi(total);
			double phibulk = Mol->GetPhiBulk();
			for (int z=1; z<=M; z++) {
				phiExc[z] -= phibulk;
			}
			int numExtr = NumExtrema(phiExc,2,M);
			Vector place(1,numExtr);
			place = PlaceExtrema(phiExc,2,M,numExtr);
			Vector value(1,numExtr);
			value = ValueExtrema(phiExc,2,M,numExtr);
			int s;
			for (s=1; s<=numExtr-1; s++) {
				Text number=Blanks(100);
				number.Putint(s);
				number = Copy(number.Strip().Frontstrip());
				if (fabs(value[s]) > Solve->GetTolerance()) {
					Out->PutReal("mol",name,"placePhiExcExt"+number,place[s]-1.5);
					Out->PutReal("mol",name,"valuePhiExcExt"+number,fabs(value[s]));
				} else {
					Out->PutReal("mol",name,"placePhiExcExt"+number,0);
					Out->PutReal("mol",name,"valuePhiExcExt"+number,0);
				}
			}
			for (s=numExtr; s<=20; s++) {
				Text number=Blanks(100);
				number.Putint(s);
				number = Copy(number.Strip().Frontstrip());
				Out->PutReal("mol",name,"placePhiExcExt"+number,0);
				Out->PutReal("mol",name,"valuePhiExcExt"+number,0);
			}
		}
	}
}
void
SF_System::GetOutputPiA(Output* Out, Lattice* Lat) const {
	int M = Lat->GetTotalNumLayers();
	Vector pressure = MolQ->GetGrandPotentialProfile();
	Out->PutReal("sys",name,"pressure2",-pressure[2]);
	pressure[2] = 0;
	Lat->MultiplyWithLatticeSites(pressure);
	double press = 0;
	for(int z=1; z<=M; z++) {
		press -= pressure[z];
	}
	Lat->DivideByLatticeSites(pressure);
	Out->PutReal("sys",name,"pressure",press);
}
void
SF_System::GetOutputBend(Output* Out, Lattice* Lat) const {
	SF_Molecule *Mol = NULL;
	int M = Lat->GetTotalNumLayers();
	for (int i=1; i<=MolQ->GetNumMolecules(); i++) {
		if (MolQ->GetMolecule(i)->GetChainLength() > 1) {
			Mol = MolQ->GetMolecule(i);
		}
	}
	if (Mol == NULL) {
		Message(literal, "Don't know how to generate isolated output, "
			"please define a polymer");
		return;

	}
	Vector phi = Mol->GetPhi(total);
	double phiBulk = Mol->GetPhiBulk();
	int z;
	Vector G(1,M);
	for (z=1; z<=M; z++) {
		G[z] = 0.5*(phi[z] - phiBulk)*(phi[z] - phiBulk);
	}
	Out->PutReal("sys",name,"phidphi1",sqrt(phi[2])*(sqrt(phi[3])-sqrt(phi[2])));
	Out->PutReal("sys",name,"phidphi2",(sqrt(phi[2]) + sqrt(phi[3]))*(sqrt(phi[3])-sqrt(phi[2])));

	double zero=0;
	double first=0;
	double second=0;
	for (z=2; z<=M; z++) {
		zero += G[z];
		first += 0.5*((z-1)*(z-1) - (z-2)*(z-2))*G[z];
		second += (1.0/3.0)*((z-1)*(z-1)*(z-1) - (z-2)*(z-2)*(z-2))*G[z];
	}
	Out->PutReal("sys",name,"0thMG",zero);
	Out->PutReal("sys",name,"1stMG",first);
	Out->PutReal("sys",name,"2ndMG",second);
	GetOutputPressure(Out,Lat);
}
void
SF_System::GetOutputMultiState(Output* Out, Lattice* Lat) const {
	SF_Molecule *Mol = NULL;
	for (int i=1; i<=MolQ->GetNumMolecules(); i++) {
		if (MolQ->GetMolecule(i)->GetChainLength() > 1) {
			Mol = MolQ->GetMolecule(i);
			Mol->StateDistribution(Out);

		}
	}
	if (Mol == NULL) {
		Message(literal, "Don't know how to generate multistate output, "
			"please define a polymer");
		return;

	}
}

