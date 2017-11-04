#include <iostream>
#include "SF_MolList2ndO.h"

//double TEMPERATURE;

SF_MolList2ndO::SF_MolList2ndO(SF_SegmentList* SegQ_,
							   Lat2ndO* Lat_,
							   Input* MyInput_)
	: SF_MoleculeList(SegQ_,MyInput_) {

	Lat = Lat_;

	numMolecules = MyInput->GetNumNames("mol");
	Array<Text> names = MyInput->GetNames("mol");
	Text name;
	SF_MolStructure* Chain;
	Text composition;
	MoleculeType MolType;
	MoleculeQ.Dim(1,numMolecules);

	int i;
	bool LatRangeStiff;
	if (MyInput->ValueSet("lat",(MyInput->GetNames("lat"))[1],"stiff_range")) {
			Text range = MyInput->GetText("lat",(MyInput->GetNames("lat"))[1],"stiff_range");
			LatRangeStiff=true;
		} else LatRangeStiff=false;

	for (i=1; i<=numMolecules; i++) {
		name = names[i];
		composition = MyInput->GetText("mol",name,"composition");
		Chain = new SF_MolStructure(composition,name,SegQ,Lat,MyInput);
		composition = Chain->GetComposition();
		MolType = Chain->GetMoleculeType();
		delete Chain;
		switch (MolType) {
			case monomer:
	//Message(literal,"SF_Monomer2ndO called for " + name);
				MoleculeQ[i] = new SF_Monomer2ndO(name,SegQ,Lat,MyInput);
				break;
			case homopolymer:
				if (LatRangeStiff) {Message(fatal,"Range in stiffness only implemented for copolymers. ");} else
				MoleculeQ[i] = new SF_Homopol2ndO(name,SegQ,Lat,MyInput); //moet gefixed worden
				//MoleculeQ[i] = new SF_Copol2ndO(name,SegQ,Lat,MyInput);
				break;
			case copolymer:
    //Message(literal,"SF_Copol2ndO called for " + name);
				if (LatRangeStiff) {MoleculeQ[i] = new SF_Copol2ndO_stiff_range(name,SegQ,Lat,MyInput);} else
				{ MoleculeQ[i] = new SF_Copol2ndO(name,SegQ,Lat,MyInput); }
				break;
			case polydisperse:
//				MoleculeQ[i] = new SF_Polydisp2ndO(name,SegQ,Lat,MyInput);
				Message(fatal,"Polydisperse not available in 'SF_MolList2ndO::SF_MolList2ndO'");
				break;
   			case dendrimer:
//				MoleculeQ[i] = new SF_Dend2ndO(name,SegQ,Lat,MyInput);
				Message(fatal,"Dend2ndO not available in 'SF_MolList2ndO::SF_MolList2ndO'");
				break;
			case asymmetric_dendrimer:
//				MoleculeQ[i] = new SF_AsymDend2ndO(name,SegQ,Lat,MyInput);
				Message(fatal,"AsymDend2ndO not available in 'SF_MolList2ndO::SF_MolList2ndO'");
				break;
			case comb:
				MoleculeQ[i] = new SF_Comb2ndO(name,SegQ,Lat,MyInput);
				//Message(fatal,"Comb2ndO not available in 'SF_MolList2ndO::SF_MolList2ndO'");
				break;
			case branched:
//				MoleculeQ[i] = new SF_Branched2ndO(name,SegQ,Lat,MyInput);
				Message(fatal,"branched not available in 'SF_MolList2ndO::SF_MolList2ndO'");
				break;
//mmm end
			default:
				Message(fatal,"Programming error in 'SF_MolList2ndO::SF_MolList2ndO'");
				break;
		}
	}
	SF_Molecule* Mol=NULL;
	SF_Molecule* Neutralizer = GetNeutralizer();
	SF_Molecule* Solv = GetSolvent();
	if (Neutralizer != NULL) {
		double bulkCharge = 0;
		Vector ChargeBoundaries(1,2*Lat->GetNumGradients());
		for (i=1; i<=numMolecules; i++) {
			Mol = MoleculeQ[i];
			if (Mol != Neutralizer && Mol != Solv) {
				bulkCharge += Mol->GetBulkCharge()*Mol->GetPhiBulk()/Mol->GetChainLength();
				if (Lat->GetNamesBulkBoundaries().Upperbound() > 0) {
					for (int j=1; j<=2*Lat->GetNumGradients(); j++) {
						ChargeBoundaries[j] += (Mol->GetPhiBulkBoundaries())[j]*Mol->GetBulkCharge()/Mol->GetChainLength();
					}
				}
			}
		}
		for (int j=1; j<=2*Lat->GetNumGradients(); j++) {
			ChargeBoundaries[j] *= -Neutralizer->GetChainLength()/Neutralizer->GetBulkCharge();
		}
		Neutralizer->SetPhiBulk(-bulkCharge*Neutralizer->GetChainLength()/Neutralizer->GetBulkCharge());
		Mol->GetMolStructure()->SetPhiBulkBoundaries(ChargeBoundaries);
		double phiBulkSolv = 1;
		Vector phiBulkSolvBoundaries(1,2*Lat->GetNumGradients());
		for (i=1; i<=2*Lat->GetNumGradients(); i++) {
			phiBulkSolvBoundaries[i] = 1;
		}
		for (i=1; i<=numMolecules; i++) {
			Mol = MoleculeQ[i];
			if (Mol != Solv) {
				phiBulkSolv -= Mol->GetPhiBulk();
				if (Lat->GetNamesBulkBoundaries().Upperbound() > 0) {
					for (int j=1; j<=2*Lat->GetNumGradients(); j++) {
						phiBulkSolvBoundaries[j] -= (Mol->GetPhiBulkBoundaries())[j];
					}
				}
			}
		}
 	 	Solv->SetPhiBulk(phiBulkSolv);
		Solv->GetMolStructure()->SetPhiBulkBoundaries(phiBulkSolvBoundaries);
	} else {
		Mol = GetSolvent();
		double phiBulkSolv = 1;
		Vector phiBulkSolvBoundaries(1,2*Lat->GetNumGradients());
		for (i=1; i<=2*Lat->GetNumGradients(); i++) {
			phiBulkSolvBoundaries[i] = 1;
		}
		for (i=1; i<=numMolecules; i++) {
			Mol = MoleculeQ[i];
			if (Mol->GetFreedom() != solvent) {
				phiBulkSolv -= Mol->GetPhiBulk();
				if (Lat->GetNamesBulkBoundaries().Upperbound() > 0) {
					for (int j=1; j<=2*Lat->GetNumGradients(); j++) {
						phiBulkSolvBoundaries[j] -= (Mol->GetPhiBulkBoundaries())[j];
					}
				}
			}
		}
		Solv->SetPhiBulk(phiBulkSolv);
		Solv->GetMolStructure()->SetPhiBulkBoundaries(phiBulkSolvBoundaries);
	}
//	if (SegQ->Charged() && Lat->GetNumGradients() > 1 && Lat->GetNum) {
//		Message(fatal,"poisson equation not implemented in 2 gradients\n"
//			"Note: this error message is generated in SF_MolList2ndO constructor");
//	}
	double maxTheta = 0;
	int M = Lat->GetTotalNumLayers();
	for (int z=1; z<=M; z++) {
		if (Lat->WithinBoundaries(z)) {
			maxTheta += Lat->GetNumLatticeSites(z);
		}
	}
	maxTheta += Lat->GetNumExtraLatticeSites();
	double thetaTot = 0;
	for (i=1; i<=numMolecules; i++) {
		Mol = MoleculeQ[i];
		thetaTot += Mol->GetTheta();
	}
	if (thetaTot >= maxTheta) {
		cout << thetaTot << "   " << maxTheta << endl;
		Message(fatal,"Values for theta given add up to an amount greater than the available latticesize");
	}
	DeleteUnusedMolecules();
	SegQ->DeleteUnusedSegments();
}
SF_MolList2ndO::~SF_MolList2ndO() {
	SF_Molecule* Mol;
	for (int i=1; i<=numMolecules; i++) {
		Mol = MoleculeQ[i];
		delete Mol;
	}
}
void
SF_MolList2ndO::GetOutput(Output* Out) const {
	SF_Molecule* Mol;
	if (SegQ->ReactionsPresent()) {
		for (int i=1; i<=numMolecules; i++) {
			Mol = MoleculeQ[i];
			if (Mol->GetFreedom() == secondGeneration) {
				Out->PutReal("mol",Mol->GetName(),"sum n FH-MU free",GetNumMolTimesChemPot(i, unconstrained));
				Out->PutReal("mol",Mol->GetName(),"sum n FH-MU restr",GetNumMolTimesChemPot(i, constrained));
				Out->PutReal("mol",Mol->GetName(),"FH-MU free state1",GetChemicalPotential(i,unconstrained));
				Out->PutReal("mol",Mol->GetName(),"FH-MU restr state1",GetChemicalPotential(i,constrained));
			} else {
				Out->PutReal("mol",Mol->GetName(),"sum n FH-MU",GetNumMolTimesChemPot(i));
				if (Mol->GetChainLength() > 1) {
					Out->PutReal("mol",Mol->GetName(),"FH-MU state1",GetChemicalPotential(i));
				} else {
					SF_MolSegment* Seg = Mol->GetMolStructure()->GetDiffSegment(1);
					int numStates = Seg->GetNumStates();
					for (int j=1; j<=numStates; j++) {
						Out->PutReal("mol", Mol->GetName(), "FH-MU - " + Seg->GetState(j)->GetName(), GetChemicalPotential(i,total,j));
					}
				}
			}
			Out->PutReal("mol",Mol->GetName(),"entropy",GetEntropy(Mol));
			Mol->GetOutput(Out);
			Mol->GetLayerAnalysis(Out);
		}
	} else {
		for (int i=1; i<=numMolecules; i++) {
			Mol = MoleculeQ[i];
			if (Mol->GetFreedom() == secondGeneration) {
				Out->PutReal("mol",Mol->GetName(),"sum n FH-MU free",GetNumMolTimesChemPot(i, unconstrained));
				Out->PutReal("mol",Mol->GetName(),"sum n FH-MU restr",GetNumMolTimesChemPot(i, constrained));
				Out->PutReal("mol",Mol->GetName(),"FH-MU free",GetChemicalPotential(i,unconstrained));
				Out->PutReal("mol",Mol->GetName(),"FH-MU restr",GetChemicalPotential(i,constrained));
			} else {
				Out->PutReal("mol",Mol->GetName(),"sum n FH-MU",GetNumMolTimesChemPot(i));
				Out->PutReal("mol",Mol->GetName(),"FH-MU",GetChemicalPotential(i));
			}
			Out->PutReal("mol",Mol->GetName(),"entropy",GetEntropy(Mol));
			Mol->GetOutput(Out);
			Mol->GetLayerAnalysis(Out);
		}
	}
	int gradients=Lat->GetNumGradients();
	for (int i=1; i<=numMolecules; i++) {
		Mol = MoleculeQ[i];
		if (Mol->GetMoleculeType() == copolymer || Mol->GetMoleculeType() == homopolymer) {
			if (gradients==1) {
				Out->PutProfile("mol",Mol->GetName(),"phi_x",Mol->GetBondOrientation(x_dir));
				Out->PutProfile("mol",Mol->GetName(),"phi_yz",Mol->GetBondOrientation(yz_dir));
			} else {
				Out->PutProfile("mol",Mol->GetName(),"phi_x",Mol->GetBondOrientation(x_dir));
				Out->PutProfile("mol",Mol->GetName(),"phi_y",Mol->GetBondOrientation(y_dir));
				Out->PutProfile("mol",Mol->GetName(),"phi_z",Mol->GetBondOrientation(z_dir));
			}
		}
	}
}

void
SF_MolList2ndO::ReComputePhi() {
	Message(fatal,"ReComputePhi not implemented in SF_MolList2ndO");
}
void
SF_MolList2ndO::ComputePhi() {
	int i;
	SF_Molecule* Mol;
	SF_Molecule* Solv = GetSolvent();
	SF_Molecule* Neutralizer = GetNeutralizer();
	if (Neutralizer != NULL) {
		double bulkCharge = 0;
		for (i=1; i<=numMolecules; i++) {
			Mol = MoleculeQ[i];
			if (Mol != Neutralizer && Mol != Solv) {
				if (Mol->MuFixed()) { // only implemented for second generation
					Mol->SetLnNormRestricted(0);
					double lnCt = Mol->GetMuFixed() - GetChemicalPotential(i,constrained);
					Mol->SetLnNormRestricted(lnCt);
				}
				Mol->ComputePhi();
				bulkCharge += Mol->GetBulkCharge()*Mol->GetPhiBulk()/Mol->GetChainLength();
			}
		}
		Neutralizer->SetPhiBulk(-bulkCharge*Neutralizer->GetChainLength()/Neutralizer->GetBulkCharge());
		Neutralizer->ComputePhi();
		double phiBulkSolv = 1;
		for (i=1; i<=numMolecules; i++) {
			Mol = MoleculeQ[i];
			if (Mol != Solv) {
				phiBulkSolv -= Mol->GetPhiBulk();
			}
		}
 	 	Solv->SetPhiBulk(phiBulkSolv);
	    Solv->ComputePhi();
	} else {
		double phiBulkSolv = 1;
		for (i=1; i<=numMolecules; i++) {
			Mol = MoleculeQ[i];
			if (Mol != Solv) {
				if (Mol->MuFixed()) { // only implemented for second generation
					Mol->SetLnNormRestricted(0);
					double lnCt = Mol->GetMuFixed() - GetChemicalPotential(i,constrained);
					Mol->SetLnNormRestricted(lnCt);
				}
				Mol->ComputePhi();
				phiBulkSolv -= Mol->GetPhiBulk();
			}
		}
 	 	Solv->SetPhiBulk(phiBulkSolv);
	    Solv->ComputePhi();
	}
	SegQ->UpdatePhiBulk();
	SegQ->UpdatePhiStates();
}
void
SF_MolList2ndO::ComputePhi(double phiBulkSolv) {
	SF_Molecule* Mol;
	Mol = GetSolvent();
    Mol->SetPhiBulk(phiBulkSolv);
	for (int i=1; i<=numMolecules; i++) {
		Mol = MoleculeQ[i];
		if (Mol->MuFixed()) { // only implemented for second generation
			Mol->SetLnNormRestricted(0);
			double lnCt = Mol->GetMuFixed() - GetChemicalPotential(i,constrained);
			Mol->SetLnNormRestricted(lnCt);
		}
		Mol->ComputePhi();
	}
	SegQ->UpdatePhiBulk();
	SegQ->UpdatePhiStates();
}
void
SF_MolList2ndO::ComputePhi(double phiBulkSolv, double phiBulkNeutralizer) {
	SF_Molecule* Mol;
	Mol = GetSolvent();
    Mol->SetPhiBulk(phiBulkSolv);
	Mol = GetNeutralizer();
    Mol->SetPhiBulk(phiBulkNeutralizer);
	for (int i=1; i<=numMolecules; i++) {
		Mol = MoleculeQ[i];
		if (Mol->MuFixed()) { // only implemented for second generation
			Mol->SetLnNormRestricted(0);
			double lnCt = Mol->GetMuFixed() - GetChemicalPotential(i,constrained);
			Mol->SetLnNormRestricted(lnCt);
		}
		Mol->ComputePhi();
	}
	SegQ->UpdatePhiBulk();
	SegQ->UpdatePhiStates();
}
Boolean
SF_MolList2ndO::CheckMolecules(double accuracy) const {
	SF_Molecule* Mol;
	SF_MolStructure* Chain;
	SF_MolSegment* MolSeg;
	int M = Lat->GetTotalNumLayers();
	double theta;
	Vector phi;
	int numMolSeg;
	double thetaPart;
	double avLength;
	double avLengthPart;
	for (int i=1; i<=numMolecules; i++) {
		Mol = MoleculeQ[i];
		Chain = Mol->GetMolStructure();
		theta = Mol->ComputeTheta();
		avLength = Chain->GetAvLength();
		numMolSeg = Chain->GetNumDiffSegments();
		for (int j=1; j<=numMolSeg; j++) {
			MolSeg = Chain->GetDiffSegment(j);
			phi = MolSeg->GetPhi(total);
			Lat->SubtractBoundaries(phi);
			Lat->MultiplyWithLatticeSites(phi);
			thetaPart = 0;
			for (int z=1; z<=M; z++) {
				thetaPart += phi[z];
			}
			Lat->DivideByLatticeSites(phi);
			Lat->RestoreBoundaries(phi);
			Lat->SetBulkBoundaries(phi,MolSeg->GetPhiBulkBoundaries());
			avLengthPart = Chain->GetAvNumSegments(MolSeg);

//cout << "avLengthPart" << avLengthPart << "theta " << theta << "thetaPart" << thetaPart << "avLength" << avLength << endl;

			if (fabs((avLengthPart*theta)/(avLength*thetaPart) - 1) > accuracy*10) {
				Message(literal,"Some problem appears: \n"
				"Molecule '" + Mol->GetName() + "' not intact");

				cout << "AvLengthPart " << avLengthPart << endl;
				cout << "theta " << theta << endl;
				cout << "avLength " << avLength << endl;
				cout << "thetaPart " << thetaPart << endl;

				cout << "relative error : " << fabs((avLengthPart*theta)/(avLength*thetaPart) - 1) << endl;
				return false;
			}
		}
		if (Mol->GetFreedom() == fixedTheta) {
			double compTheta = Mol->ComputeTheta();
			if (Lat->GetNumExtraLatticeSites() != 0) {
				compTheta += Lat->GetNumExtraLatticeSites()*Mol->GetPhiBulk();
			}
			if (fabs(Mol->GetTheta()/compTheta - 1) > accuracy*10) {
				Message(literal,"Some problem appears: \n"
				"Molecule '" + Mol->GetName() + "' does not have the desired theta");
				cout << "desired theta : " << Mol->GetTheta() << "  actual theta : " << compTheta << endl;
				return true;
			}
		}
	}
	return false;
}
double
SF_MolList2ndO::GetChemicalPotential(int molNum, DensityPart densPart, int state) const {
	double value = 0;
	int i,j;
	SF_Molecule* Mol;
	SF_Molecule* Mol2;
	Mol = MoleculeQ[molNum];
	double N = Mol->GetChainLength();
	if (N > 1 && state != 1) {
		Message(fatal,"programming error in call to SF_MolList2ndO::GetChemicalPotential");
	}
	if (Mol->GetFreedom() == secondGeneration) {
		if (densPart == constrained) {
			value = Mol->GetLnNormRestricted() + log(N) + 1;
		} else {
			value = log(Mol->GetPhiBulk()) + 1;
		}
	} else {
		if (Mol->GetMolStructure()->SomeSegmentsPinned()) {
			value = Mol->GetLnNormRestricted() + log(N) + 1;
		} else if (Mol->GetMolStructure()->SomeSegmentsGrafted()) {
			SF_MolSegment* GraftedSegment = Mol->GetMolStructure()->GetGraftedSegment();
			double phiGrafted = GraftedSegment->GetPhiGrafted();
			value = Mol->GetLnNormRestricted() + log(N) - log(phiGrafted);
		} else {
			value = log(Mol->GetPhiBulk()) + 1;
		}
	}
	for (i=1; i<=numMolecules; i++) {
		Mol2 = MoleculeQ[i];
		value -= N * (Mol2->GetPhiBulk() / Mol2->GetChainLength());
	}
	SF_MolStructure* Chain;
	Chain = Mol->GetMolStructure();
	SF_MolSegment* MolSeg;
	SF_MolState* MolState;
	double dummy = 0;
	for (i=1; i<=Chain->GetNumDiffSegments(); i++) {
		MolSeg = Chain->GetDiffSegment(i);
		for (j=1; j<=MolSeg->GetNumStates(); j++) {
			MolState = MolSeg->GetState(j);
			dummy += MolSeg->GetPhiRef() * MolState->GetAlphaBulk() * (SegQ->ChemIntBulk(MolState));
		}
	}
	value += N * (dummy - Chain->ChemIntRef());
	dummy = 0;
	SF_State* State;
	for (i=1; i<=SegQ->GetNumStates(); i++) {
		State = SegQ->GetState(i);
		dummy += State->GetPhiBulk() * SegQ->ChemIntBulk(State);
	}
	value -= 0.5*Mol->GetChainLength() * dummy;
	if (SegQ->ReactionsPresent()) {
		for (i=1; i<=Chain->GetNumDiffSegments(); i++) {
			dummy = 0;
			MolSeg = Chain->GetDiffSegment(i);
			MolState = MolSeg->GetState(state);
			dummy += log(MolState->GetAlphaBulk()) + MolState->GetInternalFreeEnergy();
			dummy -= MolSeg->GetState(1)->GetInternalFreeEnergy(); // state1 is reference
			dummy *= Chain->GetAvNumSegments(MolSeg);
			value += dummy;
		}
	}
	return value;
}
double
SF_MolList2ndO::GetChemicalPotential(const SF_Molecule* Mol, DensityPart densPart, int state) const {
	for (int i=1; i<=numMolecules; i++) {
		if (MoleculeQ[i] == Mol) return GetChemicalPotential(i,densPart,state);
	}
	Message(fatal,"error in call to SF_MolList2ndO::GetChemicalPotential"
		"(const SF_Molecule* Mol)");
	return 0; // never get here
}
double
SF_MolList2ndO::GetNumMolTimesChemPot(int molNum, DensityPart densPart) const {
	double value = 0;
	int i,j;
	SF_Molecule* Mol;
	SF_Molecule* Mol2;
	Mol = MoleculeQ[molNum];
	double N = Mol->GetChainLength();
	double theta = 0;
	if (Mol->GetFreedom() == secondGeneration) {
		if (densPart == constrained) {
			value = Mol->GetLnNormRestricted() + log(N) + 1;
			theta = Mol->GetTheta();
		} else {
			value = log(Mol->GetPhiBulk()) + 1;
			theta = Mol->ComputeTheta() - Mol->GetTheta();
		}
	} else {
		theta = Mol->ComputeTheta();
		if (Mol->GetMolStructure()->SomeSegmentsPinned()) {
			value = Mol->GetLnNormRestricted() + log(N) + 1;
		} else if (Mol->GetMolStructure()->SomeSegmentsGrafted()) {
			SF_MolSegment* GraftedSegment = Mol->GetMolStructure()->GetGraftedSegment();
			double phiGrafted = GraftedSegment->GetPhiGrafted();
			value = Mol->GetLnNormRestricted() + log(N) - log(phiGrafted);
		} else {
			value = log(Mol->GetPhiBulk()) + 1;
		}
	}
	value *= theta/N;
	double dummy = 0;
	for (i=1; i<=numMolecules; i++) {
		Mol2 = MoleculeQ[i];
		dummy -= (Mol2->GetPhiBulk() / Mol2->GetChainLength());
	}
	value += theta*dummy;
//	cout << Mol->GetName().MainC() << '\t' << value << '\t' << dummy << endl;
	SF_MolStructure* Chain;
	Chain = Mol->GetMolStructure();
	SF_MolSegment* MolSeg;
	SF_MolState* MolState;
	dummy = 0;
	for (i=1; i<=Chain->GetNumDiffSegments(); i++) {
		MolSeg = Chain->GetDiffSegment(i);
		for (j=1; j<=MolSeg->GetNumStates(); j++) {
			MolState = MolSeg->GetState(j);
			dummy += MolState->GetTheta(densPart) * (SegQ->ChemIntBulk(MolState));
		}
	}
	value += dummy - theta*Chain->ChemIntRef();
//	cout << Mol->GetName().MainC() << '\t' << value << '\t' << dummy << endl;
	dummy = 0;
	SF_State* State;
	for (i=1; i<=SegQ->GetNumStates(); i++) {
		State = SegQ->GetState(i);
		dummy += State->GetPhiBulk() * SegQ->ChemIntBulk(State);
	}
	value -= 0.5 * theta * dummy;
//	cout << Mol->GetName().MainC() << '\t' << value << '\t' << dummy << endl;
	if (SegQ->ReactionsPresent()) {
		for (i=1; i<=Chain->GetNumDiffSegments(); i++) {
			dummy = 0;
			MolSeg = Chain->GetDiffSegment(i);
			double intFreeEnergyRef = MolSeg->GetState(1)->GetInternalFreeEnergy();
			for (j=1; j<=MolSeg->GetNumStates(); j++) {
				MolState = MolSeg->GetState(j);
				dummy += MolState->GetTheta(densPart) * (log(MolState->GetAlphaBulk())
					+ MolState->GetInternalFreeEnergy() - intFreeEnergyRef);
			}
			value += dummy;
		}
	}
	return value;
}
double
SF_MolList2ndO::GetNumMolTimesChemPot(const SF_Molecule* Mol, DensityPart densPart) const {
	for (int i=1; i<=numMolecules; i++) {
		if (MoleculeQ[i] == Mol) return GetNumMolTimesChemPot(i,densPart);
	}
	Message(fatal,"error in call to SF_MolList2ndO::GetNumMolTimesChemPot"
		"(const SF_Molecule* Mol)");
	return 0; // never get here
}
double
SF_MolList2ndO::GetChemicalPotentialLayer(int molNum, int z) const {
	double value = 0;
	int i,j;
	SF_Molecule* Mol;
	SF_Molecule* Mol2;
	Mol = MoleculeQ[molNum];
	Vector phi = Mol->GetPhi(total);
	if (phi[z] == 0) return 0;
	value = log(phi[z]) + 1;
	for (i=1; i<=numMolecules; i++) {
		Mol2 = MoleculeQ[i];
		Vector phi2 = Mol2->GetPhi(total);
		value -= Mol->GetChainLength() * (phi2[z] / Mol2->GetChainLength());
	}
	SF_MolStructure* Chain;
	Chain = Mol->GetMolStructure();
	value -= Mol->GetChainLength() * Chain->ChemIntRef();
	SF_MolSegment* MolSeg;
	SF_MolState* MolState;
	double dummy = 0;
	for (i=1; i<=Chain->GetNumDiffSegments(); i++) {
		MolSeg = Chain->GetDiffSegment(i);
		for (j=1; j<=MolSeg->GetNumStates(); j++) {
			MolState = MolSeg->GetState(j);
			dummy += MolState->GetPhiRef() * SegQ->ChemInt(MolState,z);
		}
	}
	value += Mol->GetChainLength() * dummy;
	dummy = 0;
	SF_State* State;
	for (i=1; i<=SegQ->GetNumStates(); i++) {
		State = SegQ->GetState(i);
		phi = State->GetPhi();
		dummy += phi[z] * SegQ->ChemInt(State,z);
	}
	value -= 0.5 * Mol->GetChainLength() * dummy;
/*	if (SegQ->ReactionsPresent()) {
		for (i=1; i<=Chain->GetNumDiffSegments(); i++) {
			dummy = 0;
			MolSeg = Chain->GetDiffSegment(i);
			for (j=1; j<=MolSeg->GetNumStates(); j++) {
				MolState = MolSeg->GetState(j);
				dummy += MolState->GetAlphaBulk() * log(MolState->GetAlphaBulk());
			}
			dummy *= Chain->GetAvNumSegments(MolSeg);
			value += dummy;
		}
	}*/
	return value;
}
double
SF_MolList2ndO::GetChemicalPotentialLayer(const SF_Molecule* Mol, int z) const {
	for (int i=1; i<=numMolecules; i++) {
		if (MoleculeQ[i] == Mol) return GetChemicalPotentialLayer(i,z);
	}
	Message(fatal,"error in call to SF_MolList2ndO::GetChemicalPotential"
		"(const SF_Molecule* Mol)");
	return 0; // never get here
}
double
SF_MolList2ndO::GetFreeEnergy() const {
	double value = 0;
	int M = Lat->GetTotalNumLayers();
	Vector freeEnergy = GetFreeEnergyProfile();
	Lat->MultiplyWithLatticeSites(freeEnergy);
	for (int z=1; z<=M; z++) {
		value += freeEnergy[z];
	}
	if (Lat->GetNumExtraLatticeSites() != 0) {
		double extraLatticeSites = Lat->GetNumExtraLatticeSites();
		int i;
		for (i=1; i<=numMolecules; i++) {
			SF_Molecule* Mol = MoleculeQ[i];
			double phiBulk = Mol->GetPhiBulk();
			if (phiBulk > 0) {
				value += extraLatticeSites*(phiBulk/Mol->GetChainLength())*log(phiBulk);
				double chemIntRef = extraLatticeSites*Mol->GetMolStructure()->ChemIntRef();
				value -= phiBulk*chemIntRef;
			}
		}
		int numStates = SegQ->GetNumStates();
		for (i=1; i<=numStates; i++) {
			SF_State* State = SegQ->GetState(i);
			double phiBulk = State->GetPhiBulk();
			value += 0.5*extraLatticeSites*phiBulk*SegQ->ChemIntBulk(State);
			value += phiBulk*extraLatticeSites*log(State->GetAlphaBulk());
		}
	}
	return value;
}
double
SF_MolList2ndO::GetGrandPotential() const {
	Vector grandPotential = GetGrandPotentialProfile();
	Lat->MultiplyWithLatticeSites(grandPotential);
	int M = Lat->GetTotalNumLayers();
	double value = 0;
	for (int z=1; z<=M; z++) {
		value += grandPotential[z];
	}
	return value;
}

double
SF_MolList2ndO::GetGrandPotentialExcess() const {
	return 0;
}

double
SF_MolList2ndO::GetGrandPotential2() const {
	double grandPotential = GetFreeEnergy();
	for	(int i=1; i<=numMolecules; i++)	{
		SF_Molecule* Mol = MoleculeQ[i];
		if (Mol->GetFreedom() == secondGeneration) {
			grandPotential -= GetNumMolTimesChemPot(i,unconstrained);
			grandPotential -= GetNumMolTimesChemPot(i,constrained);
		} else {
			grandPotential -= GetNumMolTimesChemPot(i);
		}
	}
	return grandPotential;
}
double
SF_MolList2ndO::GetSoumiFreeEnergy() const	{
	return 0;
}
double
SF_MolList2ndO::GetExcessFreeEnergy() const	{
	Vector Fexc	= GetGrandPotentialProfile();
	Lat->MultiplyWithLatticeSites(Fexc);
	int	M =	Lat->GetTotalNumLayers();
	double value = 0;
	int z;
	for	(z=1; z<=M;	z++) {
		value += Fexc[z];
	}
	SF_Molecule* Mol;
	double thetaExcess;
	for	(int i=1; i<=numMolecules; i++)	{
		Mol	= MoleculeQ[i];
		double phiBulk = Mol->GetPhiBulk();
		if (Mol->GetFreedom() == secondGeneration) {
			int M = Lat->GetTotalNumLayers();
			Vector thetaExc(1,M);
			Vector phi = Mol->GetPhi(unconstrained);
			thetaExcess = 0;
			for (z=1; z<=M; z++) {
				if (phi[z] > 0) {
					thetaExc[z] = phi[z] - phiBulk;
				}
			}
			Lat->SubtractBoundaries(thetaExc);
			Lat->MultiplyWithLatticeSites(thetaExc);
			thetaExcess = 0;
			for (z=1; z<=M; z++) {
				thetaExcess += thetaExc[z];
			}
			value += thetaExcess*GetChemicalPotential(Mol,unconstrained)/Mol->GetChainLength();
			value += Mol->GetTheta()*GetChemicalPotential(Mol,constrained)/Mol->GetChainLength();
		} else {
			thetaExcess	= Mol->ComputeThetaExcess();
			value += thetaExcess*GetChemicalPotential(Mol)/Mol->GetChainLength();
		}
	}
	return value;
}
double
SF_MolList2ndO::GetEntropy(SF_Molecule* Mol2) const {
	double value = 0;
	int M = Lat->GetTotalNumLayers();
	Vector entropy = GetEntropyProfile(Mol2);
	Lat->MultiplyWithLatticeSites(entropy);
	for (int z=1; z<=M; z++) {
		value += entropy[z];
	}
	return value;
}
double
SF_MolList2ndO::GetContactInteractions() const {
	double value = 0;
	int M = Lat->GetTotalNumLayers();
	Vector contactInt = GetContactInteractionsProfile();
	Lat->MultiplyWithLatticeSites(contactInt);
	for (int z=1; z<=M; z++) {
		value += contactInt[z];
	}
	return value;
}
double
SF_MolList2ndO::GetInternalFreeEnergy() const {
	if (!SegQ->ReactionsPresent()) {
		return 0;
	}
	double value = 0;
	int M = Lat->GetTotalNumLayers();
	Vector intFreeEnergy = GetInternalFreeEnergyProfile();
	Lat->MultiplyWithLatticeSites(intFreeEnergy);
	for (int z=1; z<=M; z++) {
		value += intFreeEnergy[z];
	}
	return value;
}
double
SF_MolList2ndO::GetElectricInteractions() const {
	if (!SegQ->Charged()) {
		return 0;
	}
	double value = 0;
	int M = Lat->GetTotalNumLayers();
	Vector electricInt = GetElectricInteractionsProfile();
	Lat->MultiplyWithLatticeSites(electricInt);
	for (int z=1; z<=M; z++) {
		value += electricInt[z];
	}
	return value;
}
Vector
SF_MolList2ndO::GetFreeEnergyProfile() const {

	int M = Lat->GetTotalNumLayers();
	Vector freeEnergy(1,M);

	Vector entropy = GetEntropyProfile();

	Vector contactInt = GetContactInteractionsProfile();

	Vector elecInt = GetElectricInteractionsProfile();
;
	Vector intFreeEnergy = GetInternalFreeEnergyProfile();

	for (int z=1; z<=M; z++) {
		freeEnergy[z] = -entropy[z] + contactInt[z] + elecInt[z] + intFreeEnergy[z];
	}
	return freeEnergy;
}
Vector
SF_MolList2ndO::GetEntropyProfile(SF_Molecule* Mol2) const {
	SF_Molecule* Mol;
	double N;
	int i,j,z;
	int M = Lat->GetTotalNumLayers();
	Vector entropy(1,M);
	Vector phi;

	for (i=1; i<=numMolecules; i++) {
		if (Mol2 != 0) {
			Mol = Mol2;
			i = numMolecules;
		} else {
			Mol = MoleculeQ[i];
		}
		N = Mol->GetChainLength();
		if (Mol->GetFreedom() == secondGeneration) {
			if (Mol->GetPhiBulk() > 0) {
				phi = Mol->GetPhi(unconstrained);
				for (z=1; z<=M; z++) entropy[z] -= (phi[z]/N)*(log(N) + Mol->GetLnNormFree());
			}
			phi = Mol->GetPhi(constrained);
			for (z=1; z<=M; z++) entropy[z] -= (phi[z]/N)*(log(N) + Mol->GetLnNormRestricted());
		} else if (Mol->GetFreedom() == thirdGeneration) {
			if (Mol->GetPhiBulk() > 0) {
				phi = Mol->GetPhi(unconstrained);
				for (z=1; z<=M; z++) entropy[z] -= (phi[z]/N)*(log(N) + Mol->GetLnNormFree());
			}
			phi = Mol->GetPhi(constrained);
			for (z=1; z<=M; z++) entropy[z] -= (phi[z]/N)*(log(N) + Mol->GetLnNormRestricted());
			SF_MolStructure* Chain = Mol->GetMolStructure();
			if (Chain->GetNumDiffSegments() > 1) {
				Message(fatal, "programming error in SF_MolList2ndO::GetFreeEnergyProfileOld()");
			}
			for (j=1; j<=Chain->GetNumDiffSegments(); j++) {
				SF_MolSegment* MolSeg = Chain->GetDiffSegment(j);
				phi = MolSeg->GetPhi(renorm);
				Vector norm = Mol->GetRenorm(j);
				for (z=1; z<=M; z++) entropy[z] -= (phi[z]/N)*log(N*norm[(z-1)/Lat->GetNumLayers(1)+1]);
			}
		} else if (Mol->GetFreedom() == fixedTheta || Mol->GetFreedom() == rangeRestricted) {
			phi = Mol->GetPhi(total);
			for (z=1; z<=M; z++) entropy[z] -= (phi[z]/N)*(log(N) + Mol->GetLnNormRestricted());
		} else {
			phi = Mol->GetPhi(total);
			for (z=1; z<=M; z++) entropy[z] -= (phi[z]/N)*(log(N) + Mol->GetLnNormFree());
		}
		if (Mol->GetMolStructure()->SomeSegmentsGrafted()) {
			phi = Mol->GetPhi(total);
			double phiGrafted = Mol->GetMolStructure()->GetGraftedSegment()->GetPhiGrafted();
			for (z=1; z<=M; z++) entropy[z] += (phi[z]/N)*log(phiGrafted);
		}
	}
	SF_MolStructure* Chain;
	SF_MolSegment* MolSeg;
	Vector G;
	for (i=1; i<=numMolecules; i++) {
		if (Mol2 != 0) {
			Mol = Mol2;
			i = numMolecules;
		} else {
			Mol = MoleculeQ[i];
		}
		Chain = Mol->GetMolStructure();
		for (j=1; j<=Chain->GetNumDiffSegments(); j++) {
			MolSeg = Chain->GetDiffSegment(j);
			G = MolSeg->GetSWF();
			phi = MolSeg->GetPhi(total);
			for (z=1; z<=M; z++) {
				if (G[z] > 0) {
					entropy[z] -= phi[z]*log(G[z]);
				}
			}
		}
	}
	Lat->SubtractBoundaries(entropy);
	return entropy;
}

Vector
SF_MolList2ndO::GetContactInteractionsProfile() const {
	SF_Molecule* Mol;
	int i,z;
	int M = Lat->GetTotalNumLayers();
	Vector contactInt(1,M);
	Vector phi;
	SF_MolStructure* Chain;
	for (i=1; i<=numMolecules; i++) {
		Mol = MoleculeQ[i];
		Chain = Mol->GetMolStructure();
		double chemIntRef = Chain->ChemIntRef();
		phi = Mol->GetPhi(total);
		for (z=1; z<=M; z++) {
			contactInt[z] -= phi[z]*chemIntRef;
		}
	}
	SF_State* State;
	for (i=1; i<=SegQ->GetNumStates(); i++) {
		State = SegQ->GetState(i);
		phi = State->GetPhi();
		for (z=1; z<=M; z++) {
			contactInt[z] += 0.5*phi[z]*SegQ->ChemInt(State,z);
		}
	}
	Lat->SubtractBoundaries(contactInt);
	return contactInt;
}

Vector
SF_MolList2ndO::GetElectricInteractionsProfile() const {
	int M = Lat->GetTotalNumLayers();
	Vector electricInt(1,M);
	if (SegQ->Charged()) {
		Vector charge = SegQ->GetCharge();
		Lat->DivideByLatticeSites(charge);
		Vector psi = SegQ->GetElectricPotential();
		for (int z=1; z<=M; z++) {
			electricInt[z] += 0.5*psi[z]*charge[z]*ELEM_CHARGE/(BOLTZMANN*TEMPERATURE);
		}
	}
//	if (SegQ->Charged()) {
//		Vector charge = SegQ->GetCharge();
//		Lat->DivideByLatticeSites(charge);
//		Vector epsilon = SegQ->GetAverageEpsilon();
//		Vector psi = SegQ->GetElectricPotential();
//		double preFactor = (EPS0*Lat->GetSiteDistance())/(BOLTZMANN*TEMPERATURE);
//		for (int z=2; z<M; z++) {
////			electricInt[z] += 0.5*psi[z]*charge[z]*ELEM_CHARGE/(BOLTZMANN*TEMPERATURE);
////			electricInt[z] -= 0.5*preFactor*epsilon[z]*psi[z]*(psi[z-1]-2*psi[z]+psi[z+1]);
////			electricInt[z] += 0.5*preFactor*0.25*((epsilon[z-1] + epsilon[z])*(psi[z]-psi[z-1])*(psi[z]-psi[z-1]) + (epsilon[z] + epsilon[z+1])*(psi[z+1]-psi[z])*(psi[z+1]-psi[z]));
////			electricInt[z] += 0.5*preFactor*2*epsilon[z]* ( (epsilon[z-1]*(psi[z]-psi[z-1])/(epsilon[z]+epsilon[z-1]))*(epsilon[z-1]*(psi[z]-psi[z-1])/(epsilon[z]+epsilon[z-1]))
////														  + (epsilon[z+1]*(psi[z+1]-psi[z])/(epsilon[z]+epsilon[z+1]))*(epsilon[z+1]*(psi[z+1]-psi[z])/(epsilon[z]+epsilon[z+1])));
//		}
//	}
	Lat->SubtractBoundaries(electricInt);
	return electricInt;
}
Vector
SF_MolList2ndO::GetInternalFreeEnergyProfile() const {
	int M = Lat->GetTotalNumLayers();
	Vector internalFreeEnergy(1,M);
	if (SegQ->ReactionsPresent()) {
		for (int i=1; i<=numMolecules; i++) {
			SF_Molecule* Mol = MoleculeQ[i];
			SF_MolStructure* Chain = Mol->GetMolStructure();
			for (int j=1; j<=Chain->GetNumDiffSegments(); j++) {
				SF_MolSegment* MolSeg = Chain->GetDiffSegment(j);
				double intFreeEnergyRef = MolSeg->GetState(1)->GetInternalFreeEnergy();
				for (int k=1; k<=MolSeg->GetNumStates(); k++) {
					SF_MolState* MolState = MolSeg->GetState(k);
					Vector phi = MolState->GetPhi(total);
					double intFreeEnergy = MolState->GetInternalFreeEnergy();
					Vector alpha = MolSeg->GetAlpha(k);
					for (int z=1; z<=M; z++) {
						if (alpha[z] > 0) {
							internalFreeEnergy[z] += phi[z]*(log(alpha[z])+intFreeEnergy);
						}
						internalFreeEnergy[z] -= phi[z]*intFreeEnergyRef;
					}
 				}
			}
		}
	}
	Lat->SubtractBoundaries(internalFreeEnergy);
	return internalFreeEnergy;
}
Vector
SF_MolList2ndO::GetGrandPotentialProfile() const {
	// not implemented for off-equilibrium
	SF_Molecule* Mol;
	int M = Lat->GetTotalNumLayers();
	Boolean charged = SegQ->Charged();
	Vector psi;
	if (charged) {
		psi = SegQ->GetElectricPotential();
	}
	Vector grandPotential(1,M);
	int i,z;
	for (i=1; i<=numMolecules; i++) {
		Mol = MoleculeQ[i];
		double N = Mol->GetChainLength();
		double phiBulk = Mol->GetPhiBulk();
		Vector phi = Mol->GetPhi(total);
		if (!Mol->GetMolStructure()->SomeSegmentsGrafted()) {
			for (z=1; z<=M; z++) {
				if (phi[z] > 0) {
					grandPotential[z] -= (phi[z] - phiBulk)/N;
				}
			}
		}
	}
	int numStates = SegQ->GetNumStates();
	for (i=1; i<=numStates; i++) {
		SF_State* State = SegQ->GetState(i);
		if (State->GetFreedom() != frozen) {
			Vector SWF = State->GetSWF();
			Vector phi = State->GetPhi();
			double uPrime;
			double chemIntBulk = SegQ->ChemIntBulk(State);
			for (z=1; z<=M; z++) {
				uPrime = 0;
				if (SWF[z] > 0) {
					uPrime = -log(SWF[z]);
					uPrime -= SegQ->ChemInt(State,z);
					uPrime += chemIntBulk;
					if (charged) {
						uPrime -= State->GetValence()*psi[z]*
							ELEM_CHARGE/(BOLTZMANN*TEMPERATURE);
					}
					if (State->ExternalPotential()) {
						uPrime -= (State->GetExternalPotential())[z];
					}
				}
				grandPotential[z] -= uPrime*phi[z];
			}
		}
	}
	Matrix chi = SegQ->GetChiMatrix();
	for (i=1; i<=numStates; i++) {
		SF_State* State1 = SegQ->GetState(i);
		if (State1->GetFreedom() != frozen) {
			Vector phi1 = State1->GetPhi();
			for (int j=1; j<=numStates; j++) {
				SF_State* State2 = SegQ->GetState(j);
				if (State2->GetFreedom() != frozen) {
					//Vector phi2 = State2->GetPhi();
					double phiBulk1 = State1->GetPhiBulk();
					double phiBulk2 = State2->GetPhiBulk();
					for (z=1; z<=M; z++) {
						grandPotential[z] -= 0.5*chi[i][j]*phi1[z]*State2->GetRHO(z);
						//grandPotential[z] -= 0.5*chi[i][j]*phi1[z]*Lat->SideFraction(phi2,z);
						if (phi1[z] > 0) {
							grandPotential[z] += 0.5*chi[i][j]*phiBulk1*phiBulk2;
						}
					}
				}
			}
		}
	}
	if (charged) {
//		Vector charge = SegQ->GetCharge();
//		Lat->DivideByLatticeSites(charge);
		Vector electricInteraction = GetElectricInteractionsProfile();
		for (int z=1; z<=M; z++) {
			Boolean surface = false;
			SF_State* State;
			LatticeRange* LatRange;
			int numStates = SegQ->GetNumStates();
			for (int i=1; i<=numStates; i++) {
				State = SegQ->GetState(i);
				if (State->GetFreedom() == frozen) {
					LatRange = State->GetLatRange();
					if (LatRange->InRange(z)) {
						surface = true;
					}
				}
			}
			if (surface) {
				grandPotential[z] += electricInteraction[z];
			} else {
				grandPotential[z] -= electricInteraction[z];
			}
		}
	}
	Lat->SubtractBoundaries(grandPotential);
	return grandPotential;
}
Vector
SF_MolList2ndO::GetGrandPotentialProfile2() const {
	SF_Molecule* Mol;
	int M = Lat->GetTotalNumLayers();
	Vector grandPotential = GetFreeEnergyProfile();
	for (int i=1; i<=numMolecules; i++) {
		Mol = MoleculeQ[i];
		double N = Mol->GetChainLength();
		if (Mol->GetFreedom() == secondGeneration) {
			double chemPot = GetChemicalPotential(i,unconstrained);
			Vector phiFree = Mol->GetPhi(unconstrained);
			Lat->SubtractBoundaries(phiFree);
			int z;
			for (z=1; z<=M; z++) {
				grandPotential[z] -= phiFree[z]*chemPot/N;
			}
			Lat->RestoreBoundaries(phiFree);
			Lat->SetBulkBoundaries(phiFree,Mol->GetPhiBulkBoundaries());
			chemPot = GetChemicalPotential(i,constrained);
			Vector phiRestr = Mol->GetPhi(constrained);
			Lat->SubtractBoundaries(phiRestr);
			for (z=1; z<=M; z++) {
				grandPotential[z] -= phiRestr[z]*chemPot/N;
			}
			Lat->RestoreBoundaries(phiRestr);
		} else if (Mol->GetFreedom() == thirdGeneration) {
			Message(implementation,"SF_MolList2ndO::GetPressureProfile for thirdgeneration");
		} else {
			double chemPot = GetChemicalPotential(i);
			Vector phi = Mol->GetPhi(total);
			Lat->SubtractBoundaries(phi);
			for (int z=1; z<=M; z++) {
				grandPotential[z] -= phi[z]*chemPot/N;
			}
			Lat->RestoreBoundaries(phi);
			Lat->SetBulkBoundaries(phi,Mol->GetPhiBulkBoundaries());
		}
	}
	return grandPotential;
}
Vector
SF_MolList2ndO::GetExcessFreeEnergyProfile() const {
	Vector Fexc;
	int z;
	if (SegQ->ReactionsPresent()) {
		Fexc = GetGrandPotentialProfile2();
	} else {
		Fexc = GetGrandPotentialProfile();
	}
	int	M =	Lat->GetTotalNumLayers();
	SF_Molecule* Mol;
	for	(int i=1; i<=numMolecules; i++)	{
		Mol	= MoleculeQ[i];
		double phiBulk = Mol->GetPhiBulk();
		if (Mol->GetFreedom() == thirdGeneration) {
			Message(implementation, "GetExcessFreeEnergyProfile");
		} else if (Mol->GetFreedom() == secondGeneration) {
			int M = Lat->GetTotalNumLayers();
			Vector phiExc(1,M);
			Vector phiFree = Mol->GetPhi(unconstrained);
			Vector phiRestr = Mol->GetPhi(constrained);
			for (z=1; z<=M; z++) {
				if (phiFree[z] > 0) {
					phiExc[z] = phiFree[z] - phiBulk;
				}
			}
			Lat->SubtractBoundaries(phiExc);
			Lat->SubtractBoundaries(phiRestr);
			double muUnconstr = GetChemicalPotential(Mol,unconstrained);
			double muConstr = GetChemicalPotential(Mol,constrained);
			double N = Mol->GetChainLength();
			for (z=1; z<=M; z++) {
				Fexc[z] += phiExc[z]*muUnconstr/N;
				Fexc[z] += phiRestr[z]*muConstr/N;
			}
			Lat->RestoreBoundaries(phiRestr);
		} else {
			Vector phi = Mol->GetPhi(total);
			Vector phiExc(1,M);
			for (z=1; z<=M; z++) {
				if (phi[z] > 0) {
					phiExc[z] = phi[z] - phiBulk;
				}
			}
			Lat->SubtractBoundaries(phiExc);
			double mu = GetChemicalPotential(Mol);
			double N = Mol->GetChainLength();
			for (z=1; z<=M; z++) {
				Fexc[z] += phiExc[z]*mu/N;
			}
		}
	}
	return Fexc;
}




