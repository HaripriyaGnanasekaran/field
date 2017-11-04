#include "SF_Molecule.h"


SF_Molecule::SF_Molecule(Text name_,
						 SF_SegmentList* SegQ_,
						 Lattice* Lat_,
						 Input* MyInput_) {
	name = name_;
	LatBase = Lat_;
	SegQ = SegQ_;
	MyInput = MyInput_;
	muFixed = false;
	M = LatBase->GetTotalNumLayers();
	gradients = LatBase->GetNumGradients();
	Text AliasValue;


	if (MyInput->GetNumNames("alias") > 0) {
		numAlias = MyInput->GetNumNames("alias");
		AliasNames = MyInput->GetNames("alias");
		for (int i=1; i<=numAlias; i++) {
			Array<Text> parameters;
			parameters.Dim(1,1);
			parameters[1] = "value";
			MyInput->CheckParameterNames("alias",AliasNames[i],parameters);
			AliasValue = MyInput -> GetText("alias",AliasNames[i],"value");
		}
	}

	double maxTheta = 0;
	int z;
	for (z=1; z<=M; z++) {
		if (LatBase->WithinBoundaries(z)) {
			maxTheta += LatBase->GetNumLatticeSites(z);
		}
	}
	maxTheta += LatBase->GetNumExtraLatticeSites();
	Array<Text> freedoms(1,7);
	composition = MyInput->GetText("mol", name, "composition", name);
	saveMemory = MyInput->GetBoolean("mol", name, "save_memory", false);
	if (MyInput->ValueSet("mol",name,"force")){
	force = MyInput->GetReal("mol",name,"force",-1000.0,1000.0);
	} else force=0.0;

	Chain = new SF_MolStructure(composition,name,SegQ,LatBase,MyInput);
	composition = Chain->GetComposition();
	N = Chain->GetMaxLength();
	freedoms[1] = "free";
	freedoms[2] = "restricted";
	freedoms[3] = "solvent";
	freedoms[4] = "neutralizer";
	freedoms[5] = "second_generation";
	freedoms[6] = "third_generation";
	freedoms[7] = "range_restricted";
	int choice = MyInput->GetChoice("mol",name,"freedom",freedoms,1);
	Text range;
	// see what freedom the molecule has and if it is consistent with other parameters
	switch (choice) {
		case 1:
			freedom = fixedPhiBulk;
			if (Chain->SomeSegmentsPinned()) {
				Message(fatal,MyInput,"freedom cannot be set to 'free' "
					"for a molecule that contains pinned segments");
			}
			if (Chain->SomeSegmentsGrafted()) {
				Message(fatal,MyInput,"freedom cannot be set to 'free' "
					"for a molecule that contains grafted segments");
			}
			MyInput->DontCombineParam("mol",name,"freedom","theta");
			MyInput->AlwaysCombineParam("mol",name,"freedom","phibulk");
			phiBulk = MyInput->GetReal("mol",name,"phibulk",0,1);
			Chain->SetPhiBulk(phiBulk);
			theta = 0;
			LatRange = LatBase->NewLatticeRange("free");
			break;
		case 2:
			freedom = fixedTheta;
			MyInput->DontCombineParam("mol",name,"freedom","phibulk");
			//MyInput->AlwaysCombineParam("mol",name,"freedom","theta");
			theta = -1;
			if (MyInput->ValueSet("mol",name,"theta")) {
				theta = MyInput->GetReal("mol",name,"theta",0,maxTheta);
				if (MyInput->ValueSet("mol",name,"n"))
				Message(fatal,MyInput,"theta-value and n-value cannot be specified simutaneously for molecule " + name);
			}
			else {
				double Mollength = GetChainLength();
				if (MyInput->ValueSet("mol",name,"n")) theta = MyInput->GetReal("mol",name,"n",0,maxTheta/Mollength)*Mollength;
			}
			if (theta < 0) Message(fatal,MyInput,"theta-value or number of molecules n should be specified for molecule " + name);
			phiBulk = 0;
			#if DEBUG
			printf("SF_Molecule: theta=%f\n\n",theta); fflush(stdout);
			#endif
			LatRange = LatBase->NewLatticeRange("free");
			break;
		case 3:
			freedom = solvent;
			MyInput->DontCombineParam("mol",name,"freedom","phibulk");
			MyInput->DontCombineParam("mol",name,"freedom","theta");
			if (Chain->SomeSegmentsPinned()) {
				Message(fatal,MyInput,"freedom cannot be set to 'solvent' "
					"for a molecule that contains pinned segments");
			}
			if (Chain->SomeSegmentsGrafted()) {
				Message(fatal,MyInput,"freedom cannot be set to 'solvent' "
					"for a molecule that contains grafted segments");
			}
			theta = 0;
			phiBulk = 0;
			LatRange = LatBase->NewLatticeRange("free");
			break;
		case 4:
			freedom = neutralizer;
			MyInput->DontCombineParam("mol",name,"freedom","phibulk");
			MyInput->DontCombineParam("mol",name,"freedom","theta");
			MyInput->DontCombineParam("mol",name,"freedom","restricted_layer");
			MyInput->DontCombineParam("mol",name,"freedom","phi_layer");
			if (Chain->SomeSegmentsPinned()) {
				Message(fatal,MyInput,"freedom cannot be set to 'neutralizer' "
					"for a molecule that contains pinned segments");
			}
			if (Chain->SomeSegmentsGrafted()) {
				Message(fatal,MyInput,"freedom cannot be set to 'neutralizer' "
					"for a molecule that contains grafted segments");
			}
			theta = 0;
			phiBulk = 0;
			LatRange = LatBase->NewLatticeRange("free");
			break;
		case 5:
			freedom = secondGeneration;
			MyInput->DontCombineParam("mol",name,"theta","mu");
			muFixed = MyInput->ValueSet("mol",name,"mu");
			if (muFixed) {
				muFixedValue = MyInput->GetReal("mol",name,"mu",-DBL_MAX,DBL_MAX);
				theta = 0;
			} else {
				theta = MyInput->GetReal("mol",name,"theta",0,maxTheta,Chain->GetAvLength());
			}
			phiBulk = MyInput->GetReal("mol",name,"phibulk",0,1,0);
			Chain->SetPhiBulk(phiBulk);
			MyInput->AlwaysCombineParam("mol",name,"freedom","graft_range");
			range = MyInput->GetText("mol",name,"graft_range");
			LatRange = LatBase->NewLatticeRange(range);
			break;
		case 6:
			freedom = thirdGeneration;
			if (LatBase->GetNumGradients() != 2) {
				Message(fatal,"third generation only available for 2D cylinder lattice");
			}
			phiBulk = MyInput->GetReal("mol",name,"phibulk",0,1);
			Chain->SetPhiBulk(phiBulk);
			theta = MyInput->GetReal("mol",name,"theta",0,maxTheta,Chain->GetAvLength());
			MyInput->AlwaysCombineParam("mol",name,"freedom","theta_renorm");
			thetaRenorm = MyInput->GetReal("mol",name,"theta_renorm",0,M);
			MyInput->AlwaysCombineParam("mol",name,"freedom","graft_range");
			range = MyInput->GetText("mol",name,"graft_range");
			LatRange = LatBase->NewLatticeRange(range);
			MyInput->AlwaysCombineParam("mol",name,"freedom","graft_range_renorm");
 			if (LatBase->GetNumLayers(1) > 1) {
				range = MyInput->GetText("mol",name,"graft_range_renorm");
				LatRangeRenorm = LatBase->NewLatticeRange(range);
			} else {
				LatRangeRenorm = NULL;
			}
			RenormQ.Dim(1,Chain->GetNumDiffSegments());
			break;
		case 7:
			freedom = rangeRestricted;
			MyInput->DontCombineParam("mol",name,"freedom","theta");
			MyInput->DontCombineParam("mol",name,"freedom","phibulk");
			MyInput->AlwaysCombineParam("mol",name,"freedom","restricted_range");
			MyInput->AlwaysCombineParam("mol",name,"freedom","phi_range");
			range = MyInput->GetText("mol",name,"restricted_range");
			restrictedRange = LatBase->NewLatticeRange(range);
			theta = 0;
			phiBulk = 0;
			phiRange = MyInput->GetReal("mol",name,"phi_range",0,1);
			LatRange = LatBase->NewLatticeRange("free");
			break;
		default:
			Message(fatal,MyInput,"Programming error reading freedom in SF_Molecule::SF_Molecule");
			break;
	}
	// further consistency checks
	if (MyInput->ValueSet("mol",name,"phibulk") && Chain->SomeSegmentsPinned()) {
		Message(fatal,MyInput,"phibulk cannot be set for a molecule that contains pinned segments");
	}
	if (MyInput->ValueSet("mol",name,"phibulk") && Chain->SomeSegmentsGrafted()) {
		Message(fatal,MyInput,"phibulk cannot be set for a molecule that contains grafted segments");
	}
	if (MyInput->ValueSet("mol",name,"trains_range")) {
		if(Chain->GetMoleculeType()==branched) Message(fatal,"Attempt to set trains_range for molecule " + name + " Trains are not implemented for a branched molecule");
		range = MyInput->GetText("mol",name,"trains_range");
		LatRangeTrains = LatBase->NewLatticeRange(range);
		if (!MyInput->ValueSet("mol",name,"start_loops_range")) {
			Message(fatal,"when setting 'mol : " + name +
				" : trains_range', start_loops_range must "
				"also be set");
		}
	} else {
		LatRangeTrains = NULL;
	}
	if (MyInput->ValueSet("mol",name,"start_loops_range")) {
		if(Chain->GetMoleculeType()==branched) Message(fatal,"Attempt to set start_loops_range for molecule " + name + " Loops are not implemented for a branched molecule");
		range = MyInput->GetText("mol",name,"start_loops_range");
		LatRangeStartLoops = LatBase->NewLatticeRange(range);
		if (!MyInput->ValueSet("mol",name,"trains_range")) {
			Message(fatal,"when setting 'mol : " + name +
				" : start_loops_range', trains_range must "
				"also be set");
		}
	} else {
		LatRangeStartLoops = NULL;
	}
	// if passed till here, process further parameters
	Array<Text> parameters;
	Array<Text> phiBulkBoundsParam = LatBase->GetNamesBulkBoundaries();
	int numParamTemp;
	if (freedom == thirdGeneration) {
		numParamTemp = 17;
		parameters.Dim(1,numParamTemp+phiBulkBoundsParam.Upperbound());
		parameters[15] = "graft_range";
		parameters[16] = "graft_range_renorm";
		parameters[17] = "theta_renorm";
	} else if (freedom == secondGeneration) {
		numParamTemp = 16;
		parameters.Dim(1,numParamTemp+phiBulkBoundsParam.Upperbound());
		parameters[15] = "graft_range";
		parameters[16] = "mu";
	} else if (freedom == rangeRestricted) {
		numParamTemp = 16;
		parameters.Dim(1,numParamTemp+phiBulkBoundsParam.Upperbound());
		parameters[15] = "restricted_range";
		parameters[16] = "phi_range";
	} else {
		numParamTemp = 14;
		parameters.Dim(1,numParamTemp+phiBulkBoundsParam.Upperbound());
	}
	parameters[1] = "phibulk";
	parameters[2] = "theta";
	parameters[3] = "n";
	parameters[4] = "composition";
	parameters[5] = "freedom";
	parameters[6] = "trains_range";
	parameters[7] = "start_loops_range";
	parameters[8] = "save_memory";
	parameters[9] = "stiff_restricted";
	parameters[10] = "force";
	parameters[11] = "GPU";
	parameters[12] = "find_phi_value";
	parameters[13] = "include_in_TP";
	parameters[14] = "ring";
	int i;
	for (i=numParamTemp+1; i<=numParamTemp+phiBulkBoundsParam.Upperbound(); i++) {
		parameters[i] = "phibulk_" + phiBulkBoundsParam[i-numParamTemp];
	}


	//if (MyInput->ValueSet("mol",name,"force")) {
		//if (!Chain->SomeSegmentsPinned()) {
			//Message(fatal,"Molecule must contain pinned-segments when setting force.");
		//}
	//}
	MyInput->CheckParameterNames("mol",name,parameters);
	numPhi = 0;
	lnGN = 0;
	lnCb = 0;
	lnCt = 0;
	onGPU=false;
	InTP = false;
	if (MyInput->ValueSet("mol",name,"GPU"))
	onGPU = MyInput->GetBoolean("mol", name, "GPU", false);

	if (MyInput->ValueSet("mol",name,"include_in_TP")) {
		InTP = MyInput->GetBoolean("mol", name, "include_in_TP", false);
	}



	if (phiBulkBoundsParam.Upperbound() == 0) {
		phiBulkBoundaries.Dim(0,0);
	} else {
		phiBulkBoundaries.Dim(1,2*LatBase->GetNumGradients());
		// this part depends on the lattice implementation!
		if (MyInput->ValueSet("mol",name,"phibulk")) {
			for (i=1; i<=phiBulkBoundsParam.Upperbound(); i++) {
				if (!MyInput->ValueSet("mol",name,"phibulk_" + phiBulkBoundsParam[i])) {
					Message(fatal,"when setting '" + phiBulkBoundsParam[i] + "' to 'bulk' "
						"'mol : " + name + " : phibulk_" + phiBulkBoundsParam[i]
						+ "' should also be set");
				} else {
					double value = MyInput->GetReal("mol",name,"phibulk_" + phiBulkBoundsParam[i],0,1);
					if (*phiBulkBoundsParam[i] == *Copy("lowerbound")) {
						phiBulkBoundaries[1] = value;
					}
					if (*phiBulkBoundsParam[i] == *Copy("upperbound")) {
						phiBulkBoundaries[2] = value;
					}
					if (*phiBulkBoundsParam[i] == *Copy("lowerbound_x")) {
						phiBulkBoundaries[1] = value;
					}
					if (*phiBulkBoundsParam[i] == *Copy("upperbound_x")) {
						phiBulkBoundaries[2] = value;
					}
					if (*phiBulkBoundsParam[i] == *Copy("lowerbound_y")) {
						phiBulkBoundaries[3] = value;
					}
					if (*phiBulkBoundsParam[i] == *Copy("upperbound_y")) {
						phiBulkBoundaries[4] = value;
					}
				}
			}
		}
		Chain->SetPhiBulkBoundaries(phiBulkBoundaries);
	}
}
SF_Molecule::~SF_Molecule() {
	delete LatRange;
	delete Chain;
	if (MyInput->ValueSet("mol",name,"start_loops_range")) {
		delete LatRangeStartLoops;
	}
	if (MyInput->ValueSet("mol",name,"trains_range")) {
		delete LatRangeTrains;
	}
	if (freedom == thirdGeneration) {
		delete LatRangeRenorm;
	}
	if (freedom == rangeRestricted) {
		delete restrictedRange;
	}
}
Text
SF_Molecule::GetName() const {
	return name;
}

Boolean
SF_Molecule::GetIncludeTP() const {
	return InTP;
}

void
SF_Molecule::GetOutput(Output* Out) const {
	Text t;
	if (freedom == fixedPhiBulk) t = "free";
	if (freedom == fixedTheta) t = "restricted";
	if (freedom == solvent) t = "solvent";
	if (freedom == neutralizer) t = "counterion";
	if (freedom == secondGeneration) t = "secondGeneration";
	if (freedom == thirdGeneration) t = "thirdGeneration";
	if (freedom == rangeRestricted) t = "range_restricted";
	Out->PutText("mol",name,"freedom",t);
	Out->PutReal("mol",name,"phibulk",phiBulk);
	if (phiBulk>0) 	Out->PutReal("mol",name,"log(phibulk)",log10(phiBulk)); else Out->PutReal("mol",name,"log(phibulk)",0);
	double thetaComputed = ComputeTheta();
	Out->PutReal("mol",name,"phi-av",thetaComputed/LatBase->GetTotalNumLatticeSites());
	Array<Text> phiBulkBoundsParam;
	int i;
	for (i=1; i<=phiBulkBoundsParam.Upperbound(); i++) {
		if (MyInput->ValueSet("mol",name,"phibulk_" + phiBulkBoundsParam[i])) {
			if (*phiBulkBoundsParam[i] == *Copy("lowerbound")) {
				Out->PutReal("mol",name,"phibulk_" + phiBulkBoundsParam[i],phiBulkBoundaries[1]);
			}
			if (*phiBulkBoundsParam[i] == *Copy("upperbound")) {
				Out->PutReal("mol",name,"phibulk_" + phiBulkBoundsParam[i],phiBulkBoundaries[2]);
			}
			if (*phiBulkBoundsParam[i] == *Copy("lowerbound_x")) {
				Out->PutReal("mol",name,"phibulk_" + phiBulkBoundsParam[i],phiBulkBoundaries[1]);
			}
			if (*phiBulkBoundsParam[i] == *Copy("upperbound_x")) {
				Out->PutReal("mol",name,"phibulk_" + phiBulkBoundsParam[i],phiBulkBoundaries[2]);
			}
			if (*phiBulkBoundsParam[i] == *Copy("lowerbound_y")) {
				Out->PutReal("mol",name,"phibulk_" + phiBulkBoundsParam[i],phiBulkBoundaries[3]);
			}
			if (*phiBulkBoundsParam[i] == *Copy("upperbound_y")) {
				Out->PutReal("mol",name,"phibulk_" + phiBulkBoundsParam[i],phiBulkBoundaries[4]);
			}
		}
	}
	double ThetaExcess = ComputeThetaExcess();
	Out->PutReal("mol",name,"theta",thetaComputed);
	Out->PutReal("mol",name,"theta excess",ThetaExcess);
	Out->PutReal("mol",name,"n",thetaComputed/Chain->GetAvLength());
	Out->PutReal("mol",name,"n exc",ThetaExcess/Chain->GetAvLength());
	Out->PutText("mol",name,"composition",composition);
	if (MyInput->GetBoolean("mol",name,"ring",false)) {
		Out->PutBoolean("mol",name,"ring",MyInput->GetBoolean("mol",name,"ring",false));
	}
	if (int(Chain->GetAvLength()) -Chain->GetAvLength() == 0) {
		Out->PutInt("mol",name,"chainlength",int(Chain->GetAvLength()));
	} else {
		Out->PutReal("mol",name,"av. chainlength",Chain->GetAvLength());
	}
	Out->PutReal("mol",name,"ln(G(N|1))",lnGN);
	Out->PutReal("mol",name,"force",GetForce());
	if (freedom != fixedTheta && freedom != rangeRestricted) {
		Out->PutReal("mol",name,"ln(norm_free)",lnCb);
	}
	if (freedom != fixedPhiBulk && freedom != solvent && freedom != neutralizer) {
		Out->PutReal("mol",name,"ln(norm_restricted)",lnCt);
	}
	if (freedom == secondGeneration) {
		Out->PutReal("mol",name,"thetaRestricted",theta);
		Out->PutText("mol",name,"graft_range",LatRange->GetOutput());
	} else if (freedom == rangeRestricted) {
		Out->PutReal("mol",name,"phi_range",phiRange);
		Out->PutText("mol",name,"restricted_range",restrictedRange->GetOutput());
	} else {
		if (MyInput->ValueSet("mol",name,"trains_range")) {
			Out->PutText("mol",name,"trains_range",LatRangeTrains->GetOutput());
		}
		if (MyInput->ValueSet("mol",name,"start_loops_range")) {
			Out->PutText("mol",name,"start_loops_range",LatRangeStartLoops->GetOutput());
		}
	}
	if (Charged()) {
		Out->PutReal("mol",name,"charge", GetCharge());
		Out->PutReal("mol",name,"number of charges",GetNumberOfCharges());
	}
	if (LatBase->GetNumGradients() == 1) {
		Vector phiExc = GetPhi(total);
		int z;
		for (z=1; z<=M; z++) {
			phiExc[z] -= phiBulk;
		}
		LatBase->MultiplyWithLatticeSites(phiExc);
		LatticeRange* LayerZero = new LatticeRange1D(1,1,M);
		double R = LatBase->MomentUnweighted(phiExc,1,LayerZero,-0.5)/LatBase->MomentUnweighted(phiExc,0,LayerZero,-0.5);
		double R2 = LatBase->MomentUnweighted(phiExc,2,LayerZero,-0.5)/LatBase->MomentUnweighted(phiExc,0,LayerZero,-0.5);
		Out->PutReal("mol",name,"1st moment exc. n",R);
		Out->PutReal("mol",name,"2nd moment exc. n",pow(R2,0.5));

		LatBase->DivideByLatticeSites(phiExc);
		R = LatBase->MomentUnweighted(phiExc,1,LayerZero,-0.5)/LatBase->MomentUnweighted(phiExc,0,LayerZero,-0.5);
		R2 = LatBase->MomentUnweighted(phiExc,2,LayerZero,-0.5)/LatBase->MomentUnweighted(phiExc,0,LayerZero,-0.5);
		delete LayerZero;
		for (z=1; z<=M; z++) {
			phiExc[z] += phiBulk;
		}

		Out->PutReal("mol",name,"1st moment exc. phi",R);
		Out->PutReal("mol",name,"2nd moment exc. phi",pow(R2,0.5));


	}
	Vector PHI;
	Vector PHITOT(1,M);
	//SF_MolSegment* SEG;

	for (i=1; i<=numPhi; i++) {
		Vector phix;
		switch (DensityPartQ[i]) {
			case total:
				//numSegents = Chain->GetNumDiffSegments();
				//for (int j=1; j<=numSegents; j++) {
				//	SEG = Chain->GetDiffSegment(j);
				//	PHI = SEG->GetPhi(total);
				//	for (int z=1; z<=M; z++) PHITOT[z]+=PHI[z];
				//}
				//Out->PutProfile("mol",name,"phi",PHITOT);
				Out->PutProfile("mol",name,"phi",PhiQ[i]);
				break;
			case unconstrained:
				Out->PutProfile("mol",name,"phi unconstrained",PhiQ[i]);
				break;
			case constrained:
				Out->PutProfile("mol",name,"phi constrained",PhiQ[i]);
				break;
			case renorm:
				Out->PutProfile("mol",name,"phi renorm",PhiQ[i]);
				break;
			//case x_dir:
			//	Out->PutProfile("mol",name,"phi_x",PhiQ[i]);
			//	break;
			//case y_dir:
			//	Out->PutProfile("mol",name,"phi_y",PhiQ[i]);
			//	break;
			//case z_dir:
			//	Out->PutProfile("mol",name,"phi_z",PhiQ[i]);
			//	break;
			//case yz_dir:
			//	Out->PutProfile("mol",name,"phi_yz",PhiQ[i]);
			//	break;
			default:
				Message(fatal,"programming error: unknown densitypart requested in SF_Molecule::GetOutput");
				break;
		}

	}
	//alias output;
	if (MyInput->GetNumNames("alias") > 0) {
		Text AliasValue;
		int numAlias = MyInput->GetNumNames("alias");
		for (int i=1; i<=numAlias; i++) {
			AliasValue = MyInput -> GetText("alias",AliasNames[i],"value");
			if (MyInput->CheckInt(AliasValue)) {
				int aliasvalue = MyInput->GetInt("alias",AliasNames[i],"value",0,10000000);
				Out->PutInt("alias",AliasNames[i],"value",aliasvalue);
			} else {
				Out->PutText("alias",AliasNames[i],"value",AliasValue);
			}
		}
	}



	int numMolSeg = Chain->GetNumDiffSegments();
	if (numMolSeg == 1 && (Chain->GetDiffSegment(1))->GetNumStates() == 1) {
		numMolSeg = 0;
	}
	for (i=1; i<=numMolSeg; i++) {
		(Chain->GetDiffSegment(i))->GetOutput(Out);
	}
}
Vector
SF_Molecule::GetPhi(const DensityPart densityPart) const {
	int i;
	// FIXME this produces wrong results when used for branched homopolymer ??
	for (i=1; i<=numPhi; i++) {
		if (DensityPartQ[i] == densityPart) {
			LatBase->SubtractBoundaries(PhiQ[i]);
			LatBase->RestoreBoundaries(PhiQ[i]);
			return PhiQ[i];
		}
	}
	if (densityPart == total) {
		Vector phiTotal(1,M);
		Vector phi;

		for (i=1; i<=numPhi; i++) {
			phi = PhiQ[i];
			for (int z=1; z<=M; z++) {
				phiTotal[z] += phi[z];
			}
		}
		return phiTotal;
	}
	Message(fatal,"Programming error in call to SF_Molecule::"
		"GetPhi(DensityPart densityPart), Create phi before using it");
	return PhiQ[1]; // never get here
}
// allocate the Phi[z] vectors, descends to SF_MolSegemnt->SF_MolState
void
SF_Molecule::CreatePhi(const DensityPart densityPart) {
	int i;
	SF_MolSegment* MolSeg;

	int numSeg = Chain->GetNumDiffSegments();
	for (i=1; i<=numSeg; i++) {
		MolSeg = Chain->GetDiffSegment(i);
		MolSeg->CreatePhi(densityPart);

		//if (LatBase->GetMatrixApproach()==secondOrder && densityPart==total) {
		//	if (gradients==1) {MolSeg->CreatePhi(x_dir); MolSeg->CreatePhi(yz_dir);} else {
		//		MolSeg->CreatePhi(x_dir);MolSeg->CreatePhi(y_dir);MolSeg->CreatePhi(z_dir);
		//	}
		//}

	}

	if (numPhi == 0) {
	// if there is no Phi so far, create Q's of Dim(1,1):
		PhiQ.Dim(1,1);
		DensityPartQ.Dim(1,1);
		Vector phi;
		if (numSeg > 1) phi.Dim(1,LatBase->GetTotalNumLayers());
		else {
			MolSeg = Chain->GetDiffSegment(1);
			phi = MolSeg->GetPhi(densityPart);
		}
		PhiQ[1] = phi;
		DensityPartQ[1] = densityPart;
		numPhi = 1;
	} else {
	// else create Q's of Dim(1,numPhi++)
		Array<Vector> PhiQtemp(1,numPhi+1);
		Array<DensityPart> DensityPartQtemp(1,numPhi+1);
		for (i=1; i<=numPhi; i++) {
			PhiQtemp[i] = PhiQ[i];
			DensityPartQtemp[i] = DensityPartQ[i];
		}
		Vector phi;
		if (numSeg > 1) phi.Dim(1,LatBase->GetTotalNumLayers());
		else {
			MolSeg = Chain->GetDiffSegment(1);
			phi = MolSeg->GetPhi(densityPart);
		}
		numPhi++;
		PhiQtemp[numPhi] = phi;
		DensityPartQtemp[i] = densityPart;
		PhiQ = PhiQtemp;
		DensityPartQ = DensityPartQtemp;
	}
}
void
SF_Molecule::DeletePhi(const DensityPart densityPart) {
	int i;
	int number = 0;
	SF_MolSegment* MolSeg;
	for (i=1; i<=numPhi; i++) {
		if (DensityPartQ[i] == densityPart) number = i;
	}
	if (number == 0) {
		Message(fatal,"Programming error in call to SF_Molecule::"
		"DeletePhi(DensityPart densityPart), Create phi before deleting it");
	}

	int numSeg = Chain->GetNumDiffSegments();
	for (i=1; i<=numSeg; i++) {
		MolSeg = Chain->GetDiffSegment(i);
		//if (LatBase->GetMatrixApproach()==secondOrder && densityPart==total){
		//	if (gradients==1) {MolSeg->DeletePhi(x_dir); MolSeg->DeletePhi(yz_dir);} else {
		//		MolSeg->DeletePhi(x_dir);MolSeg->DeletePhi(y_dir);MolSeg->DeletePhi(z_dir);
		//	}
		//}
	}

	Array<Vector> PhiQtemp(1,numPhi-1);
	Array<DensityPart> DensityPartQtemp(1,numPhi-1);
	for (i=1; i<=numPhi; i++) {
		if (i < number) {
			PhiQtemp[i] = PhiQ[i];
			DensityPartQtemp[i] = DensityPartQ[i];
		} else if (i > number) {
			PhiQtemp[i-1] = PhiQ[i];
			DensityPartQtemp[i-1] = DensityPartQ[i];
		}
	}
	numPhi--;
	PhiQ = PhiQtemp;
	DensityPartQ = DensityPartQtemp;

	for (i=1; i<=numSeg; i++) {
		MolSeg = Chain->GetDiffSegment(i);
		MolSeg->DeletePhi(densityPart);
	}
}
Array<DensityPart>
SF_Molecule::GetDensityPartQ() const {
	return DensityPartQ;
}
SF_MolStructure*
SF_Molecule::GetMolStructure() const {
	return Chain;
}
MolFreedom
SF_Molecule::GetFreedom() const {
	return freedom;
}
void
SF_Molecule::SetFreedom(MolFreedom freedom_) {
	freedom = freedom_;
}
double
SF_Molecule::GetChainLength() const {
	return Chain->GetAvLength();
}

Boolean
SF_Molecule::GetGPU() const{
	return onGPU;
}

double
SF_Molecule::GetForce() const {
	return force;
}

//Boolean
//SF_Molecule::MayerSaupeSet() const{
	//cout << "MSset in SF_Molecule" << endl;
//	Boolean MS_found=false;
//	int numMolSeg = Chain->GetNumDiffSegments();
//	for (int i=1; i<=numMolSeg; i++) {
//	if (Chain->GetDiffSegment(i)->MayerSaupeSet()) MS_found=true;
//	}
//	return MS_found;
//}

void
SF_Molecule::SetPhiBulk(double phiBulk_) {
	phiBulk = phiBulk_;
	Chain->SetPhiBulk(phiBulk);
	SegQ->UpdatePhiBulk();
}
double
SF_Molecule::GetPhiBulk(void) const {
	return phiBulk;
}
Vector
SF_Molecule::GetPhiBulkBoundaries(void) const {
	return phiBulkBoundaries;
}
void
SF_Molecule::SetTheta(double value) {
	theta = value;
}
double
SF_Molecule::GetTheta(void) const {
	return theta;
}
double
SF_Molecule::ComputeTheta(void) const {
	Vector phiTot(1,M);
	Vector phi;
	int z;
	for (int i=1; i<=numPhi; i++) {
		phi = PhiQ[i];
		for (z=1; z<=M; z++) {
			phiTot[z] += phi[z];
		}
	}
	LatBase->SubtractBoundaries(phiTot);

	LatBase->MultiplyWithLatticeSites(phiTot);
	double calcTheta = 0;
	for (z=1; z<=M; z++) {
		calcTheta += phiTot[z];
	}
	return calcTheta;
}
double
SF_Molecule::ComputeThetaExcess(void) const {
	Vector phiTot(1,M);
	Vector phi;
	int z;
	for (int i=1; i<=numPhi; i++) {
		phi = PhiQ[i];
		for (z=1; z<=M; z++) {
			phiTot[z] += phi[z];
		}
	}
	for (z=1; z<=M; z++) {
		if (phiTot[z] > 0) {
			phiTot[z] -= phiBulk;
		}
	}
	LatBase->SubtractBoundaries(phiTot);
	LatBase->MultiplyWithLatticeSites(phiTot);
	double calcTheta = 0;
	for (z=1; z<=M; z++) {
		calcTheta += phiTot[z];
	}
	return calcTheta;
}
double
SF_Molecule::GetLnGN() const {
	return lnGN;
}
double
SF_Molecule::GetLnNormFree() const {
	return lnCb;
}
double
SF_Molecule::GetLnNormRestricted() const {
	return lnCt;
}
void
SF_Molecule::SetLnNormRestricted(double value) {
	lnCt = value;
}
double
SF_Molecule::GetThetaRenorm() const {
	if (freedom != thirdGeneration) {
		Message(fatal,"programming error in call to SF_Molecule::GetThetaRenorm");
	}
	return thetaRenorm;
}
Vector
SF_Molecule::GetRenorm(int i) const {
	return RenormQ[i];
}
void
SF_Molecule::SetRenorm(Vector renorm, int i) {
	if (freedom != thirdGeneration) {
		Message(fatal,"programming error in call to SF_Molecule::SetPhiRestricted");
	}
	RenormQ[i] = renorm;
}
Boolean
SF_Molecule::Charged() const {
	int numDiffSeg = Chain->GetNumDiffSegments();
	SF_MolSegment* MolSeg;
	SF_MolState* MolState;

	for (int i=1; i<=numDiffSeg; i++) {
		MolSeg = Chain->GetDiffSegment(i);
		int numStates = MolSeg->GetNumStates();
		for (int j=1; j<= numStates; j++) {
			MolState = MolSeg->GetState(j);
			if (MolState->GetValence() != 0) {
				return true;
			}
		}
	}
	return false;
}
double
SF_Molecule::GetBulkCharge() const {
	int numDiffSeg = Chain->GetNumDiffSegments();
	SF_MolSegment* MolSeg;
	SF_MolState* MolState;
	double charge = 0;
	for (int i=1; i<=numDiffSeg; i++) {
		MolSeg = Chain->GetDiffSegment(i);
		double number = Chain->GetAvNumSegments(MolSeg);
		int numStates = MolSeg->GetNumStates();
		for (int j=1; j<= numStates; j++) {
			MolState = MolSeg->GetState(j);
			charge += MolState->GetAlphaBulk()*MolState->GetValence()*number;
		}
	}
	return charge;
}
double
SF_Molecule::GetCharge() const {
	int z,i,j;
	int numDiffSeg, numStates;
	double q, sum, theta;
	SF_MolSegment* MolSeg;
	SF_MolState* MolState;
	Vector phi;

	numDiffSeg = Chain->GetNumDiffSegments();
	q = theta = 0;
	for (i=1; i<=numDiffSeg; i++) {
		MolSeg = Chain->GetDiffSegment(i);
		numStates = MolSeg->GetNumStates();
		for (j=1; j<= numStates; j++) {
			MolState = MolSeg->GetState(j);
			phi = MolState->GetPhi(total);
			sum = 0;
			for (z=1; z<=M; z++) {
				sum += phi[z]
					*LatBase->GetNumLatticeSites(z);
			}
			theta += sum;
			q += sum * MolState->GetValence();
		}
	}
	q *= Chain->GetAvLength()/theta;
	return q;
}
/* this function returns the center of mass of all charges in this molecule
   The variable 'sign' gives the sign of the charges to be counted. If sign
   is zero, then all charges are counted.
*/
double
SF_Molecule::GetChargeCoM(int sign) const {
	int i, j, z;
	int numDiffSeg, numStates;
	double CoM;
	double valence;
	double sum;
	double c;
	Vector charge, phi;
	SF_MolSegment* molSeg;
	SF_MolState* molState;

	sum = 0;
	charge.Dim(1,M);
	numDiffSeg = Chain->GetNumDiffSegments();
	for (i=1; i<=numDiffSeg; i++) {
		molSeg = Chain->GetDiffSegment(i);
		numStates = molSeg->GetNumStates();
		for (j=1; j<= numStates; j++) {
			molState = molSeg->GetState(j);
			valence = molState->GetValence();
			if (sign == 0 || valence*sign > 0) {
				phi = molState->GetPhi(total);
				valence = fabs(valence);
				for (z=1; z<=M; z++) {
					c = phi[z]*valence;
					charge[z] += c;
					sum += c*LatBase->GetNumLatticeSites(z);
				}
			}
		}
	}
	CoM = LatBase->Moment(charge,1,NULL,2); //y coordinate of CoM
	CoM /= sum;
	return CoM;
}
double
SF_Molecule::GetNumberOfCharges() const {
	int z,i,j;
	int numDiffSeg, numStates;
	double q, sum, theta;
	SF_MolSegment* MolSeg;
	SF_MolState* MolState;
	Vector phi;

	numDiffSeg = Chain->GetNumDiffSegments();
	q = theta = 0;
	for (i=1; i<=numDiffSeg; i++) {
		MolSeg = Chain->GetDiffSegment(i);
		numStates = MolSeg->GetNumStates();
		for (j=1; j<= numStates; j++) {
			MolState = MolSeg->GetState(j);
			phi = MolState->GetPhi(total);
			sum = 0;
			for (z=1; z<=M; z++) {
				sum += phi[z]
					*LatBase->GetNumLatticeSites(z);
			}
			theta += sum;
			q += sum * fabs(MolState->GetValence());
		}
	}
	q *= Chain->GetAvLength()/theta;
	return q;
}
Boolean
SF_Molecule::MuFixed(void) const {
	return muFixed;
}
double
SF_Molecule::GetMuFixed(void) const {
	return muFixedValue;
}

