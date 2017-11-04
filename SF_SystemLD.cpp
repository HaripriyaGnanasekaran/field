#include "SF_SystemLD.h"

SF_SystemLD::SF_SystemLD(Input* MyInput_, Text name_, Boolean compute_) {
	MyInput = MyInput_;
	name = name_;
	compute = compute_;

	Array<Text> param(1,9);
	param[1] = "calculation_type";
	param[2] = "matrix_approach";
	param[3] = "overflow_protection";
	param[4] = "temperature";
	param[5] = "warnings";
	param[6] = "LDT";
	param[7] = "mobility";
	param[8] = "free_energy_type";
	param[9] = "random_seed";

	//MayerSaupe=false;
	MyInput->CheckParameterNames("sys",name,param);
	MyInput->SetAllMessagesOff();
	if (MyInput->ValueSet("sys",name,"warnings")) {
		if (MyInput->GetBoolean("sys",name,"warnings")) {
			MyInput->SetAllMessagesOn();
		}
	}

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

	overflowProtection = MyInput->GetBoolean("sys",name,"overflow_protection",false);
	temperature = MyInput->GetReal("sys",name,"temperature",0.0001,DBL_MAX,ROOMTEMPERATURE);
	TEMPERATURE = temperature;

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
	TEMPERATURE = temperature;
	bool LD_found=false;
	for (int i=1; i<=MolQ->GetNumMolecules(); i++) {
		SF_Molecule* Mol = MolQ->GetMolecule(i);
		if (Mol->GetMolStructure()->SomeSegmentsLD()) {
			LD_found=true;
		}
	}

	if (!LD_found) {
		Message(fatal,"No molecule with pinned segments and activated LD found");
	}
	if (!MyInput->ValueSet("sys", name, "LDT") ) Message(fatal,"For LD_SCF you need to set LDT value (typically 10 000)");
	else LDT = MyInput->GetInt("sys",name,"LDT",1,1000000);
	if (!MyInput->ValueSet("sys", name, "mobility") ) {
		Message(literal,"For LD_SCF you should set 'mobility' value (default 1 is used)"); mobility=1.0;
	} else 	mobility = MyInput->GetReal("sys",name,"mobility",0.000000001,100,1.0);

	if (MyInput->ValueSet("sys",name,"free_energy_type")) {
	 	Array<Text> FChoices(1,2);
		FChoices[1] = "free_energy";
		FChoices[2] = "free_energy(po)";

		LD_F = MyInput->GetChoice("sys",name,"free_energy_type",FChoices)==1;
	} else Message(fatal,"free_energy_type should be specified when LD_SCF option is selected; set this quantity to 'free_energy', or 'free_energy(po)'");

}

double
SF_SystemLD::random_d() {
	return rand()/(RAND_MAX+1.0);
}

int
SF_SystemLD::random_int(int low, int high) {
	int range = high-low+1;
	return low+int(range*random_d());
}
double
SF_SystemLD::random(double low, double high) {
	double range = high-low;
	return low+range*random_d();
}

double
SF_SystemLD::normal(double low, double high) {
	double U1=random(low,high);
	double U2=random(low,high);
	return sqrt(-log(U1))*cos(2*PI*U2);
}

bool // checks wether two particles have the same position and returns false if this is the case;
SF_SystemLD::CheckOverlap(){
	for (int i=1; i <= numPos-1; i++) {
		for (int j=i+1; j <= numPos; j++) {
			if (r_new[i]== r_new[j]) {return false; }
		}
	}
	return true;
}


bool
SF_SystemLD::MovePos(SF_Segment** Seg){
	double step;
	double distance;
	for (int p=1; p<=numPos; p++) {
		distance = pow(rx_new[p]-rx_old[p],2)+pow(ry_new[p]-ry_old[p],2)+pow(rz_new[p]-rz_old[p],2);
		if (distance==0) distance=1;
		step=(F_new[p]-F_old[p])*(rx_new[p]-rx_old[p])/distance*mobility + sqrt(2*mobility)*normal(0,1);
		rx_new[p]=rx_old[p]+step;
		step=(F_new[p]-F_old[p])*(ry_new[p]-ry_old[p])/distance*mobility + sqrt(2*mobility)*normal(0,1);
		ry_new[p]=ry_old[p]+step;
		step=(F_new[p]-F_old[p])*(rz_new[p]-rz_old[p])/distance*mobility + sqrt(2*mobility)*normal(0,1);
		rz_new[p]=rz_old[p]+step;

		if (rx_new[p]<1) rx_new[p]=rx_new[p]+n_layers_x;
		if (ry_new[p]<1) ry_new[p]=ry_new[p]+n_layers_y;
		if (rz_new[p]<1) rz_new[p]=rz_new[p]+n_layers_z;
		if (rx_new[p]>n_layers_x+1) rx_new[p]=rx_new[p]-n_layers_x;
		if (ry_new[p]>n_layers_y+1) ry_new[p]=ry_new[p]-n_layers_y;
		if (rz_new[p]>n_layers_z+1) rz_new[p]=rz_new[p]-n_layers_z;
		x_new[p]=int(rx_new[p]); y_new[p]=int(ry_new[p]); z_new[p]=int(rz_new[p]);
		r_new[p]=jx*x_new[p]+jy*y_new[p]+z_new[p]+1;
		if (!free[r_new[p]]) {return false;}
	}
	return CheckOverlap();
}


void
SF_SystemLD::Go(Output* Out,Lattice* Lat) {
	numPos=0;
	int numSeg = SegQ->GetNumSegments();
	int numStates = SegQ->GetNumStates();
	int dim = Lat->GetNumGradients();

	if (dim != 3) {Message(fatal,"Sorry: for LD we need 3 gradients. ");}
	n_layers_x = Lat->GetNumLayers(1)-2;
	n_layers_y = Lat->GetNumLayers(2)-2;
	n_layers_z = Lat->GetNumLayers(3)-2;
	jx = (n_layers_y+2)*(n_layers_z+2);
	jy = n_layers_z+2;
	jz = 1;
	int M = (n_layers_x+2)*(n_layers_y+2)*(n_layers_z+2);
	free = new bool[M+1];
	rho_0.Dim(1,MolQ->GetNumMolecules());
	for (int i=1; i<=MolQ->GetNumMolecules(); i++) {
		rho_0[i].Dim(1,n_layers_z);
	}
	Omega_0.Dim(1,n_layers_z);
	F_0.Dim(1,n_layers_z);
	Fpo_0.Dim(1,n_layers_z);
	if (SegQ->Charged()) PSI_0.Dim(1,n_layers_z);

	LatticeRange* LatRange;
	SF_State* State;

	for (int z=0; z<=M; z++) free[z]=true;
	for (int i=1; i<=numStates; i++) {
		State = SegQ->GetState(i);
		if (State->GetFreedom() == frozen) {
		LatRange = State->GetLatRange();
			for (int z=0; z<=M; z++) {
				if (LatRange->InRange(z)) {
					free[z] = false;
				}
			}
		}
	}

	for (int i=1; i<=numSeg; i++) {
		SF_Segment* Seg = SegQ->GetSegment(i);
		if (Seg->GetLD()) {
			LatRange=Seg->GetLatRange();
			numPos += LatRange->GetNumPos();
		}
	}
	if (Solve->e_info || Solve->s_info) {cout << "LD started" << endl;}
	cout << "Num of positions " << numPos  << endl;
	//double Fpos0[numPos+1],Fpos[numPos+1];
	LatticeRange* Range[numPos+1];
	SF_Segment*  Seg[numPos+1];
	//double FposTot=0;
	x_old = new int[numPos+1];
	y_old = new int[numPos+1];
	z_old = new int[numPos+1];
	r_old = new int[numPos+1];
	x_new = new int[numPos+1];
	y_new = new int[numPos+1];
	z_new = new int[numPos+1];
	r_new = new int[numPos+1];
	r_last = new int[numPos+1];
	submask = new double[27];
	rx_old = new double[numPos+1];
	ry_old = new double[numPos+1];
	rz_old = new double[numPos+1];
	rx_new = new double[numPos+1];
	ry_new = new double[numPos+1];
	rz_new = new double[numPos+1];
	F_old = new double[numPos+1];
	F_new = new double[numPos+1];

	int k=0;
	int rr,rx,ry,rz,np;
	for (int i=1; i<=numSeg; i++) {
		SF_Segment* Segm = SegQ->GetSegment(i);
		if (Segm->GetLD()) {
			LatRange=Segm->GetLatRange();
			np = LatRange->GetNumPos();
			for (int j=1; j<=np; j++) {
				k++;
				Seg[k]=Segm;
				Range[k]=LatRange;
				r_new[k] = Range[k]->GetPos(j);
				rr=r_new[k]-1;
				if (!free[r_new[k]]) {
					cout << "pos of particle " << k << "= " <<rr << endl;
					Message(fatal,"You try to insert a particle in frozen range. ");
				}
				rz = ((rr % jx) % jy)/jz;
				ry = ((rr % jx)-rz*jz)/jy;
				rx = (rr-ry*jy-rz*jz)/jx;

				if (rx>0 && rx <n_layers_x+1) x_new[k]=rx; else {
					cout << "pos of particle " << k << "= " <<rr << endl;
					Message(fatal,"Initial x value out of bounds. "); }
				if (ry>0 && ry <n_layers_y+1) y_new[k]=ry; else {
					cout << "pos of particle " << k << "= " <<rr << endl;
					Message(fatal,"Initial y value out of bounds. "); }
				if (rz>0 && rz <n_layers_z+1) z_new[k]=rz; else {
					cout << "pos of particle " << k << "= " <<rr << endl;
					Message(fatal,"Initial z value out of bounds. "); }
			}

		}
	}


	for (int i=1; i<=numSeg; i++) {
		SF_Segment* Segm = SegQ->GetSegment(i);
		if (Segm->GetLD()) {
			Segm ->ClearAllPos();
		}
	}

	if (!CheckOverlap()) Message(fatal,"Initial positions overlap. Possibly you have two types particles moved by LD; Try new initial positions");

	for (int p=1; p<=numPos; p++) {
		rx_old[p]=rx_new[p]=1.0*x_new[p];
		ry_old[p]=ry_new[p]=1.0*y_new[p];
		rz_old[p]=rz_new[p]=1.0*z_new[p];
		int l=0;
		for (int i=-1; i<2; i++)
		for (int j=-1; j<2; j++)
		for (int k=-1; k<2; k++) {if (free[jx*(x_new[p]-i) +jy*(y_new[p]-j)+jz*(z_new[p]-k)+1]) submask[l]=1.0; else submask[l]=0.0; l++; }
		Seg[p]->UpdatePos(rx_new[p],ry_new[p],rz_new[p],submask);
	}

	Solve->SetInitialGuess(Solve);
	Solve->Iterate();

	Vector profile;
	if (LD_F) {
		profile = MolQ->GetFreeEnergyProfile();
		F_system_N=MolQ->GetFreeEnergy();
	} else {
		profile = MolQ->GetGrandPotentialProfile();
		AddChemPotGraft(profile,Lat);
		F_system_N=GetFreeEnergyPo(Lat);
	}
	for (int p=1; p<=numPos; p++) {
		r_last[p]=r_new[p];
		int rr;
		double distance=0;
		int l=0; double norm=0;
		for (int i=-1; i<2; i++)
		for (int j=-1; j<2; j++)
		for (int k=-1; k<2; k++) {
			rr=jx*(x_new[p]-i) +jy*(y_new[p]-j)+jz*(z_new[p]-k)+1;
			distance = pow((rx_new[p]-x_new[p]),2)+pow((ry_new[p]-y_new[p]),2)+pow((rz_new[p]-z_new[p]),2);
			if (free[jx*(x_new[p]-i) +jy*(y_new[p]-j)+jz*(z_new[p]-k)+1]) submask[l]=exp(-distance); else submask[l]=0.0;
			norm+=submask[l]; l++;
		}
		l=0; F_new[p]=0;
		for (int i=-1; i<2; i++)
		for (int j=-1; j<2; j++)
		for (int k=-1; k<2; k++) {
			rr=jx*(x_new[p]-i) +jy*(y_new[p]-j)+jz*(z_new[p]-k)+1;
			F_new[p] += profile[rr]*submask[l]/norm; l++;
		}
		F_old[p]=F_new[p];
	}
	ldt=0;
	cout << "LDT "<< ldt << " F = " << F_system_N << endl;
	for (int k=1; k<=numPos; k++) cout << "("<< x_new[k] << ","<< y_new[k] <<","<< z_new[k] <<")";
	cout << endl;
	SetRho(Lat,rho_0); SetF(Lat,F_0); SetFpo(Lat,Fpo_0); SetOmega(Lat,Omega_0);
	if (SegQ->Charged()) SetPsi(Lat,PSI_0);
	if (Solve->ErrorOccurred()) {F_system_N = 10000000.0;}
	bool solved = false;
	while ( ldt < LDT) {
		bool write_output = true;
		ldt++;
		F_system_O = F_system_N;
		//for (int k=1; k<=numPos; k++) {
		//	x_old[k]=x_new[k];
		//	y_old[k]=y_new[k];
		//	z_old[k]=z_new[k];
		//	r_old[k]=r_new[k];
		//	F_old[k]=F_new[k];
		//}

		//for (int k=1; k<=numPos; k++) {
		//	rx_old[k]=rx_new[k];
		//	ry_old[k]=ry_new[k];
		//	rz_old[k]=rz_new[k];
		//}

		while (!MovePos(Seg)){// If a particle has been placed on top of another particle or in a frozen range the particles are set back to their original positions
			for (int k=1; k<=numPos; k++) {
				x_new[k]=x_old[k];
				y_new[k]=y_old[k];
				z_new[k]=z_old[k];
				r_new[k]=r_old[k];
			}

			for (int k=1; k<=numPos; k++) {
				rx_new[k]=rx_old[k];
				ry_new[k]=ry_old[k];
				rz_new[k]=rz_old[k];
			}
			cout << "Overlap has occured. Try to move particles again." << endl;
		}


		Solve->MoveInitialGuess(r_new,r_last,numPos);
		Solve->SetInitialGuess(GetSolve());
		for (int i=1; i<=numSeg; i++) {
			SF_Segment* Segm = SegQ->GetSegment(i);
			if (Segm->GetLD()) {
				Segm ->ClearAllPos();
			}
		}

		for (int p=1; p<=numPos; p++) {
			int l=0;
			for (int i=-1; i<2; i++)
			for (int j=-1; j<2; j++)
			for (int k=-1; k<2; k++) {if (free[jx*(x_new[p]-i) +jy*(y_new[p]-j)+jz*(z_new[p]-k)+1]) submask[l]=1; else submask[l]=0; l++; }
			Seg[p]->UpdatePos(rx_new[p],ry_new[p],rz_new[p],submask);
		}

		Solve->SetInitialGuess(Solve);
		Solve->Iterate();
		solved = !Solve->ErrorOccurred();
		if (solved) {
			if (LD_F) {
			profile = MolQ->GetFreeEnergyProfile();
				F_system_N=MolQ->GetFreeEnergy();
			} else {
				profile = MolQ->GetGrandPotentialProfile();
				AddChemPotGraft(profile,Lat);
				F_system_N=GetFreeEnergyPo(Lat);
			}
			for (int p=1; p<=numPos; p++) {
				r_last[p]=r_new[p];
				int rr;
				double distance=0;
				int l=0; double norm=0;
				for (int i=-1; i<2; i++)
				for (int j=-1; j<2; j++)
				for (int k=-1; k<2; k++) {
					rr=jx*(x_new[p]-i) +jy*(y_new[p]-j)+jz*(z_new[p]-k)+1;
					distance = pow((rx_new[p]-x_new[p]),2)+pow((ry_new[p]-y_new[p]),2)+pow((rz_new[p]-z_new[p]),2);
					if (free[jx*(x_new[p]-i) +jy*(y_new[p]-j)+jz*(z_new[p]-k)+1]) submask[l]=exp(-distance); else submask[l]=0.0;
					norm+=submask[l]; l++;
				}
				l=0; F_new[p]=0;
				for (int i=-1; i<2; i++)
				for (int j=-1; j<2; j++)
				for (int k=-1; k<2; k++) {
					rr=jx*(x_new[p]-i) +jy*(y_new[p]-j)+jz*(z_new[p]-k)+1;
					F_new[p] += profile[rr]*submask[l]/norm; l++;
				}
			}
		}	else {
			write_output = false;
			ldt--;
			for (int k=1; k<= numPos; k++){
				x_new[k]=x_old[k];
				y_new[k]=y_old[k];
				z_new[k]=z_old[k];
				r_new[k]=r_old[k];
			}
			for (int k=1; k<= numPos; k++){
				rx_new[k]=rx_old[k];
				ry_new[k]=ry_old[k];
				rz_new[k]=rz_old[k];
			}
			F_system_N= F_system_O;
		}

		if (write_output) {
			//double FN;

			SetRho(Lat,rho_0); SetF(Lat,F_0);
			SetFpo(Lat,Fpo_0); SetOmega(Lat,Omega_0);
			if (SegQ->Charged()) SetPsi(Lat,PSI_0);

			GetOutput(Out,Lat);
			if (ldt==1) Out->WriteOutput(); else Out->WriteMCOutput(ldt);
			Out->Clear();

			//profile = MolQ->GetFreeEnergyProfile();

			//if (LD_F) {FN=MolQ->GetFreeEnergy();
			//} else {
			//	FN=GetFreeEnergyPo(Lat);
			//}

			cout << "LDT "<< ldt << " F = " << F_system_N  << endl;
			for (int k=1; k<=numPos; k++) {
					//cout << "("<< x_new[k] << ","<< y_new[k] <<","<< z_new[k] <<")";
					cout << "("<< rx_new[k] << ","<< ry_new[k] <<","<< rz_new[k] <<")";
			} cout << endl;
		}
	}
}

void
SF_SystemLD::GetOutput(Output* Out,Lattice* Lat) const {
	Out->PutText("sys",name,"inputfile",MyInput->GetFileName());
	Out->PutText("sys",name,"calculation_type","LD_SCF");
	if (approach == firstOrder) {
		Out->PutText("sys",name,"matrix_approach","first_order");
	}
	if (approach == secondOrder) {
			Out->PutText("sys",name,"matrix_approach","second_order");
	}
	Out->PutBoolean("sys",name,"overflow_protection",overflowProtection);
	Out->PutReal("sys",name,"temperature",TEMPERATURE);
	Out->PutReal("sys",name,"free energy",MolQ->GetFreeEnergy());
	Out->PutReal("sys",name,"excess free energy",MolQ->GetExcessFreeEnergy());
	if (LD_F) Out->PutText("sys",name,"free_energy_type","free_energy"); else Out->PutText("sys",name,"free_energy_type","free_energy(po)");
	Out->PutInt("sys",name,"ldt", ldt);
	Out->PutInt("sys",name,"LDT", LDT);
	Out->PutReal("sys",name,"mobility",mobility);

	if (SegQ->Charged()) {
		Out->PutProfile("sys",name,"potential",SegQ->GetElectricPotential());
		Out->PutProfile("sys",name,"charge",SegQ->GetCharge());
		Out->PutProfile("sys",name,"epsilon",SegQ->GetAverageEpsilon());
	}

	GetLDOutput(Out,Lat);

	if (MyInput->ValueSet("sys",name,"Probe_molecule")) GetEnd_to_EndOutput(Out,Lat);

	Lat1st->GetOutput(Out);
	MolQ->GetOutput(Out);
	SegQ->GetOutput(Out);
	ReactionQ->GetOutput(Out);
}

void
SF_SystemLD::GetEnd_to_EndOutput(Output* Out,Lattice* Lat) const{

	MyInput->AlwaysCombineParam("sys",name,"Probe_molecule","extra_output");
	Text ProbeMolName = MyInput->GetText("sys",name,"Probe_molecule");
	if (!MolQ->MoleculeDefined(ProbeMolName)) {Message(literal,"'sys : " + name + " : ProbeMolName : " + ProbeMolName + " is a molecule that is not defined aborting the end-to-end output"); return;}
	SF_Molecule* ProbeMol;
	ProbeMol = MolQ->GetMolecule(ProbeMolName);
	SF_MolStructure* Chain = ProbeMol->GetMolStructure();
	int N = ProbeMol->GetChainLength();
	if (N<2) {Message(literal,"Probe chain is too short: aborting the end-to-end output"); return;}
	int segnr=1; while (Chain->GetSegment(1) != Chain->GetDiffSegment(segnr)) segnr++;
	Vector phi1 = Chain->GetDiffSegment(segnr)->GetPhi(total);
	segnr=1; while (Chain->GetSegment(N) != Chain->GetDiffSegment(segnr)) segnr++;
	Vector phiN = Chain->GetDiffSegment(segnr)->GetPhi(total);
	if (Chain->GetSegment(1)->GetFreedom()==pinned) {

		//int n_layers_x = Lat->GetNumLayers(1)-2;
		//int n_layers_y = Lat->GetNumLayers(2)-2;
		//int n_layers_z = Lat->GetNumLayers(3)-2;
		//int jx = (n_layers_y+2)*(n_layers_z+2);
		//int jy = n_layers_z+2;
		int r;
		double CoMx =0.0,CoMy =0.0,CoMz =0.0;
		double tM1 =0.0;
		double tMN =0.0;
		double ete=0.0;
		double de = 0.0;
		double r4,r2;
		for (int x=1; x<=n_layers_x; x++)
		for (int y=1; y<=n_layers_y; y++)
		for (int z=1; z<=n_layers_z; z++) {
			r=jx*x+jy*y+z+1;
			tM1 +=phi1[r];
			tMN +=phiN[r];
			CoMx += x*phi1[r];
			CoMy += y*phi1[r];
			CoMz += z*phi1[r];
		}
		Out->PutReal("mol",ProbeMol->GetName(),"pos(1,x)",CoMx/tM1);
		Out->PutReal("mol",ProbeMol->GetName(),"pos(1,y)",CoMy/tM1);
		Out->PutReal("mol",ProbeMol->GetName(),"pos(1,z)",CoMz/tM1);

		for (int x=1; x<=n_layers_x; x++)
		for (int y=1; y<=n_layers_y; y++)
		for (int z=1; z<=n_layers_z; z++) {
			r=jx*x+jy*y+z+1;
			ete += (pow(x-CoMx,2)+pow(y-CoMy,2)+pow(z-CoMz,2))*phiN[r];
			de +=pow(pow(x-CoMx,2)+pow(y-CoMy,2)+pow(z-CoMz,2),2)*phiN[r];
		}
		r2=ete/tMN;
		r4=de/tMN;
		Out->PutReal("mol",ProbeMol->GetName(),"<r(ee)**2>**0.5",pow(r2,0.5));
		Out->PutReal("mol",ProbeMol->GetName(),"<r(ee)**4>**0.25",pow(r4,0.25));
		Out->PutReal("mol",ProbeMol->GetName(),"dr =(<r4>-<r2>**2)**0.5/<r2>",pow(r4-r2*r2,0.5)/r2);

	} else {Message(literal,"Probe chain is not pinned with first segment: aborting the end-to-end output"); return;}
}

void
SF_SystemLD::GetLDOutput(Output* Out,Lattice* Lat) const{
	//int n_layers_x = Lat->GetNumLayers(1)-2;
	//int n_layers_y = Lat->GetNumLayers(2)-2;
	//int n_layers_z = Lat->GetNumLayers(3)-2;
	//int jx = (n_layers_y+2)*(n_layers_z+2);
	//int jy = n_layers_z+2;
	//int jz = 1;
	int numSeg = SegQ->GetNumSegments();
	int position,pos_x,pos_y,pos_z;
	LatticeRange* LatRange;
	Text dim;
	Vector profile;
	double pro[n_layers_z+1];
	Out->PutReal("sys",name,"LD_Free_energy",F_system_N);
	for (int i=1; i<=MolQ->GetNumMolecules(); i++) {
		SF_Molecule* Mol = MolQ->GetMolecule(i);
		for (int z=1; z<=n_layers_z; z++) {
			dim = Blanks(9);
			dim.Putint(z);
			dim = Copy(dim.Strip());
			dim = Copy(dim.Frontstrip());
			Out->PutReal("mol",Mol->GetName(),"rho["+ dim + "]",rho_0[i][z]);
		}
	}
	for (int z=1; z<=n_layers_z; z++) {
		dim = Blanks(9);
		dim.Putint(z);
		dim = Copy(dim.Strip());
		dim = Copy(dim.Frontstrip());
		Out->PutReal("sys",name,"F["+ dim + "]",F_0[z]);
	}
	for (int z=1; z<=n_layers_z; z++) {
		dim = Blanks(9);
		dim.Putint(z);
		dim = Copy(dim.Strip());
		dim = Copy(dim.Frontstrip());
		Out->PutReal("sys",name,"Omega["+ dim + "]",Omega_0[z]);
	}
	for (int z=1; z<=n_layers_z; z++) {
		dim = Blanks(9);
		dim.Putint(z);
		dim = Copy(dim.Strip());
		dim = Copy(dim.Frontstrip());
		Out->PutReal("sys",name,"Fpo["+ dim + "]",Fpo_0[z]);
	}
	if (SegQ->Charged()) {
		for (int z=1; z<=n_layers_z; z++) {
			dim = Blanks(9);
			dim.Putint(z);
			dim = Copy(dim.Strip());
			dim = Copy(dim.Frontstrip());
			Out->PutReal("sys",name,"psi["+ dim + "]",PSI_0[z]);
		}
	}
	int part_pos=0;
	double mean_x=0;
	double mean_y=0;
	double mean_z=0;

	for (int i=1; i<=numSeg; i++) {
		SF_Segment* Seg = SegQ->GetSegment(i);
		if (Seg->GetLD()) {
			LatRange=Seg->GetLatRange();
			int num_seg_pos = LatRange->GetNumPos();
			for (int k=1; k<=num_seg_pos; k++){
				part_pos++;
				dim = Blanks(9);
				dim.Putint(k);
				dim = Copy(dim.Strip());
				dim = Copy(dim.Frontstrip());
				position=LatRange->GetPos(k);
				Out->PutInt("mon",Seg->GetName(),"pos["+ dim + "]",position);
				position = position -1;
				pos_z = ((position % jx) % jy)/jz; 	 mean_z +=pos_z;
				pos_y = ((position % jx)-pos_z*jz)/jy;	 mean_y +=pos_y;
				pos_x = (position-pos_y*jy-pos_z*jz)/jx;   mean_x +=pos_x;
				Out->PutInt("mon",Seg->GetName(),"pos_x["+ dim + "]",pos_x);
				Out->PutInt("mon",Seg->GetName(),"pos_y["+ dim + "]",pos_y);
				Out->PutInt("mon",Seg->GetName(),"pos_z["+ dim + "]",pos_z);
				if (part_pos < numPos+1) { //gaat fout met meerdere monomeren die een pinning stamp hebben
					Out->PutInt("mon",Seg->GetName(),"r[" + dim + "]_x",x_new[part_pos]);
					Out->PutInt("mon",Seg->GetName(),"r[" + dim + "]_y",y_new[part_pos]);
					Out->PutInt("mon",Seg->GetName(),"r[" + dim + "]_z",z_new[part_pos]);
				}
			}
		}
	}
	mean_x = mean_x/numPos; mean_y = mean_y/numPos; mean_z = mean_z/numPos; //estimate of the center of mass.
	double mx=0.0, my=0.0, mz=0.0, mmx=0.0, mmy=0.0, mmz=0.0;
	double R_RMS=0.0;
	for (int k=1; k<=numPos; k++) {
		mmx=x_new[k];	mmy=y_new[k];	mmz=z_new[k];
		if (x_new[k]-mean_x > n_layers_x/2)    mmx = x_new[k]-n_layers_x;
		if (mean_x-x_new[k] > n_layers_x/2)    mmx = x_new[k]+n_layers_x;
		if (y_new[k]-mean_y > n_layers_y/2)    mmy = y_new[k]-n_layers_y;
		if (mean_y-y_new[k] > n_layers_y/2)    mmy = y_new[k]+n_layers_y;
		if (z_new[k]-mean_z > n_layers_z/2)    mmz = z_new[k]-n_layers_z;
		if (mean_z-z_new[k] > n_layers_z/2)    mmz = z_new[k]+n_layers_z;
		mx +=mmx;
		my +=mmy;
		mz +=mmz;
	}
	mean_x = mx/numPos; mean_y = my/numPos; mean_z=mz/numPos; //improved center of mass
	Out->PutReal("sys",name,"COM_x",mean_x);
	Out->PutReal("sys",name,"COM_y",mean_y);
	Out->PutReal("sys",name,"COM_z",mean_z);
	for (int k=1; k<=numPos; k++) {
		mx = x_new[k]; my = y_new[k]; mz = z_new[k];
		if (x_new[k]-mean_x > n_layers_x/2)    mx = x_new[k]-n_layers_x;
		if (mean_x-x_new[k] > n_layers_x/2)    mx = x_new[k]+n_layers_x;
		if (y_new[k]-mean_y > n_layers_y/2)    my = y_new[k]-n_layers_y;
		if (mean_y-y_new[k] > n_layers_y/2)    my = y_new[k]+n_layers_y;
		if (z_new[k]-mean_z > n_layers_z/2)    mz = z_new[k]-n_layers_z;
		if (mean_z-z_new[k] > n_layers_z/2)    mz = z_new[k]+n_layers_z;
		R_RMS +=(mx-mean_x)*(mx-mean_x)+(my-mean_y)*(my-mean_y)+(mz-mean_z)*(mz-mean_z);
	}

	R_RMS = pow(R_RMS/numPos,0.5);// average distance of crosslink from center of mass
	Out->PutReal("sys",name,"R_RMS",R_RMS);
	for (int k=0; k<=n_layers_z; k++) pro[k]=0;
	int distance;
	for (int k=1; k<=numPos; k++) {
		mx = x_new[k]; my = y_new[k]; mz = z_new[k];
		if (x_new[k]-mean_x > n_layers_x/2)    mx = x_new[k]-n_layers_x;
		if (mean_x-x_new[k] > n_layers_x/2)    mx = x_new[k]+n_layers_x;
		if (y_new[k]-mean_y > n_layers_y/2)    my = y_new[k]-n_layers_y;
		if (mean_y-y_new[k] > n_layers_y/2)    my = y_new[k]+n_layers_y;
		if (z_new[k]-mean_z > n_layers_z/2)    mz = z_new[k]-n_layers_z;
		if (mean_z-z_new[k] > n_layers_z/2)    mz = z_new[k]+n_layers_z;
		distance = static_cast<int>(pow((mx-mean_x)*(mx-mean_x)+(my-mean_y)*(my-mean_y)+(mz-mean_z)*(mz-mean_z),0.5));
		if (distance <= n_layers_z) pro[distance]++;
	}
	for (int k=0; k<=n_layers_z; k++) {
		dim = Blanks(9);
		dim.Putint(k);
		dim = Copy(dim.Strip());
		dim = Copy(dim.Frontstrip());
		Out->PutReal("sys",name,"N_link["+ dim + "]",pro[k]); //radial cross-link profile (not normalized);

	}
}


void
SF_SystemLD::SetF(Lattice* Lat, Vector  F_0) {
	//int n_layers_x = Lat->GetNumLayers(1)-2;
	//int n_layers_y = Lat->GetNumLayers(2)-2;
	//int n_layers_z = Lat->GetNumLayers(3)-2;
	//int jx = (n_layers_y+2)*(n_layers_z+2);
	//int jy = n_layers_z+2;
	double area = n_layers_x*n_layers_y;
	Vector profile;
	profile = MolQ->GetFreeEnergyProfile();
	for (int z=1; z<=n_layers_z; z++) {
		F_0[z]=0;
		for (int x=1; x<=n_layers_x; x++){
			for (int y=1; y<=n_layers_y; y++){
				F_0[z] = F_0[z]+ profile[jx*x+jy*y+z+1];
			}
		}
		F_0[z] = F_0[z]/area;
	}
}

void
SF_SystemLD::SetFpo(Lattice* Lat, Vector Fpo_0) {
	//int n_layers_x = Lat->GetNumLayers(1)-2;
	//int n_layers_y = Lat->GetNumLayers(2)-2;
	//int n_layers_z = Lat->GetNumLayers(3)-2;
	//int jx = (n_layers_y+2)*(n_layers_z+2);
	//int jy = n_layers_z+2;
	double area = n_layers_x*n_layers_y;
	Vector profile;
	profile = MolQ->GetGrandPotentialProfile();
	AddChemPotGraft(profile,Lat);
	for (int z=1; z<=n_layers_z; z++) {
		Fpo_0[z]=0;
		for (int x=1; x<=n_layers_x; x++){
			for (int y=1; y<=n_layers_y; y++){
				Fpo_0[z] = Fpo_0[z]+ profile[jx*x+jy*y+z+1];
			}
		}
		Fpo_0[z] = Fpo_0[z]/area;
	}
}
void
SF_SystemLD::SetOmega(Lattice* Lat, Vector Omega_0) {
	//int n_layers_x = Lat->GetNumLayers(1)-2;
	//int n_layers_y = Lat->GetNumLayers(2)-2;
	//int n_layers_z = Lat->GetNumLayers(3)-2;
	//int jx = (n_layers_y+2)*(n_layers_z+2);
	//int jy = n_layers_z+2;
	double area = n_layers_x*n_layers_y;
	Vector profile;
	profile = MolQ->GetGrandPotentialProfile();
	for (int z=1; z<=n_layers_z; z++) {
		Omega_0[z]=0;
		for (int x=1; x<=n_layers_x; x++){
			for (int y=1; y<=n_layers_y; y++){
				Omega_0[z] = Omega_0[z]+ profile[jx*x+jy*y+z+1];
			}
		}
		Omega_0[z] = Omega_0[z]/area;
	}
}

void
SF_SystemLD::SetPsi(Lattice* Lat, Vector PSI_0) {
	//int n_layers_x = Lat->GetNumLayers(1)-2;
	//int n_layers_y = Lat->GetNumLayers(2)-2;
	//int n_layers_z = Lat->GetNumLayers(3)-2;
	//int jx = (n_layers_y+2)*(n_layers_z+2);
	//int jy = n_layers_z+2;
	double area = n_layers_x*n_layers_y;
	Vector profile;
	profile = SegQ->GetElectricPotential();
	for (int z=1; z<=n_layers_z; z++) {
		PSI_0[z]=0;
		for (int x=1; x<=n_layers_x; x++){
			for (int y=1; y<=n_layers_y; y++){
				PSI_0[z] = PSI_0[z]+ profile[jx*x+jy*y+z+1];
			}
		}
		PSI_0[z] = PSI_0[z]/area;
	}
}

void
SF_SystemLD::SetRho(Lattice* Lat, Array<Vector>  rho_0) const {
	//int n_layers_x = Lat->GetNumLayers(1)-2;
	//int n_layers_y = Lat->GetNumLayers(2)-2;
	//int n_layers_z = Lat->GetNumLayers(3)-2;
	//int jx = (n_layers_y+2)*(n_layers_z+2);
	//int jy = n_layers_z+2;
	double area = n_layers_x*n_layers_y;
	Vector profile;
	for (int i=1; i<=MolQ->GetNumMolecules(); i++) {
		SF_Molecule* Mol = MolQ->GetMolecule(i);
		profile = Mol->GetPhi(total);
		for (int z=1; z<=n_layers_z; z++) {
			rho_0[i][z]=0;
			for (int x=1; x<=n_layers_x; x++){
				for (int y=1; y<=n_layers_y; y++){
					rho_0[i][z] = rho_0[i][z]+ profile[jx*x+jy*y+z+1];
				}
			}
			rho_0[i][z] = rho_0[i][z]/area;
		}
	}
}
