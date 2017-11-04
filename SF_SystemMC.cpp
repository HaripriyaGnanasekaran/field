#include "SF_SystemMC.h"

SF_SystemMC::SF_SystemMC(Input* MyInput_, Text name_, Boolean compute_) {
	MyInput = MyInput_;
	name = name_;
	compute = compute_;
	Array<Text> param(1,21);
	param[1] = "calculation_type";
	param[2] = "matrix_approach";
	param[3] = "overflow_protection";
	param[4] = "temperature";
	param[5] = "warnings";
	param[6] = "MCS";
	param[7] = "move_strategy";
	param[8] = "d_pos_max";
	param[9] = "num_of_moves";
	param[10] = "lattice_moves";
	param[11] = "free_energy_type";
	param[12] = "MC_temperature";
	param[13] = "cluster_move";
	param[14] = "teleportfreq_Z"; // determines the chance that a particle is moved to a random z position instead of the normal MC move
	param[15] = "swapfreq";
	param[16] = "equilibrationsteps";
	param[17] = "extra_output";
	param[18] = "allow_overlap";
	param[19] = "probe_chain";
	param[20] = "loop_count";
	param[21] = "random_seed";

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
	mc_temperature = MyInput->GetReal("sys",name,"MC_temperature",0.0001,DBL_MAX,temperature);

	if (approach == firstOrder) {
		Text latName;
		MyInput->GetNumNames("lat",1,1);
		Array<Text> latNames = MyInput->GetNames("lat");
		latName = latNames[1];
		Lat1st = NewLat1stO(latName);
		SegQ= new SF_SegmentList(Lat1st,MyInput);
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
		SegQ= new SF_SegmentList(Lat2nd,MyInput);
		MolQ = new SF_MolList2ndO(SegQ,Lat2nd,MyInput);
		ReactionQ = new SF_ReactionList(MolQ,SegQ,MyInput);
		Solve = NewSolve(ReactionQ,MolQ,SegQ,Lat2nd);
	}

	if (changeChiWithTemperature) {
		UpdateChi();
	}
	TEMPERATURE = temperature;

	bool MC_found=false;
	for (int i=1; i<=MolQ->GetNumMolecules(); i++) {
		SF_Molecule* Mol = MolQ->GetMolecule(i);
		if (Mol->GetMolStructure()->SomeSegmentsMC()) {
			MC_found = true;
		}
	}

	if (!MC_found) {
		for (int i=1; i<=SegQ->GetNumSegments(); i++) {
		SF_Segment* Seg = SegQ->GetSegment(i);
		    if (Seg->GetMC() && Seg->GetFreedom()==frozen) MC_found = true;
		}
	}

	if (!MC_found) {
		Message(fatal,"No molecule with pinned segments with activated MC found");
	}


	if (!MyInput->ValueSet("sys", name, "MCS") ) Message(fatal,"For MC_SCF you need to set MCS value (typically 10 000)");
	else MCS = MyInput->GetInt("sys",name,"MCS",1,1000000);

	if (!MyInput->ValueSet("sys", name, "d_pos_max") ) {
		Message(literal,"For MC_SCF you should set 'd_pos_max' value (default 1 is used)"); d_max_pos = 1;
	} else d_max_pos = MyInput->GetInt("sys",name,"d_pos_max",1,100,1);

	if (!MyInput->ValueSet("sys", name, "lattice_moves") ) {
		Message(literal,"For MC_SCF you are advised to set 'lattice_moves' value (default 'true' is used)");
		lat_moves = true;
	} else lat_moves = MyInput->GetBoolean("sys",name,"lattice_moves",true);
	if (MyInput->ValueSet("sys", name, "equilibrationsteps") ) n_equi_steps = MyInput->GetInt("sys",name,"equilibrationsteps",-1,INT_MAX,-1);
	if (MyInput->ValueSet("sys",name,"move_strategy")) {
		Array<Text> strategy(1,2);
		strategy[1] = "fixed";
		strategy[2] = "variable";
		opt_mc = MyInput->GetChoice("sys",name,"move_strategy",strategy,2)==2 ;
	} else {
		opt_mc = true;
		}

	n_moves=0;
	if (MyInput->ValueSet("sys", name, "num_of_moves") ) n_moves = MyInput->GetInt("sys",name,"num_of_moves",1,9999,1);
	if ((!opt_mc || n_equi_steps==0 ) && n_moves==0) Message(fatal, "When 'move_strategy = fixed or equilibrationsteps = 0' one should set 'num_of_moves'");



	if (MyInput->ValueSet("sys",name,"free_energy_type")) {
	 	Array<Text> FChoices(1,2);
		FChoices[1] = "free_energy";
		FChoices[2] = "free_energy(po)";

		MC_F = MyInput->GetChoice("sys",name,"free_energy_type",FChoices)==1;
	} else Message(fatal,"free_energy_type should be specified when MC_SCF option is selected; set this quantity to 'free_energy', or 'free_energy(po)'");

	swapfreq = MyInput->GetReal("sys",name,"swapfreq",0.0,1.0,0.0);

	Cluster=false; neighbor=0;

	if (MyInput->ValueSet("sys",name,"cluster_move")){ Cluster=MyInput->GetBoolean("sys",name,"cluster_move",true);}
	teleportfreq_Z = MyInput->GetReal("sys",name,"teleportfreq_Z",0.0,1.0,0.0);
	if (teleportfreq_Z + swapfreq > 1.0){Message(fatal,"The sum of the swapfreq and teleportfreq_Z should not be larger than one.");}

}


double
SF_SystemMC::random_d() {
	return rand()/(RAND_MAX+1.0);
}

int
SF_SystemMC::random_int(int low, int high) {
	int range = high-low+1;
	return low+int(range*random_d());
}
double
SF_SystemMC::random(double low, double high) {
	double range = high-low;
	return low+range*random_d();
}


bool // checks wether two particles have the same position and returns false if this is the case;
SF_SystemMC::CheckOverlap(int R){
	int xi,xj,yi,yj,zi,zj,rr;
	double distance;
	bool allowoverlap = MyInput->GetBoolean("sys",name,"allow_overlap",false);
	if (!allowoverlap){
		for (int i=1; i <= numPos-1; i++) {
			for (int j=i+1; j <= numPos; j++) {
				rr=r_new[i]-1;
				zi = ((rr % jx) % jy)/jz;
				yi = ((rr % jx)-zi*jz)/jy;
				xi = (rr-yi*jy-zi*jz)/jx;
				rr=r_new[j]-1;
				zj = ((rr % jx) % jy)/jz;
				yj = ((rr % jx)-zj*jz)/jy;
				xj = (rr-yj*jy-zj*jz)/jx;
				distance = pow(pow(xi-xj,2.0)+pow(yi-yj,2.0)+pow(zi-zj,2.0),0.5);
				if (distance<R/2) {return false; }
			}
		}
	}

	// check wether there is enough space for all pinned segments

	// Step 1: Determine the total amounts of each segment type that is pinned and the number of lattice sites within the pinned range.
	int numSeg = SegQ->GetNumSegments(); // determine the number of different kinds of segments.
	int numMol = MolQ->GetNumMolecules(); // determine the number of different kinds of molecules
	int Nseginmol ;
	std::vector<int> Nposmax;// list containg the number of voxels of the pinned ranges for each segment type.
	int Npos;
	int Ntotal;
	double theta;
	double Nsegtotal;
	int NumMC =0;
	int lowlimit;
	std::vector<double> segNlist; //list containing the amount of each segment
	std::vector<int> pinnedseglist; // list containing the indexes of of the pinned segments  in SegQ
	int numpinnedseg = 0;

	for (int m=1;m <=numSeg; m++){	// loop over all segment types
		SF_Segment* Seg = SegQ->GetSegment(m);
		if (Seg->GetFreedom() == pinned){ // check whether this segment type is pinned.
			Text SegName = Seg->GetName();
			Nsegtotal = 0;
			LatticeRange* LatRange=Seg->GetLatRange();
			Npos = LatRange->GetNumPos();// determine the number of pinned sites
			Nposmax.push_back(Npos);
			if (!allowoverlap){ // check whether two monte carlo nodes of the same type overlap.
				if (Seg->GetMC()){	// If so the number of pinned sites found is less than the number of monte carlo sites times the spotsize
					NumMC++;
					if (NumMC == 1)	{
						lowlimit = 0;
					}
					else {
						lowlimit = numposseg[NumMC-2];
					}
					int spot = MyInput->GetInt(SEGMENT,SegName,"spot_size",1,27,1);
					if (Npos < spot*(numposseg[NumMC-1]-lowlimit)){
						Message(literal,"Segment "+ (SegQ->GetSegment(m))->GetName()+ " overlaps with other monte carlo nodes of the same type.");
						return false;
					}
				}
			}
			numpinnedseg++; //count the number of segment types that are pinned.
			for (int i=1; i<=numMol; i++) { //loop over all molecules
				SF_Molecule* Mol=MolQ->GetMolecule(i);
				SF_MolStructure* Chain = Mol->GetMolStructure();
				if (Chain->SomeSegmentsPinned()){
					int numDiffSeg = Chain->GetNumDiffSegments();
					for (int n=1; n<=numDiffSeg; n++){
						SF_MolSegment* MolSeg = Chain->GetDiffSegment(n);
						if (MolSeg->GetName() == SegName){
							Nseginmol = Chain->GetAvNumSegments(MolSeg); // count the number of segments of segment type Seg in the molecule
							Ntotal = Mol->GetChainLength(); // Determine the total length of the chain
							theta = Mol->GetTheta(); // get the total amount of this molecule
							Nsegtotal =+ theta*Nseginmol/Ntotal; // add the amount of the segment in this molecule to the total amount.
						}
					}
				}
			}
			segNlist.push_back(Nsegtotal);
			pinnedseglist.push_back(m);
		}
	}


	// step 2: Determine the amount of overlap between pinned and frozen ranges
	int frozenoverlap[numpinnedseg];
	int Overlaparray[numpinnedseg][numpinnedseg];
	for (int q=0;q<=numpinnedseg-1;q++){
		for (int r=0;r<=numpinnedseg-1;r++){
			Overlaparray[q][r]=0;// fill overlaparray with 0.
		}
	}

	int position1;
	int position2;
	for (int m=1; m <= numpinnedseg;m++){// loop over all pinned segment types
		SF_Segment* Seg = SegQ->GetSegment(pinnedseglist[m-1]);
		LatticeRange* LatRange1=Seg->GetLatRange();
		//Determine the total number of coordinates in the pinnend range for this segment
		frozenoverlap[m-1]=0;
		for(int n=1;n <= Nposmax[m-1];n++){// loop over all locations in the pinned range
			position1= LatRange1->GetPos(n); // get the position of the n th pinned segment

			if (!free[position1]){
				if (!allowoverlap){
					Message(literal,"Segment "+ (SegQ->GetSegment(pinnedseglist[m-1]))->GetName()+ " overlaps with frozen range.");
					return false;
				}
				frozenoverlap[m-1]++;  // check how many of these coordinates overlap with frozen ranges
			}
			for (int o=m;o <= numpinnedseg-1;o++){ // loop over all other pinned segment types
				SF_Segment* Seg = SegQ->GetSegment(pinnedseglist[o]);
				LatticeRange* LatRange2=Seg->GetLatRange();
				for(int p=1;p <= Nposmax[o];p++){ // loop over the locations of the pinned sites of other pinned segmants
					position2 = LatRange2->GetPos(p);
					if (position1 == position2){ // check whether two pinned segments overlap
						if (!allowoverlap){
							Message(literal,"Segment "+ (SegQ->GetSegment(pinnedseglist[m-1]))->GetName()+ " overlaps with pinned range of Segment " +(SegQ->GetSegment(pinnedseglist[o]))->GetName()+".");
							return false;

						}
						Overlaparray[m-1][o]++;
					}
				}
			}
		}
	}
	for (int m=0; m<=numpinnedseg-2;m++){
		for (int o=m+1;o<=numpinnedseg-1;o++){
			Overlaparray[o][m] = Overlaparray[m][o]; // Fill other half overlap array
		}
	}
	//step 4: Determine wether there is enough space for all pinned segments
	// First make a rough estimate whether there is enough space because exact calculation takes a lot of time.
	int overlapsum;
	for (int m=0; m<=numpinnedseg-1;m++){ // loop over all segment types
		overlapsum = 0;
		for (int o=0;o<=numpinnedseg-1;o++){
			overlapsum += Overlaparray[o][m]; // add all the overlapping regions together
		}
		if (Nposmax[m]-overlapsum-frozenoverlap[m]< segNlist[m]){ // check wether there is enough space left after subtracting all overlap.
			if (Nposmax[m]-frozenoverlap[m]< segNlist[m]){ // check if there is enough space after subtracting frozen range from pinned range. If not the segments will not fit.
				if (!MyInput->ValueSet(SEGMENT,SegQ->GetSegment(pinnedseglist[m])->GetName(),"spot_mix")){
					Message(literal,"Not enough space for pinned segments "+ (SegQ->GetSegment(pinnedseglist[m]))->GetName()+ ".");
					cout << " Nposmax " << Nposmax[m] << " segNlist[m] "<< endl;
					return false;
				}
			}
			Message(literal,"Possibly not enough space for pinned segments "+ (SegQ->GetSegment(pinnedseglist[m]))->GetName() + ".");
		}
	}
	return true;
}

//Randomly selects a particle and moves it by a random amount in a random direction.  n_move specifies the number of times a particle is moved;
// The function returns false if a particle has been placed on an other particle or within a frozen range. Other wise it will return true;
// When applicable, i.e., clustermoves, a closest neighbor is identified which is moved in the same way as the randomly selected particle.

bool
SF_SystemMC::MovePos(int n_move,SF_Segment** Seg, int neighbor){

	if (!lat_moves) return Off_Lattice_MovePos(n_move,Seg,neighbor);

	for (int i=1; i <= n_move; i++){

		int delta =0, delta_x=0,delta_y=0,delta_z=0,rndpart,rndpart2,neighpart=0,swapmove = 0,rstorage,direction;
		double distance,closest=pow(neighbor+1,2.0);
		rndpart = random_int(1,numPos);

		if ((double)rand() / RAND_MAX < teleportfreq_Z){

			direction = 3;
			delta = random_int(1,n_layers_z-1);//by choosing n_layers_z-1 as the highest value the new z position will always be different from the original z position.
		}
		else if ((double)rand() / RAND_MAX < (swapfreq/(1-teleportfreq_Z))){ // do a swap move

			swapmove = 1;
			direction = 1;
			bool successfullswap = false;
			while (!successfullswap){
				rndpart2 = random_int(1,numPos); // randopmly select a particle
				int rndrange1=-1;
				int rndrange2=-1;
				int j = 0;
				while ((rndrange1 < 0) || (rndrange2 <0)){ // while for one of the particles the segmnet type has not yet been determined
					if ((rndpart <= numposseg[j]) & (rndrange1 <0)) rndrange1=j; // check is the particle belongs to the first j+1 types of segment.
					if ((rndpart2 <= numposseg[j]) & (rndrange2 <0)) rndrange2=j;
					j++;
				}
				if (rndrange1 != rndrange2){ 	// check that both particles areof a different segmenttype.
					rstorage = r_new[rndpart];
					r_new[rndpart]= r_new[rndpart2];
					r_new[rndpart2]= rstorage;
					successfullswap = true;
				}
			}
		}
		else {

			direction = random_int(1,3);
			while (delta == 0){
				delta = random_int(-1*d_max_pos,d_max_pos);// Due to the poor random number generator there can be apreference for moving in one of the directions. Although this preference should be small ~ 1/10923 for dmaxpos=1
			}
		}

		if (neighbor>0) { //alleen clustermoves when neighbor distance is defined.
			while (delta == 0){
				delta = random_int(-1*d_max_pos,d_max_pos);
			}
			if ((double)rand()>0.5) { //50% van de keren clustermoves implementeren.
				for (int j=1; j<=numPos; j++) {
					if (j != rndpart) {
						distance = (x_new[j]-x_new[rndpart])*(x_new[j]-x_new[rndpart])+(y_new[j]-y_new[rndpart])*(y_new[j]-y_new[rndpart])+(z_new[j]-z_new[rndpart])*(z_new[j]-z_new[rndpart]);
						if (distance < closest) {closest = distance; neighpart=j;}
					}
				}
				if (closest > neighbor*neighbor) neighpart =0; //als closest particle niet dichterbij zit dan de neighbordistance dan is er geen sprake van een cluster.
			}

		}


		if (swapmove == 0){
			if (direction == 1 && Seg[rndpart]->MoveAllowed(direction)) delta_x= delta;
			if (direction == 2 && Seg[rndpart]->MoveAllowed(direction)) delta_y= delta;
			if (direction == 3 && Seg[rndpart]->MoveAllowed(direction)) delta_z= delta;
			x_new[rndpart]=x_new[rndpart]+delta_x;
			y_new[rndpart]=y_new[rndpart]+delta_y;
			z_new[rndpart]=z_new[rndpart]+delta_z;


			if (x_new[rndpart]<1) x_new[rndpart]=x_new[rndpart]+n_layers_x;
			if (y_new[rndpart]<1) y_new[rndpart]=y_new[rndpart]+n_layers_y;
			if (z_new[rndpart]<1) z_new[rndpart]=z_new[rndpart]+n_layers_z;
			if (x_new[rndpart]>n_layers_x) x_new[rndpart]=x_new[rndpart]-n_layers_x;
			if (y_new[rndpart]>n_layers_y) y_new[rndpart]=y_new[rndpart]-n_layers_y;
			if (z_new[rndpart]>n_layers_z) z_new[rndpart]=z_new[rndpart]-n_layers_z;

			r_new[rndpart]=jx*x_new[rndpart]+jy*y_new[rndpart]+z_new[rndpart]+1;


			if (neighpart > 0) { //als neighbor particle has been identified, verzet het deeltje dan op dezelfde manier als rndpart.
				x_new[neighpart]=x_new[neighpart]+delta_x;
				y_new[neighpart]=y_new[neighpart]+delta_y;
				z_new[neighpart]=z_new[neighpart]+delta_z;
				if (x_new[neighpart]<1) x_new[neighpart]=x_new[neighpart]+n_layers_x;
				if (y_new[neighpart]<1) y_new[neighpart]=y_new[neighpart]+n_layers_y;
				if (z_new[neighpart]<1) z_new[neighpart]=z_new[neighpart]+n_layers_z;
				if (x_new[neighpart]>n_layers_x) x_new[neighpart]=x_new[neighpart]-n_layers_x;
				if (y_new[neighpart]>n_layers_y) y_new[neighpart]=y_new[neighpart]-n_layers_y;
				if (z_new[neighpart]>n_layers_z) z_new[neighpart]=z_new[neighpart]-n_layers_z;
				r_new[neighpart]=jx*x_new[neighpart]+jy*y_new[neighpart]+z_new[neighpart]+1;
			}

			if (!free[r_new[rndpart]]) {return false;}
			if (neighpart > 0) { if (!free[r_new[neighpart]]) {return false;}}
		}
	}
	// reset mask file and positions before checking whether there is enough space for the pinned segments

	SF_Segment* Seg2[numPos+1];
	SF_Segment* Seg3[numPos+1];
	int l=-1;
	int lowlimit;
	int numSeg = SegQ->GetNumSegments();// determine the number of segment types
	for (int i=1; i<=numSeg; i++) {
		SF_Segment* Segm = SegQ->GetSegment(i);
		if (Segm->GetMC()) {// determine for each segment type wether it is a MC segment.
			Segm ->ClearAllPos(); // just to be sure remove oldpos.
		}
	}
	for (int i=1; i<=numSeg; i++) {
		SF_Segment* Segm = SegQ->GetSegment(i);
		if (Segm->GetMC()) {// determine for each segment type wether it is a MC segment.
			l++;//count the number of MC segment types
			if (l == 0)	{
				lowlimit = 1;
			}
			else {lowlimit = numposseg[l-1]+1; //determine the lowest rank number of the node of type Segm
			}
			for (int k=lowlimit ; k<=numposseg[l]; k++) {
				Seg2[k]=Segm;
				double waarde;
				if (MyInput->ValueSet("mon", Segm->GetName(), "Spot_mix") ) { // if spotmix is present we can mix two endtypes
					Text segname = MyInput->GetText("mon",Segm->GetName(),"Spot_mix");
					double waarde = MyInput->GetInt("mon",segname,"Spot_frac",0.0,1.0,1.0);
					SF_Segment* Segm2=NULL;
					for (int i=1; i<=numSeg; i++) {
						if ((SegQ->GetSegment(i))->GetName() == segname) {
							Segm2 = SegQ->GetSegment(i);
							break;
						}
					}

					Seg3[k]=Segm2;
					Seg3[k]->UpdatePos(r_new[k],waarde);
				}
				waarde = MyInput->GetInt("mon",Segm->GetName(),"Spot_frac",0.0,1.0,1.0);
				Seg2[k]->UpdatePos(r_new[k],waarde);

			}
				//If so update the the positions and mask file.


		}
	}
	return CheckOverlap(SpotSize);
}

bool
SF_SystemMC::Off_Lattice_MovePos(int n_move,SF_Segment** Seg, int neighbor){
	for (int i=1; i <= n_move; i++){
		double delta =0.0, delta_x=0.0,delta_y=0.0, delta_z=0.0;
		int rndpart,neighpart=0;
		double distance,closest=pow(neighbor+1,2.0);
		rndpart = random_int(1,numPos);
		if (neighbor>0) { //alleen clustermoves when neighbor distance is defined.
			while (delta == 0){
				delta = random(-1*d_max_pos,d_max_pos);
			}
			if (delta>0) { //50% van de keren clustermoves implementeren.
				for (int j=1; j<=numPos; j++) {
					if (j != rndpart) {
						distance = (x_new[j]-x_new[rndpart])*(x_new[j]-x_new[rndpart])+(y_new[j]-y_new[rndpart])*(y_new[j]-y_new[rndpart])+(z_new[j]-z_new[rndpart])*(z_new[j]-z_new[rndpart]);
						if (distance < closest) {closest = distance; neighpart=j;}
					}
				}
				if (closest > neighbor*neighbor) neighpart =0; //als closest particle niet dichterbij zit dan de neighbordistance dan is er geen sprake van een cluster.
			}

		}
		delta=0;
		int direction = random_int(1,3);
		while (delta == 0){
			delta = random(-1*d_max_pos,d_max_pos);
		}
		if (direction == 1 && Seg[rndpart]->MoveAllowed(direction)) delta_x= delta;
		if (direction == 2 && Seg[rndpart]->MoveAllowed(direction)) delta_y= delta;
		if (direction == 3 && Seg[rndpart]->MoveAllowed(direction)) delta_z= delta;
		rx_new[rndpart]=rx_new[rndpart]+delta_x;
		ry_new[rndpart]=ry_new[rndpart]+delta_y;
		rz_new[rndpart]=rz_new[rndpart]+delta_z;


		if (rx_new[rndpart]<1) rx_new[rndpart]=rx_new[rndpart]+n_layers_x;
		if (ry_new[rndpart]<1) ry_new[rndpart]=ry_new[rndpart]+n_layers_y;
		if (rz_new[rndpart]<1) rz_new[rndpart]=rz_new[rndpart]+n_layers_z;
		if (rx_new[rndpart]>n_layers_x+1) rx_new[rndpart]=rx_new[rndpart]-n_layers_x;
		if (ry_new[rndpart]>n_layers_y+1) ry_new[rndpart]=ry_new[rndpart]-n_layers_y;
		if (rz_new[rndpart]>n_layers_z+1) rz_new[rndpart]=rz_new[rndpart]-n_layers_z;
		x_new[rndpart]=int(rx_new[rndpart]); y_new[rndpart]=int(ry_new[rndpart]); z_new[rndpart]=int(rz_new[rndpart]);
		r_new[rndpart]=jx*x_new[rndpart]+jy*y_new[rndpart]+z_new[rndpart]+1;


		if (neighpart > 0) { //als neighbor particle has been identified, verzet het deeltje dan op dezelfde manier als rndpart.
			rx_new[neighpart]=rx_new[neighpart]+delta_x;
			ry_new[neighpart]=ry_new[neighpart]+delta_y;
			rz_new[neighpart]=rz_new[neighpart]+delta_z;
			if (rx_new[neighpart]<1) rx_new[neighpart]=rx_new[neighpart]+n_layers_x;
			if (ry_new[neighpart]<1) ry_new[neighpart]=ry_new[neighpart]+n_layers_y;
			if (rz_new[neighpart]<1) rz_new[neighpart]=rz_new[neighpart]+n_layers_z;
			if (rx_new[neighpart]>n_layers_x) rx_new[neighpart]=rx_new[neighpart]-n_layers_x;
			if (ry_new[neighpart]>n_layers_y) ry_new[neighpart]=ry_new[neighpart]-n_layers_y;
			if (rz_new[neighpart]>n_layers_z) rz_new[neighpart]=rz_new[neighpart]-n_layers_z;
			x_new[neighpart]=int(rx_new[neighpart]); y_new[neighpart]=int(ry_new[neighpart]); z_new[neighpart]=int(rz_new[neighpart]);
			r_new[neighpart]=jx*x_new[neighpart]+jy*y_new[neighpart]+z_new[neighpart]+1;

		}
		if (!free[r_new[rndpart]]) {return false;}
		if (neighpart > 0) { if (!free[r_new[neighpart]]) {return false;}}
	}
	return CheckOverlap(SpotSize);
}


void
SF_SystemMC::Go(Output* Out,Lattice* Lat) {
	bool loopcount = MyInput->GetBoolean("sys",name,"loop_count",false);
	numPos=0;
	numMCseg =0;
	int numSeg = SegQ->GetNumSegments();
	//int numStates = SegQ->GetNumStates();
	int dim = Lat->GetNumGradients();

	if (dim != 3) {Message(fatal,"Sorry: for MC we need 3 gradients. ");}
	n_layers_x = Lat->GetNumLayers(1)-2;
	n_layers_y = Lat->GetNumLayers(2)-2;
	n_layers_z = Lat->GetNumLayers(3)-2;
	jx = (n_layers_y+2)*(n_layers_z+2);
	jy = n_layers_z+2;
	jz = 1;
	int M = (n_layers_x+2)*(n_layers_y+2)*(n_layers_z+2);
	free = new bool[M+1];
	if (Cluster) neighbor= MyInput->GetInt("sys",name,"cluster_move",1,n_layers_z,3);

	rho_0.Dim(1,MolQ->GetNumMolecules());
	for (int i=1; i<=MolQ->GetNumMolecules(); i++) {
		rho_0[i].Dim(1,n_layers_z);
	}

	phi_0.Dim(1,MolQ->GetNumMolecules());
	for (int i=1; i<=MolQ->GetNumMolecules(); i++) {
		phi_0[i].Dim(1,M);
	}

	Omega_0.Dim(1,n_layers_z);
	F_0.Dim(1,n_layers_z);
	Fpo_0.Dim(1,n_layers_z);
	if (SegQ->Charged()) PSI_0.Dim(1,n_layers_z);

	LatticeRange* LatRange;
	//SF_State* State;

	for (int z=0; z<=M; z++) free[z]=true;
	for (int i=1; i<=numSeg; i++) {
		SF_Segment* Seg = SegQ->GetSegment(i);
		if (Seg->GetFreedom() == frozen && !Seg->GetMC()) {
		LatRange = Seg->GetLatRange();
			for (int z=0; z<=M; z++) {
				if (LatRange->InRange(z)) {
					free[z] = false;
				}
			}
		}
	}



	SpotSize=1; // overlap allowed
	int VarSpotSize;
	for (int i=1; i<=numSeg; i++) {
		SF_Segment* Seg = SegQ->GetSegment(i);
		if (Seg->GetMC()) {
			if (Seg->GetSpotType() ==1) { //spottype 1 = sphere, spottype 2 = node
				VarSpotSize=Seg->GetSpotSize();
				if (SpotSize>1 && VarSpotSize != SpotSize) {
					Message(fatal,"Only one SpotSize allowed (temporarily....). ");
				} else SpotSize=VarSpotSize;
			}
			numMCseg++;
			LatRange=Seg->GetLatRange();
			numPos += LatRange->GetNumPos();
			numposseg.push_back (numPos);
		}
	}

	if (numMCseg == 1 && swapfreq > 0.0){Message(fatal,"Only one type of MC node is present while swapping is only usefull when more than one type is present.");}

	if (Solve->e_info || Solve->s_info) {cout << "MC started" << endl;}
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



	int k=0;
	int rr,rx,ry,rz,np,accept;
	for (int i=1; i<=numSeg; i++) {
		SF_Segment* Segm = SegQ->GetSegment(i);
		if (Segm->GetMC()) {
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
		if (Segm->GetMC()) {
			Segm ->ClearAllPos();
		}
	}

	if (!lat_moves) {
		for (int kk=1; kk<=numPos; kk++) {
			rx_new[kk]=1.0*x_new[kk];
			ry_new[kk]=1.0*y_new[kk];
			rz_new[kk]=1.0*z_new[kk];
			int l=0;
			for (int i=-1; i<2; i++)
			for (int j=-1; j<2; j++)
			for (int k=-1; k<2; k++) {
				if (free[jx*(x_new[kk]-i) +jy*(y_new[kk]-j)+jz*(z_new[kk]-k)+1]){
					submask[l]=1.0;
				}
				else {
					submask[l]=0.0;
				}
				l++;
			}
			Seg[kk]->UpdatePos(rx_new[kk],ry_new[kk],rz_new[kk],submask);
		}
	} else 	for (int k=1; k<=numPos; k++) {Seg[k]->UpdatePos(r_new[k]); }

	if (!CheckOverlap(SpotSize)) Message(fatal,"Initial positions overlap or there is not enough space in the pinned ranges");

	Solve->SetInitialGuess(Solve);

	Solve->Iterate();

	for (int p=1; p<=numPos; p++) {r_last[p]=r_new[p];} //hoop van zegen.
	//Vector profile;
	//profile = MolQ->GetFreeEnergyProfile();
	if (MC_F) {F_system_N=MolQ->GetFreeEnergy();} else {F_system_N=GetFreeEnergyPo(Lat);}
	mcs=0;
	cout << "MCS "<< mcs << " F = " << F_system_N << endl;
	for (int k=1; k<=numPos; k++) cout << "("<< x_new[k] << ","<< y_new[k] <<","<< z_new[k] <<")";
	cout << endl;

	SetRho(Lat,rho_0); SetF(Lat,F_0); SetFpo(Lat,Fpo_0); SetOmega(Lat,Omega_0);
	DataMgr(1);
	if (SegQ->Charged()) SetPsi(Lat,PSI_0);
	if (Solve->ErrorOccurred()) {F_system_N = 10000000.0;}

	accept=0;
	double desiredaccept = 0.234;//(this defines the desired acceptance ratio for the monte carlo steps)
	int accepttotal=0;
	double MCSdamping = pow(0.5,0.5);//this sets the amount of damping for the adjustment of the size of the montecarlosteps. The higher the value the stronger the damping and thus the smaller the fluctuations. It however also takes longer to reach the optimal value.
	double n_moves_new = (double)n_moves;
	bool solved = false;
	bool nooverlap = true;
	while ( mcs < MCS) { // hoop !Solve->ErrorOccurred();
		bool write_output = true;

		mcs++;
		F_system_O = F_system_N;
		for (int k=1; k<=numPos; k++) {
			x_old[k]=x_new[k];
			y_old[k]=y_new[k];
			z_old[k]=z_new[k];
			r_old[k]=r_new[k];
			//cout << "1: old = new" << endl;
		}


		if (!lat_moves) {
			for (int k=1; k<=numPos; k++) {
				rx_old[k]=rx_new[k];
				ry_old[k]=ry_new[k];
				rz_old[k]=rz_new[k];
				//cout << "2: old = new" << endl;
			}//
		}

		if (!MovePos(n_moves,Seg,neighbor)){// If a particle has been placed on top of another particle or in a frozen range the particles are set back to their original positions
			nooverlap = false;
			for (int k=1; k<=numPos; k++) {
				x_new[k]=x_old[k];
				y_new[k]=y_old[k];
				z_new[k]=z_old[k];
				r_new[k]=r_old[k];
				//cout << "3: new = old" << endl;
				solved = true;
			}

			if (!lat_moves) {
				for (int k=1; k<=numPos; k++) {
					rx_new[k]=rx_old[k];
					ry_new[k]=ry_old[k];
					rz_new[k]=rz_old[k];
					//cout << "4: new = old" << endl;
				}
			}
		}
		else {
			nooverlap = true;
			Solve->MoveInitialGuess(r_new,r_last,numPos);
			Solve->SetInitialGuess(GetSolve());
			for (int i=1; i<=numSeg; i++) {
				SF_Segment* Segm = SegQ->GetSegment(i);
				if (Segm->GetMC()) {
					Segm ->ClearAllPos();
				}
			}

			if (!lat_moves) {
				for (int kk=1; kk<=numPos; kk++) {
					int l=0;
					for (int i=-1; i<2; i++)
					for (int j=-1; j<2; j++)
					for (int k=-1; k<2; k++) {if (free[jx*(x_new[kk]-i) +jy*(y_new[kk]-j)+jz*(z_new[kk]-k)+1]) submask[l]=1; else submask[l]=0; l++; }
					Seg[kk]->UpdatePos(rx_new[kk],ry_new[kk],rz_new[kk],submask);
				}
			} else 	{for (int k=1; k<=numPos; k++) Seg[k]->UpdatePos(r_new[k]); }

			Solve->SetInitialGuess(Solve);
			Solve->Iterate();
			solved = !Solve->ErrorOccurred();
			if (solved) {
				for (int k=1; k<=numPos; k++) {r_last[k]=r_new[k];}
				if (MC_F) {F_system_N=MolQ->GetFreeEnergy();
				} else {
					F_system_N=GetFreeEnergyPo(Lat);
				}
				if (F_system_N > F_system_O) {
					if (random_d() > exp((F_system_O-F_system_N)*temperature/mc_temperature)) {
						for (int k=1; k<= numPos; k++){
							x_new[k]=x_old[k];
							y_new[k]=y_old[k];
							z_new[k]=z_old[k];
							//Fpos[k]=Fpos0[k];
							r_new[k]=r_old[k];
							//cout << "5: new = old" << endl;
						}
						if (!lat_moves) {
							for (int k=1; k<= numPos; k++){
								rx_new[k]=rx_old[k];
								ry_new[k]=ry_old[k];
								rz_new[k]=rz_old[k];
								//cout << "6: new = old" << endl;
							}
						}
						F_system_N= F_system_O;

						for (int i=1; i<=numSeg; i++) {
							SF_Segment* Segm = SegQ->GetSegment(i);
							if (Segm->GetMC()) {
								Segm ->ClearAllPos();
							}
						}
						if (!lat_moves) {
							for (int kk=1; kk<=numPos; kk++) {
								int l=0;
								for (int i=-1; i<2; i++)
								for (int j=-1; j<2; j++)
								for (int k=-1; k<2; k++) {if (free[jx*(x_new[kk]-i) +jy*(y_new[kk]-j)+jz*(z_new[kk]-k)+1]) submask[l]=1; else submask[l]=0; l++; }
								Seg[kk]->UpdatePos(rx_new[kk],ry_new[kk],rz_new[kk],submask);
							}
						} else
						for (int k=1; k<=numPos; k++) {Seg[k]->UpdatePos(r_new[k]); }
					} else {accept = 1;}
				} else {accept = 1;}
			} else {
				write_output = false;
				mcs--;
				cout << "Failed to calculate the free energy." << endl;
				cout << "A new sety of particle positions will be generated." << endl;
				cout << "The particle positions for which the free energy could not be calculated will " << endl;
				for (int k=1; k<= numPos; k++){
					cout << x_new[k] << y_new[k] << z_new[k] << endl;
					x_new[k]=x_old[k];
					y_new[k]=y_old[k];
					z_new[k]=z_old[k];
					r_new[k]=r_old[k];

				}
				if (!lat_moves) {
					for (int k=1; k<= numPos; k++){
						rx_new[k]=rx_old[k];
						ry_new[k]=ry_old[k];
						rz_new[k]=rz_old[k];
						//cout << "8: new = old" << endl;
					}
				}

				F_system_N= F_system_O;
			}
		}
			//Here the number of steps (n_moves) is adjusted so the acceptance ratio for the MC steps is close to the desired value
			// danger if the acceptance rate is not dependant on the n_moves the value can drift to an extreme value.
		if ((n_equi_steps >= 0) && (mcs > n_equi_steps)) opt_mc = false;
		if (opt_mc && solved && nooverlap){
			n_moves_new *= 1+(2*accept-2*desiredaccept)/(pow(mcs,MCSdamping)+1);
			n_moves = (int)n_moves_new;
			if (n_moves < 1){
				n_moves = 1 ;
				n_moves_new = 1;
			}
			if (n_moves > numPos*3*10){
				n_moves = numPos*30 ;
				n_moves_new = 30.0*numPos;
			}
		}
		accepttotal += accept;
		if (write_output) {
			//double FN;
			if (accept==1) {
				SetRho(Lat,rho_0); SetF(Lat,F_0);
				SetFpo(Lat,Fpo_0); SetOmega(Lat,Omega_0);
				if (SegQ->Charged()) SetPsi(Lat,PSI_0);
			}
			DataMgr(accept);
			GetOutput(Out,Lat);
			if (loopcount){
				Loopcounting(Out,Lat,Seg,r_new,numPos,accept);
			}
			if (mcs==1) Out->WriteOutput(); else Out->WriteMCOutput(mcs);

			Out->Clear();


			//profile = MolQ->GetFreeEnergyProfile();

			//if (MC_F) {FN=MolQ->GetFreeEnergy();
			//} else {
			//	FN=GetFreeEnergyPo(Lat);
			//}

			cout << "MCS "<< mcs << setprecision(7)<<" F = " << F_system_N << setprecision(0)<< " accept " << accept << " n_moves " << n_moves << " accepttotal " << accepttotal <<endl;
			for (int k=1; k<=numPos; k++) {
					cout << "("<< x_new[k] << ","<< y_new[k] <<","<< z_new[k] <<")";
					if (!lat_moves){
						cout << "("<< rx_new[k] << ","<< ry_new[k] <<","<< rz_new[k] <<")";}
			} cout << endl;
		accept = 0;

		}
	}
}

void
SF_SystemMC::GetExtraOutput(Output* Out, Lattice* Lat) const{
	Array<Text> choices(1,1);
	choices[1] = "ParticleInBrush";

	int choice = MyInput->GetChoice("sys",name,"extra_output",choices);
	switch(choice) {
		case 1 :
			GetOutputParticleInBrush(Out,Lat);
			break;
		default :
			Message(fatal,"Error in SF_System::GetExtraOutput");
			break;
	}
}

void
SF_SystemMC:: Loopcounting(Output* Out,Lattice* Lat,SF_Segment** Seg,int* r_new,int numPos,int accept){
	int numSeg = SegQ->GetNumSegments();
	if (accept==1){
		if (!MyInput->ValueSet("sys",name,"probe_chain")) Message(fatal,"sys : probe_chain not defined");
		Text chainName = MyInput->GetText("sys",name,"probe_chain");
		if (!MolQ->MoleculeDefined(chainName)) {Message(fatal,"sys : " + name + " : probe_chain : " + chainName + " is a molecule that is not defined");}
		SF_Molecule* Mol = MolQ->GetMolecule(chainName);
		Vector nloops(1,numPos);
		int n_layers_x = Lat->GetNumLayers(1)-2;
		int n_layers_y = Lat->GetNumLayers(2)-2;
		int n_layers_z = Lat->GetNumLayers(3)-2;
		int jx = (n_layers_y+2)*(n_layers_z+2);
		int jy = n_layers_z+2;
		double theta=0;
		double meanloops=0;
		double dloop=0;
		double Length = Mol->GetChainLength();
		for (int k=1; k<=numPos; k++) {
			for (int i=1; i<=numSeg; i++) {
				SF_Segment* Segm = SegQ->GetSegment(i);
				if (Segm->GetMC()) {
					Segm->ClearAllPos();
				}
			}
			Seg[k]->UpdatePos(r_new[k]);
			Solve->ReComputePhi(true);
			Vector phi = Mol->GetPhi(total);
			theta=0;
			for (int x=1; x<=n_layers_x; x++)
			for (int y=1; y<=n_layers_y; y++)
			for (int z=1; z<=n_layers_z; z++) theta+=phi[jx*x+jy*y+z];
			nloops[k]=theta/Length;
			meanloops +=nloops[k];
		}
		meanloops /= numPos;
		for (int k=1; k<=numPos; k++) {
			dloop += pow(nloops[k]-meanloops,2);
		}
		lastvarloop = pow(dloop/numPos,0.5);
		lastaverageloop =meanloops;
	}

	Out->PutReal("sys",name,"<nloop>",lastaverageloop);
	Out->PutReal("sys",name,"d_loop",lastvarloop);
	if (accept==1){
		for (int i=1; i<=numSeg; i++) {
			SF_Segment* Segm = SegQ->GetSegment(i);
			if (Segm->GetMC()) {
				Segm->ClearAllPos();
			}
		}
		for (int k=1; k<=numPos; k++) {Seg[k]->UpdatePos(r_new[k]); }
		Solve->ReComputePhi(false);
	}
}


void
SF_SystemMC::GetOutputParticleInBrush(Output* Out,Lattice* Lat) const{
	Text Tk;
	for (int k=1; k<=numPos; k++) {
		Tk = Blanks(9);
		Tk.Putint(k);
		Tk = Copy(Tk.Strip());
		Tk = Copy(Tk.Frontstrip());
		Out ->PutInt("sys",name,"x["+Tk+"]",x_new[k]);
		Out ->PutInt("sys",name,"y["+Tk+"]",y_new[k]);
		Out ->PutInt("sys",name,"z["+Tk+"]",z_new[k]);
	}
}



void
SF_SystemMC::GetOutput(Output* Out,Lattice* Lat) const {
	Out->PutText("sys",name,"inputfile",MyInput->GetFileName());
	Out->PutText("sys",name,"calculation_type","MC_SCF");
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
	if (MC_F) Out->PutText("sys",name,"free_energy_type","free_energy"); else Out->PutText("sys",name,"free_energy_type","free_energy(po)");
	if (opt_mc) Out->PutText("sys",name,"move_strategy","variable"); else Out->PutText("sys",name,"move_strategy","fixed");
	Out->PutReal("sys",name,"MC_temperature",mc_temperature);
	Out->PutInt("sys",name,"mcs", mcs);
	Out->PutInt("sys",name,"MCS", MCS);
	Out->PutInt("sys",name,"d_max_pos",d_max_pos);
	Out->PutInt("sys",name,"num_of_moves",n_moves);
	if (Cluster) Out->PutInt("sys",name,"cluster_moves",neighbor);
	if (lat_moves) Out->PutText("sys",name,"lattice_moves","true"); else Out->PutText("sys",name,"lattice_moves","false");

	if (SegQ->Charged()) {
		Out->PutProfile("sys",name,"potential",SegQ->GetElectricPotential());
		Out->PutProfile("sys",name,"charge",SegQ->GetCharge());
		Out->PutProfile("sys",name,"epsilon",SegQ->GetAverageEpsilon());
	}

	GetMCOutput(Out,Lat);

	if (MyInput->ValueSet("sys",name,"Probe_molecule")) GetEnd_to_EndOutput(Out,Lat);
	if (MyInput->ValueSet("sys",name,"extra_output")) {	GetExtraOutput(Out,Lat);}

	Lat1st->GetOutput(Out);
	MolQ->GetOutput(Out);
	SegQ->GetOutput(Out);
	ReactionQ->GetOutput(Out);
}

void
SF_SystemMC::GetEnd_to_EndOutput(Output* Out,Lattice* Lat) const{

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
SF_SystemMC::GetMCOutput(Output* Out,Lattice* Lat) const{
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
	Out->PutReal("sys",name,"MC_Free_energy",F_system_N);
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
	double mean_x=0;
	double mean_y=0;
	double mean_z=0;
	int segnumber =0; // counts the number of MC segment types found so far
	int klow;
	for (int i=1; i<=numSeg; i++) { // loop over all segment types
		SF_Segment* Seg = SegQ->GetSegment(i);
		if (Seg->GetMC()) { // check wether segment type is moved by montcarlo
			if (segnumber == 0){  // determine the indexnumber belonging to the last monte carlo position belonging to the previous type
				klow = 0;
			}
			else {
				klow = numposseg[segnumber-1];
			}

			for (int k=klow+1 ; k<= numposseg[segnumber];k++){ // loop over all positions of the present monte carlo type
				dim = Blanks(9);
				dim.Putint(k-klow);
				dim = Copy(dim.Strip());
				dim = Copy(dim.Frontstrip());
				Out->PutInt("mon",Seg->GetName(),"r[" + dim + "]_x",x_new[k]);
				Out->PutInt("mon",Seg->GetName(),"r[" + dim + "]_y",y_new[k]);
				Out->PutInt("mon",Seg->GetName(),"r[" + dim + "]_z",z_new[k]);
			}
			segnumber++;
			LatRange=Seg->GetLatRange();
			int num_seg_pos = LatRange->GetNumPos();
			for (int k=1; k<=num_seg_pos; k++){
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
SF_SystemMC::SetF(Lattice* Lat, Vector  F_0) {
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
SF_SystemMC::SetFpo(Lattice* Lat, Vector Fpo_0) {
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
SF_SystemMC::SetOmega(Lattice* Lat, Vector Omega_0) {
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
SF_SystemMC::SetPsi(Lattice* Lat, Vector PSI_0) {
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
SF_SystemMC::SetRho(Lattice* Lat, Array<Vector>  rho_0) const {
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

void
SF_SystemMC::DataMgr(int accept) {

	int M = (n_layers_x+2)*(n_layers_y+2)*(n_layers_z+2);
	Vector profile;
	if (accept == 1) {
		for (int i=1; i<=MolQ->GetNumMolecules(); i++) {
			SF_Molecule* Mol = MolQ->GetMolecule(i);
			profile = Mol->GetPhi(total);
			for (int z=1; z<=M; z++) phi_0[i][z]=profile[z];
		}

	} else {
		for (int i=1; i<=MolQ->GetNumMolecules(); i++) {
			SF_Molecule* Mol = MolQ->GetMolecule(i);
			profile = Mol->GetPhi(total);
			for (int z=1; z<=M; z++) profile[z]=phi_0[i][z];
		}

	}

}
