#include "SF_Comb1stO.h"
#include "tools.h"

SF_Comb1stO::SF_Comb1stO(Text name_,
						   SF_SegmentList* SegQ_,
						   Lat1stO* Lat_,
						   Input* MyInput_)
: SF_Molecule(name_, SegQ_, Lat_, MyInput_)
{
	Lat = Lat_;
	//if (MayerSaupeSet()) {
	//	Message(fatal,MyInput,"In SF_Comb1stO: Mayer_Saupe not implemented for first order propagators");
	//}
	force_set = (!GetForce()==0);
	if (force_set) {
		Message(fatal,MyInput,
				"Force ensemble not implemented in SF_Comb1stO");
	}
	if (GetGPU()) {
		Message(literal,MyInput,"For comb-polymer 1st order the GPU card is not yet activated. Contact Frans Leermakers. Going classical instead");
	}

	if (Chain->GetMoleculeType() != comb) {
		Message(fatal,MyInput,
		"Programming error, trying to create a comb-chain for molecule '" +
		name + "'");
	}

	numDiffSegments = Chain->GetNumDiffSegments();
	if (freedom == secondGeneration) {
		CreatePhi(constrained);
		CreatePhi(unconstrained);
	} else {
		CreatePhi(total);
	}

	SetPhiBulk(phiBulk);

	numArms=Chain->GetNumArms(); NG=numArms;
	lin=1; den=2; a_den=3; bra=4;
	NN=0;
	side_type = Chain->GetSideType();

	symmetric = Chain->Symmetric();
	//if (symmetric) cout << "symmetric" << endl;
	//else cout << "not symmetric" << endl;

	SegBlock1 = Chain->GetSegmentBlockComb(1); N1=SegBlock1->GetMaxLength();
	SegBlock2 = Chain->GetSegmentBlockComb(2); N2=SegBlock2->GetMaxLength();
	if (side_type == lin) { SegBlock3 = Chain->GetSegmentBlockComb(3); N3=SegBlock3->GetMaxLength();}
	SegBlock4 = Chain->GetSegmentBlockComb(4); N4=SegBlock4->GetMaxLength();

	if (side_type == den) {
		N3 = 0;
		SF_SegmentBlock* SegBlock;
		numGenerations=Chain->GetNumGenerations();
		for (int i=1; i<=numGenerations; i++) {
			SegBlock = Chain->GetSegmentBlock(i);
			N3 += SegBlock->GetMaxLength();
		}
	}

	shortest=1; longest=2;
	if (side_type == a_den) {
		numGenerations=Chain->GetNumGenerations();
		cumlengths_short = new int [numGenerations+1];cumlengths_short[0]=0;
		cumlengths_long = new int [numGenerations+1]; cumlengths_long[0]=0;

		for (int i=numGenerations; i>=1; i--) {
			SegBlockA1 = Chain->GetSegmentBlock(i,shortest);
			SegBlockA2 = Chain->GetSegmentBlock(i,longest);
			cumlengths_short[numGenerations-i+1] = cumlengths_short[numGenerations-i]+SegBlockA1->GetMaxLength();
			cumlengths_long[numGenerations-i+1]  = cumlengths_long[numGenerations-i]+ SegBlockA2->GetMaxLength();
		}
		N3 = cumlengths_short[numGenerations];
		NN = cumlengths_long[numGenerations];
	}

	n_links=0;
	if (side_type ==bra) {
		SF_LinkNodeTree *tree=Chain->GetLinkNodeTree();
		n_links=tree->n_links;
		n_nodes=n_links+1;
		Links=tree->LinkList;
		Nodes=tree->NodeList;
		N3=0;
	}

	M=Lat->GetTotalNumLayers();
	NDS=numDiffSegments;
	phi_long.Dim(1,M); //backbone
	phi_short.Dim(1,M); //sides
}

SF_Comb1stO::~SF_Comb1stO() {
	if (freedom == secondGeneration) {
		DeletePhi(constrained);
		DeletePhi(unconstrained);
	} else {
		DeletePhi(total);
	}
	if (side_type == a_den) {
		delete [] cumlengths_short;
		delete [] cumlengths_long;
	}
}
void
SF_Comb1stO::ReComputePhi() {
	Message(fatal,"ReComputePhi not implemented in SF_Comb1stO");
}
void
SF_Comb1stO::ComputePhi() {
	if (Lat->GetNamesBulkBoundaries().Upperbound() > 0) {
		CreatePhi(bulk);
		MatrixBulk(bulk);
	}
	int i;
	for (i=1; i<=numDiffSegments; i++) {
		Seg = Chain->GetDiffSegment(i);
		Vector G = Seg->GetSWF();
		Lat->SetBoundaries(G);
	}
	if (freedom == thirdGeneration) {
		Message(implementation, "thirdGeneration not implemented for a dendrimer in SF_Comb1stO");
	}
	if (freedom != secondGeneration) {
		if (!saveMemory) {
			if (symmetric)
				SymMatrix1(total);
			else Matrix1(total);
		} else {
			Message(implementation, "cannot save memory for a dendrimer in SF_Comb1stO");
		}
	} else {
		if (!saveMemory) {
			Matrix2ndGen(LatRange,constrained,unconstrained);
		} else {
			Message(implementation, "cannot save memory for a dendrimer in SF_Comb1stO");
		}
	}
	if (Chain->SomeSegmentsPinned()) SetPhiBulk(0);
	if (Chain->SomeSegmentsGrafted()) SetPhiBulk(0);
	if (Lat->GetNamesBulkBoundaries().Upperbound() > 0) {
		if (freedom == thirdGeneration || freedom == secondGeneration) {
			CopyBulkBoundaries(unconstrained);
		} else {
			CopyBulkBoundaries(total);
		}
		DeletePhi(bulk);
	}
	int z;
	if (freedom == secondGeneration && numDiffSegments > 1) {
		Vector phi1Tot = GetPhi(constrained);;
		Vector phi2Tot = GetPhi(unconstrained);
		for (z=1; z<=M; z++) {
			phi1Tot[z] = 0;
			phi2Tot[z] = 0;
		}
		for (i=1; i<=numDiffSegments; i++) {
			Seg = Chain->GetDiffSegment(i);
			Vector phi1 = Seg->GetPhi(constrained);
			Vector phi2 = Seg->GetPhi(unconstrained);
			for (z=1; z<=M; z++) {
				phi1Tot[z] += phi1[z];
				phi2Tot[z] += phi2[z];
			}
		}
	} else if (numDiffSegments > 1) {
		Vector phiTot = GetPhi(total);
		for (z=1; z<=M; z++) {
			phiTot[z] = 0;
		}
		for (i=1; i<=numDiffSegments; i++) {
			Seg = Chain->GetDiffSegment(i);
			Vector phi = Seg->GetPhi(total);
			for (z=1; z<=M; z++) {
				phiTot[z] += phi[z];
			}
		}
	}
	if (freedom == secondGeneration) {
		if (muFixed) {
			Vector phi1Tot = GetPhi(constrained);;
			Lat->SubtractBoundaries(phi1Tot);
			Lat->MultiplyWithLatticeSites(phi1Tot);
			theta = 0;
			for (int z=1; z<=M; z++) {
				theta += phi1Tot[z];
			}
			Lat->DivideByLatticeSites(phi1Tot);
			Lat->RestoreBoundaries(phi1Tot);
		}
	}
}
void
SF_Comb1stO::GetLayerAnalysis(Output*) {
	if (LatRangeStartLoops == NULL || LatRangeTrains == NULL) {
		return;
	}
	Message(implementation,"SF_Comb1stO::GetLayerAnalysis");
}
Vector
SF_Comb1stO::GetLong() {
	return phi_long;
}
Vector
SF_Comb1stO::GetShort() {
	return phi_short;
}

Vector
SF_Comb1stO::GetBondOrientation( const DensityPart DensPar) {
	Message(fatal,"GetBondOrientation not implemented in SF_Comb1stO");
	Vector x;
	return x; 
};

MoleculeType
SF_Comb1stO::GetMoleculeType() const {
	return Chain->GetMoleculeType();
}

void
SF_Comb1stO::MatrixLin(Matrix gi) {
	Vector G=SegBlock3->GetSegment(N3)->GetSWF();
	for (int z=1; z<=M; z++) gi[z][1]=G[z];
	for (int s=N3-1; s>=1; s--) {
		G=SegBlock3->GetSegment(s)->GetSWF();
		Lat->PropagateG(gi,G,N3-s+1);
	}
	G=SegBlock2->GetSegment(N2)->GetSWF();
	Lat->PropagateG(gi,G,N3+1);
}
void
SF_Comb1stO::MatrixLin(Matrix gi,Vector Gi_inv,const DensityPart DensPart) {
	Vector phi,G;
	for (int s=1; s<=N3; s++){
		Seg =SegBlock3->GetSegment(s);
		G=Seg->GetSWF();
		phi = Seg->GetPhi(DensPart);
		Lat->PropagateG(Gi_inv,G);
		for (int z=1; z<=M; z++) phi[z] += Gi_inv[z]*gi[z][N3-s+1];
	}
}

void
SF_Comb1stO::MatrixDen(Matrix Gi) {
	Vector Gold,G;
	Vector Gi_inv(1,M);
	Vector Gi_add;

	SF_SegmentBlock* SegBlock;
	int s=0;
	int repeat;
	for (int i=numGenerations; i>=1; i--) {
		SegBlock = Chain->GetSegmentBlock(i);
		// descend through different homo-blocks within each generation
		for (int j=SegBlock->GetMaxLength(); j>=1; j--) {
			Gold = G;
			G = SegBlock->GetSegment(j)->GetSWF();
			s++;
			if (s == 1) {
				for (int z=1; z<=M; z++) Gi[z][1] = G[z];
			} else if (j==SegBlock->GetMaxLength() && Chain->GetNumRepeat(i+1) > 1) {
				repeat = Chain->GetNumRepeat(i+1);
				Vector Gi_temp(1,M);
				for (int z=1; z<=M; z++) Gi_temp[z] = Gi_inv[z] = Gi[z][s-1]; // Gi_inv used locally here
				for (int k=1; k<=repeat-1; k++) {
					Gi_add.Dim(1,M);
					Lat->MakeSafe(Gi_add);
					Lat->ConnectG(Gi_inv,Gi_temp,Gi_add);
					Lat->CorrectDoubleCountG(Gi_add,Gold);
					for (int z=1; z<=M; z++) Gi_temp[z] = Gi_add[z];
				}
				Lat->PropagateG(Gi_add,G);
				for (int z=1; z<=M; z++) Gi[z][s] = Gi_add[z];
			} else {
				Lat->PropagateG(Gi,G,s);
			}
		}
	}
	Gold = G;
	G = SegBlock2->GetSegment(N2)->GetSWF();
	s++;
	repeat = Chain->GetNumRepeat(1);
	if (repeat > 1) {
		Vector Gi_temp(1,M);
		for (int z=1; z<=M; z++) Gi_temp[z] = Gi_inv[z] = Gi[z][s-1]; // Gi_inv used locally here
		for (int k=1; k<=repeat-1; k++) {
			Gi_add.Dim(1,M);
			Lat->MakeSafe(Gi_add);
			Lat->ConnectG(Gi_inv,Gi_temp,Gi_add);
			Lat->CorrectDoubleCountG(Gi_add,Gold);
			for (int z=1; z<=M; z++) Gi_temp[z] = Gi_add[z];
		}
		Lat->PropagateG(Gi_add,G);
		for (int z=1; z<=M; z++) Gi[z][s] = Gi_add[z];
	} else {
		Lat->PropagateG(Gi,G,s);
	}
}

void
SF_Comb1stO::MatrixDen(Matrix Gi,Vector Gi_inv,const DensityPart DensPart) {
	int s=N3;
	Vector phi,G;
	SF_SegmentBlock* SegBlock;
	int repeat;
	double NormRepeat = 1;
	for (int i=1; i<=numGenerations; i++) {
		SegBlock = Chain->GetSegmentBlock(i);
		repeat = Chain->GetNumRepeat(i);
		for (int j=1; j<=SegBlock->GetMaxLength(); j++) {
			Seg = SegBlock->GetSegment(j);
			G = Seg->GetSWF();
			Lat->PropagateG(Gi_inv,G);
			phi = Seg->GetPhi(DensPart);
			if (j == 1) { // branchpoint
				for (int k=1; k<=repeat-1; k++) {
					Vector Gi_add(1,M);
					Lat->MakeSafe(Gi_add);
					Lat->ConnectG(Gi_inv,Gi,s,Gi_add);
					Lat->CorrectDoubleCountG(Gi_add,G);
					for (int z=1; z<=M; z++) Gi_inv[z] = Gi_add[z];
				}
				Lat->NormPhiFree(Gi_inv,NormRepeat);
				Lat->ConnectG(Gi_inv,Gi,s,phi);
				Lat->NormPhiFree(Gi_inv,1.0/NormRepeat);
				NormRepeat *= repeat;
			} else {
				Lat->NormPhiFree(Gi_inv,NormRepeat);
				Lat->ConnectG(Gi_inv,Gi,s,phi);
				Lat->NormPhiFree(Gi_inv,1.0/NormRepeat);
			}
			s--;
		}
	}
}

void
SF_Comb1stO::MatrixADen(Matrix gi,Matrix Gi) {
	int length_short=0, length_long=0;
	Vector G;
	Vector Gi_add;
		Gi_add.Dim(1,M);
		Lat->MakeSafe(Gi_add);
	Vector Gold;
	Vector Gi_inv(1,M);
	Vector Gi_temp(1,M);
	int ss=0;
	int s=0;
	for (int i=numGenerations; i>=1; i--) {
		SegBlockA2 = Chain->GetSegmentBlock(i,longest); length_long=SegBlockA2->GetMaxLength();
		SegBlockA1 = Chain->GetSegmentBlock(i,shortest); length_short=SegBlockA1->GetMaxLength();

		for (int j=length_long; j>=1; j--) {
			Gold = G;
			G = SegBlockA2->GetSegment(j)->GetSWF();
			s++;
			if (s == 1) {
				for (int z=1; z<=M; z++) Gi[z][1] = G[z];
			} else if (j==length_long) {

				for (int z=1; z<=M; z++) {
					Gi_inv[z] = Gi[z][s-1]; // Gi_inv used locally here
					Gi_temp[z] = gi[z][ss]; //ss van vorige loop
				}
				for (int z=1; z<=M; z++) Gi_add[z]=0;
				Lat->ConnectG(Gi_inv,Gi_temp,Gi_add);
				Lat->CorrectDoubleCountG(Gi_add,Gold);
				for (int z=1; z<=M; z++) {Gi_temp[z] = Gi_add[z];} //used in the other loop.

				Lat->PropagateG(Gi_add,G);
				for (int z=1; z<=M; z++) {Gi[z][s] = Gi_add[z];	}
			} else {
				Lat->PropagateG(Gi,G,s);
			}
		}
		for (int k=length_short; k>=1; k--) {
			G = SegBlockA1->GetSegment(k)->GetSWF();
			ss++;
			if (ss ==1) {for (int z=1; z<=M; z++) gi[z][1] = G[z];}
			else if (k == length_short) {
				Lat->PropagateG(Gi_temp,G);
				for (int z=1; z<=M; z++) {gi[z][ss] = Gi_temp[z];}
			} else {
				Lat->PropagateG(gi,G,ss);
			}
		}
	}
	Gold=G;
	G = SegBlock2->GetSegment(N2)->GetSWF();
	for (int z=1; z<=M; z++) if (Gold[z]>0) Gi_add[z]=Gi[z][NN]*gi[z][N3]/Gold[z];
	Lat->PropagateG(Gi_add,G);
	for (int z=1; z<=M; z++) {gi[z][N3+1] = Gi_add[z];}
}

void
SF_Comb1stO::ComputePhi(int i, Matrix Gi, Matrix gi, Vector Gb, const DensityPart DensPart) {
	if (i>numGenerations) return;
	int length_short,length_long;
	SF_SegmentBlock* SegBlock;
	Vector phi,G;
	Vector Gi_inv(1,M);
	Vector Gi_temp(1,M); for (int z=1; z<=M; z++) Gi_temp[z]=Gb[z];
	SegBlockA2 = Chain->GetSegmentBlock(i,longest); length_long=SegBlockA2->GetMaxLength();
	SegBlockA1 = Chain->GetSegmentBlock(i,shortest); length_short=SegBlockA1->GetMaxLength();
	for (int s=1; s<=length_long; s++) {
		Seg = SegBlockA2->GetSegment(s);
		phi = Seg->GetPhi(DensPart);
		G = Seg->GetSWF();
		if (s==1) {
			for (int z=1; z<=M; z++) Gi_inv[z]=Gi_temp[z]*gi[z][cumlengths_short[numGenerations-i+1]]; //omdat Gb al dedeeld is door gz hoeft dit niet meer.
		} else {
			Lat->PropagateG(Gi_inv,G);
			for (int z=1; z<=M; z++)
				phi[z]+=Gi_inv[z]*Gi[z][cumlengths_long[numGenerations-i+1]-s+1];
		}
	}
	if (i<numGenerations) {
		SegBlock = Chain->GetSegmentBlock(i+1,longest);
		Seg = SegBlock->GetSegment(1);
		G = Seg->GetSWF();
		phi = Seg->GetPhi(DensPart);
		Lat->PropagateG(Gi_inv,G);
		Lat->CorrectDoubleCountG(Gi_inv,G);

		for (int z=1; z<=M; z++) {
			Gb[z]=Gi_inv[z];
			phi[z]+=Gi[z][cumlengths_long[numGenerations-i]]*gi[z][cumlengths_short[numGenerations-i]]*Gi_inv[z];
		}
		ComputePhi(i+1,Gi,gi,Gb, DensPart);
	}
	for (int s=1; s<=length_short; s++) {
		Seg = SegBlockA1->GetSegment(s);
		phi = Seg->GetPhi(DensPart);
		G = Seg->GetSWF();
		if (s==1) {
			for (int z=1; z<=M; z++) Gi_inv[z]=Gi_temp[z]*Gi[z][cumlengths_long[numGenerations-i+1]];
		} else {
			Lat->PropagateG(Gi_inv,G);
			for (int z=1; z<=M; z++)
				phi[z]+=Gi_inv[z]*gi[z][cumlengths_short[numGenerations-i+1]-s+1];
		}
	}
	if (i<numGenerations) {
		SegBlock = Chain->GetSegmentBlock(i+1,shortest);
		Seg = SegBlock->GetSegment(1);
		G = Seg->GetSWF();
		phi = Seg->GetPhi(DensPart);
		Lat->PropagateG(Gi_inv,G);
		Lat->CorrectDoubleCountG(Gi_inv,G);

		for (int z=1; z<=M; z++) {
			Gb[z]=Gi_inv[z];
			phi[z]+=Gi[z][cumlengths_long[numGenerations-i]]*gi[z][cumlengths_short[numGenerations-i]]*Gi_inv[z];
		}
		ComputePhi(i+1,Gi,gi,Gb, DensPart);
	}
}

void SF_Comb1stO::ComputeGzsNode(const int nodenum, const int linkto) {
	int arm; // arm number in the internal numbering of the node (1..n_links)
	int linknum; // link number in the global numbering
	Vector G; // storage space for G(z) of a particular segment
	SF_Node node=Nodes[nodenum]; // to simplify the notation
	SF_Link Linkto; // the link numbered linkto
	if(linkto!=0)Linkto=Links[linkto]; // Links[0] would produce an error
	for(arm=node.n_links();arm>0;arm--) {
		linknum=node.link(arm);
		if(linknum!=linkto) {
			ComputeGzsLink(linknum,nodenum);
		}
	}
	if (linkto!=0)  PropagateGzsNode(nodenum,linkto); // connect arms and proceed to linkto
}


void SF_Comb1stO::ComputeGzsLink(const int linknum, const int nodefrom) {
	int s; // index variable
	int segnum; // internal segment number of a link
	int step; // determines the direction of propagation
	int nodeto; // node number to which the link is connected, other than nodefrom
	SF_Link Link=Links[linknum]; // to simplify the notation - the link we are workin on
	int length=Link.length(); // length of our link
	Vector G; // storage space fo G(z)
	VecArray GzsLink=Gzs[linknum]; // the Gzs matrix of the current link
	if(Link.node(1)==nodefrom) {
		segnum=length; nodeto=Link.node(2); step =-1;
	} else {
		segnum=1; nodeto=Link.node(2); step = 1;
	}
	ComputeGzsNode(nodeto,linknum); // Ask nodeto for Gzs
	for (s=2; s<=length; s++) { // s is an index to make the proper number of cycles
		segnum+=step;  // first increment segnum, the Gzs of the end segment has been provided by the node
		G = Chain->GetSegment(linknum,segnum)->GetSWF();
		Lat->PropagateG(GzsLink[segnum-step], G, GzsLink[segnum]);
	}
}

void SF_Comb1stO::PropagateGzsNode(const int nodenum, const int linkto) {
	int first_arm=1; // boolean to see if
	int arm; // link number in the internal numbering inside the node (1..n_links)
	int link; // link number in the global numbering
	int seg; // segment number of the link
	Vector Gout(1,M); // Gzs of linkto - that where we store the output
	Vector G; // G(z) of the segment
	Vector TempG(1,M); // temporary G Vector
	SF_Node Node=Nodes[nodenum]; // the node we are working on
	int n_links=Node.n_links(); // number of links of the Node
	int segto=Links[linkto].seg(nodenum); // segment number by which linkto is connected to the Node
	G=Chain->GetSegment(linkto,segto)->GetSWF(); // GetSWF() of any of the node segments
	if(n_links==1) Gout.Copy(G); // for the end-node, Gout=G
	else { // go through all the arms but linkto and connect G
		for(arm=1;arm<=n_links;arm++) {
			link=Node.link(arm);
			if(link!=linkto) {
				seg=Links[link].seg(nodenum);
				if(first_arm) {
					first_arm=0;
					TempG.Copy(Gzs[link][seg]); // for the first arm, copy Gzs[link][seg] to Gout
				}
				else {
					Gout.Dim(1,M);
					Lat->MakeSafe(Gout);
					Lat->ConnectG(Gzs[link][seg],TempG,Gout);
					Lat->CorrectDoubleCountG(Gout,G); // divide by G
					TempG.Copy(Gout); // copy the Gout to TempG
				}
			}
		}
	}
	Gzs[linkto][segto].Copy(Gout);
}

void SF_Comb1stO::ComputePhiNode(int nodenum, int linkfrom, const DensityPart DensPart) {
	int arm; // internal number of an arm of the node
	int link; // link number in the absolute numbering
	int seg; // segment number by which the link is connected to the Node
	int linkto; // link number to which we will propagate
	Vector Phi; // Phi(z) vector of the node segment
	Vector PhiNode; // locally stores computed Phi*G of the node
	Vector TempPhi; // temporary storage space for PhiNode
	Vector G; // G(z)
	SF_Node Node=Nodes[nodenum]; // the node we are working on
	int n_links=Node.n_links(); // number of links of the Node
	PhiNode.Dim(1,M); // initialize PhiNode[z] to 0.0
	Lat->MakeSafe(PhiNode);
	TempPhi.Dim(1,M); // initialize PhiNode[z] to 0.0
	Lat->MakeSafe(TempPhi);
	link=Node.link(1); // no matter which link we take, GetSegment gives the same result for the whole node
	seg=Links[link].seg(nodenum);
	Seg=Chain->GetSegment(link,seg); // it does not matter from which link of the node we call GetSWF
	G=Seg->GetSWF(); // it makes no difference to which segment of the node we call GetSWF()
	Phi = Seg->GetPhi(DensPart); // that is where we store the phi in the end
	if (n_links==1) { // it is either a start or end-node
		if(nodenum==1) link=1; // it is a start node
		else link=linkfrom; // if not, Gi_inv of linkfrom is provided by ComputePhiLink
		seg=Links[link].seg(nodenum);
		//if(nodenum==1) {
		//	gi_inv.Copy(G); // initialize the Gi_inv for the 1st node
		//}
		Lat->ConnectG(gi_inv,Gzs[link][seg],PhiNode); // first do that for the linkfrom
	} else { // node has more links
		seg=Links[linkfrom].seg(nodenum);
		TempPhi.Copy(gi_inv); // Copy the Gzs of linkfrom to TempPhi
		for(arm=1;arm<=n_links && linkfrom!=0;arm++) {// through all links but linkfrom
			link=Node.link(arm);
			if(link!=linkfrom) { // for linkfrom we have Gi_inv
				PhiNode.Dim(1,M); // initialize PhiNode to zero
				Lat->MakeSafe(PhiNode);
				seg=Links[link].seg(nodenum);
				Lat->ConnectG(TempPhi,Gzs[link][seg],PhiNode);
				TempPhi.Copy(PhiNode);
			}
		}
		for(arm=2;arm<n_links;arm++) {
			Lat->CorrectDoubleCountG(PhiNode,G);
		}
	}
	// FIXME (?) this will only work with overflow protection switched off
	if(Lat->OverflowProtection()) {
		// this is a workaround
		TempPhi.Dim(1,M); // TempPhi[z]=0.0, in logarithmic form this means G[z]=1.0
		Lat->ConnectG(TempPhi,PhiNode,Phi); // multiply PhiNode by TempPhi and add to Phi
	} else Phi += PhiNode; // simply add the computed PhiNode to the Phi of the segment
	for(arm=n_links;arm>0;arm--) {
		linkto=Node.link(arm);
		if(linkto != linkfrom) { // it is not linkfrom
			seg=Links[linkto].seg(nodenum);
			gi_inv.Copy(PhiNode); // Copy Phi to Gi_inv
			Lat->CorrectDoubleCountG(gi_inv,Gzs[linkto][seg]); // divide by Gzs of segto
			ComputePhiLink(linkto,nodenum,DensPart);
		}
	}
}

void SF_Comb1stO::AddPhi(const Vector Phi1, Vector PhiOut) {
	int z;
	for(z=1;z<=M;z++) {
		if (Phi1[z] > LOGMAXDOUBLE-1 || PhiOut[z] == LOGZERO) PhiOut[z] = Phi1[z];
		else if(PhiOut[z]<LOGMAXDOUBLE) PhiOut[z] = log1p(exp(PhiOut[z])+exp(Phi1[z]));
		// if (PhiOut >=LOGMAXDOUBLE) do not change its value
	}
}

void  SF_Comb1stO::ComputePhiLink(int linknum, int nodefrom, const DensityPart DensPart) {
	int segnum; // segment number in internal numbering of the link
	int step; // direction of propagation
	int nodeto; // number of the other node (other than nodefrom) to which the link is connected
	int s; // index variable
	VecArray GzsLink=Gzs[linknum]; // the Gzs matrix of the link
	SF_Link Link=Links[linknum]; // to simplify the notation - this link
	int length=Link.length(); // length of the Link
	Vector G; // G(z)
	Vector Phi;
	segnum=Link.seg(nodefrom); // segment by which the link is connected to the node
	nodeto=(segnum==1) ? Link.node(2) : Link.node(1); // nodeto is the other node
	step=(segnum==1) ? 1 : -1; // propagate forwards or backwards
	for(s=2;s<length;s++) {
		segnum+=step;
		Seg=Chain->GetSegment(linknum,segnum);
		G=Seg->GetSWF();
		Phi=Seg->GetPhi(DensPart);
		Lat->PropagateG(gi_inv,G,gi_inv);
		Lat->ConnectG(gi_inv,GzsLink[segnum],Phi);
	}
	segnum+=step;
	Seg=Chain->GetSegment(linknum,segnum);
	G=Seg->GetSWF();
	Lat->PropagateG(gi_inv,G,gi_inv);
	ComputePhiNode(nodeto,linknum,DensPart);
}

void
SF_Comb1stO::Matrix1(const DensityPart DensPart) {
	int z,s,ss,i;
	const int lin=1, den=2, a_den=3, bra=4;
	Vector phi,G;
	Matrix Gi(1,M,1,N1+N4+N2*NG);
	Matrix gi(1,M,1,N3+1);
	Matrix gli(1,M,1,NN);
	Vector Gb(1,M);
	Vector Gi_inv(1,M);
	Vector Gi_temp(1,M);
	Gzs.Dim(1,n_links);
	 for (int i=1;i<=n_links;i++) {
		Gzs[i].Dim(1,Links[i].length());
		for (int j=1;j<=Links[i].length();j++) Gzs[i][j].Dim(1,M);
	}

	for (i=1; i<=numDiffSegments; i++) {
		Seg = Chain->GetDiffSegment(i);
		G = Seg->GetSWF();
		Lat->MakeSafe(G);
		phi = Seg->GetPhi(DensPart);
		for (z=1; z<=M; z++) {phi[z] = 0;}
		Lat->MakeSafe(phi);
	}

	Zero(phi_long,M);
	Zero(phi_short,M);
	s = 0;
	ss = 0;

	switch (side_type) {
		case lin:
			MatrixLin(gi);
			break;
		case den:
			MatrixDen(gi);
			break;
		case a_den:
			MatrixADen(gi,gli);
			break;
		case bra:
			ComputeGzsNode(1,0);
			for (z=1; z<=M; z++) Gi_inv[z] = Gzs[1][1][z];
			G=SegBlock2->GetSegment(N2)->GetSWF();
			Lat->PropagateG(Gi_inv,G);
			 for (z=1; z<=M; z++) gi[z][1]=Gi_inv[z];
			break;
		default:
			Message(fatal,"Unknown side chain in comb.");
			break;
	}

	for (int j=1; j<=N1; j++){
		s++;
		G=SegBlock1->GetSegment(j)->GetSWF();
		if (s==1) for (z=1; z<=M; z++) Gi[z][s]=G[z];
		else {
			Lat->PropagateG(Gi,G,s);
		}
	}

	for (i=1; i<=NG; i++) {
		ss=0;
		for (int j=1; j<=N2; j++){
			s++; ss++;
			G=SegBlock2->GetSegment(ss)->GetSWF();
			if (ss==1 && i>1) {
				Lat->PropagateG(Gi_temp,G);
				for (z=1; z<=M; z++) Gi[z][s]=Gi_temp[z];
			}
			else {
				Lat->PropagateG(Gi,G,s);
			}
		}
		for (z=1; z<=M; z++) {if (G[z]>0) Gi_temp[z]=gi[z][N3+1]*Gi[z][s]/G[z];}
	}

	ss=0;
	for (int j=1; j<=N4; j++) { //dinally propagate along the last stretch of the backbone (once).
	 	ss++; s++;
	 	G=SegBlock4->GetSegment(ss)->GetSWF();
	 	if (ss==1) {
			Lat->PropagateG(Gi_temp,G);
			for (z=1; z<=M; z++) Gi[z][s]=Gi_temp[z];
		}
		else Lat->PropagateG(Gi,G,s);
	}

	s=N1+N4+N2*NG+1;
	for (int j=N4; j>=1; j--) {
		s--;
		Seg=SegBlock4->GetSegment(j);
		G=Seg->GetSWF();
		phi=Seg ->GetPhi(DensPart);
		if (j==N4) {
			for (z=1; z<=M; z++) Gi_inv[z]=G[z];
		} else {
			Lat->PropagateG(Gi_inv,G);
		}

		//double phiz,sum=0;
		//for (z=1; z<=M; z++) {phiz = Gi_inv[z]*Gi[z][s]; phi[z]+=phiz; if (G[z]>0) sum +=phiz/G[z]; }
		//cout << "s : " << s << " = " << sum << endl;
		for (z=1; z<=M; z++) phi[z] += Gi_inv[z]*Gi[z][s];
	}

	for (i=NG; i>=1; i--) {
		for (int j=N2; j>=1; j--) {
			s--;
			Seg = SegBlock2->GetSegment(j);
			G=Seg->GetSWF();
			phi = Seg ->GetPhi(DensPart);
			Lat->PropagateG(Gi_inv,G);
			if (j==N2) for (z=1; z<=M; z++) {
				if (G[z]>0) Gi_inv[z] /=G[z];
				Gi[z][s] *=Gi_inv[z]; //prepare for phi side chains
				phi[z]+=Gi[z][s]*gi[z][N3+1]; //compute density branchpoint
				Gi_inv[z] *= gi[z][N3+1]; //prepare for further propagation of backbone
			} else {
				//double phiz,sum=0;
				//for (z=1; z<=M; z++) {phiz = Gi_inv[z]*Gi[z][s]; phi[z]+=phiz; if (G[z]>0) sum +=phiz/G[z]; }
				//cout << "s : " << s << " = " << sum << endl;
				for (z=1; z<=M; z++) phi[z] += Gi_inv[z]*Gi[z][s];
			}
		}
	}

	for (int j=N1; j>=1; j--) {
		s--;
		Seg=SegBlock1->GetSegment(j);
		G=Seg->GetSWF();
		phi=Seg ->GetPhi(DensPart);
		Lat->PropagateG(Gi_inv,G);
		//double phiz,sum=0;
		//for (z=1; z<=M; z++) {phiz = Gi_inv[z]*Gi[z][s]; if (G[z]>0) sum +=phiz/G[z]; }
		//cout << "s : " << s << " = " << sum << endl;
		for (z=1; z<=M; z++) phi[z] += Gi_inv[z]*Gi[z][s];
	}

	lnGN = Lat->ComputeLnGN(Gi_inv);
	for (i=1; i<=numDiffSegments; i++) {
		Seg = Chain->GetDiffSegment(i);
		G = Seg->GetSWF();
		phi = Seg->GetPhi(DensPart);
		for (z=1; z<=M; z++) {
			if (G[z]>0) phi_long[z] +=phi[z]/G[z]; //backbone
		}
	}

	for (i=1; i<=NG; i++) {
		for (z=1; z<=M; z++) Gi_inv[z]=Gi[z][N1+i*N2];
		switch (side_type) {
			case lin:
				MatrixLin(gi,Gi_inv,DensPart);
				break;
			case den:
				MatrixDen(gi,Gi_inv,DensPart);
				break;
			case a_den:
				SF_SegmentBlock* SegBlock;
				SegBlock =  Chain->GetSegmentBlock(1,1);
				Seg = SegBlock->GetSegment(1);
				G = Seg->GetSWF();
				phi = Seg->GetPhi(DensPart);
				Lat->PropagateG(Gi_inv,G);
				for (z=1; z<=M; z++) if (G[z]>0) Gb[z]=Gi_inv[z]/G[z];
				for (z=1; z<=M; z++) phi[z] += gi[z][N3]*gli[z][NN]*Gb[z];
				ComputePhi(1,gli,gi,Gb,DensPart);
				break;
			case bra:
				gi_inv.Dim(1,M);
				for (z=1; z<=M; z++) gi_inv[z]=Gi_inv[z];
				Seg=Chain->GetSegment(1,1);
				G=Seg->GetSWF();
				Lat->PropagateG(gi_inv,G);
				ComputePhiNode(1,0,DensPart);
				break;
			default:
				Message(fatal,"Unknown side chain in comb.");
				break;
		}
	}

	for (i=1; i<=numDiffSegments; i++) {
		Seg = Chain->GetDiffSegment(i);
		G = Seg->GetSWF();
		phi = Seg->GetPhi(DensPart);
		//Lat->CorrectDoubleCountG(phi,G);
		for (int z=1; z<=M; z++) {
			if (G[z] > 0) {
				phi[z] /=G[z];
				phi_short[z] +=phi[z]; //teeth + backbone
			}
		}
	}
	for (int z=1; z<=M; z++) {
		phi_short[z] -= phi_long[z]; //teeth only
	}

	// and do the various renormalizations
	if (freedom == fixedTheta) {
		lnCt = log(theta) - lnGN - log(1.*N);
		phiBulk = theta/exp(lnGN);
		SetPhiBulk(phiBulk);
	} else if (freedom == rangeRestricted) {
		double phiPreNorm = 0;
		for (i=1; i<=numDiffSegments; i++) {
			Seg = Chain->GetDiffSegment(i);
			phi = Seg->GetPhi(DensPart);
			for (z=1; z<=M; z++) {
				if (restrictedRange->InRange(z)) {
					phiPreNorm += phi[z];
				}
			}
		}
		lnCt = log(phiRange/phiPreNorm);
	} else {
		lnCb = log(phiBulk/N);
	}
	for (i=1; i<=numDiffSegments; i++) {
		Seg = Chain->GetDiffSegment(i);
		phi = Seg->GetPhi(DensPart);
		if (freedom == fixedTheta) {
			//Lat->NormPhiRestr(phi,Gi_inv,theta/N);
			double C=theta/(N*exp(lnGN));
			for (int z=1; z<=M; z++) {
				phi[z] *= C;
				phi_long[z] *=C;
				phi_short[z]*=C;
			}
		} else if (freedom == rangeRestricted) {
			//Lat->NormPhiFree(phi,exp(lnCt));
			double C=exp(lnCt);
			for (int z=1; z<=M; z++) {
				phi[z] *= C;
				phi_long[z] *=C;
				phi_short[z]*=C;
			}
		} else {
			//Lat->NormPhiFree(phi,phiBulk/N);
			double C=phiBulk/N;
			for (int z=1; z<=M; z++) {
				phi[z] *= C;
				phi_long[z] *=C;
				phi_short[z]*=C;
			}

		}
		G = Seg->GetSWF();
		Lat->RestoreFromSafe(G);
		Lat->RestoreFromSafe(phi);
	}
	if (freedom == rangeRestricted) {
		theta = ComputeTheta();
		phiBulk = theta/exp(lnGN);
		SetPhiBulk(phiBulk);
	}
}

void
SF_Comb1stO::SymMatrix1(const DensityPart DensPart) {
	int z,s,i;
	int Nb=N1+N4+N2*NG;
	const int lin=1, den=2, a_den=3, bra=4;
	Vector phi,G;
	Matrix Gi(1,M,1,(Nb+1)/2);
	Matrix gi(1,M,1,N3+1);
	Matrix gli(1,M,1,NN);
	Vector Gb(1,M);
	Vector Gi_inv(1,M);
	Vector Gi_temp(1,M);
	Gzs.Dim(1,n_links);
	 for (int i=1;i<=n_links;i++) {
		Gzs[i].Dim(1,Links[i].length());
		for (int j=1;j<=Links[i].length();j++) Gzs[i][j].Dim(1,M);
	}

	for (i=1; i<=numDiffSegments; i++) {
		Seg = Chain->GetDiffSegment(i);
		G = Seg->GetSWF();
		Lat->MakeSafe(G);
		phi = Seg->GetPhi(DensPart);
		for (z=1; z<=M; z++) {phi[z] = 0;}
		Lat->MakeSafe(phi);
	}

	Zero(phi_long,M);
	Zero(phi_short,M);
	s = 0;
	//ss = 0;
	switch (side_type) {
		case lin:
			MatrixLin(gi);
			break;
		case den:
			MatrixDen(gi);
			break;
		case a_den:
			MatrixADen(gi,gli);
			break;
		case bra:
			ComputeGzsNode(1,0);
			for (z=1; z<=M; z++) Gi_inv[z] = Gzs[1][1][z];
			G=SegBlock2->GetSegment(N2)->GetSWF();
			Lat->PropagateG(Gi_inv,G);
			 for (z=1; z<=M; z++) gi[z][1]=Gi_inv[z];
			break;
		default:
			Message(fatal,"Unknown side chain in comb.");
			break;
	}
	for (int j=1; j<=N1; j++){
		s++;
		G=SegBlock1->GetSegment(j)->GetSWF();
		if (s==1) for (z=1; z<=M; z++) Gi[z][s]=G[z];
		else {
			Lat->PropagateG(Gi,G,s);
		}
	}
	for (i=1; i<=(NG+1)/2; i++) {
		for (int j=1; j<=N2; j++){
			s++;
			G=SegBlock2->GetSegment(j)->GetSWF();
			if (j==1 && i>1) {
				Lat->PropagateG(Gi_temp,G);
				for (z=1; z<=M; z++) Gi[z][s]=Gi_temp[z];
			}
			else {
				Lat->PropagateG(Gi,G,s);
			}
		}
		for (z=1; z<=M; z++) {if (G[z]>0) Gi_temp[z]=gi[z][N3+1]*Gi[z][s]/G[z];}
		if (NG%2==0) {
			for (int j=1; j<=N2/2; j++){
				s++;
				G=SegBlock2->GetSegment(j)->GetSWF();
				if (j==1) {
					Lat->PropagateG(Gi_temp,G);
					for (z=1; z<=M; z++) Gi[z][s]=Gi_temp[z];
				}
				else {
					Lat->PropagateG(Gi,G,s);
				}
			}
		}
	}
	if (NG%2==1) for (z=1; z<=M; z++) Gi_inv[z]=Gi[z][s-1];
	else for (z=1; z<=M; z++) Gi_inv[z]=Gi[z][s];

	if (NG%2==0 && Nb%2==1) {
		Seg = SegBlock2->GetSegment(N2/2);
		G = Seg->GetSWF();
		phi = Seg->GetPhi(DensPart);
		//Lat->PropagateG(Gi_inv,G);
		Lat->ConnectG(Gi_inv,Gi_inv,phi);
		Lat->NormPhiFree(phi,0.5);s--;
	}
	if (NG%2==0) {s++;
		for (int j=(N2+3)/2; j<=N2; j++){
			s--;
			G=SegBlock2->GetSegment(j)->GetSWF();
			phi = Seg->GetPhi(DensPart);
			Lat->PropagateG(Gi_inv,G);
			Lat->ConnectG(Gi_inv,Gi,s,phi);
		}
	}
	if (NG%2==1) s++;

	for (int i=(NG+1)/2; i>=1; i--) {
		for (int j=N2; j>=1; j--) {
			s--;
			Seg = SegBlock2->GetSegment(j);
			G=Seg->GetSWF();
			phi = Seg ->GetPhi(DensPart);
			Lat->PropagateG(Gi_inv,G);
			if (j==N2) {
				for (z=1; z<=M; z++) {
					if (G[z]>0) Gi_inv[z] /=G[z];
					Gi[z][s] *=Gi_inv[z]; //prepare for phi side chains
					phi[z]+=Gi[z][s]*gi[z][N3+1]; //compute density branchpoint
					Gi_inv[z] *= gi[z][N3+1]; //prepare for further propagation of backbone
				}
				if (NG%2==1 && i==(NG+1)/2){
					for (int ii=1; ii<=numDiffSegments; ii++) {
						Seg = Chain->GetDiffSegment(ii);
						G = Seg->GetSWF();
						phi = Seg->GetPhi(DensPart);
						if (NG%2==1) {
							Lat->NormPhiFree(phi,0.5);
						}
					}
				}
			} else {
				//double phiz,sum=0;
				//for (z=1; z<=M; z++) {phiz = Gi_inv[z]*Gi[z][s]; phi[z]+=phiz; if (G[z]>0) sum +=phiz/G[z]; }
				//cout << "s : " << s << " = " << sum << endl;
				for (z=1; z<=M; z++) phi[z] += Gi_inv[z]*Gi[z][s];
			}
		}
	}
	for (int j=N1; j>=1; j--) {
		s--;
		Seg=SegBlock1->GetSegment(j);
		G=Seg->GetSWF();
		phi=Seg ->GetPhi(DensPart);
		Lat->PropagateG(Gi_inv,G);
		//double phiz,sum=0;
		//for (z=1; z<=M; z++) {phiz = Gi_inv[z]*Gi[z][s]; if (G[z]>0) sum +=phiz/G[z]; }
		//cout << "s : " << s << " = " << sum << endl;
		for (int z=1; z<=M; z++) phi[z] += Gi_inv[z]*Gi[z][s];
	}

	lnGN = Lat->ComputeLnGN(Gi_inv);
	for (int ii=1; ii<=numDiffSegments; ii++) {
		Seg = Chain->GetDiffSegment(ii);
		G = Seg->GetSWF();
		phi = Seg->GetPhi(DensPart);
		if (NG%2==1) {
			Lat->NormPhiFree(phi,2);
		}
		for (z=1; z<=M; z++) {
			if (G[z]>0) phi_long[z] +=phi[z]/G[z]; //backbone
		}
	}
	for (i=(NG+1)/2; i>=1; i--) {
		for (z=1; z<=M; z++) Gi_inv[z]=Gi[z][N1+i*N2];
		switch (side_type) {
			case lin:
				MatrixLin(gi,Gi_inv,DensPart);
				break;
			case den:
				MatrixDen(gi,Gi_inv,DensPart);
				break;
			case a_den:
				SF_SegmentBlock* SegBlock;
				SegBlock =  Chain->GetSegmentBlock(1,1);
				Seg = SegBlock->GetSegment(1);
				G = Seg->GetSWF();
				phi = Seg->GetPhi(DensPart);
				Lat->PropagateG(Gi_inv,G);
				for (z=1; z<=M; z++) if (G[z]>0) Gb[z]=Gi_inv[z]/G[z];
				for (z=1; z<=M; z++) phi[z] += gi[z][N3]*gli[z][NN]*Gb[z];
				ComputePhi(1,gli,gi,Gb,DensPart);
				break;
			case bra:
				gi_inv.Dim(1,M);
				for (z=1; z<=M; z++) gi_inv[z]=Gi_inv[z];
				Seg=Chain->GetSegment(1,1);
				G=Seg->GetSWF();
				Lat->PropagateG(gi_inv,G);
				ComputePhiNode(1,0,DensPart);
				break;
			default:
				Message(fatal,"Unknown side chain in comb.");
				break;
		}
		if (i==(NG+1)/2 && NG%2==1) {
			for (int ii=1; ii<=numDiffSegments; ii++) {
				Seg = Chain->GetDiffSegment(ii);
				phi = Seg->GetPhi(DensPart);
				Lat->NormPhiFree(phi,0.5);
			}
		}
	}
	for (int ii=1; ii<=numDiffSegments; ii++) {
		Seg = Chain->GetDiffSegment(ii);
		G = Seg->GetSWF();
		phi = Seg->GetPhi(DensPart);
		if (NG%2==1) {
				Lat->NormPhiFree(phi,2);
			}
		//Lat->CorrectDoubleCountG(phi,G);
		for (int z=1; z<=M; z++) {
			if (G[z] > 0) {
				phi[z] /=G[z];
				phi_short[z] +=phi[z]; //teeth + backbone
			}
		}
	}
	for (z=1; z<=M; z++) {
		phi_short[z] -= phi_long[z]; //teeth only
	}

	if (NG%2==0) {
		for (int ii=1; ii<=numDiffSegments; ii++) {
			Seg = Chain->GetDiffSegment(ii);
			G = Seg->GetSWF();
			phi = Seg->GetPhi(DensPart);
			Lat->NormPhiFree(phi,2);

		}
		Lat->NormPhiFree(phi_short,2);
		Lat->NormPhiFree(phi_long,2);
	}

	// and do the various renormalizations
	if (freedom == fixedTheta) {
		lnCt = log(theta) - lnGN - log(1.*N);
		phiBulk = theta/exp(lnGN);
		SetPhiBulk(phiBulk);
	} else if (freedom == rangeRestricted) {
		double phiPreNorm = 0;
		for (i=1; i<=numDiffSegments; i++) {
			Seg = Chain->GetDiffSegment(i);
			phi = Seg->GetPhi(DensPart);
			for (z=1; z<=M; z++) {
				if (restrictedRange->InRange(z)) {
					phiPreNorm += phi[z];
				}
			}
		}
		lnCt = log(phiRange/phiPreNorm);
	} else {
		lnCb = log(phiBulk/N);
	}
	for (i=1; i<=numDiffSegments; i++) {
		Seg = Chain->GetDiffSegment(i);
		phi = Seg->GetPhi(DensPart);
		if (freedom == fixedTheta) {
			//Lat->NormPhiRestr(phi,Gi_inv,theta/N);
			double C=theta/(N*exp(lnGN));
			for (int z=1; z<=M; z++) {
				phi[z] *= C;
				phi_long[z] *=C;
				phi_short[z]*=C;
			}
		} else if (freedom == rangeRestricted) {
			//Lat->NormPhiFree(phi,exp(lnCt));
			double C=exp(lnCt);
			for (int z=1; z<=M; z++) {
				phi[z] *= C;
				phi_long[z] *=C;
				phi_short[z]*=C;
			}
		} else {
			//Lat->NormPhiFree(phi,phiBulk/N);
			double C=phiBulk/N;
			for (int z=1; z<=M; z++) {
				phi[z] *= C;
				phi_long[z] *=C;
				phi_short[z]*=C;
			}

		}
		G = Seg->GetSWF();
		Lat->RestoreFromSafe(G);
		Lat->RestoreFromSafe(phi);
	}
	if (freedom == rangeRestricted) {
		theta = ComputeTheta();
		phiBulk = theta/exp(lnGN);
		SetPhiBulk(phiBulk);
	}
}

void
SF_Comb1stO::Matrix2ndGen(const LatticeRange* Range,
							const DensityPart DensPart1,
							const DensityPart DensPart2) {

}
void
SF_Comb1stO::MatrixBulk(const DensityPart DensPart) {

}
void
SF_Comb1stO::CopyBulkBoundaries(const DensityPart DensPart) {
	for (int i=1; i<=numDiffSegments; i++) {
		Seg = Chain->GetDiffSegment(i);
		Vector phi = Seg->GetPhi(DensPart);
		Vector PhiBulk = Seg->GetPhi(bulk);
		for (int z=1; z<=M; z++) {
			if (Lat->BulkBoundary(z)) {
				phi[z] = PhiBulk[z];
			}
		}
	}
}

