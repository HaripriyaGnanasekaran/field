#include "SF_Branched1stO.h"
#include <iostream>
#ifdef CUDA
#include "Cuda_tools.h"
#endif

SF_Branched1stO::SF_Branched1stO(Text name_, SF_SegmentList* SegQ_,  Lat1stO* Lat_, Input* MyInput_)
: SF_Molecule(name_, SegQ_, Lat_, MyInput_)
{

	Vector vec; // temporary variable for the vector
	if (Chain->GetMoleculeType() != branched) {
		Message(fatal,MyInput,
		" SF_Branched1stO::SF_Branched1stO: Programming error, trying to create a branched polymer for molecule '" +
		name + "'");
	}
	symmetric = Chain->Symmetric();
	//if (MayerSaupeSet()) {
		//cout << "punt 1" << endl;
	//	Message(fatal,MyInput,"In SF_Branched1stO: Mayer_Saupe not implemented for first order propagators");
	//}
	force_set = (!GetForce()==0);
	if (force_set) {Message(fatal,MyInput,"Force ensemble not implemented in SF_Brancehd1stO");
	}

#ifdef CUDA
	if (GetGPU()) {GPUactivated = true;
		//Message(literal, "Still under construction for GPU; Outcome uncertain. In case of trouble contact frans or do not invoke GPU" );
	}else GPUactivated = false;
#else
	GPUactivated = false;
	if (GetGPU()) {
		Message(literal,"Compile with nvcc (option CUDA=1) to activate GPU... Contact Frans Leermakers. Going classical");
	}
#endif

	numDiffSegments = Chain->GetNumDiffSegments();
	Lat = Lat_;
	// Check for things that are not implemented or do not make sense for a branched molecule
	if (freedom == secondGeneration)
Message(fatal,MyInput,"SF_Branched1stO::SF_Branched1stO: Second Generation not implemented for branched copolymer");
	if (saveMemory)
		Message(fatal,MyInput,"SF_Branched1stO::SF_Branched1stO: SaveMemory not implemented for branched copolymer");
	#if !DEBUG
	//if (Lat->OverflowProtection())
	//	Message(fatal,MyInput,"SF_Branched1stO::SF_Branched1stO: OverflowProtection not implemented for branched copolymer");
	#endif
	CreatePhi(total);
	// calculate the size of the full matrix
	SetPhiBulk(phiBulk);
	// use the information from the LinkNodeTree
	SF_LinkNodeTree *tree=Chain->GetLinkNodeTree();
	n_links=tree->n_links;
	n_nodes=n_links+1;
	Links=tree->LinkList;
	Nodes=tree->NodeList;
	#if DEBUG
	tree->Dump();
	Message(debug,MyInput,"End of the constructor for SF_Branched1stO");
	#endif

}

SF_Branched1stO::~SF_Branched1stO() {
		DeletePhi(total);
}

void
SF_Branched1stO::ReComputePhi() {
	Message(fatal,"ReComputePhi not implemented in SF_Branched1stO");
}

Vector
SF_Branched1stO::GetLong(){
	Message(fatal,"GetLong not implemented in SF_Branched1stO");
	Vector x;
	return x;
}
Vector
SF_Branched1stO::GetShort(){
	Message(fatal,"GetShort not implemented in SF_Branched1stO");
	Vector x;
	return x;
}
Vector
SF_Branched1stO::GetBondOrientation( const DensityPart DensPar) {
	Message(fatal,"GetBondOrientation not implemented in SF_Branched1stO");
	Vector x;
	return x;
};

MoleculeType
SF_Branched1stO::GetMoleculeType() const {
	return Chain->GetMoleculeType();
}

void
SF_Branched1stO::ComputePhi() {
#ifdef DEBUG
	double SumPhi=0; // for debugging only
	int z;
#endif
	Text Segname=Blanks(3); // debug variable
	Vector phiTot; // the actual system storage space for PhiTotal
	Vector PhiTmp; // local temporary storage for PhiTotal
	if (Lat->GetNamesBulkBoundaries().Upperbound() > 0) {
		// descends to SF_Molecule->SF_MolSegemnt->SF_MolState
		CreatePhi(bulk);
		MatrixBulk(bulk);
	}
	int i;
	for (i=1; i<=numDiffSegments; i++) {
		Seg = Chain->GetDiffSegment(i); // return i-th member of the SegmentQ of Chain - the i-th homo-block
		Vector G = Seg->GetSWF();
		Lat->SetBoundaries(G);
	}
	// decide what kind of matrix to use for G
	if (GPUactivated) CudaMatrix1(total); else Matrix1(total);
	if (Chain->SomeSegmentsPinned()) SetPhiBulk(0);
	if (Chain->SomeSegmentsGrafted()) SetPhiBulk(0);
	if (Lat->GetNamesBulkBoundaries().Upperbound() > 0) { // Upperbound() is the size of the array
		CopyBulkBoundaries(total);
		DeletePhi(bulk);
	}
	/* WARNING: Calling GetPhi(Total) for a homopolymer returns an existing
	data object from DensityPartQ, while for a copolymer it creates a new vector
	which is filled by values form DensityPartQ. */
	phiTot = GetPhi(total);
	PhiTmp.Dim(1,M); // sets PhiTmp[z]=0.0
	for (i=1; i<=numDiffSegments; i++) {
		#if DEBUG
		SumPhi=0;
		#endif
		Seg = Chain->GetDiffSegment(i);
		Vector phi = Seg->GetPhi(total);
		#if DEBUG
		printf("\nDifferent segment number %d\n",i);
		fflush(stdout);
		#endif
		PhiTmp += phi;
		#if DEBUG
		Segname=Seg->GetName();
		Message(debug,"Segment name " + Segname);
		for (z=2; z<M; z++) SumPhi+=phi[z];
		phi._print();
		printf("SumPhi=%e\n",SumPhi);
		fflush(stdout);
		#endif
	}
	phiTot.Copy(PhiTmp); // Copy the calculated Phi from the local to system storage
	#if DEBUG
	SumPhi=0;
	for (z=2; z<M; z++) SumPhi+=phiTot[z];
	printf("\nFinal phiTot=\n");
	phiTot._print();
	printf("\nSumPhi=%1.9e\n\n",SumPhi);
	fflush(stdout);
	#endif
	return;
}

//! This function is identical to that in SF_Copol1stO
void
SF_Branched1stO::GetLayerAnalysis(Output* Out) {
	if (LatRangeStartLoops == NULL || LatRangeTrains == NULL) {
		return;
	}
	#if DEBUG
	Message(debug,MyInput,"SF_Branched1stO::GetLayerAnalysis not implemented.");
	#endif
	int i,z,s;

	Matrix Gads(1,M,1,N);
	Vector Gads_inv(1,M);
	Matrix Gfree(1,M,1,N);
	Vector Gfree_inv(1,M);
	Vector G;
	Vector phiTrainLoops;
	Vector phiTails;
	Vector phiFree;
	CreatePhi(trainLoops);
	CreatePhi(tails);
	CreatePhi(unadsorbed);
	for (i=1; i<=numDiffSegments; i++) {
		Seg = Chain->GetDiffSegment(i);
		G = Seg->GetSWF();
		Lat->MakeSafe(G);
		phiTrainLoops = Seg->GetPhi(trainLoops);
		phiTails = Seg->GetPhi(tails);
		phiFree = Seg->GetPhi(unadsorbed);
		Lat->MakeSafe(phiTrainLoops);
		Lat->MakeSafe(phiTails);
		Lat->MakeSafe(phiFree);
	}
	G = Chain->GetSegment(1)->GetSWF();
	Lat->Init2G(Gads_inv,Gfree_inv,G,LatRangeTrains);
	for (z=1; z<=M; z++) {
		Gads[z][1] = Gads_inv[z];
		Gfree[z][1] = Gfree_inv[z];
	}
	for (s=2; s<=N; s++) {
		G = Chain->GetSegment(s)->GetSWF();
		Lat->Propagate2G(Gads,Gfree,G,s,LatRangeTrains);
	}
	Seg = Chain->GetSegment(N);
	phiTrainLoops = Seg->GetPhi(trainLoops);
	phiTails = Seg->GetPhi(tails);
	phiFree = Seg->GetPhi(unadsorbed);
	G = Chain->GetSegment(N)->GetSWF();
	Lat->Init2G(Gads_inv,Gfree_inv,G,LatRangeTrains);
	Lat->ConnectG(Gads_inv,Gads,N,phiTrainLoops);
	Lat->Connect2G(Gads_inv,Gads,N,Gfree_inv,Gfree,N,phiTails);
	Lat->ConnectG(Gfree_inv,Gfree,N,phiFree);
	for (s=N-1; s>=1; s--) {
		Seg = Chain->GetSegment(s);
		G = Seg->GetSWF();
		Lat->Propagate2G(Gads_inv,Gfree_inv,G,LatRangeTrains);
		phiTrainLoops = Seg->GetPhi(trainLoops);
		phiTails = Seg->GetPhi(tails);
		phiFree = Seg->GetPhi(unadsorbed);
		Lat->ConnectG(Gads_inv,Gads,s,phiTrainLoops);
		Lat->Connect2G(Gads_inv,Gads,s,Gfree_inv,Gfree,s,phiTails);
		Lat->ConnectG(Gfree_inv,Gfree,s,phiFree);
	}
	// special hack needed here to do it right for overflow protection
	Vector GiNorm(1,M);
	if (freedom == fixedTheta) {
		Seg = Chain->GetSegment(N);
		G = Seg->GetSWF();
		Lat->MakeSafe(GiNorm);
		Lat->ConnectG(G,Gfree,N,GiNorm);
		Lat->ConnectG(G,Gads,N,GiNorm);
		Lat->CorrectDoubleCountG(GiNorm,G);
	}
	for (i=1; i<=numDiffSegments; i++) {
		Seg = Chain->GetDiffSegment(i);
		G = Seg->GetSWF();
		phiTrainLoops = Seg->GetPhi(trainLoops);
		phiTails = Seg->GetPhi(tails);
		phiFree = Seg->GetPhi(unadsorbed);
		Lat->CorrectDoubleCountG(phiTrainLoops,G);
		Lat->CorrectDoubleCountG(phiTails,G);
		Lat->CorrectDoubleCountG(phiFree,G);
		if (freedom == fixedTheta) {
			Lat->NormPhiRestr(phiTrainLoops,GiNorm,theta/N);
			Lat->NormPhiRestr(phiTails,GiNorm,theta/N);
			Lat->NormPhiRestr(phiFree,GiNorm,theta/N);
		} else {
			Lat->NormPhiFree(phiTrainLoops,exp(lnCb));
			Lat->NormPhiFree(phiTails,exp(lnCb));
			Lat->NormPhiFree(phiFree,exp(lnCb));
		}
		Lat->RestoreFromSafe(phiTrainLoops);
		Lat->RestoreFromSafe(phiTails);
		Lat->RestoreFromSafe(phiFree);
		Lat->RestoreFromSafe(G);
	}
	double GNads = exp(Lat->ComputeLnGN(Gads_inv));
	Out->PutReal("mol",name,"GNads",GNads);
	double data = 0;
	Vector phiTrainLoopsTot = GetPhi(trainLoops);
	Vector phiTailsTot = GetPhi(tails);
	Vector phiFreeTot = GetPhi(unadsorbed);
	for (z=1; z<=M; z++) {
		phiTrainLoopsTot[z] = 0;
		phiTailsTot[z] = 0;
		phiFreeTot[z] = 0;
	}
	for (i=1; i<=numDiffSegments; i++) {
		Seg = Chain->GetDiffSegment(i);
		phiTrainLoops = Seg->GetPhi(trainLoops);
		phiTails = Seg->GetPhi(tails);
		phiFree = Seg->GetPhi(unadsorbed);
		for (z=1; z<=M; z++) {
			phiTrainLoopsTot[z] += phiTrainLoops[z];
			phiTailsTot[z] += phiTails[z];
			phiFreeTot[z] += phiFree[z];
		}
		Lat->SubtractBoundaries(phiTrainLoops);
		Lat->SubtractBoundaries(phiTails);
		for (z=1; z<=M; z++) {
			data += (phiTrainLoops[z] + phiTails[z])*Lat->GetNumLatticeSites(z);
		}
		Lat->RestoreBoundaries(phiTrainLoops);
		Lat->RestoreBoundaries(phiTails);
	}
	Out->PutReal("mol",name,"theta ads",data);
	Out->PutProfile("mol",name,"phi trains+loops",phiTrainLoopsTot);
	Out->PutProfile("mol",name,"phi tails",phiTailsTot);
	Out->PutProfile("mol",name,"phi free",phiFreeTot);
	for (i=1; i<=numDiffSegments; i++) {
		Seg = Chain->GetDiffSegment(i);
		phiTrainLoops = Seg->GetPhi(trainLoops);
		phiTails = Seg->GetPhi(tails);
		phiFree = Seg->GetPhi(unadsorbed);
		Out->PutProfile("mol",name,"phi-" + Seg->GetName() + " trains+loops",phiTrainLoops);
		Out->PutProfile("mol",name,"phi-" + Seg->GetName() + " tails",phiTails);
		Out->PutProfile("mol",name,"phi-" + Seg->GetName() + " free",phiFree);
	}
	// size distribution functions, not implemented for an overflow protection
	// and for copolymers with segments that have different interactions
	Vector numTrains(1,N);
	Vector numLoops(1,N);
	Vector numTails(1,N);

	if (Chain->AllInteractionsEqual() && !Lat->OverflowProtection()) {
	Matrix Gads_inv2(1,M,1,N);
	Matrix Gfree_inv2(1,M,1,N);
	G = Chain->GetSegment(N)->GetSWF();
	Lat->Init2G(Gads_inv,Gfree_inv,G,LatRangeTrains);
	for (z=1; z<=M; z++) {
		Gads_inv2[z][1] = Gads_inv[z];
		Gfree_inv2[z][1] = Gfree_inv[z];
	}
	for (s=2; s<=N; s++) {
		Seg = Chain->GetSegment(N-s+1);
		G = Seg->GetSWF();
		Lat->Propagate2G(Gads_inv2,Gfree_inv2,G,s,LatRangeTrains);
	}
	//general version
	// trains
	Array<Vector> G2Seg(1,numDiffSegments);
	Array<Vector> GBackupSeg(1,numDiffSegments);
	for (i=1; i<=numDiffSegments; i++) {
		G2Seg[i].Dim(1,M);
		GBackupSeg[i].Dim(1,M);
	}
	Vector G2(1,M);
	Vector Gused(1,M);
	for (i=1; i<=numDiffSegments; i++) {
		G = Chain->GetDiffSegment(i)->GetSWF();
		Lat->Init2G(Gused,G2,G,LatRangeTrains);
	}
	for (z=1; z<=M; z++) {
		if (LatRangeTrains->InRange(z)) {
			G2.Dim(1,M);
			G2[z] = G[z];
			Lat->SetBoundaries(G2);
			Matrix Gt(1,M,1,N);
			for (int zq=1; zq<=M; zq++) {
				Gt[zq][1] = G2[zq];
			}
			for (s=2; s<=N; s++) {
				Lat->PropagateG(Gt,Gused,s);
			}
			for (int z2=1; z2<=M; z2++) {
				double ltrlp1 = Lat->GetLambda(z,z2);
				if (ltrlp1 > 0 && LatRangeStartLoops->InRange(z2)) {
					for (int z3=1; z3<=M; z3++) {
						if (LatRangeTrains->InRange(z3)) {
							for (int z4=1; z4<=M; z4++) {
								double ltrlp2 = Lat->GetLambda(z3,z4);
								if (ltrlp2 > 0
									&& LatRangeStartLoops->InRange(z4)) {
									for (s=1; s<=N; s++) {
										double value = 0;
										if (s<N) {
											value += (Gads[z2][N-s]
												+ Gfree[z2][N-s])/ltrlp2;
        									value += (Gads_inv2[z4][N-s]
												+ Gfree_inv2[z4][N-s])/ltrlp1;
        								} else {
        									value += 1/(ltrlp1*ltrlp2);
        								}
										for (int t=1; t<=N-s-1; t++) {
											value += (Gads[z2][t]
												+ Gfree[z2][t])*(Gads_inv2[z4][N-s-t]
												+ Gfree_inv2[z4][N-s-t]);
										}
										numTrains[s] +=
											Lat->GetNumLatticeSites(z3)*ltrlp1
											*ltrlp2*Gt[z3][s]*value/GNads;
									}
								}
							}
						}
					}
				}
			}
		}
	}
	// loops
	Lat->Init2G(G2,Gused,G,LatRangeTrains);
	for (z=1; z<=M; z++) {
		if (LatRangeStartLoops->InRange(z)) {
			G2.Dim(1,M);
			G2[z] = G[z];
			Lat->SetBoundaries(G2);
			Matrix Gt(1,M,1,N);
			for (int zq=1; zq<=M; zq++) {
				Gt[zq][1] = G2[zq];
			}
			for (s=2; s<=N-2; s++) {
				Lat->PropagateG(Gt,Gused,s);
			}
			for (int z2=1; z2<=M; z2++) {
				double llptr1 = Lat->GetLambda(z,z2);
				if (llptr1 > 0 && LatRangeTrains->InRange(z2)) {
					for (int z3=1; z3<=M; z3++) {
						if (LatRangeStartLoops->InRange(z3)) {
							for (int z4=1; z4<=M; z4++) {
								double llptr2 = Lat->GetLambda(z3,z4);
								if (llptr2 > 0 && LatRangeTrains->InRange(z4)) {
									for (s=1; s<=N-2; s++) {
										double value = 0;
										for (int t=1; t<=N-s-1; t++) {
											value += Gads[z2][t]
												*Gads_inv2[z4][N-s-t];
										}
										numLoops[s] +=
											Lat->GetNumLatticeSites(z3)*llptr1
											*llptr2*Gt[z3][s]*value/GNads;
									}
								}
							}
						}
					}
				}
			}
		}
	}
	// tails
	for (s=1; s<=N-1; s++) {
		for (int ztr=1; ztr<=M; ztr++) {
			if (LatRangeTrains->InRange(ztr)) {
				for (int ztl=1; ztl<=M; ztl++) {
					if (LatRangeStartLoops->InRange(ztl)) {
						numTails[s] += Lat->GetNumLatticeSites(ztl)
							*Lat->GetLambda(ztl,ztr)*Gfree[ztl][s]
							*Gads_inv2[ztr][N-s]/GNads;
						numTails[s] += Lat->GetNumLatticeSites(ztl)
							*Lat->GetLambda(ztl,ztr)*Gfree_inv2[ztl][s]
							*Gads[ztr][N-s]/GNads;
					}
				}
			}
		}
	}
	}
	if (Lat->OverflowProtection()) {
		Message(implementation, "train,loop and tail size distribution not"
			" implemented for overflow protection");
	}
	if (!Chain->AllInteractionsEqual()) {
		Message(implementation, "train,loop and tail size distribution not"
			" implemented for copolymers with segments that differ in interactions");
	}
	if (!Lat->OverflowProtection() && Chain->AllInteractionsEqual()) {
		double frTrains = 0;
		double frLoops = 0;
		double frTails = 0;
		double fluctFrTrains = 0;
		double fluctFrLoops = 0;
		double fluctFrTails = 0;
		for (s=1; s<=N; s++) {
			frTrains += s*numTrains[s];
			frLoops += s*numLoops[s];
			frTails += s*numTails[s];
			fluctFrTrains += s*numTrains[s]*s;
			fluctFrLoops += s*numLoops[s]*s;
			fluctFrTails += s*numTails[s]*s;
		}
		frTrains /= N;
		frLoops /= N;
		frTails /= N;
		fluctFrTrains /= N*N;
		fluctFrLoops /= N*N;
		fluctFrTails /= N*N;
		fluctFrTrains -= frTrains*frTrains;
		fluctFrLoops -= frLoops*frLoops;
		fluctFrTails -= frTails*frTails;
		double numTrainsTot = 0;
		double numLoopsTot = 0;
		double numTailsTot = 0;
		for (s=1; s<=N; s++) {
			numTrainsTot += numTrains[s];
			numLoopsTot += numLoops[s];
			numTailsTot += numTails[s];
		}
		double avLengthTrains = N*frTrains/numTrainsTot;
		double avLengthLoops = N*frLoops/numLoopsTot;
		double avLengthTails = N*frTails/numTailsTot;
		double fluctAvLengthTrains = N*N*fluctFrTrains/(numTrainsTot*numTrainsTot);
		double fluctAvLengthLoops = N*N*fluctFrLoops/(numLoopsTot*numLoopsTot);
		double fluctAvLengthTails = N*N*fluctFrTails/(numTailsTot*numTailsTot);
		Out->PutReal("mol", name, "trains fract", frTrains);
		Out->PutReal("mol", name, "fluct trains fract", fluctFrTrains);
		Out->PutReal("mol", name, "trains number", numTrainsTot);
		Out->PutReal("mol", name, "trains av length", avLengthTrains);
		Out->PutReal("mol", name, "fluct trains av length", fluctAvLengthTrains);
		Out->PutVector("mol", name, "trains fract", numTrains, N);
		Out->PutReal("mol", name, "loops fract", frLoops);
		Out->PutReal("mol", name, "fluct loops fract", fluctFrLoops);
		Out->PutReal("mol", name, "loops number", numLoopsTot);
		Out->PutReal("mol", name, "loops av length", avLengthLoops);
		Out->PutReal("mol", name, "fluct loops av length", fluctAvLengthLoops);
		Out->PutVector("mol", name, "loops fract", numLoops, N);
		Out->PutReal("mol", name, "tails fract", frTails);
		Out->PutReal("mol", name, "fluct tails fract", fluctFrTails);
		Out->PutReal("mol", name, "tails number", numTailsTot);
		Out->PutReal("mol", name, "tails av length", avLengthTails);
		Out->PutReal("mol", name, "fluct tails av length", fluctAvLengthTails);
		Out->PutVector("mol", name, "tails fract", numTails, N);
	} else {
		double frTrains = 0;
		double frLoops = 0;
		double frTails = 0;
		for (i=1; i<=numDiffSegments; i++) {
			Seg = Chain->GetDiffSegment(i);
			phiTrainLoops = Seg->GetPhi(trainLoops);
			phiTails = Seg->GetPhi(tails);
			Lat->SubtractBoundaries(phiTrainLoops);
			Lat->SubtractBoundaries(phiTails);
			for (z=1; z<=M; z++) {
				if (LatRangeTrains->InRange(z)) {
					frTrains += phiTrainLoops[z]*Lat->GetNumLatticeSites(z);
				} else {
					frLoops += phiTrainLoops[z]*Lat->GetNumLatticeSites(z);
					frTails += phiTails[z]*Lat->GetNumLatticeSites(z);
				}
			}
			Lat->RestoreBoundaries(phiTrainLoops);
			Lat->RestoreBoundaries(phiTails);
		}
		double phiads = frTrains + frLoops + frTails;
		frTrains /= phiads;
		frLoops /= phiads;
		frTails /= phiads;
		Out->PutReal("mol", name, "trains fract", frTrains);
		Out->PutReal("mol", name, "loops fract", frLoops);
		Out->PutReal("mol", name, "tails fract", frTails);
	}
	DeletePhi(trainLoops);
	DeletePhi(tails);
	DeletePhi(unadsorbed);
}
// end of GetLayerAnalysis

/*!
Calls ComputeGzsLink to compute Gzs of all links but linkto.
Connects Gzs inside the node and propagates to the next segment of the linkto.
*/
void SF_Branched1stO::ComputeGzsNode(const int nodenum, const int linkto) {
	// declare the variables and allocate space
	int arm; // arm number in the internal numbering of the node (1..n_links)
	int linknum; // link number in the global numbering
	Vector G; // storage space for G(z) of a particular segment
	SF_Node node=Nodes[nodenum]; // to simplify the notation
	SF_Link Linkto; // the link numbered linkto
	if(linkto!=0)Linkto=Links[linkto]; // Links[0] would produce an error
		#if DEBUG
		Text t1=Blanks(2); // text field to be filled for debug messages
		t1.Putint(nodenum);
		Message(debug,MyInput,"ComputeGzsNode[" + t1 + "]");
		#endif
	// end of variable declarations and space allocation

	// Compute Gzs of all links but linkto
	for(arm=node.n_links();arm>0;arm--) {
		linknum=node.link(arm);
		#if DEBUG
		if(Links[linknum].length()==1) { // by now I see no reason why this should happen
			Message(fatal,MyInput,"Programming error, link with single segment");
		}
		#endif
		if(linknum!=linkto) {
			ComputeGzsLink(linknum,nodenum);
			#if DEBUG
			Message(debug,MyInput,"return control to ComputeGzsNode[" + t1 + "]");
			#endif
		}
	}
	// now we do in PropagateGzsNode what has been done here before
	if (linkto!=0)  PropagateGzsNode(nodenum,linkto); // connect arms and proceed to linkto
	return;
}

void SF_Branched1stO::CudaComputeGzsNode(const int nodenum, const int linkto) {
#ifdef CUDA
	// declare the variables and allocate space
	int arm; // arm number in the internal numbering of the node (1..n_links)
	int linknum; // link number in the global numbering
	//Vector G; // storage space for G(z) of a particular segment
	SF_Node node=Nodes[nodenum]; // to simplify the notation
	SF_Link Linkto; // the link numbered linkto
	if (linkto!=0) Linkto=Links[linkto]; // Links[0] would produce an error
	// end of variable declarations and space allocation
	// Compute Gzs of all links but linkto
	for(arm=node.n_links();arm>0;arm--) {
		linknum=node.link(arm);
		if(linknum!=linkto) {
			CudaComputeGzsLink(linknum,nodenum);
		}
	}
	// now we do in PropagateGzsNode what has been done here before
	if (linkto!=0)  CudaPropagateGzsNode(nodenum,linkto); // connect arms and proceed to linkto
#endif
	return;
}

/*!
Calls ComputeGzsNode to provide the Gzs of the other segment than nodefrom.
After it gets Gzs of the other segment, ComputeGzsLink does the propagation from
the 2nd (next to the other node) to the last segment (the nodefrom) of the link.
*/
void SF_Branched1stO::ComputeGzsLink(const int linknum, const int nodefrom) {
	// declare the variables and allocate space
	int s; // index variable
	int segnum; // internal segment number of a link
	int step; // determines the direction of propagation
	int nodeto; // node number to which the link is connected, other than nodefrom
	SF_Link Link=Links[linknum]; // to simplify the notation - the link we are workin on
	int length=Link.length(); // length of our link
	Vector G; // storage space fo G(z)
	VecArray GzsLink=Gzs[linknum]; // the Gzs matrix of the current link
		#if DEBUG
		Text t1=Blanks(2); // text variable for debugging
		t1.Putint(linknum);
		Message(debug,MyInput,"ComputeGzsLink[" + t1 +"]");
		#endif
	// end of variable declarations and space allocation

	// if we came from node(1), the propagation will go from node(2) or the other way round
	if(Link.node(1)==nodefrom) {
		//segnum=length; // the segment number of nodeto
		//nodeto=Link.node(2); // the node number of nodeto
		//step = -1;
		segnum=length; nodeto=Link.node(2); step =-1;
	} else {
		//segnum=1;
		//nodeto=Link.node(2);
		//step = 1;
		segnum=1; nodeto=Link.node(2); step = 1;
	}
	ComputeGzsNode(nodeto,linknum); // Ask nodeto for Gzs
		#if DEBUG
		Message(debug,MyInput,"Return control to ComputeGzsLink[" + t1 +"]");
		#endif
	// now propagate through the link

	//step=(segnum==1) ? 1 : -1; // determine the direction of propagation
	// start from the 2nd segment, Gzs of the 1st segment been computed by ComputeGzsNode
	for (s=2; s<=length; s++) { // s is an index to make the proper number of cycles
		segnum+=step;  // first increment segnum, the Gzs of the end segment has been provided by the node
		G = Chain->GetSegment(linknum,segnum)->GetSWF();
		Lat->PropagateG(GzsLink[segnum-step], G, GzsLink[segnum]);
			#if DEBUG
			printf("Propagate link %d, segment %d\n",linknum,segnum);
			PrintVec(G,GzsLink[segnum-step],GzsLink[segnum],"G","GzsLink[segnum-step]","GzsLink[segnum]");
			#endif
	}
		#if DEBUG
		printf("\n"); fflush(stdout);
		#endif
	return;
}


void SF_Branched1stO::CudaComputeGzsLink(const int linknum, const int nodefrom) {
#ifdef CUDA
	// declare the variables and allocate space
	int s; // index variable
	int segnr;
	int segnum; // internal segment number of a link
	int step; // determines the direction of propagation
	int nodeto; // node number to which the link is connected, other than nodefrom
	SF_Link Link=Links[linknum]; // to simplify the notation - the link we are workin on
	int length=Link.length(); // length of our link
	//Vector G; // storage space fo G(z)

        //VecArray GzsLink=Gzs[linknum]; // the Gzs matrix of the current link

	// end of variable declarations and space allocation

	// if we came from node(1), the propagation will go from node(2) or the other way round
	if(Link.node(1)==nodefrom) {
		//segnum=length; // the segment number of nodeto
		//nodeto=Link.node(2); // the node number of nodeto
		//step = -1;
		segnum=length; nodeto=Link.node(2); step =-1;
	} else {
		//segnum=1;
		//nodeto=Link.node(2);
		//step = 1;
		segnum=1; nodeto=Link.node(2); step = 1;
	}
	CudaComputeGzsNode(nodeto,linknum); // Ask nodeto for Gzs

	// now propagate through the link
	//step=(segnum==1) ? 1 : -1; // determine the direction of propagation
	// start from the 2nd segment, Gzs of the 1st segment been computed by ComputeGzsNode
	for (s=2; s<=length; s++) { // s is an index to make the proper number of cycles
		segnum+=step;  // first increment segnum, the Gzs of the end segment has been provided by the node
		//G = Chain->GetSegment(linknum,segnum)->GetSWF();
		segnr=1; while (Chain->GetSegment(linknum,segnum) != Chain->GetDiffSegment(segnr)) segnr++;
		Propagate(Gs[linknum]+segnum-step, Gs[linknum]+segnum, M, segnr, Dgi, Dg, 0.0, Info);
		//Lat->PropagateG(GzsLink[segnum-step], G, GzsLink[segnum]);
	}
#endif
	return;
}

/*!
Multiply Gzs from all links but linkto, correct for double counting and stroe the result
in Gzs[linkto][segto] It is only called if a node has more than 1 arm.
*/
void SF_Branched1stO::PropagateGzsNode(const int nodenum, const int linkto) {
	// declare the variables and allocate space
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
		#if DEBUG
		Text t1=Blanks(2); // text variable for debugging
		t1.Putint(nodenum);
		Message(debug,MyInput,"PropagateGzsNode[" + t1 +"]");
		#endif
	// end of variable declarations and space allocation

//	Gout=Gzs[linkto][segto]; // that is where we store the result
	G=Chain->GetSegment(linkto,segto)->GetSWF(); // GetSWF() of any of the node segments
	if(n_links==1) Gout.Copy(G); // for the end-node, Gout=G
	else { // go through all the arms but linkto and connect G
		for(arm=1;arm<=n_links;arm++) {
			link=Node.link(arm);
				#if DEBUG
				t1.Putint(link);
				Message(debug,MyInput,"Connect link[" + t1 + "]");
				#endif
			if(link!=linkto) {
				seg=Links[link].seg(nodenum);
				if(first_arm) {
					first_arm=0;
					TempG.Copy(Gzs[link][seg]); // for the first arm, copy Gzs[link][seg] to Gout
						#if DEBUG
						printf("This is First arm\n\n"); fflush(stdout);
						#endif
				}
				else {
					Gout.Dim(1,M);
					Lat->MakeSafe(Gout);

					#if DEBUG
					printf("Before connection of link=%d, seg=%d:\n",link,seg);
					PrintVec(Gzs[link][seg],TempG,Gout,"Gzs[link][seg]","TempG","Gout");
					#endif

					Lat->ConnectG(Gzs[link][seg],TempG,Gout);

					#if DEBUG
					printf("Before correction of link=%d, seg=%d :",link,seg);
					PrintVec(G,Gzs[link][seg],Gout,"G","Gzs[link][seg]","Gout");
					#endif

					Lat->CorrectDoubleCountG(Gout,G); // divide by G

					#if DEBUG
					printf("After correction of link=%d, seg=%d:\n", link,seg);
					PrintVec(G,Gzs[link][seg],Gout,"G","Gzs[link][seg]","Gout");
					#endif

					TempG.Copy(Gout); // copy the Gout to TempG
				}
				#if DEBUG
				printf("After connection of link=%d, seg=%d:\n",link,seg);
				PrintVec(G,Gzs[link][seg],TempG,"G","Gzs[link][seg]","TempG");
				#endif
			}
			#if DEBUG
			else { printf("This is linkto\n\n"); fflush(stdout); }
			#endif
		}
	}
	Gzs[linkto][segto].Copy(Gout);

	#if DEBUG
	printf("Propagation of node %d finished:\n",nodenum);
	PrintVec(G,Gout,Gzs[linkto][segto],"G","Gout","Gzs[linkto][segto]");
	#endif
return;
}

void SF_Branched1stO::CudaPropagateGzsNode(const int nodenum, const int linkto) {
#ifdef CUDA
	// declare the variables and allocate space
	int first_arm=1; // boolean to see if
	int arm; // link number in the internal numbering inside the node (1..n_links)
	int link; // link number in the global numbering
	int seg; // segment number of the link
	int segnr;
	double *DGout =   (double *)AllocateMemoryOnDevice(M);
	double *DTempG =   (double *)AllocateMemoryOnDevice(M);
//	Vector Gout(1,M);  // Gzs of linkto - that where we store the output
//	Vector G;          // G(z) of the segment
//	Vector TempG(1,M); // temporary G Vector
	SF_Node Node=Nodes[nodenum]; // the node we are working on
	int n_links=Node.n_links(); // number of links of the Node
	int segto=Links[linkto].seg(nodenum); // segment number by which linkto is connected to the Node

	// end of variable declarations and space allocation

	//Gout=Gzs[linkto][segto]; // that is where we store the result; Here I use DGout
	//G=Chain->GetSegment(linkto,segto)->GetSWF(); // GetSWF() of any of the node segments
	segnr=1; while (Chain->GetSegment(linkto,segto) != Chain->GetDiffSegment(segnr)) segnr++;
	if(n_links==1) {
		//Gout.Copy(G); // for the end-node, Gout=G
		InitializeForward(1, M, segnr, DGout, Dg);
	}
	else { // go through all the arms but linkto and connect G
		for(arm=1;arm<=n_links;arm++) {
			link=Node.link(arm);

			if(link!=linkto) {
				seg=Links[link].seg(nodenum);
				if(first_arm) {
					first_arm=0;
					InitializeBackward(M,Gs[link]+seg,DTempG, Dgi);
					//TempG.Copy(Gzs[link][seg]); // for the first arm, copy Gzs[link][seg] to Gout
				}
				else {
					//Gout.Dim(1,M);
					//Lat->MakeSafe(Gout);
					Zero(M,DGout);
					//Lat->ConnectG(Gzs[link][seg],TempG,Gout);
					Composition(Gs[link]+seg, M, 1, DTempG, Dgi,DGout);
					//Lat->CorrectDoubleCountG(Gout,G); // divide by G
					CorrectDoubleCounting(M,1,segnr,DGout,Dg);
					//TempG.Copy(Gout); // copy the Gout to TempG
					InitializeForward(1, M, 1, DTempG, DGout);
				}
			}
		}
	}
	//Gzs[linkto][segto].Copy(Gout);
	InitializeForward(Gs[linkto]+segto,M,1,Dgi,DGout);

	FreeMemoryOnDevice(DGout);
	FreeMemoryOnDevice(DTempG);

#endif
return;
}

/*!
WARNING: Phi is not only the node but sum of Phi over all segments of the same
type. Moreover, we do not compute Phi but Phi*G and the correction for double counting
is performed only once per each different segment at the end.
We store Phi(Node)*G locally.

ComputePhiNode re-computes the Gi_inv for individual links by
dividing PhiNode*G by the Gzs(forward) of the corresponding link.
*/
void SF_Branched1stO::ComputePhiNode(int nodenum, int linkfrom, const DensityPart DensPart) {
	// variable declarations and space allocation
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
		#if DEBUG
		Text t1=Blanks(2); // text variable for debugging
		Text t2=Blanks(2); // text variable for debugging
		t1.Putint(nodenum);
		t2=Seg->GetName();
		Message(debug,MyInput,"ComputePhiNode[" + t1 +"] segment name " + t2);
		#endif
	// end of variable declarations and space allocation

	// Compute phi of the node itself
	PhiNode.Dim(1,M); // initialize PhiNode[z] to 0.0
	Lat->MakeSafe(PhiNode);
	TempPhi.Dim(1,M); // initialize PhiNode[z] to 0.0
	Lat->MakeSafe(TempPhi);
	link=Node.link(1); // no matter which link we take, GetSegment gives the same result for the whole node
	seg=Links[link].seg(nodenum);
	Seg=Chain->GetSegment(link,seg); // it does not matter from which link of the node we call GetSWF
	G=Seg->GetSWF(); // it makes no difference to which segment of the node we call GetSWF()
	Phi = Seg->GetPhi(DensPart); // that is where we store the phi in the end
	// end of the general part that is indepenent of node structure

	#if DEBUG
	t2=Seg->GetName();
	Message(debug,MyInput,"ComputePhiNode[" + t1 +"] segment name " + t2);
	#endif

	if (n_links==1) { // it is either a start or end-node
		if(nodenum==1) link=1; // it is a start node
		else link=linkfrom; // if not, Gi_inv of linkfrom is provided by ComputePhiLink
		seg=Links[link].seg(nodenum);
		if(nodenum==1) Gi_inv.Copy(G); // initialize the Gi_inv for the 1st node

		#if DEBUG
		printf("Before connection of linkfrom, link=%d, seg=%d :\n",link,seg);
		PrintVec(G,Gi_inv,Gzs[link][seg],PhiNode,"G","Gi_inv","Gzs[link][seg]","PhiNode");
		#endif

		Lat->ConnectG(Gi_inv,Gzs[link][seg],PhiNode); // first do that for the linkfrom

		#if DEBUG
		printf("After connection of linkfrom, link=%d, seg=%d :\n", link,seg);
		PrintVec(G,Gi_inv,Gzs[link][seg],PhiNode,"G","Gi_inv","Gzs[link][seg]","PhiNode");
		#endif

	} else { // node has more links
		seg=Links[linkfrom].seg(nodenum);
		TempPhi.Copy(Gi_inv); // Copy the Gzs of linkfrom to TempPhi

		#if DEBUG
		PrintVec(TempPhi,Gi_inv,"TempPhi","Gi_inv");
		#endif

		for(arm=1;arm<=n_links && linkfrom!=0;arm++) {// through all links but linkfrom
			link=Node.link(arm);
			if(link!=linkfrom) { // for linkfrom we have Gi_inv
				PhiNode.Dim(1,M); // initialize PhiNode to zero
				Lat->MakeSafe(PhiNode);
				seg=Links[link].seg(nodenum);

				#if DEBUG
				printf("Before connection\n");
				PrintVec(TempPhi,Gzs[link][seg],PhiNode,"TempPhi","Gzs[link][seg]","PhiNode");
				#endif

				Lat->ConnectG(TempPhi,Gzs[link][seg],PhiNode);

				#if DEBUG
				printf("After connection\n");
				PrintVec(G,Gzs[link][seg],PhiNode,"G","Gzs[link][seg]","PhiNode");
				#endif

				TempPhi.Copy(PhiNode);
			}
		}

		#if DEBUG
		printf("Correct double counting, before:\n");
		PrintVec(G,PhiNode,"G","PhiNode");
		#endif

		// perform n_links-2 corrections for double counting
		for(arm=2;arm<n_links;arm++) {
			Lat->CorrectDoubleCountG(PhiNode,G);

		#if DEBUG
		printf("After correction %d\n",arm-1);
		PrintVec(G,Gzs[link][seg],PhiNode,"G","Gzs[link][seg]","PhiNode");
		#endif
		}
	}
	// FIXME (?) this will only work with overflow protection switched off
	if(Lat->OverflowProtection()) {
		// this is a workaround
		TempPhi.Dim(1,M); // TempPhi[z]=0.0, in logarithmic form this means G[z]=1.0
		Lat->ConnectG(TempPhi,PhiNode,Phi); // multiply PhiNode by TempPhi and add to Phi
	} else Phi += PhiNode; // simply add the computed PhiNode to the Phi of the segment

	#if DEBUG
	printf("Resultant Phi - may be safe:\n");
	PrintVec(Phi,PhiNode,"Phi","PhiNode");
	#endif

	// Compute Phi of all the arms but linkfrom
	for(arm=n_links;arm>0;arm--) {
		linkto=Node.link(arm);
		if(linkto != linkfrom) { // it is not linkfrom
			seg=Links[linkto].seg(nodenum);
			Gi_inv.Copy(PhiNode); // Copy Phi to Gi_inv

			#if DEBUG
			printf("Before recalculation of Gi_inv\n");
			PrintVec(Gzs[linkto][seg],Gi_inv,PhiNode,"Gzs[linkto][seg]","Gi_inv","PhiNode");
			#endif

			Lat->CorrectDoubleCountG(Gi_inv,Gzs[linkto][seg]); // divide by Gzs of segto

			#if DEBUG
			printf("After recalculation of Gi_inv\n");
			PrintVec(Gzs[linkto][seg],Gi_inv,PhiNode,"Gzs[linkto][seg]","Gi_inv","PhiNode");
			#endif

			ComputePhiLink(linkto,nodenum,DensPart);

			#if DEBUG
			Message(debug,MyInput,"Return control to ComputePhiNode["+t1+"]");
			#endif
		}
	}
	return;
}

void SF_Branched1stO::CudaComputePhiNode(int nodenum, int linkfrom, const DensityPart DensPart) {
#ifdef CUDA
	// variable declarations and space allocation
	int arm; // internal number of an arm of the node
	int link; // link number in the absolute numbering
	int seg; // segment number by which the link is connected to the Node
	int linkto; // link number to which we will propagate
	int segnr;

	double *DPhiNode =   (double *)AllocateMemoryOnDevice(M);
	double *DTempPhi =   (double *)AllocateMemoryOnDevice(M);
	//Vector Phi; // Phi(z) vector of the node segment
	//Vector PhiNode; // locally stores computed Phi*G of the node
	//Vector TempPhi; // temporary storage space for PhiNode
	//Vector G; // G(z)
	SF_Node Node=Nodes[nodenum]; // the node we are working on
	int n_links=Node.n_links(); // number of links of the Node

	// end of variable declarations and space allocation

	// Compute phi of the node itself
	//PhiNode.Dim(1,M); // initialize PhiNode[z] to 0.0
	Zero(M,DPhiNode);

	//Lat->MakeSafe(PhiNode);
	//TempPhi.Dim(1,M); // initialize PhiNode[z] to 0.0
	Zero(M,DTempPhi); // I will use DTempPhi for TempPhi
	//Lat->MakeSafe(TempPhi);
	link=Node.link(1); // no matter which link we take, GetSegment gives the same result for the whole node

	seg=Links[link].seg(nodenum);
	//Seg=Chain->GetSegment(link,seg); // it does not matter from which link of the node we call GetSWF
	//G=Seg->GetSWF(); // it makes no difference to which segment of the node we call GetSWF()
	segnr=1; while (Chain->GetSegment(link,seg) != Chain->GetDiffSegment(segnr)) segnr++;
	//Phi = Seg->GetPhi(DensPart); // that is where we store the phi in the end
	// end of the general part that is indepenent of node structure


	if (n_links==1) { // it is either a start or end-node
		if(nodenum==1) link=1; // it is a start node
		else link=linkfrom; // if not, Gi_inv of linkfrom is provided by ComputePhiLink
		seg=Links[link].seg(nodenum);

		if(nodenum==1) {
			//Gi_inv.Copy(G); // initialize the Gi_inv for the 1st node
			InitializeBackward(M, segnr, Dgi_inv, Dg);
		}

		//Lat->ConnectG(Gi_inv,Gzs[link][seg],PhiNode); // first do that for the linkfrom
		Composition(Gs[link]+seg, M, 1, Dgi_inv, Dgi, DPhiNode);


	} else { // node has more links
		seg=Links[linkfrom].seg(nodenum);//do not know why this line is here....FL
		//TempPhi.Copy(Gi_inv); // Copy the Gzs of linkfrom to TempPhi
		InitializeBackward(M, 1, DTempPhi, Dgi_inv);

		for(arm=1;arm<=n_links && linkfrom!=0;arm++) {// through all links but linkfrom
			link=Node.link(arm);
			if(link!=linkfrom) { // for linkfrom we have Gi_inv
				//PhiNode.Dim(1,M); // initialize PhiNode to zero
				Zero(M,DPhiNode);
				//Lat->MakeSafe(PhiNode);
				seg=Links[link].seg(nodenum);

				//Lat->ConnectG(TempPhi,Gzs[link][seg],PhiNode);
				Composition(Gs[link]+seg, M, 1, DTempPhi, Dgi, DPhiNode);
				//TempPhi.Copy(PhiNode);
				InitializeBackward(M,1,DTempPhi,DPhiNode);
			}
		}


		// perform n_links-2 corrections for double counting
		for(arm=2;arm<n_links;arm++) {
			//Lat->CorrectDoubleCountG(PhiNode,G);
			CorrectDoubleCounting(M, 1, segnr, DPhiNode, Dg);

		}
	}

	//Phi += PhiNode; // simply add the computed PhiNode to the Phi of the segment
	Add(1, M, segnr, DPhiNode, Dphi);

	// Compute Phi of all the arms but linkfrom
	for(arm=n_links;arm>0;arm--) {
		linkto=Node.link(arm);
		if(linkto != linkfrom) { // it is not linkfrom
			seg=Links[linkto].seg(nodenum);
			//Gi_inv.Copy(PhiNode); // Copy Phi to Gi_inv
			InitializeBackward(M,1,Dgi_inv,DPhiNode);
			//Lat->CorrectDoubleCountG(Gi_inv,Gzs[linkto][seg]); // divide by Gzs of segto
			CorrectDoubleCounting(M, 1, Gs[linkto]+seg, Dgi_inv, Dgi);
			CudaComputePhiLink(linkto,nodenum,DensPart);
		}
	}
	FreeMemoryOnDevice(DPhiNode);
	FreeMemoryOnDevice(DTempPhi);
#endif
	return;
}

void SF_Branched1stO::AddPhi(const Vector Phi1, Vector PhiOut) {
	int z;
	for(z=1;z<=M;z++) {
		if (Phi1[z] > LOGMAXDOUBLE-1 || PhiOut[z] == LOGZERO) PhiOut[z] = Phi1[z];
		else if(PhiOut[z]<LOGMAXDOUBLE) PhiOut[z] = log1p(exp(PhiOut[z])+exp(Phi1[z]));
		// if (PhiOut >=LOGMAXDOUBLE) do not change its value
	}
	return;
}

/*! Comments from Goliath:

This pocedure computes phi of a link between two nodes. The phi of the nodes is not computed.
At the end of the procedure, the probability of the node with this link as the only end
is stored in Gi_inv to be used later in ComputePhiNode
*/
void  SF_Branched1stO::ComputePhiLink(int linknum, int nodefrom, const DensityPart DensPart) {
	// variable declarations and space allocation
	int segnum; // segment number in internal numbering of the link
	int step; // direction of propagation
	int nodeto; // number of the other node (other than nodefrom) to which the link is connected
	int s; // index variable
	VecArray GzsLink=Gzs[linknum]; // the Gzs matrix of the link
	SF_Link Link=Links[linknum]; // to simplify the notation - this link
	int length=Link.length(); // length of the Link
	Vector G; // G(z)
	Vector Phi;
		#if DEBUG
		Text t1=Blanks(1);
		t1.Putint(linknum);
		Message(debug,MyInput,"ComputePhiLink[" + t1 + "]");
		#endif
	// end of variable declarations and space allocation
	segnum=Link.seg(nodefrom); // segment by which the link is connected to the node
	nodeto=(segnum==1) ? Link.node(2) : Link.node(1); // nodeto is the other node
	step=(segnum==1) ? 1 : -1; // propagate forwards or backwards
	// Phi and Gi_inv of the first segment has been computed by ComputePhiNode
	for(s=2;s<length;s++) {
		segnum+=step;
		Seg=Chain->GetSegment(linknum,segnum);
		G=Seg->GetSWF();
		Phi=Seg->GetPhi(DensPart);
			#if DEBUG
			printf("Link %d, Segment %d\n",linknum,segnum);
			printf("Before propagation\n      Gi_inv:         G:\n");
			PrintVec(Gi_inv,G,"Gi_inv","G");
			#endif
		Lat->PropagateG(Gi_inv,G,Gi_inv);
			#if DEBUG
			printf("After propagation\n       Gzs_sym[%d][%d]:  Gi_inv:           Gzs[%d][%d]:       Phi:\n",linknum,length-segnum+1,linknum,segnum);
			PrintVec(GzsLink[length-segnum+1],Gi_inv,GzsLink[segnum],Phi,"GzsLink[length-segnum+1]","Gi_inv","GzsLink[segnum]","Phi");
			#endif
		Lat->ConnectG(Gi_inv,GzsLink[segnum],Phi);
			#if DEBUG
			printf("After connection\n       Phi:           Gi_inv:\n");
			PrintVec(Phi,Gi_inv,"Phi","Gi_inv");
			#endif
	}
	// for the last segment, just calculate Gi_inv, but do not compute phi
	segnum+=step;
	Seg=Chain->GetSegment(linknum,segnum);
	G=Seg->GetSWF();
	Lat->PropagateG(Gi_inv,G,Gi_inv);
		#if DEBUG
		printf("After last segment\n");
		PrintVec(Gi_inv,"Gi_inv");
		fflush(stdout);
		#endif
	ComputePhiNode(nodeto,linknum,DensPart);
		#if DEBUG
		Message(debug,MyInput,"Return control to ComputePhiLink[" + t1 +"]");
		#endif
	return;
}

void  SF_Branched1stO::CudaComputePhiLink(int linknum, int nodefrom, const DensityPart DensPart) {
#ifdef CUDA
	// variable declarations and space allocation
	int segnum; // segment number in internal numbering of the link
	int step; // direction of propagation
	int nodeto; // number of the other node (other than nodefrom) to which the link is connected
	int s; // index variable
	int segnr;
	//VecArray GzsLink=Gzs[linknum]; // the Gzs matrix of the link
	SF_Link Link=Links[linknum]; // to simplify the notation - this link
	int length=Link.length(); // length of the Link
	//Vector G; // G(z)
	//Vector Phi;

	// end of variable declarations and space allocation
	segnum=Link.seg(nodefrom); // segment by which the link is connected to the node
	nodeto=(segnum==1) ? Link.node(2) : Link.node(1); // nodeto is the other node
	step=(segnum==1) ? 1 : -1; // propagate forwards or backwards
	// Phi and Gi_inv of the first segment has been computed by ComputePhiNode
	for(s=2;s<length;s++) {
		segnum+=step;
		//Seg=Chain->GetSegment(linknum,segnum);
		segnr=1; while (Chain->GetSegment(linknum,segnum) != Chain->GetDiffSegment(segnr)) segnr++;
		//G=Seg->GetSWF();
		//Phi=Seg->GetPhi(DensPart);
		//Lat->PropagateG(Gi_inv,G,Gi_inv);
		Backward(M, segnr, Dgi_inv, Dg, Dgx, 0.0, Info);
		//Lat->ConnectG(Gi_inv,GzsLink[segnum],Phi);
		Composition(Gs[linknum]+segnum, M, segnr, Dgi_inv, Dgi, Dphi);

	}
	// for the last segment, just calculate Gi_inv, but do not compute phi
	segnum+=step;
	//Seg=Chain->GetSegment(linknum,segnum);
	//G=Seg->GetSWF();
	segnr=1; while (Chain->GetSegment(linknum,segnum) != Chain->GetDiffSegment(segnr)) segnr++;

	//Lat->PropagateG(Gi_inv,G,Gi_inv);
	Backward(M, segnr, Dgi_inv, Dg, Dgx, 0.0, Info);
	CudaComputePhiNode(nodeto,linknum,DensPart);

#endif
	return;
}

// This is where the propagators are computed
void
SF_Branched1stO::Matrix1(const DensityPart DensPart) {
	// variable declarations and space allocation
	int z;
	int i,j;
	Vector phi,G;
	Gi_inv.Dim(1,M); // storage space for inverse G
	fflush(stdout);
	#if DEBUG
	static int n_calls=0;
	Text t1=Blanks(3);
	printf("\n\nMATRIX1 CALLED %d times\n",++n_calls);
	fflush(stdout);
	#endif
	// Allocate space for Gzs matrix
	/* The matrix Gzs is stored as an Array of VecArrays (typedef for an Array of Vectors)
	There are no absolute segment ranking numbers!
	Each link has its own internal segment ranking from 1 to link.length()
	For each single link we create a Matrix Gzs_link that is internally represented as an
	Array of Vectors. The array of these matrices is what the full Gzs consists of
	To avoid problems creating an Array of Arrays, we define a type VecArray in the header.
	*/
	// the whole Gzs matrix
	Gzs.Dim(1,n_links); // contains n_links Vecarrays - the matrices for indovidual links
	 // initialize the matrices for all the links
	 for(i=1;i<=n_links;i++) {
		Gzs[i].Dim(1,Links[i].length());
		// for each link, initialize the matrix for all its segments
		for(j=1;j<=Links[i].length();j++) Gzs[i][j].Dim(1,M);
	}
	// end of variable declarations and space allocation

	// If necessary, convert all the G's and Phi's to safe mode
	for (i=1; i<=numDiffSegments; i++) {
		Seg = Chain->GetDiffSegment(i);
		G = Seg->GetSWF();
		Lat->MakeSafe(G);
		phi = Seg->GetPhi(DensPart);
		for (z=1; z<=M; z++) phi[z] = 0;
		Lat->MakeSafe(phi);
	}
	// The propagator starts here
	ComputeGzsNode(1,0);
	// After the propagator has finished, compute end-point distribution function for later use
	lnGN = Lat->ComputeLnGN(Gzs[1][1]);
	// Then compute the Phi
	ComputePhiNode(1,0,DensPart);
	// Correct for double counting
	for (i=1; i<=numDiffSegments; i++) {
		Seg = Chain->GetDiffSegment(i);
		G = Seg->GetSWF();
		phi = Seg->GetPhi(DensPart);
		Lat->CorrectDoubleCountG(phi,G);
	}
	// Different normalizations depending on freedom
	if (freedom == fixedTheta) {
		lnCt = log(theta) - lnGN - log(1.0*N);
		phiBulk = theta/exp(lnGN);
		SetPhiBulk(phiBulk);
	} else if (freedom == rangeRestricted) {
		double phiPreNorm = 0;
		for (i=1; i<=numDiffSegments; i++) {
			Seg = Chain->GetDiffSegment(i);
			phi = Seg->GetPhi(DensPart);
			for (z=1; z<=M; z++) {
				if (restrictedRange->InRange(z)) phiPreNorm += phi[z];
			}
		}
		lnCt = log(phiRange/phiPreNorm);
	} else lnCb = log(phiBulk/N);
	// renormalize the restricted molecule
	for (i=1; i<=numDiffSegments; i++) {
		Seg = Chain->GetDiffSegment(i);
		phi = Seg->GetPhi(DensPart);
		if (freedom == fixedTheta) {
			Lat->NormPhiRestr(phi,Gi_inv,theta/N);
		} else if (freedom == rangeRestricted) {
			Lat->NormPhiFree(phi,exp(lnCt));
		} else {
			Lat->NormPhiFree(phi,phiBulk/N);
		}
		G = Seg->GetSWF();
		Lat->RestoreFromSafe(G);
		Lat->RestoreFromSafe(phi);
	}
	//I fail to see why the following line is there.
	//if (freedom == rangeRestricted) Message(fatal,MyInput,"Programming error, rangRestricted freedom not implemented for branched copolymer yet");
}

void
SF_Branched1stO::CudaMatrix1(const DensityPart DensPart) {
#ifdef CUDA

	int z;
	int count;
	int i,k;

	Vector phi,G;
	Gi_inv.Dim(1,M);

	Gs=new int[n_links+1];
	Hgi_inv = &Gi_inv[1];  //pointers to host memory;
	Info = new int[11];	Lat->GetLatticeInfo(Info);

 // storage space for inverse G


	//Gzs.Dim(1,n_links);
	// for(i=1;i<=n_links;i++) {
	//	Gzs[i].Dim(1,Links[i].length());
	//	for(j=1;j<=Links[i].length();j++) Gzs[i][j].Dim(1,M);
	//}
	count=0;
	for(i=1;i<=n_links;i++) {
		Gs[i] = count;
		count = count + Links[i].length();
	}
	if (n!=count) Message(fatal,"program error in CudaMarix1" );
	// end of variable declarations and space allocation

	Dphi =    (double *)AllocateMemoryOnDevice(M*numDiffSegments);
	Dg =      (double *)AllocateMemoryOnDevice(M*numDiffSegments);
	Dgi =     (double *)AllocateMemoryOnDevice(M*n);
	Dgi_inv = (double *)AllocateMemoryOnDevice(M);
	Dgx =     (double *)AllocateMemoryOnDevice(M);
	Hphi =    (double *)malloc(sizeof(double)*M*numDiffSegments);
	Hg =      (double *)malloc(sizeof(double)*M*numDiffSegments);

	k=0;
	for (i=1; i<=numDiffSegments; i++) {
		Seg = Chain->GetDiffSegment(i);
		G = Seg->GetSWF();
		phi = Seg->GetPhi(DensPart);
		for (z=1; z<=M; z++) {
			phi[z] = 0;
			Hphi[k] = 0;
			Hg[k] =G[z];
			k++;
		}
	}

	TransferDataToDevice(M*numDiffSegments,Hphi, Dphi);
	TransferDataToDevice(M*numDiffSegments,Hg, Dg);

	// The propagator starts here
	CudaComputeGzsNode(1,0);
	// After the propagator has finished, compute end-point distribution function for later use


	//lnGN = Lat->ComputeLnGN(Gzs[1][1]);

        TransferDataToHost(M,Hgi_inv,Dgi);
	lnGN = Lat->ComputeLnGN(Gi_inv);

	// Then compute the Phi
	CudaComputePhiNode(1,0,DensPart);

	CorrectDoubleCounting(M*numDiffSegments,1,Dphi,Dg); //this is the way to do this for all phis
	TransferDataToHost(M*numDiffSegments,Hphi,Dphi);
	k=0;
	for (i=1; i<=numDiffSegments; i++) {
		Seg = Chain->GetDiffSegment(i);
		phi = Seg->GetPhi(DensPart);
		for (int z=1; z<=M; z++) {
			phi[z]=Hphi[k];
			k++;
		}
	}

	//// Correct for double counting
	//for (i=1; i<=numDiffSegments; i++) {
	//	Seg = Chain->GetDiffSegment(i);
	//	G = Seg->GetSWF();
	//	phi = Seg->GetPhi(DensPart);
	//	Lat->CorrectDoubleCountG(phi,G);
	//}

	// Different normalizations depending on freedom
	if (freedom == fixedTheta) {
		lnCt = log(theta) - lnGN - log(1.0*N);
		phiBulk = theta/exp(lnGN);
		SetPhiBulk(phiBulk);
	} else if (freedom == rangeRestricted) {
		double phiPreNorm = 0;
		for (i=1; i<=numDiffSegments; i++) {
			Seg = Chain->GetDiffSegment(i);
			phi = Seg->GetPhi(DensPart);
			for (z=1; z<=M; z++) {
				if (restrictedRange->InRange(z)) phiPreNorm += phi[z];
			}
		}
		lnCt = log(phiRange/phiPreNorm);
	} else lnCb = log(phiBulk/N);
	// renormalize the restricted molecule
	for (i=1; i<=numDiffSegments; i++) {
		Seg = Chain->GetDiffSegment(i);
		phi = Seg->GetPhi(DensPart);
		if (freedom == fixedTheta) {
			Lat->NormPhiRestr(phi,Gi_inv,theta/N);
		} else if (freedom == rangeRestricted) {
			Lat->NormPhiFree(phi,exp(lnCt));
		} else {
			Lat->NormPhiFree(phi,phiBulk/N);
		}
		G = Seg->GetSWF();
		Lat->RestoreFromSafe(G);
		Lat->RestoreFromSafe(phi);
	}

	FreeMemoryOnDevice(Dphi);
	FreeMemoryOnDevice(Dg);
	FreeMemoryOnDevice(Dgi);
	FreeMemoryOnDevice(Dgi_inv);
	FreeMemoryOnDevice(Dgx);
	free(Hphi); free(Hg);

	if (freedom == rangeRestricted) Message(fatal,MyInput,"Programming error, rangRestricted freedom not implemented for branched copolymer yet");
#else
	Message(fatal,MyInput,"Programming error: entered Cuda enabled routine but the compilation was not done accordingly. Compile SFBox sith Cuda=1");
#endif
}


void
SF_Branched1stO::Matrix1Long(const DensityPart DensPart) {
	Message(fatal,MyInput,"SF_Branched1stO::Matrix1Long not implemented yet.");
	return;
}

void SF_Branched1stO::Matrix2ndGen(const LatticeRange* Range, const DensityPart DensPart1, const DensityPart DensPart2) {
	Message(fatal,MyInput,"SF_Branched1stO::Matrix2ndGen not implemented yet");
	return;
}

void
SF_Branched1stO::Matrix2ndGenLong(const LatticeRange* Range, const DensityPart DensPart1, const DensityPart DensPart2) {
	Message(fatal,MyInput,"SF_Branched1stO::Matrix2ndGenLong not implemented yet");
	return;
}

//! This is identical to SF_Copol1stO but since no topology-related operations are used, it seems to be OK.
void
SF_Branched1stO::MatrixBulk(const DensityPart DensPart) {
	int i,z;
	Vector phiTot = GetPhi(DensPart);
	for (z=1; z<=M; z++) {
		if (Lat->BulkBoundary(z)) {
			phiTot[z] = 0;
		}
	}
	Seg = Chain->GetDiffSegment(1);
	Vector G = Seg->GetSWF();
	for (z=1; z<=M; z++) {
		if (Lat->BulkBoundary(z)) {
			phiTot[z] = pow(G[z],Chain->GetAvNumSegments(Seg));
		}
	}
	for (i=2; i<=numDiffSegments; i++) {
		Seg = Chain->GetDiffSegment(i);
		Vector G = Seg->GetSWF();
		for (z=1; z<=M; z++) {
			if (Lat->BulkBoundary(z)) {
				phiTot[z] *= pow(G[z],Chain->GetAvNumSegments(Seg));
			}
		}
	}
	for (z=1; z<=M; z++) {
		if (Lat->BulkBoundary(z)) {
			phiTot[z] *= phiBulk;
		}
	}
	for (i=1; i<=numDiffSegments; i++) {
		Seg = Chain->GetDiffSegment(i);
		Vector phi = Seg->GetPhi(DensPart);
		for (z=1; z<=M; z++) {
			if (Lat->BulkBoundary(z)) {
				phi[z] = phiTot[z]*Chain->GetAvNumSegments(Seg)/Chain->GetAvLength();
			}
		}
	}
}

//! This is identical to SF_Copol1stO but since no topology-related operations are used, it seems to be OK.
void
SF_Branched1stO::CopyBulkBoundaries(const DensityPart DensPart) {
	Message(fatal,MyInput,"SF_Branched1stO::CopyBulkBoundaries not implemented yet.");
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

#if DEBUG
//! Debug functions to print several vectors in columns
void SF_Branched1stO::PrintVec(const Vector V1, const Vector V2, const Vector V3, const Vector V4, const char *v1, const char *v2, const char *v3, const char *v4) {
	int dim=V1.Length()+1;
	int i;
	Lat->RestoreFromSafe(V1);
	Lat->RestoreFromSafe(V2);
	Lat->RestoreFromSafe(V3);
	Lat->RestoreFromSafe(V4);
	printf("        %9s:   %9s:   %9s:   %9s:   \n",v1,v2,v3,v4);
	for(i=1;i<dim;i++) printf("[%2d] = %.5e    %.5e    %.5e    %.5e\n",i,V1[i],V2[i],V3[i],V4[i]);
	printf("\n");
	Lat->MakeSafe(V1);
	Lat->MakeSafe(V2);
	Lat->MakeSafe(V3);
	Lat->MakeSafe(V4);
	fflush(stdout);
	return;
}
void SF_Branched1stO::PrintVec(const Vector V1, const Vector V2, const Vector V3, const char *v1, const char *v2, const char *v3) {
	int dim=V1.Length()+1;
	int i;
	Lat->RestoreFromSafe(V1);
	Lat->RestoreFromSafe(V2);
	Lat->RestoreFromSafe(V3);
	printf("        %9s:   %9s:   %9s:   \n",v1,v2,v3);
	for(i=1;i<dim;i++) printf("[%2d] = %.5e    %.5e    %.5e\n",i,V1[i],V2[i],V3[i]);
	printf("\n");
	Lat->MakeSafe(V1);
	Lat->MakeSafe(V2);
	Lat->MakeSafe(V3);
	fflush(stdout);
	return;
}
void SF_Branched1stO::PrintVec(const Vector V1, const Vector V2, const char *v1, const char *v2) {
	int dim=V1.Length()+1;
	int i;
	Lat->RestoreFromSafe(V1);
	Lat->RestoreFromSafe(V2);
	printf("        %9s:   %9s:   \n",v1,v2);
	for(i=1;i<dim;i++) printf("[%2d] = %.5e    %.5e\n",i,V1[i],V2[i]);
	printf("\n");
	Lat->MakeSafe(V1);
	Lat->MakeSafe(V2);
	fflush(stdout);
	return;
}
void SF_Branched1stO::PrintVec(const Vector V1, const char *v1) {
	int dim=V1.Length()+1;
	int i;
	Lat->RestoreFromSafe(V1);
	printf("        %9s:   \n",v1);
	for(i=1;i<dim;i++) printf("[%2d] = %.5e\n",i,V1[i]);
	printf("\n");
	Lat->MakeSafe(V1);
	fflush(stdout);
	return;
}
#endif
