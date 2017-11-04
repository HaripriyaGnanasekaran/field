#include "SF_Copol1stO.h"
#ifdef CUDA
#include "Cuda_tools.h"
#endif
#include <iostream>

SF_Copol1stO::SF_Copol1stO(Text name_,
						   SF_SegmentList* SegQ_,
						   Lat1stO* Lat_,
						   Input* MyInput_)
: SF_Molecule(name_, SegQ_, Lat_, MyInput_)
{
	if (Chain->GetMoleculeType() != copolymer) {
		Message(fatal,MyInput,
		"Programming error, trying to create a copolymer for molecule '" +
		name + "'");
	}
	symmetric = Chain->Symmetric();
	//if (MayerSaupeSet()) {
		//cout << "punt 1" << endl;
	//	Message(fatal,MyInput,"In SF_Copol1stO: Mayer_Saupe not implemented for first order propagators");
	//}
	force_set = (!GetForce()==0);
	if (force_set) {symmetric = false;
		Message(literal,MyInput,"Symmetry rejected in constant force ensemble.");
	}

#ifdef CUDA
	if (GetGPU()) GPUactivated = true; else GPUactivated = false;
#else
	GPUactivated = false;
	if (GetGPU()) {
		Message(literal,MyInput,"Compile with nvcc (option CUDA=1) to activate GPU... Contact Frans Leermakers. Going classical");
	}
#endif

	numDiffSegments = Chain->GetNumDiffSegments();
	Lat = Lat_;

	// allocate the Phi[z] vectors
	// descends through SF_Molecule->SF_MolSegemnt->SF_MolState
	if (freedom == secondGeneration) {
		CreatePhi(constrained);
		CreatePhi(unconstrained);
	} else {
		CreatePhi(total);
	}
	// calculate the size of the memory-saving matrix
	if (symmetric) {
		n = int(pow(N*3,1.0/3)+0.5);
		if (2*N < 120) n++;
		if (2*N < 60) n++;
	} else {
		n = int(pow(N*6,1.0/3)+0.5);
		if (N < 120) n++;
		if (N < 60) n++;
	}
	SetPhiBulk(phiBulk);
}
SF_Copol1stO::~SF_Copol1stO() {
	if (freedom == secondGeneration) {
		DeletePhi(constrained);
		DeletePhi(unconstrained);
	} else {
		DeletePhi(total);
	}
}

Vector
SF_Copol1stO::GetLong(){
	Message(fatal,"GetLong not implemented in SF_Copol1stO");
	Vector x;
	return x;
}
Vector
SF_Copol1stO::GetShort(){
	Message(fatal,"GetShort not implemented in SF_Copol1stO");
	Vector x;
	return x;
}
Vector
SF_Copol1stO::GetBondOrientation( const DensityPart DensPar) {
	Message(fatal,"GetBondOrientation not implemented in SF_Copol1stO");
	Vector x;
	return x; 
}; 
MoleculeType
SF_Copol1stO::GetMoleculeType() const {
	return Chain->GetMoleculeType();
}

void
SF_Copol1stO::ReComputePhi() {
	int z,i;
	//for (i=1; i<=numDiffSegments; i++) {
	//	Seg = Chain->GetDiffSegment(i);
	//		Vector phi = Seg->GetPhi(total);
	//		for (z=1; z<=M; z++) {
	//			phi[z]=0;
	//		}
	//}

	Matrix2(total);
	Vector phiTot = GetPhi(total);
	for (z=1; z<=M; z++) { phiTot[z] = 0;}
		for (i=1; i<=numDiffSegments; i++) {
			Seg = Chain->GetDiffSegment(i);
				Vector phi = Seg->GetPhi(total);
				for (z=1; z<=M; z++) {
					phiTot[z] += phi[z];
				}
		}
}

void
SF_Copol1stO::ComputePhi() {
	// allocate the Phi[z] vectors
	// descends to SF_Molecule->SF_MolSegemnt->SF_MolState
	if (Lat->GetNamesBulkBoundaries().Upperbound() > 0) {
		CreatePhi(bulk);
		MatrixBulk(bulk);
	}
	int i;
	for (i=1; i<=numDiffSegments; i++) {
		Seg = Chain->GetDiffSegment(i); // return i-th member of the SegmentQ of Chain - the i-th homo-block
		Vector G = Seg->GetSWF();
		Lat->SetBoundaries(G);
	}
	int z;
	// decide what kind of matrix to use for G
	if (symmetric) {
		if (freedom != secondGeneration) {
			if (!saveMemory) {
				if (GPUactivated) {
					CudaSymMatrix1(total);
				} else SymMatrix1(total);
			}
			else {
				if (GPUactivated) {
					CudaSymMatrix1Long(total);
				} else SymMatrix1Long(total);

			}
		} else {
			if (GPUactivated) {Message(literal,MyInput,"GPU-activation declined for Second Generation. Contact Frans Leermakers. Going classical");}
			if (!saveMemory) SymMatrix2ndGen(LatRange,constrained,unconstrained);
			else {
				SymMatrix2ndGenLong(LatRange,constrained,unconstrained);
			}
		}
	} else {
		if (freedom != secondGeneration) {
			if (!saveMemory) {
				if (GPUactivated) {
					CudaMatrix1(total);
				} else Matrix1(total);
			}
			else {
				if (GPUactivated) {
					CudaMatrix1Long(total);
				} else Matrix1Long(total);
			}
		} else {
			if (GPUactivated) {Message(literal,MyInput,"GPU-activation declined for Second Generation. Contact Frans Leermakers. Going classical");}
			if (!saveMemory) Matrix2ndGen(LatRange,constrained,unconstrained);
			else Matrix2ndGenLong(LatRange,constrained,unconstrained);
		}
	}
	if (Chain->SomeSegmentsPinned()) SetPhiBulk(0);
	if (Chain->SomeSegmentsGrafted()) SetPhiBulk(0);
	// Upperbound() is the size of the array (NOT Length()!!!)
	if (Lat->GetNamesBulkBoundaries().Upperbound() > 0) {
		if (freedom == thirdGeneration || freedom == secondGeneration) {
			CopyBulkBoundaries(unconstrained);
		} else {
			CopyBulkBoundaries(total);
		}
		DeletePhi(bulk);
	}
	if (freedom == secondGeneration) {
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
		if (muFixed) {
			Lat->SubtractBoundaries(phi1Tot);
			Lat->MultiplyWithLatticeSites(phi1Tot);
			theta = 0;
			for (int z=1; z<=M; z++) {
				theta += phi1Tot[z];
			}
			Lat->DivideByLatticeSites(phi1Tot);
			Lat->RestoreBoundaries(phi1Tot);
		}
	} else {
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
}
void
SF_Copol1stO::GetLayerAnalysis(Output* Out) {
	if (LatRangeStartLoops == NULL || LatRangeTrains == NULL) {
		return;
	}
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
		if (force_set) {
			Lat->Propagate2G(Gads,Gfree,G,s,LatRangeTrains,GetForce());
		} else {
		Lat->Propagate2G(Gads,Gfree,G,s,LatRangeTrains);
		}
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
		if (force_set) {
			Lat->Propagate2G(Gads_inv,Gfree_inv,G,LatRangeTrains,-GetForce());
		} else {
			Lat->Propagate2G(Gads_inv,Gfree_inv,G,LatRangeTrains);
		}
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
		if (force_set) {
			Lat->Propagate2G(Gads_inv2,Gfree_inv2,G,s,LatRangeTrains,GetForce());
		} else
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
				if (force_set) {
					Lat->PropagateG(Gt,Gused,s,GetForce());
				} else
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
				if (force_set) {
					Lat->PropagateG(Gt,Gused,s,GetForce());
				} else
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

		double numTrainsTot = 0;
		double numLoopsTot = 0;
		double numTailsTot = 0;
		for (s=1; s<=N; s++) {
			numTrainsTot += numTrains[s];
			numLoopsTot += numLoops[s];
			numTailsTot += numTails[s];
		}
		fluctFrTrains =fluctFrTrains/numTrainsTot - frTrains*frTrains/numTrainsTot/numTrainsTot;
		fluctFrLoops = fluctFrLoops/numLoopsTot - frLoops*frLoops/numLoopsTot/numLoopsTot;
		fluctFrTails = fluctFrTails/numTailsTot - frTails*frTails/numTailsTot/numTailsTot;
		double avLengthTrains = N*frTrains/numTrainsTot;
		double avLengthLoops = N*frLoops/numLoopsTot;
		double avLengthTails = N*frTails/numTailsTot;
		double fluctAvLengthTrains = N*N*fluctFrTrains;
		double fluctAvLengthLoops = N*N*fluctFrLoops;
		double fluctAvLengthTails = N*N*fluctFrTails;
		Out->PutReal("mol", name, "trains fract", frTrains);
		Out->PutReal("mol", name, "fluct trains fract", fluctFrTrains);
		Out->PutReal("mol", name, "trains number", numTrainsTot);
		Out->PutReal("mol", name, "trains av length", avLengthTrains);
		Out->PutReal("mol", name, "fluct trains av length", fluctAvLengthTrains);
		Out->PutVector("mol", name, "num trains", numTrains, N);
		Out->PutReal("mol", name, "loops fract", frLoops);
		Out->PutReal("mol", name, "fluct loops fract", fluctFrLoops);
		Out->PutReal("mol", name, "loops number", numLoopsTot);
		Out->PutReal("mol", name, "loops av length", avLengthLoops);
		Out->PutReal("mol", name, "fluct loops av length", fluctAvLengthLoops);
		Out->PutVector("mol", name, "num loops", numLoops, N);
		Out->PutReal("mol", name, "tails fract", frTails);
		Out->PutReal("mol", name, "fluct tails fract", fluctFrTails);
		Out->PutReal("mol", name, "tails number", numTailsTot);
		Out->PutReal("mol", name, "tails av length", avLengthTails);
		Out->PutReal("mol", name, "fluct tails av length", fluctAvLengthTails);
		Out->PutVector("mol", name, "num tails", numTails, N);
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
void
SF_Copol1stO::StateDistribution(Output* Out) const {
	// function to calculate the dissociation distribution function
	int z,s;
	int numSegWithTwoStates = 0;
	SF_MolSegment* Seg;
	for (s=1; s<=numDiffSegments; s++) {
 		Seg = Chain->GetDiffSegment(s);
 		if (Seg->GetNumStates() > 2) {
 			Message(literal, "Unable to generate multistate "
 			"output for segments with more than 2 states");
 			return;
 		} else if (Seg->GetNumStates() == 2) {
 			numSegWithTwoStates += (int) Chain->GetAvNumSegments(Seg);
 		}
	}
	Array<Vector> EndPoints(0,numSegWithTwoStates);
	Vector Gi,G;
	for (s=0; s<=numSegWithTwoStates; s++) {
		EndPoints[s].Dim(0,M);
		Lat->MakeSafe(EndPoints[s]);
	}
	for (s=1; s<=numDiffSegments; s++) {
 		Seg = Chain->GetDiffSegment(s);
 		if (Seg->GetNumStates() == 2) {
			G = Seg->GetState(1)->GetSWF();
			double alphaBulk = Seg->GetState(1)->GetAlphaBulk();
			for (z=1; z<=M; z++) {
				G[z] *= alphaBulk;
			}
			Lat->MakeSafe(G);
			G = Seg->GetState(2)->GetSWF();
			alphaBulk = Seg->GetState(2)->GetAlphaBulk();
			for (z=1; z<=M; z++) {
				G[z] *= alphaBulk;
			}
			Lat->MakeSafe(G);
 		} else {
 			G = Seg->GetSWF();
			Lat->MakeSafe(G);
 		}
	}
	int numSegWithTwoStatesCurrent = 0;
	Seg = Chain->GetSegment(1);
	if (Seg->GetNumStates() == 2) {
		numSegWithTwoStatesCurrent++;
		G = Seg->GetState(1)->GetSWF();
		Gi = EndPoints[0];
		for (z=1; z<=M; z++) {
			Gi[z] = G[z];
		}
		G = Seg->GetState(2)->GetSWF();
		Gi = EndPoints[1];
		for (z=1; z<=M; z++) {
			Gi[z] = G[z];
		}
	} else {
		G = Seg->GetSWF();
		Gi = EndPoints[0];
		for (z=1; z<=M; z++) {
			Gi[z] = G[z];
		}
	}
	/* GiDummy is needed to add two vectors later on */
	Vector GiDummy(1,M);
	for (z=1; z<=M; z++) {
		GiDummy[z] = 1;
	}
	Lat->MakeSafe(GiDummy);
	for (s=2; s<=N; s++) {
		Seg = Chain->GetSegment(s);
		if (Seg->GetNumStates() == 2) {
			for (int t=numSegWithTwoStatesCurrent; t>=0; t--) {
				Vector GiTemp(1,M);
				Gi = EndPoints[t];
				for (z=1; z<=M; z++) {
					GiTemp[z] = Gi[z];
				}
				G = Seg->GetState(1)->GetSWF();
				if (force_set) {
					Lat->PropagateG(Gi,G,GetForce());
				} else
				Lat->PropagateG(Gi,G);
				G = Seg->GetState(2)->GetSWF();
				if (force_set) {
					Lat->PropagateG(GiTemp,G,GetForce());
				} else
				Lat->PropagateG(GiTemp,G);
				/* Add GiTemp to EndPoints[t+1] */
				Lat->ConnectG(GiTemp,GiDummy,EndPoints[t+1]);
			}
			numSegWithTwoStatesCurrent++;
		} else {
			G = Seg->GetSWF();
			for (int t=numSegWithTwoStatesCurrent; t>=0; t--) {
				Gi = EndPoints[t];
				if (force_set) {
					Lat->PropagateG(Gi,G,GetForce());
				} else
				Lat->PropagateG(Gi,G);
			}
		}

	}
	double total1 = LOGZERO;
	for (s=0; s<=numSegWithTwoStates; s++) {
		Gi = EndPoints[s];
		Boolean compute = false;
		for (z=1; z<=M; z++) {
			if (Lat->OverflowProtection()) {
				if (Gi[z] > LOGZERO) {
					compute = true;
				}
			} else {
				if (Gi[z] > 0) {
					compute = true;
				}
			}
		}
		double Gtot = Lat->ComputeLnGN(Gi);
		if (Gtot>LOGZERO && compute) {
			double x = Gtot-total1;
			if (x > LOGMAXDOUBLE-1 || total1 == LOGZERO) {
				total1 = Gtot;
			} else {
				total1 += log1p(exp(x));
			}
		}
	}
	Vector Pm(0,numSegWithTwoStates);
	for (s=0; s<=numSegWithTwoStates; s++) {
		Gi = EndPoints[s];
		for (z=1; z<=M; z++) {
			if (Lat->OverflowProtection()) {
				if (Gi[z] > LOGZERO) {
					Gi[z] -= total1;
				}
			} else {
				Gi[z] /= exp(total1);
			}
		}
		Pm[s] = exp(Lat->ComputeLnGN(Gi));
		Lat->RestoreFromSafe(Gi);
		Text number=Blanks(100);
		number.Putint(s);
		number = Copy(number.Strip().Frontstrip());
		Out->PutProfile("mol",name,"EPst"+number,Gi);
	}
	Out->PutVector("mol",name,"P(states)",Pm,numSegWithTwoStates,0);
	double mAv = 0;
	double mFluct = 0;
	for (s=0; s<=numSegWithTwoStates; s++) {
		mAv += s*Pm[s];
 		mFluct += s*s*Pm[s];
	}
	Out->PutReal("mol",name,"states_av",mAv);
	Out->PutReal("mol",name,"states_fluct",mFluct);
	int numExtr = NumExtrema(Pm,0,numSegWithTwoStates);
	Vector place(1,numExtr);
	place = PlaceExtrema(Pm,0,numSegWithTwoStates,numExtr);
	Vector value(1,numExtr);
	value = ValueExtrema(Pm,0,numSegWithTwoStates,numExtr);
	for (s=1; s<=numExtr; s++) {
		Text number=Blanks(100);
		number.Putint(s);
		number = Copy(number.Strip().Frontstrip());
		Out->PutReal("mol",name,"placePstatesExt"+number,place[s]);
		Out->PutReal("mol",name,"valuePstatesExt"+number,value[s]);
	}
	for (s=1; s<=numDiffSegments; s++) {
 		Seg = Chain->GetDiffSegment(s);
 		if (Seg->GetNumStates() == 2) {
			G = Seg->GetState(1)->GetSWF();
			double alphaBulk = Seg->GetState(1)->GetAlphaBulk();
			Lat->RestoreFromSafe(G);
			for (z=1; z<=M; z++) {
				G[z] /= alphaBulk;
			}
			G = Seg->GetState(2)->GetSWF();
			alphaBulk = Seg->GetState(2)->GetAlphaBulk();
			Lat->RestoreFromSafe(G);
			for (z=1; z<=M; z++) {
				G[z] /= alphaBulk;
			}
 		} else {
 			G = Seg->GetSWF();
			Lat->RestoreFromSafe(G);
 		}
	}
}

void
SF_Copol1stO::ContactNumberDistr(Output* Out, Boolean partial) const {
	// function to calculate the endpoint distribution function as
	// a function of the number of surface contacts.
	if (partial) { // fast routine used to determine spinodals
 	 	Array<Vector> EndPoints(0,2);
		int z,s;
		Vector Gi,Gi2,G;
		for (s=0; s<=2; s++) {
			EndPoints[s].Dim(0,M);
			Lat->MakeSafe(EndPoints[s]);
		}
		for (s=1; s<=numDiffSegments; s++) {
			G = Chain->GetDiffSegment(s)->GetSWF();
			Lat->MakeSafe(G);
		}
		G = Chain->GetSegment(1)->GetSWF();
		Gi = EndPoints[0];
		Gi2 = EndPoints[1];
		Lat->Init2G(Gi2,Gi,G,LatRangeTrains);
		for (s=2; s<=N; s++) {
			G = Chain->GetSegment(s)->GetSWF();
			for (int t=2; t>=0; t--) {
				Gi = EndPoints[t];
				if (force_set) {
					Lat->PropagateG(Gi,G,GetForce());
				} else
				Lat->PropagateG(Gi,G);
				if (t < 2) {
					Gi2 = EndPoints[t+1];
				}
				for (z=1; z<=M; z++) {
					if (LatRangeTrains->InRange(z)) {
						if (t<2) {
							Gi2[z] = Gi[z];
						}
						if (t==0) {
							if (Lat->OverflowProtection()) {
								Gi[z] = LOGZERO;
							} else {
								Gi[z] = 0;
							}
						}
					}
				}
			}
		}
		double total1 = LOGZERO;
		for (s=0; s<=2; s++) {
			Gi = EndPoints[s];
			Boolean compute = false;
			for (z=1; z<=M; z++) {
				if (Lat->OverflowProtection()) {
					if (Gi[z] > LOGZERO) {
						compute = true;
					}
				} else {
					if (Gi[z] > 0) {
						compute = true;
					}
				}
			}
			double Gtot = Lat->ComputeLnGN(Gi);
			if (Gtot>LOGZERO && compute) {
				double x = Gtot-total1;
				if (x > LOGMAXDOUBLE-1 || total1 == LOGZERO) {
					total1 = Gtot;
				} else {
					total1 += log1p(exp(x));
				}
			}
		}
		Vector Pm(0,2);
		for (s=0; s<=2; s++) {
			Gi = EndPoints[s];
			for (z=1; z<=M; z++) {
				if (Lat->OverflowProtection()) {
					if (Gi[z] > LOGZERO) {
						Gi[z] -= total1;
					}
				} else {
					Gi[z] /= exp(total1);
				}
			}
			Pm[s] = exp(Lat->ComputeLnGN(Gi));
			Lat->RestoreFromSafe(Gi);
		}
		Out->PutReal("mol",name,"P(0)",Pm[0]);
		Out->PutReal("mol",name,"P(1)",Pm[1]);
		Out->PutReal("mol",name,"P(2)",Pm[2]);
		for (s=1; s<=numDiffSegments; s++) {
			G = Chain->GetDiffSegment(s)->GetSWF();
			Lat->RestoreFromSafe(G);
		}
		return;
	}
 	Array<Vector> EndPoints(0,N);
	int z,s;
	Vector Gi,Gi2,G;
	for (s=0; s<=N; s++) {
		EndPoints[s].Dim(0,M);
		Lat->MakeSafe(EndPoints[s]);
	}
	for (s=1; s<=numDiffSegments; s++) {
		G = Chain->GetDiffSegment(s)->GetSWF();
		Lat->MakeSafe(G);
	}
	G = Chain->GetSegment(1)->GetSWF();
	Gi = EndPoints[0];
	Gi2 = EndPoints[1];
	Lat->Init2G(Gi2,Gi,G,LatRangeTrains);
	for (s=2; s<=N; s++) {
		G = Chain->GetSegment(s)->GetSWF();
		for (int t=s-1; t>=0; t--) {
			Gi = EndPoints[t];
			if (force_set) {
				Lat->PropagateG(Gi,G,GetForce());
			} else
			Lat->PropagateG(Gi,G);
			Gi2 = EndPoints[t+1];
			for (z=1; z<=M; z++) {
				if (LatRangeTrains->InRange(z)) {
					Gi2[z] = Gi[z];
					if (t==0) {
						if (Lat->OverflowProtection()) {
							Gi[z] = LOGZERO;
						} else {
							Gi[z] = 0;
						}
					}
				}
			}
		}
	}
	double total1 = LOGZERO;
	for (s=0; s<=N; s++) {
		Gi = EndPoints[s];
		Boolean compute = false;
		for (z=1; z<=M; z++) {
			if (Lat->OverflowProtection()) {
				if (Gi[z] > LOGZERO) {
					compute = true;
				}
			} else {
				if (Gi[z] > 0) {
					compute = true;
				}
			}
		}
		double Gtot = Lat->ComputeLnGN(Gi);
		if (Gtot>LOGZERO && compute) {
			double x = Gtot-total1;
			if (x > LOGMAXDOUBLE-1 || total1 == LOGZERO) {
				total1 = Gtot;
			} else {
				total1 += log1p(exp(x));
			}
		}
	}
	Vector Pm(0,N);
	for (s=0; s<=N; s++) {
		Gi = EndPoints[s];
		for (z=1; z<=M; z++) {
			if (Lat->OverflowProtection()) {
				if (Gi[z] > LOGZERO) {
					Gi[z] -= total1;
				}
			} else {
				Gi[z] /= exp(total1);
			}
		}
		Pm[s] = exp(Lat->ComputeLnGN(Gi));
		Lat->RestoreFromSafe(Gi);
		Text number=Blanks(100);
		number.Putint(s);
		number = Copy(number.Strip().Frontstrip());
		Out->PutProfile("mol",name,"EPm="+number,Gi);
	}
	Out->PutVector("mol",name,"P(m)",Pm,N,0);
	double mAv = 0;
	double mFluct = 0;
	double OP_Av = 0;
	double OP_Fluct = 0;
	for (s=0; s<=N; s++) {
		mAv += s*Pm[s];
 		mFluct += s*s*Pm[s];
		double orderParam = (double(2*s-N))/double(N);
		OP_Av += orderParam*Pm[s];
 		OP_Fluct += orderParam*orderParam*Pm[s];
	}
	Out->PutReal("mol",name,"m_av",mAv);
	Out->PutReal("mol",name,"m_fluct",mFluct);
	Out->PutReal("mol",name,"avOrderParam",OP_Av);
	Out->PutReal("mol",name,"flOrderParam",OP_Fluct);
	int numExtr = NumExtrema(Pm,0,N);
	Vector place(1,numExtr);
	place = PlaceExtrema(Pm,0,N,numExtr);
	Vector value(1,numExtr);
	value = ValueExtrema(Pm,0,N,numExtr);
	for (s=1; s<=numExtr; s++) {
		Text number=Blanks(100);
		number.Putint(s);
		number = Copy(number.Strip().Frontstrip());
		Out->PutReal("mol",name,"placePmExt"+number,place[s]);
		Out->PutReal("mol",name,"valuePmExt"+number,value[s]);
	}
	Vector dPm(1,N);
	for (s=1; s<=N; s++) {
		dPm[s] = Pm[s] - Pm[s-1];
	}
	Out->PutVector("mol",name,"dP(m)",dPm,N,1);
	numExtr = NumExtrema(dPm,0,N);
	place.Dim(1,numExtr);
	place = PlaceExtrema(dPm,0,N,numExtr);
	value.Dim(1,numExtr);
	value = ValueExtrema(dPm,0,N,numExtr);
	for (s=1; s<=numExtr; s++) {
		Text number=Blanks(100);
		number.Putint(s);
		number = Copy(number.Strip().Frontstrip());
		Out->PutReal("mol",name,"placedPmExt"+number,place[s]);
		Out->PutReal("mol",name,"valuedPmExt"+number,value[s]);
	}
	Vector endpoint = Chain->GetSegment(N)->GetPhi(total);
	numExtr = NumExtrema(endpoint,1,M);
	place.Dim(1,numExtr);
	place = PlaceExtrema(endpoint,1,M,numExtr);
	value.Dim(1,numExtr);
	value = ValueExtrema(endpoint,1,M,numExtr);
	for (s=1; s<=numExtr; s++) {
		Text number=Blanks(100);
		number.Putint(s);
		number = Copy(number.Strip().Frontstrip());
		Out->PutReal("mol",name,"placeEPExt"+number,place[s]);
		Out->PutReal("mol",name,"valueEPExt"+number,value[s]);
	}
	Vector dendpoint(2,M);
	for (s=2; s<=M; s++) {
		dendpoint[s] = endpoint[s] - endpoint[s-1];
	}
	Out->PutVector("mol",name,"dEP",dendpoint,N,2);
	numExtr = NumExtrema(dendpoint,1,M);
	place.Dim(1,numExtr);
	place = PlaceExtrema(dendpoint,1,M,numExtr);
	value.Dim(1,numExtr);
	value = ValueExtrema(dendpoint,1,M,numExtr);
	for (s=1; s<=numExtr; s++) {
		Text number=Blanks(100);
		number.Putint(s);
		number = Copy(number.Strip().Frontstrip());
		Out->PutReal("mol",name,"placedEPExt"+number,place[s]);
		Out->PutReal("mol",name,"valuedEPExt"+number,value[s]);
	}
	for (s=1; s<=numDiffSegments; s++) {
		G = Chain->GetDiffSegment(s)->GetSWF();
		Lat->RestoreFromSafe(G);
	}
	/* old code to compute probability for walk towards the surface
	G = Chain->GetSegment(N)->GetSWF();
	Vector Gcopy(1,M);
	int zTrainMax;
	for (z=1; z<=M; z++) {
		if (LatRangeTrains->InRange(z)) {
			Gcopy[z] = G[z];
			zTrainMax = z;
		}
	}
	Matrix G3(1,M,1,N);
	G = Chain->GetSegment(1)->GetSWF();
	for (z=1; z<=M; z++) {
		G3[z][1] = G[z];
	}
	for (s=2; s<=N; s++) {
		if (force_set) {
			Lat->PropagateG(G3,Gcopy,s,GetForce());
		} else
		Lat->PropagateG(G3,Gcopy,s);
	}
	Vector Pzn(1,N);
	for (s=1; s<=N; s++) {
		Pzn[s] = G3[zTrainMax][s];
	}
	Out->PutVector("mol",name,"P(z,n)",Pzn,N,1);
	*/

}

void
SF_Copol1stO::Matrix2(const DensityPart DensPart) {
	int z,s,i;
	Vector phi,G;
	Matrix Gi(1,M,1,N); // here we need the full MxN matrix
	Vector Gi_inv(1,M);

	for (i=1; i<=numDiffSegments; i++) {
		Seg = Chain->GetDiffSegment(i); // Seg=first segment of i-th homogeneous block of the Chain - the i-th block in the SegBlockQ
		G = Seg->GetSWF();
		Lat->MakeSafe(G);
		phi = Seg->GetPhi(DensPart);
		for (z=1; z<=M; z++) {
			phi[z] = 0;
		}
		Lat->MakeSafe(phi);
	}
	// PropagateG from 1 to N
	G = Chain->GetSegment(1)->GetSWF();
	for (z=1; z<=M; z++) {
		Gi[z][1] = G[z];
	}
	for (s=2; s<=N; s++) {
		G = Chain->GetSegment(s)->GetSWF();
		if (force_set) {
			Lat->PropagateG(Gi,G,s,GetForce());
		} else
		Lat->PropagateG(Gi,G,s);
	}
	Seg = Chain->GetSegment(N);
	G = Seg->GetSWF();
	// and from N to 1, connecting G at the same time and forgetting the actual Gi_inv
	for (z=1; z<=M; z++) {
		Gi_inv[z] = G[z];
	}
	phi = Seg->GetPhi(DensPart);
	Lat->ConnectG(Gi_inv,Gi,N,phi);
	for (s=N-1; s>=1; s--) {
		Seg = Chain->GetSegment(s);
		G = Seg->GetSWF();
		phi = Seg->GetPhi(DensPart);
		if (force_set) {
			Lat->PropagateG(Gi_inv,G,-GetForce());
		} else
		Lat->PropagateG(Gi_inv,G);
		Lat->ConnectG(Gi_inv,Gi,s,phi);
	}
	// now we have Gs and can compute phi
	// lnGN = Lat->ComputeLnGN(Gi_inv);
	for (i=1; i<=numDiffSegments; i++) {
		Seg = Chain->GetDiffSegment(i);
		G = Seg->GetSWF();
		phi = Seg->GetPhi(DensPart);
		Lat->CorrectDoubleCountG(phi,G);
	}

	if (freedom == fixedTheta) {
		lnCt = log(theta) - lnGN - log(1.*N);

		phiBulk = theta/exp(lnGN);
		//SetPhiBulk(phiBulk);
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
			Lat->NormPhiFree(phi,exp(lnCt));
		} else if (freedom == rangeRestricted) {
			Lat->NormPhiFree(phi,exp(lnCt));
		} else {
			Lat->NormPhiFree(phi,phiBulk/N);

		}
		G = Seg->GetSWF();
		Lat->RestoreFromSafe(G);
		Lat->RestoreFromSafe(phi);
	}


}


// This is the classical (simple) version of Matrix1

void
SF_Copol1stO::Matrix1(const DensityPart DensPart) {
	int z,s,i;
	Vector phi,G;
	Matrix Gi(1,M,1,N); // here we need the full MxN matrix
	Vector Gi_inv(1,M);

	for (i=1; i<=numDiffSegments; i++) {
		Seg = Chain->GetDiffSegment(i); // Seg=first segment of i-th homogeneous block of the Chain - the i-th block in the SegBlockQ
		G = Seg->GetSWF();
		Lat->MakeSafe(G);
		phi = Seg->GetPhi(DensPart);
		for (z=1; z<=M; z++) {
			phi[z] = 0;
		}
		Lat->MakeSafe(phi);
	}
	// PropagateG from 1 to N
	G = Chain->GetSegment(1)->GetSWF();
	for (z=1; z<=M; z++) {
		Gi[z][1] = G[z];
	}
	for (s=2; s<=N; s++) {
		G = Chain->GetSegment(s)->GetSWF();
		if (force_set) {
			Lat->PropagateG(Gi,G,s,GetForce());
		} else
		Lat->PropagateG(Gi,G,s);
	}
	Seg = Chain->GetSegment(N);
	G = Seg->GetSWF();
	// and from N to 1, connecting G at the same time and forgetting the actual Gi_inv
	for (z=1; z<=M; z++) {
		Gi_inv[z] = G[z];
	}
	phi = Seg->GetPhi(DensPart);
	Lat->ConnectG(Gi_inv,Gi,N,phi);
	for (s=N-1; s>=1; s--) {
		Seg = Chain->GetSegment(s);
		G = Seg->GetSWF();
		phi = Seg->GetPhi(DensPart);
		if (force_set) {
			Lat->PropagateG(Gi_inv,G,-GetForce());
		} else
		Lat->PropagateG(Gi_inv,G);
		Lat->ConnectG(Gi_inv,Gi,s,phi);
	}
	// now we have Gs and can compute phi
	lnGN = Lat->ComputeLnGN(Gi_inv);
	for (i=1; i<=numDiffSegments; i++) {
		Seg = Chain->GetDiffSegment(i);
		G = Seg->GetSWF();
		phi = Seg->GetPhi(DensPart);
		Lat->CorrectDoubleCountG(phi,G);
	}
	// different normalizations depending on freedom
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
	if (freedom == rangeRestricted) {
		theta = ComputeTheta();
		phiBulk = theta/exp(lnGN);
		SetPhiBulk(phiBulk);
	}
}

void
SF_Copol1stO::Matrix1Long(const DensityPart DensPart) {
	int z,s,s0,t0,v0,t,rs1,i;
	Vector phi,G;
	Matrix Gi(1,M,1,n);
	Vector Gi_inv(1,M), Gs(1,M);
	for (i=1; i<=numDiffSegments; i++) {
		Seg = Chain->GetDiffSegment(i);
		G = Seg->GetSWF();
		Lat->MakeSafe(G);
		phi = Seg->GetPhi(DensPart);
		for (z=1; z<=M; z++) {
			phi[z] = 0;
		}
		Lat->MakeSafe(phi);
	}
	// GetSegment(i) returns the i-th segment of the 1st SegmentBlock in the SegBlockQ of Chain - there is only one block in SEgBlockQ for a linear polymer
	// get the first segment
	G = Chain->GetSegment(1)->GetSWF();
	for (z=1; z<=M; z++) {
		Gi[z][1] = Gs[z] = G[z];
	}
	t=1;
	v0=t0=s0 = 0;
	// get the rest of 2..N segments
	for (s=2; s<=N; s++) {
		t++;
		// n is the numbr of matrix elements
		if (t>n) {
			t0++;
			if (t0 == n) t0 = ++v0;
			t = t0 + 1;
			s0 = s - t0 - 1;
		}
		G = Chain->GetSegment(s)->GetSWF();
		if (force_set) {
			Lat->PropagateG(Gs,G,GetForce());
		} else
		Lat->PropagateG(Gs,G);
		if ((t == t0+1 && t0 == v0)
		   || (t == t0+1 && ((n-t0)*(n-t0+1) >= N-1-2*(s0+t0)))
		   || (2*(n-t+s) >= N-1))
			for (z=1; z<=M; z++) Gi[z][t] = Gs[z];
	}
	for (s=N; s>=1; s--) {
		G = Chain->GetSegment(s)->GetSWF();
		if (s == N) {
			for (z=1; z<=M; z++) {
				Gi_inv[z] = G[z];
			}
		} else {
			if (force_set) {
				Lat->PropagateG(Gi_inv,G,-GetForce());
			} else
			Lat->PropagateG(Gi_inv,G);
		}
		t = s - s0;
		if (t == t0) {
			s0 += - n + t0;
			if (t0 == v0 ) {
				s0 -= ((n - t0)*(n - t0 + 1))/2;
			}
			t0 --;
			if (t0 < v0) {
				v0 = t0;
			}
			for (z=1; z<=M; z++) {
				Gs[z] = Gi[z][t];
			}
			for (rs1=s0+t0+2; rs1<=s; rs1++) {
				t++;
				G = Chain->GetSegment(rs1)->GetSWF();
				if (force_set) {
					Lat->PropagateG(Gs,G,GetForce());
				} else
				Lat->PropagateG(Gs,G);
				if (t == t0+1 || s0+n == s) {
					for (z=1; z<=M; z++) {
						Gi[z][t] = Gs[z];
					}
				}
				if (t == n && s0+n < s) {
					t  = ++t0;
					s0 += n - t0;
				}
			}
			t = n;
		}
		Seg = Chain->GetSegment(s);
		phi = Seg->GetPhi(DensPart);
		Lat->ConnectG(Gi_inv,Gi,t,phi);
	}
	lnGN = Lat->ComputeLnGN(Gi_inv);
	for (i=1; i<=numDiffSegments; i++) {
		Seg = Chain->GetDiffSegment(i);
		G = Seg->GetSWF();
		phi = Seg->GetPhi(DensPart);
		Lat->CorrectDoubleCountG(phi,G);
	}
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
	if (freedom == rangeRestricted) {
		theta = ComputeTheta();
		phiBulk = theta/exp(lnGN);
		SetPhiBulk(phiBulk);
	}
}
void SF_Copol1stO::Matrix2ndGen(const LatticeRange* Range,
								const DensityPart DensPart1,
								const DensityPart DensPart2) {
	int z,s,i;
	Vector phi1,phi2,G;
	Matrix Gi1(1,M,1,N);
	Vector Gi_inv1(1,M);
	Matrix Gi2(1,M,1,N);
	Vector Gi_inv2(1,M);
	for (i=1; i<=numDiffSegments; i++) {
		Seg = Chain->GetDiffSegment(i);
		G = Seg->GetSWF();
		Lat->MakeSafe(G);
		phi1 = Seg->GetPhi(DensPart1);
		phi2 = Seg->GetPhi(DensPart2);
		for (z=1; z<=M; z++) {
			phi1[z] = 0;
			phi2[z] = 0;
		}
		Lat->MakeSafe(phi1);
		Lat->MakeSafe(phi2);
	}
	G = Chain->GetSegment(1)->GetSWF();
	Lat->Init2G(Gi_inv1,Gi_inv2,G,Range); // Gi_inv used here to get values
	for (z=1; z<=M; z++) {
		Gi1[z][1] = Gi_inv1[z];
		Gi2[z][1] = Gi_inv2[z];
	}
	for (s=2; s<=N; s++) {
		G = Chain->GetSegment(s)->GetSWF();
		if (force_set) {
			Lat->Propagate2G(Gi1,Gi2,G,s,Range,GetForce());
		} else
		Lat->Propagate2G(Gi1,Gi2,G,s,Range);
	}
	Seg = Chain->GetSegment(N);
	phi1 = Seg->GetPhi(DensPart1);
	phi2 = Seg->GetPhi(DensPart2);
 	G = Seg->GetSWF();
	Lat->Init2G(Gi_inv1,Gi_inv2,G,Range);
	Lat->ConnectG(Gi_inv1,Gi1,N,phi1);
	Lat->Connect2G(Gi_inv1,Gi1,N,Gi_inv2,Gi2,N,phi1);
	Lat->ConnectG(Gi_inv2,Gi2,N,phi2);
	for (s=N-1; s>=1; s--) {
		Seg = Chain->GetSegment(s);
		G = Seg->GetSWF();
		phi1 = Seg->GetPhi(DensPart1);
		phi2 = Seg->GetPhi(DensPart2);
		if (force_set) {
			Lat->Propagate2G(Gi_inv1,Gi_inv2,G,Range,-GetForce());
		} else
		Lat->Propagate2G(Gi_inv1,Gi_inv2,G,Range);
		Lat->ConnectG(Gi_inv1,Gi1,s,phi1); //'loops'
		Lat->Connect2G(Gi_inv1,Gi1,s,Gi_inv2,Gi2,s,phi1); //'tails'
		Lat->ConnectG(Gi_inv2,Gi2,s,phi2); //'free'
	}
	lnGN = Lat->ComputeLnGN(Gi_inv1);
	if (!muFixed) {
		lnCt = log(theta) - lnGN - log(1.*N);
	}
	if (phiBulk != 0) {
		lnCb = log(phiBulk/N);
	}
	for (i=1; i<=numDiffSegments; i++) {
		Seg = Chain->GetDiffSegment(i);
		G = Seg->GetSWF();
		phi1 = Seg->GetPhi(DensPart1);
		phi2 = Seg->GetPhi(DensPart2);
		Lat->CorrectDoubleCountG(phi1,G);
		Lat->CorrectDoubleCountG(phi2,G);
		if (muFixed) {
			Lat->NormPhiFree(phi1,exp(lnCt));
		} else {
			Lat->NormPhiRestr(phi1,Gi_inv1,theta/N);
		}
		Lat->NormPhiFree(phi2,phiBulk/N);
		Lat->RestoreFromSafe(G);
		Lat->RestoreFromSafe(phi1);
		Lat->RestoreFromSafe(phi2);
	}
}
void
SF_Copol1stO::Matrix2ndGenLong(const LatticeRange* Range,
							   const DensityPart DensPart1,
							   const DensityPart DensPart2) {
	int z,s,s0,t0,v0,t,rs1,i;
	Vector phi1,phi2,G;
	Matrix Gi1(1,M,1,N);
	Vector Gi_inv1(1,M), Gs1(1,M);
	Matrix Gi2(1,M,1,N);
	Vector Gi_inv2(1,M), Gs2(1,M);
	for (i=1; i<=numDiffSegments; i++) {
		Seg = Chain->GetDiffSegment(i);
		G = Seg->GetSWF();
		Lat->MakeSafe(G);
		phi1 = Seg->GetPhi(DensPart1);
		phi2 = Seg->GetPhi(DensPart2);
		for (z=1; z<=M; z++) {
			phi1[z] = 0;
			phi2[z] = 0;
		}
		Lat->MakeSafe(phi1);
		Lat->MakeSafe(phi2);
	}
	G = Chain->GetSegment(1)->GetSWF();
	Lat->Init2G(Gs1,Gs2,G,Range);
	for (z=1; z<=M; z++) {
		Gi1[z][1] = Gs1[z];
		Gi2[z][1] = Gs2[z];
	}
	t=1;
	v0=t0=s0 = 0;
	for (s=2; s<=N; s++) {
		t++;
		if (t>n) {
			t0++;
			if (t0 == n) t0 = ++v0;
			t = t0 + 1;
			s0 = s - t0 - 1;
		}
		G = Chain->GetSegment(s)->GetSWF();
		if (force_set) {
			Lat->Propagate2G(Gs1,Gs2,G,Range,GetForce());
		} else
		Lat->Propagate2G(Gs1,Gs2,G,Range);
		if ((t == t0+1 && t0 == v0)
		   || (t == t0+1 && ((n-t0)*(n-t0+1) >= N-1-2*(s0+t0)))
		   || (2*(n-t+s) >= N-1)) {
			for (z=1; z<=M; z++) {
				Gi1[z][t] = Gs1[z];
				Gi2[z][t] = Gs2[z];
			}
		}
	}
	for (s=N; s>=1; s--) {
		G = Chain->GetSegment(s)->GetSWF();
		if (s == N) {
			Lat->Init2G(Gi_inv1,Gi_inv2,G,Range);
		} else {
			if (force_set) {
				Lat->Propagate2G(Gi_inv1,Gi_inv2,G,Range,-GetForce());
			} else
			Lat->Propagate2G(Gi_inv1,Gi_inv2,G,Range);
		}
		t = s - s0;
		if (t == t0) {
			s0 += - n + t0;
			if (t0 == v0 ) {
				s0 -= ((n - t0)*(n - t0 + 1))/2;
			}
			t0 --;
			if (t0 < v0) {
				v0 = t0;
			}
			for (z=1; z<=M; z++) {
				Gs1[z] = Gi1[z][t];
				Gs2[z] = Gi2[z][t];
			}
			for (rs1=s0+t0+2; rs1<=s; rs1++) {
				t++;
				G = Chain->GetSegment(rs1)->GetSWF();
				if (force_set) {
					Lat->Propagate2G(Gs1,Gs2,G,Range,-GetForce());
				} else
				Lat->Propagate2G(Gs1,Gs2,G,Range);
				if (t == t0+1 || s0+n == s) {
					for (z=1; z<=M; z++) {
						Gi1[z][t] = Gs1[z];
						Gi2[z][t] = Gs2[z];
					}
				}
				if (t == n && s0+n < s) {
					t  = ++t0;
					s0 += n - t0;
				}
			}
			t = n;
		}
		Seg = Chain->GetSegment(s);
		phi1 = Seg->GetPhi(DensPart1);
		phi2 = Seg->GetPhi(DensPart2);
		Lat->ConnectG(Gi_inv1,Gi1,t,phi1); //'loops'
		Lat->Connect2G(Gi_inv1,Gi1,t,Gi_inv2,Gi2,t,phi1); //'tails'
		Lat->ConnectG(Gi_inv2,Gi2,t,phi2); //'free'
	}
	lnGN = Lat->ComputeLnGN(Gi_inv1);
	if (!muFixed) {
		lnCt = log(theta) - lnGN - log(1.*N);
	}
	lnCb = log(phiBulk/N);
	for (i=1; i<=numDiffSegments; i++) {
		Seg = Chain->GetDiffSegment(i);
		G = Seg->GetSWF();
		phi1 = Seg->GetPhi(DensPart1);
		phi2 = Seg->GetPhi(DensPart2);
		Lat->CorrectDoubleCountG(phi1,G);
		Lat->CorrectDoubleCountG(phi2,G);
		if (muFixed) {
			Lat->NormPhiFree(phi1,exp(lnCt));
		} else {
			Lat->NormPhiRestr(phi1,Gi_inv1,theta/N);
		}
		Lat->NormPhiFree(phi2,phiBulk/N);
		Lat->RestoreFromSafe(G);
		Lat->RestoreFromSafe(phi1);
		Lat->RestoreFromSafe(phi2);
	}
}

void
SF_Copol1stO::SymMatrix1(const DensityPart DensPart) {
	int z,s,s_inv,i;
	Vector phi,G;
	Matrix Gi(1,M,1,N/2);
	Vector Gi_inv(1,M);
	for (i=1; i<=numDiffSegments; i++) {
		Seg = Chain->GetDiffSegment(i);
		G = Seg->GetSWF();
		Lat->MakeSafe(G);
		phi = Seg->GetPhi(DensPart);
		for (z=1; z<=M; z++) {
			phi[z] = 0;
		}
		Lat->MakeSafe(phi);
	}
	G = Chain->GetSegment(1)->GetSWF();
	for (z=1; z<=M; z++) {
		Gi[z][1] = G[z];
	}
	for (s=2; s<=N/2; s++) {
		G = Chain->GetSegment(s)->GetSWF();
		if (force_set) {
			Lat->PropagateG(Gi,G,s,GetForce());
		} else
		Lat->PropagateG(Gi,G,s);
	}
	s = N/2;
	for (z=1; z<=M; z++) {
		Gi_inv[z] = Gi[z][s];
	}
	if (N%2 == 1) {
		s = N/2 + 1;
		Seg = Chain->GetSegment(s);
		G = Seg->GetSWF();
		phi = Seg->GetPhi(DensPart);
		if (force_set) {
			Lat->PropagateG(Gi_inv,G,GetForce());
		} else	Lat->PropagateG(Gi_inv,G);
		Lat->ConnectG(Gi_inv,Gi_inv,phi);
		Lat->NormPhiFree(phi,0.5); // correction needed for optimization
	}
	for (s=(N+3)/2; s<=N; s++) {
		s_inv = N-s+1;
		Seg = Chain->GetSegment(s);
		G = Seg->GetSWF();
		if (force_set) {
			Lat->PropagateG(Gi_inv,G,GetForce());
		} else Lat->PropagateG(Gi_inv,G);
		phi = Seg->GetPhi(total);
		Lat->ConnectG(Gi_inv,Gi,s_inv,phi); // optimization, chain is symmetric
			// Connect only half the chain into phi
	}
	lnGN = Lat->ComputeLnGN(Gi_inv);

	for (i=1; i<=numDiffSegments; i++) {
		Seg = Chain->GetDiffSegment(i);
		G = Seg->GetSWF();
		phi = Seg->GetPhi(DensPart);
		Lat->CorrectDoubleCountG(phi,G);
	}
	if (freedom == fixedTheta) {
		lnCt = log(theta) - lnGN -log(1.*N);
		phiBulk = theta/exp(lnGN);
		Chain->SetPhiBulk(phiBulk);
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
		lnCt = log(phiRange/(2*phiPreNorm)); // factor 2 is correction needed for optimization
	} else {
		lnCb = log(phiBulk/N);
	}
	for (i=1; i<=numDiffSegments; i++) {
		Seg = Chain->GetDiffSegment(i);
		phi = Seg->GetPhi(DensPart);
		if (freedom == fixedTheta) {
			Lat->NormPhiRestr(phi,Gi_inv,2*theta/N); // factor 2 is correction needed for optimization
		} else if (freedom == rangeRestricted) {
			Lat->NormPhiFree(phi,2*exp(lnCt));// factor 2 is correction needed for optimization
		} else {
			Lat->NormPhiFree(phi,2*phiBulk/N);// factor 2 is correction needed for optimization
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
SF_Copol1stO::SymMatrix1Long(const DensityPart DensPart) {
	int z,s;
	//int s_inv;
	int s0,t0,v0,t,rs1,i;
	Vector phi,G;
	Matrix Gi(1,M,1,n);
	Vector Gi_inv(1,M), Gs(1,M);
	for (i=1; i<=numDiffSegments; i++) {
		Seg = Chain->GetDiffSegment(i);
		G = Seg->GetSWF();
		Lat->MakeSafe(G);
		phi = Seg->GetPhi(DensPart);
		for (z=1; z<=M; z++) {
			phi[z] = 0;
		}
		Lat->MakeSafe(phi);
	}
	G = Chain->GetSegment(1)->GetSWF();
	for (z=1; z<=M; z++) {
		Gi[z][1] = Gs[z] = G[z];
	}
	t=1;
	v0=t0=s0 = 0;
	for (s=2; s<=N/2; s++) {
		t++;
		if (t>n) {
			t0++;
			if (t0 == n) t0 = ++v0;
			t = t0 + 1;
			s0 = s - t0 - 1;
		}
		G = Chain->GetSegment(s)->GetSWF();
		if (force_set) {
			Lat->PropagateG(Gs,G,-GetForce());
		} else
		Lat->PropagateG(Gs,G);
		if ((t == t0+1 && t0 == v0)
		   || (t == t0+1 && ((n-t0)*(n-t0+1) >= N-1-2*(s0+t0)))
		   || (2*(n-t+s) >= N-1)) {
			for (z=1; z<=M; z++) {
				Gi[z][t] = Gs[z];
			}
		}
	}
	for (z=1; z<=M; z++) {
		Gi_inv[z] = Gs[z];
	}
	if (N%2 == 1) {
		s = N/2 + 1;
		Seg = Chain->GetSegment(s);
		G = Seg->GetSWF();
		if (force_set) {
			Lat->PropagateG(Gi_inv,G,GetForce());
		} else	Lat->PropagateG(Gi_inv,G);
		phi = Seg->GetPhi(DensPart);
		Lat->ConnectG(Gi_inv,Gi_inv,phi);
		Lat->NormPhiFree(phi,0.5);
	}
	for (s=(N+3)/2; s<=N; s++) {
		G = Chain->GetSegment(s)->GetSWF();
		if (force_set) {
			Lat->PropagateG(Gi_inv,G,GetForce());
		} else	Lat->PropagateG(Gi_inv,G);
		//s_inv = N-s+1;
		t = N - s + 1 - s0;
		if (t == t0) {
			s0 += - n + t0;
			if (t0 == v0)
				s0 -= ((n - t0)*(n - t0 + 1))/2;
			t0 --;
			if (t0 < v0)
				v0 = t0;
			for (z=1; z<=M; z++) {
				Gs[z] = Gi[z][t];
			}
			for (rs1=s0+t0+2; rs1<=(N-s+1); rs1++) {
				t++;
				G = Chain->GetSegment(rs1)->GetSWF();
				if (force_set) {
					Lat->PropagateG(Gs,G,GetForce());
				} else	Lat->PropagateG(Gs,G);
				if (t == t0+1 || s0+n == N-s+1) {
					for (z=1; z<=M; z++) {
						Gi[z][t] = Gs[z];
					}
				}
				if (t == n && s0+n < N-s+1) {
					t  = ++t0;
					s0 += n - t0;
				}
			}
			t = n;
		}
		Seg = Chain->GetSegment(s);
		phi = Seg->GetPhi(DensPart);
		Lat->ConnectG(Gi_inv,Gi,t,phi);

	}
	lnGN = Lat->ComputeLnGN(Gi_inv);
	for (i=1; i<=numDiffSegments; i++) {
		Seg = Chain->GetDiffSegment(i);
		G = Seg->GetSWF();
		phi = Seg->GetPhi(DensPart);
		Lat->CorrectDoubleCountG(phi,G);
	}
	if (freedom == fixedTheta) {
		lnCt = log(theta) - lnGN -log(1.*N);
		phiBulk = theta/exp(lnGN);
		Chain->SetPhiBulk(phiBulk);
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
		lnCt = log(phiRange/(2*phiPreNorm)); // factor 2 is correction needed for optimization
	} else {
		lnCb = log(phiBulk/N);
	}
	for (i=1; i<=numDiffSegments; i++) {
		Seg = Chain->GetDiffSegment(i);
		phi = Seg->GetPhi(DensPart);
		if (freedom == fixedTheta) {
			Lat->NormPhiRestr(phi,Gi_inv,2*theta/N); // factor 2 is correction needed for optimization
		} else if (freedom == rangeRestricted) {
			Lat->NormPhiFree(phi,2*exp(lnCt));// factor 2 is correction needed for optimization
		} else {
			Lat->NormPhiFree(phi,2*phiBulk/N);// factor 2 is correction needed for optimization
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
void SF_Copol1stO::SymMatrix2ndGen(const LatticeRange* Range,
								   const DensityPart DensPart1,
								   const DensityPart DensPart2) {
	int z,s,s_inv,i;
	Vector phi1, phi2, G;
	Matrix Gi1(1,M,1,N/2);
	Vector Gi_inv1(1,M);
	Matrix Gi2(1,M,1,N/2);
	Vector Gi_inv2(1,M);
	for (i=1; i<=numDiffSegments; i++) {
		Seg = Chain->GetDiffSegment(i);
		G = Seg->GetSWF();
		Lat->MakeSafe(G);
		phi1 = Seg->GetPhi(DensPart1);
		phi2 = Seg->GetPhi(DensPart2);
		for (z=1; z<=M; z++) {
			phi1[z] = 0;
			phi2[z] = 0;
		}
		Lat->MakeSafe(phi1);
		Lat->MakeSafe(phi2);
	}
	G = Chain->GetSegment(1)->GetSWF();
	Lat->Init2G(Gi_inv1,Gi_inv2,G,Range); // Gi_inv used here to get values
	for (z=1; z<=M; z++) {
		Gi1[z][1] = Gi_inv1[z];
		Gi2[z][1] = Gi_inv2[z];
	}
	for (s=2; s<=N/2; s++){
		G = Chain->GetSegment(s)->GetSWF();
		if (force_set) {
			Lat->Propagate2G(Gi1,Gi2,G,s,Range,GetForce());
		} else
		Lat->Propagate2G(Gi1,Gi2,G,s,Range);
	}
	s = N/2;
	for (z=1; z<=M; z++) {
		Gi_inv1[z] = Gi1[z][s];
		Gi_inv2[z] = Gi2[z][s];
	}
	if (N%2 == 1) {
		s = N/2 + 1;
		Seg = Chain->GetSegment(s);
		G = Seg->GetSWF();
		phi1 = Seg->GetPhi(DensPart1);
		phi2 = Seg->GetPhi(DensPart2);
		if (force_set) {
			Lat->Propagate2G(Gi_inv1,Gi_inv2,G,Range,GetForce());
		} else
		Lat->Propagate2G(Gi_inv1,Gi_inv2,G,Range);
		Lat->ConnectG(Gi_inv1,Gi_inv1,phi1); //'loops'
		Lat->Connect2G(Gi_inv1,Gi_inv1,Gi_inv2,Gi_inv2,phi1); //'tails'
		Lat->ConnectG(Gi_inv2,Gi_inv2,phi2); //'free'
		Lat->NormPhiFree(phi1,0.5);
		Lat->NormPhiFree(phi2,0.5);
	}
	for (s=(N+3)/2; s<=N; s++) {
		s_inv = N-s+1;
		Seg = Chain->GetSegment(s);
		G = Seg->GetSWF();
		phi1 = Seg->GetPhi(DensPart1);
		phi2 = Seg->GetPhi(DensPart2);
		if (force_set) {
			Lat->Propagate2G(Gi_inv1,Gi_inv2,G,Range,GetForce());
		} else
		Lat->Propagate2G(Gi_inv1,Gi_inv2,G,Range);
		Lat->ConnectG(Gi_inv1,Gi1,s_inv,phi1); //'loops'
		Lat->Connect2G(Gi_inv1,Gi1,s_inv,Gi_inv2,Gi2,s_inv,phi1); //'tails'
		Lat->ConnectG(Gi_inv2,Gi2,s_inv,phi2); //'free'
	}
	lnGN = Lat->ComputeLnGN(Gi_inv1);
	if (!muFixed) {
		lnCt = log(theta) - lnGN - log(1.*N);
	}
	lnCb = log(phiBulk/N);
	for (i=1; i<=numDiffSegments; i++) {
		Seg = Chain->GetDiffSegment(i);
		G = Seg->GetSWF();
		phi1 = Seg->GetPhi(DensPart1);
		phi2 = Seg->GetPhi(DensPart2);
		Lat->CorrectDoubleCountG(phi1,G);
		Lat->CorrectDoubleCountG(phi2,G);
		if (muFixed) {
			Lat->NormPhiFree(phi1,2*exp(lnCt));
		} else {
			Lat->NormPhiRestr(phi1,Gi_inv1,2*theta/N);
		}
		Lat->NormPhiFree(phi2,2*phiBulk/N);
		Lat->RestoreFromSafe(G);
		Lat->RestoreFromSafe(phi1);
		Lat->RestoreFromSafe(phi2);
	}
}
void
SF_Copol1stO::SymMatrix2ndGenLong(const LatticeRange* Range,
								  const DensityPart DensPart1,
								  const DensityPart DensPart2) {
	int z,s;
	//int s_inv;
	int s0,t0,v0,t,rs1,i;
	Vector phi1,phi2,G;
	Matrix Gi1(1,M,1,n);
	Vector Gi_inv1(1,M), Gs1(1,M);
	Matrix Gi2(1,M,1,n);
	Vector Gi_inv2(1,M), Gs2(1,M);
	for (i=1; i<=numDiffSegments; i++) {
		Seg = Chain->GetDiffSegment(i);
		G = Seg->GetSWF();
		Lat->MakeSafe(G);
		phi1 = Seg->GetPhi(DensPart1);
		phi2 = Seg->GetPhi(DensPart2);
		for (z=1; z<=M; z++) {
			phi1[z] = 0;
			phi2[z] = 0;
		}
		Lat->MakeSafe(phi1);
		Lat->MakeSafe(phi2);
	}
	G = Chain->GetSegment(1)->GetSWF();
	Lat->Init2G(Gs1,Gs2,G,Range);
	for (z=1; z<=M; z++) {
		Gi1[z][1] = Gs1[z];
		Gi2[z][1] = Gs2[z];
	}
	t=1;
	v0=t0=s0 = 0;
	for (s=2; s<=N/2; s++) {
		t++;
		if (t>n) {
			t0++;
			if (t0 == n) t0 = ++v0;
			t = t0 + 1;
			s0 = s - t0 - 1;
		}
		G = Chain->GetSegment(s)->GetSWF();
		if (force_set) {
			Lat->Propagate2G(Gs1,Gs2,G,Range,GetForce());
		} else
		Lat->Propagate2G(Gs1,Gs2,G,Range);
		if ((t == t0+1 && t0 == v0)
		   || (t == t0+1 && ((n-t0)*(n-t0+1) >= N-1-2*(s0+t0)))
		   || (2*(n-t+s) >= N-1)) {
			for (z=1; z<=M; z++) {
				Gi1[z][t] = Gs1[z];
				Gi2[z][t] = Gs2[z];
			}
		}
	}
	for (z=1; z<=M; z++) {
		Gi_inv1[z] = Gs1[z];
		Gi_inv2[z] = Gs2[z];
	}
	if (N%2 == 1) {
		s = N/2 + 1;
		Seg = Chain->GetSegment(s);
		G = Seg->GetSWF();
		if (force_set) {
			Lat->Propagate2G(Gi_inv1,Gi_inv2,G,Range,-GetForce());
		} else
		Lat->Propagate2G(Gi_inv1,Gi_inv2,G,Range);
		phi1 = Seg->GetPhi(DensPart1);
		phi2 = Seg->GetPhi(DensPart2);
		Lat->ConnectG(Gi_inv1,Gi_inv1,phi1); //'loops'
		Lat->Connect2G(Gi_inv1,Gi_inv1,Gi_inv2,Gi_inv2,phi1); //'tails'
		Lat->ConnectG(Gi_inv2,Gi_inv2,phi2); //'free'
		Lat->NormPhiFree(phi1,0.5);
		Lat->NormPhiFree(phi2,0.5);
	}
	for (s=(N+3)/2; s<=N; s++) {
		G = Chain->GetSegment(s)->GetSWF();
		if (force_set) {
			Lat->Propagate2G(Gi_inv1,Gi_inv2,G,Range,-GetForce());
		} else
		Lat->Propagate2G(Gi_inv1,Gi_inv2,G,Range);
		//s_inv = N-s+1;
		t = N - s + 1 - s0;
		if (t == t0) {
			s0 += - n + t0;
			if (t0 == v0)
				s0 -= ((n - t0)*(n - t0 + 1))/2;
			t0 --;
			if (t0 < v0)
				v0 = t0;
			for (z=1; z<=M; z++) {
				Gs1[z] = Gi1[z][t];
				Gs2[z] = Gi2[z][t];
			}
			for (rs1=s0+t0+2; rs1<=(N-s+1); rs1++) {
				t++;
				G = Chain->GetSegment(rs1)->GetSWF();
				if (force_set) {
					Lat->Propagate2G(Gs1,Gs2,G,Range,-GetForce());
				} else
				Lat->Propagate2G(Gs1,Gs2,G,Range);
				if (t == t0+1 || s0+n == N-s+1) {
					for (z=1; z<=M; z++) {
						Gi1[z][t] = Gs1[z];
						Gi2[z][t] = Gs2[z];
					}
				}
				if (t == n && s0+n < N-s+1) {
					t  = ++t0;
					s0 += n - t0;
				}
			}
			t = n;
		}
		Seg = Chain->GetSegment(s);
		phi1 = Seg->GetPhi(DensPart1);
		phi2 = Seg->GetPhi(DensPart2);
		Lat->ConnectG(Gi_inv1,Gi1,t,phi1); //'loops'
		Lat->Connect2G(Gi_inv1,Gi1,t,Gi_inv2,Gi2,t,phi1); //'tails'
		Lat->ConnectG(Gi_inv2,Gi2,t,phi2); //'free'
	}
	lnGN = Lat->ComputeLnGN(Gi_inv1);
	if (!muFixed) {
		lnCt = log(theta) - lnGN - log(1.*N);
	}
	lnCb = log(phiBulk/N);
	for (i=1; i<=numDiffSegments; i++) {
		Seg = Chain->GetDiffSegment(i);
		G = Seg->GetSWF();
		phi1 = Seg->GetPhi(DensPart1);
		phi2 = Seg->GetPhi(DensPart2);
		Lat->CorrectDoubleCountG(phi1,G);
		Lat->CorrectDoubleCountG(phi2,G);
		if (muFixed) {
			Lat->NormPhiFree(phi1,2*exp(lnCt));
		} else {
			Lat->NormPhiRestr(phi1,Gi_inv1,2*theta/N);
		}
		Lat->NormPhiFree(phi2,2*phiBulk/N);
		Lat->RestoreFromSafe(G);
		Lat->RestoreFromSafe(phi1);
		Lat->RestoreFromSafe(phi2);
	}
}
void
SF_Copol1stO::MatrixBulk(const DensityPart DensPart) {
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
void
SF_Copol1stO::CopyBulkBoundaries(const DensityPart DensPart) {
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

// Here are the propagators that are using the GPU.

void
SF_Copol1stO::CudaMatrix1(const DensityPart DensPart) {
#ifdef CUDA
	int Info[11];
	int z,s,i,k;
	int segnr=0;
	Vector phi,G;
	Lat->GetLatticeInfo(Info);
	//cout << "lattice dim " << Info[0] << endl;
	//cout << "lattice type " << Info[1] << endl;
	//cout << "lattice Mx " << Info[2] << endl;
	//cout << "lattice My " << Info[3] << endl;
	//cout << "lattice Mz " << Info[4] << endl;
	//cout << "lattice bx1 " << Info[5] << endl;
	//cout << "lattice bxm " << Info[6] << endl;
	//cout << "lattice by1 " << Info[7] << endl;
	//cout << "lattice bym " << Info[8] << endl;
	//cout << "lattice bz1 " << Info[9] << endl;
	//cout << "lattice bzm " << Info[10] << endl;

	Vector Gi_inv(1,M);

	double *Hphi, *Hg;
	double *Hgi_inv=&Gi_inv[1];  //pointers to host memory;
	double *Dphi, *Dg, *Dgi, *Dgi_inv, *Dgx; //pointers to device memory;

	Dphi =    (double *)AllocateMemoryOnDevice(M*numDiffSegments);
	Dg =      (double *)AllocateMemoryOnDevice(M*numDiffSegments);
	Dgi =     (double *)AllocateMemoryOnDevice(M*N);
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

	segnr=1; while (Chain->GetSegment(1) != Chain->GetDiffSegment(segnr)) segnr++;

	InitializeForward(1, M, segnr, Dgi, Dg);

	for (s=2; s<=N; s++) {
		segnr=1; while (Chain->GetSegment(s) != Chain->GetDiffSegment(segnr)) segnr++;
		Forward(s,M,segnr,Dgi,Dg,GetForce(),Info);
	}

	InitializeBackward(M, segnr, Dgi_inv, Dg);

	Composition(N,M,segnr,Dgi_inv,Dgi,Dphi);
	for (s=N-1; s>=1; s--) {
		segnr=1; while (Chain ->GetSegment(s) != Chain -> GetDiffSegment(segnr)) segnr++;
		Backward(M,segnr,Dgi_inv,Dg,Dgx,-GetForce(),Info);
		Composition(s,M,segnr,Dgi_inv,Dgi,Dphi);
	}

        TransferDataToHost(M,Hgi_inv,Dgi_inv);
	lnGN = Lat->ComputeLnGN(Gi_inv);

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
	if (freedom == rangeRestricted) {
		theta = ComputeTheta();
		phiBulk = theta/exp(lnGN);
		SetPhiBulk(phiBulk);
	}

	FreeMemoryOnDevice(Dphi);
	FreeMemoryOnDevice(Dg);
	FreeMemoryOnDevice(Dgi);
	FreeMemoryOnDevice(Dgi_inv);
	FreeMemoryOnDevice(Dgx);
	free(Hphi); free(Hg);
#else
	Message(fatal,"Programming error: entered Cuda enabled routine but the compilation was not done accordingly. Compile SFBox sith Cuda=1");

#endif
}

void
SF_Copol1stO::CudaSymMatrix1(const DensityPart DensPart) {
#ifdef CUDA
	int z,s,s_inv,i,k,Ndiv2;
	double Chalf=0.5;
	int Info[11];
	int segnr=0;
	Vector phi,G;
	Vector Gi_inv(1,M);
	Lat->GetLatticeInfo(Info);
	Ndiv2=N/2;

	double *Hphi, *Hg;
	double *Hgi_inv = &Gi_inv[1];  //pointers to host memory;
	double *Dphi, *Dg, *Dgi, *Dgi_inv, *Dgx; //pointers to device memory;
	Dphi =    (double *)AllocateMemoryOnDevice(M*numDiffSegments);
	Dg =      (double *)AllocateMemoryOnDevice(M*numDiffSegments);
	Dgi =     (double *)AllocateMemoryOnDevice(M*Ndiv2);
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
	segnr=1; while (Chain->GetSegment(1) != Chain->GetDiffSegment(segnr)) segnr++;
	InitializeForward(1, M, segnr, Dgi, Dg);

	for (s=2; s<=N/2; s++) {
		segnr=1; while (Chain->GetSegment(s) != Chain->GetDiffSegment(segnr)) segnr++;
		Forward(s,M,segnr,Dgi,Dg,GetForce(),Info);
	}

	s = N/2;

	InitializeBackward(M,s,Dgi_inv,Dgi); //truck; s en segnr worden in .cu file op dezelfde manier behandeld.

	if (N%2 == 1) {
		s = N/2 + 1;
		segnr=1; while (Chain->GetSegment(s) != Chain->GetDiffSegment(segnr)) segnr++;
		Backward(M,segnr,Dgi_inv,Dg,Dgx,-GetForce(),Info);
		Composition(1,M,segnr,Dgi_inv,Dgi_inv,Dphi); //ook een truck. Hier vul ik voor eerste argument een 1 in om het programma om de tuin te leiden...
		NormPhi(M,segnr,Dphi,Chalf); // correction needed for optimization
	}
	for (s=(N+3)/2; s<=N; s++) {
		s_inv = N-s+1;
		segnr=1; while (Chain->GetSegment(s) != Chain->GetDiffSegment(segnr)) segnr++;
		Backward(M,segnr,Dgi_inv,Dg,Dgx,-GetForce(),Info);
		Composition(s_inv,M,segnr,Dgi_inv,Dgi,Dphi);
	}
	TransferDataToHost(M,Hgi_inv,Dgi_inv);
	lnGN = Lat->ComputeLnGN(Gi_inv);

	CorrectDoubleCounting(M*numDiffSegments,1,Dphi,Dg);
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

	if (freedom == fixedTheta) {
		lnCt = log(theta) - lnGN -log(1.*N);
		phiBulk = theta/exp(lnGN);
		Chain->SetPhiBulk(phiBulk);
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
		lnCt = log(phiRange/(2*phiPreNorm)); // factor 2 is correction needed for optimization
	} else {
		lnCb = log(phiBulk/N);
	}
	for (i=1; i<=numDiffSegments; i++) {
		Seg = Chain->GetDiffSegment(i);
		phi = Seg->GetPhi(DensPart);
		if (freedom == fixedTheta) {
			Lat->NormPhiRestr(phi,Gi_inv,2*theta/N); // factor 2 is correction needed for optimization
		} else if (freedom == rangeRestricted) {
			Lat->NormPhiFree(phi,2*exp(lnCt));// factor 2 is correction needed for optimization
		} else {
			Lat->NormPhiFree(phi,2*phiBulk/N);// factor 2 is correction needed for optimization
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
	FreeMemoryOnDevice(Dphi);
	FreeMemoryOnDevice(Dg);
	FreeMemoryOnDevice(Dgi);
	FreeMemoryOnDevice(Dgi_inv);
	FreeMemoryOnDevice(Dgx);
	free(Hphi); free(Hg);
#else
	Message(fatal,MyInput,"Programming error: entered Cuda enabled routine but the compilation was not done accordingly. Compile SFBox sith Cuda=1");
#endif
}

void
SF_Copol1stO::CudaMatrix1Long(const DensityPart DensPart) {
#ifdef CUDA
	int z,s,s0,t0,v0,t,rs1,i,k,segnr,segnr1;
	int Info[11];
	Vector phi,G;
	//Matrix Gi(1,M,1,n);
	Vector Gi_inv(1,M); //Gs(1,M);

	double *Hphi, *Hg;
	double *Hgi_inv=&Gi_inv[1];  //pointers to host memory;
	double *Dphi, *Dg, *Dgi, *Dgi_inv, *Dgx, *Dgs; //pointers to device memory;

	Lat ->GetLatticeInfo(Info);
	Dphi =    (double *)AllocateMemoryOnDevice(M*numDiffSegments);
	Dg =      (double *)AllocateMemoryOnDevice(M*numDiffSegments);
	Dgi =     (double *)AllocateMemoryOnDevice(M*n);
	Dgi_inv = (double *)AllocateMemoryOnDevice(M);
	Dgx =     (double *)AllocateMemoryOnDevice(M);
	Dgs =     (double *)AllocateMemoryOnDevice(M);
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

	segnr=1; while (Chain->GetSegment(1) != Chain->GetDiffSegment(segnr)) segnr++;
	InitializeForward(1, M, segnr, Dgi, Dg);
	InitializeForward(1, M, segnr, Dgs, Dg);

	//G = Chain->GetSegment(1)->GetSWF();
	//for (z=1; z<=M; z++) {
	//	Gi[z][1] = Gs[z] = G[z];
	//}

	t=1;
	v0=t0=s0 = 0;
	// get the rest of 2..N segments
	for (s=2; s<=N; s++) {
		t++;
		// n is the numbr of matrix elements
		if (t>n) {
			t0++;
			if (t0 == n) t0 = ++v0;
			t = t0 + 1;
			s0 = s - t0 - 1;
		}
		segnr=1; while (Chain->GetSegment(s) != Chain->GetDiffSegment(segnr)) segnr++;
		Backward(M,segnr,Dgs,Dg,Dgx,GetForce(),Info); //this is a forward, but using backword notation.

		//G = Chain->GetSegment(s)->GetSWF();
		//if (force_set) {
		//	Lat->PropagateG(Gs,G,GetForce());
		//} else
		//Lat->PropagateG(Gs,G);

		if ((t == t0+1 && t0 == v0)
		   || (t == t0+1 && ((n-t0)*(n-t0+1) >= N-1-2*(s0+t0)))
		   || (2*(n-t+s) >= N-1))
			//for (z=1; z<=M; z++) Gi[z][t] = Gs[z];
			InitializeForward(t, M, 1, Dgi, Dgs);
	}

	for (s=N; s>=1; s--) {
		segnr=1; while (Chain->GetSegment(s) != Chain->GetDiffSegment(segnr)) segnr++;
		//G = Chain->GetSegment(s)->GetSWF();
		if (s == N) {
			InitializeBackward(M, segnr, Dgi_inv, Dg);
			//for (z=1; z<=M; z++) {Gi_inv[z] = G[z];}

		} else {
			Backward(M,segnr,Dgi_inv,Dg,Dgx,-GetForce(),Info);
			//if (force_set) {
			//	Lat->PropagateG(Gi_inv,G,-GetForce());
			//} else
			//Lat->PropagateG(Gi_inv,G);

		}
		t = s - s0;
		if (t == t0) {
			s0 += - n + t0;
			if (t0 == v0 ) {s0 -= ((n - t0)*(n - t0 + 1))/2;}
			t0 --;
			if (t0 < v0) {v0 = t0;}
			InitializeBackward(M, t, Dgs, Dgi);
			//for (z=1; z<=M; z++) {Gs[z] = Gi[z][t];}
			for (rs1=s0+t0+2; rs1<=s; rs1++) {
				segnr1=1; while (Chain->GetSegment(rs1) != Chain->GetDiffSegment(segnr1)) segnr1++;
				t++;
				Backward(M,segnr1,Dgs,Dg,Dgx,GetForce(),Info);
				//G = Chain->GetSegment(rs1)->GetSWF();
				//if (force_set) {
				//	Lat->PropagateG(Gs,G,GetForce());
				//} else
				//Lat->PropagateG(Gs,G);
				if (t == t0+1 || s0+n == s) {
					InitializeForward(t, M, 1, Dgi, Dgs);
					//for (z=1; z<=M; z++) {Gi[z][t] = Gs[z];}
				}
				if (t == n && s0+n < s) {
					t  = ++t0;
					s0 += n - t0;
				}
			}
			t = n;
		}
		Composition(t,M,segnr,Dgi_inv,Dgi,Dphi);

		//Seg = Chain->GetSegment(s);
		//phi = Seg->GetPhi(DensPart);
		//Lat->ConnectG(Gi_inv,Gi,t,phi);
	}
	TransferDataToHost(M,Hgi_inv,Dgi_inv);
	lnGN = Lat->ComputeLnGN(Gi_inv);
	CorrectDoubleCounting(M*numDiffSegments,1,Dphi,Dg);
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
	if (freedom == rangeRestricted) {
		theta = ComputeTheta();
		phiBulk = theta/exp(lnGN);
		SetPhiBulk(phiBulk);
	}
	FreeMemoryOnDevice(Dphi);
	FreeMemoryOnDevice(Dg);
	FreeMemoryOnDevice(Dgi);
	FreeMemoryOnDevice(Dgi_inv);
	FreeMemoryOnDevice(Dgx);
	FreeMemoryOnDevice(Dgs);
	free(Hphi); free(Hg);
#else
	Message(fatal,MyInput,"Programming error: entered Cuda enabled routine but the compilation was not done accordingly. Compile SFBox sith Cuda=1");
#endif
}


void
SF_Copol1stO::CudaSymMatrix1Long(const DensityPart DensPart) {
#ifdef CUDA
	int z,s,s_inv,s0,t0,v0,t,rs1,i,k,segnr,segnr1;
	double Chalf=0.5;
	int Info[11];
	Vector phi,G;
	//Matrix Gi(1,M,1,n);
	Vector Gi_inv(1,M);
	//Gs(1,M);

	double *Hphi, *Hg;
	double *Hgi_inv=&Gi_inv[1];  //pointers to host memory;
	double *Dphi, *Dg, *Dgi, *Dgi_inv, *Dgx, *Dgs; //pointers to device memory;

	Lat ->GetLatticeInfo(Info);
	Dphi =    (double *)AllocateMemoryOnDevice(M*numDiffSegments);
	Dg =      (double *)AllocateMemoryOnDevice(M*numDiffSegments);
	Dgi =     (double *)AllocateMemoryOnDevice(M*n);
	Dgi_inv = (double *)AllocateMemoryOnDevice(M);
	Dgx =     (double *)AllocateMemoryOnDevice(M);
	Dgs =     (double *)AllocateMemoryOnDevice(M);
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

	segnr=1; while (Chain->GetSegment(1) != Chain->GetDiffSegment(segnr)) segnr++;
	InitializeForward(1, M, segnr, Dgi, Dg);
	InitializeForward(1, M, segnr, Dgs, Dg);

	//G = Chain->GetSegment(1)->GetSWF();
	//for (z=1; z<=M; z++) {
	//	Gi[z][1] = Gs[z] = G[z];
	//}
	t=1;
	v0=t0=s0 = 0;
	for (s=2; s<=N/2; s++) {
		segnr=1; while (Chain->GetSegment(s) != Chain->GetDiffSegment(segnr)) segnr++;
		t++;
		if (t>n) {
			t0++;
			if (t0 == n) t0 = ++v0;
			t = t0 + 1;
			s0 = s - t0 - 1;
		}
		Backward(M,segnr,Dgs,Dg,Dgx,GetForce(),Info); //toch voorwaards, maar met backward procedure
		//G = Chain->GetSegment(s)->GetSWF();
		//if (force_set) {
		//	Lat->PropagateG(Gs,G,GetForce());
		//} else
		//Lat->PropagateG(Gs,G);
		if ((t == t0+1 && t0 == v0)
		   || (t == t0+1 && ((n-t0)*(n-t0+1) >= N-1-2*(s0+t0)))
		   || (2*(n-t+s) >= N-1)) {
			InitializeForward(t, M, 1, Dgi, Dgs); //the way to store....sorry...
			//for (z=1; z<=M; z++) {Gi[z][t] = Gs[z];}
		}
	}

	InitializeBackward(M,1,Dgi_inv,Dgs);
	//for (z=1; z<=M; z++) {Gi_inv[z] = Gs[z];}
	if (N%2 == 1) {
		s = N/2 + 1;
		segnr=1; while (Chain->GetSegment(s) != Chain->GetDiffSegment(segnr)) segnr++;
		Backward(M,segnr,Dgi_inv,Dg,Dgx,-GetForce(),Info);
		Composition(1,M,segnr,Dgi_inv,Dgi_inv,Dphi);
		NormPhi(M,segnr,Dphi,Chalf);
		//Seg = Chain->GetSegment(s);
		//G = Seg->GetSWF();
		//if (force_set) {
		//	Lat->PropagateG(Gi_inv,G,GetForce());
		//} else	Lat->PropagateG(Gi_inv,G);
		//phi = Seg->GetPhi(DensPart);
		//Lat->ConnectG(Gi_inv,Gi_inv,phi);
		//Lat->NormPhiFree(phi,0.5);
	}
	for (s=(N+3)/2; s<=N; s++) {
		segnr=1; while (Chain->GetSegment(s) != Chain->GetDiffSegment(segnr)) segnr++;
		Backward(M,segnr,Dgi_inv,Dg,Dgx,-GetForce(),Info);
		//G = Chain->GetSegment(s)->GetSWF();
		//if (force_set) {
		//	Lat->PropagateG(Gi_inv,G,GetForce());
		//} else	Lat->PropagateG(Gi_inv,G);
		s_inv = N-s+1;
		t = N - s + 1 - s0;
		if (t == t0) {
			s0 += - n + t0;
			if (t0 == v0) s0 -= ((n - t0)*(n - t0 + 1))/2;
			t0 --;
			if (t0 < v0)	v0 = t0;
			InitializeBackward(M,t,Dgs,Dgi);
			//for (z=1; z<=M; z++) {Gs[z] = Gi[z][t];}
			for (rs1=s0+t0+2; rs1<=(N-s+1); rs1++) {
				segnr1=1; while (Chain->GetSegment(rs1) != Chain->GetDiffSegment(segnr1)) segnr1++;
				t++;
				Backward(M,segnr1,Dgs,Dg,Dgx,-GetForce(),Info);
				//G = Chain->GetSegment(rs1)->GetSWF();
				//if (force_set) {
				//	Lat->PropagateG(Gs,G,GetForce());
				//} else	Lat->PropagateG(Gs,G);
				if (t == t0+1 || s0+n == N-s+1) {
					InitializeForward(t, M, 1, Dgi, Dgs);
					//for (z=1; z<=M; z++) {Gi[z][t] = Gs[z];}
				}
				if (t == n && s0+n < N-s+1) {
					t  = ++t0;
					s0 += n - t0;
				}
			}
			t = n;
		}
		Composition(t,M,segnr,Dgi_inv,Dgi,Dphi);
		//Seg = Chain->GetSegment(s);
		//phi = Seg->GetPhi(DensPart);
		//Lat->ConnectG(Gi_inv,Gi,t,phi);

	}
	TransferDataToHost(M,Hgi_inv,Dgi_inv);
	lnGN = Lat->ComputeLnGN(Gi_inv);

	CorrectDoubleCounting(M*numDiffSegments,1,Dphi,Dg);
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
	if (freedom == fixedTheta) {
		lnCt = log(theta) - lnGN -log(1.*N);
		phiBulk = theta/exp(lnGN);
		Chain->SetPhiBulk(phiBulk);
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
		lnCt = log(phiRange/(2*phiPreNorm)); // factor 2 is correction needed for optimization
	} else {
		lnCb = log(phiBulk/N);
	}
	for (i=1; i<=numDiffSegments; i++) {
		Seg = Chain->GetDiffSegment(i);
		phi = Seg->GetPhi(DensPart);
		if (freedom == fixedTheta) {
			Lat->NormPhiRestr(phi,Gi_inv,2*theta/N); // factor 2 is correction needed for optimization
		} else if (freedom == rangeRestricted) {
			Lat->NormPhiFree(phi,2*exp(lnCt));// factor 2 is correction needed for optimization
		} else {
			Lat->NormPhiFree(phi,2*phiBulk/N);// factor 2 is correction needed for optimization
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
	FreeMemoryOnDevice(Dphi);
	FreeMemoryOnDevice(Dg);
	FreeMemoryOnDevice(Dgi);
	FreeMemoryOnDevice(Dgi_inv);
	FreeMemoryOnDevice(Dgx);
	FreeMemoryOnDevice(Dgs);
	free(Hphi); free(Hg);

#else
	Message(fatal,MyInput,"Programming error: entered Cuda enabled routine but the compilation was not done accordingly. Compile SFBox sith Cuda=1");
#endif

}


