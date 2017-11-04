#include "SF_Homopol1stO.h"
//#include <iostream>
#ifdef CUDA
#include "Cuda_tools.h"
#endif

SF_Homopol1stO::SF_Homopol1stO(Text name_,
								   SF_SegmentList* SegQ_,
								   Lat1stO* Lat_,
								   Input* MyInput_)
: SF_Molecule(name_, SegQ_, Lat_, MyInput_)
{
	Lat = Lat_;
	Seg = Chain->GetSegment(1);

	force_set = (!GetForce()==0);

	//if (MayerSaupeSet()) {
		//cout << "punt 1" << endl;
	//	Message(fatal,MyInput,"In SF_homopol1stO: Mayer_Saupe not implemented for first order propagators");
	//}


	if (force_set) {
		Message(fatal,MyInput,
					"Force ensemble not implemented in Homopol1stO; do not know how to deal with sign of force when inversion symmetry is imposed.");
	}

#ifdef CUDA
	if (GetGPU()) {GPUactivated = true;
	}else GPUactivated = false;
#else
	GPUactivated = false;
	if (GetGPU()) {
		Message(literal,"Compile with nvcc (option CUDA=1) to activate GPU... Contact Frans Leermakers. Going classical");
	}
#endif

	if (freedom == thirdGeneration) {
		CreatePhi(constrained);
		CreatePhi(unconstrained);
		CreatePhi(renorm);
	} else if (freedom == secondGeneration) {
		CreatePhi(constrained);
		CreatePhi(unconstrained);
	} else {
		CreatePhi(total);
	}
	n = int(pow(N*3,1.0/3)+0.5);
	if (2*N < 60) n++;
	if (2*N < 30) n++;
	SetPhiBulk(phiBulk);
}

SF_Homopol1stO::~SF_Homopol1stO() {
	if (freedom == thirdGeneration) {
		DeletePhi(constrained);
		DeletePhi(unconstrained);
		DeletePhi(renorm);
	} else if (freedom == secondGeneration) {
		DeletePhi(constrained);
		DeletePhi(unconstrained);
	} else {
		DeletePhi(total);
	}
}

void
SF_Homopol1stO::ReComputePhi() {
	Message(fatal,"ReComputePhi not implemented in SF_Homopol1stO");
}

Vector
SF_Homopol1stO::GetLong(){
	Message(fatal,"GetLong not implemented in SF_Homopol1stO");
	Vector x;
	return x;
}
Vector
SF_Homopol1stO::GetShort(){
	Message(fatal,"GetShort not implemented in SF_Homopol1stO");
	Vector x;
	return x;
}
Vector
SF_Homopol1stO::GetBondOrientation( const DensityPart DensPar) {
	Message(fatal,"GetBondOrientation not implemented in SF_Homopol1stO");
	Vector x;
	return x; 
}; 

MoleculeType
SF_Homopol1stO::GetMoleculeType() const {
	return Chain->GetMoleculeType();
}



void
SF_Homopol1stO::ComputePhi() {
	if (Lat->GetNamesBulkBoundaries().Upperbound() > 0) {
		CreatePhi(bulk);
		MatrixBulk(bulk);
	}
	Vector G = Seg->GetSWF();
	Lat->SetBoundaries(G);
	if (freedom == thirdGeneration) {
		if (GPUactivated) {Message(literal,"GPU declined for third generation; going classical");}
		if (!saveMemory) Matrix3rdGen(LatRange,constrained,unconstrained,renorm,
			LatRangeRenorm);
		else  Matrix3rdGenLong(LatRange,constrained,unconstrained,renorm,
			LatRangeRenorm);
	} else if (freedom == secondGeneration) {
		if (GPUactivated) {Message(literal,"GPU declined for second generation; going classical");}
		if (!saveMemory) Matrix2ndGen(LatRange,constrained,unconstrained);
		else Matrix2ndGenLong(LatRange,constrained,unconstrained);
	} else {
		if (!saveMemory) {
			if (GPUactivated) CudaMatrix1(total); else Matrix1(total);
		}
		else { if (GPUactivated) CudaMatrix1Long(total); else Matrix1Long(total);
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
}
void
SF_Homopol1stO::GetLayerAnalysis(Output* Out) {
	if (LatRangeTrains == NULL || LatRangeStartLoops == NULL) {
		return;
	}
	int z,s;

	Vector G = Seg->GetSWF();
	Matrix Gads(1,M,1,N);
	Vector Gads_inv(1,M);
	Matrix Gfree(1,M,1,N);
	Vector Gfree_inv(1,M);
	Vector phiTrainLoops(1,M);
	Vector phiTails(1,M);
	Vector phiFree(1,M);
	Lat->MakeSafe(G);
	Lat->Init2G(Gads_inv,Gfree_inv,G,LatRangeTrains);
	for (z=1; z<=M; z++) {
		Gads[z][1] = Gads_inv[z];
		Gfree[z][1] = Gfree_inv[z];
	}
	Lat->MakeSafe(phiTrainLoops);
	Lat->MakeSafe(phiTails);
	Lat->MakeSafe(phiFree);
	for (s=2; s<=N; s++) {
		if (force_set) {
			Lat->Propagate2G(Gads,Gfree,G,s,LatRangeTrains,GetForce());
		} else
		Lat->Propagate2G(Gads,Gfree,G,s,LatRangeTrains);
	}
	Lat->ConnectG(Gads_inv,Gads,N,phiTrainLoops);
	Lat->Connect2G(Gads_inv,Gads,N,Gfree_inv,Gfree,N,phiTails);
	Lat->ConnectG(Gfree_inv,Gfree,N,phiFree);

	for (s=N-1; s>=1; s--) {
		if (force_set) {
			Lat->Propagate2G(Gads_inv,Gfree_inv,G,LatRangeTrains,GetForce());
		} else
		Lat->Propagate2G(Gads_inv,Gfree_inv,G,LatRangeTrains);
		Lat->ConnectG(Gads_inv,Gads,s,phiTrainLoops);
		Lat->Connect2G(Gads_inv,Gads,s,Gfree_inv,Gfree,s,phiTails);
		Lat->ConnectG(Gfree_inv,Gfree,s,phiFree);
	}
	//return;
	Lat->CorrectDoubleCountG(phiTrainLoops,G);
	Lat->CorrectDoubleCountG(phiTails,G);
	Lat->CorrectDoubleCountG(phiFree,G);
	if (freedom == fixedTheta) {
		// special hack needed here to do it right for overflow protection
		Vector GiNorm(1,M);
		Lat->MakeSafe(GiNorm);
		Lat->ConnectG(G,Gfree,N,GiNorm);
		Lat->ConnectG(G,Gads,N,GiNorm);
		Lat->CorrectDoubleCountG(GiNorm,G);
		Lat->NormPhiRestr(phiTrainLoops,GiNorm,theta/N);
		Lat->NormPhiRestr(phiTails,GiNorm,theta/N);
		Lat->NormPhiRestr(phiFree,GiNorm,theta/N);
	} else {
		Lat->NormPhiFree(phiTrainLoops,phiBulk/N);
		Lat->NormPhiFree(phiTails,phiBulk/N);
		Lat->NormPhiFree(phiFree,phiBulk/N);
	}
	Lat->RestoreFromSafe(phiTrainLoops);
	Lat->RestoreFromSafe(phiTails);
	Lat->RestoreFromSafe(phiFree);
	Lat->RestoreFromSafe(G);
	double GNads = exp(Lat->ComputeLnGN(Gads_inv));
	Out->PutReal("mol",name,"GNads",GNads);
	double data = 0;
	Lat->SubtractBoundaries(phiTrainLoops);
	Lat->SubtractBoundaries(phiTails);

	for (z=1; z<=M; z++) {
		data += (phiTrainLoops[z] + phiTails[z])*Lat->GetNumLatticeSites(z);
	}
	Out->PutReal("mol",name,"theta ads",data);
	Lat->RestoreBoundaries(phiTrainLoops);
	Lat->RestoreBoundaries(phiTails);
	z = 2;
	while (phiTrainLoops[z] > phiTails[z] && z != M) {
		z++;
	}
	double RC1 = phiTrainLoops[z] - phiTrainLoops[z-1];
	double RC2 = phiTails[z] - phiTails[z-1];
	double a1 = phiTrainLoops[z-1] - RC1*(z-1);
	double a2 = phiTails[z-1] - RC2*(z-1);
	data = (a1 - a2)/(RC2 - RC1) - 1.5;
	Out->PutReal("mol",name,"z star",data);
	Out->PutProfile("mol",name,"phi trains+loops",phiTrainLoops);
	Out->PutProfile("mol",name,"phi tails",phiTails);
	Out->PutProfile("mol",name,"phi free",phiFree);
	// size distribution functions, not implemented for a overflow protection
	// two variants here, first one general, second one for 1D adsorption
	Vector numTrains(1,N);
	Vector numLoops(1,N);
	Vector numTails(1,N);
	Boolean general = false;
	Boolean found = false;
	for (z=1; z<=M; z++) {
		if (LatRangeTrains->InRange(z)) {
			if (found) {
				general = true;
			}
			found = true;
		}
	}

	found = false;
	for (z=1; z<=M; z++) {
		if (LatRangeStartLoops->InRange(z)) {
			if (found) {
				general = true;
			}
			found = true;
		}
	}
	if ((Lat->GetNumGradients() > 1 && !Lat->OverflowProtection()) || general) {
	//general version
	// trains
	Vector G2(1,M);
	Vector GBackup(1,M);
	for (z=1; z<=M; z++) {
		GBackup[z] = G[z]; // copy G to GBackup
	}
	Lat->Init2G(G,G2,GBackup,LatRangeTrains);
	for (z=1; z<=M; z++) {
		if (LatRangeTrains->InRange(z)) {
			G2.Dim(1,M);
			G2[z] = GBackup[z];
			Lat->SetBoundaries(G2);
			Matrix Gt(1,M,1,N);
			for (int zq=1; zq<=M; zq++) {
				Gt[zq][1] = G2[zq];
			}
			for (s=2; s<=N; s++) {
				if (force_set) {
					Lat->PropagateG(Gt,G,s,GetForce());
				} else
				Lat->PropagateG(Gt,G,s);
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
        									value += (Gads[z4][N-s]
												+ Gfree[z4][N-s])/ltrlp1;
        								} else {
        									value += 1/(ltrlp1*ltrlp2);
        								}
										for (int t=1; t<=N-s-1; t++) {
											value += (Gads[z2][t]
												+ Gfree[z2][t])*(Gads[z4][N-s-t]
												+ Gfree[z4][N-s-t]);
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
	Lat->Init2G(G2,G,GBackup,LatRangeTrains);
	for (z=1; z<=M; z++) {
		if (LatRangeStartLoops->InRange(z)) {
			G2.Dim(1,M);
			G2[z] = GBackup[z];
			Lat->SetBoundaries(G2);
			Matrix Gt(1,M,1,N);
			for (int zq=1; zq<=M; zq++) {
				Gt[zq][1] = G2[zq];
			}
			for (s=2; s<=N-2; s++) {
				if (force_set) {
					Lat->PropagateG(Gt,G,s,GetForce());
				} else
				Lat->PropagateG(Gt,G,s);
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
												*Gads[z4][N-s-t];
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
	// give G its original value back
	for (z=1; z<=M; z++) {
		G[z] = GBackup[z];
	}
	// tails
	for (s=1; s<=N-1; s++) {
		for (int ztr=1; ztr<=M; ztr++) {
			if (LatRangeTrains->InRange(ztr)) {
				for (int ztl=1; ztl<=M; ztl++) {
					if (LatRangeStartLoops->InRange(ztl)) {
						numTails[s] += 2*Lat->GetNumLatticeSites(ztl)
							*Lat->GetLambda(ztl,ztr)*Gfree[ztl][s]
							*Gads[ztr][N-s]/GNads;
					}
				}
			}
		}
	}

	} else if (!Lat->OverflowProtection()) { // routines optimized for 1D
		int ztr=0; // only one layer with trains
		int zlp=0; // only one layer with loops adjoining trains
		Boolean found = false;
		for (z=1; z<=M; z++) {
			if (LatRangeTrains->InRange(z)) {
				ztr = z;
				if (found) Message(implementation,"your trains_range contains more than "
					"1 lattice layer");
				found = true;
			}
		}
		found = false;
		for (z=1; z<=M; z++) {
			if (LatRangeStartLoops->InRange(z)) {
				zlp = z;
				if (found) Message(implementation,"your start_loops_range contains more "
					"than 1 lattice layer");
				found = true;
			}
		}
		double ltrlp = Lat->GetLambda(ztr,zlp);
		double Ltr = Lat->GetNumLatticeSites(ztr);
		double llptr = Lat->GetLambda(zlp,ztr);
		double Llp = Lat->GetNumLatticeSites(zlp);
		if (ltrlp <=0) {
			Message(fatal,"your start_loops_range and trains_range are not "
				"adjoining");
		}
		if (ztr == zlp) {
			Message(fatal,"your start_loops_range and trains_range are "
				"overlapping");
		}
		// trains
        for (s=1; s<=N; s++) {
        	if (s<N) {
        		numTrains[s] += 2*(Gads[zlp][N-s] + Gfree[zlp][N-s])/ltrlp;
        	} else {
        		numTrains[s] += 1/(ltrlp*ltrlp);
        	}
        	for (int t=1; t<=N-s-1; t++) {
        		numTrains[s] += ((Gads[zlp][t] + Gfree[zlp][t])
					*(Gads[zlp][N-s-t] + Gfree[zlp][N-s-t]));
        	}
        	numTrains[s] *= Ltr*(ltrlp*ltrlp
				*pow(Lat->GetLambda(ztr,ztr)*G[ztr],s-1)*G[ztr]/GNads);
        }
        // loops
		Vector Gt(1,M);

		Lat->Init2G(Gads_inv,Gfree_inv,G,LatRangeStartLoops); // Gads_inv and Gads_free used here to safe memory
		for (z=1; z<=M; z++) {
		 	Gt[z] = Gads_inv[z];
		  	Gfree_inv[z] = G[z]; // Gads_free used here to temporary safe info from G
		}
		Lat->Init2G(Gads_inv,G,Gfree_inv,LatRangeTrains); // G temporary overwritten
		for (s=1; s<=N-2; s++) {
			if (s > 1) {
				if (force_set) {
					Lat->PropagateG(Gt,G,GetForce());
				} else
				Lat->PropagateG(Gt,G);
			}
			for (int t=1; t<=N-s-1; t++) {
				numLoops[s] += (Gads[ztr][t]*Gads[ztr][N-s-t]);
			}
			numLoops[s] *= Llp*(llptr*llptr*Gt[zlp]/GNads);
        }
        for (z=1; z<=M; z++) {
        	G[z] = Gfree_inv[z]; // give G its original value again
        }
		// tails
        for (s=1; s<=N-1; s++) {
        	numTails[s] = 2*Llp*llptr*Gfree[zlp][s]*Gads[ztr][N-s]/GNads;
        }
 	}
 	if (Lat->OverflowProtection()) {
		Message(implementation, "train,loop and tail size distribution not"
			" implemented for overflow protection");
	}
	if (!Lat->OverflowProtection()) {
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
			fluctFrTrains += s*s*numTrains[s];
			fluctFrLoops += s*s*numLoops[s];
			fluctFrTails += s*s*numTails[s];
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

		double phiads = frTrains + frLoops + frTails;
		Lat->RestoreBoundaries(phiTrainLoops);
		Lat->RestoreBoundaries(phiTails);
		frTrains /= phiads;
		frLoops /= phiads;
		frTails /= phiads;
		Out->PutReal("mol", name, "trains fract", frTrains);
		Out->PutReal("mol", name, "loops fract", frLoops);
		Out->PutReal("mol", name, "tails fract", frTails);
	}
}/*
void
SF_Homopol1stO::GetSegmentAnalysis(Output* Out) {
	int z,s;
	Vector G;	// the segment weighting factors for the one segment type
	Matrix Gs(1,M,1,N);	// the propagator. Only one is needed, since
			//			Gs[z,s] = Gs_inv[z,N-s+1]
	Matrix phi(1,M,1,N);	// the densities for all segments: phi[z,s]

	G = Seg->GetSWF();
	for (z=1; z<=M; z++) {
		Gs(z,1)=G(z);
	}

}*/

// This is the classical (simple) version of Matrix1
void
SF_Homopol1stO::Matrix1(const DensityPart DensPart) {
	int z,s,s_inv;
	Matrix Gi(1,M,1,N/2); // the homopolymer is always symmetric, Mx(N/2) matrix is sufficient
	// from Jan Van Male thesis: Gi[z][s] is the SWF of segment s in molecule i in layer z
	Vector Gi_inv(1,M);
	// in constructor: Seg=Chain->GetSegment(1)
	// i.e. Seg is the 1st segment of the SF_SegmentBlock *Chain
	Vector phi = Seg->GetPhi(DensPart);
	Vector G = Seg->GetSWF();
	Lat->MakeSafe(G); // does nothing! if no safe mode required
	for (z=1; z<=M; z++) {
		phi[z] = 0;
		Gi[z][1] = G[z];
	}
	Lat->MakeSafe(phi);
	// PropragateG until N/2
	for (s=2; s<=N/2; s++) {
		if (force_set) {
			Lat->PropagateG(Gi,G,s,GetForce());
		} else
		Lat->PropagateG(Gi,G,s);
	}
	s = N/2; // correct fot the fact that in the end, s=N/2+1
	for (z=1; z<=M; z++) {
		Gi_inv[z] = Gi[z][s];
	}
	// if N is odd, do it defferently for the central segment
	if (N%2 == 1) {
		if (force_set) {
			Lat->PropagateG(Gi_inv,G,GetForce());
		} else
		Lat->PropagateG(Gi_inv,G);
		Lat->ConnectG(Gi_inv,Gi_inv,phi);
		Lat->NormPhiFree(phi,0.5); // correction needed for optimization
	}
	for (s=(N+3)/2; s<=N; s++) {
		s_inv = N-s+1;
		if (force_set) {
			Lat->PropagateG(Gi_inv,G,-GetForce());
		} else
		Lat->PropagateG(Gi_inv,G);
		Lat->ConnectG(Gi_inv,Gi,s_inv,phi); // optimization, chain is symmetric
			// Connect only half the chain into phi
	}
	Lat->CorrectDoubleCountG(phi,G);
	lnGN = Lat->ComputeLnGN(Gi_inv);
	if (freedom == fixedTheta) {
		lnCt = log(theta)-lnGN-log(1.0*N);
		Lat->NormPhiRestr(phi,Gi_inv,2*theta/N); // factor 2 optimization correction
		phiBulk = theta/exp(lnGN);
		SetPhiBulk(phiBulk);
	} else if (freedom == rangeRestricted) {
		double phiPreNorm = 0;
		for (z=1; z<=M; z++) {
			if (restrictedRange->InRange(z)) {
				phiPreNorm += phi[z];
			}
		}
		lnCt = log(phiRange/(2*phiPreNorm));// factor 2 optimization correction
		Lat->NormPhiFree(phi,phiRange/phiPreNorm);
		theta = ComputeTheta();
		phiBulk = exp(lnCt)*N;
		SetPhiBulk(phiBulk);
	} else {
		Lat->NormPhiFree(phi,2*phiBulk/N);// factor 2 optimization correction
		lnCb = log(phiBulk/N);
	}
	Lat->RestoreFromSafe(G);
	Lat->RestoreFromSafe(phi);
}

void
SF_Homopol1stO::CudaMatrix1(const DensityPart DensPart) {
#ifdef CUDA
	int z,s,s_inv,Ndiv2;
	double Chalf=0.5;
	int Info[11];
	int segnr=0;
	Vector phi,G;
	Vector Gi_inv(1,M);
	Lat->GetLatticeInfo(Info);

	Ndiv2=N/2;

	double *Hphi, *Hg, *Hgi_inv;  //pointers to host memory;
	double *Dphi, *Dg, *Dgi, *Dgi_inv, *Dgx; //pointers to device memory;
	Dphi =    (double *)AllocateMemoryOnDevice(M);
	Dg =      (double *)AllocateMemoryOnDevice(M);
	Dgi =     (double *)AllocateMemoryOnDevice(M*Ndiv2);
	Dgi_inv = (double *)AllocateMemoryOnDevice(M);
	Dgx =     (double *)AllocateMemoryOnDevice(M);
	Hphi =    (double *)malloc(sizeof(double)*M);
	Hg =      (double *)malloc(sizeof(double)*M);
	Hgi_inv = (double *)malloc(sizeof(double)*M);


	Seg = Chain->GetSegment(1);
	G = Seg->GetSWF();
	phi = Seg->GetPhi(DensPart);
	for (z=1; z<=M; z++) {
		phi[z] = 0;
		Hphi[z-1] = 0;
		Hg[z-1] =G[z];
	}

	TransferDataToDevice(M,Hphi, Dphi);
	TransferDataToDevice(M,Hg, Dg);
	segnr=1;
	InitializeForward(1, M, segnr, Dgi, Dg);
cout.precision(15);
cout << "norm1=" << Norm2(M,1,Dgi) << endl;
	for (s=2; s<=N/2; s++) {
		Forward(s,M,segnr,Dgi,Dg,GetForce(),Info);


	}
	s = N/2;

	InitializeBackward(M,s,Dgi_inv,Dgi); //truck; s en segnr worden in .cu file op dezelfde manier behandeld.
	if (N%2 == 1) {
		s = N/2 + 1;
		Backward(M,segnr,Dgi_inv,Dg,Dgx,-GetForce(),Info);
		Composition(1,M,segnr,Dgi_inv,Dgi_inv,Dphi); //ook een truck. Hier vul ik voor eerste argument een 1 in om het programma om de tuin te leiden...
		NormPhi(M,segnr,Dphi,Chalf); // correction needed for optimization
	}
	for (s=(N+3)/2; s<=N; s++) {
		s_inv = N-s+1;
		Backward(M,segnr,Dgi_inv,Dg,Dgx,-GetForce(),Info);
		Composition(s_inv,M,segnr,Dgi_inv,Dgi,Dphi);
	}
cout << "normN=" << Norm2(M,1,Dgi_inv) << endl;
	TransferDataToHost(M,Hgi_inv,Dgi_inv);
	for (int z=1; z<=M; z++) Gi_inv[z]=Hgi_inv[z-1];
	lnGN = Lat->ComputeLnGN(Gi_inv);
	CorrectDoubleCounting(M,1,Dphi,Dg);
	TransferDataToHost(M,Hphi,Dphi);
	phi = Seg->GetPhi(DensPart);
	for (int z=1; z<=M; z++) {
		phi[z]=Hphi[z-1];
	}

	if (freedom == fixedTheta) {
		lnCt = log(theta)-lnGN-log(1.0*N);
		Lat->NormPhiRestr(phi,Gi_inv,2*theta/N); // factor 2 optimization correction
		phiBulk = theta/exp(lnGN);
		SetPhiBulk(phiBulk);
	} else if (freedom == rangeRestricted) {
		double phiPreNorm = 0;
		for (z=1; z<=M; z++) {
			if (restrictedRange->InRange(z)) {
				phiPreNorm += phi[z];
			}
		}
		lnCt = log(phiRange/(2*phiPreNorm));// factor 2 optimization correction
		Lat->NormPhiFree(phi,phiRange/phiPreNorm);
		theta = ComputeTheta();
		phiBulk = exp(lnCt)*N;
		SetPhiBulk(phiBulk);
	} else {
		Lat->NormPhiFree(phi,2*phiBulk/N);// factor 2 optimization correction
		lnCb = log(phiBulk/N);
	}
	Lat->RestoreFromSafe(G);
	Lat->RestoreFromSafe(phi);
	FreeMemoryOnDevice(Dphi);
	FreeMemoryOnDevice(Dg);
	FreeMemoryOnDevice(Dgi);
	FreeMemoryOnDevice(Dgi_inv);
	FreeMemoryOnDevice(Dgx);
	free(Hphi); free(Hg); free(Hgi_inv);
#else
	Message(fatal,MyInput,"Programming error: entered Cuda enabled routine but the compilation was not done accordingly. Compile SFBox sith Cuda=1");
#endif
}


// This is the memory-saving long version of Matrix1
// Forget about it for now and study Matrix1() first
void
SF_Homopol1stO::Matrix1Long(const DensityPart DensPart) {
	int z,s,s0,t0,v0,t,rs1;
	Matrix Gi(1,M,1,n);
	Vector Gi_inv(1,M), Gs(1,M);
	// in constructor: Seg=Chain->GetSegment(1)
	Vector phi = Seg->GetPhi(DensPart);
	Vector G = Seg->GetSWF();
	Lat->MakeSafe(G);
	// G and phi for the first segment
	for (z=1; z<=M; z++) {
		Gi[z][1] = Gs[z] = G[z];
		phi[z] = 0;
	}
    Lat->MakeSafe(phi);
	t=1;
	v0=t0=s0 = 0;
	// and for the other 2..N
	for (s=2; s<=N/2; s++) {
		t++;
		// n=(3*N)^{1/3}+0.5 is the matrix size
		if (t>n) {    //(++t>n)
			t0++;
			if (t0 == n) t0 = ++v0;
			t = t0 + 1;
			s0 = s - t0 - 1;
		}
		if (force_set) {
			Lat->PropagateG(Gs,G,GetForce());
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
		if (force_set) {
			Lat->PropagateG(Gi_inv,G,GetForce());
		} else
		Lat->PropagateG(Gi_inv,G);
		Lat->ConnectG(Gi_inv,Gi_inv,phi);
		Lat->NormPhiFree(phi,0.5); // correction needed for optimization
	}
	for (s=(N+3)/2; s<=N; s++) {
		if (force_set) {
			Lat->PropagateG(Gi_inv,G,GetForce());
		} else
		Lat->PropagateG(Gi_inv,G);
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
				if (force_set) {
					Lat->PropagateG(Gs,G,GetForce());
				} else
				Lat->PropagateG(Gs,G);
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
		Lat->ConnectG(Gi_inv,Gi,t,phi); // optimization, chain is symmetric
	}
	Lat->CorrectDoubleCountG(phi,G);
	lnGN = Lat->ComputeLnGN(Gi_inv);
	if (freedom == fixedTheta) {
		lnCt = log(theta)-lnGN-log(1.0*N);
		Lat->NormPhiRestr(phi,Gi_inv,2*theta/N); // factor 2 optimization correction
		phiBulk = theta/exp(lnGN);
		SetPhiBulk(phiBulk);
	} else if (freedom == rangeRestricted) {
		double phiPreNorm = 0;
		for (z=1; z<=M; z++) {
			if (restrictedRange->InRange(z)) {
				phiPreNorm += phi[z];
			}
		}
		lnCt = log(phiRange/(2*phiPreNorm));// factor 2 optimization correction
		Lat->NormPhiFree(phi,2*exp(lnCt));// factor 2 optimization correction
		theta = ComputeTheta();
		phiBulk = theta/exp(lnGN);
		SetPhiBulk(phiBulk);
	} else {
		Lat->NormPhiFree(phi,2*phiBulk/N);// factor 2 optimization correction
		lnCb = log(phiBulk/N);
	}
	Lat->RestoreFromSafe(G);
	Lat->RestoreFromSafe(phi);
}

void
SF_Homopol1stO::CudaMatrix1Long(const DensityPart DensPart) {
//cout << "CudaMatrix1Long "  << endl;
#ifdef CUDA
	int z,s,s_inv,s0,t0,v0,t,rs1,k,segnr;
	double Chalf=0.5;
	int Info[11];
	//int Div2[N+1]; for (s=0; s<=N; s++) Div2[s]=0;
	Vector phi,G;
	//Matrix Gi(1,M,1,n);
	Vector Gi_inv(1,M);
	//Gs(1,M);

	double *Hphi, *Hg;
	double *Hgi_inv=&Gi_inv[1];  //pointers to host memory;
	double *Dphi, *Dg, *Dgi, *Dgi_inv, *Dgx, *Dgs; //pointers to device memory;

	Lat ->GetLatticeInfo(Info);
	Dphi =    (double *)AllocateMemoryOnDevice(M);
	Dg =      (double *)AllocateMemoryOnDevice(M);
	Dgi =     (double *)AllocateMemoryOnDevice(M*n);
	Dgi_inv = (double *)AllocateMemoryOnDevice(M);
	Dgx =     (double *)AllocateMemoryOnDevice(M);
	Dgs =     (double *)AllocateMemoryOnDevice(M);
	Hphi =    (double *)malloc(sizeof(double)*M);
	Hg =      (double *)malloc(sizeof(double)*M);

	k=0;
	Seg->GetPhi(DensPart);
	G = Seg->GetSWF();
	phi = Seg->GetPhi(DensPart);
	for (z=1; z<=M; z++) {
		phi[z] = 0;
		Hphi[k] = 0;
		Hg[k] =G[z];
		k++;
	}

	TransferDataToDevice(M,Hphi, Dphi);
	TransferDataToDevice(M,Hg, Dg);

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
				t++;
				Backward(M,segnr,Dgs,Dg,Dgx,-GetForce(),Info);
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

	CorrectDoubleCounting(M,1,Dphi,Dg);
	TransferDataToHost(M,Hphi,Dphi);


	k=0;
	phi = Seg->GetPhi(DensPart);
	for (int z=1; z<=M; z++) {
		phi[z]=Hphi[k];
		k++;
	}

	if (freedom == fixedTheta) {
		lnCt = log(theta)-lnGN-log(1.0*N);
		Lat->NormPhiRestr(phi,Gi_inv,2*theta/N); // factor 2 optimization correction
		phiBulk = theta/exp(lnGN);
		SetPhiBulk(phiBulk);
	} else if (freedom == rangeRestricted) {
		double phiPreNorm = 0;
		for (z=1; z<=M; z++) {
			if (restrictedRange->InRange(z)) {
				phiPreNorm += phi[z];
			}
		}
		lnCt = log(phiRange/(2*phiPreNorm));// factor 2 optimization correction
		Lat->NormPhiFree(phi,2*exp(lnCt));// factor 2 optimization correction
		theta = ComputeTheta();
		phiBulk = theta/exp(lnGN);
		SetPhiBulk(phiBulk);
	} else {
		Lat->NormPhiFree(phi,2*phiBulk/N);// factor 2 optimization correction
		lnCb = log(phiBulk/N);
	}
	Lat->RestoreFromSafe(G);
	Lat->RestoreFromSafe(phi);


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
SF_Homopol1stO::Matrix2ndGen(const LatticeRange* Range,
							 const DensityPart DensPart1,
							 const DensityPart DensPart2) {
	int z,s,s_inv;
	Matrix Gi1(1,M,1,N/2);
	Vector Gi_inv1(1,M);
	Matrix Gi2(1,M,1,N/2);
	Vector Gi_inv2(1,M);
	Vector G = Seg->GetSWF();
	Lat->MakeSafe(G);
	Lat->Init2G(Gi_inv1,Gi_inv2,G,Range); // Gi_inv used here to get values
	Vector phi1 = GetPhi(DensPart1);
	Vector phi2 = GetPhi(DensPart2);
	for (z=1; z<=M; z++) {
		Gi1[z][1] = Gi_inv1[z];
		Gi2[z][1] = Gi_inv2[z];
		phi1[z] = 0;
		phi2[z] = 0;
	}
	Lat->MakeSafe(phi1);
	Lat->MakeSafe(phi2);
	for (s=2; s<=N/2; s++)
		if (force_set) {
			Lat->Propagate2G(Gi1,Gi2,G,s,Range,GetForce());
		} else
			Lat->Propagate2G(Gi1,Gi2,G,s,Range);
	s = N/2;
	for (z=1; z<=M; z++) {
		Gi_inv1[z] = Gi1[z][s];
		Gi_inv2[z] = Gi2[z][s];
	}
	if (N%2 == 1) {
		s = N/2 + 1;
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
		if (force_set) {
			Lat->Propagate2G(Gi_inv1,Gi_inv2,G,Range,GetForce());
		} else
		Lat->Propagate2G(Gi_inv1,Gi_inv2,G,Range);
		Lat->ConnectG(Gi_inv1,Gi1,s_inv,phi1); //'loops'
		Lat->Connect2G(Gi_inv1,Gi1,s_inv,Gi_inv2,Gi2,s_inv,phi1); //'tails'
		Lat->ConnectG(Gi_inv2,Gi2,s_inv,phi2); //'free'
	}
	Lat->CorrectDoubleCountG(phi1,G);
	Lat->CorrectDoubleCountG(phi2,G);
	lnGN = Lat->ComputeLnGN(Gi_inv1);
	if (muFixed) {
		Lat->NormPhiFree(phi1,2*exp(lnCt));
		Lat->SubtractBoundaries(phi1);
		Lat->MultiplyWithLatticeSites(phi1);
		theta = 0;
		for (int z=1; z<=M; z++) {
			theta += phi1[z];
		}
		Lat->DivideByLatticeSites(phi1);
		Lat->RestoreBoundaries(phi1);
	} else {
		lnCt = log(theta)-lnGN-log(1.0*N);
		Lat->NormPhiRestr(phi1,Gi_inv1,2*theta/N);
	}
	Lat->NormPhiFree(phi2,2*phiBulk/N);
	lnCb = log(phiBulk/N);
	Lat->RestoreFromSafe(G);
	Lat->RestoreFromSafe(phi1);
	Lat->RestoreFromSafe(phi2);
}
void
SF_Homopol1stO::Matrix2ndGenLong(const LatticeRange* Range,
								 const DensityPart DensPart1,
								 const DensityPart DensPart2) {
	int z,s,s0,t0,v0,t,rs1;
	Matrix Gi1(1,M,1,n);
	Vector Gi_inv1(1,M), Gs1(1,M);
	Matrix Gi2(1,M,1,n);
	Vector Gi_inv2(1,M), Gs2(1,M);
	Vector G = Seg->GetSWF();
	Lat->MakeSafe(G);
	Lat->Init2G(Gs1,Gs2,G,Range);
	Vector phi1 = GetPhi(DensPart1);
	Vector phi2 = GetPhi(DensPart2);
	for (z=1; z<=M; z++) {
		Gi1[z][1] = Gs1[z];
		Gi2[z][1] = Gs2[z];
		phi1[z] = 0;
		phi2[z] = 0;
	}
	Lat->MakeSafe(phi1);
	Lat->MakeSafe(phi2);
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
		if (force_set) {
			Lat->Propagate2G(Gs1,Gs2,G,Range,GetForce());
		} else
		Lat->Propagate2G(Gs1,Gs2,G,Range);
		if ((t == t0+1 && t0 == v0)
		   || (t == t0+1 && ((n-t0)*(n-t0+1) >= N-1-2*(s0+t0)))
		   || (2*(n-t+s) >= N-1))
			for (z=1; z<=M; z++) {
				Gi1[z][t] = Gs1[z];
				Gi2[z][t] = Gs2[z];
			}
	}
	s = N/2;
	for (z=1; z<=M; z++) {
		Gi_inv1[z] = Gs1[z];
		Gi_inv2[z] = Gs2[z];
	}
	if (N%2 == 1) {
		s = N/2 + 1;
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
		if (force_set) {
			Lat->Propagate2G(Gi_inv1,Gi_inv2,G,Range,GetForce());
		} else
		Lat->Propagate2G(Gi_inv1,Gi_inv2,G,Range);
		t = N - s + 1 - s0;
		if (t == t0) {
			s0 += - n + t0;
			if (t0 == v0 )
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
				if (force_set) {
					Lat->Propagate2G(Gs1,Gs2,G,Range,GetForce());
				} else
				Lat->Propagate2G(Gs1,Gs2,G,Range);
				if (t == t0+1 || s0+n == N-s+1)
					for (z=1; z<=M; z++) {
						Gi1[z][t] = Gs1[z];
						Gi2[z][t] = Gs2[z];
					}

				if (t == n && s0+n < N-s+1) {
					t  = ++t0;
					s0 += n - t0;
				}
			}
			t = n;
		}
		Lat->ConnectG(Gi_inv1,Gi1,t,phi1); //'loops'
		Lat->Connect2G(Gi_inv1,Gi1,t,Gi_inv2,Gi2,t,phi1); //'tails'
		Lat->ConnectG(Gi_inv2,Gi2,t,phi2); //'free'
	}
	Lat->CorrectDoubleCountG(phi1,G);
	Lat->CorrectDoubleCountG(phi2,G);
	lnGN = Lat->ComputeLnGN(Gi_inv1);
	if (muFixed) {
		Lat->NormPhiFree(phi1,2*exp(lnCt));
		Lat->SubtractBoundaries(phi1);
		Lat->MultiplyWithLatticeSites(phi1);
		theta = 0;
		for (int z=1; z<=M; z++) {
			theta += phi1[z];
		}
		Lat->DivideByLatticeSites(phi1);
		Lat->RestoreBoundaries(phi1);
	} else {
		lnCt = log(theta)-lnGN-log(1.0*N);
		Lat->NormPhiRestr(phi1,Gi_inv1,2*theta/N);
	}
	Lat->NormPhiFree(phi2,2*phiBulk/N);
	lnCb = log(phiBulk/N);
	Lat->RestoreFromSafe(G);
	Lat->RestoreFromSafe(phi1);
	Lat->RestoreFromSafe(phi2);
}
void
SF_Homopol1stO::Matrix3rdGen(const LatticeRange* LatRangeConstr,
							 const DensityPart DensPart1,
							 const DensityPart DensPart2,
							 const DensityPart DensPart3,
							 const LatticeRange* LatRangeRenormed) {
	int z;
	Matrix2ndGen(LatRangeConstr,DensPart1,DensPart2);
	if (LatBase->GetNumLayers(1) > 1) {
		Vector G = Seg->GetSWF();
		Vector GBackup(1,M);
		for (z=1; z<=M; z++) {
			GBackup[z] = G[z];
		}
		Vector Dummy(1,M);
		Lat->Init2G(Dummy,G,GBackup,LatRange);
		Matrix2ndGen(LatRangeRenormed,DensPart3,DensPart2);
		for (z=1; z<=M; z++) {
			G[z] = GBackup[z];
		}
	}
}
void
SF_Homopol1stO::Matrix3rdGenLong(const LatticeRange* LatRangeConstr,
								 const DensityPart constr,
								 const DensityPart unconstr,
								 const DensityPart renormed,
								 const LatticeRange* LatRangeRenormed) {
	int z;
	Matrix2ndGenLong(LatRangeConstr,constr,unconstr);
	if (LatBase->GetNumLayers(1) > 1) {
		Vector G = Seg->GetSWF();
		Vector GBackup(1,M);
		for (z=1; z<=M; z++) {
			GBackup[z] = G[z];
		}
		Vector Dummy(1,M);
		Lat->Init2G(Dummy,G,GBackup,LatRange);
		Matrix2ndGenLong(LatRangeRenormed,renormed,unconstr);
		for (z=1; z<=M; z++) {
			G[z] = GBackup[z];
		}
	}
}
void
SF_Homopol1stO::MatrixBulk(const DensityPart DensPart) {
	Vector G = Seg->GetSWF();
	Vector phi = GetPhi(DensPart);
	for (int z=1; z<=M; z++) {
		if (Lat->BulkBoundary(z)) {
			phi[z] = phiBulk*pow(G[z],N);
		}
	}
}
void
SF_Homopol1stO::CopyBulkBoundaries(const DensityPart DensPart) {
	Vector phi = GetPhi(DensPart);
	Vector PhiBulk = GetPhi(bulk);
	for (int z=1; z<=M; z++) {
		if (Lat->BulkBoundary(z)) {
			phi[z] = PhiBulk[z];
		}
	}
}


