#include "SF_Copol2ndO_stiff_range.h"
#include "tools.h"

SF_Copol2ndO_stiff_range::SF_Copol2ndO_stiff_range(Text name_, SF_SegmentList* SegQ_,
						   Lat2ndO* Lat_, Input* MyInput_)
: SF_Molecule(name_, SegQ_, Lat_, MyInput_)
{
	if (Chain->GetMoleculeType() != copolymer) {
		Message(fatal,MyInput,
		"Programming error, trying to create a copolymer for molecule '" +
		name + "'");
	}
	symmetric = Chain->Symmetric();
	force_set = (!GetForce()==0);
	//if (MayerSaupeSet()) {
		//cout << "punt 1" << endl;
	//	Message(fatal,MyInput,"In SF_Copol2ndO_stiff_range: Mayer_Saupe not implemented: will be implemented upon request....");
	//}
	if (force_set) {
		Message(fatal,MyInput,
					"Force ensemble not implemented in Copol2stO_stiff_range");
	}
	if (GetGPU()) {
		Message(literal,MyInput,"For second order stiff range, the GPU card is not yet activated. Contact Frans Leermakers. Going classical instead");
	}
	numDiffSegments = Chain->GetNumDiffSegments();
	Lat = Lat_;
	if (freedom == secondGeneration) {
		CreatePhi(constrained);
		CreatePhi(unconstrained);
	} else {
		CreatePhi(total);
	}
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
SF_Copol2ndO_stiff_range::~SF_Copol2ndO_stiff_range() {
	if (freedom == secondGeneration) {
		DeletePhi(constrained);
		DeletePhi(unconstrained);
	} else {
		DeletePhi(total);
	}
}
void
SF_Copol2ndO_stiff_range::ReComputePhi() {
	Message(fatal,"ReComputePhi not implemented in SF_Copol2ndO_stiff_range");
}

Vector
SF_Copol2ndO_stiff_range::GetLong(){
	Message(fatal,"GetLong not implemented in SF_Copol2ndO_stiff_range");
	Vector x;
	return x;
}
Vector
SF_Copol2ndO_stiff_range::GetShort(){
	Message(fatal,"Getshort not implemented in SF_Copol2ndO_stiff_range");
	Vector x;
	return x;
}

Vector
SF_Copol2ndO_stiff_range::GetBondOrientation( const DensityPart DensPar) {
	Message(fatal,"GetBondOrientation not implemented in SF_Copol2ndO_stiff_range");
	Vector x;
	return x; 
};


MoleculeType
SF_Copol2ndO_stiff_range::GetMoleculeType() const {
	return Chain->GetMoleculeType();
}
void
SF_Copol2ndO_stiff_range::ComputePhi() {
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
	int z;
	if (symmetric) {
		if (freedom != secondGeneration) {
			if (!saveMemory) SymMatrix1(total);
			else SymMatrix1Long(total);
		} else {
			if (!saveMemory) SymMatrix2ndGen(LatRange,constrained,unconstrained);
			else SymMatrix2ndGenLong(LatRange,constrained,unconstrained);
		}
	} else {
		if (freedom != secondGeneration) {
			if (!saveMemory) Matrix1(total);
			else Matrix1Long(total);
		} else {
			if (!saveMemory) Matrix2ndGen(LatRange,constrained,unconstrained);
			else Matrix2ndGenLong(LatRange,constrained,unconstrained);
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
SF_Copol2ndO_stiff_range::GetLayerAnalysis(Output* Out) {
	//Message(literal,"getLayerAnalysis not implemented in 3d 2nd order");
}
void
SF_Copol2ndO_stiff_range::StateDistribution(Output* Out) const {
	// function to calculate the dissociation distribution function
 			Message(literal, "Unable to generate multistate "
 			"output in 3D 2nd order");
}

void
SF_Copol2ndO_stiff_range::ContactNumberDistr(Output* Out, Boolean partial) const {
	// function to calculate the endpoint distribution function as
	// a function of the number of surface contacts.
     Message(literal, "ContactNumberDistr output not generated in 3D 2nd order");
}


double SF_Copol2ndO_stiff_range::GetPf(int s1,int s2) {
	int NumPfs = SegQ->GetPf_count();
	Vector Pf = SegQ->GetPfVector();
	Array <Text> Pfx = SegQ->GetPfx();
	Array <Text> Pfy = SegQ->GetPfy();
	Text name1 = Chain->GetSegment(s1)->GetName();
	Text name2 = Chain->GetSegment(s2)->GetName();

	for (int i=1; i<=NumPfs; i++){
		if ((Pfx[i]==name1 && Pfy[i] == name2) || (Pfy[i]==name1 && Pfx[i] == name2)) return Pf[i];
	}
	return 1.0/5.0;
}

void
SF_Copol2ndO_stiff_range::Matrix1(const DensityPart DensPart) {
	int z,s,i,k;
	Vector phi,G;
	Matrix Gi(1,M*6,1,N);
	Vector Gi_inv(1,M*6);
	for (i=1; i<=numDiffSegments; i++) {
		Seg = Chain->GetDiffSegment(i);
		G = Seg->GetSWF();
		Lat->MakeSafe(G);
		phi = Seg->GetPhi(DensPart);
		Zero(phi,M);
		Lat->MakeSafe(phi);
	}
	G = Chain->GetSegment(1)->GetSWF();
	for (k=0; k<6; k++) {double *gi=&Gi[1+k*M][1], *g=&G[1]; cp(gi,g,M); }

	for (s=2; s<=N; s++) {
		G = Chain->GetSegment(s)->GetSWF();
		if (s==N) {Lat->PropagateF(Gi,G,s);} else {Lat->PropagateF(Gi,G,s,GetPf(s-1,s+1),true);}
	}
	Seg = Chain->GetSegment(N);
	G = Seg->GetSWF();
	for (k=0; k<6; k++) {double *gi_inv=&Gi_inv[1+k*M], *g=&G[1]; cp(gi_inv,g,M); }

	phi = Seg->GetPhi(DensPart);
	Lat->ConnectG(Gi_inv,Gi,N,phi);
	for (s=N-1; s>=1; s--) {
		Seg = Chain->GetSegment(s);
		G = Seg->GetSWF();
		phi = Seg->GetPhi(DensPart);
		if (s==N-1) {Lat->PropagateB(Gi_inv,G);} else {Lat->PropagateB(Gi_inv,G,GetPf(s,s+2),true);}
		Lat->ConnectG(Gi_inv,Gi,s,phi);
	}
	lnGN = Lat->ComputeLnGN(Gi_inv);
	for (i=1; i<=numDiffSegments; i++) {
		Seg = Chain->GetDiffSegment(i);
		G = Seg->GetSWF();
		phi = Seg->GetPhi(DensPart);
		Lat->CorrectDoubleCountG(phi,G);
	}
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
			Lat->NormPhiFree(phi,phiBulk/(6.0*N));
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
SF_Copol2ndO_stiff_range::Matrix1Long(const DensityPart DensPart) {
	int z,s,s0,t0,v0,t,rs1,i,k;

	Vector phi,G;
	Matrix Gi(1,M*6,1,n);
	Vector Gi_inv(1,M*6), Gs(1,M*6);

	for (i=1; i<=numDiffSegments; i++) {
		Seg = Chain->GetDiffSegment(i);
		G = Seg->GetSWF();
		Lat->MakeSafe(G);
		phi = Seg->GetPhi(DensPart);
		Zero(phi,M);
		Lat->MakeSafe(phi);
	}
	G = Chain->GetSegment(1)->GetSWF();
	for (k=0; k<6; k++) {double *gi=&Gi[1+k*M][1], *g=&G[1], *gs=&Gs[1+k*M]; cp(gi,g,M); cp(gs,g,M); }
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
		if (s==N) {Lat->PropagateF(Gs,G);} else {Lat->PropagateF(Gs,G,GetPf(s-1,s+1),true);}
		if ((t == t0+1 && t0 == v0)
		   || (t == t0+1 && ((n-t0)*(n-t0+1) >= N-1-2*(s0+t0)))
		   || (2*(n-t+s) >= N-1))
		  {double *gi=&Gi[1][t], *gs=&Gs[1]; cp(gi,gs,6*M);}
	}
	for (s=N; s>=1; s--) {
		G = Chain->GetSegment(s)->GetSWF();
		if (s == N) {
			for (k=0; k<6; k++) {double *gi_inv=&Gi_inv[1+k*M],*g=&G[1]; cp(gi_inv,g,M);}
		} else {
			if (s==N-1) {Lat->PropagateB(Gi_inv,G);} else {Lat->PropagateB(Gi_inv,G,GetPf(s,s+2),true);}
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
			double *gi=&Gi[1][t], *gs=&Gs[1]; cp(gs,gi,6*M);

			for (rs1=s0+t0+2; rs1<=s; rs1++) {
				t++;
				G = Chain->GetSegment(rs1)->GetSWF();
				if (rs1==N) {Lat->PropagateF(Gs,G);} else {Lat->PropagateF(Gs,G,GetPf(s-1,s+1),true);}
				if (t == t0+1 || s0+n == s) {
					double *gi=&Gi[1][t], *gs=&Gs[1]; cp(gi,gs,6*M);
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
			Lat->NormPhiFree(phi,phiBulk/(6.0*N));
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
void SF_Copol2ndO_stiff_range::Matrix2ndGen(const LatticeRange* Range,
								const DensityPart DensPart1,
								const DensityPart DensPart2) {
	Message(literal, "Matrix2nd Gen not implemented in 3D 2nd order");
}
void
SF_Copol2ndO_stiff_range::Matrix2ndGenLong(const LatticeRange* Range,
							   const DensityPart DensPart1,
							   const DensityPart DensPart2) {
	Message(literal, "Matrix2ndGenLong not implemented in 3D 2nd order");
}
void
SF_Copol2ndO_stiff_range::SymMatrix1(const DensityPart DensPart) {
	int z,s,s_inv,i,k;

	Vector phi,G;
	Matrix Gi(1,M*6,1,N/2);
	Vector Gi_inv(1,M*6);
	Vector Gc(1,M*6);
	for (i=1; i<=numDiffSegments; i++) {
		Seg = Chain->GetDiffSegment(i);
		G = Seg->GetSWF();
		Lat->MakeSafe(G);
		phi = Seg->GetPhi(DensPart);
		Zero(phi,M);
		Lat->MakeSafe(phi);
	}
	G = Chain->GetSegment(1)->GetSWF();
	for (k=0; k<6; k++) {double *gi=&Gi[1+k*M][1], *g=&G[1]; cp(gi,g,M); }

	for (s=2; s<=N/2; s++) {
		G = Chain->GetSegment(s)->GetSWF();
		Lat->PropagateF(Gi,G,s,GetPf(s-1,s+1),true);
	}
	s = N/2;
	double *gi_inv=&Gi_inv[1], *gi=&Gi[1][s], *gc=&Gc[1]; cp(gi_inv,gi,6*M); cp(gc,gi,6*M);
	if (N%2 == 1) {
		s = N/2 + 1;
		Seg = Chain->GetSegment(s);
		G = Seg->GetSWF();
		Lat->PropagateF(Gc,G,GetPf(s-1,s+1),true);
		phi = Seg->GetPhi(DensPart);
		Lat->PropagateB(Gi_inv,G,GetPf(s,s+2),true);
		Lat->ConnectG(Gi_inv,Gc,phi);
		//Lat->ConnectG(Gi_inv,Gi_inv,phi);
		Lat->NormPhiFree(phi,0.5); // correction needed for optimization
	}
	for (s=(N+3)/2; s<=N; s++) {
		s_inv = N-s+1;
		Seg = Chain->GetSegment(s);
		G = Seg->GetSWF();
		Lat->PropagateB(Gi_inv,G,GetPf(s,s+2),true);
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
		lnCt = log(theta) - lnGN -log(1.0*N);
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
			Lat->NormPhiFree(phi,2*phiBulk/(6.0*N));// factor 2 is correction needed for optimization
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
SF_Copol2ndO_stiff_range::SymMatrix1Long(const DensityPart DensPart) {
	int z,s;
	//int s_inv;
	int s0,t0,v0,t,rs1,i,k;

	Vector phi,G;
	Matrix Gi(1,M*6,1,n);
	Vector Gi_inv(1,M*6), Gs(1,M*6);
	Vector Gc(1,M*6);
	for (i=1; i<=numDiffSegments; i++) {
		Seg = Chain->GetDiffSegment(i);
		G = Seg->GetSWF();
		Lat->MakeSafe(G);
		phi = Seg->GetPhi(DensPart);
		Zero(phi,M);
		Lat->MakeSafe(phi);
	}
	G = Chain->GetSegment(1)->GetSWF();
	for (k=0; k<6; k++){ double *gi=&Gi[1+k*M][1], *gs=&Gs[1+k*M], *g=&G[1]; cp(gi,g,M); cp(gs,g,M);}

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
		Lat->PropagateF(Gs,G,GetPf(s-1,s+1),true);
		if ((t == t0+1 && t0 == v0)
		   || (t == t0+1 && ((n-t0)*(n-t0+1) >= N-1-2*(s0+t0)))
		   || (2*(n-t+s) >= N-1)) {
		    double *gi=&Gi[1][t], *gs=&Gs[1]; cp(gi,gs,6*M);
		}
	}
	double *gi_inv=&Gi_inv[1], *gs=&Gs[1], *gc=&Gc[1]; cp(gi_inv,gs,6*M); cp(gc,gs,6*M);
	if (N%2 == 1) {
		s = N/2 + 1;
		Seg = Chain->GetSegment(s);
		G = Seg->GetSWF();
		Lat->PropagateF(Gc,G,GetPf(s-1,s+1),true);
		phi = Seg->GetPhi(DensPart);
		Lat->PropagateB(Gi_inv,G,GetPf(s,s+2),true);
		Lat->ConnectG(Gi_inv,Gc,phi);
		Lat->NormPhiFree(phi,0.5);
	}
	for (s=(N+3)/2; s<=N; s++) {
		G = Chain->GetSegment(s)->GetSWF();
		Lat->PropagateB(Gi_inv,G,GetPf(s,s+2),true);
		//s_inv = N-s+1;
		t = N - s + 1 - s0;
		if (t == t0) {
			s0 += - n + t0;
			if (t0 == v0)
				s0 -= ((n - t0)*(n - t0 + 1))/2;
			t0 --;
			if (t0 < v0)
				v0 = t0;
			double *gs=&Gs[1], *gi=&Gi[1][t]; cp(gs,gi,6*M);

			for (rs1=s0+t0+2; rs1<=(N-s+1); rs1++) {
				t++;
				G = Chain->GetSegment(rs1)->GetSWF();
				Lat->PropagateF(Gs,G,GetPf(s-1,s+1),true);
				if (t == t0+1 || s0+n == N-s+1) {
					double *gs=&Gs[1], *gi=&Gi[1][t]; cp(gi,gs,6*M);
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
			Lat->NormPhiFree(phi,2*phiBulk/(6.0*N));// factor 2 is correction needed for optimization
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
void SF_Copol2ndO_stiff_range::SymMatrix2ndGen(const LatticeRange* Range,
								   const DensityPart DensPart1,
								   const DensityPart DensPart2) {
	Message(literal,"SymMatrix2ndGen not implemented in 3d 2nd order");
}
void
SF_Copol2ndO_stiff_range::SymMatrix2ndGenLong(const LatticeRange* Range,
								  const DensityPart DensPart1,
								  const DensityPart DensPart2) {
	Message(literal,"SymMatrix2ndGenLong not implemented in 3d 2nd order");
}
void
SF_Copol2ndO_stiff_range::MatrixBulk(const DensityPart DensPart) {
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
SF_Copol2ndO_stiff_range::CopyBulkBoundaries(const DensityPart DensPart) {
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

