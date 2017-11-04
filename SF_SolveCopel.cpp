#include "SF_SolveCopel.h"

SF_SolveCopel::SF_SolveCopel(Boolean compute_,
							 SF_ReactionList* ReactionQ_,
							 SF_MoleculeList* MolQ_,
							 SF_SegmentList* SegQ_,
							 Lattice* Lat_,
							 Input* MyInput_) 
	: SF_Solve(compute_,ReactionQ_,MolQ_,SegQ_,Lat_,MyInput_) {
	
}
SF_SolveCopel::~SF_SolveCopel() {
}
void
SF_SolveCopel::GetOutput(Output* Out) const {
	Out->PutText("newton",name,"iteration solver","copel");
	Out->PutInt("newton",name,"iterations",getiterations());
	Out->PutInt("newton",name,"iterationlimit",getiterationlimit());
	Out->PutBoolean("newton",name,"error occurred",errorOccurred);
	if (errorOccurred) {
		Out->SetErrorOccurred();
	}
	Out->PutReal("newton",name,"accuracy",getaccuracy());
	Out->PutReal("newton",name,"tolerance",gettolerance());
	Out->PutReal("newton",name,"deltamin",getdeltamin());
	Out->PutReal("newton",name,"deltamax",getdeltamax());
	Out->PutBoolean("newton",name,"pseudohessian",pseudohessian);
	Out->PutInt("newton",name,"number iter var",numItVar);
	Out->PutVector("newton",name,"iter var",x,numItVar);
	Vector funct(1,numItVar);
	GetResiduals(funct);
	Out->PutVector("newton",name,"residuals",funct,numItVar);
}
void
SF_SolveCopel::GetResiduals(Vector f) const {
	Vector phiTotal = SegQ->GetPhiTotal();
	Vector psi;
	if (charged) {
		if (potentialItVar) {
			psi = SegQ->ComputeElectricPotential(psi0);
		} else {
			psi = SegQ->ComputeElectricPotential();
		}
	}
	int j=1;
	SF_State* State;
	for (int z=1; z<=M; z++) {
		if (PotItVar(z)) {
			f[j++] = psi0[z] - psi[z] + phiTotal[z] - 1/phiTotal[z];
		}
		if (SegItVar(z)) {
			f[j++] = phiTotal[z] - 1/phiTotal[z];
			int count=NumStatesForItVar[1];
			for (int i=2; i<=numDiffItVar; i++) {
				count+=NumStatesForItVar[i];
				State = StateQ[count];
				f[j] = phiTotal[z] - 1/phiTotal[z];
				f[j] += SegQ->ChemInt(State,z) - SegQ->ChemIntBulk(State) + x[j];
				f[j] -= SegQ->ChemInt(StateQ[1],z) - SegQ->ChemIntBulk(StateQ[1]);
				j++;
			}
		}
	}
}
void 
SF_SolveCopel::residuals(double *const ff, double *const) {
	Vector f(ff, 1, numItVar);
	ComputePhi();
	GetResiduals(f);
}
void
SF_SolveCopel::UpdateSWF() {
	int i,j,k,z;
	j=1;
	for (z=1; z<=M; z++) {
		if (PotItVar(z)) {
			psi0[z] = x[j++];
		}
		if (SegItVar(z)) {
			int count = 0;
			for (k=1; k<=NumStatesForItVar[1]; k++) {
				count++;
				SWF[k][z] = exp(x[j]);
				if (charged) {
					SWF[k][z] *= exp(-(StateQ[k])->GetValence()*psi0[z]
					*ELEM_CHARGE/(BOLTZMANN*TEMPERATURE));
				}
			}
			j++;
			for (i=2; i<=numDiffItVar; i++) {
				for (k=1; k<=NumStatesForItVar[i]; k++) {
					count++;
					SWF[count][z] = exp(x[j-i+1]+x[j]);
					if (charged) {
						SWF[count][z] *= exp(-(StateQ[count])->GetValence()*psi0[z]
							*ELEM_CHARGE/(BOLTZMANN*TEMPERATURE));
					}
				}
				j++;
			}
		}
	}
	if (charged) {
		Lat->SetBoundaries(psi0);
	}
	SegQ->UpdatePhiBulk();
	SegQ->UpdateSWF();
}
void
SF_SolveCopel::UpdateItVar(Boolean){
	int i,j,z;
	j=1;
	SF_State* State;
	for (z=1; z<=M; z++) {
		if (PotItVar(z)) {
			x[j++] = psi0[z];
		}
		if (SegItVar(z)) {
			int count=NumStatesForItVar[1];
			State = StateQ[1];
			if (SWF[1][z] > 0) {
				x[j] = log(SWF[1][z]);
			}
			if (charged) {
				x[j] += State->GetValence()
				*psi0[z]*ELEM_CHARGE/(BOLTZMANN*TEMPERATURE);
			}
			j++;
			for (i=2; i<=numDiffItVar; i++) {
				count += NumStatesForItVar[i];
				State = StateQ[count];
				x[j] = -SegQ->ChemInt(State,z) + SegQ->ChemIntBulk(State) + SegQ->ChemInt(StateQ[1],z) - SegQ->ChemIntBulk(StateQ[1]);
				j++;
			}
		}
	}
	if (neutralizerPresent) {
		phiBulkNeutralizer = (MolQ->GetNeutralizer())->GetPhiBulk();
	}
	phiBulkSolvent = (MolQ->GetSolvent())->GetPhiBulk();
}
