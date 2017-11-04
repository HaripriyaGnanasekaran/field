#include <iostream>
#include "SF_Solve_Johan.h"
SF_Solve_Johan::SF_Solve_Johan(Boolean compute_,
				   SF_ReactionList* ReactionQ_,
				   SF_MoleculeList* MolQ_,
				   SF_SegmentList* SegQ_,
				   Lattice* Lat_,
				   Input* MyInput_)
				   : SF_Solve(compute_,ReactionQ_,MolQ_,SegQ_,Lat_,MyInput_)
				   {

}
SF_Solve_Johan::~SF_Solve_Johan() {
}



void
SF_Solve_Johan::GetResiduals(Vector f) const {
	int i,j,z,count;
	//int count2;
	SF_State* State;
//	SF_State* State2;
	Vector phiTotal = SegQ->GetPhiTotal();
	Vector uPrime(1,numDiffItVar);
	double uPrimeAv;
	Vector psi;
	if (charged) {
		if (potentialItVar) {
			psi = SegQ->ComputeElectricPotential(psi0);
		} else {
			psi = SegQ->ComputeElectricPotential();
		}
	}
	count=0;
	//count2=0;
	for (i=1; i<=numLiquidStates; i++) {
	//for (i=1; i<=numDiffItVar; i++) {
		//count+=NumStatesForItVar[i];
		//State = StateQ[count];
		State = StateQ[i];
		SegQ->PutSides(State);
	}

        j=1;
	for (z=1; z<=M; z++) {
		if (PotItVar(z)) {
			f[j++] = psi0[z] - psi[z];
		}
		if (SegItVar(z)) {
			count=0; //count2=0;
			uPrimeAv=0;
			for (i=1; i<=numDiffItVar; i++) {
				count+=NumStatesForItVar[i];
				State = StateQ[count];
				uPrime[i] = -log(SWF[count][z]);
				uPrime[i] -= SegQ->ChemInt(State,z)/phiTotal[z];
				uPrime[i] += SegQ->ChemIntBulk(State);
				if (charged) {
					uPrime[i] -= State->GetValence()*psi0[z]*
						ELEM_CHARGE/(BOLTZMANN*TEMPERATURE);
					if (diffEpsilon) {
						uPrime[i] += 0.5*State->GetEpsilon()*ESquared[z];
					}
				}
				if (State->ExternalPotential()) {
					uPrime[i] -= (State->GetExternalPotential())[z];
				}
				uPrimeAv += uPrime[i];

			}
			uPrimeAv /= numDiffItVar;

			for (i=1; i<=numDiffItVar; i++){
				f[j++] = - 1.0 + 1.0/phiTotal[z]+ uPrime[i] - uPrimeAv;
				//f[j++] = 0.43*(1-phiTotal[z]) + (uPrime[i] - uPrimeAv); //perhaps better in terms of convergence
			}
		}
	}
}

