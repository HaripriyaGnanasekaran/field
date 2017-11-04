#include "SF_MoleculeList.h"
SF_MoleculeList::SF_MoleculeList(SF_SegmentList* SegQ_, Input* MyInput_) {
	SegQ = SegQ_;
	MyInput = MyInput_;
}

SF_MoleculeList::~SF_MoleculeList() {
}
void
SF_MoleculeList::DeleteUnusedMolecules() {
	SF_Molecule* Mol;
	int i,j;
	for (i=1; i<=numMolecules; i++) {
		Mol = MoleculeQ[i];
		if (Mol->GetPhiBulk() == 0
		&& Mol->GetTheta() == 0
		&& Mol->GetFreedom() != neutralizer
		&& Mol->GetFreedom() != solvent
		&& Mol->GetFreedom() != rangeRestricted) {
			Message(warning,"Molecule '" + Mol->GetName() +
			"' is neither a solvent nor a neutralizer and it has neither a phibulk "
			"nor a theta defined, furthermore the freedom is not set to layer_restricted.\n"
			"It will be deleted for further calculation");
			delete Mol;
			Array<SF_Molecule*> MolQTemp(1,numMolecules-1);
			for (j=1; j<=numMolecules; j++) {
				if (j<i) MolQTemp[j] = MoleculeQ[j];
				if (j>i) MolQTemp[j-1] = MoleculeQ[j];
			}
			numMolecules--;
			i--;
			MoleculeQ = MolQTemp;
		}
	}
}
int
SF_MoleculeList::GetNumMolecules() const {
	return numMolecules;
}
Boolean
SF_MoleculeList::MoleculeDefined(Text name) const {
	for (int i=1; i<=numMolecules; i++) {
		if (*Copy((MoleculeQ[i])->GetName()) == *Copy(name)) {
			return true;
		}
	}
	return false;
}
SF_Molecule*
SF_MoleculeList::GetMolecule(int number) const {
	return MoleculeQ[number];
}
SF_Molecule*
SF_MoleculeList::GetMolecule(Text name) const {
	for (int i=1; i<=numMolecules; i++) {
		if (*Copy((MoleculeQ[i])->GetName()) == *Copy(name)) {
			return MoleculeQ[i];
		}
	}
	return NULL; // never get here
}
SF_Molecule*
SF_MoleculeList::GetSolvent() const {
	int i;
	SF_Molecule* Mol;
	SF_Molecule* Solvent=NULL;
	Boolean solventFound = false;
	for (i=1; i<=numMolecules; i++) {
		Mol = MoleculeQ[i];
		if (Mol->GetFreedom() == solvent) {
			if (solventFound) Message(fatal,MyInput,"Two solvent molecules defined");
			solventFound = true;
			Solvent = Mol;
		}
	}
	if (solventFound) {
		return Solvent;
	} else {
		Message(fatal,MyInput,"Set freedom to 'solvent' for one molecule");
		return NULL; // never get here;
	}
}
SF_Molecule*
SF_MoleculeList::GetNeutralizer() const {
	int i;
	SF_Molecule* Mol;
	SF_Molecule* Neutralizer=NULL;
	Boolean neutralizerFound = false;
	for (i=1; i<=numMolecules; i++) {
		Mol = MoleculeQ[i];
		if (Mol->GetFreedom() == neutralizer) {
			if (neutralizerFound) Message(fatal,MyInput,"Two molecules defined as neutralizer");
			neutralizerFound = true;
			Neutralizer = Mol;
		}
	}
	if (neutralizerFound) return Neutralizer;
	else return NULL;
}
