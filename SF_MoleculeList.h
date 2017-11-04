#ifndef SF_MOLECULELISTxH
#define SF_MOLECULELISTxH

#include "SF_Molecule.h"

///
class SF_MoleculeList {
  public:
///
	SF_MoleculeList(SF_SegmentList*, Input*);
///
	virtual ~SF_MoleculeList(void);
///
	virtual void GetOutput(Output*) const = 0;
///
	void DeleteUnusedMolecules(void);
///
	int GetNumMolecules(void) const;
///
	Boolean MoleculeDefined(Text) const;
///
	SF_Molecule* GetMolecule(int) const;
///
	SF_Molecule* GetMolecule(Text) const;
///
	SF_Molecule* GetSolvent(void) const;
///
	SF_Molecule* GetNeutralizer(void) const; /// returns NULL if no neutralizer is defined
	virtual Boolean CheckMolecules(double) const = 0;
///
	virtual double GetChemicalPotential(int, DensityPart = total, int = 1) const = 0;
///
	virtual double GetChemicalPotential(const SF_Molecule*, DensityPart = total, int = 1) const = 0;
///
	virtual double GetNumMolTimesChemPot(int, DensityPart = total) const = 0;
///
	virtual double GetNumMolTimesChemPot(const SF_Molecule*, DensityPart = total) const = 0;
///
	virtual double GetChemicalPotentialLayer(int,int) const = 0;
///
	virtual double GetChemicalPotentialLayer(const SF_Molecule*,int) const = 0;
///
	virtual double GetFreeEnergy(void) const = 0;
///
	virtual	double GetSoumiFreeEnergy(void) const = 0;
///
	virtual double GetGrandPotential(void) const = 0;
	virtual double GetGrandPotentialExcess(void) const = 0;
///
	virtual double GetGrandPotential2(void) const = 0;// less accurate, needed for reactions
///
	virtual double GetExcessFreeEnergy(void) const = 0;
///
	virtual double GetEntropy(SF_Molecule* = 0) const = 0;
///
	virtual double GetContactInteractions(void) const = 0;
///
	virtual double GetInternalFreeEnergy(void) const = 0;
///
	virtual double GetElectricInteractions(void) const = 0;
///
	virtual Vector GetFreeEnergyProfile(void) const = 0;
///
	virtual Vector GetGrandPotentialProfile(void) const = 0;
///
	virtual Vector GetGrandPotentialProfile2(void) const = 0;// less accurate, needed for reactions
///
	virtual Vector GetExcessFreeEnergyProfile(void) const = 0;
///
	virtual Vector GetEntropyProfile(SF_Molecule* = 0) const = 0;
///
	virtual Vector GetContactInteractionsProfile(void) const = 0;
///
	virtual Vector GetInternalFreeEnergyProfile(void) const = 0;
///
	virtual Vector GetElectricInteractionsProfile(void) const = 0;
///
	virtual void ComputePhi(void) = 0;
	virtual void ReComputePhi(void) = 0;
///
	virtual void ComputePhi(double) = 0;
///
	virtual void ComputePhi(double, double) = 0;
  protected:
///
	int numMolecules;
///
	Array<SF_Molecule*> MoleculeQ;
///
	SF_SegmentList* SegQ;
///
	Input *MyInput;
};

#endif
