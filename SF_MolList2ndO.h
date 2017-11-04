#ifndef SF_MOLLIST2NDOxH
#define SF_MOLLIST2NDOxH

#include "SF_MoleculeList.h"
#include "Lat2ndO.h"
#include "SF_Monomer2ndO.h"
#include "SF_Homopol2ndO.h"
#include "SF_Copol2ndO.h"
#include "SF_Copol2ndO_stiff_range.h"
#include "SF_Comb2ndO.h"

//#include "SF_Dend2ndO.h"
//include "SF Polydisp2ndO.h"
///
class SF_MolList2ndO: public SF_MoleculeList {
///
	Lat2ndO* Lat;
  public:
///
	SF_MolList2ndO(SF_SegmentList*, Lat2ndO*, Input*);
///
	~SF_MolList2ndO(void);
///
	void GetOutput(Output*) const;
///
	void ComputePhi(void);
	void ReComputePhi(void);
	Vector GetLong(void);
	Vector GetShort(void);
	MoleculeType GetMoleculeType() const;
///
	void ComputePhi(double);
///
	void ComputePhi(double,double);
///
	Boolean CheckMolecules(double) const;
///
	double GetChemicalPotential(int, DensityPart = total, int = 1) const;
///
	double GetChemicalPotential(const SF_Molecule*, DensityPart = total, int = 1) const;
///
	double GetNumMolTimesChemPot(int, DensityPart = total) const;
///
	double GetNumMolTimesChemPot(const SF_Molecule*, DensityPart = total) const;
///
	double GetChemicalPotentialLayer(int,int) const;
///
	double GetChemicalPotentialLayer(const SF_Molecule*,int) const;
///
	double GetFreeEnergy(void) const;
///
	double GetGrandPotential(void) const;
	double GetGrandPotentialExcess(void) const;
///
	double GetGrandPotential2(void) const; // less accurate, needed for reactions
///
	double GetExcessFreeEnergy(void) const;
	double GetSoumiFreeEnergy(void) const;
///
	double GetEntropy(SF_Molecule* = 0) const;
///
	double GetContactInteractions(void) const;
///
	double GetInternalFreeEnergy(void) const;
///
	double GetElectricInteractions(void) const;
///
	Vector GetFreeEnergyProfile(void) const;
///
	Vector GetGrandPotentialProfile(void) const;
///
	Vector GetGrandPotentialProfile2(void) const; // less accurate, needed for reactions
///
	Vector GetExcessFreeEnergyProfile(void) const;
///
	Vector GetEntropyProfile(SF_Molecule* = 0) const;
///
	Vector GetContactInteractionsProfile(void) const;
///
	Vector GetInternalFreeEnergyProfile(void) const;
///
	Vector GetElectricInteractionsProfile(void) const;
};

#endif

