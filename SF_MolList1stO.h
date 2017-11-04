#ifndef SF_MOLLIST1stOxH
#define SF_MOLLIST1stOxH

#include "SF_MoleculeList.h"
#include "Lat1stO.h"
#include "SF_Monomer1stO.h"
#include "SF_Homopol1stO.h"
#include "SF_HomoRing1stO.h"
#include "SF_Copol1stO.h"
#include "SF_Dend1stO.h"
#include "SF_AsymDend1stO.h"
#include "SF_Branched1stO.h"
#include "SF_Comb1stO.h"
//include "SF Polydisp1stO.h"
///
class SF_MolList1stO: public SF_MoleculeList {
///
	Lat1stO* Lat;
  public:
///
	SF_MolList1stO(SF_SegmentList*, Lat1stO*, Input*);
///
	~SF_MolList1stO(void);
///
	void GetOutput(Output*) const;
///
	void ComputePhi(void);
	void ReComputePhi(void);
	Vector GetLong(void);
	Vector GetShort(void);
	Vector GetBondOrientation(const DensityPart);
	MoleculeType GetMoleculeType() const;
	bool ring;

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
	double GetSoumiFreeEnergy(void) const;
///
	double GetGrandPotential(void) const;
	double GetGrandPotentialExcess(void) const;
///
	double GetGrandPotential2(void) const; // less accurate, needed for reactions
///
	double GetExcessFreeEnergy(void) const;
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

