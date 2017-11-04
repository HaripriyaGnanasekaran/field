#ifndef SF_MOLECULExH
#define SF_MOLECULExH

#include "SF_SegmentList.h"
#include "SF_MolStructure.h"

///
class SF_Molecule {
  public:
///
	SF_Molecule(Text, SF_SegmentList*, Lattice*, Input*);
///
	virtual ~SF_Molecule();

///
	virtual void ComputePhi(void) = 0;
	virtual void ReComputePhi(void) = 0;
	virtual Vector GetLong(void)=0;
	virtual Vector GetShort(void)=0;
	virtual Vector GetBondOrientation(const DensityPart)=0;
	virtual MoleculeType GetMoleculeType() const=0;
///
	Text GetName(void) const;
///
	virtual void GetOutput(Output*) const;
///
	virtual void GetLayerAnalysis(Output*) = 0;
///
	virtual void ContactNumberDistr(Output*,Boolean = false) const {};
///
	virtual void StateDistribution(Output*) const {};
///
	Vector GetPhi(const DensityPart) const;
///
	void CreatePhi(const DensityPart);
///
	void DeletePhi(const DensityPart);
///
	Array<DensityPart> GetDensityPartQ() const;
///
	SF_MolStructure* GetMolStructure() const;
///
	MolFreedom GetFreedom(void) const;
///
	void SetFreedom(MolFreedom);
///
	double GetChainLength(void) const;
///
	void SetPhiBulk(double);
///
	double GetPhiBulk(void) const;
///
	Vector GetPhiBulkBoundaries(void) const;
///
	Vector ComputePhiBoundaries(void) const;
///
	void SetTheta(double);
///
	double GetTheta(void) const;
///
	double ComputeTheta(void) const;
///
	double ComputeThetaExcess(void) const;
///
	double GetLnGN(void) const; /// !! returns log(GN)
	//Boolean MayerSaupeSet() const;
	double GetLnNormFree(void) const; /// !! returns log(Cb)
	double GetLnNormRestricted(void) const; /// !! returns log(Ct)
	void SetLnNormRestricted(double); /// !! set log(Ct)
	Vector GetRenorm(int) const;
///
	void SetRenorm(Vector, int);
///
	double GetThetaRenorm(void) const;
///
	Boolean Charged(void) const;
///
	double GetBulkCharge(void) const;
///
	double GetCharge(void) const;
///
	double GetChargeCoM(int sign) const;
///
	double GetNumberOfCharges(void) const;
///
	Boolean MuFixed(void) const;
///
	double GetMuFixed(void) const;
	double GetForce(void) const;
	Boolean GetGPU(void) const;
	Boolean GetIncludeTP(void) const;
	Array<Text> AliasNames;
	int numAlias;


///
  protected:
///
	Text name;
	double force;
///
	Input* MyInput;
///
	Lattice* LatBase;
///
	SF_SegmentList* SegQ;
///
	Text composition;
///
	double phiBulk;
///
	Vector phiBulkBoundaries;
///
	double theta;

///
	Boolean muFixed;
///
	double muFixedValue;
///
	double thetaRenorm;
///
	LatticeRange* restrictedRange;
////
	LatticeRange* MolRangeStiff;
///
	double phiRange;
///
	MolFreedom freedom;
///
	double lnGN; /// ! logarithm
	double lnCt; /// ! logarithm
	double lnCb; /// ! logarithm
	SF_MolStructure* Chain;


///
	LatticeRange* LatRange;
///
	LatticeRange* LatRangeRenorm;
///
	LatticeRange* LatRangeTrains;
///
	LatticeRange* LatRangeStartLoops;
///
	int numPhi;
///
	Array<Vector> PhiQ;
///
	Array<DensityPart> DensityPartQ;
///
	Array<Vector> RenormQ;
///
	int M;
	int gradients;
///
	int N;
///
	Boolean saveMemory,onGPU;
	Boolean InTP;
};

#endif
