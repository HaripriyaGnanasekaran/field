#ifndef SF_MONOMER2NDOxH
#define SF_MONOMER2NDOxH

#include "SF_Molecule.h"
#include "Lat2ndO.h"

///
class SF_Monomer2ndO: public SF_Molecule {
  public:
///
	Boolean force_set;
	//Boolean MayerSaupe;
	SF_Monomer2ndO(Text, SF_SegmentList*, Lat2ndO*, Input*);
///
	~SF_Monomer2ndO();
///
	void ComputePhi(void);
	void ReComputePhi(void);
	Vector GetLong(void);
	Vector GetShort(void);
	Vector GetBondOrientation( const DensityPart);
	MoleculeType GetMoleculeType() const;

///
	void GetLayerAnalysis(Output*);
	int gradients;
	int size;
  private:
///
	Lat2ndO* Lat;
///
	SF_MolSegment* Seg;
///
	void Matrix1();
///
	void Matrix2ndGen(const LatticeRange*,
					  const DensityPart,
					  const DensityPart);
//	void Matrix3thGen(const LatticeRange*,
//					  const DensityPart,
//					  const DensityPart,
//					  const DensityPart,
//					  const LatticeRange*,
///					  const DensityPart);
	void MatrixBulk(const DensityPart);
///
	void CopyBulkBoundaries(const DensityPart);
};

#endif
