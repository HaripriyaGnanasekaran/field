#ifndef SF_MONOMER1STOxH
#define SF_MONOMER1STOxH

#include "SF_Molecule.h"
#include "Lat1stO.h"

///
class SF_Monomer1stO: public SF_Molecule {
  public:
///
	//Boolean force_set;
	SF_Monomer1stO(Text, SF_SegmentList*, Lat1stO*, Input*);
///
	~SF_Monomer1stO();
///
	void ComputePhi(void);
	void ReComputePhi(void);
	Vector GetLong(void);
	Vector GetShort(void);
	Vector GetBondOrientation( const DensityPart); 
	MoleculeType GetMoleculeType() const;
///
	void GetLayerAnalysis(Output*);
  private:
///
	Lat1stO* Lat;
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
