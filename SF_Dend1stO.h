#ifndef SF_DEND1STOxH
#define SF_DEND1STOxH

#include "SF_Molecule.h"
#include "Lat1stO.h"

///
class SF_Dend1stO: public SF_Molecule {
  public:
///
	SF_Dend1stO(Text, SF_SegmentList*, Lat1stO*, Input*);
///
	~SF_Dend1stO();
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
	Boolean force_set;
///
	int n;
///
	int numDiffSegments;
///
	void Matrix1(const DensityPart);
///
	void Matrix2ndGen(const LatticeRange*,
			  const DensityPart,
			  const DensityPart);
///
	void MatrixBulk(const DensityPart);
///
	void CopyBulkBoundaries(const DensityPart);
};

#endif
