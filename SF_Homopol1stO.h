#ifndef SF_HOMOPOL1STOxH
#define SF_HOMOPOL1STOxH

#include "SF_Molecule.h"
#include "Lat1stO.h"

///
class SF_Homopol1stO: public SF_Molecule {
  public:
///
	SF_Homopol1stO(Text, SF_SegmentList*, Lat1stO*, Input*);
///
	~SF_Homopol1stO();
///
	void ComputePhi(void);
	void ReComputePhi(void);
	Vector GetLong(void);
	Vector GetShort(void);
	Vector GetBondOrientation( const DensityPart);
	MoleculeType GetMoleculeType() const;
///
	void GetLayerAnalysis(Output*);
///
//	void GetSegmentAnalysis(Output*);
  private:
///
	Lat1stO* Lat;
///
	SF_MolSegment* Seg;
///
	int n;
	Boolean force_set,GPUactivated;
///
	void Matrix1(const DensityPart);
	void CudaMatrix1(const DensityPart);
///
	void Matrix1Long(const DensityPart);
	void CudaMatrix1Long(const DensityPart);
///
	void Matrix2ndGen(const LatticeRange*,
					  const DensityPart,
					  const DensityPart);
///
	void Matrix2ndGenLong(const LatticeRange*,
						  const DensityPart,
						  const DensityPart);
///
	void Matrix3rdGen(const LatticeRange*,
					  const DensityPart,
					  const DensityPart,
					  const DensityPart,
					  const LatticeRange*);
///
	void Matrix3rdGenLong(const LatticeRange*,
						  const DensityPart,
						  const DensityPart,
						  const DensityPart,
						  const LatticeRange*);
///
	void MatrixBulk(const DensityPart);
///
	void CopyBulkBoundaries(const DensityPart);
};

#endif
