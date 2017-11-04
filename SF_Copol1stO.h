#ifndef SF_COPOL1STOxH
#define SF_COPOL1STOxH

#include "SF_Molecule.h"
#include "Lat1stO.h"
#include "misc.h"

///
class SF_Copol1stO: public SF_Molecule {
  public:
///
	SF_Copol1stO(Text, SF_SegmentList*, Lat1stO*, Input*);
///
	~SF_Copol1stO();
///
	void ComputePhi(void);
	void ReComputePhi(void);
	Vector GetLong(void);
	Vector GetShort(void);
	Vector GetBondOrientation( const DensityPart);
	MoleculeType GetMoleculeType() const;
///
	void GetLayerAnalysis(Output*);
	void ContactNumberDistr(Output*,Boolean = false) const;
	void StateDistribution(Output*) const;
  private:
///
	Lat1stO* Lat;
///
	SF_MolSegment* Seg;
///
	Boolean symmetric;
	Boolean force_set,GPUactivated;
///
	int n;
///
	int numDiffSegments;
///
	void Matrix1(const DensityPart);
	void Matrix2(const DensityPart);
       void CudaMatrix1(const DensityPart);
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
	void Matrix1Long(const DensityPart);
	void CudaMatrix1Long(const DensityPart);
///
	void Matrix2ndGenLong(const LatticeRange*,
						  const DensityPart,
						  const DensityPart);
//	void Matrix3thGenLong(const LatticeRange*,
//						  const DensityPart,
//						  const DensityPart,
//						  const DensityPart,
//						  const LatticeRange*,
///						  const DensityPart);
	void SymMatrix1(const DensityPart);
	void CudaSymMatrix1(const DensityPart);
///
	void SymMatrix2ndGen(const LatticeRange*,
						 const DensityPart,
						 const DensityPart);
//	void SymMatrix3thGen(const LatticeRange*,
//						 const DensityPart,
//						 const DensityPart,
//						 const DensityPart,
//						 const LatticeRange*,
///						 const DensityPart);
	void SymMatrix1Long(const DensityPart);
	void CudaSymMatrix1Long(const DensityPart);

///
	void SymMatrix2ndGenLong(const LatticeRange*,
							 const DensityPart,
							 const DensityPart);
//	void SymMatrix3thGenLong(const LatticeRange*,
//							 const DensityPart,
//							 const DensityPart,
//							 const DensityPart,
//							 const LatticeRange*,
///							 const DensityPart);
	void MatrixBulk(const DensityPart);
///
	void CopyBulkBoundaries(const DensityPart);
};

#endif
