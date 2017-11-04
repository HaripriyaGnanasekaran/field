#ifndef SF_COPOL2NDO_STIFF_RANGExH
#define SF_COPOL2NDO_STIFF_RANGExH

#include "SF_Molecule.h"
#include "Lat2ndO.h"
#include "misc.h"


///
class SF_Copol2ndO_stiff_range: public SF_Molecule {
  public:
///
	SF_Copol2ndO_stiff_range(Text, SF_SegmentList*, Lat2ndO*, Input*);
///
	~SF_Copol2ndO_stiff_range();
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
	Lat2ndO* Lat;

///
	SF_MolSegment* Seg;
///
	Boolean symmetric;
	Boolean force_set;
///
	int n;
///
	int numDiffSegments;
	void SetPf();
	double GetPf(int, int);
///
	void Matrix1(const DensityPart);
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
