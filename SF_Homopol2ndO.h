#ifndef SF_HOMOPOL2NDOxH
#define SF_HOMOPOL2NDOxH

#include "SF_Molecule.h"
#include "Lat2ndO.h"

///
class SF_Homopol2ndO: public SF_Molecule {
  public:
///
	SF_Homopol2ndO(Text, SF_SegmentList*, Lat2ndO*, Input*);
///
	~SF_Homopol2ndO();
///
	void ComputePhi(void);
	void ReComputePhi(void);
	Vector GetLong(void);
	Vector GetShort(void);
	Vector GetBondOrientation( const DensityPart);
	MoleculeType GetMoleculeType() const;
	//SF_SegmentList* SegQ;
///
	void GetLayerAnalysis(Output*);
	double GetPf() ;

///
//	void GetSegmentAnalysis(Output*);
  private:

///
	Lat2ndO* Lat;
///
	SF_MolSegment* Seg;

///
	int n;
	int size;
	int Mphi;
	double *PHI;
	Boolean force_set;
	//Boolean MayerSaupe;
///
	void Matrix1(const DensityPart);
///
	void Matrix1Long(const DensityPart);
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
