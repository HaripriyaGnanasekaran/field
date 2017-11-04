#ifndef SF_ASYMDEND1STOxH
#define SF_ASYMDEND1STOxH

#include "SF_Molecule.h"
#include "Lat1stO.h"

///
class SF_AsymDend1stO: public SF_Molecule {
  public:
///
	SF_AsymDend1stO(Text, SF_SegmentList*, Lat1stO*, Input*);
///
	~SF_AsymDend1stO();
///
	void ComputePhi(void);
	void ReComputePhi(void);
	Vector GetLong(void);
	Vector GetShort(void);
	Vector GetBondOrientation( const DensityPart);
	MoleculeType GetMoleculeType() const;
///
	void GetLayerAnalysis(Output*);
	int numGenerations, NG;
	int *lengths;
	int *cumlengths_long, *cumlengths_short;
	int shortest;
	int longest;
	SF_SegmentBlock* SegBlock1;
	SF_SegmentBlock* SegBlock2;
	SF_SegmentBlock* SegBlock;




  private:
///
	Lat1stO* Lat;
///
	SF_MolSegment* Seg;
	Boolean force_set;
	Vector phi_short;
	Vector phi_long;
///
	int n;
///
	int numDiffSegments, NDS;
	double *PHI;
	double *PHI_;
///
	void Matrix1(const DensityPart);
	//void ComputePhi(int , Matrix ,Matrix, Vector , Vector , const DensityPart);
///
	void ComputePhi(int , Matrix ,Matrix, Vector , Vector , bool, bool, int, int, const DensityPart);
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
