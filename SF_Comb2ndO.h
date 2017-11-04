#ifndef SF_COMB2NDOxH
#define SF_COMB2NDOxH

#include "SF_Molecule.h"
#include "Lat2ndO.h"

///
class SF_Comb2ndO: public SF_Molecule {
  public:
///
	SF_Comb2ndO(Text, SF_SegmentList*, Lat2ndO*, Input*);
///
	~SF_Comb2ndO();
///
	void ComputePhi(void);
	void ReComputePhi(void);
	Vector GetLong(void);
	Vector GetShort(void);
	Vector GetBondOrientation( const DensityPart);
	MoleculeType GetMoleculeType() const;
///
	void GetLayerAnalysis(Output*);
	int numGenerations, NG, numArms;
	int N1,N2,N3,N4;

	SF_SegmentBlock* SegBlock1;
	SF_SegmentBlock* SegBlock2;
	SF_SegmentBlock* SegBlock3;
	SF_SegmentBlock* SegBlock4;

  private:
///
	Lat2ndO* Lat;
///
	SF_MolSegment* Seg;
	SF_MolSegment* Seg_old;
	Boolean force_set;
	Vector phi_short;
	Vector phi_long;
	Vector labda;
	double GetPf(SF_SegmentBlock*, int,int);
	double GetPf(SF_SegmentBlock*,SF_SegmentBlock*, int,int);
	double GetPf(SF_MolSegment*,SF_MolSegment*);
///
	int n;
	int size;
	int Mphi;
	double *PHI;
///
	int numDiffSegments, NDS;
	int side_type;

///
	void Matrix1(const DensityPart);
	//void ComputePhi(int , Matrix ,Matrix, Vector , Vector , const DensityPart);
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
