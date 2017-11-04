#ifndef SF_HOMORING1STOxH
#define SF_HOMORING1STOxH

#include "SF_Molecule.h"
#include "Lat1stO.h"

///
class SF_HomoRing1stO: public SF_Molecule {
  public:
	SF_HomoRing1stO(Text, SF_SegmentList*, Lat1stO*, Input*);
	~SF_HomoRing1stO();
	void ComputePhi(void);
	void ReComputePhi(void);
	Vector GetLong(void);
	Vector GetShort(void);
	Vector GetBondOrientation( const DensityPart);
	MoleculeType GetMoleculeType() const;
	void GetLayerAnalysis(Output*);
	double norm_ref;

  private:

	Lat1stO* Lat;
	SF_MolSegment* Seg;

	int n;
	Boolean force_set,GPUactivated;

	void Matrix1(const DensityPart);
	void CudaMatrix1(const DensityPart);
	void Matrix1Long(const DensityPart);
	void CudaMatrix1Long(const DensityPart);


};

#endif
