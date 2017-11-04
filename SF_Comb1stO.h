#ifndef SF_COMB1STOxH
#define SF_COMB1STOxH

#include "SF_Molecule.h"
#include "Lat1stO.h"


typedef Array<Vector> VecArray;
///
class SF_Comb1stO: public SF_Molecule {
  public:
///
	SF_Comb1stO(Text, SF_SegmentList*, Lat1stO*, Input*);
///
	~SF_Comb1stO();
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
	int N1,N2,N3,N4,NN;

	SF_SegmentBlock* SegBlock1;
	SF_SegmentBlock* SegBlock2;
	SF_SegmentBlock* SegBlock3;
	SF_SegmentBlock* SegBlock4;
	SF_SegmentBlock* SegBlockA1;
	SF_SegmentBlock* SegBlockA2;
	int *cumlengths_long, *cumlengths_short;
	int shortest,longest;

  private:
  	Array<SF_Link> Links;
  	Array<SF_Node> Nodes;
  	int n_links;
  	int n_nodes;
///
	Lat1stO* Lat;
///
	SF_MolSegment* Seg;

	Boolean force_set;
	Vector phi_short;
	Vector phi_long;

	Array<VecArray> Gzs;
	Vector gi_inv;
///
	int n;
///
	int numDiffSegments, NDS;
	int side_type;
	int lin,den,a_den,bra,com;
	Boolean symmetric;

///
	void MatrixLin(Matrix);
	void MatrixLin(Matrix,Vector,const DensityPart);
	void MatrixDen(Matrix);
	void MatrixDen(Matrix,Vector,const DensityPart);
	void MatrixADen(Matrix,Matrix);
	void ComputePhi(int, Matrix, Matrix, Vector, const DensityPart);

	void Matrix1(const DensityPart);
	void SymMatrix1(const DensityPart);
	//void ComputePhi(int , Matrix ,Matrix, Vector , Vector , const DensityPart);
///
	void Matrix2ndGen(const LatticeRange*,
			  const DensityPart,
			  const DensityPart);
///
	void MatrixBulk(const DensityPart);
///
	void CopyBulkBoundaries(const DensityPart);
	void ComputeGzsNode(const int nodenum,const int linkto);
	void ComputeGzsLink(const int linknum,const int nodefrom);
	void ComputePhiNode(int nodenum,int linkfrom, const DensityPart);
	void ComputePhiLink(int linknum,int nodefrom, const DensityPart);
	void PropagateGzsNode(const int nodenum, const int linkto);
	void AddPhi(const Vector Phi1, Vector PhiOut);

};

#endif
