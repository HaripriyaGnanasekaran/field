#ifndef SF_Branched1stOxH
#define SF_Branched1stOxH

#include <fenk.h>
#include "SF_Molecule.h"
#include "Lat1stO.h"
#include "misc.h"
#include "SF_SegmentBlock.h"

// a data type for an Array of Vectors
typedef Array<Vector> VecArray;

///
class SF_Branched1stO: public SF_Molecule {
  public:
///
//! The constructor for branched polymer
/*!
\param tree provides the topology of the polymer via links and nodes
*/
	SF_Branched1stO(Text name_, SF_SegmentList* SegQ_, Lat1stO* Lat_, Input* MyInput_, SF_LinkNodeTree tree);
//! just a testing constructor, will not work in practice
	SF_Branched1stO(Text name_, SF_SegmentList* SegQ_, Lat1stO* Lat_, Input* MyInput_);
///
//	the old constructor: SF_Branched1stO(Text, SF_SegmentList*, Lat1stO*, Input*);
///
	//! The dstrutor for SF_Branched1stO
	~SF_Branched1stO();
///
//! The computation of Phi and Gz
/*! Calls an appropriate version of Matrix*, depending on other parameters set by the user
*/
	void ComputePhi(void);
	void ReComputePhi(void);
	Vector GetLong(void);
	Vector GetShort(void);
	Vector GetBondOrientation(const DensityPart);
	MoleculeType GetMoleculeType() const;
///
	void GetLayerAnalysis(Output*);

// These functions have not been implemented:
//	void ContactNumberDistr(Output*,Boolean = false) const;
//	void StateDistribution(Output*) const;
  private:
///
	//! pointer to the beginning of the array of links provided to the constructor
	Array<SF_Link> Links;

//
	//! the size of the links array
	int n_links;
	Boolean force_set,GPUactivated;
//
	//! pointer to the beginning of the array of links provided to the constructor
	Array<SF_Node> Nodes;
//
	//! the size of THE NODES array
	int n_nodes;
//
	Lat1stO* Lat;
///
	SF_MolSegment* Seg;

///
	int n;
///
	int numDiffSegments;
//! Obsolete, will be removed as soon as possible.
/*! just to avoid compile errors while not all the functions that come form
the SF_Copol1stO have been re-implemented for the branched. In the Branched, Array<Vector> Gzs
is used instead.
*/
	Matrix Gi;
//! Matrix to store all the Gzs - for the Branched implemented as an Array of Vectors

	Array<VecArray> Gzs;

//! Vector for the inverse Gi
	Vector Gi_inv;
	int *Gs;
	double *Hphi, *Hg;
	double *Hgi_inv;  //pointers to host memory;
	double *Dphi, *Dg, *Dgi, *Dgi_inv, *Dgx; //pointers to device memory;
	int *Info;
///
	//! just for backward compatibility, every copolymer is apriori treated as asymmetric
	Boolean symmetric;
///
//! This is the classical (simple) version of Matrix1
/*! Uses the full Matrix which is demanding for memory but the algorithm is
much simpler than \ref Matrix1Long
*/
	void Matrix1(const DensityPart);
	void CudaMatrix1(const DensityPart);
///
	void Matrix2ndGen(const LatticeRange*,
					  const DensityPart,
					  const DensityPart);
///
//! This is the long version of \ref Matrix1
/*! Uses the reduced  Matrix which requires less memory but the algorithm is
much more complicated than \ref Matrix1
*/
	void Matrix1Long(const DensityPart);
///
	void Matrix2ndGenLong(const LatticeRange*,
						  const DensityPart,
						  const DensityPart);
	void MatrixBulk(const DensityPart);
///
	void CopyBulkBoundaries(const DensityPart);

/*!
Functions stolen from the simula code of sfbranch.sim
*/

	//! Computes Gz(s) of a node
	/*! Compute end probabilities of all arms but 'linkto'.
	\param nodenum is the number of the node for which Gz(s) is to be computed
	\param linkto it the number of link from which it was called
	*/
	void ComputeGzsNode(const int nodenum,const int linkto);
	void CudaComputeGzsNode(const int nodenum,const int linkto);
	//! Computes Gz(s) of a link
	/*!
	\param linknum is the number of the link for which Gz(s) is to be computed
	\param nodefrom is the number of the node from which it was called
	*/
	void ComputeGzsLink(const int linknum,const int nodefrom);
	void CudaComputeGzsLink(const int linknum,const int nodefrom);
	//! Computes Phi(s) of a node
	/*! Compute end probabilities of all arms but 'linkto'.
	\param nodenum is the number of the node for which Gz(s) is to be computed
	\param linkfrom is the number of link from which it was called
	*/
	void ComputePhiNode(int nodenum,int linkfrom, const DensityPart);
	void CudaComputePhiNode(int nodenum,int linkfrom, const DensityPart);
	//! Computes Phi(s) of a link
	/*!
	\param link is the number of the link for which Gz(s) is to be computed
	\param node is the number of the node from which it was called
	*/
	void ComputePhiLink(int linknum,int nodefrom, const DensityPart);
	void CudaComputePhiLink(int linknum,int nodefrom, const DensityPart);

	//! Propagator for the node connects G's of the different node segments and propagates
	void PropagateGzsNode(const int nodenum, const int linkto);
	void CudaPropagateGzsNode(const int nodenum, const int linkto);

	//! Add Phi1 to PhiOut, takes care of the overflow protection
	void AddPhi(const Vector Phi1, Vector PhiOut);
	/*! Debug functions follow */
	//! Print n Vectors, each in one column
	void PrintVec(const Vector V1,const char *v1);
	void PrintVec(const Vector V1, const Vector V2, const char *v1, const char *v2);
	void PrintVec(const Vector V1, const Vector V2, const Vector V3, const char *v1, const char *v2, const char *v3);
	void PrintVec(const Vector V1, const Vector V2, const Vector V3, const Vector V4, const char *v1, const char *v2, const char *v3, const char *v4);
};


#endif

