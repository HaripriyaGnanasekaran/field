#ifndef SF_MOLSTATExH
#define SF_MOLSTATExH

#include "SF_Box.h"
#include "Lattice.h"

///
class SF_MolState: public Link {
  public:
    ///
    SF_MolState(Text,
    			Text,
    			double,
    			double,
    			double,
    			SegmentFreedom,
    			LatticeRange*,
    			Lattice*);
    ///
    virtual ~SF_MolState();

    ///
    Text GetName(void) const;
    ///
    Text GetMolName(void) const;
    ///
    double GetValence(void) const;
    ///
    double GetEpsilon(void) const;
    ///
    double GetInternalFreeEnergy(void) const;
    ///
    void SetInternalFreeEnergy(double);
    ///
    double GetAlphaBulk(void) const;
    ///
    void SetAlphaBulk(double);
 ///
 	Vector GetAlphaBulkBoundaries(void) const;
///
	void SetAlphaBulkBoundaries(Vector);
    ///
    double GetPhiBulk(void) const;
    ///
    void SetPhiBulk(double);
 ///
 	Vector GetPhiBulkBoundaries(void) const;
///
	void SetPhiBulkBoundaries(Vector);
///
	double GetPhiRef(void) const;
    ///
    void SetPhiRef(double);
    ///
    SegmentFreedom GetFreedom(void) const;
    ///
    LatticeRange* GetLatRange(void) const;
    ///
    void SetSWF(Vector);
    ///
    Vector GetSWF(void) const;
    ///
    void DelPos(int,int,int);
    void ClearAllPos(void);
    void UpdatePos(int,int,int);
    void UpdatePos(double,double,double,double*);
    void UpdateSWF(void);


    ///
    void SetExternalPotential(Vector);
    ///
    Vector GetExternalPotential(void);
    ///
    Boolean ExternalPotential(void);
    ///
    Vector GetPhi(DensityPart) const;
 ///
	double GetTheta(DensityPart) const;
   ///
    void CreatePhi(DensityPart);
    ///
    void DeletePhi(DensityPart);
///
	Array<DensityPart> GetDensityPartQ() const;
  private:
    ///
    Text name;
    ///
    Text molName;
    ///
    SegmentFreedom freedom;
    ///
    LatticeRange* LatRange;
    ///
    double valence;
    ///
    double epsilon;
    ///
    double internalFreeEnergy;
    ///
    double alphaBulk;
    ///
    Vector alphaBulkBoundaries;
    ///
    double phiBulk;
    ///
    Vector phiBulkBoundaries;
    ///
    double phiRef;
    ///
    Vector alpha;
    ///
    int numPhi;
    ///
    Array<Vector> PhiQ;
    ///
    Array<DensityPart> DensityPartQ;
    ///
    Vector swf;
    ///
    Vector swfFull;
    ///
    Boolean extPot;
    ///
    Vector externalPotential;
    ///
    Lattice* Lat;
};

#endif
