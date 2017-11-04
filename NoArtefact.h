#ifndef NOARTEFACTxH
#define NOARTEFACTxH

#include "Lat2DFlat1stOSten.h"
#include "Lat1DSphere1stO.h"
#include "Lat2DCyl1stO.h"
#include "Lat1DSphere1stOS.h"
#include "Lat2DCyl1stOS.h"
#include "Lat2DCyl1stOSten.h"
#include "SF_MolList1stO.h"
#include "SF_ReactionList.h"
#include "SF_SolveCopel.h"

///
class NoArtefact {
  public:
  ///
	NoArtefact(Lattice*, SF_MoleculeList*, SF_SegmentList*,SF_Solve*, double);
  ///
	void Go(void);
  ///
  protected:
///
	Lattice* Lat;
///
	SF_MoleculeList* MolQ;
///
	SF_SegmentList* SegQ;
///
	SF_Solve* Solve;
/// 
	double tol;
///
	double Delta(double);
///
	void Jump(int);
};
#endif
