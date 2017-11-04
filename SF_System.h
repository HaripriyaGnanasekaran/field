#ifndef SF_SYSTEMxH
#define SF_SYSTEMxH

#include "Lat2DFlat1stOSten.h"
#include "Lat1DSphere1stO.h"
#include "Lat2DCyl1stO.h"
#include "Lat1DSphere1stOS.h"
#include "Lat2DCyl1stOS.h"
#include "Lat2DCyl1stOSten.h"
#include "Lat2DCyl1stOStenS.h"
#include "Lat1DCyl2ndO.h"
#include "Lat2DCyl2ndO.h"
#include "Lat1DSphere2ndO.h"
#include "Lat3D1stO.h"
#include "Lat3DFCC1stO.h"
#include "Lat3DHEX1stO.h"
#include "Lat3D2ndO.h"
#include "Lat2D2ndO.h"
#include "Lat1D2ndO.h"
#include "Lat3D.h"
#include "SF_MolList1stO.h"
#include "SF_MolList2ndO.h"
#include "SF_ReactionList.h"
#include "SF_SolveCopel.h"
#include "SF_SolvePikar.h"
#include "SF_SolveCG.h"
#include "SF_SolvePhiU.h"
#include "SF_Solve_Johan.h"
#include "NoArtefact.h"
#include "misc.h"

///
class SF_System : public SFNewton {
  public:

///
	static const double ROOMTEMPERATURE;

///
	SF_System() {};
	//Boolean MayerSaupe;
///
	SF_System(Input*, Text, Boolean);
///
	virtual ~SF_System(void);
///
	virtual void Go(Output*, Lattice* );
///
	virtual void SetInitialGuess(SF_System*, Lattice*);
///
	Lattice* GetLattice(void) const;
///
	SF_SegmentList* GetSegmentList(void) const;
///
	SF_MoleculeList* GetMoleculeList(void) const;
///

///
	virtual SF_Solve* GetSolve(void) const;
  protected:
  ///


	bool Cuda_enabled;
  	Text name;
  ///
	Input* MyInput;
  ///
  	Lat1stO* Lat1st;
	Lat2ndO* Lat2nd;
  ///
  	SF_SegmentList* SegQ;
  ///
  	SF_MoleculeList* MolQ;
  ///
  	SF_ReactionList* ReactionQ;
  ///
  	SF_Solve* Solve;
  ///
  	Boolean compute;
  ///
	MatrixApproach approach;
  ///
  	Boolean overflowProtection;
  ///
  	Boolean changeChiWithTemperature;
  ///
  	double temperature;
  ///
  	Boolean graftedMolecules;
  ///
  	Boolean superIterate;
  ///
  	SF_Molecule* SuperMol;
  ///
  	enum supFunction {grandPotential, phiBulk, chemPot, theta, GPE, ETS, GibbsExc, P_Laplace};
  ///
  	supFunction superFunction;
  ///
  	double superFunctionValue;
  ///
  	Boolean iterateLatticeArtefact;
  ///
  	double artefactTolerance;
  ///
  	NoArtefact* ArtefactIter;
  ///
	void residuals(double *const, double *const);
  ///
  	Lat1stO* NewLat1stO(Text) const;
  	double GetGibbsExcess() ;

	Lat2ndO* NewLat2ndO(Text) const;
  ///
	SF_Solve* NewSolve(SF_ReactionList*,
					   SF_MoleculeList*,
					   SF_SegmentList*,
					   Lattice*) const;
///
	void UpdateChi(void);
///
	void AddChemPotGraft(Vector,Lattice*) const;
	bool AddChemPotGraftForHelmholtz(Vector,Lattice*) const;
///
	virtual void GetOutput(Output*,Lattice*) const;
///
	void GetExtraOutput(Output*,Lattice*) const;
///
	void GetOutputJoanne(Output*,Lattice*) const;
///
	void GetOutputSecondGen(Output*,Lattice*) const;
///
	void GetOutputPressure(Output*,Lattice*) const;
///
	void GetOutputHairy(Output*,Lattice*) const;
///
	void GetOutputVesicle(Output*,Lattice*) const;
///
	void GetOutputTwoPhase(Output*, Text,Lattice*) const;
///
	double FindGibbsDevPlane(Vector, Text, Lattice*) const;
///
	double CalcGammaGibbs(double, Text,Lattice*) const;
///
	void GetOutputEmulsion(Output*,Lattice*) const;
///
	void GetOutputIsolated(Output*,Lattice*) const;
	void GetOutputSpinodal(Output*,Lattice*) const;
	void GetOutputDepletion(Output*,Lattice*) const;
	void GetOutputPiA(Output*,Lattice*) const;
	void GetOutputBend(Output*,Lattice*) const;
	void GetOutputMultiState(Output*,Lattice*) const;
	void GetCeciliaOutput(Output*,Lattice*) const;
	void GetEmiliaOutput(Output*,Lattice*) const;
	void GetSergioOutput(Output*,Lattice*) const;
	void GetDSMOutput(Output*,Lattice*) const;
	void GetEgorovOutput(Output*,Lattice*) const;
	void GetEscapeOutput(Output*,Lattice*) const;
	void GetBrush3GOutput(Output*,Lattice*) const;
	void GetModuliOutput(Output*,Lattice*) const;
	void GetPoreOutput(Output*,Lattice*) const;
	void GetMoments2GOutput(Output*,Lattice*) const;
	void GetMoments3GOutput(Output*,Lattice*) const;
	void GetProbeSizeOutput(Output*,Lattice*) const;
	void GetSoumiOutput(Output*,Lattice*) const;
	void GetLorealOutput(Output*,Lattice*) const;
	void GetFransOutput(Output*,Lattice*) const;
	void GetPhiZ(Output*,Lattice*) const;
	double GetFreeEnergyPo(Lattice*);
	void GetDendronOverlap(Output* ,Lattice* ) const;
	void GrandPotentialLeftRight(Output* ,Lattice* ) const;
	void Katya(Output* ,Lattice* ) const;
	void Sabine(Output* ,Lattice* ) const;
	void Johan(Output* ,Lattice* ) const;
	void Boris(Output* ,Lattice* ) const;
	void Zeta(Output* ,Lattice* ) const;
	void CoExistBilayers(Output* ,Lattice* ) const;
	void InterfacialWidth(Output* ,Lattice* ) const;
	void LineTension(Output* ,Lattice* ) const;
	int numMCseg;

	//void GetEnd_to_EndOutput(Output*,Lattice*) const;
	//void SetRho(Lattice*,Array<Vector>)const;
	//void SetF(Lattice*,Vector);
	//void SetFpo(Lattice*,Vector);
	//void SetOmega(Lattice*,Vector);
	//void SetPsi(Lattice*,Vector);
	//void GetMCOutput(Output*,Lattice*) const;
	//bool CheckOverlap();
	//bool MovePos(int,SF_Segment**,int);
	//bool Off_Lattice_MovePos(int,SF_Segment**,int);
	//Boolean MCtrain;
	//int MCS;
	//bool simple,risky,opt_mc,lat_moves;
	//int mcs;
  	//int d_max_pos;
  	//double mc_temperature;
    //double F_system_O,F_system_N;
    //int *x_old, *y_old, *z_old,*r_old, *x_new, *y_new, *z_new, *r_new, *r_last;
	//double *rx_old, *ry_old, *rz_old, *rx_new, *ry_new, *rz_new;
	//double *submask;
	//bool *free;
	//Array<Vector> rho_0;
	//Vector F_0;
	//Vector Omega_0;
	//Vector Fpo_0;
	//Vector PSI_0;
	//int numPos,spot,n_moves;
	//int n_layers_x,n_layers_y,n_layers_z,jx,jy,jz;
};

#endif
