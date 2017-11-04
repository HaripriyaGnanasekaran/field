#ifndef SF_SYSTEMMCxH
#define SF_SYSTEMMCxH

#include "SF_System.h"

///
class SF_SystemMC : public SF_System {
  public:
///
	SF_SystemMC() {};
///
	SF_SystemMC(Input*, Text, Boolean);

///
	double random_d(void);
	//Boolean MayerSaupe;
	int random_int(int, int);
	double  random(double, double);
	bool CheckOverlap(int);
	bool MovePos(int,SF_Segment**,int);
	bool Off_Lattice_MovePos(int,SF_Segment**,int);
	double lastaverageloop;
	double lastvarloop;
	void Loopcounting(Output*,Lattice*,SF_Segment**,int*,int,int);

	void Go(Output*,Lattice*);
///
	//void SetInitialGuess(SF_System*,Lattice*);
///
	void GetOutput(Output*,Lattice*) const;
///
	//SF_Solve* GetSolve(void) const;
  protected:
///
	//SF_SegmentList* SegMCQ;
  ///
  	//SF_MoleculeList* MolMCQ;
  ///
  	//SF_ReactionList* ReactionMCQ;
///
	//void UpdateChi(void);

	//double timeNow;
	//double timeStep;
	//double timeLimit;
	//double timeStepIncr;
	//double error;
	//double maxError;
	//double minError;
	//double stopCriterion;
	//double stop;
	//int outputTimerLimit;

	int *x_old, *y_old, *z_old,*r_old, *x_new, *y_new, *z_new, *r_new, *r_last;
	double *rx_old, *ry_old, *rz_old, *rx_new, *ry_new, *rz_new;
	double *submask;
	bool *free;
	Array<Vector> rho_0;
	Vector F_0;
	Vector Omega_0;
	Vector Fpo_0;
	Vector PSI_0;
	vector<int> numposseg;
	Array<Vector> phi_0;

	int n_equi_steps;
	int numPos,spot,n_moves;
	Boolean MCtrain;
	int MCS;
	bool opt_mc,lat_moves,MC_F,Cluster;
	int mcs;
	int neighbor;
  	int d_max_pos;
  	double mc_temperature;
  	int SpotSize;
	double teleportfreq_Z,swapfreq;
  	double F_system_O,F_system_N;
  	int n_layers_x,n_layers_y,n_layers_z,jx,jy,jz;

	void GetMCOutput(Output*,Lattice*) const;
	void GetEnd_to_EndOutput(Output*,Lattice*) const;

	void SetRho(Lattice*,Array<Vector>)const;
	void SetF(Lattice*,Vector);
	void SetFpo(Lattice*,Vector);
	void SetOmega(Lattice*,Vector);
	void SetPsi(Lattice*,Vector);
	void DataMgr(int);
	void GetExtraOutput(Output*,Lattice*) const;
	void GetOutputParticleInBrush(Output*,Lattice*) const;


};

#endif










