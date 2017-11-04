#ifndef SF_SYSTEMLDxH
#define SF_SYSTEMLDxH

#include "SF_System.h"

///
class SF_SystemLD : public SF_System {
  public:
///
	SF_SystemLD() {};
///
	SF_SystemLD(Input*, Text, Boolean);

///
	//Boolean MayerSaupe;
	double random_d(void);
	int random_int(int, int);
	double  random(double, double);
	double  normal(double, double);
	bool CheckOverlap(void);
	bool MovePos(SF_Segment**);

	void Go(Output*,Lattice*);
///
	void GetOutput(Output*,Lattice*) const;

  protected:
///
	//SF_SegmentList* SegNewQ;
  ///
//SF_MoleculeList* MolNewQ;
  ///
  	//SF_ReactionList* ReactionNewQ;

	int *x_old, *y_old, *z_old,*r_old, *x_new, *y_new, *z_new, *r_new, *r_last;
	double *rx_old, *ry_old, *rz_old, *rx_new, *ry_new, *rz_new, *F_old, *F_new;
	double *submask;
	bool *free;
	Array<Vector> rho_0;
	Vector F_0;
	Vector Omega_0;
	Vector Fpo_0;
	Vector PSI_0;
	bool LD_F;
	int numPos,spot;
	int LDT;
	int ldt;
  	double mobility;
  	double F_system_O,F_system_N;
  	int n_layers_x,n_layers_y,n_layers_z,jx,jy,jz;

	void GetLDOutput(Output*,Lattice*) const;
	void GetEnd_to_EndOutput(Output*,Lattice*) const;

	void SetRho(Lattice*,Array<Vector>)const;
	void SetF(Lattice*,Vector);
	void SetFpo(Lattice*,Vector);
	void SetOmega(Lattice*,Vector);
	void SetPsi(Lattice*,Vector);

};

#endif










