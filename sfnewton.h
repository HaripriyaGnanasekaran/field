#ifndef SFNEWTONxH
#define SFNEWTONxH

#define _GNU_SOURCE 1
#include <iostream>
#include <iomanip>
#include <string.h>
#include <math.h>
#include <limits>
#include <time.h>
#ifdef USE_FENV
#include <fenv.h>
#endif
#ifdef __SVR4
#include <ieeefp.h>
#elif _WIN32
#define finite(x) (_finite(x))
#ifndef isnan
#define isnan(x) (_isnan(x))
#endif
#endif

#include "newttool.h"
#include <deque>
#include <vector>
#include <fenk/array.h>
//#include <time.h>
using namespace std;

class SFNewton {

/*      class for
        unconstrained minimization,
        systems of nonlinear algebraic equations,
        curve fitting.

References:

 (1)    "The State of the Art in Numerical Analysis"
        Edited by D.Jacobs. (1977) Academic Press.
 (2)    "Numerical Methods for Unconstrained Optimization"
        Edited by W.Murray. (1972) Academic Press.
 (3)    "Numerical Methods for Constrained Optimization"
        Edited by P.E.Gill and W.Murray. (1974) Academic Press.
 (4)    "Numerical Methods for Non-linear Algebraic Equations"
        Edited by P.Rabinowitz. (1970) Gordon and Breach.
 (5)    "Minimization Subject to Bounds on the Variables"
        P.E.Gill and W.Murray. NPL Report NAC 72 (1976).
 (6)    "Safeguarded Steplength Algorithms for Optimization Using
        Descent Methods" P.E.Gill and W.Murray. NPL Report NAC 37
        (1974).
 (7)    "The Implementation of Two Modified Newton Methods for
        Unconstrained Optimization"
        P.E.Gill, W.Murray, and S.M.Picken. NPL Rpt. NAC 24  (1972)
 (8)    "The Implementation of Two Quasi-Newton Methods for
        Unconstrained Optimization"
        P.E.Gill, W.Murray, and R.A.Pitfield. NPL Rpt NAC 11 (1972)
 (9)    "Conjugate Gradient Methods with Inexact Searches"
        D.F.Shanno. Math. of Operations Res. 3 (1978) 244-256

Author:
Jan Scheutjens (1947-1992), Wageningen Agricultural University, NL.

C Copyright (1980) (1981-1989) Wageningen Agricultural University, NL.

C++ translation:
Peter Barneveld, Wageningen Agricultural University, NL.

C Copyright (1993) Wageningen Agricultural University, NL.

 *NO PART OF THIS WORK MAY BE REPRODUCED, EITHER ELECTRONICALLY OF OTHERWISE*

*/

public:
	bool pseudohessian,samehessian,LM_BFGS,BFGS_damped,CG,Scheutjens,alpha_cg,No_LS,SD,TN,secant;
	bool d_info,e_info,i_info,g_info,h_info,s_info,v_info,x_info,picard,CG_F,BRR,DIIS,Pure_DIIS;
	bool ignore_newton_direction;
	bool haltonFPE;
	double max_accuracy_for_hessian_scaling,ys_div_yy,ss_div_ys,epsilon,normgk,normgk_1,alpha_TN;
	SFNewton();
	virtual ~SFNewton();

	virtual void residuals(double *const,double *const);
	virtual void residuals(Vector);
	virtual void inneriteration(double *const,double *const,double);
	virtual void iterate(double*,int);

	void vector_iterate(double *const, int);

	bool use_vector_iterations;
	int amount_of_stored_vectors,n_iterations_CG,reset_CG;
	double lowest_error;
	double reset_limit_gradient_inversion;
	bool print_improvement_info;
	bool print_exportable_info;
	int print_iteration_info;
	bool print_common_info;
	bool print_verbose_info;
	bool no_stars;

	void resethessian();
	void settolerance(double);
	void setlinetolerance(double);
	void setdeltamin(double);
	void setdeltamax(double);
	void setiterationlimit(int);
	void setDIIS(int);
	void setlinesearchlimit(int);
	void setoutput(std::ostream &);
	void setopenline(const char *const);
	void setcloseline(const char *const);
	int getiterations() const;
	int getiterationlimit() const;
	int getlinesearchlimit() const;
	double getaccuracy() const;
	double getdeltamax() const;
	double getdeltamin() const;
	double getalpha() const;
	double getminimum() const;
	double gettolerance() const;
	bool getnewtondirection() const;
	int getfunctioncalls() const;
	int getNumStoredArrays() const;

private:
	bool noresiduals, newtondirection;

	char bell;
	const char *openline,*closeline;

	int
	it,iterations,itmax,iterationlimit,lineiterations,linesearchlimit,it_CG,lm,it_TN,diis,
	nvar, functioncalls, trouble, resetiteration, nbits, n_reset,n_ignore, n_reset_hessians, m, n_LMBFGS,m_in_q,m_active;

	double
	minimum,tolerance,accuracy,linetolerance,PG,T_error,
	deltamin,deltamax,trustregion,trustfactor,mean_error,fluctuations,alphaMax,alphaMin;

	double *x; //ptr to vector supplied to this class
	double *x0, *p, *p0, *g, *g0, *xb, *r, *dg, *s, *d, *Mr, *sk, *yk, *CC, *DD, *EE, *tk, *I_DTC, *RR, *DTg, *Sigma, *VT;

	double *Nv, *Cv;
	double *Aij, *Ui, *Ci, *Ii, *Apij, *Rho, *Alphad, *X_X0, *G_G0, *XR;


	float *h; // hessian matrix
	//float *C;

	ostream *out;

	// output functions
	void outline(const char *);
	void warning(const char *);
	void outnumber(int,const char *,const char *);
	void iterationinfo();
	void matrixinfo();
	void statusinfo();
	void vectorinfo();
	void preconditioned_conjugated_gradient(double*,int);
	void BRR_(double* x, int nv);
	void BRR_step();
	void conjugated_gradient(double*,int);
	double Dot(double*,double*,int);
	void PlainLBFGSDIIS(double*,int);
	// calculation functions
	void decomposition();
	double newfunction();
	void newgradient(double *const);
	void direction(Vector);
	bool zero(double,bool);
	double linecriterion();
	double stepchange();
	void newtrustregion();
	double residue();
	void initialize_iteration();
	void terminate_iteration();
	double linesearch(double,bool);
	void newdirection(Vector,int&);
	void numhessian();
	void newhessian();
	void startderivatives();
	void findhessian();
	void StoreData(Vector, Vector, double* , double* , double* , double* );
	bool Hp(double*, double*, double, double);
	double Pg(double*, double);
	void PCG(Vector, bool);
//	void iterate_picard(void);
	Vector L_BFGS(Vector, double*, bool);
//	void itinfo(void);

	Vector rv;
	Vector pv;
	Vector zv;

	double* x_addr;
	long gradient_call_count;
	double relative_improvement;

	Vector x_prev;
	Vector r_prev;

	deque<Vector> x_diff;
	deque<Vector> r_diff;
	//deque<Vector> Gd;
	//deque<Vector> Xd;

	deque<double> rho;
	deque<double> abs_error;
	deque<double> ys_d_yy;

	deque<Vector> x_x0;
	deque<Vector> xR;

	double alpha;
	double beta;

	Vector computeZ(Vector&,double);
	Vector computeGradient(Vector);
	double computeBeta(Vector&);
	double computeAlpha(Vector&);

	enum ZMethod { NEGATIVE_GRADIENT, LIMITED_STORAGE_PSEUDO_HESSIAN };
	enum AlphaMethod { ALPHA_ALWAYS1, ALWAYS_INV_ITVARS, R_DEPENDENT, BISECTIONS, PARABOLA, RESET_DEPENDENT, RR_DIV_PP};
	enum BetaMethod { BETA_ALWAYS1, BETA_ALWAYS0, RVDxG_OVER_RVDxPV, R_OVER_R_PREV, CONVERGENCE_DEPENDENT, Polak_Ribiere};

	ZMethod z_method;
	AlphaMethod alpha_method;
	BetaMethod beta_method;

	clock_t tick_timer, algorithm_ticks, gradient_ticks;

};



#endif
