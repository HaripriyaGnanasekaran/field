#include "sfnewton.h"
#ifdef  HAS_LAPACK
#include "f2c.h"
#include "clapack.h"
#endif
#include <iostream> 


// PUBLIC FUNCTIONS
SFNewton::SFNewton() {
	nbits = numeric_limits<double>::digits;
	it = iterations = itmax = iterationlimit = lineiterations = it_CG = linesearchlimit = nvar = lm = 0;
	diis=1;
	functioncalls = 0;
	trouble = resetiteration = 0;
	reset_CG = 0;

	minimum = tolerance = accuracy = alpha = beta = T_error = linetolerance = PG =
	deltamin = trustregion = trustfactor = deltamax  = mean_error=alphaMax=alphaMin=0.0;

	noresiduals = pseudohessian = TN =
	samehessian = d_info = e_info = i_info = g_info = h_info = s_info = x_info =  alpha_cg  = No_LS = Scheutjens = SD = secant =
	newtondirection = ignore_newton_direction = LM_BFGS  = BFGS_damped = CG = picard  = CG_F = BRR = DIIS = no_stars = Pure_DIIS = false ;
	max_accuracy_for_hessian_scaling = 0.1;


	bell = '\a'; openline=""; closeline="";
	itmax = iterationlimit = 100;
	tolerance =  pow(10.0,-nbits/8); //alternative: DBL_EPSILON
	deltamax = pow(2.0,nbits); //alternative: 1/DBL_EPSILON or DBL_MAX
	linesearchlimit = 20;
	linetolerance = 9e-1;
	out = &cout;
	h = NULL;
	haltonFPE = false;
	amount_of_stored_vectors = m = m_in_q = m_active=0;
	epsilon = 0.1/pow(2.0,nbits/2);

}

SFNewton::~SFNewton() {
	if (h) {
		delete [] h;
	}
}


void SFNewton::inneriteration(double *const,double *const,double) {}

void SFNewton::residuals(double *const,double *const) {
	noresiduals = true;
}
void SFNewton::residuals(Vector) {
	noresiduals = true;
}

#ifdef  HAS_LAPACK
//#include "LapacTools.h"

void S_AxXpB(double* A, int Ma, int Na, double* X, double* Y, double* S)
{
	integer M = (integer)Ma;
	integer N = (integer)Na;
	integer LDA=(integer)Ma;
	integer INCX=1;
	integer INCY=1;

	char TRANS='N';
	double ALPHA = 1.0;
	double BETA = 1.0;

	memcpy(S, Y, sizeof(*Y)*Ma);
	dgemv_(&TRANS,&M,&N,&ALPHA,A,&LDA,X,&INCX,&BETA,S,&INCY);
}

void S_AxXmB(double* A, int Ma, int Na, double* X, double* Y, double* S)
//okay
//A(Ma,Na))
//X(Na)
//Y(Ma)
//S(Ma)
{
	integer M = (integer)Ma;
	integer N = (integer)Na;
	integer LDA=(integer)Ma;
	integer INCX=1;
	integer INCY=1;

	char TRANS='N';
	double ALPHA = 1.0;
	double BETA = -1.0;

	memcpy(S, Y, sizeof(*Y)*Ma);
	dgemv_(&TRANS,&M,&N,&ALPHA,A,&LDA,X,&INCX,&BETA,S,&INCY);
}

void Y_AxX(double* A, int Ma, int Na, double* X, double*Y)
//okay
//A(Ma,Na))
//X(Na)
//Y(Ma)
{
	integer M = (integer)Ma;
	integer N = (integer)Na;
	integer LDA=(integer)Ma;
	integer INCX=1;
	integer INCY=1;

	char TRANS='N';
	double ALPHA = 1.0;
	double BETA = 0.0;

	dgemv_(&TRANS,&M,&N,&ALPHA,A,&LDA,X,&INCX,&BETA,Y,&INCY);
}

void Y_ATxX(double* A, int Ma, int Na, double* X, double*Y)
//okay
//A(Ma,Na))
//X(Ma)
//Y(Na)
{
	integer M = (integer)Ma;
	integer N = (integer)Na;
	integer LDA=(integer)Ma;
	integer INCX=1;
	integer INCY=1;

	char TRANS='T';
	double ALPHA = 1.0;
	double BETA = 0.0;

	dgemv_(&TRANS,&M,&N,&ALPHA,A,&LDA,X,&INCX,&BETA,Y,&INCY);
}

void Y_mATxX(double* A, int Ma, int Na, double* X, double*Y)
//okay
{
	integer M = (integer)Ma;
	integer N = (integer)Na;
	integer LDA=(integer)Ma;
	integer INCX=1;
	integer INCY=1;

	char TRANS='T';
	double ALPHA = -1.0;
	double BETA = 0.0;

	dgemv_(&TRANS,&M,&N,&ALPHA,A,&LDA,X,&INCX,&BETA,Y,&INCY);
}


void C_AxBT(double* C, double* A, int Ma, int Na, double* B, int Mb, int Nb)
//okay
//A(Ma,Na)=A(M,K)
//B(Mb,Nb)=B(M,K) en dus BT(K,M)
//C(Ma,Mb)=C(M,M)
//dus Na=Nb
{
	integer M = (integer)Ma;
	integer K = (integer)Na;
	integer N = (integer)Mb;
	integer LDA=(integer)Ma;
	integer LDB=(integer)Mb;
	integer LDC=(integer)Ma;
	char TRANSA='N';
	char TRANSB='T';
	double ALPHA = 1.0;
	double BETA = 0.0;
	dgemm_(&TRANSA,&TRANSB,&M,&N,&K,&ALPHA,A,&LDA,B,&LDB,&BETA,C,&LDC);
}
void C_ATxB(double* C, double* A, int Ma, int Na, double* B, int Mb, int Nb)
//okay
//A(Ma,Na)=A(K,Na) en dus AT(N,K)
//B(Mb,Nb)=B(K,Nb)
//C(Na,Nb)
//dus Na=Nb
{
	integer M = (integer)Na;
	integer K = (integer)Ma;
	integer N = (integer)Nb;
	integer LDA=(integer)Ma;
	integer LDB=(integer)Mb;
	integer LDC=(integer)Na;
	char TRANSA='T';
	char TRANSB='N';
	double ALPHA = 1.0;
	double BETA = 0.0;
	dgemm_(&TRANSA,&TRANSB,&M,&N,&K,&ALPHA,A,&LDA,B,&LDB,&BETA,C,&LDC);
}

void C_mATxBpI(double* C, double* A, int Ma, int Na, double* B, int Mb, int Nb)
//okay
//A(Ma,Na)=A(K,Na) en dus AT(N,K)
//B(Mb,Nb)=B(K,Nb)
//C(Na,Nb)
//dus Na=Nb
{
	integer M = (integer)Na;
	integer K = (integer)Ma;
	integer N = (integer)Nb;
	integer LDA= (integer)Ma;
	integer LDB= (integer)Mb;
	integer LDC= (integer)Na;
	char TRANSA='T';
	char TRANSB='N';
	double ALPHA = -1.0;
	double BETA = 0.0;
	dgemm_(&TRANSA,&TRANSB,&M,&N,&K,&ALPHA,A,&LDA,B,&LDB,&BETA,C,&LDC);
	for (int i=0; i<Na; i++) C[i*Na+i]++;
}

void C_AxB(double* C, double* A, int Ma, int Na, double* B, int Mb, int Nb)
//A(Ma,Na)=A(M,K)
//B(Mb,Nb)=B(K,N)
//C(Ma,Nb)=C(M,N)
//dus Na=Mb
{
	integer M = (integer)Ma;
	integer K = (integer)Na;
	integer N = (integer)Nb;
	integer LDA= (integer)Ma;
	integer LDB= (integer)Mb;
	integer LDC= (integer)Ma;
	char TRANSA='N';
	char TRANSB='N';
	double ALPHA = 1.0;
	double BETA = 0.0;
	dgemm_(&TRANSA,&TRANSB,&M,&N,&K,&ALPHA,A,&LDA,B,&LDB,&BETA,C,&LDC);
}


void inverse(double* A, int N)
{
	integer *IPIV = new integer[N+1];
	integer LWORK;
	double *WORK;
	integer INFO;
	integer NN = (integer)N;
	int lwork;
	double WKOPT;
	dgetrf_(&NN,&NN,A,&NN,IPIV,&INFO);
	if (INFO >0) {//error message genereren
	};
	LWORK = -1;
	dgetri_(&NN,A,&NN,IPIV,&WKOPT,&LWORK,&INFO);
	lwork = (int)WKOPT;
	WORK = (double*)malloc( lwork*sizeof(double) );
	LWORK = (integer)lwork;
	dgetri_(&NN,A,&NN,IPIV,WORK,&LWORK,&INFO);
	if (INFO >0) {//error message genereren
	};
	delete IPIV;
	delete WORK;
}

void svd(double* A, int M, int N, double* U, double* S, double* VT)
//double A[LDA*N]=A[M,N]
//double S[N]
//double U[LDU*N]=U[M,N]
//double VT[LDVT*N] = VT[N,N]
//M>>N
//LDA =M
//LDU = M
//LDVT = N
{
	integer MM = (integer) M;
	integer NN = (integer) N;
	integer LDA=MM;
	integer LDU=MM;
	integer LDVT=NN;
	integer INFO;
	integer LWORK;
	int lwork;
	double WKOPT;
	double* WORK;
	char JOBU='S'; //'S' is nodig om alleen de eerste N colommen in U te schrijven.
	char JOBVT='A';

	LWORK = -1; //grootte hulpgeheugen aanvragen
	dgesvd_( &JOBU, &JOBVT, &MM, &NN, A, &LDA, S, U, &LDU, VT, &LDVT, &WKOPT, &LWORK, &INFO );
	lwork = (int)WKOPT;
	WORK = (double*)malloc( lwork*sizeof(double) );
	LWORK = (integer)lwork; //nu uitrekenen.
	dgesvd_( &JOBU, &JOBVT, &MM, &NN, A, &LDA, S, U, &LDU, VT, &LDVT, WORK, &LWORK, &INFO );
	if (INFO >0) { //error message genereren
	};
	delete WORK;
}

bool Ax_B(double* A, double* X, double* B, int N){
	integer NN = (integer) N;
	integer NRHS = 1;
	integer LDA=NN;
	integer LDAF=NN;
	integer LDB =NN;
	integer LDX = NN;
	double RCOND;
	char FACT = 'N';
	char UPLO = 'U';
	char EQUED = 'N';
	integer INFO;

	double* AF;
	AF = (double*)malloc( LDAF*NN*sizeof(double) );
	double* S;
	S = (double*)malloc( NN*sizeof(double) );
	double* FERR;
	FERR = (double*)malloc( sizeof(double) );
	double* BERR;
	BERR = (double*)malloc( sizeof(double) );
	double* WORK;
	WORK = (double*)malloc( 3*NN*sizeof(double) );
	integer* IWORK;
	IWORK = (integer*)malloc(NN*sizeof(integer) );

	dposvx_(&FACT,&UPLO,&NN,&NRHS,A,&LDA,AF,&LDAF,&EQUED,S,B,&LDB,X,&LDX,&RCOND,FERR,BERR,WORK,IWORK,&INFO);

	delete AF;
	delete S;
	delete FERR;
	delete BERR;
	delete WORK;
	return INFO==0;
}

void qr(double* A, int M, int N, double* R)
//double A[LDA*N]=A[M,N]
//double R[N,N]
//LDA =M
{
	integer MM = (integer) M;
	integer NN = (integer) N;
	integer LDA=MM;
	integer K=NN;
	integer INFO;
	integer LWORK;
	int lwork;
	double WKOPT;
	double* WORK;
	double* TAU = new double[N];

	LWORK = -1; //grootte hulpgeheugen aanvragen
	dgeqrf_( &MM, &NN, A, &LDA, TAU, &WKOPT, &LWORK, &INFO );
	lwork = (int)WKOPT;
	WORK = (double*)malloc( lwork*sizeof(double) );
	LWORK = (integer)lwork; //nu uitrekenen.
	dgeqrf_( &MM, &NN, A, &LDA, TAU, WORK, &LWORK, &INFO );
	if (INFO >0) { //error message genereren
	};
	for (int i=0; i<N; i++) {
		for (int j=0; j<N; j++) {
			if (j<=i) R[i*N+j]=A[i*M+j]; else R[i*N+j]=0;
		}
	}
	LWORK = -1;
	dorgqr_( &MM, &NN, &K, A, &LDA, TAU, &WKOPT, &LWORK, &INFO );
	lwork = (int)WKOPT;
	delete WORK;
	WORK = (double*)malloc( lwork*sizeof(double) );
	LWORK = (integer)lwork;
	dorgqr_( &MM, &NN, &K, A, &LDA, TAU, WORK, &LWORK, &INFO );
	delete WORK;
}

void SFNewton::BRR_step() {
	double normsk,xy;
	for (int i=0; i<nvar; i++) sk[i]=x[i]-x0[i]; normsk=norm2(sk,nvar);
	for (int i=0; i<nvar; i++) yk[i]=g[i]-g0[i];
	qr(DD, nvar, m, RR); //v
	C_AxBT(EE, CC, nvar, m, RR, m, m); //v
	svd(EE, nvar, m, CC, Sigma, VT); //vi
	int pos=0;
	for (int k=0; k<m; k++) { for (int i=0; i<nvar; i++) CC[pos++]*=Sigma[k];	} //vi
	C_AxBT(EE, DD, nvar, m, VT, m, m); //vi
	memcpy(DD, EE, sizeof(*EE)*nvar*m); //vi
	int j=iterations-1;
	if (j > m-1) {
		j=m-1;
		xy=0;
		for (int i=0; i<nvar; i++) xy=DD[j*nvar+i]*sk[i];
		for (int i=0; i<nvar; i++) CC[j*nvar+i] = (g[i]+ xy*CC[j*nvar+i])/normsk;
	} else {
		for (int i=0; i<nvar; i++) CC[j*nvar+i] = g[i]/normsk;
	}

	for (int i=0; i<nvar; i++) DD[j*nvar+i] =sk[i]/normsk;
	C_mATxBpI(I_DTC, DD, nvar, m, CC, nvar, m);
	inverse(I_DTC,m);
	Y_ATxX(DD, nvar, m, g, DTg);
	Y_AxX(I_DTC, m, m, DTg, tk); //i
	S_AxXmB(CC, nvar, m, tk, g, p);
}



void SFNewton::BRR_(double *x, int nv) {

	double normsk,normg,alphabound,alphamax,xy;
	if ( e_info ) {
		outline("Broyden Rank Reduction has been notified.");
		outline("Residual error:");
		*out << openline << "Your guess:" << closeline;
	}
	initialize_iteration();
	if (!alpha_cg && !No_LS ) Scheutjens=true;
	iterations=0;
	newgradient(g);
	for (int i=0; i<nv; i++) p[i]=-g[i];
	normg=norm2(g,nv);  minimum=pow(normg,2); 	accuracy = residue();

	while ((tolerance < accuracy || tolerance*10<normg) && iterations<iterationlimit && accuracy == fabs(accuracy) ) {
		iterationinfo(); iterations++; lineiterations = 0;
		C_mATxBpI(I_DTC, DD, nv, m, CC, nv, m); //Sherman Morrison formula for inversion of mxn matrix=  (B0+CDT)^-1= B0^_1-B0^-1 C(I+DTB0^-1C)^-1 D^T B0^-1 reduces to a nxn matrix inversion
		inverse(I_DTC,m);
		//Y_mATxX(DD, nv, m, g, DTg);	//Dit is dus fout (in artikel) maar 't stond dus goed in het proefschrift
		Y_ATxX(DD, nv, m, g, DTg);
		Y_AxX(I_DTC, m, m, DTg, tk); //i
		S_AxXmB(CC, nv, m, tk, g, p); //ii) //Het minteken hoort er dus. Was fout in het proefschrift en goed in artikel.
		memcpy(x0, x, sizeof(*x0)*nv);
		memcpy(g0, g, sizeof(*g0)*nv);
		if (Scheutjens) {
			newtrustregion();
			alphabound = alphamax = trustregion/(norm2(p,nvar)+1/pow(2.0,nbits));
			alpha = linesearch(alphabound,true); //changes x and g //iii
			trustfactor *= stepchange();
			trustfactor *= alpha/alphabound;
		}
		if (alpha_cg) {
			bool valid;
			alpha= alpha*norm2(p0,nvar)/norm2(p,nvar); memcpy(p0, p, sizeof(*g0)*nv);
			valid = Hp(s,p,0,alpha);
			if (valid) alpha = -deltamax*innerproduct(g,p,nvar)/innerproduct(s,p,nvar);
			for (int i=0; i<nvar; i++) x[i]=x0[i] + alpha*p[i];
			newgradient(g);
		}

		if (secant){
			double delta=1;
			double alpha_old=0;
			double alpha_current = 0;
			double alpha_new=0.00001;
			double pg_zero=innerproduct(g,p,nvar);
			if (i_info) {
				(*out).precision(5);
				cout << "   k = " << lineiterations << " pg = " << pg_zero << " alpha = " << alpha_new << endl;
			}
			double pg_old=0;
			double pg_current=pg_zero; if (pg_zero < 0) pg_zero = -1.0*pg_zero;
			double pg_tst; if (pg_current < 0) pg_tst = -pg_current; else pg_tst = pg_current;
			memcpy(x0, x, sizeof(*x0)*nvar);
			memcpy(g0, g, sizeof(*g0)*nvar);


			lineiterations=1;
			while (lineiterations <= linesearchlimit && pg_tst > linetolerance*pg_zero  && alpha_new >0) {
				alpha_old=alpha_current;
				alpha_current = alpha_new;
				pg_old = pg_current;
				pg_current=Pg(p,alpha_current);
				//if (pg_current < pg_old ) delta = delta/10; else delta = delta*2; if (delta>1) delta = 1;
				delta=1;
				alpha_new = alpha_current - delta*pg_current*(alpha_current-alpha_old)/(pg_current-pg_old);
				if (pg_current < 0) pg_tst = -pg_current; else pg_tst = pg_current;
				if (i_info) cout << "   k = " << lineiterations << " pg = " << pg_current << " alpha = " << alpha_new << endl;
				lineiterations++;
			}
			if (alpha_new < 0) {
				if (it==1) alpha = 1.0/norm2(g,nvar); else alpha = (ys_div_yy+ss_div_ys)/2.0;
				if (alpha < alphaMin) alpha=alphaMin;
				if (alpha > alphaMax) alpha=alphaMax;
				for (int i=0; i<nvar; i++) p[i]=-g[i];

			} else { alpha = alpha_new;	}
			alpha=alpha*deltamax;
			for (int i=0; i<nvar; i++) x[i]=x0[i] + alpha*p[i];
			newgradient(g);
			if (i_info) {
				cout << "   k = " << lineiterations << " pg = " << innerproduct(g,p,nvar) << " alpha = " << alpha << endl;
			}
			normg=sqrt(minimum);
			normgk=normg;
			if (lowest_error > normg) lowest_error = normg;
		}

		if (No_LS) {
			if (iterations==1) alpha = 1.0/normg;
			else {
				double ys=innerproduct(yk,sk,nvar);
				double yy=innerproduct(yk,yk,nvar);
				double ss=innerproduct(sk,sk,nvar);
				alpha = deltamax*(ys/yy+ss/ys)/2.0;///pow(log(normg),2);
			}
			for (int i=0; i<nvar; i++) x[i]=x0[i] + alpha*p[i];
			newgradient(g);
		}

		normg=norm2(g,nv); minimum=pow(normg,2); accuracy = residue();
		for (int i=0; i<nv; i++) sk[i]=x[i]-x0[i]; normsk=norm2(sk,nv);
		for (int i=0; i<nv; i++) yk[i]=g[i]-g0[i];
		qr(DD, nv, m, RR); //v
		C_AxBT(EE, CC, nv, m, RR, m, m); //v
		svd(EE, nv, m, CC, Sigma, VT); //vi
		int pos=0;
		for (int k=0; k<m; k++) { for (int i=0; i<nv; i++) CC[pos++]*=Sigma[k];	} //vi
		C_AxBT(EE, DD, nv, m, VT, m, m); //vi
		memcpy(DD, EE, sizeof(*EE)*nv*m); //vi
		int j=iterations-1;
		if (j > m-1) {
			j=m-1;
			xy=0;
			for (int i=0; i<nv; i++) xy=DD[j*nv+i]*sk[i];
			for (int i=0; i<nv; i++) CC[j*nv+i] = (g[i]+ xy*CC[j*nv+i])/normsk;
		} else {
			for (int i=0; i<nv; i++) CC[j*nv+i] = g[i]/normsk;
		}
		//for (int i=0; i<nv; i++) CC[j*nv+i] = yk[i] + sk[i];
		//for (int k=0; k<j; k++) {
		//	xy=0;
		//	for (int i=0; i<nv; i++) xy-=DD[k*nv+i]*sk[i];
		//	for (int i=0; i<nv; i++) CC[j*nv+i] += xy*CC[k*nv+i];
		//}
		//for (int i=0; i<nv; i++) CC[j*nv+i] = CC[j*nv+i]/normsk;

		//for (int i=0; i<nv; i++) CC[j*nv+i] = g[i]/normsk;
		for (int i=0; i<nv; i++) DD[j*nv+i] =sk[i]/normsk;
	}

	terminate_iteration();
}

void Ax(double* A, double* X, int N){//From Ax_B; below B is not used: it is assumed to contain a row of unities.
	double *U = new double[N*N];
	double *S = new double[N];
	double *VT = new double[N*N];
	integer MM = (integer)N, NN = (integer)N;
	integer LDA=MM, LDU=MM, LDVT=NN, INFO, LWORK;
	int lwork;
	double WKOPT;
	double* WORK;
	char JOBU='S'; //'S' is nodig om alleen de eerste N colommen in U te schrijven.
	char JOBVT='A';

	LWORK = -1; //grootte hulpgeheugen aanvragen
	dgesvd_( &JOBU, &JOBVT, &MM, &NN, A, &LDA, S, U, &LDU, VT, &LDVT, &WKOPT, &LWORK, &INFO );
	lwork = (int)WKOPT;
	WORK = (double*)malloc( lwork*sizeof(double) );
	LWORK = (integer)lwork; //nu uitrekenen.
	dgesvd_( &JOBU, &JOBVT, &MM, &NN, A, &LDA, S, U, &LDU, VT, &LDVT, WORK, &LWORK, &INFO );
	if (INFO >0) { //error message genereren
	};
	delete WORK;
	for (int i=0; i<N; i++) X[i]=0;
	for (int i=0; i<N; i++) for (int j=0; j<N; j++) X[i] += U[i*N + j];//*B[j];
	for (int i=0; i<N; i++) {S[i] = X[i]/S[i]; X[i]=0;} //S is use decause it is no longer needed.
	for (int i=0; i<N; i++) for (int j=0; j<N; j++) X[i] += VT[i*N + j]*S[j];
	delete U; delete S; delete VT;
}
#endif

double SFNewton::Dot(double *X, double *Y, int M){
	double result=0;
	for (int i=0; i<M; i++) result +=X[i]*Y[i];
	return result;
}

void SFNewton::PlainLBFGSDIIS(double *x, int iv){

	double eta=deltamax;
	//double ys_div_yy=1;
	double error;
	//double betad;
	double normC;
	//double Ddot;
	int k=0;
	int it=0;
	//int mj=0,up,
	int k_diis=1;
	int posi;
	initialize_iteration();

	newgradient(g);
	memcpy(x0, x, sizeof(*x0)*iv);
	for (int i=0; i<iv; i++) x[i]=-eta*g[i];
	for (int i=0; i<iv; i++) X_X0[i]=x[i]-x0[i];
	memcpy(XR, x, sizeof(*x0)*iv);

	error = norm2(g,iv);
	if ( e_info ) {
		outline("DIIS notified");
		outline("Residual error:");
		*out << openline << "Your guess:" << closeline;
	}
	cout << "i = " << it << " |g| = "<< error << endl;

	while (tolerance < error && it<iterationlimit) {
		it++;
		memcpy(x0, x, sizeof(*x0)*iv);
		//memcpy(g0, g, sizeof(*x0)*iv);
		newgradient(g);

		//if (it<0) {
		//	for (int i=0; i<iv; i++) G_G0[i+k*iv]=g[i]-g0[i];
		//	Ddot=Dot(X_X0+k*iv,G_G0+k*iv,iv);
		//	if (Ddot != 0) Rho[k] =1/Ddot;  else Rho[k]= Rho[(it-1)%m];
		//	//cout << "Ddot = " << Ddot << endl;
		//	ys_div_yy=1/Rho[k]/Dot(G_G0+k*iv,G_G0+k*iv,iv);
		//	up=k_diis-1;
		//	memcpy(p, g, sizeof(*x0)*iv);
		//	for (int mi = up; mi >= 0; mi--) {
		//		mj=k-up+mi; if (mj<0) mj +=k_diis;
		//		Alphad[mj] = Rho[mj]*Dot(X_X0+mj*iv,g,iv);
		//		for (int i=0; i<iv; i++) p[i] -=G_G0[i+mj*iv]*Alphad[mj];
		//	}
		//	for (int i=0; i<iv; i++) p[i] = p[i]*ys_div_yy;
//
		//	for (int mi = 0; mi <= up; mi++) {
		//		mj = k-up+mi; if (mj<0) mj +=k_diis;
		//		betad = Rho[mj]*Dot(G_G0+mj*iv,p,iv);
		//		for (int i=0; i<iv; i++) p[i] +=X_X0[i+mj*iv]*(Alphad[mj]-betad);
		//	}
		//	for (int i=0; i<iv; i++) x[i] -= p[i]*eta;
		//} else

		for (int i=0; i<iv; i++) x[i] -=g[i]*eta;

		k=it % m;  k_diis++;
		memcpy(XR+k*iv, x, sizeof(*x0)*iv);
		for (int i=0; i<iv; i++) X_X0[i+k*iv]=x[i]-x0[i];

		normC=0;
		if (k_diis>m) { k_diis =m;
			for (int i=1; i<m; i++) for (int j=1; j<m; j++)
			Aij[m*(i-1)+j-1]=Aij[m*i+j]; //remove oldest elements
		}
		for (int i=0; i<k_diis; i++) {
			posi = k-k_diis+1+i; if (posi<0) posi +=m;
			Aij[i+m*(k_diis-1)] = Aij[k_diis-1+m*i] = Dot(X_X0+posi*iv, X_X0+k*iv,iv);
		}
		for (int i=0; i<k_diis; i++) for (int j=0; j<k_diis; j++) {
			Apij[j+k_diis*i] = Aij[j+m*i];
		}

		Ax(Apij,Ci,k_diis);
		for (int i=0; i<k_diis; i++) normC +=Ci[i];
		for (int i=0; i<k_diis; i++) {Ci[i] =Ci[i]/normC; }

		posi = k-k_diis+1; if (posi<0) posi +=m;
		for (int i=0; i<iv; i++) x[i]=Ci[0]*XR[i+posi*iv];
		for (int i=1; i<k_diis; i++) {
			posi = k-k_diis+1+i; if (posi<0) posi +=m;
			for (int j=0; j<iv; j++) x[j]+=Ci[i]*XR[j+posi*iv];
		}
		error = norm2(g,iv);

		if ( e_info ) {
			if (print_iteration_info ==0) print_iteration_info=1;
			if (iterations % print_iteration_info == 0)
			cout << "i = " << it << " |g| = "<< error << endl;
		}

	}
	terminate_iteration();

}

void SFNewton::preconditioned_conjugated_gradient(double *x, int nv){
	double alphamax=0,alphabound=0;
	alphaMax=deltamax;
	alphaMin=deltamin;
	nvar = nv; 	initialize_iteration();
	Vector pv(1,nvar);
	Vector d0(1,nvar); for (int i=0; i<nvar; i++) d0[i+1]=0;
	int j=0;
	int j_max=10;
	double delta_new=0;
	double delta_mid;
	double delta_old;
	double sigma_0=deltamax;
	double eta_prev,eta;
	double delta_d=0;
	double inner_err;
	//double alpha_old;

	double* r = p0;
	double* d = p;
	double* s = xb;
	bool proceed; 	LM_BFGS=true;
	int k=0;
	if ( e_info ) {
		outline("Preconditioned nonlinear conjugate gradients with secant and Polak-Ribiere has been notified.");
		outline("Residual error:");
		*out << openline << "Your guess:" << closeline;
	}
	iterations=0;
	newgradient(g);
	newdirection(pv,k);
	newtrustregion();
	alphabound = alphamax = trustregion/(norm2(p,nvar)+1/pow(2.0,nbits));
	memcpy(x0, x, sizeof(*x0)*nvar);
	memcpy(g0, g, sizeof(*g0)*nvar);
	alpha = linesearch(alphabound,true); //changes x and g
	trustfactor *= stepchange();
	trustfactor *= alpha/alphabound;
	StoreData(pv,d0,x,x0,g,g0);

	for (int z=0; z<nv; z++) r[z] = -g[z];
	pv = L_BFGS(pv, r, true);
	//alpha_old=sigma_0;
	for (int z=0; z<nv; z++) {s[z]=pv[z+1]; d[z]=s[z]; delta_new += r[z]*d[z];}
	accuracy=pow(delta_new,0.5);
	cout << "i = " << iterations << " |g| = "<< accuracy << endl;
	while (tolerance < accuracy && iterations<iterationlimit) {
		j=0;
		delta_d=0;
		memcpy(x0, x, sizeof(*x0)*nvar);
		memcpy(g0, g, sizeof(*g0)*nvar);
		for (int z=0; z<nv; z++) {delta_d +=d[z]*d[z];}
		alpha = -sigma_0;
		//cout << "alpha " << alpha_old << endl;
		for (int z=0; z<nv; z++) x[z]+=sigma_0*d[z];
		newgradient(dg);
		eta_prev =0;
		memcpy(x, x0, sizeof(*x0)*nvar);
		for (int z=0; z<nv; z++) {eta_prev +=dg[z]*d[z];}
		inner_err=pow(delta_d,0.5);
		proceed=true;
		//alpha_old=0;
		while (proceed) {
			eta=0;
			for (int z=0; z<nv; z++) eta+=g[z]*d[z];
			alpha = alpha*eta/(eta_prev-eta);
			//alpha_old +=alpha;
			for (int z=0; z<nv; z++) x[z] += alpha*d[z];
			newgradient(g);
			eta_prev=eta;
			j++;
			proceed = (j<j_max && alpha*inner_err > 0.1*tolerance);
		}
		for (int z=0; z<nv; z++) r[z]=-g[z];
		delta_old = delta_new;
		delta_mid =0;
		for (int z=0; z<nv; z++) delta_mid+=r[z]*s[z];
		iterations++;
		StoreData(pv,d0,x,x0,g,g0);
		pv=L_BFGS(pv, r, true);
		delta_new=0;
		for (int z=0; z<nv; z++) {s[z]=pv[z+1]; delta_new += r[z]*s[z];}
		beta = (delta_new-delta_mid)/delta_old;
		k++;
		if (k == nv || beta < 0) {
			for (int z=0; z<nv; z++) d[z]=s[z]; k=0; beta=0;
		} else {
			for (int z=0; z<nv; z++) d[z]=s[z]+beta*d[z];
		}

		accuracy = pow(pow(delta_new,2),0.25);
		if ( e_info ) {
			if (print_iteration_info ==0) print_iteration_info=1;
			if (iterations % print_iteration_info == 0)
			cout << "i = " << iterations << " |g| = "<< accuracy << "  alpha("<<j<<") = " << alpha << "  beta = " << beta << endl;
		}
	}
	if (s_info) statusinfo();
	terminate_iteration();
}

void
SFNewton::conjugated_gradient(double *x, int nv) {
// CG with Newton-Raphson and Fletcher-Reeves
	int  k=0, j=0, j_max =5;
	double inner_err=0, delta_new=0, delta_old=0, delta_d=0;
	double teller=0, noemer=0;
	double rd=0;
	nvar = nv; 	initialize_iteration();
	double* r = xb;
	double* d = p;
	double* Hd =p0;
	bool proceed;
	//tolerance=1e-7;

	if ( e_info ) {
		outline("Nonlinear conjugate gradients with Newton-Raphson and Fletcher-Reeves has been notified.");
		outline("Residual error:");
		*out << openline << "Your guess:" << closeline;
	}
	iterations=0;
	newgradient(g);

	for (int z=0; z<nv; z++) {
		r[z] = -g[z];
		d[z]=r[z];
		delta_new += r[z]*r[z];
	}
	accuracy=pow(delta_new,0.5);
	cout << "i = " << iterations << " |g| = "<< accuracy << endl;
	while (tolerance < accuracy && iterations<iterationlimit) {
		j=0;
		delta_d=0;
		for (int z=0; z<nv; z++) delta_d += d[z]*d[z];
		inner_err = pow(delta_d,0.5);
		proceed=true;
		while (proceed) {
			teller=0;
			for (int z=0; z<nv; z++) {
				teller -= g[z]*d[z];
				x0[z]=x[z];
			}
			Hp(Hd,d,0,0);
			noemer=0;
			for (int z=0; z<nv; z++) noemer +=Hd[z]*d[z];
			alpha=teller/noemer;
			for (int z=0; z<nv; z++) {x[z]=x0[z]+alpha*d[z];}
			newgradient(g);
			j++;
			proceed =(j<j_max && alpha*inner_err > tolerance);
		}

		for (int z=0; z<nv; z++) r[z]=-g[z];
		delta_old=delta_new;
		delta_new=0;
		for (int z=0; z<nv; z++) delta_new += r[z]*r[z];
		beta=delta_new/delta_old;
		for (int z=0; z<nv; z++) d[z]=r[z]+beta*d[z];
		k++;
		rd=0;
		for (int z=0; z<nv; z++) rd +=r[z]*d[z];
		if (k == nv || rd<0) {
			k=0; beta=0;
			for (int z=0; z<nv; z++) d[z]=r[z];
		}
		iterations++;
		accuracy = pow(delta_new,0.5);
		if ( e_info ) {
			if (print_iteration_info ==0) print_iteration_info=1;
			if (iterations % print_iteration_info == 0)
			cout << "i = " << iterations << " |g| = "<< accuracy << "  alpha("<<j<<") = " << alpha << "  beta = " << beta << endl;
		}
	}
	if (s_info) statusinfo();
	terminate_iteration();
}


void SFNewton::iterate(double *x_, int nv) {
	double normg;
	if (nv < 1) return;
	nvar = nv;
	x=x_;
	if (CG_F) {
		if (amount_of_stored_vectors==1) {
			conjugated_gradient(x,nv);
			return;
		} else {
			preconditioned_conjugated_gradient(x,nv); return;
		}
	}
	bool no_lapack=false;
#ifdef HAS_LAPACK
	if (BRR) { BRR_(x,nv); return;}
	if (Pure_DIIS) {

			PlainLBFGSDIIS(x,nvar); return;}
#endif
	if (BRR) {
		no_lapack=true;
		cout << "BRR not available. Please compile with LAPACK=1. Proceed with LBFGS..." <<endl;
		LM_BFGS=true;
		BRR=false;
	}
	if (Pure_DIIS) {
		no_lapack=true;
		cout << "LBFGS + DIIS not available. Please compile with LAPACK=1. Proceed with LBFGS..." <<endl;
		LM_BFGS=true;
		Pure_DIIS=false;
	}

	Vector pv(1,nvar); for (int i=0; i<nvar; i++) pv[i+1]=0;
	Vector d0(1,nvar); for (int i=0; i<nvar; i++) d0[i+1]=0;

	if (!alpha_cg && !No_LS ) Scheutjens=true;


	int k=0;
	int k_diis=0;
	double alphamax=0,alphabound=0;
	alphaMax=deltamax;
	alphaMin=deltamin;
	initialize_iteration();
	T_error=0;it_CG=0;
	if (DIIS) {
			if (no_lapack) { cout << "DIIS not available. Please compile with LAPACK=1." <<endl;
			DIIS=false;
		} else {
			if (secant) { Scheutjens=true;alpha_cg=false;secant=false;No_LS=false;
				 cout << "DIIS  cannot be combined with secant linesearch; Scheutjens line search used instead" <<endl;
			}

			//if (linesearchlimit>1) {
			//	linesearchlimit=1;
			//	cout << "DIIS does not allow for a line search (for the time being): linesearchlimit is set to 1." <<endl;
			//}
		}
		x_x0.clear(); xR.clear(); k_diis=0;
	}

	if ( e_info ) {
		if (SD) 		{outline("Steepest decent has been notified.");}
		else if (CG) 		{outline("Conjugate gradients has been notified.");}
		else if (TN) 		{outline("Truncated Newton has been notified.");}
		else if (LM_BFGS) 	{outline("L-BFGS has been notified.");}
		else 			{outline("NEWTON has been notified.");}
		if (DIIS) {outline("DIIS activated");}
		outline("Residual error:");
		*out << openline << "Your guess:" << closeline;
	}

	newgradient(g);
	minimum = newfunction();
	if (TN || LM_BFGS) resethessian();
	newdirection(pv,k);
	normg=sqrt(minimum); normgk=normg; normgk_1 = normgk;
	lowest_error=normg;
	//double oldnormg = normg;

	while ((tolerance < accuracy || tolerance*10<normg) && iterations<iterationlimit && accuracy == fabs(accuracy) ) {
		//double maxerror=0;
		//double Evalue;
		//if (iterations%10==0) {
			//for (int i=0; i<nvar; i++) { Evalue = g[i]; if (maxerror<Evalue*Evalue) maxerror=Evalue*Evalue;
			//}
			//cout << "Max error: " << sqrt(maxerror) << endl;
		//}
		iterationinfo();
		if ( x_info || (d_info && alpha>0) || g_info ) {
			vectorinfo();
		}
		k = it = iterations = iterations+1;
		lineiterations = 0;
		if (Scheutjens) { //old methods.
			newtrustregion();
			// alternative for 1/pow(2.0,nbits) DBL_EPSILON
			// maybe DBL_MIN is better...
			alphabound = alphamax = trustregion/(norm2(p,nvar)+1/pow(2.0,nbits));

			memcpy(x0, x, sizeof(*x0)*nvar);
			memcpy(g0, g, sizeof(*g0)*nvar);

			alpha = linesearch(alphabound,!DIIS); //with DIIS we postpone gradient calculation.
			if (!DIIS) {
				trustfactor *= stepchange(); // alpha is modified as well!
				trustfactor *= alpha/alphabound;
			}
		};

		if (alpha_cg) {
			bool valid = true;
			memcpy(x0, x, sizeof(*x0)*nvar);
			memcpy(g0, g, sizeof(*g0)*nvar);
			normgk_1=normgk;


			if (alpha_TN >0) {alpha = alpha_TN; alpha_TN = 0; }
			else {
				if (alpha > 0) alpha= alpha*norm2(p0,nvar)/norm2(p,nvar);
				valid = Hp(s,p,0,alpha);
				if (valid) alpha = -deltamax*innerproduct(g,p,nvar)/innerproduct(s,p,nvar);
			}
			if (alpha <1e-6 || !valid) {
				alpha = 1.0/norm2(g,nvar); if (alpha > 1) alpha = 1;
				for (int i=0; i<nvar; i++) p[i] = - g[i];
			}
			while (!zero(alpha,!DIIS)) {
				memcpy(x, xb, sizeof(*xb)*nvar); newgradient(g); memcpy(x0, x, sizeof(*x0)*nvar);
				for (int i=0; i<nvar; i++) p[i]=-(1.0-0.5*rand()/RAND_MAX)*g[i]; resethessian();
				alpha *=0.01;
			}
			normg=sqrt(minimum);
			normgk=normg;
			//if (lowest_error > normg) lowest_error = normg;
		}
		if (secant){
			double delta=1;
			double alpha_old=0;
			double alpha_current = 0;
			double alpha_new=0.00001;
			double pg_zero=innerproduct(g,p,nvar);
			if (i_info) {
				(*out).precision(5);
				cout << "   k = " << lineiterations << " pg = " << pg_zero << " alpha = " << alpha_new << endl;
			}
			double pg_old=0;
			double pg_current=pg_zero; if (pg_zero < 0) pg_zero = -1.0*pg_zero;
			double pg_tst; if (pg_current < 0) pg_tst = -pg_current; else pg_tst = pg_current;
			memcpy(x0, x, sizeof(*x0)*nvar);
			memcpy(g0, g, sizeof(*g0)*nvar);


			lineiterations=1;
			while (lineiterations <= linesearchlimit && pg_tst > linetolerance*pg_zero  && alpha_new >0) {
				alpha_old=alpha_current;
				alpha_current = alpha_new;
				pg_old = pg_current;
				pg_current=Pg(p,alpha_current);
				//if (pg_current < pg_old ) delta = delta/10; else delta = delta*2; if (delta>1) delta = 1;
				delta=1;
				alpha_new = alpha_current - delta*pg_current*(alpha_current-alpha_old)/(pg_current-pg_old);
				if (pg_current < 0) pg_tst = -pg_current; else pg_tst = pg_current;
				if (i_info) cout << "   k = " << lineiterations << " pg = " << pg_current << " alpha = " << alpha_new << endl;
				lineiterations++;
			}
			if (alpha_new < 0) {
				if (it==1) alpha = 1.0/norm2(g,nvar); else alpha = (ys_div_yy+ss_div_ys)/2.0;
				if (alpha < alphaMin) alpha=alphaMin;
				if (alpha > alphaMax) alpha=alphaMax;
				for (int i=0; i<nvar; i++) p[i]=-g[i];

			} else { alpha = alpha_new;	}
			alpha=alpha*deltamax;
			for (int i=0; i<nvar; i++) x[i]=x0[i] + alpha*p[i];
			newgradient(g);
			if (i_info) {
				cout << "   k = " << lineiterations << " pg = " << innerproduct(g,p,nvar) << " alpha = " << alpha << endl;
			}
			normg=sqrt(minimum);
			normgk=normg;
			if (lowest_error > normg) lowest_error = normg;
		}
		if (No_LS) {

			if (it==1) alpha = 1.0/norm2(g,nvar); else alpha = deltamax*(ys_div_yy+ss_div_ys)/2.0;

			if (alpha < alphaMin) alpha=alphaMin;
			if (alpha > 10*alphaMax) alpha=10*alphaMax;


			memcpy(x0, x, sizeof(*x0)*nvar);
			memcpy(g0, g, sizeof(*g0)*nvar);

			while (!zero(alpha,!DIIS)) {
				memcpy(x, xb, sizeof(*xb)*nvar); newgradient(g); memcpy(x0, x, sizeof(*x0)*nvar);
				for (int i=0; i<nvar; i++) p[i]=-(1.0-0.5*rand()/RAND_MAX)*g[i]; resethessian();
				alpha *=0.01;
			}
		}

		if (DIIS) {
			k_diis++;

			if (k_diis>m) {
				k_diis=m;
				int jump =m+1;
				int pos;
				x_x0.pop_front(); xR.pop_front(); //remove storage from deque
				for (int i=1; i<m; i++)
				for (int j=1; j<m; j++) {pos=j+m*i; Aij[pos-jump]=Aij[pos];} //move elements of A up one and one to the left: dismiss oldest A elements.
			}

			for (int i=0; i< nvar; i++) pv[i+1] = -g[i];//+alpha*p[i];
				x_x0.push_back(pv+d0); //add to deque
			for (int i=0; i< nvar; i++) pv[i+1] = x0[i]-g[i]; //+alpha*p[i];
				xR.push_back(pv+d0); //add to deque

			//add elements to expanded A matrix
			for (int i=0; i<k_diis; i++) {
				Aij[i+m*(k_diis-1)] = Aij[k_diis-1+m*i] = x_x0[i].dotproduct(x_x0[k_diis-1]);
			}
			// write to (compressed) matrix Apij
			for (int i=0; i<k_diis; i++)
			for (int j=0; j<k_diis; j++) {
				Apij[j+k_diis*i] = Aij[j+m*i];
			}

			//find Ci
			double normC=0;
#ifdef HAS_LAPACK
			//if (Ax_B(Apij,Ci,Ui,k_diis)){ //Ui[i]=1 for all i; Ci should contain the solution.
				Ax(Apij,Ci,k_diis);
				if (iterations%diis==0) {
					for (int i=0; i<k_diis; i++) normC +=Ci[i];
					for (int i=0; i<k_diis; i++) {Ci[i] =Ci[i]/normC; }
					pv = Ci[0]*xR[0];
					for (int i=1; i<k_diis; i++) pv += Ci[i]*xR[i];
					for (int i=0; i<nvar; i++) x[i]=pv[i+1];
					//cout <<"!";
				}
			//} else {
			//	if (e_info) cout << "Unable to solve for Ci in DIIS. DIIS restarted." << endl;
			//	k_diis=0;
			//	x_x0.clear(); xR.clear();
			//}
#endif


			newgradient(g);
			bool valid=true;
			for (int i=0; i<nvar && valid; i++) {
				if (!finite(g[i])) {
					valid = false;
					warning("invalid numbers in gradient");
					g[i] = 1;
				}
			}
			minimum = newfunction();
			if (Scheutjens) {//Continue Scheutjens linesearch
				trustfactor *= stepchange(); // alpha is modified as well!
				trustfactor *= alpha/alphabound;
			}

		}

		if (LM_BFGS) StoreData(pv,d0,x,x0,g,g0);

		if (accuracy < lowest_error) {
			lowest_error=accuracy;
			memcpy(xb, x, sizeof(*xb)*nvar);
		}

		if ( !samehessian && !LM_BFGS) {
			newhessian();
		}

		newdirection(pv,k);
		normg=sqrt(minimum);
	}

	iterationinfo();
	if ( s_info&& !e_info) statusinfo();
	if ( iterations==iterationlimit ) {
		warning("[NEWTON:Maximum number of iterations reached *** problem terminated]");
		outline("Newton terminated.");
	} else if ( e_info ) {
		if (iterations==0) {
			outline("You hit the nail on the head.");
			outline("Newton done.");
		} else if (accuracy<.01*tolerance) {
			outline("That will do.");
			outline("Newton done.");
		} else {
			outline("Newton done.");
		}
	}
	terminate_iteration();
}

void SFNewton::StoreData(Vector pv, Vector d0, double* x, double* x0, double* g, double* g0){
	double r_dot_x;
	double noemer;
	if (iterations%nvar == 0) {resethessian();}
	if (m_in_q == m) {
		x_diff.pop_front();
		r_diff.pop_front();
		rho.pop_front();
		abs_error.pop_front();
		ys_d_yy.pop_front();
		m_in_q--;

	}
	for (int i=0; i< nvar; i++) {pv[i+1] = x[i]-x0[i]; } x_diff.push_back(pv+d0);

	//double	normpv =0;
	//for (int kk=0; kk<nvar; kk++) normpv +=pv[kk+1];
	//cout << "Norm pv x" << normpv << endl;

	for (int i=0; i< nvar; i++) {pv[i+1] = g[i]-g0[i];}  r_diff.push_back(pv+d0);

	//normpv =0;
	//for (int kk=0; kk<nvar; kk++) normpv +=pv[kk+1];
	//cout << "Norm pv r" << normpv << endl;

	m_in_q++;
	if (m_active<m) {m_active++;}
	r_dot_x =r_diff.back().dotproduct(x_diff.back());
	if (r_dot_x==0) {cout<<" can't compute rho" << endl; rho.push_back(1.0); r_dot_x=1;}
	else rho.push_back(1/r_dot_x);
	abs_error.push_back(norm2(g,nvar));
	noemer =r_diff.back().dotproduct(r_diff.back());
	if (noemer == 0) {cout << " can't compute ys_div_yy" << endl; ys_div_yy=1;
	}
	else ys_div_yy = r_dot_x/noemer;
	ys_d_yy.push_back(ys_div_yy);

	ss_div_ys = x_diff.back().dotproduct(x_diff.back())/r_dot_x;

}

void SFNewton::terminate_iteration(){
	delete [] x0;
	delete [] xb;
	delete [] p;
	delete [] g;
	delete [] g0;
	delete [] p0;

	if (alpha_cg || secant || TN ||  CG_F){
		delete [] dg;
	}
	if (alpha_cg || TN){
		delete [] s;
	}
	//if (LM_BFGS && !SD && !CG && !TN) {
		//delete [] Cv;
		//delete [] Nv;
	//}
	if (TN){
		delete [] r;
		delete [] d;
	}
	if (BRR){
		delete [] sk;
		delete [] yk;
		delete [] CC;
		delete [] DD;
		delete [] EE;
		delete [] tk;
		delete [] I_DTC;
		delete [] DTg;
		delete [] RR;
		delete [] VT;
		delete [] Sigma;
	}
	if (DIIS ) {
		delete [] Aij;
		delete [] Ci;
		delete [] Ui;
		delete [] Apij;
	}
	if (Pure_DIIS) {
		delete [] Aij;
		delete [] Ci;
		delete [] Ui;
		delete [] Apij;
		//delete [] G_G0;
		delete [] X_X0;
		delete [] XR;
		//delete [] Alphad;
		//delete [] Rho;
	}
}

void SFNewton::resethessian() {

	trouble = 0;
	if (LM_BFGS) {
		m_active = 0;
		n_reset++;
		n_reset_hessians++;
		//
		x_diff.clear(); r_diff.clear(); rho.clear(); abs_error.clear(); ys_d_yy.clear();
		m_in_q=0;
		//
	} else {
		startderivatives();
	}
	resetiteration = iterations;
}

void SFNewton::settolerance(double t) {
	if (t > 0) {
		tolerance = t;
	}
}
void SFNewton::setlinetolerance(double t) {
	if (t > 0) {
		linetolerance = t;
	}
}
void SFNewton::setdeltamin(double d) {
	if (d > 0 && d <= deltamax) {
		deltamin = d;
	}
}
void SFNewton::setdeltamax(double d) {
	if (d > deltamin) {
		deltamax = d;
	}
}
void SFNewton::setDIIS(int i) {
	if (i>0) {DIIS=true;
		diis = i;
	}
}

void SFNewton::setiterationlimit(int i) {
	if (i >= 0) {
		iterationlimit = i;
	}
}
void SFNewton::setlinesearchlimit(int i) {
	if (i >= 0) {
		linesearchlimit = i;
	}
}
void SFNewton::setoutput(ostream &o) {
	out = &o; // like this, out will never be NULL
}
void SFNewton::setopenline(const char *const o) {if(o){openline = o;}}
void SFNewton::setcloseline(const char *const c) {if(c){closeline = c;}}
bool SFNewton::getnewtondirection() const {return newtondirection;}
int SFNewton::getiterations() const {return iterations;}
int SFNewton::getiterationlimit() const {return iterationlimit;}
int SFNewton::getlinesearchlimit() const {return linesearchlimit;}
int SFNewton::getfunctioncalls() const {return functioncalls;}
double SFNewton::gettolerance() const {return tolerance;}
double SFNewton::getaccuracy() const {return accuracy;}
double SFNewton::getdeltamax() const {return deltamax;}
double SFNewton::getdeltamin() const {return deltamin;}
double SFNewton::getalpha() const {return alpha;}
double SFNewton::getminimum() const {return minimum;}
int SFNewton::getNumStoredArrays() const {return amount_of_stored_vectors;}

// PRIVATE FUNCTIONS
void SFNewton::initialize_iteration() {
	int i;
	it = iterations = functioncalls = 0;
	alpha = 1;
	ys_div_yy=ss_div_ys=1.0;
	n_reset = n_LMBFGS = n_reset_hessians = n_ignore = 0;
	trustregion = deltamax;
	trustfactor = 1;
	m=amount_of_stored_vectors;
	alpha_TN =0;
	it_TN=-m;

	if (reset_CG == 0) reset_CG = 10;
	x0 = new double[nvar];
	xb = new double[nvar];
	p = new double[nvar];
	g = new double[nvar];
	g0 = new double[nvar];
	p0 = new double[nvar];

	for (i=0; i<nvar; i++) {
		x0[i]=0;
		xb[i]=0;
		p[i]=0;
		g[i]=0;
		g0[i]=0;
		p0[i]=0;
	}

	if (BRR) {
		sk = new double[nvar];
		yk = new double[nvar];
		CC = new double[nvar*m];
		DD = new double[nvar*m];
		EE = new double[nvar*m];
		tk = new double[m];
		I_DTC = new double[m*m];
		RR = new double[m*m];
		DTg = new double[m];
		VT = new double[m*m];
		Sigma = new double[m];
		for (i=0; i<nvar; i++) {
			sk[i]=0;
			yk[i]=0;
		}
		for (i=0; i<m; i++) {
			tk[i]=0;
			DTg[i]=0;
			Sigma[i]=0;
		}
		for (i=nvar*m-1; i>-1; i--) {
			CC[i]=0;
			DD[i]=0;
			EE[i]=0;
		}
		for (i=m*m-1; i>-1; i--) {
			I_DTC[i]=0; RR[i]=0; VT[i]=0;
		}

	}
	if (alpha_cg || secant || TN || CG_F){
		dg = new double[nvar];
		for (i=0; i<nvar; i++) {
			dg[i]=0;
		}
	}
	if (alpha_cg || TN){
		s = new double[nvar];
		for (i=0; i<nvar; i++) s[i]=0;
	}

	if (TN){
		r = new double[nvar];
		d = new double[nvar];
		for (i=0; i<nvar; i++) {
			r[i]=0;
			d[i]=0;
		}
	}
	if (DIIS ) {
		Aij= new double[m*m]; for (i=0; i<m*m; i++) Aij[i]=0;
		Ci = new double[m]; for (i=0; i<m; i++) Ci[i]=0;
		Ui = new double[m]; for (i=0; i<m; i++) Ui[i]=1.0;
		Apij = new double[m*m]; for (i=0; i<m*m; i++) Apij[i]=0;

	}
	if (Pure_DIIS) {


		//Rho = new double[m];
		//Alphad = new double[m];
		XR = new double[m*nvar];
		X_X0 = new double[m*nvar];
		//G_G0 = new double[m*nvar];

		Aij= new double[m*m]; for (i=0; i<m*m; i++) Aij[i]=0;
		Ci = new double[m]; for (i=0; i<m; i++) Ci[i]=0;
		Ui = new double[m]; for (i=0; i<m; i++) Ui[i]=1.0;
		Apij = new double[m*m]; for (i=0; i<m*m; i++) Apij[i]=0;

	}
	(*out).precision(7);
	(*out).setf(ios::scientific,ios::floatfield);
}


void SFNewton::findhessian() {
	if ( !samehessian ) {
		if ( iterations==0 ) {
			resethessian();
		}
		numhessian(); // passes through residuals so check pseudohessian
		if (!pseudohessian) {
			if ( h_info ) {
				matrixinfo();
			}
			decomposition();
		}
	}
}

void SFNewton::newhessian() {
	if (LM_BFGS) {return;}
	int i=0;
	double dmin=0,sum=0,theta=0,php=0,dg=0,gg=0,g2=0,py=0,y2=0;
	dmin = 1/pow(2.0,nbits); // alternative: DBL_EPSILON or DBL_MIN
	if (!pseudohessian) { findhessian();
	} else if (!samehessian && alpha!=0 &&	iterations!=0) {
		double *y = new double[nvar];
		double *hp = new double[nvar];
		py = php = y2 = gg = g2 = 0;
		for (i=0; i<nvar; i++) {
			y[i] = dg = g[i]-g0[i];
			py += p[i]*dg;
			y2 += pow(y[i],2);
			gg += g[i]*g0[i];
			g2 += pow(g[i],2);
		}

		if ( !newtondirection ) {
			multiply(hp,1,h,p,nvar);
		} else {
			for (i=0; i<nvar; i++) {
				hp[i] = -g0[i];
			}

		}

		php = innerproduct(p,hp,nvar);
		theta = py/(10*dmin+alpha*php);
		if ( nvar>=1 && theta>0	&& iterations==resetiteration+1 && accuracy > max_accuracy_for_hessian_scaling) {
			if (e_info ) {
				(*(out)).precision(2);
				*(out) << openline
					<< "hessian scaling: " << theta
					<< closeline << endl;
			}
			//COMMENT See (9) page 248;
			alpha *= theta;
			py /= theta;
			php /= theta;
			for (i=0; i<nvar; i++) {
				p[i] /= theta;
				h[i+nvar*i] *= theta;
			}
		}
		if (nvar>=1) {
			trustfactor *= (4/(pow(theta-1,2)+1)+0.5);
		}
		if ( nvar>1 ) {
			sum = alpha*pow(norm2(p,nvar),2);
			//COMMENT See (4) page 132;
			theta = fabs(py/(alpha*php));
			if ( theta<.01 ) {
				sum /= 0.8;
			} else if ( theta>100 ) {
				sum *= theta/50;
			}
			for (i=0; i<nvar; i++) {
				y[i] -= alpha*hp[i];
			}
			updatpos(h,y,p,nvar,1.0/sum);
			trouble -= signdeterminant(h,nvar);
			if ( trouble<0 ) {
				trouble = 0;
			} else if ( trouble>=3 ) {
				resethessian();
			}
		} else if ( nvar>=1 && py>0 ) {
			trouble = 0;
			theta = py>0.2*alpha*php ? 1
			: 0.8*alpha*php/(alpha*php-py);
			if ( theta<1 ) {
				py = 0;
				for (i=0; i<nvar; i++) {
					y[i] = theta*y[i]+(1-theta)*alpha*hp[i];
					py += p[i]*y[i];
				}
			}
			updatpos(h,y,y,nvar,1.0/(alpha*py));
			updateneg(h,hp,nvar,-1.0/php);
		}

		delete [] y;
		delete [] hp;
	} else if ( !samehessian ) {
		resethessian();
	}

}

void SFNewton::startderivatives() {
	int i=0;
	float diagonal = 1+norm2(g,nvar);
	memset(h,0,nvar*nvar*sizeof(*h));
	for (i=0; i<nvar; i++) {
		h[i+nvar*i] = diagonal;
	}
}

void SFNewton::decomposition() {
	//COMMENT See (2) pages 64-67;
	int i=0,ntr=0;
	decompos(h,nvar,ntr);

	if (e_info) {
		if (iterations==0) {
			*out << openline;
			if (ntr==0) {
				*out << "Not too bad.";
			} else {
				*out << " sign might be wrong.";
			}
		} else if (ntr>0 && trouble==0) {
			if (ntr==1) {
				*out << "Wait a sec.";
			} else {
				*out << "some TROUBLES appear.";
			}
		} else if (ntr>trouble+1) {
			for (i=1; i<= ntr; i++) {
				*out << 'O';
			}
			*out << "H!";
		} else if (trouble>0) {
			if (ntr==0) {
				*out << "Here we go.";
			} else if (ntr<trouble-1) {
				*out << "Hold on.";
			} else if (ntr < trouble) {
				*out << "There is some hope for you.";
			} else if (ntr == trouble) {
				*out << "no Progress.";
			} else if (ntr > 4) {
				*out << "We won't fix it.";
			} else {
				*out << "More troubles.";
			}
		}
		if (iterations==0 || trouble>0 || ntr>0) {
			*out << closeline << endl;
		}
	}
	trouble = ntr;
}

double SFNewton::newfunction() {
	return pow(norm2(g,nvar),2);
}

void SFNewton::newgradient(double *const g1) {
#if USE_FENV
	int e;
#endif
	if ( !noresiduals ) {
		// TODO: add code for catching num. errors.
#if USE_FENV
		feclearexcept (FE_ALL_EXCEPT);
		if (haltonFPE) {
			feenableexcept(FE_OVERFLOW|FE_INVALID|FE_DIVBYZERO
				|FE_UNDERFLOW);
		}
#endif
		residuals(g1,x);

#if USE_FENV
		e = fetestexcept(FE_ALL_EXCEPT);
		e = fetestexcept(FE_ALL_EXCEPT);
		feclearexcept (FE_ALL_EXCEPT);
		fedisableexcept(FE_OVERFLOW|FE_INVALID|FE_DIVBYZERO|FE_UNDERFLOW);
		if (e && e != FE_INEXACT) {
			*out << "errors in residuals: ";
        		if (e & FE_DIVBYZERO) {
				*out << "FE_DIVBYZERO ";
			}
			if (e & FE_INVALID) {
				*out << "FE_INVALID ";
			}
			if (e & FE_OVERFLOW) {
				*out << "FE_OVERFLOW ";
			}
			if (e & FE_UNDERFLOW) {
				*out << "FE_UNDERFLOW ";
			}
			*out << endl;
		}
#endif
	}
	if ( !noresiduals ) {
		functioncalls++;
	} else {
		warning("[NEWTON:Iteration impossible: procedure residuals not found]");
	}
}

void SFNewton::numhessian() {
	int i=0,j=0;
	double dmax2=0,dmax3=0,di=0;
	double *g1, *xt;
	g1 = new double[nvar];
	xt = new double[nvar];
	dmax2 = pow(2.0,nbits/2); //alternative 2*pow(DBL_EPSILON,-0.5)?
	dmax3 = pow(2.0,nbits/3); //alternative 2*pow(DBL_EPSILON,-1.0/3)?
	for (i=0; i<nvar; i++) {
		xt[i] = x[i];
		di = (1/(dmax3*dmax3*fabs(h[i+nvar*i])+dmax3+fabs(g[i]))
			+1/dmax2)*(1+fabs(x[i]));
		if ( di<deltamin ) {
			di = deltamin;
		}
		x[i] += di;
		newgradient(g1);
		x[i] = xt[i];
		for (j=0; j<nvar; j++ ) {
			h[j+nvar*i] = (g1[j]-g[j])/di;
		}
	}
	delete [] g1;
	delete [] xt;
	newgradient(g);
}

double SFNewton::Pg(double *p,double alpha) {
	double pg=0;
	for (int i=0; i<nvar; i++) x[i] = x0[i] + alpha*p[i];
	newgradient(dg);
	for (int i=0; i<nvar; i++) pg +=dg[i]*p[i];
	return pg;
}

bool SFNewton::Hp(double *H_q, double *q, double normx,double delta) {
	bool valid;
	if (delta==0 ){
		double epsilon = 2e-8; //double precision. Machine error =10^-16; epsilon = 2 sqrt(precision)
		double normq = norm2(q,nvar);
		if (normx==0) normx = norm2(x0,nvar);
		//if ( (1+normx)/normq > 1000 || (1+normx)/normq < 0.001 ) {
			//delta = epsilon;
		//} else {
			delta = epsilon* (1+normx)/normq;
		//}
		//cout << "delta =  "<<delta <<" normx =" << normx << "normq = "<< normq << endl;
	}
	for (int i=0; i<nvar; i++) x[i] = x0[i] + delta*q[i];
	newgradient(dg);
	valid = true;
	for (int i=0; i<nvar && valid; i++) {
		if (!finite(dg[i])) {
			valid = false;
			warning("invalid numbers in gradient");
			dg[i] = 1;
		}
	}

	for (int i=0; i<nvar; i++) H_q[i] = (dg[i]-g[i])/delta;
	return valid;
}


void SFNewton::PCG(Vector pv, bool pre) {
	double alpha,beta,delta_old,delta_new,rd;
	//double q_old, q_new;
	double T_old=0;

	it_CG = 0;

	int it_limit;
	double err;
	//double rd;
	if (it>1) err=abs_error.back(); else err=10.1;
	it_limit = -5*log10(err) +3;
	if (it_limit < 3) it_limit=3;
	if (it_limit > linesearchlimit) it_limit = linesearchlimit;
	//cout << "it_limit: " << it_limit << endl;
	memcpy(x0, x, sizeof(*x0)*nvar);
	double normx = norm2(x0,nvar);
	Hp(s,p,normx,0);
	for (int i=0; i<nvar; i++) r[i] =  -g[i] - s[i];
	if (pre) {
		pv=L_BFGS(pv,r,true);
		for (int i=0; i<nvar; i++) d[i] = pv[i+1];
		delta_new = innerproduct(r,d,nvar);
	} else {
		for (int i=0; i<nvar; i++) d[i] = r[i];
		delta_new = innerproduct(r,r,nvar);
	}
	delta_old = 1.1*delta_new;
	if (delta_new>0) {T_error = sqrt(delta_new);} else  T_error = err*1.1;
	//rd = delta_new;
	T_old=T_error;
	//q_old=0;

	while ((it_CG < it_limit && T_error > 0.1*err && T_old > 0.5*T_error) || it_CG<1 ) { //&&  rd > -0.1*err*err && T_error > 0.01*err) {
		T_old=T_error;
		//q_old=q_new;
		Hp(s,d,normx,0);
		alpha = delta_new/innerproduct(s,d,nvar);
		for (int i=0; i<nvar; i++) {
			p0[i]=p[i];
			p[i]=p[i]+alpha*d[i];
		}

		for (int i=0; i<nvar; i++) { r[i]=r[i]-alpha*s[i]; }

		delta_old = delta_new;
		if (pre) {
			pv=L_BFGS(pv,r,true);
			for (int i=0; i<nvar; i++) s[i] = pv[i+1];
			delta_new = innerproduct(r,s,nvar);
		} else {
			delta_new = innerproduct(r,r,nvar);
		}
		beta = delta_new/delta_old;
		rd = innerproduct(r,d,nvar);
		if (rd<0) beta=0;
		for (int i=0; i<nvar; i++) d[i] = s[i] + beta*d[i];
		if (delta_new>0) T_error = sqrt(delta_new); else T_error = sqrt(-delta_new);
		//q_new=0.5*(rd+innerproduct(g,d,nvar));
		//cout << "1-q_o/q_n=" <<1-q_old/q_new << endl;
		it_CG++;
	}
	//if (it_CG == it_limit)     cout << "it limit"<< endl;
	//if (T_error < err)        cout << "convergence" << endl;
	//if (T_old < 0.5*T_error) cout << "walking backwards"<< endl;
	//memcpy(x, x0, sizeof(*x0)*nvar);
	//if (rd < 0) {
	//	memcpy(p, p0, sizeof(*p0)*nvar);
	//} else {
	//	alpha_TN = -deltamax*innerproduct(g,p,nvar)/innerproduct(s,p,nvar);
	//}
	if (T_error>10*err) {
		for (int i=0; i<nvar; i++) p[i] = -g[i];
		memcpy(x, x0, sizeof(*x0)*nvar);
	}

	alpha_TN = -deltamax*innerproduct(g,p,nvar)/innerproduct(s,p,nvar);
}

Vector SFNewton::L_BFGS(Vector pv, double *g, bool pv_g) {
	int upbound,lowbound;
	double vk,v1,v2;
	deque<double> alphad;
	deque<double> betad;

	upbound = m_in_q - 1;
	lowbound=m_in_q-m_active;

	if (pv_g ) {for (int i=0; i<nvar; i++) pv[i+1] = g[i]; }

	int mj=0;
	for (int mi = upbound; mi >= lowbound; --mi) {
		mj = mi; //+increment now zero
		alphad[mi] = rho[mj] * x_diff[mj].dotproduct(pv);
		pv = pv - alphad[mi] * r_diff[mj];
	}
	if (lowbound>upbound) vk=1; else
	{
		v1=ys_d_yy[lowbound];
		v2=ys_d_yy[upbound];
		if (v2>v1) vk=v2; else vk=v1;
	}

	//pv *= ys_div_yy;
	pv *=vk;


	for (int mi = lowbound; mi <= upbound; ++mi) {

		betad[mj] = rho[mj] * r_diff[mj].dotproduct(pv);
		pv = pv + x_diff[mj] * (alphad[mi] - betad[mi]);
	}

		//Note that -p is returned (als g gegeven is).

	return pv;
}

void SFNewton::direction(Vector pv) {

	if (CG) {
		//double theta;
		double yy=0;
		//double pg=0;
		double gg=0;
		double yg=0;
		double dg=0;
		double dy=0;
		double g0g=0;
		//double gs=0;
		//double ss=0;
		//double sy = 0;
		//double gdp = 0;
		//double dg0=0;// only needed for CG method of Liu and Storey or Fletcher.
		double g0g0=0;
		double upbound = m_in_q - 1;
		//double previous=upbound-1;
		//double beta_;
		//double y_y=0;
		//double d_y_=0;



		if (it>1) {
			for (int i=0; i<nvar; i++) {

					//gdp += g[i]*(p[i]-p0[i]);
				yg += g[i]  * r_diff[upbound][i+1];
				yy += r_diff[upbound][i+1] * r_diff[upbound][i+1];
					//y_y += r_diff[previous][i+1] * r_diff[upbound][i+1];
				dg += p[i] * g[i];

					//dg0 += p[i] * g0[i];// only needed for CG method of Liu and Storey or Fletcher.

				dy += p[i] * r_diff[upbound][i+1];
					//d_y_ += p0[i] * r_diff[previous][i+1];

				gg += g[i]*g[i]; //only needed for cg method of Dai and Yuan and Fletcher

				g0g += g0[i]*g[i]; //not needed for the method used now.
				g0g0 += g0[i]*g0[i]; //not needed for the method used now.
					//gs += g[i] * x_diff[upbound][i+1]; //not needed for the method used now.
					//sy +=  x_diff[upbound][i+1] *  r_diff[upbound][i+1];
					//ss +=  x_diff[upbound][i+1] *  x_diff[upbound][i+1];
			}
			beta = (yg - 2*dg*yy/dy)/dy; //Hanger and Zhang 2005.

				//theta = ss/sy;

				//beta = gdp/dg; //inspired by buckley AG 1978

				//beta = (theta*yg-gs)/sy;

				//  t = 2 * yy/sy;
				//  beta_DL = yg/dy; if (beta_DL<0) beta_DL=0;
				//  beta_DL -= t*gs/dy;
				//  mu = dg /dy;

				//if (sqrt(g0g*g0g) > 0.2*gg) beta=0;
				//if (beta<0) beta=gg/dy; //Dai and Yuan;
				//beta = -yg/dg0; //Liu and Storey

				//beta = -gg/dg0; //Fletcher


				//beta = gg/g0g0;
				//if (abs_error.back()>10*lowest_error) {beta=0; lowest_error = abs_error.back(); }

				//if (g0g < -0.1*normgk*normgk_1) {beta = 0; }

			pv=L_BFGS(pv,g,true);
			if (iterations%reset_CG ==0) {beta=0;
				for (int i=0; i<nvar; i++) p[i]=-pv[i+1];
			}
			else {
				//for (int i=0; i<nvar; i++) {
				//	  += p[i]*g[i];
				//}
				if (g0g+0.99*sqrt(gg*g0g0) < 0) {beta=0;}

				for (int i=0; i<nvar; i++) p[i]=-pv[i+1] + beta*p[i];
			}
			// for (int i=0; i<nvar; i++) p[i]=-g[i] + beta_DL*p[i]-mu*(r_diff.back()[i+1]-t*x_diff.back()[i+1]);
			//  beta = beta_DL;

		} else {
			beta=0;
			pv=L_BFGS(pv,g,true);
			for (int i=0; i<nvar; i++) p[i]=-pv[i+1];
		}
	} else
	if (SD) {
		for (int i=0; i<nvar; i++) p[i]=-g[i];
	} else
	if (TN) {
		pv=L_BFGS(pv,g,true);
		for (int i=0; i<nvar; i++) p[i]=-pv[i+1];
		//if (sqrt(minimum)<1) {
		if (((accuracy <100 && it % 10 == 0) || sqrt(minimum)<10) && (it>1)) {
			PCG(pv, true);
		}
	} else
	if (LM_BFGS) {
		//cout << "In BFGS" << endl;
		//for (int i=0; i<nvar; i++) pv[i+1] = -g[i]; // set starting value for pv (for LM BFGS)
		pv=L_BFGS(pv,g,true);
		for (int i=0; i<nvar; i++) p[i] = -pv[i+1];
	} else {
		//finally, standard newton....
		int i;

		if (iterations==0 && !samehessian) {
			if (h) delete [] h;
			h = new float[nvar*nvar];
			newhessian();
		}
		newtondirection = true;
		gausa(h,p,g,nvar);
		gausb(h,p,nvar);
		if (ignore_newton_direction) {
			newtondirection = true;
		} else {
			newtondirection = signdeterminant(h,nvar)>0;
		}
		if ( !newtondirection ) {
			for (i=0; i<nvar; i++) {
				p[i] *= -1;
			}
			if ( e_info && !no_stars) {
				*out << '*';
			}
		}
	}
}

bool SFNewton::zero(double alpha_,bool compute_g) {
	bool valid, timedep;
	int i=0;
	lineiterations++;
	if ( (lineiterations==5)  != TN) {
		memcpy(x, x0, sizeof(*x)*nvar);
		newgradient(g);
		valid = true;
		timedep = false;
		for (i=0; i<nvar && valid && !timedep; i++) {
			if ( g[i]!=g0[i] && !timedep) {
				*out << closeline << endl;
				warning("[NEWTON:ERROR?: your functions are time dependent!]");
				*out << endl;
				timedep = true;
			} else if (!finite(g[i]) && valid) {
				warning("invalid numbers in gradient, reversing search direction");
				valid = false;
				alpha_ *= -1; // reverse search direction
			}
		}
	}
	// adjust iteration variables
	for (i=0; i<nvar; i++) x[i] = x0[i]+alpha_*p[i];
	valid = true;

	if (compute_g) {
		newgradient(g);
		for (i=0; i<nvar && valid; i++) {
			if (!finite(g[i])) {
				valid = false;
				warning("invalid numbers in gradient");
				g[i] = 1;
			}
		}
		minimum = newfunction();
	}


	return valid;
}

double SFNewton::linesearch(double b, bool compute_g) {
	//COMMENT See (6);
	double newalpha = b<1 ? b : 1;
	zero(newalpha,compute_g);
	return newalpha;
}

double SFNewton::linecriterion() {
	double normg,gg0;
	normg = norm2(g0,nvar);
	gg0 = innerproduct(g,g0,nvar)/normg/normg;
	normg = pow(norm2(g,nvar)/normg,2);
	if ( (gg0>1 || normg>1) && normg-gg0*fabs(gg0)<0.2 ) {
		normg = 1.5*normg;
	}
	if (gg0<0 && normg<2) {
		return 1;
	} else if (normg>10) {
		return .01;
	} else {
		return 0.4*(1+0.75*normg)/(normg-gg0*fabs(gg0)+0.1);
	}
}

double SFNewton::stepchange() {
	double change, crit;
	change = crit = linecriterion();
	while ( crit<0.35 && lineiterations<linesearchlimit ) {
		alpha /= 4;
		zero(alpha,true);
		crit = linecriterion();
		change = 1;
	}
	return change;
}

void SFNewton::newtrustregion() {
	double normp0 = norm2(p0,nvar);

	if ( normp0>0 && trustregion>2*alpha*normp0 ) {
		trustregion = 2*alpha*normp0;
	}
	trustregion *= trustfactor;
	trustfactor = 1.0;
	if ( trustregion>deltamax ) {
		trustregion = deltamax;
	}
	if ( trustregion<deltamin ) {
		trustregion = deltamin;
	}
}

double SFNewton::residue() {
	return sqrt(norm2(p,nvar)*norm2(g,nvar)/(1+norm2(x,nvar)));
}

void SFNewton::newdirection(Vector pv, int& k) {

	inneriteration(g,x,accuracy);
	memcpy(p0, p, sizeof(*p0)*nvar);
	direction(pv);
	accuracy = residue();
	if ( it!=k ) {
		iterations = it;
	}
	if ( itmax!=100 ) {
		iterationlimit = itmax;
	}
}

// OUTPUT FUNCTIONS
void SFNewton::outline(const char *t) {
	*out << openline << t << closeline << endl;
}

void SFNewton::warning(const char *t) {
	//*out << bell;
	outline(t);
}

void SFNewton::outnumber(int n,const char *t1,const char *t2) {
	if ( n==0 ) {
		*out << "no";
	} else {
		*out << n;
	}
	*out << t1;
	if ( !(n==1 || t1==NULL || *t1 == '\0') ) {
		*out << 's';
	}
	*out << t2;
}

void SFNewton::statusinfo() {
	*out << openline;
	outnumber(iterations," iteration","");
	if ( functioncalls>0 ) {
		*out << "  ";
		outnumber(functioncalls, " residual","");
	}
	*out << closeline << endl;
}

void SFNewton::iterationinfo() {
	if ( !e_info ) return;
	if (print_iteration_info ==0) print_iteration_info=1;
	if (iterations % print_iteration_info == 0) {
		(*out).precision(5);
		*out << openline <<"i= "<< iterations <<"  " << accuracy << "  |g|=";
		if (minimum < 0) {
			(*out).precision(4);
		} else {
			(*out).precision(5);
		}
		if (minimum > 0) *out << sqrt(minimum); else *out << minimum;
		if ( iterations>0 && alpha!=2 )	{
			(*out).precision(2);
			*out << "  alpha=";
			*out << alpha;
		}
		if (BRR) {
			*out << "  S[m] = " << Sigma[m-1];
		}
		if (TN && T_error != 0) {
			*out << "  TN"<<it_CG << "="  << T_error;
			T_error = 0; it_CG = 0;
		}
		if (CG || (TN && beta != 0)) {
		    *out << "  beta=" << beta;
		}

		if ( iterations>0 && lineiterations>1 )	{
			*out << "  ";lineiterations--;
			outnumber(lineiterations," step","");
		}
		if (n_reset > 0)	*out << "  *"<<n_reset <<"*"; n_reset=0;
		if (n_ignore>0)      *out << "  !"<<n_ignore<<"!"; n_ignore=0;
		if (n_LMBFGS> 0)	*out << "  ("<<n_LMBFGS <<")"; n_LMBFGS=0;
		*out <<closeline << endl;
		if ((iterations>0 && (iterations%10)==0)){
			//|| iterations==iterationlimit || accuracy<tolerance) {
			if (s_info) statusinfo();
		}
	}
	//if ((iterations==iterationlimit || accuracy < tolerance) && !(iterations % print_iteration_info == 0)) {statusinfo();}
}

void SFNewton::vectorinfo() {
	int i=0;
	*out << closeline << endl << openline;
	if ( x_info ) {
		*out << "      variables";
	}
	if ( d_info && alpha>0 ) {
		*out << "      last step";
	}
	if ( g_info ) {
		*out << "      gradient";
	}
	*out << " (1:" << nvar << ')';
	for (i=0; i<nvar; i++) {
		if ( x_info || ((d_info && alpha>0) && p[i]!=0) ||
		(g_info && g[i]!=0)) {
			*out << closeline << endl << openline << '(' << i
				<< ')';
			if ( x_info ) {
				*out << setprecision(6) << x[i];
			}
			if ( (d_info && alpha>0 &&
				p[i]!=0) || (g_info && g[i]!=0)) {
				*out << setprecision(6) << alpha*p[i];
			}
			if ( g_info && g[i]!=0) {
				*out << setprecision(6) << g[i];
			}
		}
	}
	*out << closeline << endl << endl;
}

void SFNewton::matrixinfo() {
	int i=0,j=0,jj=0;
	*out << closeline << endl << openline << "hessian (1:" << nvar
		<<",1:" << nvar << ')';
	for (j=0; j<nvar; j++) {
		for (i=0;i<nvar; i++) {
			if (h[i+nvar*j]!=0 ) {
				if (jj!=j ) {
					jj = j;
					*out << endl;
				}
				*out << closeline << endl << openline;
				if (i==j) {
					*out << "   *(";
				} else {
					*out << "(";
				}
				*out << i << ',' << j << ')'
					<< setprecision(12) << h[i+nvar*j];
			}
		}
	}
	*out << closeline << endl << endl;
}

