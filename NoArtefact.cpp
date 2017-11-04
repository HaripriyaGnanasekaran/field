#include "NoArtefact.h"
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

NoArtefact::NoArtefact(Lattice* Lat_,
					   SF_MoleculeList* MolQ_,
					   SF_SegmentList* SegQ_,
					   SF_Solve* Solve_,
					   double tolerance_) {
	Lat = Lat_;
	MolQ = MolQ_;
	SegQ = SegQ_;
	Solve = Solve_;
	tol = tolerance_;
	int numMol = MolQ->GetNumMolecules();
	for (int i=1; i<=numMol; i++) {
		SF_Molecule* Mol = MolQ->GetMolecule(i);
		if (Mol->GetFreedom() == secondGeneration && Mol->GetPhiBulk() > 0) {
			Message(fatal,"lattice artefact iteration not possible"
			" for second generation molecules with a non-zero phibulk");
		}
		if (Mol->GetFreedom() == thirdGeneration) {
			Message(fatal,"lattice artefact iteration not possible"
			" for third generation molecules");
		}
		if (Mol->GetFreedom() == rangeRestricted) {
			Message(fatal,"lattice artefact iteration not possible"
			" for range restricted molecules");
		}
	}
}
void
NoArtefact::Go(void) {
	double X = MolQ->GetFreeEnergy();
	if (Solve->e_info || Solve->s_info) {
		cout << "LATTICE ARTEFACT ITERATION layer adjustment: " 
			<< Lat->GetLayerAdjustment() << " delta: " << X << endl;
	}
	double layerAdjust = 1.1;
	double X2 = Delta(layerAdjust);
	if (X2 > X) {
		layerAdjust = 0.9;
		double X3 = Delta(layerAdjust);
		if (X3 > X) { // minimum found
			layerAdjust = 1;
		} else {
			while (X3 < X) {
				layerAdjust -= 0.1;
				if (layerAdjust < 0.7) {
					layerAdjust++;
					Solve->Jump(1);
				}
				X = X3;
				X3 = Delta(layerAdjust);
			}
			layerAdjust += 0.1;
			if (layerAdjust > 1.8) {
				layerAdjust--;
				Solve->Jump(-1);
			}
		}
	} else {
		layerAdjust = 1.2;
		double X3 = Delta(layerAdjust);
		if (X3 > X2) { // minimum found
			layerAdjust = 1.1;
		} else {
			while (X3 < X2) {
				layerAdjust += 0.1;
				if (layerAdjust > 1.8) {
					layerAdjust--;
					Solve->Jump(-1);
				}
				X2 = X3;
				X3 = Delta(layerAdjust);
			}
			layerAdjust -= 0.1;
			if (layerAdjust < 0.7) {
				layerAdjust++;
				Solve->Jump(+1);
			}
		}
	}
	double ax,bx,cx,af,bf,cf;
	bx = layerAdjust;
	bf = Delta(layerAdjust);
	ax = layerAdjust-0.1;
	af = Delta(layerAdjust-0.1);
	cx = layerAdjust+0.1;
	cf = Delta(layerAdjust+0.1);
	if (af < bf || cf < bf) {
		Message(literal,"LATTICE ARTEFACT ITERATION Unable to do "
			"artefact iteration, history in system. It could be that your"
			" system is too large to handle or is there more than one "
			"phase boundary?");
		Solve->SetErrorOccurred(true);
		return;
	} else {
		if (Solve->e_info || Solve->s_info) {
			cout << "LATTICE ARTEFACT ITERATION environment of "
				"minimum found, trying to zoom in" << endl;
		}
	}
	// Brents algorithm from numerical recipies
	const int ITMAX=100;
	const double CGOLD=0.3819660;
	const double ZEPS=1.0e-10;
	int iter;
	double a,b,d=0,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
	double e=0.0;
	a=(ax < cx ? ax : cx);
	b=(ax > cx ? ax : cx);
	x=w=v=bx;
	fw=fv=fx=bf;
	for (iter=1;iter<=ITMAX;iter++) {
		xm=0.5*(a+b);
		tol2=2.0*(tol1=tol*fabs(x)+ZEPS);
		if (fabs(x-xm) <= (tol2-0.5*(b-a))) {
			if (Solve->e_info || Solve->s_info) {
				cout << endl << "LATTICE ARTEFACT ITERATION finished"
					<< endl << endl;
			}
			return;
		}
		if (fabs(e) > tol1) {
			r=(x-w)*(fx-fv);
			q=(x-v)*(fx-fw);
			p=(x-v)*q-(x-w)*r;
			q=2.0*(q-r);
			if (q > 0.0) p = -p;
			q=fabs(q);
			etemp=e;
			e=d;
			if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
				d=CGOLD*(e=(x >= xm ? a-x : b-x));
			else {
				d=p/q;
				u=x+d;
				if (u-a < tol2 || b-u < tol2) 
					d=SIGN(tol1,xm-x);
			}
		} else {
			d=CGOLD*(e=(x >= xm ? a-x : b-x));
		}
		u=(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
		fu=Delta(u);
		if (fu <= fx) {
			if (u >= x) a=x; else b=x;
			SHFT(v,w,x,u)
			SHFT(fv,fw,fx,fu)
		} else {
			if (u < x) a=u; else b=u;
			if (fu <= fw || w == x) {
				v=w;
				w=u;
				fv=fw;
				fw=fu;
			} else if (fu <= fv || v == x || v == w) {
				v=u;
				fv=fu;
			}
		}
	}
	Message(warning,"Too many iterations in brent");
	return;
}
double
NoArtefact::Delta(double layerAdjustment) {
	Lat->SetLayerAdjustment(layerAdjustment);
	Solve->samehessian = true;
	Solve->Iterate();
	double X = MolQ->GetFreeEnergy();
	cout.precision(-(int)log10(tol));
	if (Solve->e_info || Solve->s_info) {
		cout << "LATTICE ARTEFACT ITERATION layer adjustment: " 
			<< layerAdjustment << " delta: " << X << endl;
	}
	return X;
}
#undef SHFT
#undef SIGN
