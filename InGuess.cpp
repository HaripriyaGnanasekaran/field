#include <math.h>
#include <cmath>
#include <iostream>
#include <fenk.h>
#include <fenk/newton.h>
Outfile *Out;
Infile *In;
Text Filename;
int M,NF,N,Z;
double k,alpha,Ginf,G1,G2,G3;
bool B1,B2,B3;
int mode;
Vector Nv;
Matrix OutIv;
Vector NewIv;

class SFsystem: public Newton {
	public:
	Vector iv;


	SFsystem() {
		iv.Dim(1,1);

		//jacobian=true;
		iterationlimit=1000;
		pseudohessian=false;
		tolerance=1e-8;
		deltamax=0.0001;
		e_info = false;

	}

	void gradient(Vector g, Vector iv, int dim) {
		Ginf = iv[1];
		double SumX=0,SumY=0,SumX2=0,SumXY=0;
		double y1;
		if (OutIv[1][Z]>Ginf) y1= log(OutIv[1][Z]-Ginf); else  y1= log(Ginf-OutIv[1][Z]);
		double y2;
		if (OutIv[2][Z]>Ginf) y2= log(OutIv[2][Z]-Ginf); else  y2= log(Ginf-OutIv[2][Z]);
		double y3;
		if (OutIv[3][Z]>Ginf) y3= log(OutIv[3][Z]-Ginf); else  y3= log(Ginf-OutIv[3][Z]);
		double x1= log(Nv[1]);
		double x2= log(Nv[2]);
		double x3= log(Nv[3]);
		//Sysout().Outtext("Ginf= "); Sysout().Outreal(Ginf,8,16); Sysout().Outimage();
		//Sysout().Outtext("Y1= "); Sysout().Outreal(y1,8,16); Sysout().Outimage();
		//Sysout().Outtext("Y2= "); Sysout().Outreal(y2,8,16); Sysout().Outimage();
		//Sysout().Outtext("Y3= "); Sysout().Outreal(y3,8,16); Sysout().Outimage();
		//Sysout().Outtext("X1= "); Sysout().Outreal(x1,8,16); Sysout().Outimage();
		//Sysout().Outtext("X2= "); Sysout().Outreal(x2,8,16); Sysout().Outimage();
		//Sysout().Outtext("X3= "); Sysout().Outreal(x3,8,16); Sysout().Outimage();
		//Sysout().Outtext("iv1= "); Sysout().Outreal(OutIv[1][Z],8,16); Sysout().Outimage();
		//Sysout().Outtext("iv2= "); Sysout().Outreal(OutIv[2][Z],8,16); Sysout().Outimage();
		//Sysout().Outtext("iv3= "); Sysout().Outreal(OutIv[3][Z],8,16); Sysout().Outimage();
		SumX = x1+x2+x3;
		SumY = y1+y2+y3;
		SumXY = x1*y1+x2*y2+x3*y3;
		SumX2 = x1*x1+x2*x2+x3*x3;
		alpha = (SumX*SumY-3*SumXY)/(SumX*SumX-3*SumX2);
		k = (SumY-alpha*SumX)/3;
		g[1] =pow(1-(alpha*x1+k)/y1,2);
		g[1]+=pow(1-(alpha*x2+k)/y2,2);
		g[1]+=pow(1-(alpha*x3+k)/y3,2);
		g[1] = sqrt(g[1]);
	}

	void iters() {
		for (int z=1; z<=M; z++) {
			Sysout().Outtext("z="); Sysout().Outint(z,4); Sysout().Outimage();
			Z=z;
			G1 = OutIv[1][Z];
			G2 = OutIv[2][Z];
			G3 = OutIv[3][Z];
			//iv[1]=G3*(G3/G2);
			iv[1]=G3+(G3-G2)/(log(Nv[3])-log(Nv[2]))*(log(N)-log(Nv[3]));
			if (G1 > G2 && G2 > G3) {
				mode = 1;
				iterate(iv,1,1);
				if (iterations < iterationlimit && G3>Ginf  && Ginf < G3*1.1) NewIv[Z] = Ginf+exp(k)*pow(N,alpha); else NewIv[Z]=G3;
			} else if (G1 < G2 && G2 < G3) {
				mode = 2;
				iterate(iv,1,1);
				if (iterations < iterationlimit && G3<Ginf && Ginf>G1*0.9) NewIv[Z] = Ginf+exp(k)*pow(N,alpha); else NewIv[Z]=G3;
			} else {
				double X1=log(Nv[1]); double X2=log(Nv[2]); double X3 = log(Nv[3]);
				double G12=G1-G2;
				double G13=G1-G3;
				double XX12=X1*X1-X2*X2;
				double XX13=X1*X1-X3*X3;
				double X12=X1-X2;
				double X13=X1-X3;
				double a= (G12*X13-G13*X12)/(XX12*X13-XX13*X12);
				double b= (G12-a*XX12)/X12;
				double c=G1-a*X1*X1-b*X1;
				mode = 3;
				NewIv[Z]= a*log(N)*log(N)+b*log(N)+c;
			}

			if (NewIv[Z]==NewIv[Z]) {NewIv[Z]=NewIv[Z];} else NewIv[Z] = OutIv[3][Z];
			if (Z>60) NewIv[Z]=OutIv[3][Z];

			Sysout().Outreal(OutIv[1][Z],8,16);
			Sysout().Outreal(OutIv[2][Z],8,16);
			Sysout().Outreal(OutIv[3][Z],8,16);
			Sysout().Outreal(NewIv[Z],8,16);
			Sysout().Outint(mode,5);
			Sysout().Outimage();
		}
	}
};

int main() {
	SFsystem* Iv;
	Sysout().Outtext("Extrapolate initial guesses.");                Sysout().Outimage();
	Sysout().Outimage();
	NF = 3; M=1500;
	//Sysout().Outtext("Give M:  ");                                Sysout().Breakoutimage(); Sysin().Inimage(); M = Sysin().Inint();


	Nv.Dim(1,NF);
	OutIv.Dim(1,NF,1,M);
	NewIv.Dim(1,M);

	for (int i=1; i<=NF; i++) {
		//Sysout().Outtext("Give input filename: ");
		//Sysout().Breakoutimage(); Sysin().Inimage();
		//Filename=Copy(Sysin().Image.Strip()==notext ? "Nofile" :Sysin().Image.Strip());
		//Sysout().Outtext("Give corresponding N:  ");              Sysout().Breakoutimage(); Sysin().Inimage(); Nv[i] = Sysin().Inint();
if (i==1) {In = new Infile("NOC6_100000.outi"); Nv[1]=20000;}
if (i==2) {In = new Infile("NOC6_200000.outi"); Nv[2]=50000;}
if (i==3) {In = new Infile("NOC6_300000.outi"); Nv[3]=100000;}
		//In= new Infile(Filename);
	    In->Open(Blanks(1000));
	    In->Inimage();
	    In->Inimage();
	    In->Inimage();
	    In->Inimage();
	    In->Inimage();
	    In->Inimage();
		In->Inimage();
		In->Inimage();
	    for (int z=1; z<=M; z++) {
			In->Inimage();
			OutIv[i][z]=In->Inreal();
		}

	    In->Close();
	}
	//Sysout().Outtext("Give extrapolation N:  ");              Sysout().Breakoutimage(); Sysin().Inimage(); N = Sysin().Inint();
N=500000;
	//Sysout().Outtext("Give output filename:  ");
	//Sysout().Breakoutimage(); Sysin().Inimage();
	//Filename=Copy(Sysin().Image.Strip()==notext ? "Nofile" :Sysin().Image.Strip());

	Iv = new SFsystem();
	Iv -> iters();

	//Out= new Outfile(Filename);
	Out =new Outfile("NOC6__500000.outi");
	Out->Open(Blanks(1000));
	Out->Outtext("gradients"); Out->Outimage();
	Out->Outtext("1");         Out->Outimage();
	Out->Outtext("1502");      Out->Outimage();
	Out->Outtext("molecule");  Out->Outimage();
	Out->Outtext("all");       Out->Outimage();
	Out->Outtext("state");     Out->Outimage();
	Out->Outtext("W");         Out->Outimage();
	Out->Outtext("0.0");       Out->Outimage();
	for (int z=1; z<=M; z++) {
		Out->Outreal(NewIv[z],12,20); Out->Outimage();
	}
	Out->Outreal(NewIv[1500],12,20); Out->Outimage();
	Out->Outtext("molecule");  Out->Outimage();
	Out->Outtext("all");       Out->Outimage();
	Out->Outtext("state");     Out->Outimage();
	Out->Outtext("A");         Out->Outimage();
	Out->Outtext("0.0");       Out->Outimage();
	Out->Outreal(exp(log(NewIv[1])+1.0),12,20); Out->Outimage();
	for (int z=2; z<=M; z++) {
		Out->Outreal(NewIv[z],12,20); Out->Outimage();
	}
	Out->Outreal(NewIv[1500],12,20); Out->Outimage();
	Out->Outtext("phibulk solvent");   Out->Outimage();
	Out->Outreal(1.0-1.0/N,12,20);      Out->Outimage();
	Out->Outtext("alphabulk");         Out->Outimage();
	Out->Outtext("W");                 Out->Outimage();
	Out->Outtext("1.0");               Out->Outimage();
	Out->Outtext("alphabulk");         Out->Outimage();
	Out->Outtext("A");                 Out->Outimage();
	Out->Outtext("1.0");               Out->Outimage();
	Out->Close();
};
