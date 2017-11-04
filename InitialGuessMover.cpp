#include <math.h>
#include <fenk.h>
#include <fenk/newton.h>
Outfile *Out;
Infile *In;
Text Filename;
Text Filename_in;
Text Filename_out;
int Increment;
Text molecule;
Text segname;



int M;
int counter;
int jx,jy,jz;
int Mx,My,Mz;
int xx,yy,zz;

int main() {

	Sysout().Outtext("Modifies the 3-G initial guess file by one layer in x y and z"); Sysout().Outimage();
	Sysout().Outimage();

	Sysout().Outtext("Give filename that contains initial guess \\Gyroid43.outi\\");
	Sysout().Breakoutimage(); Sysin().Inimage();
	Filename=Copy(Sysin().Image.Strip()==notext ? "Gyroid43.outi" : Sysin().Image.Strip());
	Sysout().Outtext("Increment: ");
	Sysout().Breakoutimage(); Sysin().Inimage();
	Increment = Sysin().Inint();
	if (Increment>0) Increment = 1;
	if (Increment<0) Increment = -1;
	if (Increment==0) Increment = 0;
	Sysout().Outtext("Number of iv per coordinate: ");
	Sysout().Breakoutimage(); Sysin().Inimage();
	counter = Sysin().Inint();


	Filename_in = Filename;
	Filename_out = Filename+"1";
	In = new Infile(Filename_in);
	In -> Open(Blanks(1000));
	Out= new Outfile(Filename_out);
	Out -> Open(Blanks(1000));

	In->Inimage(); //gradients
	Out->Outtext("gradients"); Out-> Outimage();
	In->Inimage(); //3
	Out->Outtext("3"); Out-> Outimage();
	In->Inimage(); Mx =In->Inint();
	Out->Outint(Mx+Increment,5); Out->Outimage();
	In->Inimage(); My =In->Inint();
	Out->Outint(My+Increment,5); Out->Outimage();
	In->Inimage(); Mz =In->Inint();
	Out->Outint(Mz+Increment,5); Out->Outimage();


	jx=My*Mz;
	jy=Mz;
	jz = 1;

	for (int k=1; k<=counter; k++) {
		In->Inimage();
		Out->Outtext("molecule"); Out-> Outimage();
		In->Inimage(); //all
		Out->Outtext("all"); Out-> Outimage();
		In->Inimage(); //state
		Out->Outtext("state"); Out-> Outimage();
		In->Inimage(); segname=Copy(In->Image.Strip());
		Out->Outtext(segname); Out->Outimage();
		Vector iv;
	    iv.Dim(1,Mx*My*Mz);
	    for (int i=1; i<=Mx*My*Mz; i++) {
			In->Inimage();
			iv[i]=In->Inreal();
		}
	    for (int x=1; x<=Mx+Increment; x++)
	    for (int y=1; y<=My+Increment; y++)
	    for (int z=1; z<=Mz+Increment; z++){
			if (x>Mx) xx=Mx-1; else xx = x-1;
			if (y>My) yy=My-1;else yy = y-1;
			if (z>Mz) zz=Mz-1;else zz = z-1;
			Out->Outreal(iv[jx*xx+jy*yy+jz*zz+1],20,26); Out->Outimage();
		}
	}

	In ->Inimage(); Out->Outtext(Copy(In->Image.Strip())); Out->Outimage();
	In ->Inimage(); Out->Outtext(Copy(In->Image.Strip())); Out->Outimage();
	for (int k=1; k<=counter; k++) {
		In ->Inimage(); Out->Outtext(Copy(In->Image.Strip())); Out->Outimage();
		In ->Inimage(); Out->Outtext(Copy(In->Image.Strip())); Out->Outimage();
		In ->Inimage(); Out->Outtext(Copy(In->Image.Strip())); Out->Outimage();
	}

	In->Close();
	Out -> Close();
};
