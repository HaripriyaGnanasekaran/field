#include <math.h>
#include <fenk.h>
#include <fenk/newton.h>
Outfile *Out;
Infile *In;
Text Filename;
Text Filename_in;
Text Filename_out;

int Mx,My,Mz,mirror;
int nx,ny,nz;
int jx,jy,jz;

int main() {

	Sysout().Outtext("Extends vtk file"); Sysout().Outimage();
	Sysout().Outimage();

	Sysout().Outtext("Give vtk filename (no extension) \\G60\\");
	Sysout().Breakoutimage(); Sysin().Inimage();
	Filename=Copy(Sysin().Image.Strip()==notext ? "G60" :Sysin().Image.Strip());
	Sysout().Outtext("Give box size in x-direction: ");
	Sysout().Breakoutimage(); Sysin().Inimage();
	Mx =Sysin().Inint();
	Sysout().Outtext("Give box size in y-direction: ");
	Sysout().Breakoutimage(); Sysin().Inimage();
	My =Sysin().Inint();
	Sysout().Outtext("Give box size in z-direction: ");
	Sysout().Breakoutimage(); Sysin().Inimage();
	Mz =Sysin().Inint();
	Sysout().Outtext("1: periodic; 2:mirror1; Please select \\1\\ : ");
	Sysout().Breakoutimage(); Sysin().Inimage();
	if (Sysin().Image.Strip()==notext)  mirror = 1; else mirror = Sysin().Inint();
	if (mirror>1) mirror=2;
	if (mirror<2) mirror=1;
	Sysout().Outtext("nx \\2\\ : ");
	Sysout().Breakoutimage(); Sysin().Inimage();
	if (Sysin().Image.Strip()==notext)  nx = 2; else nx = Sysin().Inint();
	Sysout().Outtext("ny \\2\\ : ");
	Sysout().Breakoutimage(); Sysin().Inimage();
	if (Sysin().Image.Strip()==notext)  ny = 2; else ny = Sysin().Inint();
	Sysout().Outtext("nz \\2\\ : ");
	Sysout().Breakoutimage(); Sysin().Inimage();
	if (Sysin().Image.Strip()==notext)  nz = 2; else nz = Sysin().Inint();
	//Vector phi;
	//phi.Dim(1,Mx*My*Mz);
	Vector phi(1,Mx*My*Mz);


	Filename_in = Filename+".vtk";
	In= new Infile(Filename_in);
	In -> Open(Blanks(1000));
	In->Inimage(); //# vtk DataFile Version 3.0
	In->Inimage(); //vtk output
	In->Inimage(); //ASCII
	In->Inimage(); //DATASET STRUCTURED_POINTS
	In->Inimage(); //DIMENSIONS M M M
	In->Inimage(); //SPACING 1 1 1
	In->Inimage(); //ORIGIN 0 0 0
	In->Inimage(); //POINT_DATA M^3
	In->Inimage(); //SCALARS SFBox_profile double
	In->Inimage(); //LOOKUP_TABLE default

	jx=My*Mz;
	jy=Mz;
	jz = 1;

	for (int z=0; z<Mz; z++)
	for (int y=0; y<My; y++)
	for (int x=0; x<Mx; x++) {
		In->Inimage();
		phi[jx*(Mx-x-1)+jy*(My-y-1)+jz*z+1]=In->Inreal(); //Spiegeling in x en y for sabine's project.
	}

	In->Close();

	Filename_out= Filename+"_t.vtk";
	Out= new Outfile(Filename_out);
	Out -> Open(Blanks(1000));
	Out->Outtext("# vtk DataFile Version 3.0"); Out-> Outimage();
	Out->Outtext("vtk output"); Out-> Outimage();
	Out->Outtext("ASCII"); Out-> Outimage();
	Out->Outtext("DATASET STRUCTURED_POINTS"); Out-> Outimage();
	Out->Outtext("DIMENSIONS " ); Out -> Outint(nx*Mx,4); Out -> Outint(ny*My,4); Out -> Outint(nz*Mz,4);  Out-> Outimage();
	Out->Outtext("SPACING 1 1 1"); Out-> Outimage();
	Out->Outtext("ORIGIN 0 0 0"); Out-> Outimage();
	Out->Outtext("POINT_DATA "); Out -> Outint(nx*ny*nz*Mx*My*Mz,7); Out-> Outimage();
	Out->Outtext("SCALARS SFBox_profile double"); Out-> Outimage();
	Out->Outtext("LOOKUP_TABLE default "); Out-> Outimage();

	int xx,yy,zz;
	int counterx, countery, counterz;
	for (int z=0; z < nz*Mz; z++)
	for (int y=0; y < ny*My; y++)
	for (int x=0; x < nx*Mx; x++) {
		xx=x; yy=y; zz=z; counterx=0; countery=0; counterz=0;
		while (xx>Mx-1) {counterx++; xx=xx-Mx;}
		while (yy>My-1) {countery++; yy=yy-My;}
		while (zz>Mz-1) {counterz++; zz=zz-Mz;}
		if (mirror==1) {} //{xx--; yy--; zz--;} //moet gechecked.
		else {
			if (counterx % 2 == 0) {} //xx--; //moet gechecked in algemeen
			else xx = Mx-xx-1;
			if (countery % 2 == 0) {} //yy--;
			else yy = My-yy-1;
			if (counterz % 2 == 0) {} //zz--;
			else zz = Mz-zz-1;
		}

		Out -> Outreal(phi[jz*zz+jy*yy+jx*xx+1],8,16); Out->Outimage();
	}
	Out -> Close();
};
