#include <math.h>
#include <fenk.h>
Outfile *Out;
Text Filename;

double A1,A2,B1,B2;
int n_layers_x=0;
int n_layers_y=0;
int n_layers_z=0;
double R;
double D;
int M;
int Np;
int Npt;
int H;
Vector coordinate;
Vector pos_x,pos_y,pos_z;

int main(int argc, char** argv) {
	
	Sysout().Outtext("Give filename \\holes.mask\\");
	Sysout().Breakoutimage(); Sysin().Inimage();
	Filename=Copy(Sysin().Image.Strip()==notext ? "holes.mask" : Lowcase(Sysin().Image.Strip()));

	Sysout().Outtext("Diameter of Particles ");  Sysout().Breakoutimage(); Sysin().Inimage();
	R=Sysin().Inreal()/2;

	Sysout().Outtext("Center-center distance in a row (even number please)");  Sysout().Breakoutimage(); Sysin().Inimage();
	D=Sysin().Inreal();

	Sysout().Outtext("Number of Particles in each row");  Sysout().Breakoutimage(); Sysin().Inimage();
	Np=Sysin().Inint()+1;
	Npt=(2*Np)*(Np);

	Sysout().Outtext("Cut-off height ");  Sysout().Breakoutimage(); Sysin().Inimage();
	H=Sysin().Inint();

	int i=0;
	double ypos, xpos;
	pos_x.Dim(1,Npt);
	pos_y.Dim(1,Npt);
	pos_z.Dim(1,Npt);

	xpos=1/2; ypos=1/2;
	for (int x=1; x<=Np; x++) 
	for (int y=1; y<=Np; y++) {i++;
		pos_x[i] = xpos+(x-1)*D; 
		pos_y[i] = ypos+(y-1)*D;
		pos_z[i] = R;
	}  

	xpos +=D/2; ypos +=D/2;
	for (int x=1; x<=Np; x++) 
	for (int y=1; y<=Np; y++) {i++;
		pos_x[i] = xpos+(x-1)*D; 
		pos_y[i] = ypos+(y-1)*D;
		pos_z[i] = R;
	}  
	n_layers_x = (Np-1)*D-D/2;
	n_layers_y=n_layers_x;
	n_layers_z=n_layers_x;
	M=(n_layers_x+2)*(n_layers_y+2)*(n_layers_z+2);
	Text dim;
	dim = Blanks(9);
	dim.Putint(n_layers_z);
	dim = Copy(dim.Strip());
	dim = Copy(dim.Frontstrip());
	Sysout().Outtext("n_layers = " + dim); Sysout().Outimage(); 
	
	coordinate.Dim(1,M);

	int j=0;
	for (int x=0; x<=n_layers_x+1; x++)
	for (int y=0; y<=n_layers_y+1; y++)
	for (int z=0; z<=n_layers_z+1; z++){
		j++;
		coordinate[j]=0;
		for (int i=1; i<=Npt; i++){
			if ((x-pos_x[i])*(x-pos_x[i])+ (y-pos_y[i])*(y-pos_y[i])+(z-pos_z[i])*(z-pos_z[i])<=R*R) coordinate[j]=1;
		}
	}

	j=0;
	for (int x=0; x<=n_layers_x+1; x++)
	for (int y=0; y<=n_layers_y+1; y++)
	for (int z=0; z<=n_layers_z+1; z++){
		j++;
		if (coordinate[j]==1) coordinate[j]=0; else coordinate[j]=1;
		if (z>H) coordinate[j] = 0;
	}

	
	Out= new Outfile(Filename);
	Out -> Open(Blanks(1000));
	j=0;
	int k=0;
	for (int x=0; x<=n_layers_x+1; x++)
	for (int y=0; y<=n_layers_y+1; y++)
	for (int z=0; z<=n_layers_z+1; z++){
		j++;
		if (coordinate[j]==1) k=1; else k=0;
		if (x<1 || x>n_layers_x || y<1 || y>n_layers_y || z<1 || z>n_layers_z) k=0;
		Out -> Outint(k,1); Out -> Outimage();
	}
	Out -> Close();
	
};
