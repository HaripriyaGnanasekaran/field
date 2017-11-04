#include <math.h>
#include <fenk.h>
Outfile *Out;
Text Filename_low,Filename_high;

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
	
	Filename_low = "low.mask";
	Filename_high = "high.mask";

	Sysout().Outtext("Diameter of Particles ");  Sysout().Breakoutimage(); Sysin().Inimage();
	R=Sysin().Inreal()/2;

	Sysout().Outtext("Center-center distance in a row (even number please)");  Sysout().Breakoutimage(); Sysin().Inimage();
	D=Sysin().Inreal();

	Sysout().Outtext("Number of Particles in each row");  Sysout().Breakoutimage(); Sysin().Inimage();
	Np=Sysin().Inint()+1;
	Npt=(2*Np);

	Sysout().Outtext("Cut-off height ");  Sysout().Breakoutimage(); Sysin().Inimage();
	H=Sysin().Inint();

	int i=0;
	double ypos, xpos;
	pos_x.Dim(1,Npt);
	pos_y.Dim(1,Npt);
	pos_z.Dim(1,Npt);

	xpos=1/2; ypos=1/2;
	int x=1;
	for (int y=1; y<=Np; y++) {i++;
		pos_x[i] = xpos+(x-1)*D; 
		pos_y[i] = ypos+(y-1)*D;
		pos_z[i] = R;
	}  

	xpos +=D/2; ypos +=D/2;
	for (int y=1; y<=Np; y++) {i++;
		pos_x[i] = xpos+(x-1)*D; 
		pos_y[i] = ypos+(y-1)*D;
		pos_z[i] = R;
	}  
	n_layers_x = D/2;
	n_layers_y= Np*D-3*D/2;
	n_layers_z= 3*D;
	M=(n_layers_x+2)*(n_layers_y+2)*(n_layers_z+2);
	Text dim_x;
	Text dim_y;
	Text dim_z;
	dim_x = Blanks(9);
	dim_y = Blanks(9);
	dim_z = Blanks(9);
	dim_x.Putint(n_layers_x);
	dim_y.Putint(n_layers_y);
	dim_z.Putint(n_layers_z);
	dim_x = Copy(dim_x.Strip());
	dim_y = Copy(dim_y.Strip());
	dim_z = Copy(dim_z.Strip());
	dim_x = Copy(dim_x.Frontstrip());
	dim_y = Copy(dim_y.Frontstrip());
	dim_z = Copy(dim_z.Frontstrip());
	Sysout().Outtext("n_layers_x = " + dim_x); Sysout().Outimage(); 
	Sysout().Outtext("n_layers_y = " + dim_y); Sysout().Outimage(); 
	Sysout().Outtext("n_layers_z = " + dim_z); Sysout().Outimage(); 
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

	
	Out= new Outfile(Filename_low);
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
	for (int i=1; i<=Npt; i++) {
		pos_x[i] =n_layers_x+1-pos_x[i];
		pos_y[i] =n_layers_y+1-pos_y[i];
		pos_z[i] =n_layers_z+1-pos_z[i];
	}
	
	j=0;
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
		if (z<n_layers_z+1-H) coordinate[j] = 0;
	}

	
	Out= new Outfile(Filename_high);
	Out -> Open(Blanks(1000));
	j=0;
	k=0;
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
