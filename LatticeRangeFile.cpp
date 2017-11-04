#include "LatticeRangeFile.h"
#include <math.h>
#include <cstdlib>
#include <iostream>



LatticeRangeFile::LatticeRangeFile() {
};

LatticeRangeFile::~LatticeRangeFile() {
};

double random_d2() {
	return rand()/(RAND_MAX+1.0);
}

int random_int2(int low, int high) {
	int range = high-low+1;
	return low+int(range*random_d2());
}

LatticeRangeFile::LatticeRangeFile(Array <bool> _mask, int _numlayers) {
	NumLayers = _numlayers;
	Mask.Dim(1, NumLayers);
	MaskValue.Dim(1,NumLayers);
	for (int z=1; z<=NumLayers; z++) { Mask[z] = _mask[z];
		if (Mask[z]) MaskValue[z]=1.0; else MaskValue[z] = 0.0;
	}
}

LatticeRangeFile::LatticeRangeFile(Text rangetext, int mx, int my, int mz) {

	NumLayers = (mx+2)*(my+2)*(mz+2);
	jx = (my+2)*(mz+2);
	jy = (mz+2);
	jz = 1;
	dim =3;
	char test = '?';
	n_layers_x=mx;
	n_layers_y=my;
	n_layers_z=mz;
	bool randomposx, randomposy, randomposz;
	bool input_entered=false;
	int x_value,y_value,z_value,D_value=0,x2_value=0,y2_value=0,z2_value=0;
	Text xValue, yValue, zValue,DValue;
	int posnr=0;
	if (rangetext.Getchar()=='{') {//pore in z-direction and cylindrical particles ={x1,y1,z1}[z2,R] and  ={x1,y1,z1}[z2,-R] respectively.
		input_entered=true;
		Mask.Dim(1,NumLayers);
		MaskValue.Dim(1,NumLayers);
		for (int z=1; z<=NumLayers; z++) {Mask[z] =false; MaskValue[z]=0.0; }

		if (!rangetext.More()) Message(fatal,"use ={x1,y1,z1}[z2,R] or  ={x1,y1,z1}[z2,-R] for cylinder or pore, respectively" );
		xValue = Copy(rangetext.Scanto(','));
		xValue = Copy(xValue.Strip());
		xValue = Copy(xValue.Frontstrip());
		if (!rangetext.More()) Message(fatal,"use ={x1,y1,z1}[z2,R] or  ={x1,y1,z1}[z2,-R] for cylinder or pore, respectively" );
		yValue = Copy(rangetext.Scanto(','));
		yValue = Copy(yValue.Strip());
		yValue = Copy(yValue.Frontstrip());
		if (!rangetext.More()) Message(fatal,"use ={x1,y1,z1}[z2,R] or  ={x1,y1,z1}[z2,-R] for cylinder or pore, respectively" );
		zValue = Copy(rangetext.Scanto('}'));
		zValue = Copy(zValue.Strip());
		zValue = Copy(zValue.Frontstrip());
		int length = xValue.Length();
		for (int i=1; i<= length; i++) {
			test = xValue.Getchar();
			if (!isdigit(test)) Message(fatal,"Pos of pore/cyl error: x-value not a digit. Use ={x1,y1,z1}[z2,R] or  ={x1,y1,z1}[z2,-R] for cylinder or pore, respectively ");
		}
		xValue.Setpos(1);
		x_value = xValue.Getint();
		if (x_value < 1) Message(fatal,"Pos of pore/cyl error: x-value below lowerbound ");
		if (x_value > mx) Message(fatal,"Pos of pore/cyl error: x-value above upperbound ");

		length = yValue.Length();
		for (int i=1; i<= length; i++) {
			test = yValue.Getchar();
			if (!isdigit(test)) Message(fatal,"Pos of pore/cyl error: y-value not a digit. Use ={x1,y1,z1}[z2,R] or  ={x1,y1,z1}[z2,-R] for cylinder or pore, respectively ");
		}
		yValue.Setpos(1);
		y_value = yValue.Getint();
		if (y_value < 1) Message(fatal,"Pos of pore/cyl error: y-value below lowerbound ");
		if (y_value > my) Message(fatal,"Pos of pore/cyl error: y-value above upperbound ");

		length = zValue.Length();
		for (int i=1; i<= length; i++) {
			test = zValue.Getchar();
			if (!isdigit(test)) Message(fatal,"Pos of pore/cyl error: z-value not a digit. Use ={x1,y1,z1}[z2,R] or  ={x1,y1,z1}[z2,-R] for cylinder or pore, respectively ");
		}
		zValue.Setpos(1);
		z_value = zValue.Getint();
		//Message(literal, "x = " + xValue + " y = " + yValue + " z = " +zValue);
		if (z_value < 1) Message(fatal,"Pos of pore/cyl error: z-value below lowerbound ");
		if (z_value > mz) Message(fatal,"Pos of pore/cyl error: z-value above upperbound ");
		if (!rangetext.More()) Message(fatal,"use ={x1,y1,z1}[z2,R] or  ={x1,y1,z1}[z2,-R] for cylinder or pore, respectively" );
		if (!rangetext.Getchar()=='[') Message(fatal,"use ={x1,y1,z1}[z2,R] or  ={x1,y1,z1}[z2,-R] for cylinder or pore, respectively" );
		if (!rangetext.More()) Message(fatal,"use ={x1,y1,z1}[z2,R] or  ={x1,y1,z1}[z2,-R] for cylinder or pore, respectively" );
		zValue = Copy(rangetext.Scanto(','));
		zValue = Copy(zValue.Strip());
		zValue = Copy(zValue.Frontstrip());
		length = zValue.Length();
		for (int i=1; i<= length; i++) {
			test = zValue.Getchar();
			if (!isdigit(test)) Message(fatal,"Pos of pore/cyl error: second z-value not a digit ");
		}
		zValue.Setpos(1);
		z2_value = zValue.Getint();
		if (z2_value < z_value) Message(fatal,"use ={x1,y1,z1}[z2,R] or  ={x1,y1,z1}[z2,-R] for cylinder or pore, respectively. z2 must be larger or equal to z1" );
		if (!rangetext.More()) Message(fatal,"use ={x1,y1,z1}[z2,R] or  ={x1,y1,z1}[z2,-R] for cylinder or pore, respectively" );
		DValue = Copy(rangetext.Scanto(']'));
		DValue = Copy(DValue.Strip());
		DValue = Copy(DValue.Frontstrip());
		length = DValue.Length();
		for (int i=1; i<= length; i++) {
			test = DValue.Getchar();
			if (!(isdigit(test) || test == '-')) Message(fatal,"Pos of pore/cyl error: R-value not a digit ");
		}
		DValue.Setpos(1);
		D_value = DValue.Getint();
		bool pore = D_value < 0;


		if (pore) D_value = -D_value;

		for (int z=z_value; z<=z2_value; z++) {
			for (int x=1; x<=mx; x++) {
				for (int y=1; y<=my; y++){
					double phi=0;
					if ( sqrt((x-x_value) *(x-x_value)+ (y-y_value)*(y-y_value)) <= D_value-1) {phi=1.0;
					} else if ( sqrt((x-x_value) *(x-x_value)+ (y-y_value)*(y-y_value)) <= D_value+1) {
						int inside=0;
						double xx,yy;
						for (int ex=1; ex<=100; ex++) for (int ey=1; ey<=100; ey++) {
							xx=1.0*x+1.0*ex/100.0-0.555;
							yy=1.0*y+1.0*ey/100.0-0.555;
							if ( sqrt((xx-1.0*x_value) *(xx-1.0*x_value)+ (yy-1.0*y_value)*(yy-1.0*y_value)) <= D_value)  {inside++;}
						}
						phi = 1.0*inside/10000.0;
					}

					if (pore) phi = 1-phi;
					if (phi>0) {
						if (!SetPos(jx*x+jy*y+jz*z+1,1.0,phi)) Message(fatal,"error: duplicate positions.??");
					}
				}
			}
		}
	} else {
		rangetext.Setpos(1);
		if (rangetext.Getchar()=='(') { // particle coordinates: points, spherical and ellipsoidal particles
			Mask.Dim(1,NumLayers);
			MaskValue.Dim(1,NumLayers);
			for (int z=1; z<=NumLayers; z++) {Mask[z] =false; MaskValue[z]=0.0; }
			input_entered=true;

			while (rangetext.More()) {
				input_entered=true;
				randomposx = false;
				randomposy = false;
				randomposz = false;
				posnr++;
				Text Posnumber = Blanks(5);
				Posnumber.Putint(posnr);
				Posnumber = Copy(Posnumber.Strip().Frontstrip());
				if (!rangetext.More()) Message(fatal,"Pos of particle " + Posnumber + " incomplete only x-value or ',' missing. For multi-pos entry: =(x1,y1,z1)(x2,y2,z2) etc.; For spheres entry: =(x1,y1,z1)[D] etc; For ellipsoids =(x1,y1,z1);(x2,y2,z2)[D] ");
				xValue = Copy(rangetext.Scanto(','));
				xValue = Copy(xValue.Strip());
				xValue = Copy(xValue.Frontstrip());
				if (!rangetext.More()) Message(fatal,"Pos of particle " + Posnumber + " incomplete only y-value or ',' missing.  For multi-pos entry: =(x1,y1,z1)(x2,y2,z2) etc.; For spheres entry: =(x1,y1,z1)[D] etc; For ellipsoids =(x1,y1,z1);(x2,y2,z2)[D] ");
				yValue = Copy(rangetext.Scanto(','));
				yValue = Copy(yValue.Strip());
				yValue = Copy(yValue.Frontstrip());
				if (!rangetext.More()) Message(fatal,"Pos of particle " + Posnumber + " incomplete missing z-value or ',' missing. For multi-pos entry: =(x1,y1,z1)(x2,y2,z2) etc.; For spheres entry: =(x1,y1,z1)[D] etc; For ellipsoids =(x1,y1,z1);(x2,y2,z2)[D] ");
				zValue = Copy(rangetext.Scanto(')'));
				zValue = Copy(zValue.Strip());
				zValue = Copy(zValue.Frontstrip());
				if (xValue.Length() == 0) Message(fatal,"Pos of particle " + Posnumber + " incomplete: x-value empty.  For multi-pos entry: =(x1,y1,z1)(x2,y2,z2) etc.; For spheres entry: =(x1,y1,z1)[D] etc; For ellipsoids =(x1,y1,z1);(x2,y2,z2)[D] ");
				if (yValue.Length() == 0) Message(fatal,"Pos of particle " + Posnumber + " incomplete: y-value empty.  For multi-pos entry: =(x1,y1,z1)(x2,y2,z2) etc.; For spheres entry: =(x1,y1,z1)[D] etc; For ellipsoids =(x1,y1,z1);(x2,y2,z2)[D] ");
				if (zValue.Length() == 0) Message(fatal,"Pos of particle " + Posnumber + " incomplete: z-value empty.  For multi-pos entry: =(x1,y1,z1)(x2,y2,z2) etc.; For spheres entry: =(x1,y1,z1)[D] etc; For ellipsoids =(x1,y1,z1);(x2,y2,z2)[D] ");
				int length = xValue.Length();
				for (int i=1; i<= length; i++) {
					test = xValue.Getchar();
					if (test == '*'){randomposx = true;}
					else {if (!isdigit(test)) Message(fatal,"Pos of particle " + Posnumber + " error: x-value not a digit. For multi-pos entry: =(x1,y1,z1)(x2,y2,z2) etc.; For spheres entry: =(x1,y1,z1)[D] etc; For ellipsoids =(x1,y1,z1);(x2,y2,z2)[D] ");
					}
				}
				xValue.Setpos(1);
				if (test == '*'){x_value = random_int2(1,mx);}
				else {x_value = xValue.Getint();}
				if (x_value < 1) Message(fatal,"Pos of particle " + Posnumber + " error: x-value below lowerbound.  For multi-pos entry: =(x1,y1,z1)(x2,y2,z2) etc.; For spheres entry: =(x1,y1,z1)[D] etc; For ellipsoids =(x1,y1,z1);(x2,y2,z2)[D] ");
				if (x_value > mx) Message(fatal,"Pos of particle " + Posnumber + " error: x-value above upperbound.  For multi-pos entry: =(x1,y1,z1)(x2,y2,z2) etc.; For spheres entry: =(x1,y1,z1)[D] etc; For ellipsoids =(x1,y1,z1);(x2,y2,z2)[D] ");

				length = yValue.Length();
				for (int i=1; i<= length; i++) {
					test = yValue.Getchar();
					if (test == '*'){randomposy = true;}
					else {if (!isdigit(test)) Message(fatal,"Pos of particle " + Posnumber + " error: y-value not a digit.  For multi-pos entry: =(x1,y1,z1)(x2,y2,z2) etc.; For spheres entry: =(x1,y1,z1)[D] etc; For ellipsoids =(x1,y1,z1);(x2,y2,z2)[D] ");
					}
				}
				yValue.Setpos(1);
				if (test == '*'){y_value = random_int2(1,my);}
				else {y_value = yValue.Getint();}
				if (y_value < 1) Message(fatal,"Pos of particle " + Posnumber + " error: y-value below lowerbound.  For multi-pos entry: =(x1,y1,z1)(x2,y2,z2) etc.; For spheres entry: =(x1,y1,z1)[D] etc; For ellipsoids =(x1,y1,z1);(x2,y2,z2)[D] ");
				if (y_value > my) Message(fatal,"Pos of particle " + Posnumber + " error: y-value above upperbound.  For multi-pos entry: =(x1,y1,z1)(x2,y2,z2) etc.; For spheres entry: =(x1,y1,z1)[D] etc; For ellipsoids =(x1,y1,z1);(x2,y2,z2)[D] ");

				length = zValue.Length();
				for (int i=1; i<= length; i++) {
					test = zValue.Getchar();
					if (test == '*'){randomposz = true;}
					else {if (!isdigit(test)) Message(fatal,"Pos of particle " + Posnumber + " error: z-value not a digit.  For multi-pos entry: =(x1,y1,z1)(x2,y2,z2) etc.; For spheres entry: =(x1,y1,z1)[D] etc; For ellipsoids =(x1,y1,z1);(x2,y2,z2)[D] ");
					}
				}
				zValue.Setpos(1);
				if (test == '*'){z_value = random_int2(1,mz);}
				else {z_value = zValue.Getint();}
				if (z_value < 1) Message(fatal,"Pos of particle " + Posnumber + " error: z-value below lowerbound ");
				if (z_value > mz) Message(fatal,"Pos of particle " + Posnumber + "error: z-value above upperbound ");
				if (rangetext.More()) {
					test = rangetext.Getchar();
					if (test=='(') D_value=0;
					if (test==';') {
						if (rangetext.Getchar()=='(') {
							if (!rangetext.More()) Message(fatal,"Pos of particle " + Posnumber + " incomplete only x2-value or ',' missing. For ellipsods =(x1,y1,z1);(x2,y2,z2)[D]  ");
							xValue = Copy(rangetext.Scanto(','));
							xValue = Copy(xValue.Strip());
							xValue = Copy(xValue.Frontstrip());
							if (!rangetext.More()) Message(fatal,"Pos of particle " + Posnumber + " incomplete only y2-value or ',' missing. For ellipsods =(x1,y1,z1);(x2,y2,z2)[D]  ");
							yValue = Copy(rangetext.Scanto(','));
							yValue = Copy(yValue.Strip());
							yValue = Copy(yValue.Frontstrip());
							if (!rangetext.More()) Message(fatal,"Pos of particle " + Posnumber + " incomplete missing z2-value or ',' missing. For ellipsods =(x1,y1,z1);(x2,y2,z2)[D]  ");
							zValue = Copy(rangetext.Scanto(')'));
							zValue = Copy(zValue.Strip());
							zValue = Copy(zValue.Frontstrip());
							if (xValue.Length() == 0) Message(fatal,"Pos of particle " + Posnumber + " incomplete: x2-value empty. For ellipsods =(x1,y1,z1);(x2,y2,z2)[D]  ");
							if (yValue.Length() == 0) Message(fatal,"Pos of particle " + Posnumber + " incomplete: y2-value empty. For ellipsods =(x1,y1,z1);(x2,y2,z2)[D]  ");
							if (zValue.Length() == 0) Message(fatal,"Pos of particle " + Posnumber + " incomplete: z2-value empty. For ellipsods =(x1,y1,z1);(x2,y2,z2)[D]  ");
							int length = xValue.Length();
							for (int i=1; i<= length; i++) {
								test = xValue.Getchar();
								if (test == '*'){randomposx = true;}
								else {if (!isdigit(test)) Message(fatal,"Pos of particle " + Posnumber + " error: x2-value not a digit. For ellipsods =(x1,y1,z1);(x2,y2,z2)[D]  ");
								}
							}
							xValue.Setpos(1);
							if (test == '*'){x2_value = random_int2(1,mx);}
							else {x2_value = xValue.Getint();}
							if (x2_value < 1) Message(fatal,"Pos of particle " + Posnumber + " error: x2-value below lowerbound. For ellipsods =(x1,y1,z1);(x2,y2,z2)[D]  ");
							if (x2_value > mx) Message(fatal,"Pos of particle " + Posnumber + " error: x2-value above upperbound. For ellipsods =(x1,y1,z1);(x2,y2,z2)[D]  ");

							length = yValue.Length();
							for (int i=1; i<= length; i++) {
								test = yValue.Getchar();
								if (test == '*'){randomposy = true;}
								else {if (!isdigit(test)) Message(fatal,"Pos of particle " + Posnumber + " error: y2-value not a digit. For ellipsods =(x1,y1,z1);(x2,y2,z2)[D]  ");
								}
							}
							yValue.Setpos(1);
							if (test == '*'){y2_value = random_int2(1,my);}
							else {y2_value = yValue.Getint();}
							if (y2_value < 1) Message(fatal,"Pos of particle " + Posnumber + " error: y2-value below lowerbound. For ellipsods =(x1,y1,z1);(x2,y2,z2)[D]  ");
							if (y2_value > my) Message(fatal,"Pos of particle " + Posnumber + " error: y2-value above upperbound. For ellipsods =(x1,y1,z1);(x2,y2,z2)[D]  ");

							length = zValue.Length();
							for (int i=1; i<= length; i++) {
								test = zValue.Getchar();
								if (test == '*'){randomposz = true;}
								else {if (!isdigit(test)) Message(fatal,"Pos of particle " + Posnumber + " error: z2-value not a digit. For ellipsods =(x1,y1,z1);(x2,y2,z2)[D]  ");
								}
							}
							zValue.Setpos(1);
							if (test == '*'){z2_value = random_int2(1,mz);}
							else {z2_value = zValue.Getint();}
							if (z2_value < 1) Message(fatal,"Pos of particle " + Posnumber + " error: z2-value below lowerbound. For ellipsods =(x1,y1,z1);(x2,y2,z2)[D]  ");
							if (z2_value > mz) Message(fatal,"Pos of particle " + Posnumber + "error: z2-value above upperbound. For ellipsods =(x1,y1,z1);(x2,y2,z2)[D]  ");
							test = rangetext.Getchar();
						}
					}
					if (test=='['){
						DValue = Copy(rangetext.Scanto(']')); if (rangetext.More()) rangetext.Scanto('(');
						length = DValue.Length();
						for (int i=1; i<= length; i++) {
							test = DValue.Getchar();
							if (!(isdigit(test) || test == '-')) Message(fatal,"Diameter of particle " + Posnumber + " error: D-value not a digit. For spheres entry: =(x1,y1,z1)[R] etc; For ellipsods (x1,y1,z1);(x2,y2,z2)[D]  ");
						}
						DValue.Setpos(1);
						D_value = DValue.Getint();
					}
				}
				bool cavity = D_value<0;
				if (cavity) D_value =-D_value;

				if (x2_value ==0) x2_value = x_value;
				if (y2_value ==0) y2_value = y_value;
				if (z2_value ==0) z2_value = z_value;

				if (D_value ==0) {
					if (x_value < 1 || x_value > mx || y_value < 1 || y_value > my || z_value < 1 || z_value > mz) Message(fatal,"Part of particle " + Posnumber + " outside box. ");
					while (!SetPos(jx*x_value+jy*y_value+jz*z_value+1,1.0,1.0))
						{if (randomposx == false && randomposy == false && randomposz == false) Message(fatal,"Pos of particle " + Posnumber + "error: duplicate positions");
						else {	if (randomposx == true){x_value = random_int2(1,mx);}
							if (randomposy == true){y_value = random_int2(1,my);}
							if (randomposz == true){z_value = random_int2(1,mz);}
						}
					}
				} else {
					for (int x=x_value-D_value;  x<= x_value+D_value; x++)
					for (int y=y_value-D_value;  y<= y_value+D_value; y++)
					for (int z=z_value-D_value;  z<= z_value+D_value; z++) {
						if ( sqrt((x-x_value) *(x-x_value)+ (y-y_value)*(y-y_value) +  (z-z_value)*(z-z_value))+
				     	    	 sqrt((x-x2_value)*(x-x2_value)+(y-y2_value)*(y-y2_value)+ (z-z2_value)*(z-z2_value)) <= D_value-1) {
							if (x < 1 || x > mx || y < 1 || y > my || z < 1 || z > mz) Message(fatal,"Part of particle " + Posnumber + " outside box. ");

							if (!cavity) { if (!SetPos(jx*x+jy*y+jz*z+1,1.0,1.0))
									{Message(fatal,"Pos of particle " + Posnumber + "error: duplicate positions");
								}
							}

						} else
						if ( sqrt((x-x_value) *(x-x_value)+ (y-y_value)*(y-y_value) +  (z-z_value)*(z-z_value))+
				     	    	 sqrt((x-x2_value)*(x-x2_value)+(y-y2_value)*(y-y2_value)+ (z-z2_value)*(z-z2_value)) <= D_value+1) {

							int inside=0;
							double xx,yy,zz;
							for (int ex=1; ex<=10; ex++) for (int ey=1; ey<=10; ey++) for (int ez=1; ez<=10; ez++) {
								xx=1.0*x+1.0*ex/10.0-0.55;
								yy=1.0*y+1.0*ey/10.0-0.55;
								zz=1.0*z+1.0*ez/10.0-0.55;
								if ( sqrt((xx-1.0*x_value) *(xx-1.0*x_value)+ (yy-1.0*y_value)*(yy-1.0*y_value) +  (zz-1.0*z_value)*(zz-1.0*z_value))+
				     	           	   	  sqrt((xx-1.0*x2_value)*(xx-1.0*x2_value)+(yy-1.0*y2_value)*(yy-1.0*y2_value)+ (zz-1.0*z2_value)*(z-1.0*z2_value)) <= D_value)  inside++;
							}
							if (inside > 0) {
								double phi = 1.0*inside/1000.0;
								if (x < 1 || x > mx || y < 1 || y > my || z < 1 || z > mz) Message(fatal,"Part of particle " + Posnumber + " outside box. ");
								if (cavity) phi  = 1-phi;
								if (!SetPos(jx*x+jy*y+jz*z+1,1,phi)) Message(fatal,"Pos of particle " + Posnumber + "error: duplicate positions");
 							} else {
								if (cavity) SetPos(jx*x+jy*y+jz*z+1,1.0,1.0);
							}

						} else {
							if (cavity) SetPos(jx*x+jy*y+jz*z+1,1.0,1.0);
						}
					}
				}
				x2_value = y2_value = z2_value = 0.0;
				//if (!SetPos(jx*x_value+jy*y_value+jz*z_value+1)) Message(fatal,"Pos of particle " + Posnumber + "error: duplicate positions");
				//if (rangetext.More()) rangetext.Scanto('(');
			}
		} else {
			rangetext.Setpos(1);
			char TValue;
			TValue = rangetext.Getchar();
			if ((TValue == '$') || (TValue == '<') || ( TValue =='>')) {
				// particle coordinates: points, can be non-integer coordinates...$ is a heck for putting a sphere at arbitray place. \ and / are for quadupole particles
				int QValue;
				if (TValue =='$') {QValue = 0;}
				if (TValue =='<') {QValue = 1;}
				if (TValue =='>') {QValue = -1;}
				Mask.Dim(1,NumLayers);
				MaskValue.Dim(1,NumLayers);
				for (int z=1; z<=NumLayers; z++) {Mask[z] =false; MaskValue[z]=0.0; }
				input_entered=true;
				double x_valueReal;
				double y_valueReal;
				double z_valueReal;

				while (rangetext.More()) {
					input_entered=true;
					posnr++;
					Text Posnumber = Blanks(5);
					Posnumber.Putint(posnr);
					Posnumber = Copy(Posnumber.Strip().Frontstrip());
					if (!rangetext.More()) Message(fatal,"Pos of particle " + Posnumber + " incomplete only x-value or ',' missing. For multi-pos entry: =(x1,y1,z1)(x2,y2,z2) etc.; For spheres entry: =(x1,y1,z1)[D] etc; For ellipsoids =(x1,y1,z1);(x2,y2,z2)[D] ");
					xValue = Copy(rangetext.Scanto(','));
					xValue = Copy(xValue.Strip());
					xValue = Copy(xValue.Frontstrip());
					if (!rangetext.More()) Message(fatal,"Pos of particle " + Posnumber + " incomplete only y-value or ',' missing.  For multi-pos entry: =(x1,y1,z1)(x2,y2,z2) etc.; For spheres entry: =(x1,y1,z1)[D] etc; For ellipsoids =(x1,y1,z1);(x2,y2,z2)[D] ");
					yValue = Copy(rangetext.Scanto(','));
					yValue = Copy(yValue.Strip());
					yValue = Copy(yValue.Frontstrip());
					if (!rangetext.More()) Message(fatal,"Pos of particle " + Posnumber + " incomplete missing z-value or ',' missing. For multi-pos entry: =(x1,y1,z1)(x2,y2,z2) etc.; For spheres entry: =(x1,y1,z1)[D] etc; For ellipsoids =(x1,y1,z1);(x2,y2,z2)[D] ");
					if (TValue =='$') zValue = Copy(rangetext.Scanto('$'));
					if (TValue =='<') zValue = Copy(rangetext.Scanto('<'));
					if (TValue =='>') zValue = Copy(rangetext.Scanto('>'));
					zValue = Copy(zValue.Strip());
					zValue = Copy(zValue.Frontstrip());
					if (xValue.Length() == 0) Message(fatal,"Pos of particle " + Posnumber + " incomplete: x-value empty.  For multi-pos entry: =(x1,y1,z1)(x2,y2,z2) etc.; For spheres entry: =(x1,y1,z1)[D] etc; For ellipsoids =(x1,y1,z1);(x2,y2,z2)[D] ");
					if (yValue.Length() == 0) Message(fatal,"Pos of particle " + Posnumber + " incomplete: y-value empty.  For multi-pos entry: =(x1,y1,z1)(x2,y2,z2) etc.; For spheres entry: =(x1,y1,z1)[D] etc; For ellipsoids =(x1,y1,z1);(x2,y2,z2)[D] ");
					if (zValue.Length() == 0) Message(fatal,"Pos of particle " + Posnumber + " incomplete: z-value empty.  For multi-pos entry: =(x1,y1,z1)(x2,y2,z2) etc.; For spheres entry: =(x1,y1,z1)[D] etc; For ellipsoids =(x1,y1,z1);(x2,y2,z2)[D] ");
					int length = xValue.Length();
					//for (int i=1; i<= length; i++) {
						//test = xValue.Getchar();
						//if (!isdigit(test)) Message(fatal,"Pos of particle " + Posnumber + " error: x-value not a digit. For multi-pos entry: =(x1,y1,z1)(x2,y2,z2) etc.; For spheres entry: =(x1,y1,z1)[D] etc; For ellipsoids =(x1,y1,z1);(x2,y2,z2)[D] ");
					//}
					xValue.Setpos(1);
					x_valueReal = xValue.Getreal();
					if (x_valueReal < 1) Message(fatal,"Pos of particle " + Posnumber + " error: x-value below lowerbound.  For multi-pos entry: =(x1,y1,z1)(x2,y2,z2) etc.; For spheres entry: =(x1,y1,z1)[D] etc; For ellipsoids =(x1,y1,z1);(x2,y2,z2)[D] ");
					if (x_valueReal > mx) Message(fatal,"Pos of particle " + Posnumber + " error: x-value above upperbound.  For multi-pos entry: =(x1,y1,z1)(x2,y2,z2) etc.; For spheres entry: =(x1,y1,z1)[D] etc; For ellipsoids =(x1,y1,z1);(x2,y2,z2)[D] ");

					length = yValue.Length();
					//for (int i=1; i<= length; i++) {
						//test = yValue.Getchar();
						//if (!isdigit(test)) Message(fatal,"Pos of particle " + Posnumber + " error: y-value not a digit.  For multi-pos entry: =(x1,y1,z1)(x2,y2,z2) etc.; For spheres entry: =(x1,y1,z1)[D] etc; For ellipsoids =(x1,y1,z1);(x2,y2,z2)[D] ");
					//}
					yValue.Setpos(1);
					y_valueReal = yValue.Getreal();
					if (y_valueReal < 1) Message(fatal,"Pos of particle " + Posnumber + " error: y-value below lowerbound.  For multi-pos entry: =(x1,y1,z1)(x2,y2,z2) etc.; For spheres entry: =(x1,y1,z1)[D] etc; For ellipsoids =(x1,y1,z1);(x2,y2,z2)[D] ");
					if (y_valueReal > my) Message(fatal,"Pos of particle " + Posnumber + " error: y-value above upperbound.  For multi-pos entry: =(x1,y1,z1)(x2,y2,z2) etc.; For spheres entry: =(x1,y1,z1)[D] etc; For ellipsoids =(x1,y1,z1);(x2,y2,z2)[D] ");

					length = zValue.Length();
					//for (int i=1; i<= length; i++) {
						//test = zValue.Getchar();
						//if (!isdigit(test)) Message(fatal,"Pos of particle " + Posnumber + " error: z-value not a digit.  For multi-pos entry: =(x1,y1,z1)(x2,y2,z2) etc.; For spheres entry: =(x1,y1,z1)[D] etc; For ellipsoids =(x1,y1,z1);(x2,y2,z2)[D] ");
					//}
					zValue.Setpos(1);
					z_valueReal = zValue.Getreal();
					if (z_valueReal < 1) Message(fatal,"Pos of particle " + Posnumber + " error: z-value below lowerbound ");
					if (z_valueReal > mz) Message(fatal,"Pos of particle " + Posnumber + "error: z-value above upperbound ");
					if (rangetext.More()) {
						test = rangetext.Getchar();
						if (test=='['){
							DValue = Copy(rangetext.Scanto(']')); if (rangetext.More()) rangetext.Scanto('(');
							length = DValue.Length();
							for (int i=1; i<= length; i++) {
								test = DValue.Getchar();
								if (!isdigit(test)) Message(fatal,"Diameter of particle " + Posnumber + " error: D-value not a digit. For spheres entry: =(x1,y1,z1)[R] etc; For ellipsods (x1,y1,z1);(x2,y2,z2)[D]  ");
							}
							DValue.Setpos(1);
							D_value = DValue.Getint();
						}
					}
//new
					double R_value=1.0*D_value/2.0;
					double phi=0;
					double phix=0;
					double phit=0;
					double phiP=0;
					double phiM=0;
					double dx,dy,dz;
					for (int x=x_valueReal-R_value-1.0; x<=x_valueReal+R_value+1.0; x++)
					for (int y=y_valueReal-R_value-1.0; y<=y_valueReal+R_value+1.0; y++)
					for (int z=z_valueReal-R_value-1.0; z<=z_valueReal+R_value+1.0; z++) {
						if (x < 1 || x > mx || y < 1 || y > my || z < 1 || z > mz) { //outside system; do not write pos.
						} else {
							phi=0; phiP=0; phiM=0;
							dx=1.0*x-x_valueReal;
							dy=1.0*y-y_valueReal;
							dz=1.0*z-z_valueReal;
							//double Qvalue=QValue*dx*dy*dz;
							if (sqrt(dx*dx+dy*dy+dz*dz)>R_value+1.5) { phi=0;//outside particle; do not write pos.
							} else {
								double xx=0,yy=0,zz=0;
								for (int ex=1; ex<=10; ex++) for (int ey=1; ey<=10; ey++) for (int ez=1; ez<=10; ez++) {
									xx=1.0*x+1.0*ex/10.0-0.55;
									yy=1.0*y+1.0*ey/10.0-0.55;
									zz=1.0*z+1.0*ez/10.0-0.55;
									if ( sqrt((xx-x_valueReal) *(xx-x_valueReal)+ (yy-y_valueReal)*(yy-y_valueReal) +  (zz-z_valueReal)*(zz-z_valueReal)) <= R_value)  {
										phi +=0.001;
										dx=xx-x_valueReal;
										dy=yy-y_valueReal;
										dz=zz-z_valueReal;
										if (dx*dy*dz>0) phiP +=0.001; else phiM +=0.001;
									}
								}

								switch (QValue) {
									case (0):
										break;
									case (1):
										phi = phiP*0.99;
										break;
									case (-1):
										phi = phiM*0.99;
										break;
								}
								phit +=phi;
								phix +=phi*x;
								if (!SetPos(jx*x+jy*y+jz*z+1,1,phi)) Message(fatal,"Pos of particle " + Posnumber + "error: duplicate positions");
							} //outside particle
						} //outside system
					} // for x , for y, for z
					Text ltext = Blanks(15);
					ltext.Putreal(phit,5);
					Message(literal,ltext);

//end new

//					for (int x=x_valueReal-D_value-1;  x<= x_valueReal+D_value+1; x++)
//					for (int y=y_valueReal-D_value-1;  y<= y_valueReal+D_value+1; y++)
//					for (int z=z_valueReal-D_value-1;  z<= z_valueReal+D_value+1; z++)
//					{
//						if ( sqrt((1.0*x-x_valueReal) *(1.0*x-x_valueReal)+ (1.0*y-y_valueReal)*(1.0*y-y_valueReal) +  (1.0*z-z_valueReal)*(1.0*z-z_valueReal))*2<= D_value-1) {
//							if (x < 1 || x > mx || y < 1 || y > my || z < 1 || z > mz) {} else {
//
//								if (!SetPos(jx*x+jy*y+jz*z+1,1.0,1.0))	{Message(fatal,"Pos of particle " + Posnumber + "error: duplicate positions");}
//							}
//
//						} else {
//							if ( sqrt((1.0*x-x_valueReal) *(1.0*x-x_valueReal)+ (1.0*y-y_valueReal)*(y-y_valueReal) +  (1.0*z-z_valueReal)*(1.0*z-z_valueReal))*2 <= D_value+1) {
//
//								int inside=0;
//								double xx,yy,zz;
//								for (int ex=1; ex<=10; ex++) for (int ey=1; ey<=10; ey++) for (int ez=1; ez<=10; ez++) {
//									xx=1.0*x+1.0*ex/10.0-0.55;
//									yy=1.0*y+1.0*ey/10.0-0.55;
//									zz=1.0*z+1.0*ez/10.0-0.55;
//									if ( sqrt((xx-x_valueReal) *(xx-x_valueReal)+ (yy-y_valueReal)*(yy-y_valueReal) +  (zz-z_valueReal)*(zz-z_valueReal))*2 <= D_value)  inside++;
//								}
//								if (inside > 0) {
//									double phi = 1.0*inside/1000.0;
//									if (x < 1 || x > mx || y < 1 || y > my || z < 1 || z > mz) {} else {
//										if (!SetPos(jx*x+jy*y+jz*z+1,1,phi)) Message(fatal,"Pos of particle " + Posnumber + "error: duplicate positions");
//									}
//								}
//							}
//						}

//					}
				}
			}
		}
	}

	if (!input_entered) {
		Message(fatal,"Valid inputs:  ={x1,y1,z1}[z2,R] or  ={x1,y1,z1}[z2,-R] for one cylinder or one pore, respectively; or For multi-pos entry: =(x1,y1,z1)(x2,y2,z2), etc.; For spheres: =(x1,y1,z1)[D], etc; For ellipsoids =(x1,y1,z1);(x2,y2,z2)[D], etc ");
	}
}



LatticeRangeFile::LatticeRangeFile(Text rangetext, int mx, int my) {
	bool on_edge = false; 
	NumLayers = (mx+2)*(my+2);
	jx = (my+2);
	jy = 1;
	jz = 0;
	int posnr=0;
	dim =2;
	n_layers_x=mx;
	n_layers_y=my;
	n_layers_z=0;
	if (rangetext.Getchar()=='(') {
		Mask.Dim(1,NumLayers);
		MaskValue.Dim(1,NumLayers);
		for (int z=1; z<=NumLayers; z++) {Mask[z] =false; MaskValue[z]=0.0; }
		int x_value,y_value,R_value=0;
		Text xValue, yValue,RValue;
		char test;
		while (rangetext.More()) {
			posnr++;
			Text Posnumber = Blanks(5);
			Posnumber.Putint(posnr);
			Posnumber = Copy(Posnumber.Strip().Frontstrip());
			if (!rangetext.More()) Message(fatal,"Pos of particle " + Posnumber + " incomplete only x-value or ',' missing ");
			xValue = Copy(rangetext.Scanto(','));
			xValue = Copy(xValue.Strip());
			xValue = Copy(xValue.Frontstrip());
			if (!rangetext.More()) Message(fatal,"Pos of particle " + Posnumber + " incomplete missing y-value or ',' missing ");
			yValue = Copy(rangetext.Scanto(')'));
			yValue = Copy(yValue.Strip());
			yValue = Copy(yValue.Frontstrip());
			if (xValue.Length() == 0) Message(fatal,"Pos of particle " + Posnumber + " incomplete: x-value empty ");
			if (yValue.Length() == 0) Message(fatal,"Pos of particle " + Posnumber + " incomplete: y-value empty ");

			int length = xValue.Length();
			for (int i=1; i<= length; i++) {
				test = xValue.Getchar();
				if (!isdigit(test)) Message(fatal,"Pos of particle " + Posnumber + " error: x-value not a digit ");
			}
			xValue.Setpos(1);
			x_value = xValue.Getint();
			if (x_value < 1) {
				if (x_value == 0) on_edge = true; 
				else Message(fatal,"Pos of particle " + Posnumber + " error: x-value below lowerbound ");
			}
			if (x_value > mx) Message(fatal,"Pos of particle " + Posnumber + " error: x-value above upperbound ");

			length = yValue.Length();
			for (int i=1; i<= length; i++) {
				test = yValue.Getchar();
				if (!isdigit(test)) Message(fatal,"Pos of particle " + Posnumber + " error: y-value not a digit ");
			}
			yValue.Setpos(1);
			y_value = yValue.Getint();
			if (y_value < 1) Message(fatal,"Pos of particle " + Posnumber + " error: y-value below lowerbound ");
			if (y_value > my) Message(fatal,"Pos of particle " + Posnumber + "error: y-value above upperbound ");
			if (rangetext.More()) {
				test = rangetext.Getchar();
				if (test=='(') R_value=0;
					else {
					RValue = Copy(rangetext.Scanto(']')); if (rangetext.More()) rangetext.Scanto('(');
					length = RValue.Length();
					for (int i=1; i<= length; i++) {
						test = RValue.Getchar();
						if (!isdigit(test)) Message(fatal,"Radius of particle " + Posnumber + " error: R-value not a digit ");
					}
					RValue.Setpos(1);
					R_value = RValue.Getint();
				}
			}
			if (on_edge) {
				x_value =1;  
				for (int x=1;  x<= x_value+R_value; x++)
				for (int y=y_value-R_value;  y<= y_value+R_value; y++) {
					if (x > mx || y < 1 || y > my) Message(fatal,"Part of particle " + Posnumber + " outside box. ");
					if ( (1.0*x-1.0*x_value)*(1.0*x-1.0*x_value) + (1.0*y-1.0*y_value)*(1.0*y-1.0*y_value) <= (1.0*R_value+0.5)*(1.0*R_value+0.5)) {
						int inside=0;
						double xx,yy;
						for (int ex=1; ex<=10; ex++) for (int ey=1; ey<=10; ey++) {
							xx=1.0*x+1.0*ex/10.0-0.55;
							yy=1.0*y+1.0*ey/10.0-0.55;
							if ( (xx-1.0*x_value) *(xx-1.0*x_value)+ (yy-1.0*y_value)*(yy-1.0*y_value) <= R_value*R_value)  inside++;
						}
						if (inside > 0) {
							double phi = 1.0*inside/100.0;
							if (!SetPos(jx*x+jy*y+1,1.0,phi)) Message(fatal,"Pos of particle " + Posnumber + "error: duplicate positions");
						}
					}
				}
			} else {
				for (int x=x_value-R_value;  x<= x_value+R_value; x++)
				for (int y=y_value-R_value;  y<= y_value+R_value; y++) {
					if (x < 1 || x > mx || y < 1 || y > my) Message(fatal,"Part of particle " + Posnumber + " outside box. ");
					if ( (1.0*x-1.0*x_value)*(1.0*x-1.0*x_value) + (1.0*y-1.0*y_value)*(1.0*y-1.0*y_value) <= (1.0*R_value+0.5)*(1.0*R_value+0.5)) {
						int inside=0;
						double xx,yy;
						for (int ex=1; ex<=100; ex++) for (int ey=1; ey<=100; ey++) {
							xx=1.0*x+1.0*ex/100.0-0.505;
							yy=1.0*y+1.0*ey/100.0-0.505;
							if ( (xx-1.0*x_value) *(xx-1.0*x_value)+ (yy-1.0*y_value)*(yy-1.0*y_value) <= R_value*R_value)  inside++;
						}
						if (inside > 0) {
							double phi = 1.0*inside/10000.0;
							if (!SetPos(jx*x+jy*y+1,1.0,phi)) Message(fatal,"Pos of particle " + Posnumber + "error: duplicate positions");
						}
					}
				}
			}
		}
	}
};




LatticeRangeFile::LatticeRangeFile(Text _filename, int _numLayers) {
	Infile in(_filename);
	NumLayers = _numLayers;
	if (!in.Open(Blanks(100))) {
		Message(fatal,"The range file '" + _filename + "' cannot be opened");
	};
	in.Inimage();
	if (in.Endfile()) {
		Message(fatal, "The range file '" + _filename + "' is empty");
	};
	Mask.Dim(1,NumLayers);
	MaskValue.Dim(1,NumLayers);
	Text line;
	for (int i=1; i<=NumLayers; i++){
		if (in.Endfile()) Message(fatal, "The range file '" + _filename + "' is shorter than expected");
		else line = Copy(in.Image.Strip());
		Mask[i] = ((line.Getint()==0)?0:1);
		if (Mask[i]) MaskValue[i] = 1.0; else MaskValue[i] = 0.0;
		in.Inimage();
	};
	in.Close();
};



double
LatticeRangeFile::GetVolPos()  {
	double vol=0;
	for (int z=1; z<=NumLayers; z++) vol = vol + MaskValue[z];
	return vol;
}

int
LatticeRangeFile::GetNumPos()  {
	int aantal=0;
	for (int z=1; z<=NumLayers; z++)
		if (Mask[z]) aantal++;
	return aantal;
}

bool
LatticeRangeFile::ClearAllPos()  {
	for (int z=1; z<=NumLayers; z++) {
		Mask[z]=false; MaskValue[z]=0.0;
	}
	return true;
}

bool
LatticeRangeFile::SetPosLocal(int pos, double waarde)  {
    if (waarde == 0 ) {Mask[pos] = false; MaskValue[pos] = waarde;  return true; }
    if (Mask[pos]) { return false; } //pos already taken
    else {
    	Mask[pos]=true; MaskValue[pos] = waarde; return true;
    }
}

void
LatticeRangeFile::PastePos(int pos, double waarde)  {
   Mask[pos]=true; MaskValue[pos] = waarde;
}

bool
LatticeRangeFile::SetPos(double rx,double ry,double rz,double* submask) {
	int x,y,z,posr;
	int xx,yy,zz;
	double distance=0,sum=0;

	x=int(rx);
	y=int(ry);
	z=int(rz);
	int l=0;
	for (int i=-1; i<2; i++)
	for (int j=-1; j<2; j++)
	for (int k=-1; k<2; k++) {
		distance=pow(rx-(x+i),2)+ pow(ry-(y+j),2)+pow(rz-(z+k),2);
		submask[l]*=exp(-distance); sum+=submask[l];
		l++;
	}
	l=0;
	for (int i=-1; i<2; i++)
	for (int j=-1; j<2; j++)
	for (int k=-1; k<2; k++) {
		submask[l]/=sum;
		xx=x+i;  if (xx>n_layers_x) xx = xx-n_layers_x; if (xx<1) xx = xx+n_layers_x;
		yy=y+j;  if (yy>n_layers_y) yy = yy-n_layers_y; if (yy<1) yy = yy+n_layers_y;
		zz=z+k;  if (zz>n_layers_z) zz = zz-n_layers_z; if (zz<1) zz = zz+n_layers_z;
		posr=jx*xx+jy*yy+jz*zz+1;
		SetPosLocal(posr,GetRangeValue(posr)+submask[l]);
		if (GetRangeValue(posr)>1) {
			Text Posnumber = Blanks(10);
			Posnumber.Putint(posr);
			Posnumber = Copy(Posnumber.Strip().Frontstrip());
			Message(fatal, "At position " + Posnumber + " the Range Value (in LatticeRangefile) is larger than unity. Contact Frans.");
		}
		l++;
	}
	return true;
}


bool
LatticeRangeFile::SetPos(int pos, int SpotType, int SpotSize, double waarde)  {

	if (SpotType ==2) return SetPos(pos, SpotSize, waarde);
	int position,pos_x,pos_y,pos_z;
	int xm,ym,zm;
	double r;
	position = pos -1;
	pos_z = ((position % jx) % jy)/jz;
	pos_y = ((position % jx)-pos_z*jz)/jy;
	pos_x = (position-pos_y*jy-pos_z*jz)/jx;



	for (int x=pos_x-SpotSize;  x<= pos_x+SpotSize; x++)
	for (int y=pos_y-SpotSize;  y<= pos_y+SpotSize; y++)
	for (int z=pos_z-SpotSize;  z<= pos_z+SpotSize; z++) {
		r = 2*sqrt(  (x-pos_x)*(x-pos_x)+ (y-pos_y)*(y-pos_y) +  (z-pos_z)*(z-pos_z)  );
		if ( r <= SpotSize-1) {
			if (x<1) xm=x+n_layers_x; else xm=x;
			if (x>n_layers_x) xm = x-n_layers_x; else xm=x;
			if (y<1) ym=y+n_layers_y; else ym=y;
			if (y>n_layers_y) ym = y-n_layers_y; else ym=y;
			if (z<1) zm=z+n_layers_z; else zm=z;
			if (z>n_layers_z) zm = z-n_layers_z; else zm=z;
			PastePos(jx*xm+jy*ym+jz*zm+1,1.0*waarde);
		} else if ( r <= SpotSize+1) {
			int inside=0;
			double xx,yy,zz;
			for (int ex=1; ex<=10; ex++) for (int ey=1; ey<=10; ey++) for (int ez=1; ez<=10; ez++) {
				xx=1.0*x+1.0*ex/10.0-0.55;
				yy=1.0*y+1.0*ey/10.0-0.55;
				zz=1.0*z+1.0*ez/10.0-0.55;
				if ( 2*sqrt((xx-pos_x) *(xx-pos_x)+ (yy-pos_y)*(yy-pos_y) +  (zz-pos_z)*(zz-pos_z)) <= SpotSize)  inside++;
			}
			if (inside > 0) {
				double phi = 1.0*inside/1000.0;
				if (x<1) xm=x+n_layers_x;else xm=x;
				if (x>n_layers_x) xm = x-n_layers_x; else xm=x;
				if (y<1) ym=y+n_layers_y;else ym=y;
				if (y>n_layers_y) ym = y-n_layers_y;else ym=y;
				if (z<1) zm=z+n_layers_z;else zm=z;
				if (z>n_layers_z) zm = z-n_layers_z;else zm=z;

				PastePos(jx*xm+jy*ym+jz*zm+1,1.0*phi);
 			}
		}
	}
    return true;
}


bool
LatticeRangeFile::SetPos(int pos, int SpotSize, double waarde)  {
    if (dim ==2) {
		int position,pos_x,pos_y,xx=0,yy=0;
		position = pos -1;
		pos_y = ((position % jx))/jy;
		pos_x = (position-pos_y*jy)/jx;
		SetPosLocal(pos,waarde);
		if (SpotSize>1) {xx = pos_x+1; if (xx>n_layers_x) xx = xx-n_layers_x; SetPosLocal(jx*xx+jy*pos_y+1,waarde); }// assume periodic boundary conditions
		if (SpotSize>2) {yy = pos_y+1; if (yy>n_layers_y) yy = yy-n_layers_y; SetPosLocal(jx*pos_x+jy*yy+1,waarde); }
		if (SpotSize>3) {xx = pos_x-1; if (xx<1) xx = xx+n_layers_x; SetPosLocal(jx*xx+jy*pos_y+1,waarde); }
    } else {
	int position,pos_x,pos_y,pos_z,xx=0,yy=0,zz=0;
	position = pos -1;
	pos_z = ((position % jx) % jy)/jz;
	pos_y = ((position % jx)-pos_z*jz)/jy;
	pos_x = (position-pos_y*jy-pos_z*jz)/jx;

	switch (SpotSize){
		case 1:
			SetPosLocal(pos,waarde);
			break;
		case 2:
			SetPosLocal(pos,waarde);
			xx = pos_x+1; if (xx>n_layers_x) xx = xx-n_layers_x; SetPosLocal(jx*xx+jy*pos_y+jz*pos_z+1,waarde);
			break;
		case 4:
			for (int x=0; x<=1; x++)
			for (int y=0; y<=1; y++) {
				xx = pos_x+x; if (xx>n_layers_x) xx = xx-n_layers_x;
				yy = pos_y+y; if (yy>n_layers_y) yy = yy-n_layers_y;
				SetPosLocal(jx*xx+jy*yy+jz*pos_z+1,waarde);
			}
			break;
		case 8:
			for (int x=0; x<=1; x++)
			for (int y=0; y<=1; y++)
			for (int z=0; z<=1; z++){
				xx = pos_x+x; if (xx>n_layers_x) xx = xx-n_layers_x;
				yy = pos_y+y; if (yy>n_layers_y) yy = yy-n_layers_y;
				zz = pos_z+z; if (zz>n_layers_z) zz = zz-n_layers_z;
				SetPosLocal(jx*xx+jy*yy+jz*zz+1,waarde);
			}
			//printf("pos=%d\n", pos);
			break;
		case 27:
			for (int x=0; x<=2; x++)
			for (int y=0; y<=2; y++)
			for (int z=0; z<=2; z++){
				xx = pos_x+x; if (xx>n_layers_x) xx = xx-n_layers_x;
				yy = pos_y+y; if (yy>n_layers_y) yy = yy-n_layers_y;
				zz = pos_z+z; if (zz>n_layers_z) zz = zz-n_layers_z;
				SetPosLocal(jx*xx+jy*yy+jz*zz+1,waarde);
			}
			break;
		default:
			Message(fatal,"Spotsize not supported in LatticeRangeFile");
			break;
	}

	//if (SpotSize>1) {xx = pos_x+1; if (xx>n_layers_x) xx = xx-n_layers_x; SetPosLocal(jx*xx+jy*pos_y+jz*pos_z+1,waarde); } //assume periodic boundary conditions
	//if (SpotSize>2) {yy = pos_y+1; if (yy>n_layers_y) yy = yy-n_layers_y; SetPosLocal(jx*pos_x+jy*yy+jz*pos_z+1,waarde); }
	//if (SpotSize>3) {zz = pos_z+1; if (zz>n_layers_z) zz = zz-n_layers_z; SetPosLocal(jx*pos_x+jy*pos_y+jz*zz+1,waarde); }
	//if (SpotSize>4) {xx = pos_x-1; if (xx<1) xx = xx+n_layers_x; SetPosLocal(jx*xx+jy*pos_y+jz*pos_z+1,waarde); }
	//if (SpotSize>5) {yy = pos_y-1; if (yy<1) yy = yy+n_layers_y; SetPosLocal(jx*pos_x+jy*yy+jz*pos_z+1,waarde); }
	//if (SpotSize>6) {zz = pos_z-1; if (zz<1) zz = zz+n_layers_z; SetPosLocal(jx*pos_x+jy*pos_y+jz*zz+1,waarde); }

    }
    return true;
}

int
LatticeRangeFile::GetPos(int pos_number)  {
    int aantal=0;
    int z=0;
    while (z<=NumLayers && aantal < pos_number) {
    	z++;
    	if (Mask[z]) aantal++;
    }
    return z;
}

bool
LatticeRangeFile::ChangePos(int pos_new, int pos_old){
	Message(fatal,"in LatticeRangeFile; ChangePos used; I thought it was only used in MC...; contact FL.");
	if ((pos_new<=NumLayers) && (pos_new>0) && (pos_old<=NumLayers) && (pos_old>0)) {
		Mask[pos_old] = false ; Mask[pos_new] = true;
		MaskValue[pos_old] = 0.0; MaskValue[pos_new] = 1.0;
		return true;
	}
	else return false;
}

void
LatticeRangeFile::MapVector(Vector out, Vector in) const {
	for (int i=1; i<=NumLayers; i++) {
		out[i] = Mask[i] * in[i];
	}
};

bool
LatticeRangeFile::InRange(int z) const {
	if ((z<=NumLayers) && (z>0)) {
		return Mask[z];
	}
	else printf("out of range z=%d\n", z);
	return false;
};

double
LatticeRangeFile::GetRangeValue(int z) const {
	if (Mask[z]) return MaskValue[z];
	return 0.0;
}

int
LatticeRangeFile::GetNumLayers() const {
	int n = 0;
	for (int z=1; z<=NumLayers; z++) n=n+Mask[z];
	printf("numlayers in range n=%d\n", n);
	return n;
};
Text
LatticeRangeFile::GetOutput() const {
	return "Mask";
};

