#ifndef LATTICERANGEGENERICxH
#define LATTICERANGEGENERICxH

#include "LatticeRange.h"

///
class LatticeRangeFile : public LatticeRange {
  public:
	LatticeRangeFile(void);
	LatticeRangeFile(Text _filename, int mx, int my, int mz);
	LatticeRangeFile(Text _filename, int my, int mz);
	LatticeRangeFile(Text _filename, int _numlayers);
	LatticeRangeFile(Array <bool> _mask, int _numlayers);
	virtual ~LatticeRangeFile(void);
       void MapVector(Vector, Vector) const;
	bool InRange(int) const;
	double GetRangeValue(int) const;
	int GetNumPos(void) ;
	double GetVolPos(void) ;
	int GetPos(int);
	bool SetPos(double,double,double,double*);
	bool SetPosLocal(int,double);
	bool SetPos(int,int,double);
	bool SetPos(int,int,int, double);
	void PastePos(int, double);

	bool ChangePos(int,int);
	bool ClearAllPos();
	int GetNumLayers(void) const;
	Text GetOutput(void) const;
	int dim,n_layers_x,n_layers_y,n_layers_z,jx,jy,jz;
protected:
	Text Range;
	Array <bool> Mask;
	Array <double> MaskValue;
	int NumLayers;
};

#endif
