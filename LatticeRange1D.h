#ifndef LATTICERANGE1DxH
#define LATTICERANGE1DxH
#include "LatticeRange.h"

///
class LatticeRange1D : public LatticeRange {
  public:
///
	LatticeRange1D(void);
///
	LatticeRange1D(int,int,int);
///
	virtual ~LatticeRange1D(void);
///
	int Dimensions(void) const;
    ///
    void MapVector(Vector, Vector) const;
///
	bool InRange(int) const;
	double GetRangeValue(int) const;
///
	int GetNumLayers(void) const;
///
	Text GetOutput(void) const;
	int GetNumLayers_x(void) const;
	int GetNumLayers_y(void) const;
	int GetNumLayers_z(void) const;
	int Get_jx(void) const;
	int Get_jy(void) const;
	int Get_jz(void) const;
	int Get_x0(void) const;
	int Get_y0(void) const;
	int Get_z0(void) const;
	int Get_xm(void) const;
	int Get_ym(void) const;
	int Get_zm(void) const;
	int Get_M(void) const;
	int GetNumPos(void);
	double GetVolPos(void);
	int GetPos(int);
	bool SetPosLocal(int,double);
	bool SetPos(double,double,double,double*);
	bool SetPos(int,int,int, double);
	bool SetPos(int,int,double);
	bool ChangePos(int, int);
	bool ClearAllPos();
  private:
///
	int lowX;
///
	int highX;
///
	int maxX;
///
	Text Range;
};

#endif
