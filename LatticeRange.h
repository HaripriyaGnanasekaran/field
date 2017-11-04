#ifndef LATTICERANGExH
#define LATTICERANGExH
#include "Message.h"
#include <stdlib.h>
#include <iostream>

///
class LatticeRange {

	public:
	virtual ~LatticeRange();
///
    virtual void MapVector(Vector, Vector) const = 0;
///
	virtual bool InRange(int) const = 0;
	virtual double GetRangeValue(int) const = 0;
///
	virtual int GetNumLayers(void) const = 0;
///
	virtual Text GetOutput(void) const = 0;
	virtual int GetNumPos(void) = 0;
	virtual double GetVolPos(void) = 0;
	virtual int GetPos(int) = 0;
	virtual bool SetPosLocal(int,double) = 0;
	virtual bool SetPos(int,int,int,double) = 0;
	virtual bool SetPos(int,int,double) = 0;
	virtual bool SetPos(double,double,double,double*) = 0;
	virtual bool ChangePos(int, int) = 0;
	virtual bool ClearAllPos() = 0;
//
//	virtual int GetNumLayers_x(void) const= 0;
//	virtual int GetNumLayers_y(void) const= 0;
//	virtual int GetNumLayers_z(void) const= 0;
//	virtual int Get_jx(void) const= 0;
//	virtual int Get_jy(void) const= 0;
//	virtual int Get_jz(void) const= 0;
//	virtual int Get_x0(void) const= 0;
//	virtual int Get_y0(void) const= 0;
//	virtual int Get_z0(void) const= 0;
//	virtual int Get_xm(void) const= 0;
//	virtual int Get_ym(void) const= 0;
//	virtual int Get_zm(void) const= 0;
//	virtual int Get_M(void) const= 0;
//

};

#endif
