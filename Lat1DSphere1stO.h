#ifndef LAT1DSPHERE1STOxH
#define LAT1DSPHERE1STOxH

#include "Lat1DSphere.h"
#include "Lat1DCyl1stO.h"
#include "Input.h"

///
class Lat1DSphere1stO : virtual public Lat1DCyl1stO, public Lat1DSphere {
  public:
///
	Lat1DSphere1stO(Input*, Text);
///
	virtual ~Lat1DSphere1stO();
};

#endif
