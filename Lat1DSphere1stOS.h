#ifndef LAT1DSPHERE1STOSxH
#define LAT1DSPHERE1STOSxH

#include "Lat1DSphere1stO.h"
#include "Lat1DCyl1stOS.h"

///
class Lat1DSphere1stOS : public virtual Lat1DSphere, virtual public Lat1DCyl1stOS {
  public:
///
	Lat1DSphere1stOS(Input*, Text);
///
	virtual ~Lat1DSphere1stOS();
};

#endif
