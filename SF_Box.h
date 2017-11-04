#ifndef SF_BOXxH
#define SF_BOXxH

#include <time.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#ifdef USE_FENV
#include <pthread.h>
#endif
#include "Input.h"
#include "Output.h"
#include "misc.h"

///
enum SegmentFreedom {loose, pinned, grafted, frozen};
///
enum MoleculeType {monomer, homopolymer, copolymer, polydisperse, dendrimer, asymmetric_dendrimer, branched, comb};
///
enum MolFreedom {fixedPhiBulk, fixedTheta,
	solvent, neutralizer, secondGeneration, thirdGeneration, rangeRestricted};
///
//enum DensityPart {total, unconstrained,  constrained, renorm, bulk, trainLoops, tails, unadsorbed, z_dir, y_dir, x_dir, yz_dir, longpath, shortpath};
enum DensityPart {total, unconstrained,  constrained, renorm, bulk, trainLoops, tails, unadsorbed, z_dir, y_dir, x_dir, yz_dir, longpath, shortpath};

///
extern Text SEGMENT;
///
extern Text STATE;
///
const double LOGMAXDOUBLE = log(DBL_MAX);
///
const double LOGMINDOUBLE = log(DBL_MIN);
///
const double PI = 4.0*atan(1.0);
///
const double BOLTZMANN = 1.38065812e-23;
///
extern double TEMPERATURE;
///
const double EPS0 = 8.85418785e-12;
///
const double ELEM_CHARGE = 1.6021773349e-19;
///
const double LOGCUTOFF = 100;
///
const double LOGTOLERANCE = 50;
///
const double LOGMAXPREFACTOR = LOGMAXDOUBLE - LOGCUTOFF - LOGTOLERANCE;
///
const double LOGZERO = -DBL_MAX/2;
///
const double CUTOFF = exp(LOGCUTOFF);
//@Include: *.h
///

void fpeon();
void fpeoff();
#endif
