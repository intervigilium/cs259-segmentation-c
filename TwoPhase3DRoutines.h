//
//#####################################################################
//                     TwoPhase3DRoutines.h
//#####################################################################
//
//              Active Contours without Edges - Two Phase in 3D,
//                                Routine Function Prototypes 
// 
//#####################################################################
// Igor Yanovsky (C) UCLA                               April 17, 2005
//#####################################################################
//
//#include "DoubleArray3D.h"
#include "GridParameters.h"
#include "RunParameters.h"


void advanceTimeStepEXPLICIT(double phi[M][N][P], const double u0[M][N][P],
			     const GridParameters * gridData,
			     const RunParameters * runData,
			     TwoPhase3D * AC);

void advanceTimeStepIMPLICIT(double phi[M][N][P], const double u0[M][N][P],
			     const GridParameters * gridData,
			     const RunParameters * runData,
			     TwoPhase3D * AC);

double getCourantLimit(const GridParameters * gridData,
		       const TwoPhase3D * AC);

void initializePHI(double phi[M][N][P], const GridParameters * gridData);
