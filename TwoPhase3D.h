//
//#####################################################################
//                       TwoPhase3D.h 
//
//      TwoPhase3D Class for Active Contours without Edges
//
//#####################################################################
//
//
#ifndef __TwoPhase3D__
#define __TwoPhase3D__

#include "stdio.h"
#include "ctype.h"

//#include "DoubleArray3D.h"

#define M 256
#define N 61
#define P 256

#define PPHI(i, j, k) (phi[i][j][k])

typedef struct {

	double mu;
	double nu;
	double lambda1;
	double lambda2;

	double c1;		// mean intensities values
	double c2;

	char method;		// method identifier -- Explicit or Semi-Implicit
	char grad;		// gradient identifier -- Grad or DiracDelta
} TwoPhase3D;

#endif
