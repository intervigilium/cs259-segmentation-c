//
//#####################################################################
//                       RunParameters.h 
//
//         RunParameters Class For Two Phase 2D CV Model
//#####################################################################
//
//
#ifndef __RunParameters__
#define __RunParameters__

typedef struct {
    char scheme;		// Difference scheme identifier
    double dt;			// Time step size
    int rkOrder;		// Order of Runge-Kutta method

    double epsilon;
} RunParameters;


#endif
