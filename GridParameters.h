//
//#####################################################################
//                       GridParameters.h 
//
//             GridParameters Class for Level Sets
//#####################################################################
//
//
#ifndef __GridParameters__
#define __GridParameters__


typedef struct {
    double xMin;		// Computational Region is [xMin,xMax]X,[yMin,yMax]
    double xMax;
    double yMin;
    double yMax;
    double zMin;
    double zMax;
    long m;			// Number of Points in the x direction
    long n;			// Number of Points in the y direction
    long p;			// Number of Points in the z direction
    double dx;			// mesh width in x direction
    double dy;			// mesh widht in y direction
    double dz;			// mesh widht in z direction
} GridParameters;

#endif
