#ifndef __IORoutines__
#define __IORoutines__

//#include "DoubleArray3D.h"



//void getDAT3DInfo( string filename, long& m, long& n, long& p );
void getDAT3DInfo(const char *filename, int *m, int *n, int *p);

//void readDAT3D( DoubleArray3D& A, string filename );
void readDAT3D(double A[M][N][P], const char *filename);


//void writeDAT3D( const DoubleArray3D& A, string filename );



#endif
