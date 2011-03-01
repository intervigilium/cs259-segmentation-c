//#include "DoubleArray3D.h"

#include  "TwoPhase3D.h"

//#include <iostream>
//#include <sstream>
//#include <iomanip>
//using namespace std;

//#include <fstream>

//
//#####################################################################
//                                                GET IMAGE INFO
//#####################################################################
//
void getDAT3DInfo(const char *filename, int *m, int *n, int *p)
{
	FILE *fp;
	if (fp = fopen(filename, "rb")) {
		fscanf(fp, "%d %d %d\n", m, n, p);
		fclose(fp);
	}

}

//
//#####################################################################
//                                                      READ DAT
//#####################################################################
//
void readDAT3D(double A[M][N][P], const char *filename)
{
	FILE *fp;
	long m, n, p;
	long i, j, k;

	if (fp = fopen(filename, "rb")) {
		fscanf(fp, "%ld %ld %ld\n", &m, &n, &p);

		for (k = 0; k < p; k++) {
			for (j = 0; j < n; j++) {
				for (i = 0; i < m; i++) {
					fscanf(fp, "%lf ", &A[i][j][k]);
				}
			}
		}

		fclose(fp);
	} else {
		printf("Error in reading file.\n");
		getchar();
	}
}

//
//#####################################################################
//                                                      WRITE DAT
//#####################################################################
//
//void writeDAT3D( const DoubleArray3D& A, string filename )
//{
//      ofstream outfile( filename.c_str() );
//      if(!outfile)
//      {
//              cout << "Error opening  " << filename << "!!!" << endl;
//              exit(-1);
//      }
//
//      outfile << A.getIndex1Size() << " " << A.getIndex2Size() << " " << A.getIndex3Size() << endl;
//
//    long i; long j; long k;
//      for(k = A.getIndex3Begin(); k <= A.getIndex3End(); k++)
//      {
//              for(j = A.getIndex2Begin(); j <= A.getIndex2End(); j++)
//              {
//                      for(i = A.getIndex1Begin(); i <= A.getIndex1End(); i++)
//                      {
//                              outfile <<  setw(5) << A(i,j,k) << " ";
//                      }
//                      outfile << endl;
//              }
//    }
//      outfile.close();
//}
