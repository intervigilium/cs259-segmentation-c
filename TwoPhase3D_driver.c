//
//#####################################################################
//                                      TwoPhase3D_driver.cpp
//#####################################################################
//
//              Implementation of Active Contours - Two Phase in 3D
//                                              by Chan and Vese
//
//#####################################################################
// Igor Yanovsky (C) UCLA                               April 16, 2005
//#####################################################################
//


//#include <iostream>
//#include <sstream>
//#include <iomanip>
//using namespace std;
//
//#include <fstream>
//
#include <math.h>
#include <getopt.h>
#include "GridParameters.h"	// Parameter Classes
#include "RunParameters.h"
#include "TwoPhase3D.h"
#include "TwoPhase3DRoutines.h"	// Convection Routines
#include "IORoutines.h"		// IO for images
#include "papi.h"

#define PI 4.0*atan2(1.0,1.0)
#define IMG_FILE "8_4_nor_tal_nuc_ss.dat"

// Timing:
double realtime(void);
//
//#####################################################################
//                Local Function Declarations
//#####################################################################
//
//

void outputData(double phi[M][N][P], const GridParameters gridData,
		long stepCount);

double u0[M][N][P];

double phi[M][N][P];

unsigned batch_id = 0;

void help()
{
    printf("run the program\n");
}

int configure(int argc, char *argv[])
{
    char option_def[] = "?vh:b:";
    int ch;

    for (;;) {
	ch = getopt(argc, argv, option_def);
	if (ch == -1)
	    break;

	switch (ch) {
	case '?':
	case 'v':
	case 'h':
	    help();
	    return -1;
	case 'b':
	    batch_id = atoi(optarg);
	    printf("papi batch id is : %d \n", batch_id);
	    break;
	default:
	    printf("Unknown option detection\n");
	    return 0;
	}
    }

    return 0;
}

//
//#####################################################################
//                           Main Program
//#####################################################################
//
int main(int argc, char *argv[])
{
    //
    //*******************************************************
    //      Load an image:
    //*******************************************************
    //

    /*
       //
       // Test Image:
       //
       long m = 100;  long n = 80;  long p = 50;

       DoubleArray3D u0(m,n,p);

       long i;  long j;  long k;

       for(k=0; k < p; k++)
       {
       for(j=0; j < n; j++)
       {
       for(i=0; i < m; i++)
       {
       if( i>42 && i<58)
       if(j>32 && j<48)
       if(k>17 && k<33)
       u0(i,j,k) = 200.0;
       }
       }
       } */

    int m;
    int n;
    int p;
    GridParameters gridData;
    TwoPhase3D AC;
    double timeStart, timeTaken;
    RunParameters runData;

    long TimeSteps = 0;
    long outputCount = 0;

    double cont;
    double totalTime = 0.0;
    long kk = 0;
    char stemp[1024];

    if (configure(argc, argv) == -1)
	return -1;

    getDAT3DInfo(IMG_FILE, &m, &n, &p);	// get image dimensions

    //cout << endl << "Dimensions: " << m << " " << n << " " << p << endl;
    printf("\n Dimensions: %d %d %d \n", m, n, p);

    readDAT3D(u0, IMG_FILE);

    //
    //*******************************************************
    //  Specify the computational domain and the grid size
    //*******************************************************
    //

    gridData.m = m;		// m = the number of pixels in x-direction
    gridData.n = n;		// n = the number of pixels in y-direction
    gridData.p = p;

    gridData.dx = 1.0;
    gridData.dy = 1.0;
    gridData.dz = 1.0;

    gridData.xMin = 0.0;
    gridData.xMax = (m - 1) * gridData.dx;
    gridData.yMin = 0.0;
    gridData.yMax = (n - 1) * gridData.dy;
    gridData.zMin = 0.0;
    gridData.zMax = (p - 1) * gridData.dz;

    //
    //*******************************************************
    //  Declare Computational Variables and Initialize
    //*******************************************************
    //

    // AC.mu = 0.2 * 255 * 255;
    AC.mu = 0.18 * 255 * 255;
    AC.nu = 0.0;
    AC.lambda1 = 1.0;
    AC.lambda2 = 1.0;

    //cout << "Enter the Method        " << endl;
    //cout << "e =  Explicit           " << endl;
    //cout << "i =  Semi-Implicit      " << endl;
    //cin >> AC.method;

    printf("Enter the Method        \n");
    printf("e =  Explicit           \n");
    printf("i =  Semi-Implicit      \n");
    AC.method = 'e';
    printf("%c\n", AC.method);

    //cout << "Enter the Approximation to the Gradient        " << endl;
    //cout << "g =  Grad                                        " << endl;
    //cout << "d =  Dirac Delta                               " << endl;

    //cin >> AC.grad;

    printf("Enter the Approximation to the Gradient        \n");
    printf("g =  Grad  		                                 \n");
    printf("d =  Dirac Delta                               \n");
    AC.grad = 'g';
    printf("%c\n", AC.grad);

    initializePHI(phi, &gridData);

    //
    //*******************************************************
    //  Specify run and output parameters
    //*******************************************************
    //


    runData.epsilon = 10e-10;

    //cout << "Number of TimeSteps:     " << endl;
    //cin >> TimeSteps;
    printf("Number of TimeSteps:     \n");
    TimeSteps = 10;
    printf("%d\n", TimeSteps);

    //cout << "Output Every nth Step:   " << endl;
    //cin >> outputCount;
    printf("Output Every nth Step:   \n");
    outputCount = 1;
    printf("%d\n", outputCount);

    runData.rkOrder = 1;

    runData.dt = getCourantLimit(&gridData, &AC);

    if (AC.method == 'i') {
	runData.dt = 10e7 * runData.dt;
    }
    //
    //*******************************************************
    //  Output Parameters and Initial Conditions
    //*******************************************************
    //

    //cout << endl;
    ////cout << "Image: " << img_name << endl;
    //cout << "M  = " << gridData.m <<  " N =  " << gridData.n << " P =  " << gridData.p << endl;
    //cout << "dx = " << gridData.dx << " dy = " << gridData.dy << " dz = " << gridData.dz << endl;
    //cout << endl;
    //cout << "Time Stepping Parameters " << endl;
    //cout << "Number of Steps: " << TimeSteps << endl;
    //cout << "dt = " << runData.dt << " RK order : " << runData.rkOrder << endl;
    //cout << endl;
    //cout << "Active Contours CV 3D Parameters" << endl;
    //cout << "mu = " << AC.mu / (255*255) << " * 255 * 255, nu = " << AC.nu << ", lambda1 = " << AC.lambda1 << ", lambda2 = " << AC.lambda2 << endl;
    //cout << endl;
    //cout << "Enter any character to continue..." << endl;
    //cin >> cont;



    //  Plot initial data
    outputData(phi, gridData, 0);

    //
    //*******************************************************
    //  Initialize Time
    //*******************************************************
    //


    //cout << "Time : " << totalTime << endl;
    printf("Time : %lf \n", totalTime);


    //
    //*******************************************************
    //  Main Time Stepping Loop
    //*******************************************************
    //

    timeStart = realtime();

    int Events[5];
    u_long_long papi_values[5];
    util_start_papi(batch_id, Events);
    for (kk = 1; kk <= TimeSteps; kk++) {
	if (AC.method == 'e') {
	    advanceTimeStepEXPLICIT(phi, u0, &gridData, &runData, &AC);
	} else if (AC.method == 'i') {
	    advanceTimeStepIMPLICIT(phi, u0, &gridData, &runData, &AC);
	}
	//cout << "Step " << kk << endl;
	printf("Step %d\n", kk);
	totalTime += runData.dt;

	if ((kk % outputCount) == 0) {
	    //outputData(phi,gridData,kk);
	    outputData(phi, gridData, kk);
	}
	//cout << endl;
	//cout << "Time : " << totalTime << endl;
	printf("\nTime : %lf \n", totalTime);
    }
    util_stop_papi(batch_id, papi_values);
    util_print_papi(batch_id, papi_values, (batch_id == 0));
    timeTaken = realtime() - timeStart;

    //cout << "Filtering time  : " << timeTaken << " milliseconds.  (timeb.h)" << endl;
    //cout << "Filtering time  : " << timeTaken/1000 << " seconds.      (timeb.h)" << endl;
    printf("Filtering time : %lf milliseconds (timeb.h)\n", timeTaken);

    return 0;
}

void outputData(double phi[M][N][P], const GridParameters gridData,
		long stepCount)
{
    //
    //  Create ostringstream for construction output file names
    // 
    //ostringstream outs;
    char fileName[256];

    double sum = 0;

    long m = gridData.m;
    long n = gridData.n;
    long p = gridData.p;
    long i;
    long j;
    long k;

    //
    //  Open and then write to a file
    //
    FILE *dataFile;

    // 
    // Create file name of form uXXX.dat for step XXX 
    //
    //outs.str("");

    //outs << "phi" << stepCount << ".dat";

    //strcpy(fileName,(outs.str()).c_str());

    sprintf(fileName, "phi%d.dat", stepCount);


    if ((dataFile = fopen(fileName, "w+")) == NULL) {
	printf("The file %s could not be  opened\n", fileName);
	return;
    }
    //
    //  Output the data.
    //

    // x is the fastest-running variable:
    //
    for (k = 0; k < p; k++) {
	for (j = 0; j < n; j++) {
	    for (i = 0; i < m; i++) {
		fprintf(dataFile, "%-10.5e ", PPHI(i, j, k));
	    }
	    fprintf(dataFile, "\n");
	}
	fprintf(dataFile, "\n");
    }
    fclose(dataFile);

    for (k = 0; k < p; k++) {
	for (j = 0; j < n; j++) {
	    for (i = 0; i < m; i++) {
		sum += PPHI(i, j, k);
	    }
	}
    }
    printf("sum = %lf. \n", sum);


}
