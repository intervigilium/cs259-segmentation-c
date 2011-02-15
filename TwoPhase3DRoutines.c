//
//#####################################################################
//                     TwoPhase3DRoutines.cpp
//#####################################################################
//
//		  Active Contours without Edges - Two Phase in 3D
//
//#####################################################################
// Igor Yanovsky (C) UCLA                               April 16, 2005
//#####################################################################
//
#include <math.h>
//#include <iostream>
//using namespace std;

//#include "DoubleArray1D.h"
//#include "DoubleArray3D.h"

#include "GridParameters.h"
#include "RunParameters.h"
#include "TwoPhase3D.h"

#define PI 4.0*atan2(1.0,1.0)
#define min(x,y) ((x<y)?x:y)

double curvature_motion_part[M][N][P];


// 
//#####################################################################
// 	           GET COURANT LIMIT
//#####################################################################
//
double getCourantLimit( const GridParameters * gridData, const TwoPhase3D * AC )
{
	//
	// Computes the Courant number for Active Contours CV
	//

	double mu = AC->mu;

	double dx = gridData->dx;
	double dy = gridData->dy;
	double dz = gridData->dz;

	double term1 = 2.0*mu /(dx*dx);
	double term2 = 2.0*mu /(dy*dy);
	double term3 = 2.0*mu /(dz*dz);

	double sum = term1 + term2 + term3;

	double dt = 100.0;

	if( sum > 1.0e-15 )
		dt = 0.9/sum;		// 90% CFL

	return dt;
}
//
//#####################################################################
//                     INITIALIZE PHI
//#####################################################################
//
void initializePHI(double phi[M][N][P], const GridParameters * gridData)
{
//
//  This routine inputs the initial condition for phi.
//
	double epsilon = 10e-10;

	long i; 
	long j;
	long k;

	double xcent = (gridData->m-1) / 2.0;
	double ycent = (gridData->n-1) / 2.0;
	double zcent = (gridData->p-1) / 2.0;
	double r = (min( min(gridData->m,gridData->n), gridData->p )) / 3.0;

	double xx; double yy; double zz;

	for(i=0; i < gridData->m; i++)
	{
		for(j=0; j < gridData->n; j++)
		{
			for(k=0; k < gridData->p; k++)
			{
				xx = (double)(i)*gridData->dx + gridData->xMin;
				yy = (double)(j)*gridData->dy + gridData->yMin;
				zz = (double)(k)*gridData->dz + gridData->zMin;
				PPHI(i,j,k) = sqrt( (xx-xcent)*(xx-xcent) + (yy-ycent)*(yy-ycent) + (zz-zcent)*(zz-zcent)) - r;
			}
		}
	}
}
// 
//#####################################################################
// 	           NEUMANN BC
//#####################################################################
//
void NeumannBC( double phi[M][N][P] )
{
	long m = M; //A.getIndex1Size();
	long n = N; //A.getIndex2Size();
	long p = P; //A.getIndex3Size();

	long i;  long j;  long k;

	for(j=0; j < n; j++)
	{
		for(k=0; k < p; k++)
		{
			PPHI( 0, j,k) = PPHI( 1, j,k);
			PPHI(m-1,j,k) = PPHI(m-2,j,k);
		}
	}

	for(i=0; i < m; i++)
	{
		for(k=0; k < p; k++)
		{
			PPHI(i, 0, k) = PPHI(i, 1, k);
			PPHI(i,n-1,k) = PPHI(i,n-2,k);
		}
	}

	for(i=0; i < m; i++)
	{
		for(j=0; j < n; j++)
		{
			PPHI(i,j, 0 ) = PPHI(i,j, 1 );
			PPHI(i,j,p-1) = PPHI(i,j,p-2);
		}
	}

	PPHI( 0 , 0 , 0 ) = PPHI( 1 , 1 , 1 );
	PPHI(m-1, 0 , 0 ) = PPHI(m-2, 1 , 1 );
	PPHI( 0 ,n-1, 0 ) = PPHI( 1 ,n-2, 1 );
	PPHI( 0 , 0 ,p-1) = PPHI( 1 , 1 ,p-2);
	PPHI(m-1,n-1, 0 ) = PPHI(m-2,n-2, 1 );
	PPHI(m-1, 0 ,p-1) = PPHI(m-2, 1 ,p-2);
	PPHI( 0 ,n-1,p-1) = PPHI( 1 ,n-2,p-2);
	PPHI(m-1,n-1,p-1) = PPHI(m-2,n-2,p-2);
}

// 
//#####################################################################
// 	           DIRAC DELTA
//#####################################################################
//

double DiracDelta( double a )
{
	// delta_2 approximation to delta, as in Chan-Vese, 2001

	double value;

	double eps = 1.0;

	// as in Chan-Vese, 2001
	value  =  eps / ( (double)(PI) * ( a*a + eps*eps ) );

	return value;
}
// 
//#####################################################################
//					GET MEAN INTENSITIES
//#####################################################################
//
void getMeanIntensities1( const double u0[M][N][P], const double phi[M][N][P], TwoPhase3D * AC )
{
	long m = M;//phi.getIndex1Size();
	long n = N;//phi.getIndex2Size();
	long p = P;//phi.getIndex3Size();

	double num1 = 0.0;
	double num2 = 0.0;
	double den1 = 0.0;
	double den2 = 0.0;

	// c1 is the average(u0) INSIDE  the zero contour
	// Need to sum up the values of u0 INSIDE the zero countour
	// and divide by the number of grid points INSIDE the zero contour.
	//
	// c2 is the average(u0) OUTSIDE the zero contour
	// Need to sum up the values of u0 OUTSIDE the zero countour
	// and divide by the number of grid points OUTSIDE the zero contour.

	long i;  long j;  long k;

	for(i=0; i < m; i++)
	{
		for(j=0; j < n; j++)
		{
			for(k=0; k < p; k++)
			{
				if(PPHI(i,j,k) < 0)
				{
					num1 = num1 + u0[i][j][k]; //(i,j,k);
					den1 = den1 + 1;
				}
				else if(PPHI(i,j,k) > 0)
				{
					num2 = num2 + u0[i][j][k]; //u0(i,j,k);
					den2 = den2 + 1;
				}
			}
		}
	}
	AC->c1 = num1/den1;
	AC->c2 = num2/den2;
}
// 
//#####################################################################
//			ADVANCE Chan-Vese Two Phase 3D IMPLICIT
//#####################################################################
//
void advanceTimeStepIMPLICIT( double phi[M][N][P], const double u0[M][N][P],
									   const GridParameters * gridData,
									   const RunParameters * runData, 
									   TwoPhase3D * AC )
{
	double eps = 10e-7;
	double eps_sqrd = eps*eps;

	double mu = AC->mu;
	double nu = AC->nu;
	double lambda1 = AC->lambda1;
	double lambda2 = AC->lambda2;

	double dt = runData->dt;

	long   m  = gridData->m;
	long   n  = gridData->n;
	long   p  = gridData->p;
	double dx = gridData->dx;
	double dy = gridData->dy;
	double dz = gridData->dz;

	double dx2 = dx*2.0;
	double dy2 = dy*2.0;
	double dz2 = dz*2.0;


	double c1;// = AC->c1;	// mean intensities
	double c2;// = AC->c2;

	long i; long j; long k;

	double Grad;

	double Dx_p;
	double Dx_m;
	double Dy_p;
	double Dy_m;
	double Dz_p;
	double Dz_m;
	double Dx_0;
	double Dy_0;
	double Dz_0;

	double C1x, C2x, C3y, C4y, C5z, C6z;
	double MM, CC, C1x_2x, C3y_4y, C5z_6z;


	getMeanIntensities1( u0, phi, AC );

	c1 = AC->c1;	// mean intensities
	c2 = AC->c2;


	for(i=1; i < m-1; i++)
	{
		for(j=1; j < n-1; j++)
		{
			for(k=1; k < p-1; k++)
			{
				Dx_p = (PPHI(i+1,j,k) - PPHI(i,j,k))/dx;
				Dx_m = (PPHI(i,j,k) - PPHI(i-1,j,k))/dx;
				Dy_p = (PPHI(i,j+1,k) - PPHI(i,j,k))/dy;
				Dy_m = (PPHI(i,j,k) - PPHI(i,j-1,k))/dy;
				Dz_p = (PPHI(i,j,k+1) - PPHI(i,j,k))/dz;
				Dz_m = (PPHI(i,j,k) - PPHI(i,j,k-1))/dz;

				Dx_0 = (PPHI(i+1,j,k) - PPHI(i-1,j,k))/dx2;
				Dy_0 = (PPHI(i,j+1,k) - PPHI(i,j-1,k))/dy2;
				Dz_0 = (PPHI(i,j,k+1) - PPHI(i,j,k-1))/dz2;

				C1x  =  1.0 / sqrt( Dx_p*Dx_p + Dy_0*Dy_0 + Dz_0*Dz_0  + eps_sqrd );
				C2x  =  1.0 / sqrt( Dx_m*Dx_m + Dy_0*Dy_0 + Dz_0*Dz_0  + eps_sqrd );
				C3y  =  1.0 / sqrt( Dx_0*Dx_0 + Dy_p*Dy_p + Dz_0*Dz_0  + eps_sqrd );
				C4y  =  1.0 / sqrt( Dx_0*Dx_0 + Dy_m*Dy_m + Dz_0*Dz_0  + eps_sqrd );
				C5z  =  1.0 / sqrt( Dx_0*Dx_0 + Dy_0*Dy_0 + Dz_p*Dz_p  + eps_sqrd);
				C6z  =  1.0 / sqrt( Dx_0*Dx_0 + Dy_0*Dy_0 + Dz_m*Dz_m  + eps_sqrd);
	        
				if( AC->grad == 'g' )
				{
					Grad = sqrt(Dx_0*Dx_0 + Dy_0*Dy_0 + Dz_0*Dz_0);
				}
				else if( AC->grad == 'd' )
				{
					Grad = DiracDelta( PPHI(i,j,k) );
				}

				MM  =  (dt/(dx*dy)) * Grad * mu;
				CC  =  1 + M*(C1x + C2x + C3y + C4y + C5z + C6z);
	        
				C1x_2x  =  C1x*PPHI(i+1,j,k) + C2x*PPHI(i-1,j,k);
				C3y_4y  =  C3y*PPHI(i,j+1,k) + C4y*PPHI(i,j-1,k);
				C5z_6z  =  C5z*PPHI(i,j,k+1) + C6z*PPHI(i,j,k-1);
	        
				curvature_motion_part[i][j][k] = (1.0 / CC) * ( PPHI(i,j,k) + MM*( C1x_2x + C3y_4y + C5z_6z ) 
										+ (dt*Grad)* ( (lambda1*(u0[i][j][k] - c1)*(u0[i][j][k] - c1)) 
													 - (lambda2*(u0[i][j][k] - c2)*(u0[i][j][k] - c2)) ) );
			}
		}
	}

	NeumannBC( curvature_motion_part );

	for (i=0;i<M;i++)
		for (j=0;j<N;j++)
		{
			for (k=0;k<P;k++)
			{
				PPHI(i,j,k) = curvature_motion_part[i][j][k];
			}
		}
}



// 
//#####################################################################
//			Evaluate Chan-Vese Two Phase 3D Op EXPLICIT
//#####################################################################
//
void evaluateTwoPhase3DopExplicit( double phi[M][N][P], const double u0[M][N][P],
											const GridParameters * gridData,
											const RunParameters * runData, 
											TwoPhase3D * AC , double dt)
{
	double mu = AC->mu;
	double nu = AC->nu;
	double lambda1 = AC->lambda1;
	double lambda2 = AC->lambda2;

	long   m  = gridData->m;
	long   n  = gridData->n;
	long   p  = gridData->p;
	double dx = gridData->dx;
	double dy = gridData->dy;
	double dz = gridData->dz;

	double dx2 = dx*2.0;
	double dy2 = dy*2.0;
	double dz2 = dz*2.0;

	long i; long j; long k;



	double Grad;

	double K;			// curvature

	double Dx_p;
	double Dx_m;
	double Dy_p;
	double Dy_m;
	double Dz_p;
	double Dz_m;
	double Dx_0;
	double Dy_0;
	double Dz_0;

	double Dxx, Dyy, Dzz;
	double Dxy, Dxz, Dyz;


	double c1;// = AC.c1;	// mean intensities
	double c2;// = AC.c2;


	getMeanIntensities1( u0, phi, AC );

	c1 = AC->c1;	// mean intensities
	c2 = AC->c2;


	for(i=1; i < m-1; i++)
	{
		for(j=1; j < n-1; j++)
		{
			for(k=1; k < p-1; k++)
			{

				Dx_p = (PPHI(i+1,j,k) - PPHI(i,j,k))/dx;
				Dx_m = (PPHI(i,j,k) - PPHI(i-1,j,k))/dx;
				Dy_p = (PPHI(i,j+1,k) - PPHI(i,j,k))/dy;
				Dy_m = (PPHI(i,j,k) - PPHI(i,j-1,k))/dy;
				Dz_p = (PPHI(i,j,k+1) - PPHI(i,j,k))/dz;
				Dz_m = (PPHI(i,j,k) - PPHI(i,j,k-1))/dz;

				Dx_0 = (PPHI(i+1,j,k) - PPHI(i-1,j,k))/dx2;
				Dy_0 = (PPHI(i,j+1,k) - PPHI(i,j-1,k))/dy2;
				Dz_0 = (PPHI(i,j,k+1) - PPHI(i,j,k-1))/dz2;

				Dxx = (Dx_p - Dx_m) / dx;
				Dyy = (Dy_p - Dy_m) / dy;
				Dzz = (Dz_p - Dz_m) / dz;

				Dxy = ( PPHI(i+1,j+1,k) - PPHI(i+1,j-1,k) - PPHI(i-1,j+1,k) + PPHI(i-1,j-1,k) ) / (4*dx*dy);
				Dxz = ( PPHI(i+1,j,k+1) - PPHI(i+1,j,k-1) - PPHI(i-1,j,k+1) + PPHI(i-1,j,k-1) ) / (4*dx*dz);
				Dyz = ( PPHI(i,j+1,k+1) - PPHI(i,j+1,k-1) - PPHI(i,j-1,k+1) + PPHI(i,j-1,k-1) ) / (4*dy*dz);

				K  = (  Dx_0*Dx_0*Dyy - 2.0*Dx_0*Dy_0*Dxy + Dy_0*Dy_0*Dxx
					+ Dx_0*Dx_0*Dzz - 2.0*Dx_0*Dz_0*Dxz + Dz_0*Dz_0*Dxx
					+ Dy_0*Dy_0*Dzz - 2.0*Dy_0*Dz_0*Dyz + Dz_0*Dz_0*Dyy )
					/ ( pow(Dx_0*Dx_0 + Dy_0*Dy_0 + Dz_0*Dz_0, 1.5) + runData->epsilon );

				if( AC->grad == 'g' )
				{
					Grad = sqrt(Dx_0*Dx_0 + Dy_0*Dy_0 + Dz_0*Dz_0);
				}
				else if( AC->grad == 'd' )
				{
					Grad = DiracDelta( PPHI(i,j,k) );
				}

				curvature_motion_part[i][j][k] = 
					Grad * ( mu * K + lambda1*(u0[i][j][k]-c1)*(u0[i][j][k]-c1) 
									- lambda2*(u0[i][j][k]-c2)*(u0[i][j][k]-c2) );
			}
		}
	}

	NeumannBC( curvature_motion_part );


	for (i=0;i<M;i++)
		for (j=0;j<N;j++)
		{
			for (k=0;k<P;k++)
			{
				PPHI(i,j,k) += curvature_motion_part[i][j][k] * dt;
			}
		}

	return ;
}
// 
//#####################################################################
// 	              Advance Time Step EXPLICIT
//#####################################################################
//

void advanceTimeStepEXPLICIT( double phi[M][N][P], const double u0[M][N][P],
					  const GridParameters * gridData, const RunParameters * runData,
					  TwoPhase3D * AC )
{
   double dt       = runData->dt;
   long rkOrder    = runData->rkOrder;
//
// The method is a variable order Runge-Kutta
// time-stepping scheme, the order of which is selected by rkOrder.
//
// Obtain Appropriate  Runge-Kutta Coefficients
//
	//DoubleArray1D  rk(4);
	//DoubleArray1D rkk(4);
	//switch(rkOrder)
	//	{
	//	case 1:  rk(0)  = 0.0;
	//		     rkk(0) = 1.0;
	//	break;
	//	case 2:  rk(0)  = 0.0; rk(1)  = 1.0;
	//		     rkk(0) = 0.5; rkk(1) = 0.5;
 //       break;
	//	case 3:  rk(0)  = 0.0;   rk(1)  =  1.0/3.0;  rk(2)  = 2.0/3.0;
 //                rkk(0) =.25;    rkk(1) =  0.0;     rkk(2) =.75;
	//	break;
 //       case 4:  rk(0)  = 0.0;    rk(1)  = 1.0/2.0;rk(2)  = 1.0/2.0;rk(3) = 1.0;
	//		    rkk(0) = 1.0/6.0;rkk(1)  = 1.0/3.0;rkk(2) = 1.0/3.0;rkk(3)= 1.0/6.0;
 //       break;
	//	default: cout << "Incorrect Specification of RK order " << endl;
 //            exit(1);
	//	break;
 //       } 
//
//  Declare and Initialize temporary variables
//
	//long m = gridData->m;
	//long n = gridData->n;
	//long p = gridData->p;

	//DoubleArray3D phin(m,n,p);
	//DoubleArray3D phinp1(m,n,p);
	//DoubleArray3D phistar(m,n,p);
	//DoubleArray3D FPHI(m,n,p);

	//phin     = phi;
	//phinp1   = phi;
	//phistar  = phi;
//
// RK loop - loop over stages
//
	//long i;
	//for (i=0; i < rkOrder; i++)
	//{
	//	FPHI    = evaluateTwoPhase3DopExplicit( phistar, u0, gridData, runData, AC );
	//	phinp1  = phinp1 + dt*rkk(i)*FPHI;
	//	//if((i+1) < rkOrder)
	//	//{
	//	//	phistar  = phin + dt*rk(i+1)*FPHI;
	//	//}
	//}
	//phi = phinp1;

	evaluateTwoPhase3DopExplicit( phi, u0, gridData, runData, AC, dt );
}

