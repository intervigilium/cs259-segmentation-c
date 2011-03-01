/*
 * Segmentation kernel for FPGA implementation
 */

#include <stdlib.h>
#include <math.h>


#define PI 4.0*atan2(1.0,1.0)


double DiracDelta(double a)
{
    /* delta_2 approximation to delta, as in Chan-Vese, 2001 */
    double value;
    double eps = 1.0;
    /* as in Chan-Vese, 2001 */
    value = eps / ((double) (PI) * (a * a + eps * eps));
    return value;
}


void NeumannBC(double phi[M][N][P])
{

}


void getMeanIntensities(const double u0[M][N][P],
             const double phi[M][N][P], TwoPhase3D * AC)
{

}


void evaluateTwoPhase3DopExplicit(double phi[M][N][P],
                  const double u0[M][N][P],
                  const GridParameters * gridData,
                  const RunParameters * runData,
                  TwoPhase3D * AC, double dt,
                  double c1, double c2)
{
    /*
     * bandwidth dependant, stencil it
     */
    double mu = AC->mu;
    double nu = AC->nu;
    double lambda1 = AC->lambda1;
    double lambda2 = AC->lambda2;

    long m = gridData->m;
    long n = gridData->n;
    long p = gridData->p;
    double dx = gridData->dx;
    double dy = gridData->dy;
    double dz = gridData->dz;

    double dx2 = dx * 2.0;
    double dy2 = dy * 2.0;
    double dz2 = dz * 2.0;

    long i;
    long j;
    long k;

    double Grad;
    double K;           // curvature

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

    double curvature_motion_part[M][N][P];
    double stencil[3][3][3];

    for (i = 1; i < m - 1; i++) {
    for (j = 1; j < n - 1; j++) {
        for (k = 1; k < p - 1; k++) {

        /* stencil code goes here */
        stencil[0][0][0] = stencil[0][0][1];
        stencil[0][1][0] = stencil[0][1][1];
        stencil[0][2][0] = stencil[0][2][1];

        stencil[0][0][1] = stencil[0][0][2];
        stencil[0][1][1] = stencil[0][1][2];
        stencil[0][2][1] = stencil[0][2][2];
        
        stencil[0][0][2] = phi[i - 1][j - 1][k + 1];
        stencil[0][1][2] = phi[i - 1][j][k + 1];
        stencil[0][2][2] = phi[i - 1][j + 1][k + 1];

        stencil[1][0][0] = stencil[1][0][1];
        stencil[1][1][0] = stencil[1][2][1];
        stencil[1][2][0] = stencil[1][2][1];

        stencil[1][0][1] = stencil[1][0][2];
        stencil[1][1][1] = stencil[1][1][2];
        stencil[1][2][1] = stencil[1][2][2];

        stencil[1][0][2] = phi[i][j - 1][k + 1];
        stencil[1][1][2] = phi[i][j][k + 1];
        stencil[1][2][2] = phi[i][j + 1][k + 1];

        stencil[2][0][0] = stencil[2][0][1];
        stencil[2][1][0] = stencil[2][1][1];
        stencil[2][2][0] = stencil[2][2][1];

        stencil[2][0][1] = stencil[2][0][2];
        stencil[2][1][1] = stencil[2][1][2];
        stencil[2][2][1] = stencil[2][2][2];

        stencil[2][0][2] = phi[i + 1][j - 1][k + 1];
        stencil[2][1][2] = phi[i + 1][j][k + 1];
        stencil[2][2][2] = phi[i + 1][j + 1][k + 1];

        /* regular calculation here */
        Dx_p = (stencil[2][1][1] - stencil[1][1][1]) / dx;
        Dx_m = (stencil[1][1][1] - stencil[0][1][1]) / dx;
        Dy_p = (stencil[1][2][1] - stencil[1][1][1]) / dy;
        Dy_m = (stencil[1][1][1] - stencil[1][0][1]) / dy;
        Dz_p = (stencil[1][1][2] - stencil[1][1][1]) / dz;
        Dz_m = (stencil[1][1][1] - stencil[1][1][0]) / dz;

        Dx_0 = (stencil[2][1][1] - stencil[0][1][1]) / dx2;
        Dy_0 = (stencil[1][2][1] - stencil[1][0][1]) / dy2;
        Dz_0 = (stencil[1][1][2] - stencil[1][1][0]) / dz2;

        Dxx = (Dx_p - Dx_m) / dx;
        Dyy = (Dy_p - Dy_m) / dy;
        Dzz = (Dz_p - Dz_m) / dz;

        Dxy =
            (stencil[2][2][1] - stencil[2][0][1] -
             stencil[0][2][1] - stencil[0][0][1]) / (4 * dx * dy);
        Dxz =
            (stencil[2][1][2] - stencil[2][1][0] - 
             stencil[0][1][2] + stencil[0][1][0]) / (4 * dx * dz);
        Dyz =
            (stencil[1][2][2] - stencil[1][2][0] -
             stencil[1][0][2] + stencil[1][0][0]) / (4 * dy * dz);

        K = (Dx_0 * Dx_0 * Dyy - 2.0 * Dx_0 * Dy_0 * Dxy +
             Dy_0 * Dy_0 * Dxx + Dx_0 * Dx_0 * Dzz -
             2.0 * Dx_0 * Dz_0 * Dxz + Dz_0 * Dz_0 * Dxx +
             Dy_0 * Dy_0 * Dzz - 2.0 * Dy_0 * Dz_0 * Dyz +
             Dz_0 * Dz_0 * Dyy)
            / (pow(Dx_0 * Dx_0 + Dy_0 * Dy_0 + Dz_0 * Dz_0, 1.5) +
               runData->epsilon);

        if (AC->grad == 'g') {
            Grad = sqrt(Dx_0 * Dx_0 + Dy_0 * Dy_0 + Dz_0 * Dz_0);
        } else if (AC->grad == 'd') {
            Grad = DiracDelta(stencil[1][1][1]);
        }

        curvature_motion_part[i][j][k] =
            Grad * (mu * K +
                lambda1 * (u0[i][j][k] - c1) * (u0[i][j][k] -
                                c1)
                - lambda2 * (u0[i][j][k] - c2) * (u0[i][j][k] -
                                  c2));
        }
    }
    }

    /* inline this */
    NeumannBC(curvature_motion_part);

    /* pipeline this */
    for (i = 0; i < M; i++) {
        for (j = 0; j < N; j++) {
            for (k = 0; k < P; k++) {
                phi[i][j][k] += curvature_motion_part[i][j][k] * dt;
            }
        }
    }
}
