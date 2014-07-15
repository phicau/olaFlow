/*---------------------------------------------------------------------------*\
| Work developed at IH Cantabria       IIIII H   H FFFFF OOOOO AAAAA M   M    |
|                                        I   H   H F     O   O A   A MM MM    |
|   Coder: Pablo Higuera Caubilla        I   HHHHH FFFF  O   O AAAAA M M M    |
|   Bug reports: higuerap@unican.es      I   H   H F     O   O A   A M   M    |
|                                      IIIII H   H F     OOOOO A   A M   M    |
|   -----------------------------------------------------------------------   |
| References:                                                                 |
|                                                                             |
| - Realistic wave generation and active wave absorption for Navier-Stokes    |
|    models: Application to OpenFOAM.                                         |
|    Higuera, P., Lara, J.L. and Losada, I.J. (2013)                          |
|    Coastal Engineering, Vol. 71, 102-118.                                   |
|    http://dx.doi.org/10.1016/j.coastaleng.2012.07.002                       |
|                                                                             |
| - Simulating coastal engineering processes with OpenFOAM                    |
|    Higuera, P., Lara, J.L. and Losada, I.J. (2013)                          |
|    Coastal Engineering, Vol. 71, 119-134.                                   |
|    http://dx.doi.org/10.1016/j.coastaleng.2012.06.002                       |
|                                                                             |
\*---------------------------------------------------------------------------*/

#ifndef waveFun_H
#define waveFun_H

namespace otherFun
{
    double interpolation (double x1, double x2, double y1, double y2, double xInt);
}

namespace StokesIFun
{
    double waveLength (double h, double T);
    double eta (double H, double Kx, double x, double Ky, double y, double omega, double t, double phase);
    double U (double H, double h, double Kx, double x, double Ky, double y, double omega, double t, double phase, double z);
    double W (double H, double h, double Kx, double x, double Ky, double y, double omega, double t, double phase, double z);
}

namespace StokesIIFun
{
    double eta (double H, double h, double Kx, double x, double Ky, double y, double omega, double t, double phase);
    double U (double H, double h, double Kx, double x, double Ky, double y, double omega, double t, double phase, double z);
    double W (double H, double h, double Kx, double x, double Ky, double y, double omega, double t, double phase, double z);
    double timeLag (double H, double h, double Kx, double x, double Ky, double y, double T, double phase);
}

namespace Elliptic
{
    void ellipticIntegralsKE (double m, double* K, double* E);
    double JacobiAmp (double u, double m);
    void JacobiSnCnDn (double u, double m, double* Sn, double* Cn, double* Dn);
}

namespace cnoidalFun
{
    double eta (double H, double m, double kx, double ky, double T, double x, double y, double t);
    double etaCnoidal1D (double H, double m, double t, double T);
    double etaMeanSq (double H, double m, double T);
    double d1EtaDx (double H, double m, double uCnoidal, double L);
    double d2EtaDx (double H, double m, double uCnoidal, double L);
    double d3EtaDx (double H, double m, double uCnoidal, double L);
    double timeLag (double H, double m, double kx, double ky, double T, double x, double y, double phase);
    int calculations (double H, double d, double T, double* mOut, double* LOut);
    double U (double H, double h, double m, double kx, double ky, double T, double x, double y, double t, double z);
    double W (double H, double h, double m, double kx, double ky, double T, double x, double y, double t, double z);
}

namespace stokesVFun
{
    double A11 (double h, double k);
    double A13 (double h, double k);
    double A15 (double h, double k);
    double A22 (double h, double k);
    double A24 (double h, double k);
    double A33 (double h, double k);
    double A35 (double h, double k);
    double A44 (double h, double k);
    double A55 (double h, double k);
    double B22 (double h, double k);
    double B24 (double h, double k);
    double B33 (double h, double k);
    double B33k (double h, double k);
    double B35 (double h, double k);
    double B35k (double h, double k);
    double B44 (double h, double k);
    double B55 (double h, double k);
    double B55k (double h, double k);
    double C1 (double h, double k);
    double C1k (double h, double k);
    double C2 (double h, double k);
    double C2k (double h, double k);
    double C3 (double h, double k);
    double C4 (double h, double k);
    int StokesVNR (double H, double d, double T, double* kOut, double* LambdaOut, double* f1Out, double* f2Out );
    int StokesExtendedSolver (double H, double d, double T, double* kOut, double* LambdaOut, double* LambdaErrOut );
    double eta (double h, double kx, double ky, double lambda, double T, double x, double y, double t, double phase);
    double timeLag (double h, double kx, double ky, double lambda, double T, double x, double y, double phase);
    double phaseLag (double d, double k, double lambda, double T);
    double U (double d, double kx, double ky, double lambda, double T, double x, double y, double t, double phase, double z);
    double V (double d, double kx, double ky, double lambda, double T, double x, double y, double t, double phase, double z);
}

namespace stokesVFentonFun
{
    double C0 (double h, double k);
    double C2 (double h, double k);
    double C4 (double h, double k);
    double D2 (double h, double k);
    double D4 (double h, double k);
    double error (double H, double h, double k, double T);
    int StokesSolver (double H, double d, double T, double* kOut, double* errorOut );
    double B22 (double h, double k);
    double B31 (double h, double k);
    double B42 (double h, double k);
    double B44 (double h, double k);
    double B53 (double h, double k);
    double B55 (double h, double k);
    double eta (double H, double d, double kx, double ky, double T, double x, double y, double t, double phase);
}

namespace secondOrderFun
{

    double C (double sigma1, double sigma2, double alphaSO1, double alphaSO2);
    double E (double a1, double a2, double sigma1, double sigma2, double alphaSO1, double alphaSO2);
    double etaSO (double H1, double H2, double sigma1, double sigma2, double phase1, double phase2, double kx1, double ky1, double kx2, double ky2, double x, double y, double t, double h);
    double uSO (double H1, double H2, double sigma1, double sigma2, double phase1, double phase2, double kx1, double ky1, double kx2, double ky2, double x, double y, double t, double h, double z);
    double wSO (double H1, double H2, double sigma1, double sigma2, double phase1, double phase2, double kx1, double ky1, double kx2, double ky2, double x, double y, double t, double h, double z);
}

namespace BoussinesqFun
{
    double eta (double H, double h, double x, double y, double theta, double t, double X0);
    double Deta1 (double H, double h, double x, double y, double theta, double t, double X0);
    double Deta2 (double H, double h, double x, double y, double theta, double t, double X0);
    double Deta3 (double H, double h, double x, double y, double theta, double t, double X0);
    double U (double H, double h, double x, double y, double theta, double t, double X0, double z);
    double W (double H, double h, double x, double y, double theta, double t, double X0, double z);
}

#endif
