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

#include "waveFun.H"
#include <math.h>
#include <limits>

namespace otherFun
{
    #define PII 3.1415926535897932384626433832795028
    #define grav 9.81

    double interpolation (double x1, double x2, double y1, double y2, double xInt)
    {

        double yInt = y1 + (y2-y1)/(x2-x1)*(xInt-x1);

        return yInt;
    }
}

namespace StokesIFun
{
    #define PII 3.1415926535897932384626433832795028
    #define G 9.81

    double waveLength (double h, double T)
    {
        double L0 = G*T*T/(2.0*PII);
        double L = L0;

        for(int i=1; i<=100; i++)
        {
            L = L0*tanh(2.0*PII*h/L);
        }

        return L;
    }
    
    double eta (double H, double Kx, double x, double Ky, double y, double omega, double t, double phase)
    {
        double faseTot = Kx*x + Ky*y - omega*t + phase;
        
        double sup = H*0.5*cos(faseTot);
        
        return sup;
    }

    double U (double H, double h, double Kx, double x, double Ky, double y, double omega, double t, double phase, double z)
    {
        double k = sqrt(Kx*Kx + Ky*Ky);
        double faseTot = Kx*x + Ky*y - omega*t + phase;

        double velocity = H*0.5*omega*cos(faseTot)*cosh(k*z)/sinh(k*h);

        return velocity;
    }
    
    double W (double H, double h, double Kx, double x, double Ky, double y, double omega, double t, double phase, double z)
    {
        double k = sqrt(Kx*Kx + Ky*Ky);
        double faseTot = Kx*x + Ky*y - omega*t + phase;

        double velocity = H*0.5*omega*sin(faseTot)*sinh(k*z)/sinh(k*h);

        return velocity;
    }   
}

namespace StokesIIFun
{
    #define PII 3.1415926535897932384626433832795028
    #define G 9.81 
    
    double eta (double H, double h, double Kx, double x, double Ky, double y, double omega, double t, double phase)
    {
        double k = sqrt(Kx*Kx + Ky*Ky);
        double sigma = tanh(k*h);
        double faseTot = Kx*x + Ky*y - omega*t + phase;
        
        double sup = H*0.5*cos(faseTot) + k*H*H/4.0*(3.0-sigma*sigma)/(4.0*pow(sigma,3))*cos(2.0*faseTot);
        
        return sup;
    }
    
    double U (double H, double h, double Kx, double x, double Ky, double y, double omega, double t, double phase, double z)
    {
        double k = sqrt(Kx*Kx + Ky*Ky);
        double faseTot = Kx*x + Ky*y - omega*t + phase;

        double velocity = H*0.5*omega*cos(faseTot)*cosh(k*z)/sinh(k*h) + 3.0/4.0*H*H/4.0*omega*k*cosh(2.0*k*z)/pow(sinh(k*h),4)*cos(2.0*faseTot);

        return velocity;
    }
    
    double W (double H, double h, double Kx, double x, double Ky, double y, double omega, double t, double phase, double z)
    {
        double k = sqrt(Kx*Kx + Ky*Ky);
        double faseTot = Kx*x + Ky*y - omega*t + phase;

        double velocity = H*0.5*omega*sin(faseTot)*sinh(k*z)/sinh(k*h) + 3.0/4.0*H*H/4.0*omega*k*sinh(2.0*k*z)/pow(sinh(k*h),4)*sin(2.0*faseTot);

        return velocity;
    }   
    
    double timeLag (double H, double h, double Kx, double x, double Ky, double y, double T, double phase)
    {
        double omega = 2.0*PII/T;
        
        double startTime = 0;
        double lagIncrement = 0.001;

        double lag = 0;
        double lag0 = 0; // phase = 0
        double lag1 = 0; // phase = PI/2
        double lag2 = 0; // phase = PI
        double lag3 = 0; // phase = 3*PI/2
        double lag4 = 0; // phase = 2*PI

        double oldEta = 0;
        double newEta = 0;
        
        while (phase >= 2.0*PII)
        {
            phase -= 2.0*PII;
        }

        // Calculate time lag to start at phase = 0 (max)
        startTime = -T/4.0;
        while (true)
        {
            oldEta = eta(H, h, Kx, x, Ky, y, omega, startTime, phase);

            startTime += lagIncrement;

            newEta = eta(H, h, Kx, x, Ky, y, omega, startTime, phase);

            if ( oldEta > newEta )
            {
                lag0 = startTime;
                break;
            }
        }

        // Calculate time lag to start at phase = PI/2
        startTime = lag0;
        while (true)
        {
            oldEta = eta(H, h, Kx, x, Ky, y, omega, startTime, phase);

            startTime += lagIncrement;

            newEta = eta(H, h, Kx, x, Ky, y, omega, startTime, phase);

            if ( oldEta == 0.0 && newEta < oldEta )
            {
                lag1 = startTime - lagIncrement;
                break;
            }
            else if ( oldEta > 0.0 && newEta < 0.0 )
            {
                lag1 = startTime;
                break;
            }
        }

        // Calculate time lag to start at phase = PI
        startTime = lag1;
        while (true)
        {
            oldEta = eta(H, h, Kx, x, Ky, y, omega, startTime, phase);

            startTime += lagIncrement;

            newEta = eta(H, h, Kx, x, Ky, y, omega, startTime, phase);

            if ( oldEta < newEta )
            {
                lag2 = startTime;
                break;
            }
        }

        // Calculate time lag to start at phase = 3*PI/2
        startTime = lag2;
        while (true)
        {
            oldEta = eta(H, h, Kx, x, Ky, y, omega, startTime, phase);

            startTime += lagIncrement;

            newEta = eta(H, h, Kx, x, Ky, y, omega, startTime, phase);

            if ( oldEta == 0.0 && newEta > oldEta )
            {
                lag3 = startTime - lagIncrement;
                break;
            }
            else if ( oldEta < 0.0 && newEta > 0.0 )
            {
                lag3 = startTime;
                break;
            }
        }

        // Calculate time lag to start at phase = 2*PI
        startTime = lag3;
        while (true)
        {
            oldEta = eta(H, h, Kx, x, Ky, y, omega, startTime, phase);

            startTime += lagIncrement;

            newEta = eta(H, h, Kx, x, Ky, y, omega, startTime, phase);

            if ( oldEta > newEta )
            {
                lag4 = startTime;
                break;
            }
        }

        if ( phase >= 0.0 && phase < PII/2.0 )
        {
            lag = otherFun::interpolation(0.0, PII/2.0, lag0, lag1, phase);
        }
        else if ( phase >= PII/2.0 && phase < PII )
        {
            lag = otherFun::interpolation(PII/2.0, PII, lag1, lag2, phase);
        }
        else if ( phase >= PII && phase < 3.0*PII/2.0 )
        {
            lag = otherFun::interpolation(PII, 3.0*PII/2.0, lag2, lag3, phase);
        }
        else if ( phase >= 3.0*PII/2.0 && phase < 2.0*PII )
        {
            lag = otherFun::interpolation(3.0*PII/2.0, 2.0*PII, lag3, lag4, phase);
        }

        return lag;
    }
}

namespace Elliptic
{
    #define PII 3.1415926535897932384626433832795028
    #define ITER 50
    #define TOL std::numeric_limits<long double>::digits10

    void ellipticIntegralsKE (double m, double* K, double* E)
    {
        long double a,aOld,g,gOld,aux,sum;

        if ( m == 0.0 ) 
        {
            *K = PII/2.;
            *E = PII/2.;
            return;
        }

        a = 1.0L;
        g = sqrt(1.0L - m);
        aux = 1.0L;
        sum = 2.0L - m;

        for (int n = 0; n < ITER; n++)
        {
            gOld = g;
            aOld = a;
            a = 0.5L * (gOld + aOld);
            g = gOld * aOld;
            aux += aux;
            sum -= aux * (a * a - g);

            if ( fabsl(aOld - gOld) <= (aOld * pow(10.,-TOL)) )
            { 
                break;
            }

            g = sqrt(g);
        }

        *K = (double) (PII/2. / a);
        *E = (double) (PII/4. / a * sum);
        return;
    }

    double JacobiAmp (double u, double m)
    {
        long double a[ITER+1], g[ITER+1], c[ITER+1];
        long double aux, amp;
        int n;

        m = fabsl(m);

        if ( m == 0.0 ) 
        {
            return u;
        }
        
        if ( m == 1. ) 
        {
            return 2. * atan( exp(u) ) - PII/2.;
        }

        a[0] = 1.0L;
        g[0] = sqrtl(1.0L - m);
        c[0] = sqrtl(m);

        aux = 1.0L;

        for (n = 0; n < ITER; n++)
        {
            if ( fabsl(a[n] - g[n]) < (a[n] * pow(10.,-TOL)) ) 
            {
                break;
            }

            aux += aux;
            a[n+1] = 0.5L * (a[n] + g[n]);
            g[n+1] = sqrtl(a[n] * g[n]);
            c[n+1] = 0.5L * (a[n] - g[n]);
        }

        amp = aux * a[n] * u;

        for (; n > 0; n--) 
        {
            amp = 0.5L * ( amp + asinl( c[n] * sinl(amp) / a[n]) );
        }

        return (double) amp; 
    }

    void JacobiSnCnDn (double u, double m, double* Sn, double* Cn, double* Dn)
    {
        double amp = Elliptic::JacobiAmp( u, m );

        *Sn = sin( amp );
        *Cn = cos( amp );
        *Dn = sqrt(1.0 - m * sin( amp ) * sin( amp ));
        return;
    }

}

namespace cnoidalFun
{
    #define PII 3.1415926535897932384626433832795028
    #define G 9.81

    double eta (double H, double m, double kx, double ky, double T, double x, double y, double t)
    {
        double K, E;
        Elliptic::ellipticIntegralsKE(m, &K, &E);

        double uCnoidal = K/PII*(kx*x + ky*y - 2.0*PII*t/T);

        double sn, cn, dn;
        Elliptic::JacobiSnCnDn(uCnoidal, m, &sn, &cn, &dn);

        double etaCnoidal = H*((1.0-E/K)/m - 1.0 + pow(cn,2));

        return etaCnoidal;
    }

    double etaCnoidal1D (double H, double m, double t, double T)
    {
        double K, E;
        Elliptic::ellipticIntegralsKE(m, &K, &E);

        double uCnoidal = -2.0*K*(t/T); // for cn Jacobian elliptic function

        double sn, cn, dn;
        Elliptic::JacobiSnCnDn(uCnoidal, m, &sn, &cn, &dn);

        double etaCnoidal = H*((1.0-E/K)/m - 1.0 + pow(cn,2));

        return etaCnoidal;
    }

    double etaMeanSq (double H, double m, double T)
    {
        double eta = 0;
        double etaSumSq = 0;

        for(int i=0; i<1000; i++)
        {
            eta = etaCnoidal1D(H, m, i*T/(1000.0), T);
            etaSumSq += eta*eta;
        }

        etaSumSq /= 1000.0;
        return etaSumSq;
    }

    double d1EtaDx (double H, double m, double uCnoidal, double L)
    {
        double K, E;
        Elliptic::ellipticIntegralsKE(m, &K, &E);

        double dudx = 0;
        dudx = 2.0*K/L;

        double sn, cn, dn;
        Elliptic::JacobiSnCnDn(uCnoidal, m, &sn, &cn, &dn);

        double deriv = -2.0*H*cn*dn*sn*dudx;

        return deriv;
    }

    double d2EtaDx (double H, double m, double uCnoidal, double L)
    {
        double K, E;
        Elliptic::ellipticIntegralsKE(m, &K, &E);

        double dudxx = 0;
        dudxx = 4.0*K*K/L/L;

        double sn, cn, dn;
        Elliptic::JacobiSnCnDn(uCnoidal, m, &sn, &cn, &dn);

        double deriv = 2.0*H*(dn*dn*sn*sn - cn*cn*dn*dn + m*cn*cn*sn*sn)*dudxx;

        return deriv;
    }

    double d3EtaDx (double H, double m, double uCnoidal, double L)
    {
        double K, E;
        Elliptic::ellipticIntegralsKE(m, &K, &E);

        double dudxxx = 0;
        dudxxx = 8.0*K*K*K/L/L/L;

        double sn, cn, dn;
        Elliptic::JacobiSnCnDn(uCnoidal, m, &sn, &cn, &dn);

        double deriv = 8.0*H*cn*dn*sn*(m*cn*cn + dn*dn - m*sn*sn)*dudxxx;

        return deriv;
    }

    double timeLag (double H, double m, double kx, double ky, double T, double x, double y, double phase)
    {
        double startTime = 0;
        double lagIncrement = 0.001;

        double lag = 0;
        double lag0 = 0; // phase = 0
        double lag1 = 0; // phase = PI/2
        double lag2 = 0; // phase = PI
        double lag3 = 0; // phase = 3*PI/2
        double lag4 = 0; // phase = 2*PI

        double oldEta = 0;
        double newEta = 0;
        
        while (phase >= 2.0*PII)
        {
            phase -= 2.0*PII;
        }

        // Calculate time lag to start at phase = 0 (max)
        startTime = -T/4.0;
        while (true)
        {
            oldEta = eta(H, m, kx, ky, T, x, y, startTime);

            startTime += lagIncrement;

            newEta = eta(H, m, kx, ky, T, x, y, startTime);

            if ( oldEta > newEta )
            {
                lag0 = startTime;
                break;
            }
        }

        // Calculate time lag to start at phase = PI/2
        startTime = lag0;
        while (true)
        {
            oldEta = eta(H, m, kx, ky, T, x, y, startTime);

            startTime += lagIncrement;

            newEta = eta(H, m, kx, ky, T, x, y, startTime);

            if ( oldEta == 0.0 && newEta < oldEta )
            {
                lag1 = startTime - lagIncrement;
                break;
            }
            else if ( oldEta > 0.0 && newEta < 0.0 )
            {
                lag1 = startTime;
                break;
            }
        }

        // Calculate time lag to start at phase = PI
        startTime = lag1;
        while (true)
        {
            oldEta = eta(H, m, kx, ky, T, x, y, startTime);

            startTime += lagIncrement;

            newEta = eta(H, m, kx, ky, T, x, y, startTime);

            if ( oldEta < newEta )
            {
                lag2 = startTime;
                break;
            }
        }

        // Calculate time lag to start at phase = 3*PI/2
        startTime = lag2;
        while (true)
        {
            oldEta = eta(H, m, kx, ky, T, x, y, startTime);

            startTime += lagIncrement;

            newEta = eta(H, m, kx, ky, T, x, y, startTime);

            if ( oldEta == 0.0 && newEta > oldEta )
            {
                lag3 = startTime - lagIncrement;
                break;
            }
            else if ( oldEta < 0.0 && newEta > 0.0 )
            {
                lag3 = startTime;
                break;
            }
        }

        // Calculate time lag to start at phase = 2*PI
        startTime = lag3;
        while (true)
        {
            oldEta = eta(H, m, kx, ky, T, x, y, startTime);

            startTime += lagIncrement;

            newEta = eta(H, m, kx, ky, T, x, y, startTime);

            if ( oldEta > newEta )
            {
                lag4 = startTime;
                break;
            }
        }

        if ( phase >= 0.0 && phase < PII/2.0 )
        {
            lag = otherFun::interpolation(0.0, PII/2.0, lag0, lag1, phase);
        }
        else if ( phase >= PII/2.0 && phase < PII )
        {
            lag = otherFun::interpolation(PII/2.0, PII, lag1, lag2, phase);
        }
        else if ( phase >= PII && phase < 3.0*PII/2.0 )
        {
            lag = otherFun::interpolation(PII, 3.0*PII/2.0, lag2, lag3, phase);
        }
        else if ( phase >= 3.0*PII/2.0 && phase < 2.0*PII )
        {
            lag = otherFun::interpolation(3.0*PII/2.0, 2.0*PII, lag3, lag4, phase);
        }

        return lag;
    }

    int calculations (double H, double d, double T, double* mOut, double* LOut)
    {
        double mTolerance = 0.0001;
        double mElliptic = 0.5;
        double LElliptic = 0;
        double KElliptic = 0;
        double EElliptic = 0;
        double phaseSpeed = 0;

        double mError = 0.0;
        double mMinError = 999;

        while (mElliptic < 1.0)
        {
            // Calculate K & E
            double KElliptic, EElliptic;
            Elliptic::ellipticIntegralsKE(mElliptic, &KElliptic, &EElliptic);

            LElliptic = KElliptic*sqrt(16.0*pow(d,3)*mElliptic/(3.0*H));        // (Svendsen p404)

            phaseSpeed = sqrt(G*d)*(1.0 - H/d/2.0 + H/d/mElliptic*(1.0-3.0/2.0*EElliptic/KElliptic));        // (Svendsen p405)

            mError = fabs(T-LElliptic/phaseSpeed);

            if (mError <= mMinError)
            {
                *mOut = mElliptic;
                *LOut = LElliptic;
                mMinError = mError;
            }

            mElliptic += mTolerance;
        }

        return 0;
    }

    double U (double H, double h, double m, double kx, double ky, double T, double x, double y, double t, double z)
    {
        double K, E;
        Elliptic::ellipticIntegralsKE(m, &K, &E);

        double uCnoidal = K/PII*(kx*x + ky*y - 2.0*PII*t/T);
        double k = sqrt(kx*kx + ky*ky);
        double L = 2.0*PII/k;
        double c = L/T;

        double etaCN = eta(H, m, kx, ky, T, x, y, t);
        double etaXX = d2EtaDx(H, m, uCnoidal, L);
        double etaMS = etaMeanSq(H, m, T);

        double velocity = c*etaCN/h - c*(etaCN*etaCN/h/h + etaMS/h/h) + 1.0/2.0*c*h*(1.0/3.0 - z*z/h/h)*etaXX;

        return velocity;
    }

    double W (double H, double h, double m, double kx, double ky, double T, double x, double y, double t, double z)
    {
        double K, E;
        Elliptic::ellipticIntegralsKE(m, &K, &E);

        double uCnoidal = K/PII*(kx*x + ky*y - 2.0*PII*t/T);
        double k = sqrt(kx*kx + ky*ky);
        double L = 2.0*PII/k;
        double c = L/T;

        double etaCN = eta(H, m, kx, ky, T, x, y, t);
        double etaX = d1EtaDx(H, m, uCnoidal, L);
        double etaXXX = d3EtaDx(H, m, uCnoidal, L);

        double velocity = -c*z*( etaX/h*(1.0-2.0*etaCN/h) + 1.0/6.0*h*(1.0-z*z/h/h)*etaXXX ) ;

        return velocity;
    }
}

namespace stokesVFun
{
    // Skjelbreia 1960
    
    #define PII 3.1415926535897932384626433832795028
    #define G 9.81

    double A11 (double h, double k)
    {
        double s = sinh(k*h);

        double A = 1.0/s;

        return A;
    }

    double A13 (double h, double k)
    {
        double s = sinh(k*h);
        double c = cosh(k*h);

        double A = 
            -pow(c,2)*(5.0*pow(c,2)+1.0)
                /
            (8.0*pow(s,5));

        return A;
    }
 
    double A15 (double h, double k)
    {
        double s = sinh(k*h);
        double c = cosh(k*h);

        double A = 
            -(1184.0*pow(c,10)-1440.0*pow(c,8)-1992.0*pow(c,6)+2641.0*pow(c,4)-249.0*pow(c,2)+18)
                /
            (1536.0*pow(s,11));

        return A;
    }   

    double A22 (double h, double k)
    {
        double s = sinh(k*h);

        double A = 
            3.0
                /
            (8.0*pow(s,4));

        return A;
    }

    double A24 (double h, double k)
    {
        double s = sinh(k*h);
        double c = cosh(k*h);

        double A = 
            (192.0*pow(c,8)-424.0*pow(c,6)-312.0*pow(c,4)+480.0*pow(c,2)-17)
                /
            (768.0*pow(s,10));

        return A;
    }

    double A33 (double h, double k)
    {
        double s = sinh(k*h);
        double c = cosh(k*h);

        double A = 
            (13.0-4.0*pow(c,2))
                /
            (64.0*pow(s,7));

        return A;
    }

    double A35 (double h, double k)
    {
        double s = sinh(k*h);
        double c = cosh(k*h);

        double A = 
            (512.0*pow(c,12)+4224.0*pow(c,10)-6800.0*pow(c,8)-12808.0*pow(c,6)+16704.0*pow(c,4)-3154.0*pow(c,2)+107.0)
                /
            (4096.0*pow(s,13)*(6.0*pow(c,2)-1.0));

        return A;
    }

    double A44 (double h, double k)
    {
        double s = sinh(k*h);
        double c = cosh(k*h);

        double A = 
            (80.0*pow(c,6)-816.0*pow(c,4)+1338.0*pow(c,2)-197.0)
                /
            (1536.0*pow(s,10)*(6.0*pow(c,2)-1.0));

        return A;
    }

    double A55 (double h, double k)
    {
        double s = sinh(k*h);
        double c = cosh(k*h);

        double A = 
            -(2880.0*pow(c,10)-72480.0*pow(c,8)+324000.0*pow(c,6)-432000.0*pow(c,4)+163470.0*pow(c,2)-16245.0)
                /
            (61440.0*pow(s,11)*(6.0*pow(c,2)-1.0)*(8.0*pow(c,4)-11.0*pow(c,2)+3.0));

        return A;
    }

    double B22 (double h, double k)
    {
        double s = sinh(k*h);
        double c = cosh(k*h);

        double B = 
            (2.0*pow(c,2)+1)*c
                /
            (4.0*pow(s,3));

        return B;
    }

    double B24 (double h, double k)
    {
        double s = sinh(k*h);
        double c = cosh(k*h);

        double B = 
            (272.0*pow(c,8)-504.0*pow(c,6)-192.0*pow(c,4)+322.0*pow(c,2)+21.0)*c
                /
            (384.0*pow(s,9));

        return B;
    }

    double B33 (double h, double k)
    {
        double s = sinh(k*h);
        double c = cosh(k*h);

        double B = 
            (8.0*pow(c,6)+1.0)*3.0
                /
            (64.0*pow(s,6));

        return B;
    }

    double B33k (double h, double k) // d B33 / d k
    {
        double s = sinh(k*h);
        double c = cosh(k*h);

        double sk = h*s;
        double ck = h*c;

        double B = 
            9.0*pow(c,5)*ck/(4.0*pow(s,6))- 
            (9.0*(8.0*pow(c,6)+1.0))/(32.0*pow(s,7))*sk;

        return B;
    }

    double B35 (double h, double k)
    {
        double s = sinh(k*h);
        double c = cosh(k*h);

        double B = 
            (88128.0*pow(c,14)-208224.0*pow(c,12)+70848.0*pow(c,10)+54000.0*pow(c,8)-21816.0*pow(c,6)+6264.0*pow(c,4)-54.0*pow(c,2)-81)
                /
            (12288.0*pow(s,12)*(6.0*pow(c,2)-1.0));

        return B;
    }

    double B35k (double h, double k) // d B35 / d k
    {
        double s = sinh(k*h);
        double c = cosh(k*h);

        double sk = h*s;
        double ck = h*c;

        double B = 
            (14.0*88128.0*pow(c,13)*ck-12.0*208224.0*pow(c,11)*ck
            +10.0*70848.0*pow(c,9)*ck+8.0*54000.0*pow(c,7)*ck
            -6.0*21816.0*pow(c,5)*ck+4.0*6264.0*pow(c,3)*ck
            -2.0*54.0*c*ck)
            /(12288.0*pow(s,12)*(6.0*pow(c,2)-1.0))
            -(88128.0*pow(c,14)-208224.0*pow(c,12)+70848.0*pow(c,10)+54000.0*pow(c,8) 
            -21816.0*pow(c,6)+6264.0*pow(c,4)-54.0*pow(c,2)-81.0)*12.0 
            /(12288.0*pow(s,13)*(6.0*pow(c,2)-1.0))*sk 
            -(88128.0*pow(c,14)-208224.0*pow(c,12)+70848.0*pow(c,10)+54000.0*pow(c,8) 
            -21816.0*pow(c,6)+6264.0*pow(c,4)-54.0*pow(c,2)-81.0)*12.0*c*ck 
            /(12288.0*pow(s,12)*pow((6.0*pow(c,2)-1.0),2));

        return B;
    }

    double B44 (double h, double k)
    {
        double s = sinh(k*h);
        double c = cosh(k*h);

        double B = 
            (768.0*pow(c,10)-448.0*pow(c,8)-48.0*pow(c,6)+48.0*pow(c,4)+106.0*pow(c,2)-21.0)*c
                /
            (384.0*pow(s,9)*(6.0*pow(c,2)-1.0));

        return B;
    }

    double B55 (double h, double k)
    {
        double s = sinh(k*h);
        double c = cosh(k*h);

        double B = 
            (192000.0*pow(c,16)-262720.0*pow(c,14)+83680.0*pow(c,12)+20160.0*pow(c,10)-7280.0*pow(c,8)+7160.0*pow(c,6)-1800.0*pow(c,4)-1050.0*pow(c,2)+225.0)
                /
            (12288.0*pow(s,10)*(6.0*pow(c,2)-1.0)*(8.0*pow(c,4)-11.0*pow(c,2)+3.0));

        return B;
    }

    double B55k (double h, double k) // d B55 / d k
    {
        double s = sinh(k*h);
        double c = cosh(k*h);

        double sk = h*s;
        double ck = h*c;

        double B = 
            (16.0*192000.0*pow(c,15)*ck-14.0*262720.0*pow(c,13)*ck
            +12.0*83680.0*pow(c,11)*ck+10.0*20160.0*pow(c,9)*ck
            -8.0*7280.0*pow(c,7)*ck+6.0*7160.0*pow(c,5)*ck
            -4.0*1800.0*pow(c,3)*ck-2.0*1050.0*pow(c,1)*ck)
            /(12288.0*pow(s,10)*(6.0*pow(c,2)-1.0)*(8.0*pow(c,4)-11.0*pow(c,2)+3.0))
            -(192000.0*pow(c,16)-262720.0*pow(c,14)+83680.0*pow(c,12)+20160.0*pow(c,10) 
            -7280.0*pow(c,8)+7160.0*pow(c,6)-1800.0*pow(c,4)-1050.0*pow(c,2)+225.0)*10.0
            /(12288.0*pow(s,11)*(6.0*pow(c,2)-1.0)*(8.0*pow(c,4)-11.0*pow(c,2)+3.0))*sk
            -(192000.0*pow(c,16)-262720.0*pow(c,14)+83680.0*pow(c,12)+20160.0*pow(c,10) 
            -7280.0*pow(c,8)+7160.0*pow(c,6)-1800.0*pow(c,4)-1050.0*pow(c,2)+225.0)
            *12.0*c*ck 
            /(12288.0*pow(s,10)*pow((6.0*pow(c,2)-1.0),2)*(8.0*pow(c,4)-11.0*pow(c,2)+3.0))
            -(192000.0*pow(c,16)-262720.0*pow(c,14)+83680.0*pow(c,12)+20160.0*pow(c,10) 
            -7280.0*pow(c,8)+7160.0*pow(c,6)-1800.0*pow(c,4)-1050.0*pow(c,2)+225.0)
            *(32.0*pow(c,3)-22.0*c)*ck
            /(12288.0*pow(s,10)*(6.0*pow(c,2)-1.0)*pow((8.0*pow(c,4)-11.0*pow(c,2)+3.0),2));

        return B;
    }

    double C1 (double h, double k)
    {
        double s = sinh(k*h);
        double c = cosh(k*h);

        double C = 
            (8.0*pow(c,4)-8.0*pow(c,2)+9.0)
                /
            (8.0*pow(s,4));

        return C;
    }

    double C1k (double h, double k)
    {
        double s = sinh(k*h);
        double c = cosh(k*h);

        double sk = h*s;
        double ck = h*c;

        double C = 
            (4.0*8.0*pow(c,3)*ck-2.0*8.0*c*ck)/(8.0*pow(s,4))
            -(8.0*pow(c,4)-8.0*pow(c,2)+9.0)*4.0*sk/(8.0*pow(s,5));

        return C;
    }

    double C2 (double h, double k)
    {
        double s = sinh(k*h);
        double c = cosh(k*h);

        double C = 
            (3840.0*pow(c,12)-4096.0*pow(c,10) + 2592.0*pow(c,8)-1008.0*pow(c,6)+5944.0*pow(c,4)-1830.0*pow(c,2)+147.0) // - 2592
                /
            (512.0*pow(s,10)*(6.0*pow(c,2)-1.0));

        return C;
    }

    double C2k (double h, double k)
    {
        double s = sinh(k*h);
        double c = cosh(k*h);

        double sk = h*s;
        double ck = h*c;

        double C = 
            (12.0*3840.0*pow(c,11)*ck-10.0*4096.0*pow(c,9)*ck
            +8.0*2592.0*pow(c,7)*ck-6.0*1008.0*pow(c,5)*ck+
            4.0*5944.0*pow(c,3)*ck-2.0*1830.0*c*ck)
            /(512.0*pow(s,10)*(6.0*pow(c,2)-1.0))
            -(3840.0*pow(c,12)-4096.0*pow(c,10)+2592.0*pow(c,8)-1008.0*pow(c,6)+
            5944.0*pow(c,4)-1830.0*pow(c,2)+147.0)*10.0*sk/
            (512.0*pow(s,11)*(6.0*pow(c,2)-1.0))
            -(3840.0*pow(c,12)-4096.0*pow(c,10)+2592.0*pow(c,8)-1008.0*pow(c,6)+
            5944.0*pow(c,4)-1830.0*pow(c,2)+147.0)*12.0*c*ck
            /(512.0*pow(s,10)*pow((6.0*pow(c,2)-1.0),2));

        return C;
    }

    double C3 (double h, double k)
    {
        double s = sinh(k*h);
        double c = cosh(k*h);

        double C = 
            (-1.0)
                /
            (4.0*s*c);

        return C;
    }

    double C4 (double h, double k)
    {
        double s = sinh(k*h);
        double c = cosh(k*h);

        double C = 
            (12.0*pow(c,8)+36.0*pow(c,6)-162.0*pow(c,4)+141.0*pow(c,2)-27.0)
                /
            (192.0*c*pow(s,9));

        return C;
    }

    int StokesVNR (double H, double d, double T, double* kOut, double* LambdaOut, double* f1Out, double* f2Out )
    {
        // Solves L and Lambda using Newton-Raphson

        // Two nonlinear algebriac equations (nonlinear wave
        // Lambdalitude equation (f1) and dispersion relationship (f2)).

        // f1=PII*H/d-2*PII/(k*d)*(Lambda+Lambda**3*b33+Lambda**5*(b35+b55))
        // f2=(2*PII*d)/(-gy*xxt**2)
        //	-k*d/(2*PII)*dtanh(k*d)*(1+Lambda**2*c1+Lambda**4*c2)

        double f1 = 1;
        double f2 = 1;

        double k = 2.0*PII/(sqrt(G*d)*T);
        double Lambda = H/2.0*k;

        double Bmat11, Bmat12, Bmat21, Bmat22;
        double b33, b35, b55, b33k, b35k, b55k;
        double c1, c2, c1k, c2k;
        double kPr, LambdaPr;

        int n = 0;

        //while ( (abs(f1)>1.0e-12 || abs(f2)>1.0e-12) && n<1000 )
        while ( (fabs(f1)>1.0e-12 || fabs(f2)>1.0e-12) && n<10000 )
        {
            b33 = B33(d, k);
            b35 = B35(d, k);
            b55 = B55(d, k);
            c1 = C1(d, k);
            c2 = C2(d, k);

            // Derivatives
            b33k = B33k(d, k);
            b35k = B35k(d, k);
            b55k = B55k(d, k);
            c1k = C1k(d, k);
            c2k = C2k(d, k);

            // d f1 / d k
            Bmat11 = 2.0*PII/(pow(k,2)*d)*(Lambda+pow(Lambda,3)*b33+pow(Lambda,5)*(b35+b55))
                    -2.0*PII/(k*d)*(pow(Lambda,3)*b33k+pow(Lambda,5)*(b35k+b55k));
            // d f1 / d Lambda
            Bmat12 = -2.0*PII/(k*d)*
                     (1.0+3.0*pow(Lambda,2)*b33+5.0*pow(Lambda,4)*(b35+b55));
            // d f2 / d k
            Bmat21 = -d/(2.0*PII)*tanh(k*d)*(1.0+pow(Lambda,2)*c1+pow(Lambda,4)*c2)
                    -k*d/(2.0*PII)*(1.0-pow((tanh(k*d)),2))*d*
                    (1.0+pow(Lambda,2)*c1+pow(Lambda,4)*c2)
                    -k*d/(2.0*PII)*tanh(k*d)
                    *(pow(Lambda,2)*c1k+pow(Lambda,4)*c2k);
            // d f2 / d Lambda
            Bmat22 = -k*d/(2.0*PII)*tanh(k*d)
                     *(2.0*Lambda*c1+4.0*pow(Lambda,3)*c2);

            // Calculate f1 and f2
            f1 = PII*H/d-2.0*PII/(k*d)*(Lambda+pow(Lambda,3)*b33+pow(Lambda,5)*(b35+b55));

            f2 = (2.0*PII*d)/(G*pow(T,2))-k*d/(2.0*PII)*tanh(k*d)
                *(1.0+pow(Lambda,2)*c1+pow(Lambda,4)*c2);
            // Solve the two-equation system Bx=-F
            // b(1,1)*L + b(1,2)*La = -f1
            // b(2,1)*L + b(2,2)*La = -f2

            LambdaPr = (f1*Bmat21-f2*Bmat11)/(Bmat11*Bmat22-Bmat12*Bmat21);
            kPr = (f2*Bmat12-f1*Bmat22)/(Bmat11*Bmat22-Bmat12*Bmat21);

            Lambda += LambdaPr;
            k += kPr;

            n++;
        }

        *kOut = k;
        *LambdaOut = Lambda;

        *f1Out = fabs(f1);
        *f2Out = fabs(f2);

        return 1;
    }

    int StokesExtendedSolver (double H, double d, double T, double* kOut, double* LambdaOut, double* LambdaErrOut )
    {
        double b33, b35, b55;
        double c1, c2;
        double Aeq, Beq, Ceq;

        double L0 = G*T*T/(2.0*PII); // L deep water
        double Lred = L0*PII/10.0; // L limit shallow water 

        double k = 0.0;
        double Lambda = 0.0;

        double Ltest = 0.0;
        double ktest = 0.0;
        double LambdaTest = 0.0;
        double LambdaCalc = 0.0;
        double LambdaErr = 0.0;
        double LambdaMAXErr = 999999.0;

        for(int i = 0; i<=1000; i++)
        {
            Ltest = L0 - i*(L0-Lred)/1000;
            ktest = 2.0*PII/Ltest;

            b33 = B33(d, ktest);
            b35 = B35(d, ktest);
            b55 = B55(d, ktest);
            c1 = C1(d, ktest);
            c2 = C2(d, ktest);
            
            Aeq = c2;
            Beq = c1;
            Ceq = 1.0-Ltest/(L0*tanh(ktest*d));

            // Check if solvable: X
            if ( (Beq*Beq-4.0*Aeq*Ceq) >= 0.0 )
            {
                // Check if solvable: Lambda ( +sqrt, since Beq > 0 always )
                if ( (-Beq + sqrt(Beq*Beq-4.0*Aeq*Ceq))/(2.0*Aeq) >= 0.0 )
                {
                    LambdaTest = sqrt( (-Beq + sqrt(Beq*Beq-4.0*Aeq*Ceq))/(2.0*Aeq) );

                    // Second equation
                    LambdaCalc = PII*H/(Ltest*(1 + b33*Lambda*Lambda + (b35+b55)*pow(Lambda, 4) ));

                    // Calculate difference (error)
                    LambdaErr = fabs(LambdaTest-LambdaCalc);

                    if ( LambdaErr <= LambdaMAXErr && LambdaTest != 0 )
                    {
                        LambdaMAXErr = LambdaErr;
                        Lambda = LambdaTest;
                        k = ktest;
                    }
                }
            }
        }

        *kOut = k;
        *LambdaOut = Lambda;
        *LambdaErrOut = LambdaMAXErr;

        return 1;
    }

    double eta (double h, double kx, double ky, double lambda, double T, double x, double y, double t, double phase)
    {
        // Eta Stokes V
        double k = sqrt(kx*kx + ky*ky);

        double b22 = B22(h, k);
        double b24 = B24(h, k);
        double b33 = B33(h, k);
        double b35 = B35(h, k);
        double b44 = B44(h, k);
        double b55 = B55(h, k);

        double amp1 = lambda/k;
        double amp2 = (b22*pow(lambda,2)+b24*pow(lambda,4))/k;
        double amp3 = (b33*pow(lambda,3)+b35*pow(lambda,5))/k;
        double amp4 = b44*pow(lambda,4)/k;
        double amp5 = b55*pow(lambda,5)/k;

        double theta = kx*x + ky*y - 2.0*PII/T*t + phase;

        double C = 
            amp1*cos(theta) + amp2*cos(2*theta)
            + amp3*cos(3*theta) + amp4*cos(4*theta)
            + amp5*cos(5*theta);

        return C;
    }
    
    double timeLag (double h, double kx, double ky, double lambda, double T, double x, double y, double phase)
    {
        double startTime = 0;
        double lagIncrement = 0.001;

        double lag = 0;
        double lag0 = 0; // phase = 0
        double lag1 = 0; // phase = PI/2
        double lag2 = 0; // phase = PI
        double lag3 = 0; // phase = 3*PI/2
        double lag4 = 0; // phase = 2*PI

        double oldEta = 0;
        double newEta = 0;
        
        while (phase >= 2.0*PII)
        {
            phase -= 2.0*PII;
        }

        // Calculate time lag to start at phase = 0 (max)
        startTime = -T/4.0;
        while (true)
        {
            oldEta = eta(h, kx, ky, lambda, T, x, y, startTime, phase);

            startTime += lagIncrement;

            newEta = eta(h, kx, ky, lambda, T, x, y, startTime, phase);

            if ( oldEta > newEta )
            {
                lag0 = startTime;
                break;
            }
        }

        // Calculate time lag to start at phase = PI/2
        startTime = lag0;
        while (true)
        {
            oldEta = eta(h, kx, ky, lambda, T, x, y, startTime, phase);

            startTime += lagIncrement;

            newEta = eta(h, kx, ky, lambda, T, x, y, startTime, phase);

            if ( oldEta == 0.0 && newEta < oldEta )
            {
                lag1 = startTime - lagIncrement;
                break;
            }
            else if ( oldEta > 0.0 && newEta < 0.0 )
            {
                lag1 = startTime;
                break;
            }
        }

        // Calculate time lag to start at phase = PI
        startTime = lag1;
        while (true)
        {
            oldEta = eta(h, kx, ky, lambda, T, x, y, startTime, phase);

            startTime += lagIncrement;

            newEta = eta(h, kx, ky, lambda, T, x, y, startTime, phase);

            if ( oldEta < newEta )
            {
                lag2 = startTime;
                break;
            }
        }

        // Calculate time lag to start at phase = 3*PI/2
        startTime = lag2;
        while (true)
        {
            oldEta = eta(h, kx, ky, lambda, T, x, y, startTime, phase);

            startTime += lagIncrement;

            newEta = eta(h, kx, ky, lambda, T, x, y, startTime, phase);

            if ( oldEta == 0.0 && newEta > oldEta )
            {
                lag3 = startTime - lagIncrement;
                break;
            }
            else if ( oldEta < 0.0 && newEta > 0.0 )
            {
                lag3 = startTime;
                break;
            }
        }

        // Calculate time lag to start at phase = 2*PI
        startTime = lag3;
        while (true)
        {
            oldEta = eta(h, kx, ky, lambda, T, x, y, startTime, phase);

            startTime += lagIncrement;

            newEta = eta(h, kx, ky, lambda, T, x, y, startTime, phase);

            if ( oldEta > newEta )
            {
                lag4 = startTime;
                break;
            }
        }

        if ( phase >= 0.0 && phase < PII/2.0 )
        {
            lag = otherFun::interpolation(0.0, PII/2.0, lag0, lag1, phase);
        }
        else if ( phase >= PII/2.0 && phase < PII )
        {
            lag = otherFun::interpolation(PII/2.0, PII, lag1, lag2, phase);
        }
        else if ( phase >= PII && phase < 3.0*PII/2.0 )
        {
            lag = otherFun::interpolation(PII, 3.0*PII/2.0, lag2, lag3, phase);
        }
        else if ( phase >= 3.0*PII/2.0 && phase < 2.0*PII )
        {
            lag = otherFun::interpolation(3.0*PII/2.0, 2.0*PII, lag3, lag4, phase);
        }

        return lag;
    }
    
    double phaseLag (double d, double k, double lambda, double T)
    {
        // Lag Stokes V

        double phase = 0.0;
        double phaseIncrement = 0.001;
        double oldEta = 0.0;
        double newEta = 0.0;
        
        // Calculate phase lag to start at waterDepth_ and descending
        while (true)
        {
            oldEta = eta(d, k, 0, lambda, T, 0.0, 0.0, 0.0, phase);

            phase += phaseIncrement;

            newEta = eta(d, k, 0, lambda, T, 0.0, 0.0, 0.0, phase);

            if ( oldEta == 0.0 && newEta>oldEta )
            {
                break;
            }
            else if ( oldEta<0.0 && newEta>0.0 )
            {
                break;
            }
        }
        return phase;
    }

    double U (double d, double kx, double ky, double lambda, double T, double x, double y, double t, double phase, double z)
    {
        // Horizontal velocity Stokes V
        double k = sqrt(kx*kx + ky*ky);

        double a11 = A11(d, k);
        double a13 = A13(d, k);
        double a15 = A15(d, k);
        double a22 = A22(d, k);
        double a24 = A24(d, k);
        double a33 = A33(d, k);
        double a35 = A35(d, k);
        double a44 = A44(d, k);
        double a55 = A55(d, k);

        double a1u=2.0*PII/T/k*(lambda*a11+pow(lambda,3)*a13+pow(lambda,5)*a15);
        double a2u=2.0*2.0*PII/T/k*(pow(lambda,2)*a22+pow(lambda,4)*a24);
        double a3u=3.0*2.0*PII/T/k*(pow(lambda,3)*a33+pow(lambda,5)*a35);
        double a4u=4.0*2.0*PII/T/k*(pow(lambda,4)*a44);
        double a5u=5.0*2.0*PII/T/k*(pow(lambda,5)*a55);

        double velU = 0;

        double theta = kx*x + ky*y - 2.0*PII/T*t + phase;

        velU = a1u*cosh(k*z)*cos(theta)
                 + a2u*cosh(2.0*k*z)*cos(2.0*(theta))
                 + a3u*cosh(3.0*k*z)*cos(3.0*(theta))
                 + a4u*cosh(4.0*k*z)*cos(4.0*(theta))
                 + a5u*cosh(5.0*k*z)*cos(5.0*(theta));

        return velU;
    }

    double V (double d, double kx, double ky, double lambda, double T, double x, double y, double t, double phase, double z)
    {
        // Horizontal velocity Stokes V
        double k = sqrt(kx*kx + ky*ky);

        double a11 = A11(d, k);
        double a13 = A13(d, k);
        double a15 = A15(d, k);
        double a22 = A22(d, k);
        double a24 = A24(d, k);
        double a33 = A33(d, k);
        double a35 = A35(d, k);
        double a44 = A44(d, k);
        double a55 = A55(d, k);

        double a1u=2.0*PII/T/k*(lambda*a11+pow(lambda,3)*a13+pow(lambda,5)*a15);
        double a2u=2.0*2.0*PII/T/k*(pow(lambda,2)*a22+pow(lambda,4)*a24);
        double a3u=3.0*2.0*PII/T/k*(pow(lambda,3)*a33+pow(lambda,5)*a35);
        double a4u=4.0*2.0*PII/T/k*(pow(lambda,4)*a44);
        double a5u=5.0*2.0*PII/T/k*(pow(lambda,5)*a55);

        double velV = 0;

        double theta = kx*x+ky*y-2.0*PII/T*t+phase;

        velV = a1u*sinh(k*z)*sin(theta)
                + a2u*sinh(2.0*k*z)*sin(2.0*(theta))
                + a3u*sinh(3.0*k*z)*sin(3.0*(theta))
                + a4u*sinh(4.0*k*z)*sin(4.0*(theta))
                + a5u*sinh(5.0*k*z)*sin(5.0*(theta));

        return velV;
    }
}

namespace stokesVFentonFun
{
    // Fenton 1985
    
    #define PII 3.1415926535897932384626433832795028
    #define G 9.81

    double C0 (double h, double k)
    {
        double C = sqrt( tanh(k*h) );

        return C;
    }

    double C2 (double h, double k)
    {
        double s = 1.0/cosh(2.0*k*h);

        double C = C0(h,k)*(2.0 + 7.0*s*s)/(4.0*(1.0-s)*(1.0-s));

        return C;
    }

    double C4 (double h, double k)
    {
        double s = 1.0/cosh(2.0*k*h);

        double C = 
            C0(h,k)*(4.0 + 32.0*s - 116.0*s*s - 400.0*pow(s,3) - 71.0*pow(s,4) + 146.0*pow(s,5))
            /(32.0*pow((1.0-s),5));

        return C;
    }

    double D2 (double h, double k)
    {
        double CTH = 1.0/tanh(k*h);

        double D = -sqrt(CTH)/2.0;

        return D;
    }

    double D4 (double h, double k)
    {
        double CTH = 1.0/tanh(k*h);
        double s = 1.0/cosh(2.0*k*h);

        double D = sqrt(CTH)*(2.0 + 4.0*s + s*s + 2.0*pow(s,3))/(8.0*pow((1.0-s),3));

        return D;
    }

    double error (double H, double h, double k, double T)
    {
        double c0 = C0(h, k);
        double c2 = C2(h, k);
        double c4 = C4(h, k);
        double d2 = D2(h, k);
        double d4 = D4(h, k);

        double eps = k*H/2.0;
        double kd = k*h;

        double err = -2.0*PII/T/sqrt(G*k) + c0 + eps*eps*(c2+d2/kd) + pow(eps, 4)*(c4+d4/kd);

        return err;
    }

    int StokesSolver (double H, double d, double T, double* kOut, double* errorOut )
    {
        double tolerance = 0.0000000001;
        double ErrI = 9999.0;
        double ErrM = 9999.0;

        double KI = 2.0*PII/StokesIFun::waveLength(d, T);
        double KD = 40*KI;
        double Kmid = (KI+KD)/2.0;

        while ( fabs(ErrM) > tolerance )
        {
            ErrI = error(H, d, KI, T);
            Kmid = (KI+KD)/2.0;
            ErrM = error(H, d, Kmid, T);

            if ( ErrI*ErrM >= 0 ) // same sign
            {
                KI = Kmid;
                ErrI = ErrM;
            }
            else
            {
                KD = Kmid;
            }
        }

        *kOut = Kmid;
        *errorOut = ErrM;

        return 1;
    }

    double B22 (double h, double k)
    {
        double CTH = 1.0/tanh(k*h);
        double s = 1.0/cosh(2.0*k*h);

        double B = CTH*(1.0 + 2.0*s)/(2.0*(1.0-s));

        return B;
    }

    double B31 (double h, double k)
    {
        double s = 1.0/cosh(2.0*k*h);

        double B = 
            -3.0*(1.0 + 3.0*s + 3.0*s*s + 2.0*pow(s,3))
            /(8.0*pow((1.0-s), 3));

        return B;
    }

    double B42 (double h, double k)
    {
        double CTH = 1.0/tanh(k*h);
        double s = 1.0/cosh(2.0*k*h);

        double B = 
            CTH*(6.0 - 26.0*s - 182.0*s*s - 204.0*pow(s,3) - 25.0*pow(s,4) + 26.0*pow(s,5))
            /(6.0*(3.0+2.0*s)*pow((1.0-s),4));

        return B;
    }

    double B44 (double h, double k)
    {
        double CTH = 1.0/tanh(k*h);
        double s = 1.0/cosh(2.0*k*h);

        double B = 
            CTH*(24.0 + 92.0*s + 122.0*s*s + 66.0*pow(s,3) + 67.0*pow(s,4) + 34.0*pow(s,5))
            /(24.0*(3.0+2.0*s)*pow((1.0-s),4));

        return B;
    }

    double B53 (double h, double k)
    {
        double s = 1.0/cosh(2.0*k*h);

        double B = 
            9.0*(132.0 + 17.0*s - 2216.0*s*s - 5897.0*pow(s,3) - 6292.0*pow(s,4) - 2687.0*pow(s,5) 
            + 194.0*pow(s,6) + 467.0*pow(s,7) + 82.0*pow(s,8))
            /(128.0*(3.0+2.0*s)*(4.0+s)*pow((1.0-s),6));

        return B;
    }

    double B55 (double h, double k)
    {
        double s = 1.0/cosh(2.0*k*h);

        double B = 
            5.0*(300.0 + 1579.0*s + 3176.0*s*s + 2949.0*pow(s,3) + 1188.0*pow(s,4) + 675.0*pow(s,5) 
            + 1326.0*pow(s,6) + 827.0*pow(s,7) + 130.0*pow(s,8))
            /(384.0*(3.0+2.0*s)*(4.0+s)*pow((1.0-s),6));

        return B;
    }

    double eta (double H, double d, double kx, double ky, double T, double x, double y, double t, double phase)
    {
        // Eta Stokes V: Fenton 85
        double k = sqrt(kx*kx + ky*ky);
        double eps = k*H/2.0;

        double b22 = B22(d, k);
        double b31 = B31(d, k);
        double b42 = B42(d, k);
        double b44 = B44(d, k);
        double b53 = B53(d, k);
        double b55 = B55(d, k);

        double theta = kx*x+ky*y-2.0*PII/T*t+phase;

        double freeS = 0;        
        
        freeS = eps*cos(theta) + pow(eps,2)*b22*cos(2.0*theta) + pow(eps,3)*b31*(cos(theta)-cos(3.0*theta))
                + pow(eps,4)*(b42*cos(2.0*theta)+b44*cos(4.0*theta))
                + pow(eps,5)*(-(b53+b55)*cos(theta) + b53*cos(3.0*theta) + b55*cos(5.0*theta));
        freeS = freeS/k;

        return freeS;
    }

}


namespace secondOrderFun
{
    #define PII 3.1415926535897932384626433832795028
    #define grav 9.81

    double C (double sigma1, double sigma2, double alphaSO1, double alphaSO2)
    {
        double denom = pow(sigma1, 2)*(pow(alphaSO1, 2)-1.0)-2.0*sigma1*sigma2*(alphaSO1*alphaSO2-1.0) + pow(sigma2, 2)*(pow(alphaSO2, 2)-1.0);
        
        if(denom == 0.0)
        {
            // Both components with 1.0/tanh(kh) = 1, very deep water
            return 0.0;
        }
        else
        {
            double value = (  ( 2.0*sigma1*sigma2*(sigma1-sigma2)*(1.0+alphaSO1*alphaSO2) + pow(sigma1, 3)*(pow(alphaSO1, 2)-1.0) - pow(sigma2, 3)*(pow(alphaSO2, 2)-1.0) )*(sigma1-sigma2)*(alphaSO1*alphaSO2-1.0)  )
                / denom + pow(sigma1, 2) + pow(sigma2, 2) - sigma1*sigma2*(alphaSO1*alphaSO2+1.0);
            return value;
        }

//( (2.0*sigma1*sigma2*(sigma1-sigma2)*(1.0+alphaSO1*alphaSO2) +
 //            pow(sigma1,3)*(pow(alphaSO1,2)-1.0)-pow(sigma2,3)*(pow(alphaSO2,2)-1.0) )*(sigma1-sigma2)*
  //           (alphaSO1*alphaSO2-1.0))/( pow(sigma1,2)*(pow(alphaSO1,2)-1.0)-2.0*sigma1*sigma2*
   //          (alphaSO1*alphaSO2-1.0)+pow(sigma2,2)*(pow(alphaSO2,2)-1.0)) + (pow(sigma1,2)+pow(sigma2,2))-
    //         sigma1*sigma2*(alphaSO1*alphaSO2+1.0);

        return 0.0;
    }

    double E (double a1, double a2, double sigma1, double sigma2, double alphaSO1, double alphaSO2)
    {
         double value = -0.5*a1*a2*( 2.0*sigma1*sigma2*(sigma1-sigma2)*(1.0+alphaSO1*alphaSO2) + pow(sigma1, 3)*(pow(alphaSO1, 2)-1.0) - pow(sigma2, 3)*(pow(alphaSO2, 2)-1.0) );
//-0.5*a1*a2*(2*sigma1*sigma2*(sigma1-sigma2)*(1.0+alphaSO1*alphaSO2)+
 //          pow(sigma1,3)*(pow(alphaSO1,2)-1.0)-pow(sigma2,3)*(pow(alphaSO2,2)-1.0));

        return value;
    }

    double etaSO (double H1, double H2, double sigma1, double sigma2, double phase1, double phase2, double kx1, double ky1, double kx2, double ky2, double x, double y, double t, double h)
    {
        double k1 = sqrt( pow(kx1,2) + pow(ky1,2) );
        double k2 = sqrt( pow(kx2,2) + pow(ky2,2) );

        double alphaSO1 = 1.0/tanh(k1*h);
        double alphaSO2 = 1.0/tanh(k2*h);

        double fullPhase1 =  kx1*x + ky1*y - sigma1*t + phase1;
        double fullPhase2 =  kx2*x + ky2*y - sigma2*t + phase2;

        double CSO = C(sigma1, sigma2, alphaSO1, alphaSO2);

        double eta = H1/2.0*H2/2.0/(2.0*grav)*CSO*cos(fullPhase1-fullPhase2);

        return eta;
    }

    double uSO (double H1, double H2, double sigma1, double sigma2, double phase1, double phase2, double kx1, double ky1, double kx2, double ky2, double x, double y, double t, double h, double z)
    {
        double k1 = sqrt( pow(kx1,2) + pow(ky1,2) );
        double k2 = sqrt( pow(kx2,2) + pow(ky2,2) );

        double alphaSO1 = 1.0/tanh(k1*h);
        double alphaSO2 = 1.0/tanh(k2*h);

        double fullPhase1 =  kx1*x + ky1*y - sigma1*t + phase1;
        double fullPhase2 =  kx2*x + ky2*y - sigma2*t + phase2;

        double ESO = E(H1/2.0, H2/2.0, sigma1, sigma2, alphaSO1, alphaSO2);

        double diff = sigma1-sigma2;

        double U = ESO*cosh((k1-k2)*(z))*cos(fullPhase1-fullPhase2)*(k1-k2)/(grav*(k1-k2)*sinh((k1-k2)*h)-pow(diff,2)*cosh((k1-k2)*h));

        return U;
    }

    double wSO (double H1, double H2, double sigma1, double sigma2, double phase1, double phase2, double kx1, double ky1, double kx2, double ky2, double x, double y, double t, double h, double z)
    {
        double k1 = sqrt( pow(kx1,2) + pow(ky1,2) );
        double k2 = sqrt( pow(kx2,2) + pow(ky2,2) );

        double alphaSO1 = 1.0/tanh(k1*h);
        double alphaSO2 = 1.0/tanh(k2*h);

        double fullPhase1 =  kx1*x + ky1*y - sigma1*t + phase1;
        double fullPhase2 =  kx2*x + ky2*y - sigma2*t + phase2;

        double ESO = E(H1/2.0, H2/2.0, sigma1, sigma2, alphaSO1, alphaSO2);

        double diff = sigma1-sigma2;

        double W = ESO*sinh((k1-k2)*(z))*sin(fullPhase1-fullPhase2)*(k1-k2)/(grav*(k1-k2)*sinh((k1-k2)*h)-pow(diff,2)*cosh((k1-k2)*h));

        return W;
    }
}

namespace BoussinesqFun
{
    #define PII 3.1415926535897932384626433832795028
    #define G 9.81 
    
    double eta (double H, double h, double x, double y, double theta, double t, double X0)
    {
        double C = sqrt(G*(H+h));
        double ts = 3.5*h/sqrt(H/h);
        double aux = sqrt(3.0*H/(4.0*h))/h;
        double Xa = -C * t + ts - X0 + x * cos(theta) + y * sin(theta);

        double sup = H * 1.0/pow(cosh( aux * Xa ),2);

        return sup;
    }

    double Deta1 (double H, double h, double x, double y, double theta, double t, double X0)
    {
        double C = sqrt(G*(H+h));
        double ts = 3.5*h/sqrt(H/h);
        double aux = sqrt(3.0*H/(4.0*h))/h;
        double Xa = -C * t + ts - X0 + x * cos(theta) + y * sin(theta);

        double deta = 8.0*aux*h * exp(2.0*aux*Xa) * (1.0-exp(2.0*aux*Xa)) 
                        / pow(1.0+exp(2.0*aux*Xa),3);

        return deta;
    }

    double Deta2 (double H, double h, double x, double y, double theta, double t, double X0)
    {
        double C = sqrt(G*(H+h));
        double ts = 3.5*h/sqrt(H/h);
        double aux = sqrt(3.0*H/(4.0*h))/h;
        double Xa = -C * t + ts - X0 + x * cos(theta) + y * sin(theta);

        double deta = 16.0*pow(aux,2)*h * exp(2.0*aux*Xa) * (exp(4.0*aux*Xa)
                        - 4.0*exp(2.0*aux*Xa)+1.0) / pow(1.0+exp(2.0*aux*Xa),4);

        return deta;
    }

    double Deta3 (double H, double h, double x, double y, double theta, double t, double X0)
    {
        double C = sqrt(G*(H+h));
        double ts = 3.5*h/sqrt(H/h);
        double aux = sqrt(3.0*H/(4.0*h))/h;
        double Xa = -C * t + ts - X0 + x * cos(theta) + y * sin(theta);

        double deta = -32.0*pow(aux,3)*h * exp(2.0*aux*Xa) * (exp(6.0*aux*Xa)
                        - 11.0*exp(4.0*aux*Xa) + 11.0*exp(2.0*aux*Xa)-1.0) / pow(1.0+exp(2.0*aux*Xa),5);

        return deta;
    }

    double U (double H, double h, double x, double y, double theta, double t, double X0, double z)
    {
        double C = sqrt(G*(H+h));
        double etaSolit = eta( H, h, x, y, theta, t, X0);
        double Detas2 = Deta2( H, h, x, y, theta, t, X0);

        double vel = C*etaSolit/h
                *(1.0 - etaSolit/(4.0*h) +
                pow(h,2)/(3.0*etaSolit)
                *(1.0 - 3.0/2.0*pow(z/h,2))
                *Detas2);

        return vel;
    }

    double W (double H, double h, double x, double y, double theta, double t, double X0, double z)
    {
        double C = sqrt(G*(H+h));
        double etaSolit = eta( H, h, x, y, theta, t, X0);
        double Detas1 = Deta1( H, h, x, y, theta, t, X0);
        double Detas3 = Deta3( H, h, x, y, theta, t, X0);

        double vel = -C*z/h
                *((1.0 - etaSolit/(2.0*h))*Detas1 +
                pow(h,2)/3.0
                *(1.0 - 1.0/2.0*pow(z/h,2))
                *Detas3);

        return vel;
    }
}

