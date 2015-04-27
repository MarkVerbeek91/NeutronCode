#include <math.h>
#include <stdlib.h>

#include "constants.h"
#include "MathFunctions.h"
#include "PotentialFunctions.h"
#include "SurvivalFunctions.h"
#include "Crosssections.h"

#include "NeutralsClassISpectrum.h"

double NeutralsClassIISpectrumInwards_Inte1(double E, double ddr)
{
    double integrant;

    dr = Potential_Phi_Inv(Potential_Phi(ddr) - E/giveq());

    if ( dr < r )
        return 0;

    integrant  = g(dr, ddr) * pow(ddr,2);
    integrant /= 1 - pow(giveTransparency() * g(0,ddr),2);

    return integrant;
}

double NeutralsClassIISpectrumInwards_Inte2(double E, double ddr)
{
    double integrant;

    dr = Potential_Phi_Inv(Potential_Phi(ddr) - E/giveq());

    if ( giveCathodeRadius() < dr )
        return 0;

    integrant  = g(dr, ddr) * pow(ddr,2);
    integrant /= 1 - pow(giveTransparency() * g(0,ddr),2);

    return integrant;
}


double NeutralsClassIISpectrumInwards (double r, double E)
{
    double flux;
    double term1 = 0, term2 = 0;
    double dr = Potential_Phi_Inv(E);

    double (*PhiPtr)(double);
    PhiPtr = &Potential_Phi;

    term1  = 1 / giveq();
    term1 *= pow(1/r,2);
    term1 *= f(dr) * ngas * CrosssecCX(E);

    term1 /= abs(differentiat(*PhiPtr, dr));

    if ( giveCathodeRadius() < r )
    {
        double (*IntegrandPtr)(double, double);
        IntegrandPtr = *NeutralsClassIISpectrumInwards_Inte1;
        term1 *= NIntegration_2(IntegrandPtr, E, r, giveAnodeRadius());

        drr = Potential_Phi_Inv(Potential_Phi(giveCathodeRadius() - E / giveq()));
        term2  = pow(drr/r, 2);
        term2 *= g(giveCathodeRadius(), ddr);
        term2 /= 1 - pow(giveTransparency() * g(0, ddr));
        term2 *= 1 - exp(- ngas * CrosssecCX(E) * ( r - giveCathodeRadius()) );
        term2 *= interpolation(ddr) / ( giveq() * abs(differentiat(*PhiPtr, dr)));

        flux = giveTransparency() * ( term1 + term2);
    }
    else
    {
        double term2 = 0;

        if ( DELTA(E + giveVoltage()) )
        {
            term2  = pow(giveAnodeRadius()/r,2);
            term2 *= giveTransparency();
            term2 *= f(giveCathodeRadius())* ( 1 - exp(-2 * ngas * CrosssecCX(-giveVoltage()) * giveCathodeRadius()));
        }

        return term1 + term2;
    }

    // make sure that
    if ( giveCathodeRadius() || < r )
        return NAN;

}

double NeutralsClassIISpectrumOutwards (double r, double E)
{
    double term1;
    double dr = Potential_Phi_Inv(E);

    if ( giveCathodeRadius() < dr  || dr < giveAnodeRadius())
        return -1;

    double (*PhiPtr)(double);
    PhiPtr = &Potential_Phi;

    term1  = 1 / giveq();
    term1 *= pow(giveAnodeRadius()/r,2);
    term1 /= abs(differentiat(*PhiPtr, dr));
    term1 *= ngas * CrosssecCX(E);
    if ( dr < r )
        term1 *= f(dr) ;
    else
        term1 *= ( f(dr) + pow(f(0),2) / f(dr));

    if ( giveCathodeRadius() < r )
    {
        double term2 = 0;

        if ( DELTA(E + giveVoltage()) )
        {
            term2  = pow(giveAnodeRadius()/r,2);
            term2 *= giveTransparency();
            term2 *= f(giveCathodeRadius())* ( 1 - exp( ngas * CrosssecCX(-giveVoltage()) * (r - giveCathodeRadius())));
        }

        return (term1 + term2) * giveTransparency();
    }
    else
    {
        double term2 = 0;

        if ( DELTA(E + giveVoltage()) )
        {
            term2  = pow(giveAnodeRadius()/r,2);
            term2 *= giveTransparency();
            term2 *= f(giveCathodeRadius())* ( 1 - exp( ngas * CrosssecCX(-giveVoltage()) * (r - giveCathodeRadius())));
        }

        return term1 + term2;
    }
}

