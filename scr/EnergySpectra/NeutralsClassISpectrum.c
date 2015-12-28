#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "constants.h"
#include "MathFunctions.h"
#include "PotentialFunctions.h"
#include "SurvivalFunctions.h"
#include "Crosssections.h"

#include "NeutralsClassISpectrum.h"

// equation 40 and 47
double NeutralsClassISpectrumInwards(double r, double E)
{
    double flux;
    double term1 = 0, term2 = 0;
    double dr = Potential_Phi_Inv(E);

    // make sure that only physical numbers are calculated
    if ( giveCathodeRadius() < dr || r < dr || dr < giveAnodeRadius()) // @@@@@@ r < dr ?=> dr < r
    {
        printf("NeutralsClassISpectrumInwards error: a < dr, r < dr or dr < a\n");
        return NAN;
    }


    double (*PhiPtr)(double);
    PhiPtr = &Potential_Phi;

    term1  = 1 / giveq();
    term1 *= pow(giveAnodeRadius()/r,2);
    term1 *= ngas * CrosssecCX(E);
    term1 *= EdgeIonFlux;
    term1 /= abs(differentiat(*PhiPtr, dr));
    term1 *= f(dr);

    if ( giveCathodeRadius() < r )
    {
        flux = term1;
    }
    else
    {
        if ( DELTA(E + giveVoltage()) )
        {
            term2  = pow(giveAnodeRadius()/r,2) * EdgeIonFlux;
            term2 *= f(giveCathodeRadius())* ( 1 - exp(-2 * ngas * CrosssecCX(-giveVoltage()) * giveCathodeRadius()));
        }

        flux = giveTransparency() * (term1 + term2);
    }

    return flux;
}

// equation 51 and 55
double NeutralsClassISpectrumOutwards(double r, double E)
{
    double flux;
    double term1 = 0, term2 = 0;
    double dr = Potential_Phi_Inv(E);

    if (  giveCathodeRadius() < dr || r < dr || dr < giveAnodeRadius() ) // @@@@ r < dr ?=> dr < r
    {
        printf("NeutralsClassISpectrumOutwards error: a < dr, r < dr or dr < a\n");
        return NAN;
    }

    double (*PhiPtr)(double);
    PhiPtr = &Potential_Phi;

    term1  = 1 / giveq();
    term1 *= pow(giveAnodeRadius()/r,2);
    term1 *= ngas * CrosssecCX(E);
    term1 /= abs(differentiat(*PhiPtr, dr));

    if ( giveCathodeRadius() > r )
    {
        term1 *= f(dr);

        if ( DELTA(E + giveVoltage()) )
        {
            term2  = pow(giveAnodeRadius()/r,2);
            term2 *= giveTransparency();
            term2 *= f(giveCathodeRadius())* ( 1 - exp( ngas * CrosssecCX(-giveVoltage()) * (r - giveCathodeRadius())));
        }

        flux = giveTransparency() * (term1 + term2);
    }
    else
    {
        term1 *= f(dr) + pow(f(0),2) / f(dr);

        if ( DELTA(E + giveVoltage()) )
        {
            term2  = pow(giveAnodeRadius()/r,2);
            term2 *= giveTransparency();
            term2 *= f(giveCathodeRadius())* ( 1 - exp( ngas * CrosssecCX(-giveVoltage()) * giveCathodeRadius()));
        }

        flux = pow(giveTransparency(),2) * (term1 + term2);
    }

    return flux;
}








