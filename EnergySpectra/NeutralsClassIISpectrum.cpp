#include <math.h>
#include <stdlib.h>

#include "constants.h"
#include "MathFunctions.h"
#include "PotentialFunctions.h"
#include "SurvivalFunctions.h"
#include "Crosssections.h"

#include "NeutralsClassISpectrum.h"

double NeutralsClassIISpectrumInwards (double r, double E)
{
    double term1;
    double dr = r_shell(r, E);

    double (*PhiPtr)(double);
    PhiPtr = &Potential_Phi;

    term1  = 1 / giveq();
    term1 *= pow(1/r,2);
    term1 *= f(dr) * ngas * CrosssecCX(E);


    term1 /= abs(differentiat(*PhiPtr, dr));

    if ( giveCathodeRadius() < r )
        return term1;
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
    double dr = r_shell(r, E);

    if ( dr < giveCathodeRadius() | dr > giveAnodeRadius())
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

