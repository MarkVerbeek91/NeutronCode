/**
 * /brief These functions calculate the neutral particle flux on a given
 * radius (r) for a given energy (E). This are two main functions. One for
 * inward traveling ions and another for outward traving ions.
 */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "constants.h"
#include "MathFunctions.h"
#include "PotentialFunctions.h"
#include "SurvivalFunctions.h"
#include "CrossSections.h"

#include "NeutralsClassISpectrum.h"

// equation 40 and 47
double NeutralsClassISpectrumInwards(double r, double E)
{
    double flux, dr;
    double term1 = 0, term2 = 0;
	  double (*PhiPtr)(double);
    PhiPtr = &Potential_Phi;

    dr = Potential_Phi_Inv((E0 - E)/giveq());

    if ( dr < r )
        dr = roundf(dr * 10000) / 10000 + 1e-5;

    term1  = 1 / giveq();
    term1 *= pow(giveAnodeRadius()/r,2);
    term1 *= ngas * CrosssecCX(E);
    term1 *= EdgeIonFlux;
    term1 /= abs(differentiat(*PhiPtr, dr));
    term1 *= f(dr);

    if ( giveCathodeRadius() < r )
    {
            // make sure that only physical numbers are calculated
        if ( dr < r || dr > giveAnodeRadius() )
        {
            printf("NeutralsClassISpectrumInwards error: r < dr or dr < a\n");
            return NAN;
        }

        flux = term1;
    }
    else
    {
            // make sure that only physical numbers are calculated (47)
        if ( dr < giveCathodeRadius() )
        {
            printf("NeutralsClassISpectrumInwards error: r < dr or dr < a\n");
            return NAN;
        }


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
    double flux, dr;
    double term1 = 0, term2 = 0;
    double (*PhiPtr)(double);
    PhiPtr = &Potential_Phi;

    dr = Potential_Phi_Inv((E0 - E)/giveq());

    if ( dr < r )
        dr = roundf(dr * 10000) / 10000 + 1e-5;

    term1  = 1 / giveq();
    term1 *= pow(giveAnodeRadius()/r,2);
    term1 *= ngas * CrosssecCX(E);
    term1 /= abs(differentiat(*PhiPtr, dr));

    if ( r < giveCathodeRadius() )
    {
        if ( dr < giveCathodeRadius() )
        {
            fprintf(stderr,"NeutralsClassISpectrumOutwards error: dr < a:\n\tr = %E\n \tdr = %E\n \tE = %E\n", r, dr, E);
            return NAN;
        }


        term1 *= f(dr);

        if ( DELTA(E + giveVoltage()) )
        {
            term2  = pow(giveAnodeRadius()/r,2) * EdgeIonFlux;
            term2 *= f(giveCathodeRadius())* ( 1 - exp( ngas * CrosssecCX(-giveVoltage()) * (r - giveCathodeRadius())));
        }

        flux = giveTransparency() * (term1 + term2);
    }
    else
    {
        if ( dr < giveCathodeRadius() || dr < r )
        {
            fprintf(stderr,"NeutralsClassISpectrumOutwards error: dr < a || dr < r:\n\tr = %E\n \tdr = %E\n \tE = %E\n", r, dr, E);
            return NAN;
        }

        term1 *= f(dr) + pow(f(0),2) / f(dr);

        if ( DELTA(E + giveVoltage()) )
        {
            term2  = pow(giveAnodeRadius()/r,2) * EdgeIonFlux;
            term2 *= f(giveCathodeRadius())* ( 1 - exp( ngas * CrosssecCX(-giveVoltage()) * giveCathodeRadius()));
        }

        flux = pow(giveTransparency(),2) * (term1 + term2);
    }

    return flux;
}
