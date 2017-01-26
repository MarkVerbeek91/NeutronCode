/**
 * /brief These functions calculate the neutral particle flux on a given
 * radius (r) for a given energy (E). The neutrals come from the Class II
 * ions. This are two main functions. One for inward traveling ions and
 * another for outward traving ions.
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

/**
 *         INWARD TRAVELING NEUTRALS FROM CLASS II IONS
 */

// intregral of neutral from claas II ions _outside_ the cathode _inwards_
double NeutralsClassIISpectrumInwards_Inte1(double r, double E, double ddr)
{
    double integrant, dr;
    double (*PhiPtr)(double);
    PhiPtr = &Potential_Phi;

    dr = Potential_Phi_Inv(Potential_Phi(ddr) - E/giveq());

    if ( dr < r || ddr < dr )  // last part is a bit hacky
        return 0;

    integrant  = interpolation(Table->S, ddr);
    integrant /= abs(differentiat(*PhiPtr, dr));
    integrant *= g(dr, ddr);
    integrant /= 1 - pow(giveTransparency() * g(0,ddr),2);
    integrant *= pow(ddr,2);

    return integrant;
}

// intregral of neutral from claas II ions _inside_ the cathode _inwards_
double NeutralsClassIISpectrumInwards_Inte2(double E, double ddr)
{
    double integrant, dr;
    double (*PhiPtr)(double);
    PhiPtr = &Potential_Phi;

    // from Eq. 56
    dr = Potential_Phi_Inv(Potential_Phi(ddr) - E/giveq() );

    if ( dr < fusor->a || ddr < dr ) // last part is a bit hacky
        return 0;

    integrant  = interpolation(Table->S, ddr);
    integrant /= abs(differentiat(*PhiPtr, dr));
    integrant *= g(dr, ddr);
    integrant /= 1 - pow(giveTransparency() * g(0,ddr),2);
    integrant *= pow(ddr,2);

    return integrant;
}

// equation 58 and 65
double NeutralsClassIISpectrumInwards (double r, double E)
{
    double flux;
    double term1 = 0, term2 = 0;

    double (*PhiPtr)(double);
    PhiPtr = &Potential_Phi;

    term1  = 1 / giveq();
    term1 /= pow(r,2);
    term1 *= ngas * CrosssecCX(E);

    if ( fusor->a < r)
    {
        double (*FunctPtr)(double, double, double);
        FunctPtr = &NeutralsClassIISpectrumInwards_Inte1;

        term1 *= NIntegration_3(*FunctPtr, r, E, r, fusor->b);

        flux = term1;
    }
    else
    {
        // inside the cathode region
        if ( false )
        {
            printf("NeutralsClassIISpectrumInwards error: r < dr or dr < a\n");
            return NAN;
        }

        double ddr, dr;
        double (*FunctPtr)(double, double);
        FunctPtr = &NeutralsClassIISpectrumInwards_Inte2;

        term1 *= NIntegration_2(FunctPtr, E, fusor->a, fusor->b);

        ddr = Potential_Phi_Inv(Potential_Phi(fusor->a - E / giveq()));
        dr = Potential_Phi_Inv(Potential_Phi(ddr) - E/giveq());

        term2  = pow(ddr/r, 2);
        term2 *= g(fusor->a, ddr);
        term2 /= 1 - pow(giveTransparency() * g(0, ddr),2);
        term2 *= 1 - exp(ngas * CrosssecCX(E) * ( r - fusor->a) );
        term2 *= interpolation(Table->S, ddr) / ( giveq() * abs(differentiat(*PhiPtr, dr)));

        flux = giveTransparency() * ( term1 + term2 );
    }

    return flux;
}

/**
        OUTWARD TRAVELING NEUTRALS FROM CLASS II IONS
*/

// intregral of neutrals from claas II ions _outside_ the cathode _outward_
double NeutralsClassIISpectrumOutwards_Inte1(double E, double ddr)
{
    double integrant, dr;
    double (*PhiPtr)(double);
    PhiPtr = &Potential_Phi;

    dr = Potential_Phi_Inv(Potential_Phi(ddr) - E/giveq());

    if ( dr < fusor->a || ddr < fusor->a ) // last part is a bit hacky
        return 0;

    integrant  = interpolation(Table->S, ddr);
    integrant /= abs(differentiat(*PhiPtr, dr));
    integrant *= g(fusor->a, ddr);
    integrant /= 1 - pow(giveTransparency() * g(0,ddr),2);
    integrant *= pow(ddr,2);

    return integrant;
}

// intregral of neutral from claas II ions _inside_ the cathode _outwards_
double NeutralsClassIISpectrumOutwards_Inte2(double r, double E, double ddr)
{
    double integrant = 0, dr;
    double (*PhiPtr)(double);
    PhiPtr = &Potential_Phi;

    dr = Potential_Phi_Inv(Potential_Phi(ddr) - E/giveq());

    if ( dr < fusor->a )
        integrant = g(dr, ddr);

    if ( r < dr )
        integrant += pow(g(0,ddr),2) / g(dr, ddr);

    if ( integrant == 0)
        return 0;

    integrant  = interpolation(Table->S, ddr);
    integrant /= abs(differentiat(*PhiPtr, dr));
    integrant /= 1 - pow(giveTransparency() * g(0,ddr),2);
    integrant *= pow(ddr,2);

    return integrant;
}

double NeutralsClassIISpectrumOutwards (double r, double E)
{
    double flux, dr, ddr;
    double term1 = 0, term2 = 0;

    double (*PhiPtr)(double);
    PhiPtr = &Potential_Phi;

    term1  = 1 / giveq();
    term1 /= pow(r,2);
    term1 *= ngas * CrosssecCX(E);

    if ( r < fusor->a )
    {
        // inside the cathode region
        double (*FunctPtr)(double, double);
        FunctPtr = &NeutralsClassIISpectrumOutwards_Inte1;

        term1 *= NIntegration_2(*FunctPtr, E, fusor->a, fusor->b);

        ddr = Potential_Phi_Inv(Potential_Phi(fusor->a - E / giveq()));
        dr = Potential_Phi_Inv(Potential_Phi(ddr) - E/giveq());

        term2  = pow(ddr/r, 2);
        term2 *= g(fusor->a, ddr);
        term2 /= 1 - pow(giveTransparency() * g(0, ddr),2);
        term2 *= 1 - exp(- ngas * CrosssecCX(E) * ( r + fusor->a) );
        term2 *= interpolation(Table->S, ddr) / ( giveq() * abs(differentiat(*PhiPtr, dr)));

        flux = giveTransparency() * ( term1 + term2);
    }
    else
    {
        // outside the cathode region
        double (*FunctPtr)(double, double, double);
        FunctPtr = &NeutralsClassIISpectrumOutwards_Inte2;

        term1 *= NIntegration_3(FunctPtr, r, E, fusor->a, fusor->b);

        ddr = Potential_Phi_Inv(Potential_Phi(fusor->a - E / giveq()));
        dr = Potential_Phi_Inv(Potential_Phi(ddr) - E/giveq());

        term2  = pow(ddr/r, 2);
        term2 *= g(fusor->a, ddr);
        term2 /= 1 - pow(giveTransparency() * g(0, ddr),2);
        term2 *= 1 - exp(-2 * ngas * CrosssecCX(E) * fusor->a);
        term2 *= interpolation(Table->S, ddr) / ( giveq() * abs(differentiat(*PhiPtr, dr)));

        flux = pow(giveTransparency(),2) * ( term1 + term2);
    }

    return flux;
}
