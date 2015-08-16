
#include <math.h>
#include <stdlib.h>

#include "constants.h"
#include "MathFunctions.h"
#include "PotentialFunctions.h"
#include "SurvivalFunctions.h"
#include "Crosssections.h"

#include "NeutralsClassISpectrum.h"

/**
        INWARD TRAVELING NEUTRALS FROM CLASS II IONS
*/

// intregral of neutral from claas II ions _outside_ the cathode _inwards_
double NeutralsClassIISpectrumInwards_Inte1(double E, double ddr)
{
    dr = Potential_Phi_Inv(Potential_Phi(ddr) - E/giveq());

    if ( dr < r )
    {
        printf("NeutralsClassIISpectrumInwards_Inte1 error: dr < r \n");
        return NAN;
    }
    double integrant;

    integrant  = interpolation(ddr);
    integrant /= abs(differentiat(*PhiPtr, dr));
    integrant *= g(dr, ddr);
    integrant /= 1 - pow(giveTransparency() * g(0,ddr),2);
    integrant *= pow(ddr,2);

    return integrant;
}

// intregral of neutral from claas II ions _inside_ the cathode _inwards_
double NeutralsClassIISpectrumInwards_Inte2(double E, double ddr)
{
    dr = Potential_Phi_Inv(Potential_Phi(ddr) - E/giveq());

    if ( giveCathodeRadius() < dr )
    {
        printf("NeutralsClassIISPectrumInwards error: d < dr \n");
        return NAN;
    }


    double integrant;

    integrant  = interpolation(ddr);
    integrant /= abs(differentiat(*PhiPtr, dr));
    integrant *= g(dr, ddr);
    integrant /= 1 - pow(giveTransparency() * g(0,ddr),2);
    integrant *= pow(ddr,2);

    return integrant;
}

double NeutralsClassIISpectrumInwards (double r, double E)
{
    double flux;
    double term1 = 0, term2 = 0;

    double (*PhiPtr)(double);
    PhiPtr = &Potential_Phi;

    term1  = 1 / giveq();
    term1 /= pow(r,2);
    term1 *= ngas * CrosssecCX(E);

    if ( giveCathodeRadius() < r)
    {
        // outside the cathode region
        double (*FunctPtr)(double, double);
        FunctPtr = &NeutralsClassIISpectrumInwardsInte1;

        term1 *= NIntegration_2(*FunctPtr, E, r, giveAnodeRadius());

        flux = term1;
    }
    else
    {
        // inside the cathode region
        double (*FunctPtr)(double, double);
        FunctPtr = &NeutralsClassIISpectrumInwardsInte2;

        term1 *= NIntegration_2(FunctPtr, E, giveCathodeRadius(), giveAnodeRadius());

        dr = Potential_Phi_Inv(Potential_Phi(ddr) - E/giveq());
        drr = Potential_Phi_Inv(Potential_Phi(giveCathodeRadius() - E / giveq()));

        term2  = pow(drr/r, 2);
        term2 *= g(giveCathodeRadius(), ddr);
        term2 /= 1 - pow(giveTransparency() * g(0, ddr));
        term2 *= 1 - exp(ngas * CrosssecCX(E) * ( r - giveCathodeRadius()) );
        term2 *= interpolation(ddr) / ( giveq() * abs(differentiat(*PhiPtr, dr)));

        flux = giveTransparency() * ( term1 + term2);
    }

    return flux;
}

/**
        OUTWARD TRAVELING NEUTRALS FROM CLASS II IONS
*/

// intregral of neutrals from claas II ions _outside_ the cathode _outward_
double NeutralsClassIISpectrumOutwards_Inte1(double E, double ddr)
{
    dr = Potential_Phi_Inv(Potential_Phi(ddr) - E/giveq());

    if ( dr < r )
    {
        printf("NeutralsClassIISpectrumOutwards_Inte1 error: dr < r \n");
        return NAN;
    }


    double integrant;

    integrant  = interpolation(ddr);
    integrant /= abs(differentiat(*PhiPtr, dr));
    integrant *= g(giveCathodeRadius(), ddr);
    integrant /= 1 - pow(giveTransparency() * g(0,ddr),2);
    integrant *= pow(ddr,2);

    return integrant;
}

// intregral of neutral from claas II ions _inside_ the cathode _outwards_
double NeutralsClassIISpectrumOutwards_Inte2(double E, double r, double ddr)
{
    dr = Potential_Phi_Inv(Potential_Phi(ddr) - E/giveq());

    double intergrant = 0;
    if ( giveCathodeRadius() < dr )
        intergrant = g(dr, ddr);

    if ( dr < r )
        intergrant += pow(g(0,ddr),2) / g(dr, ddr);

    if ( intergrant == 0)
        return 0;

    integrant  = interpolation(ddr);
    integrant /= abs(differentiat(*PhiPtr, dr));
    integrant /= 1 - pow(giveTransparency() * g(0,ddr),2);
    integrant *= pow(ddr,2);

    return integrant;
}

double NeutralsClassIISpectrumOutwards (double r, double E)
{
    double flux;
    double term1 = 0, term2 = 0;

    double (*PhiPtr)(double);
    PhiPtr = &Potential_Phi;

    term1  = 1 / giveq();
    term1 /= pow(r,2);
    term1 *= ngas * CrosssecCX(E);

    if ( r < giveCathodeRadius() )
    {
        // inside the cathode region

        double (*FunctPtr)(double, double);
        FunctPtr = &NeutralsClassIISpectrumOutwards_Inte2;

        term1 *= NIntegration_2(*FunctPtr, E, r, giveAnodeRadius());

        dr = Potential_Phi_Inv(Potential_Phi(ddr) - E/giveq());
        drr = Potential_Phi_Inv(Potential_Phi(giveCathodeRadius() - E / giveq()));

        term2  = pow(drr/r, 2);
        term2 *= g(giveCathodeRadius(), ddr);
        term2 /= 1 - pow(giveTransparency() * g(0, ddr));
        term2 *= 1 - exp(- ngas * CrosssecCX(E) * ( r + giveCathodeRadius()) );
        term2 *= interpolation(ddr) / ( giveq() * abs(differentiat(*PhiPtr, dr)));

        flux = giveTransparency() * ( term1 + term2);
    }
    else
    {
        // outside the cathode region
        double (*FunctPtr)(double, double);
        FunctPtr = &NeutralsClassIISpectrumOutwards_Inte1;

        term1 *= NIntegration_3(FunctPtr, E, r, giveCathodeRadius(), giveAnodeRadius());

        dr = Potential_Phi_Inv(Potential_Phi(ddr) - E/giveq());
        drr = Potential_Phi_Inv(Potential_Phi(giveCathodeRadius() - E / giveq()));

        term2  = pow(drr/r, 2);
        term2 *= g(giveCathodeRadius(), ddr);
        term2 /= 1 - pow(giveTransparency() * g(0, ddr));
        term2 *= 1 - exp(-2 * ngas * CrosssecCX(E) * giveCathodeRadius());
        term2 *= interpolation(ddr) / ( giveq() * abs(differentiat(*PhiPtr, dr)));

        flux = pow(giveTransparency(),2) * ( term1 + term2);
    }

    return flux;
}
