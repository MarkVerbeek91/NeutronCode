#include <math.h>
#include <stdlib.h>

#include "constants.h"
#include "MathFunctions.h"
#include "PotentialFunctions.h"
#include "SurvivalFunctions.h"
#include "CrossSections.h"

#include "NeutronProductionIons.h"

/** \brief Neutron production per radii produced by ions in the fusor.

*/

/**
        NEUTRONS FROM IONS GOING INWARDS
*/
// outside the cathod
double NeutronsIonFluxInwards_Inte1(double r, double dr)
{
    double E, term1;

    E = giveq() * ( Potential_Phi(r) - Potential_Phi(dr));

    term1  = CrosssecFusion(E);
    term1 *= pow(dr/r,2);
    term1 *= interpolation(dr);
    term1 *= 1 / ( 1 - pow(giveTransparency() * g(0,dr),2));
    term1 *= g(r,dr);

    return term1;
}

// inside the cathod
double NeutronsIonFluxInwards_Inte2(double r, double dr)
{
    double E, term1;

    E = giveq() * ( Potential_Phi(r) - Potential_Phi(dr));

    term1  = CrosssecFusion(E);
    term1 *= pow(dr/r,2);
    term1 *= interpolation(dr);
    term1 *= 1 / ( 1 - pow(giveTransparency() * g(0,dr),2));
    term1 *= g(giveCathodeRadius(),dr) * exp(ngas * CrosssecCX(ParticleEnergy2(giveCathodeRadius(),dr)) * (r - giveCathodeRadius()));

    return r;
}

double NeutronsIonFluxInwards(double r)
{
    double NeutronFlux;
    double term1 = 0, term2 = 0;

    if ( giveCathodeRadius() < r )
    {
        double (*FunctPtr)(double, double);
        FunctPtr = &NeutronsIonFluxInwards_Inte1;

        term1 = NIntegration_2(FunctPtr, r, r, giveAnodeRadius());


        term2  = pow(giveAnodeRadius()/r,2) * EdgeIonFlux * f(r);
        term2 *= CrosssecFusion(giveq() * Potential_Phi(r));

        NeutronFlux = term1 + term2;
    }
    else
    {
        double (*FunctPtr)(double, double);
        FunctPtr = &NeutronsIonFluxInwards_Inte2;

        term1 = NIntegration_2(FunctPtr, r, r, giveAnodeRadius());

        term2  = pow(giveAnodeRadius()/r,2) * EdgeIonFlux ;
        term2 *= f(giveCathodeRadius()) * exp(ngas * CrosssecCX(ParticleEnergy1(giveCathodeRadius())) * (r - giveCathodeRadius()));
        term2 *= CrosssecFusion(giveq() * Potential_Phi(r));

        NeutronFlux = giveTransparency() * (term1 + term2);
    }

    return NeutronFlux;
}

/**
        NEUTRONS FROM IONS GOING OUTWARDS
*/

// inside the cathod
double NeutronsIonFluxOutwards_Inte1(double r, double dr)
{
    double E, term1;

    E = giveq() * ( Potential_Phi(r) - Potential_Phi(dr));

    term1  = CrosssecFusion(E);
    term1 *= pow(dr/r,2);
    term1 *= interpolation(dr);
    term1 *= 1 / ( 1 - pow(giveTransparency() * g(0,dr),2));
    term1 /= g(giveCathodeRadius(), dr) * exp(ngas * CrosssecCX(ParticleEnergy2(giveCathodeRadius(),dr)) * ( r - giveCathodeRadius()));

    return term1;}

// outside the cathod
double NeutronsIonFluxOutwards_Inte2(double r, double dr)
{
    double E, term1;

    E = giveq() * ( Potential_Phi(r) - Potential_Phi(dr));

    term1  = CrosssecFusion(E);
    term1 *= pow(dr/r,2);
    term1 *= interpolation(dr);
    term1 *= 1 / ( 1 - pow(giveTransparency() * g(0,dr),2));
    term1 /= g(r, dr);

    return term1;
}

double NeutronsIonFluxOutwards(double r)
{
    double NeutronFlux;
    double term1 = 0, term2 = 0;

    if ( giveCathodeRadius() < r )
    {
        double (*FunctPtr)(double, double);
        FunctPtr = &NeutronsIonFluxOutwards_Inte2;

        term1 = NIntegration_2(FunctPtr, r, r, giveAnodeRadius());

        term2  = pow(giveAnodeRadius() * f(0)/r,2) * EdgeIonFlux;
        term2 /= f(giveCathodeRadius()) * exp(ngas * CrosssecCX(ParticleEnergy1(giveCathodeRadius())) * (r - giveCathodeRadius()));
        term2 *= CrosssecFusion(giveq() * Potential_Phi(r));

        NeutronFlux = giveTransparency() * ( term1 + term2);
    }
    else
    {
        double (*FunctPtr)(double, double);
        FunctPtr = &NeutronsIonFluxOutwards_Inte1;

        term1 = NIntegration_2(FunctPtr, r, r, giveAnodeRadius());

        term2 =  pow(giveAnodeRadius() * f(0)/r,2) * EdgeIonFlux / f(r);
        term2 *= CrosssecFusion(giveq() * Potential_Phi(r));

        NeutronFlux = pow(giveTransparency(),2) * ( term1 + term2);
    }


    return NeutronFlux;
}
