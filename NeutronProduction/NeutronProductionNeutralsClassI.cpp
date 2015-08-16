#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "constants.h"
#include "MathFunctions.h"
#include "PotentialFunctions.h"
#include "SurvivalFunctions.h"
#include "CrossSections.h"

#include "NeutronProductionNeutralsClassI.h"

/**
        NEUTRONS FROM NEUTRALS FROM CLASS I IONS GOING INWARDS
*/
double NeutronsNeutralsClassIFluxInwards_Inte1(double r, double dr)
{
    double E, term1;

    if ( r < dr )
    {
        printf("NeutronsNeutralsClassIFluxInwards_Inte1 error: r < dr \n");
        return NAN;
    }

    E = giveq() * ( Potential_Phi(r) - Potential_Phi(dr));

    term1  = CrosssecFusion(E);
    term1 *= pow(giveAnodeRadius()/r,2);
    term1 *= ngas * CrosssecCX(E);
    term1 *= EdgeIonFlux;
    term1 *= f(dr);

    return term1;
}

double NeutronsNeutralsClassIFluxInwards_Inte2(double r, double dr)
{
    double E, term1;

    E = giveq() * ( Potential_Phi(r) - Potential_Phi(dr));

    term1  = CrosssecFusion(E);
    term1 *= pow(giveAnodeRadius()/r,2);
    term1 *= ngas * CrosssecCX(E);
    term1 *= EdgeIonFlux;
    term1 *= f(dr);

    return term1;
}

double NeutronsNeutralsClassIFluxInwards(double r)
{
    double NeutronFlux;
    double term1 = 0, term2 = 0;

    if ( giveCathodeRadius() < r )
    {
        double (*FunctPtr)(double, double);
        FunctPtr = &NeutronsNeutralsClassIFluxInwards_Inte1;

        term1 = NIntegration_2(FunctPtr, r, r, giveAnodeRadius());

        NeutronFlux = term1;
    }
    else
    {
        double (*FunctPtr)(double, double);
        FunctPtr = &NeutronsNeutralsClassIFluxInwards_Inte2;

        term1 = NIntegration_2(FunctPtr, r, r, giveAnodeRadius());

        term2  = pow(giveAnodeRadius()/r,2) * EdgeIonFlux ;
        term2 *= f(giveCathodeRadius()) * ( 1 - exp(-2 * ngas * CrosssecCX(ParticleEnergy1(giveCathodeRadius())) * giveCathodeRadius()));
        term2 *= CrosssecFusion(giveq() * Potential_Phi(r));

        NeutronFlux = giveTransparency() * (term1 + term2);
    }

    return NeutronFlux;
}

/**
        NEUTRONS FROM NEUTRALS FROM CLASS I IONS GOING INWARDS
*/
// inside the cathode
double NeutronsNeutralsClassIFluxOutwards_Inte1(double r, double dr)
{
    double E, term1;

    E = giveq() * ( Potential_Phi(r) - Potential_Phi(dr));

    term1  = CrosssecFusion(E);
    term1 *= pow(giveAnodeRadius()/r,2);
    term1 *= ngas * CrosssecCX(E);
    term1 *= EdgeIonFlux;
    term1 *= f(dr);

    return term1;
}

// outside the cathode
double NeutronsNeutralsClassIFluxOutwards_Inte2(double r, double dr)
{
    double E, term1;

    E = giveq() * ( Potential_Phi(r) - Potential_Phi(dr));

    term1 *= f(dr);

    if ( dr < r)
        term1 += pow(f(0),2)/f(dr);

    term1  = CrosssecFusion(E);
    term1 *= pow(giveAnodeRadius()/r,2);
    term1 *= ngas * CrosssecCX(E);
    term1 *= EdgeIonFlux;

    return term1;
}

double NeutronsNeutralsClassIFluxOutwards(double r)
{
    double NeutronFlux;
    double term1 = 0, term2 = 0;

    if ( giveCathodeRadius() < r )
    {
        double (*FunctPtr)(double, double);
        FunctPtr = &NeutronsNeutralsClassIFluxOutwards_Inte1;

        term1 = NIntegration_2(FunctPtr, r, r, giveAnodeRadius());

        term2  = pow(giveAnodeRadius()/r,2) * EdgeIonFlux ;
        term2 *= f(giveCathodeRadius()) * ( 1 - exp(2 * ngas * CrosssecCX(ParticleEnergy1(giveCathodeRadius())) * (r - giveCathodeRadius())));
        term2 *= CrosssecFusion(giveq() * Potential_Phi(r));

        NeutronFlux = giveTransparency() * (term1 + term2);
    }
    else
    {
        double (*FunctPtr)(double, double);
        FunctPtr = &NeutronsNeutralsClassIFluxOutwards_Inte1;

        term1 = NIntegration_2(FunctPtr, r, r, giveAnodeRadius());

        term2  = pow(giveAnodeRadius()/r,2) * EdgeIonFlux ;
        term2 *= f(giveCathodeRadius()) * ( 1 - exp(-2 * ngas * CrosssecCX(ParticleEnergy1(giveCathodeRadius())) * giveCathodeRadius()));
        term2 *= CrosssecFusion(giveq() * Potential_Phi(r));

        NeutronFlux = giveTransparency() * (term1 + term2);
    }

    return NeutronFlux;
}
