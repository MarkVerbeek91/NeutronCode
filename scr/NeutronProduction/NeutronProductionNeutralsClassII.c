#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "constants.h"
#include "MathFunctions.h"
#include "PotentialFunctions.h"
#include "SurvivalFunctions.h"
#include "CrossSections.h"

#include "NeutronProductionNeutralsClassII.h"

/**
        NEUTRONS FROM NEUTRALS FROM CLASS II IONS GOING INWARDS
*/

// intregral of neutral from claas II ions _outside_ the cathode _inwards_
double NeutralsClassIISpectrumInwards_Inte1_Inte(double E, double r, double ddr)
{
    double dr = Potential_Phi_Inv(Potential_Phi(ddr) - E/giveq());

    if ( dr > r )
        return 0;

    double integrant;

    integrant  = interpolation(Table->S, ddr);
    integrant *= g(dr, ddr);
    integrant /= 1 - pow(giveTransparency() * g(0,ddr),2);
    integrant *= pow(ddr,2);

    return integrant;
}

// outside the cathode
double NeutronsNeutralsClassIIFluxInwards_Inte1(double r, double dr)
{
    double E, term1;

    if ( r < dr )
    {
        printf("NeutronsNeutralsClassIIFluxInwards_Inte1 error: r < dr\n");
        return NAN;
    }

    E = ParticleEnergy1(dr);

    term1  = CrosssecFusion(E);
    term1 /= pow(r,2);
    term1 *= CrosssecCX(E);

    double (*FunctPtr)(double, double, double);
    FunctPtr = &NeutralsClassIISpectrumInwards_Inte1_Inte;

    term1 *= NIntegration_3(*FunctPtr, E, r, r, giveAnodeRadius());

    return term1;
}

// intregral of neutral from claas II ions _inside_ the cathode _inwards_
double NeutralsClassIISpectrumInwards_Inte2_Inte(double E, double ddr)
{
    double dr = Potential_Phi_Inv(Potential_Phi(ddr) - E/giveq());

    if ( giveCathodeRadius() > dr )
    {
        return 0;
    }

    double integrant;

    integrant  = interpolation(Table->S, ddr);
    integrant *= g(dr, ddr);
    integrant /= 1 - pow(giveTransparency() * g(0,ddr),2);
    integrant *= pow(ddr,2);

    return integrant;
}

// inside the cathode
double NeutronsNeutralsClassIIFluxInwards_Inte2(double r, double dr)
{
    double E, term1;

    if ( r > giveCathodeRadius() )
    {
        printf("NeutralsClassIISpectrumInwards_Inte2 error: r > d\n");
        return NAN;
    }


    E = ParticleEnergy1(dr);

    term1  = CrosssecFusion(E);
    term1 /= pow(r,2);
    term1 *= CrosssecCX(E);

    double (*FunctPtr)(double, double);
    FunctPtr = &NeutralsClassIISpectrumInwards_Inte2_Inte;

    term1 *= NIntegration_2(*FunctPtr, E, giveCathodeRadius(), giveAnodeRadius());


    return term1;
}

double NeutronsNeutralsClassIIFluxInwards(double r)
{
    double NeutronFlux;
    double term1 = 0, term2 = 0;

    if ( giveCathodeRadius() < r )
    {
        double (*FunctPtr)(double, double);
        FunctPtr = &NeutronsNeutralsClassIIFluxInwards_Inte1;

        term1 = NIntegration_2(FunctPtr, r, r, giveAnodeRadius());

        NeutronFlux = term1;
    }
    else
    {
        double (*FunctPtr)(double, double);
        FunctPtr = &NeutronsNeutralsClassIIFluxInwards_Inte2;

        term1 = NIntegration_2(FunctPtr, r, r, giveAnodeRadius());

        term2  = pow(giveAnodeRadius()/r,2) * EdgeIonFlux ;
        term2 *= f(giveCathodeRadius()) * ( 1 - exp(-2 * ngas * CrosssecCX(ParticleEnergy1(giveCathodeRadius())) * giveCathodeRadius()));
        term2 *= CrosssecFusion(ParticleEnergy1(r));

        NeutronFlux = giveTransparency() * (term1 + term2);
    }

    return ngas * NeutronFlux;
}

/**
        NEUTRONS FROM NEUTRALS FROM CLASS II IONS GOING OUTWARDS
*/
// intregral of neutrals from claas II ions _outside_ the cathode _outward_
double NeutralsClassIISpectrumOutwards_Inte1_Int(double E, double ddr)
{
    double dr = Potential_Phi_Inv(Potential_Phi(ddr) - E/giveq());

    if ( dr < giveCathodeRadius() )
    {
        printf("NeutralsClassIISpectrumOutwards_Inte1_Int error: dr < b\n");
        return NAN;
    }

    double integrant;

    integrant  = interpolation(Table->S, ddr);
    integrant *= g(giveCathodeRadius(), ddr);
    integrant /= 1 - pow(giveTransparency() * g(0,ddr),2);
    integrant *= pow(ddr,2);

    return integrant;
}

// inside the cathode
double NeutronsNeutralsClassIIFluxOutwards_Inte1(double r, double dr)
{
    double E, term1;

    E = ParticleEnergy2(r, dr);

    term1  = pow(giveAnodeRadius()/r,2);
    term1 *= ngas * CrosssecCX(E);
    term1 *= EdgeIonFlux;
    term1 *= f(dr);

    double (*FunctPtr)(double, double);
    FunctPtr = &NeutralsClassIISpectrumOutwards_Inte1_Int;

    term1 *= NIntegration_2(*FunctPtr, E, giveCathodeRadius(), giveAnodeRadius());

    return term1;
}

// intregral of neutral from claas II ions _inside_ the cathode _outwards_
double NeutralsClassIISpectrumOutwards_Inte2_Int(double E, double r, double ddr)
{
    double dr = Potential_Phi_Inv(Potential_Phi(ddr) - E/giveq());

    double integrant = 0;
    if ( giveCathodeRadius() < dr )
        integrant = g(dr, ddr);

    if ( dr < r )
        integrant += pow(g(0,ddr),2) / g(dr, ddr);

    if ( integrant == 0)
        return 0;

    integrant  = interpolation(Table->S, ddr);
    integrant /= 1 - pow(giveTransparency() * g(0,ddr),2);
    integrant *= pow(ddr,2);

    return integrant;
}


// outside the cathode
double NeutronsNeutralsClassIIFluxOutwards_Inte2(double r, double dr)
{
    double E, term1;

    E = ParticleEnergy2(r, dr);

    term1 = f(dr);

    if ( dr < r)
        term1 += pow(f(0),2)/f(dr);

    term1 *= pow(giveAnodeRadius()/r,2);
    term1 *= ngas * CrosssecCX(E);
    term1 *= EdgeIonFlux;

    double (*FunctPtr)(double, double, double);
    FunctPtr = &NeutralsClassIISpectrumOutwards_Inte2_Int;

    term1 *= NIntegration_3(FunctPtr, E, r, giveCathodeRadius(), giveAnodeRadius());

    return term1;
}

double NeutronsNeutralsClassIIFluxOutwards(double r)
{
    double NeutronFlux;
    double term1 = 0, term2 = 0;

    if ( giveCathodeRadius() < r )
    {
        double (*FunctPtr)(double, double);
        FunctPtr = &NeutronsNeutralsClassIIFluxOutwards_Inte2;

        term1 = NIntegration_2(FunctPtr, r, r, giveAnodeRadius());

        term2  = 0; // TODO

        NeutronFlux = giveTransparency() * (term1 + term2);
    }
    else
    {
        double (*FunctPtr)(double, double);
        FunctPtr = &NeutronsNeutralsClassIIFluxOutwards_Inte1;

        term1 = NIntegration_2(FunctPtr, r, r, giveAnodeRadius());

        term2  = 0; // TODO

        NeutronFlux = giveTransparency() * (term1 + term2);
    }

    return NeutronFlux;
}
