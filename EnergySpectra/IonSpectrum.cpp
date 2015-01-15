
#include <math.h>
#include <stdlib.h>

#include "constants.h"
#include "MathFunctions.h"
#include "PotentialFunctions.h"
#include "SurvivalFunctions.h"
#include "Crosssections.h"

#include "IonSpectrum.h"

// equations 28 and 30;
double IonSpectrumInwards(double r, double E)
{
    double flux, term2 = 0;
    double dr = r_shell(r, E);

    if ( dr < giveCathodeRadius() | dr > giveAnodeRadius())
        return -1;

    double (*PhiPtr)(double);
    PhiPtr = &Potential_Phi;

    flux  = 1/giveq();
    flux *= pow(dr/r,2);
    flux *= interpolation(dr) / abs(differentiat(*PhiPtr, dr));
    flux *= 1 / ( 1 - pow(giveTransparency() * g(0,dr),2) );

    bool delta = DELTA(E - ParticleEnergy1(r));

    if ( r < giveCathodeRadius() )
    {
        flux *= giveTransparency() * g(giveCathodeRadius(),dr) * exp(ngas * CrosssecCX(ParticleEnergy2(giveCathodeRadius(),dr)) * (r - giveCathodeRadius()));

        if ( delta )
        {
            term2  = pow(giveCathodeRadius()/r,2) * EdgeIonFlux ;
            term2 *= f(giveCathodeRadius()) * exp(ngas * CrosssecCX(ParticleEnergy1(giveCathodeRadius())) * (r - giveCathodeRadius()));
        }
    }
    else
    {
        flux  *= g(r,dr);

        if ( delta )
        {
            term2  = pow(giveCathodeRadius()/r,2) * EdgeIonFlux ;
            term2 *= f(r);
        }
    }

    return flux + term2;
}

// equations 29 and 31
double IonSpectrumOutwards(double r, double E)
{
    double flux, term2 = 0;
    double dr = r_shell(r, E);

    if ( dr < giveCathodeRadius() | dr > giveAnodeRadius())
        return -1;

    double (*PhiPtr)(double);
    PhiPtr = &Potential_Phi;

    flux  = 1/giveq();
    flux *= pow(dr/r,2);
    flux *= interpolation(dr) / abs(differentiat(*PhiPtr, dr));


    flux *= giveTransparency();
    flux *= pow(g(0, dr),2) / ( 1 - pow(giveTransparency() * g(0,dr),2) );

    bool delta = DELTA(E - ParticleEnergy1(r));

    if ( r > giveCathodeRadius() )
    {
        flux *= giveTransparency();

        if ( delta )
        {
            flux *= 1 / g(r, dr);
            flux +=  pow(giveCathodeRadius() * f(0)/r,2) * EdgeIonFlux / f(r);
        }
    }
    else
    {
        flux *= 1 / (g(giveCathodeRadius(), dr) * exp(ngas * CrosssecCX(ParticleEnergy2(giveCathodeRadius(),dr)) * ( r - giveCathodeRadius())));

        if ( delta )
        {
            term2  = pow(giveCathodeRadius()/r,2) * EdgeIonFlux ;
            term2 /= f(giveCathodeRadius()) * exp(ngas * CrosssecCX(ParticleEnergy1(giveCathodeRadius())) * (r - giveCathodeRadius()));
        }
    }

    return flux + term2;
}
