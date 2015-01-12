

#include "constants.h"
#include "MathFunctions.h"
#include "PotentialFunctions.h"
#include "SurvivalFunctions.h"

#include "IonSpectrum.h"


double IonSpectrumInwards(double r, double E)
{
    double f;
    double dr = r_shell(r, E);

    double (*PhiPtr)(double);
    PhiPtr = &Potential_Phi;

    f  = 1/giveq();
    f *= interpolation(dr) / differentiat(*PhiPtr, dr);
    f *= g(r, dr) / ( 1 - pow(giveTransparency() * g(0,dr),2) );

    bool delta = DELTA(E - ParticleEnergy1(r));

    if ( delta )
        f += pow(giveCathodeRadius()/r,) * EdgeIonFlux * f(r);

    if ( r < giveCathodeRadius() )
        f *= giveTransparency();

    return f;
}

double IonSpectrumOutwards(double r, double E)
{
    double f;
    double dr = r_shell(r, E);

    double (*PhiPtr)(double);
    PhiPtr = &Potential_Phi;

    f  = 1/giveq();
    f *= interpolation(dr) / differentiat(*PhiPtr, dr);



    f *= g(r, dr) / ( 1 - pow(giveTransparency() * g(0,dr),2) );
    f *= 1 / g(r, dr);

    bool delta = DELTA(E - ParticleEnergy1(r));

    if ( delta )
    {
        f +=  pow(giveCathodeRadius() * f(0)/r,) * EdgeIonFlux / f(r);
    }

    f *= giveTransparency();
    if ( r > giveCathodeRadius() )
        f *= giveTransparency();

    return f;
}
