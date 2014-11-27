
#include <math.h>

#include "constants.h"
#include "SurvivalFunctions.h"
#include "MathFunctions.h"
#include "EnergySpectrumIons.h"

// inward ions.

double f_min(double r, double E)
{
    double result;

    if ( r > r_shell(r, E))
        return -1;

    result = (1/giveq()) * pow(r_shell(r,E)/r,2) * interpolation(r_shell()) * dPhi_dr(r) * g(r,dr)/( 1 - pow(giveTransparency()*g(0,r_shell()),2) )
    result *= 1/g(r,r_shell());

    // a delta term should be included here
    return result;
}

double r_shell(double r, double E)
{
    double radius;

    radius =  pow(giveVoltage(),2) * giveq() * giveAnodeRadius() * giveCathodeRadius() * r;
    radius /= E * giveCathodeRadius() * r - E * giveAnodeRadius() * r + giveq() * giveAnodeRadius() * giveCathodeRadius() * pow(giveVoltage(),2);

    return radius;
}

double dPhi_dr(double r)
{
    double phi;

    phi =

    return phi;
}
