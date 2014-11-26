

#include "constants.h"
#include "SurvivalFunctions.h"
-


// inward ions.

double f_min(double r, double E)
{
    double result;

    result = (1/giveq()) * pow(r,2) * interpolation(r) * dPhi_dr(r) * g(r,dr)/( 1 - pow(giveTransparency()*g(0,dr),2) )


    return result;
}

double r_shell(double r, double E)
{
    double radius;

    radius =  pow(giveVoltage(),2) * giveq() * giveAnodeRadius() * giveCathodeRadius() * r;
    radius /= E * giveCathodeRadius() * r - E * giveAnodeRadius() * r + giveq() * giveAnodeRadius() * giveCathodeRadius() * pow(giveVoltage(),2);

    return radius;
}
