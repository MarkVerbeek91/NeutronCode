

#include "constants.h"
#include "SurvivalFunctions.h"



// inward ions.

double f_min(double r, double E)
{
    double result;

    result = (1/giveq()) * pow(r,2) * interpolation(r) * dPhi_dr(r) * g(r,dr)/( 1 - pow(giveTransparency()*g(0,dr),2) )


    return result;
}



