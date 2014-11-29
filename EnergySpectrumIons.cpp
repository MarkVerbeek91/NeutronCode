
//#include <stdio.h>
#include <math.h>

#include "constants.h"
#include "PotentialFunctions.h"
#include "SurvivalFunctions.h"
#include "MathFunctions.h"
#include "EnergySpectrumIons.h"

// inward ions.

double f_min(double r, double E)
{
    double result;
    double dr = r_shell(r, E);

//    printf("\t E = %E, r = %E, dr = %E\n", E, r, dr);

    if ( r > dr)
        return -1;

    double (*PhiPtr)(double);
    PhiPtr = &Potential_Phi;

    double dPhi_dr = differentiat(*PhiPtr, dr);

    result = (1/giveq()) * pow(dr/r,2) * interpolation(dr) * dPhi_dr * g(r,dr)/( 1 - pow(giveTransparency()*g(0,dr),2) );
    result *= 1/g(r,dr);

    // a delta term should be included here
    return result;
}

double r_shell(double r, double E)
{
    double radius;

    radius =  -giveVoltage() * giveAnodeRadius() * giveCathodeRadius() * r;
    double teller = (giveAnodeRadius() * giveCathodeRadius() * -giveVoltage() + (giveCathodeRadius() - giveAnodeRadius()) * E * r );
    radius /= teller;

    return radius;
}
