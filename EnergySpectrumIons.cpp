
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

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

//    printf("E = %E, r = %E, dr = %E\n", E, r, dr);

    if ( r > dr || dr > giveAnodeRadius())
        return -1;

    double (*PhiPtr)(double);
    PhiPtr = &Potential_Phi;

    double dPhi_dr = differentiat(*PhiPtr, dr);

//    printf("\t S = %E, dPhi/dr = %E\n", dr, dPhi_dr);

    result = pow(dr/r,2) * interpolation(dr) * abs(dPhi_dr) * g(r,dr)/( 1 - pow(giveTransparency()*g(0,dr),2) );

    // a delta term should be included here
    return result;
}

double f_plus(double r, double E)
{
    double result;
    double dr = r_shell(r, E);

//    printf("E = %E, r = %E, dr = %E\n", E, r, dr);

    if ( r > dr || dr > giveAnodeRadius())
        return -1;

    double (*PhiPtr)(double);
    PhiPtr = &Potential_Phi;

    double dPhi_dr = differentiat(*PhiPtr, dr);

//    printf("\t S = %E, dPhi/dr = %E\n", dr, dPhi_dr);

    result = pow(dr/r,2) * interpolation(dr) * abs(dPhi_dr) * pow(giveTransparency() * g(r,dr),2)/( 1 - pow(giveTransparency()*g(0,dr),2) );
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
