#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "constants.h"
#include "PotentialFunctions.h"
#include "SurvivalFunctions.h"
#include "MathFunctions.h"
#include "EnergySpectrumIons.h"

// inward ions.

double f_min(double r, double E)
{
    double result;
    double dr = Potential_Phi_Inv(E);

    if ( r > dr || dr > giveAnodeRadius())
    {
        printf("f_min function error: r > dr or dr > a\n");
        return NAN;
    }


    double (*PhiPtr)(double);
    PhiPtr = &Potential_Phi;

    double dPhi_dr = differentiat(*PhiPtr, dr);

    result = pow(dr/r,2) * (interpolation(Table.S, dr) / abs(dPhi_dr)) * g(r,dr)/( 1 - pow(giveTransparency()*g(0,dr),2) );

    // a delta term should be included here
    double delta = 0;
    if ( E - Potential_Phi(r) < 0.001 )
        delta = 1;

    result += pow(giveAnodeRadius()/r,2) * EdgeIonFlux * f(r) * delta;

    if ( r < giveCathodeRadius())
        result *= giveTransparency();

    return result;
}

double f_plus(double r, double E)
{
    double result;
    double dr = r_shell2(r, E);

//    printf("E = %E, r = %E, dr = %E\n", E, r, dr);

    if ( r > dr || dr > giveAnodeRadius())
    {
        printf("f_plus function error: r > dr or dr > Anode Radius\n");
        return NAN;
    }


    double (*PhiPtr)(double);
    PhiPtr = &Potential_Phi;

    double dPhi_dr = differentiat(*PhiPtr, dr);

    if ( abs(dPhi_dr) < 0.0001 )
    {
        printf("f_plus error: dPhi to small\n");
        return NAN;
    }


    result = pow(dr/r,2) * (interpolation(Table.S, dr) / abs(dPhi_dr)) * pow(giveTransparency() * g(0,dr),2)/( 1 - pow(giveTransparency()*g(0,dr),2) );
    result *= 1/g(r,dr);

    double delta = 0;
    if ( E - abs(Potential_Phi(r)) < 0.001 )
        delta = 1;

    result += pow(giveAnodeRadius() * giveTransparency() * f(0) / r,2) * EdgeIonFlux * delta / f(r);

    if ( r < giveCathodeRadius())
        result *= giveTransparency();

    return result;
}

double r_shell2(double r, double E)
{
    double radius;

    radius =  -giveVoltage() * giveAnodeRadius() * giveCathodeRadius() * r;
    double teller = (giveAnodeRadius() * giveCathodeRadius() * -giveVoltage() + (giveCathodeRadius() - giveAnodeRadius()) * E * r );
    radius /= teller;

    return radius;
}
