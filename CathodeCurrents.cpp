/** the next section will calculate currents and such:

*/

#include <stdio.h>
#include <math.h>

#include "constants.h"
#include "SurvivalFunctions.h"
#include "CrossSections.h"
#include "PotentialFunctions.h"
#include "MathFunctions.h"

#include "CathodeCurrents.h"

/**
    In E in eV
*/
double SIIEE(double energy)
{
    double tmp;
    tmp = 1.5 * pow((1.15*(energy/97861)),-0.667)*(1-exp(-1.8*pow(energy/97891,1.2)));
    return tmp;
}

double I_c1(void)
{
    double current;

    current = 4 * 3.141529 * pow(giveAnodeRadius(),2) * giveq() * (1.0 - giveTransparency()) * (f(giveCathodeRadius()) + (giveTransparency() * pow(f(0.0),2))/f(giveCathodeRadius()))*(1 + SIIEE(-giveVoltage()));

    return current;
}

double I_c2inte(double dr)
{
    double  fac  = interpolation(dr) / ( 1 - pow(giveTransparency(),2)*g(0,dr));
            fac *= (g(giveCathodeRadius(),dr) + giveTransparency() * g(0,dr)/g(giveCathodeRadius(),dr));
            fac *= (1 + SIIEE(ParticleEnergy2(0,dr))) * pow(dr,2);

    return fac;
}

double I_c2(void)
{
    double current = 0;
    double sum = 0;

    double (*I_c2intePtr)(double);
    I_c2intePtr = &I_c2inte;

    sum = NIntegration(*I_c2intePtr, giveCathodeRadius(), giveAnodeRadius());

    current = 4 * 3.141529 * q * (1 - giveTransparency()) * sum;

    return current;
}

double I_c3(void)
{
    double current;

    current  = 4 * 3.141529 * pow(giveAnodeRadius(),2) * q * (CrosssecTot(-giveVoltage())/CrosssecCX(-giveVoltage()));
    current *= giveTransparency() * f(giveCathodeRadius()) * ( 1 - exp( -2 * ngas * giveCathodeRadius() * CrosssecCX(-giveVoltage())));

    return current;
}

double I_c4inte(double dr)
{
    double fac;

    fac  = pow(dr,2);
    fac *= CrosssecTot(ParticleEnergy2(giveCathodeRadius(),dr))/CrosssecCX(ParticleEnergy2(giveCathodeRadius(),dr));
    fac *= (interpolation(dr) * g(giveCathodeRadius(),dr))/( 1 - pow(giveTransparency(),2) * g(0,dr));
    fac *= ( 1 - exp( -2 * ngas * giveCathodeRadius() * CrosssecCX(ParticleEnergy2(giveCathodeRadius(),dr))));

    return fac;
}

double I_c4(void)
{
    double current = 0;
    double sum = 0;

    double (*I_c4intePtr)(double);
    I_c4intePtr = &I_c4inte;

    sum = NIntegration(*I_c4intePtr, giveCathodeRadius(), giveAnodeRadius());

    current = 4 * 3.141529 * pow(giveAnodeRadius(),2) * q * sum;

    return current;
}


