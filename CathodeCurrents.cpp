/** the next section will calculate currents and such:

*/

#include <stdio.h>
#include <math.h>

#include "SurvivalFunctions.h"
#include "CrossSections.hpp"
#include "PotentialFunctions.h"
#include "MathFunctions.h"
#include "constants.hpp"
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

    current = 4 * 3.141529 * pow(fusor.b,2) * q * (1.0 - fusor.Tc) * (f(fusor.a) + (fusor.Tc * pow(f(0.0),2))/f(fusor.a))*(1 + SIIEE(-fusor.V0));

    return current;
}

double I_c2inte(double dr)
{
    double  fac  = interpolation(dr) / ( 1 - pow(fusor.Tc,2)*g(0,dr));
            fac *= (g(fusor.a,dr) + fusor.Tc * g(0,dr)/g(fusor.a,dr));
            fac *= (1 + SIIEE(ParticleEnergy2(0,dr))) * pow(dr,2);

    return fac;
}

double I_c2(void)
{
    double current = 0;
    double sum = 0;

    double (*I_c2intePtr)(double);
    I_c2intePtr = &I_c2inte;

    sum = NIntegration(*I_c2intePtr, fusor.a, fusor.b);

    current = 4 * 3.141529 * q * (1 - fusor.Tc) * sum;

    return current;
}

double I_c3(void)
{
    double current;

    current  = 4 * 3.141529 * pow(fusor.b,2) * q * (CrosssecTot(-fusor.V0)/CrosssecCX(-fusor.V0));
    current *= fusor.Tc * f(fusor.a) * ( 1 - exp( -2 * ngas * fusor.a * CrosssecCX(-fusor.V0)));

    return current;
}

double I_c4inte(double dr)
{
    double fac;

    fac  = pow(dr,2);
    fac *= CrosssecTot(ParticleEnergy2(fusor.a,dr))/CrosssecCX(ParticleEnergy2(fusor.a,dr));
    fac *= (interpolation(dr) * g(fusor.a,dr))/( 1 - pow(fusor.Tc,2) * g(0,dr));
    fac *= ( 1 - exp( -2 * ngas * fusor.a * CrosssecCX(ParticleEnergy2(fusor.a,dr))));

    return fac;
}

double I_c4(void)
{
    double current = 0;
    double sum = 0;

    double (*I_c4intePtr)(double);
    I_c4intePtr = &I_c4inte;

    sum = NIntegration(*I_c4intePtr, fusor.a, fusor.b);

    current = 4 * 3.141529 * pow(fusor.b,2) * q * sum;

    return current;
}


