/**
    The next function compute the chance that a particle survise to that radius
*/

#include <stdio.h>
#include <math.h>

#include "constants.h"
#include "MathFunctions.h"
#include "CrossSections.h"
#include "PotentialFunctions.h"

#include "SurvivalFunctions.h"

double f_inte(double r)
{
    return CrosssecCX(ParticleEnergy1(r));
}

double f(double r)
{
    double (*f_intePtr)(double);
    f_intePtr = &f_inte;

    double sum = NIntegration(*f_intePtr, r, giveAnodeRadius());

    sum *= ngas;
    sum = exp(-1 * sum);

    return sum;
}

double A(double r)
{
    return ngas * CrosssecTot(ParticleEnergy1(r)) * pow(giveAnodeRadius()/r,2) * (f(r) + pow(giveTransparency() * f(0),2)/f(r));
}

double g_inte(double dr, double ddr)
{
    return CrosssecCX(ParticleEnergy2(ddr,dr));
}

double g(double r, double dr)
{
    if ( dr < r)
    {
        printf("G function error: dr < r \n");
        return NAN;
    }

    double (*g_intePtr)(double, double);
    g_intePtr = &g_inte;

    double sum = NIntegration_2(*g_intePtr, dr, r, dr);

    sum *= ngas;
    sum = exp(-sum);

    return sum;
}

double kernel(double r, double dr)
{
    double tmp;
    if ( r <= dr )
    {
        tmp  = ngas * CrosssecTot(ParticleEnergy2(r,dr));
        tmp *= pow(dr/r,2) * ((g(r,dr) + (pow(giveTransparency()*g(0,dr),2)/g(r,dr)))/(1.0-pow(giveTransparency()*g(0,dr),2)));
    }
    else
        tmp = 0;

    return tmp;
}

void CalculateTables(void)
{
    double step = (giveAnodeRadius() - giveCathodeRadius())/(N_TABLE-1);

	fprintf(stdout,"\n# - Filling tables ----------------------------------- #\n");

    printf("#");
    for ( int i = 0; i < N_TABLE; i++)
    {
        for ( int j = 0; j < N_TABLE; j++)
        {
            Table->K[i][j] = kernel(i*step+giveCathodeRadius(),j*step+giveCathodeRadius());
        }

        if ( i *(N_TABLE / 50) % 20 == 0)
            printf(".");

        Table->A[i] = A(i*step+giveCathodeRadius());

        Table->R[i] = i*step+giveCathodeRadius();
    }

    printf("\n");

    return;
}
