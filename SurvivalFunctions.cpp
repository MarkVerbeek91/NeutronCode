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

double Intergrant(double r)
{
    double tmp;

    tmp = ngas * CrosssecCX(ParticleEnergy1(r));

    return tmp;

}

double gamma(double r)
{
    double Gamma = -1;
    Gamma = pow(giveAnodeRadius()/r,2) * (f(r) + pow(giveTransparency() * f(0),2)/f(r));
    return Gamma;
}

double A(double r)
{
    double tmp = -1;
    tmp = ngas * CrosssecTot(ParticleEnergy1(r)) * gamma(r);

    if (isnan(tmp))
    {
        printf("went to the end\n");
        tmp = 0;
    }
    return tmp;
}

double g_inte(double dr, double ddr)
{
    return CrosssecCX(ParticleEnergy2(ddr,dr));
}

double g(double r, double dr)
{
    if ( dr < r)
    {
        printf("error: dr < r \n");
        return -2;
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
    if ( r < r1 )
    {
        tmp = pow(dr/r,2) * ((g(r,dr) + (pow(giveTransparency()*g(0,dr),2)/g(r,dr)))/(1.0-pow(giveTransparency()*g(0,dr),2)));
        tmp = tmp * ngas * CrosssecTot(ParticleEnergy2(r,dr));
    }
    else
        tmp = 0;

    return tmp;
}

void kernel_to_table(void)
{
    double step = (giveAnodeRadius() - giveCathodeRadius())/(N_TABLE-1);

    printf("#");
    for ( int i = 0; i < N_TABLE; i++)
    {
        for ( int j = 0; j < N_TABLE; j++)
        {
            Table.K[i][j] = kernel(i*step+giveCathodeRadius(),j*step+giveCathodeRadius());
        }

        if ( i *(N_TABLE / 50) % 20 == 0)
            printf(".");

        Table.A[i] = A(i*step+giveCathodeRadius());

        Table.R[i] = i*step+giveCathodeRadius();
    }

    // to do, change this.
//    Table.A[N_TABLE-1] = 43.9944;

    printf("\n");

    return;
}
