/**
 * \brief The functions in this file calculate the differend currents that
 * exist inside the cathode.

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
 * \brief Calculation of the Secondary Ion Induced Electron Emision coeffiecent
 * \parameter E for energy in eV
 * \return gamma, dimentionless
 */
double SIIEE(double E)
{
    double gamma;
    gamma = 1.5 * pow((1.15*(E/97861)),-0.667)*(1-exp(-1.8*pow(E/97891,1.2)));
    return gamma;
}

double I_c1(void)
{
    double current;
	// C * m^2
    current  = 4 * 3.141529 * Q_ELECTRON * (1.0 - fusor->Tc) * pow(fusor->b,2);
    current *= f(fusor->a) + (fusor->Tc * pow(f(0.0),2))/f(fusor->a);
    current *= 1 + SIIEE(-fusor->V0);

    return current;
}

double I_c2inte(double dr)
{
    double  fac  = interpolation(Table->S, dr) / ( 1 - pow(fusor->Tc*g(0,dr),2));
            fac *= (g(fusor->a,dr) + fusor->Tc * pow(g(0,dr),2)/g(fusor->a,dr));
            fac *= (1 + SIIEE(ParticleEnergy2(0,dr))) * pow(dr,2);

    return fac;
}

double I_c2(void)
{
    double term;
    double (*FuncPtr)(double);
    FuncPtr = &I_c2inte;

    term = NIntegration(*FuncPtr, fusor->a, fusor->b);

    term *= 4 * 3.141529 * Q_ELECTRON * (1 - fusor->Tc) * term;

    return term;
}

double I_c3(void)
{
    double current;
    double Emax = ParticleEnergy1(fusor->a);

    current  = 4 * 3.141529 * Q_ELECTRON * fusor->Tc * (CrosssecTot(Emax)/CrosssecCX(Emax));
    current *= pow(fusor->b,2) * f(fusor->a) * ( 1.0 - exp( -2 * fusor->ngas* CrosssecCX(Emax) * fusor->a));

    return current;
}

double I_c4inte(double dr)
{
    double fac;
	double E = ParticleEnergy2(fusor->a,dr);

    fac  = CrosssecTot(E)/CrosssecCX(E);
    fac *= interpolation(Table->S, dr) * g(fusor->a,dr) / ( 1.0 - pow(fusor->Tc * g(0,dr),2));
    fac *= 1 - exp( -2 * fusor->ngas* CrosssecCX(E) * fusor->a);
	fac *= pow(dr,2);

    return fac;
}

double I_c4(void)
{
    double term;
    double (*FuncPtr)(double);
    FuncPtr = &I_c4inte;

    term = NIntegration(*FuncPtr, fusor->a, fusor->b);

    term *= 4 * 3.141529 * Q_ELECTRON * term;

    return term;
}


