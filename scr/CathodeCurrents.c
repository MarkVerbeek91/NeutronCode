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
    current  = 4 * 3.141529 * Q_ELECTRON * (1.0 - giveTransparency()) * pow(giveAnodeRadius(),2);
    current *= f(giveCathodeRadius()) + (giveTransparency() * pow(f(0.0),2))/f(giveCathodeRadius());
    current *= 1 + SIIEE(-giveVoltage());

    return current;
}

double I_c2inte(double dr)
{
    double  fac  = interpolation(Table->S, dr) / ( 1 - pow(giveTransparency()*g(0,dr),2));
            fac *= (g(giveCathodeRadius(),dr) + giveTransparency() * pow(g(0,dr),2)/g(giveCathodeRadius(),dr));
            fac *= (1 + SIIEE(ParticleEnergy2(0,dr))) * pow(dr,2);

    return fac;
}

double I_c2(void)
{
    double term;
    double (*FuncPtr)(double);
    FuncPtr = &I_c2inte;

    term = NIntegration(*FuncPtr, giveCathodeRadius(), giveAnodeRadius());

    term *= 4 * 3.141529 * Q_ELECTRON * (1 - giveTransparency()) * term;

    return term;
}

double I_c3(void)
{
    double current;
    double Emax = ParticleEnergy1(giveCathodeRadius());

    current  = 4 * 3.141529 * Q_ELECTRON * giveTransparency() * (CrosssecTot(Emax)/CrosssecCX(Emax));
    current *= pow(giveAnodeRadius(),2) * f(giveCathodeRadius()) * ( 1.0 - exp( -2 * ngas * CrosssecCX(Emax) * giveCathodeRadius()));

    return current;
}

double I_c4inte(double dr)
{
    double fac;
	double E = ParticleEnergy2(giveCathodeRadius(),dr);
	
    fac  = CrosssecTot(E)/CrosssecCX(E);
    fac *= interpolation(Table->S, dr) * g(giveCathodeRadius(),dr) / ( 1.0 - pow(giveTransparency() * g(0,dr),2));
    fac *= 1 - exp( -2 * ngas * CrosssecCX(E) * giveCathodeRadius());
	fac *= pow(dr,2);
	
    return fac;
}

double I_c4(void)
{
    double term;
    double (*FuncPtr)(double);
    FuncPtr = &I_c4inte;

    term = NIntegration(*FuncPtr, giveCathodeRadius(), giveAnodeRadius());

    term *= 4 * 3.141529 * Q_ELECTRON * term; 

    return term;
}


