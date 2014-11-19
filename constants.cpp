
#include "constants.h"

/**
 * Declaration of constants and tables. Here the properties of a fusor are
 * descripted. Some phicical constants defined like charge of a proton and
 * all parameters needed for the cross section calculation.
 *
 * The tables will contain the values that take a long time to calculate.
 *
 */

struct Fusor fusor;

fusor.a = 0.05;
fusor.b = 0.25;
fusor.V0 = -40000;
fusor.wire_diameter = 0.005;
fusor.Tc = 0.95;

// defining of some variables
double pressure = 0.5;  // Pa
double Tgas = 400; // K
double ngas = 9.05401e19; //6.022e23 * pressure / (8.314 * Tgas);
double E0 = 0.0001;          // reducing errors
double Itot = 0.1;

double EdgeIonFlux;

Tables Table;

double giveCathodeRadius(void)
{
    return fusor.a;
}

double giveAnodeRadius(void)
{
    return fusor.b;
}

double giveVoltage(void)
{
    return fusor.b;
}

double giveTransparency(void)
{
    return fusor.Tc;
}

double giveq(void)
{
    return q;
}

double giveTgas(void)
{
    return Tgas;
}

double giveNgas(void)
{
    return ngas;
}

double giveItot(void)
{
    return Itot;
}

double giveEdgeIonFlux(void)
{
    return EdgeIonFlux;
}



