
/**
    Four function are needed to calculate the neutron production at given radius.

    a function for ingoing ions in the cathode
    a function for outgoing ions in the cathode
    a function for ingoing ions outside the cathode
    a function for outgoing ions outside the cathode
*/

#include <math.h>

#include "CrossSections.h"
#include "PotentialFunctions.h"
#include "constants.h"
#include "MathFunctions.h"
#include "SurvivalFunctions.h"

#include "ParticleFlux.h"

double Sfi_OutMinInte(double r, double dr)
{
    double fac  = CrosssecFusion(ParticleEnergy2(r,dr)) * pow(dr,2) * interpolation(r);
           fac *= g(r,dr)/( 1 - pow(fusor.Tc*g(0,dr),2) );

    return fac;
}

// outside cathode inward.
double Sfi_OutMin(double r)
{
    double S;
    double term1, term2;

    double (*Sfi_OutMinIntePtr)(double, double);
    Sfi_OutMinIntePtr = &Sfi_OutMinInte;

    term1 = NIntegration_2(*Sfi_OutMinIntePtr, r, r, fusor.b);
    term1 *= ngas * EdgeIonFlux;

    term2 = ngas * pow(fusor.b/r,2) * EdgeIonFlux * f(r) * CrosssecFusion(ParticleEnergy1(r));

    S = term1 + term2;

    return S;
}

// outside cathode outward.
double Sfi_OutPlusInte(double r, double dr)
{
    double fac  = CrosssecFusion(ParticleEnergy2(r,dr)) * pow(dr,2) * interpolation(r);
           fac *= (pow(fusor.Tc * g(0,dr),2)/( 1 - pow(fusor.Tc*g(0,dr),2) )) * 1 / g(r,dr);

    return fac;
}

double Sfi_OutPlus(double r)
{
    double S;
    double term1, term2;

    double (*Sfi_OutPlusIntePtr)(double, double);
    Sfi_OutPlusIntePtr = &Sfi_OutPlusInte;

    term1 = NIntegration_2(*Sfi_OutPlusIntePtr, r, r, fusor.b);
    term1 *= ngas * EdgeIonFlux;

    term2 = ngas * pow(fusor.b/r,2) * EdgeIonFlux * (pow(f(0),2) / f(r)) * CrosssecFusion(ParticleEnergy1(r));

    S = term1 + term2;

    return S;
}

double Sfi_InMinInte(double r, double dr)
{
    double fac  = CrosssecFusion(ParticleEnergy2(r,dr)) * pow(dr,2) * interpolation(r);
           fac *= fusor.Tc * g(r,dr) / ( 1 - pow(fusor.Tc*g(0,dr),2) ) ;
           fac *= exp(ngas * CrosssecCX(ParticleEnergy2(fusor.a,dr)) * ( r - fusor.a));

    return fac;
}

double Sfi_InMin(double r)
{
    double S;
    double term1, term2;

    double (*Sfi_InMinIntePtr)(double, double);
    Sfi_InMinIntePtr = &Sfi_InMinInte;

    term1 = NIntegration_2(*Sfi_InMinIntePtr, r, fusor.a, fusor.b);
    term1 *= ngas * EdgeIonFlux;

    term2  = ngas * pow(fusor.b/r,2) * EdgeIonFlux * f(fusor.a) * CrosssecFusion(ParticleEnergy1(r));
    term2 *= exp(ngas * CrosssecCX(-fusor.V0) * ( r - fusor.a));

    S = term1 + term2;

    return S;
}

double Sfi_InPlusInte(double r, double dr)
{
    double fac = CrosssecFusion(ParticleEnergy2(r,dr)) * pow(dr,2) * interpolation(r);
           fac *= pow(fusor.Tc * g(0,dr),2) / ( 1 - pow(fusor.Tc*g(0,dr),2) ) * 1 / g(fusor.a,dr) ;
           fac *= exp( -1 * ngas * CrosssecCX(ParticleEnergy2(fusor.a,dr)) * ( r - fusor.a));

    return fac;
}

double Sfi_InPlus(double r)
{
    double S;
    double term1, term2;

    double (*Sfi_InPlusPtr)(double, double);
    Sfi_InPlusPtr = &Sfi_InPlusInte;

    term1 = NIntegration_2(*Sfi_InPlusPtr, r, fusor.a, fusor.b);
    term1 *= ngas * EdgeIonFlux;

    term2  = ngas * pow(fusor.b/r,2) * EdgeIonFlux * pow(f(fusor.a),2) * CrosssecFusion(ParticleEnergy1(r)) / f(fusor.a);
    term2 *= exp( -1 * ngas * CrosssecCX(-fusor.V0) * ( r - fusor.a));

    S = term1 + term2;

    return S;
}

