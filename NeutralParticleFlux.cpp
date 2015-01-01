
#include "math.h"

#include "constants.h"
#include "CrossSections.h"
#include "SurvivalFunctions.h"
#include "PotentialFunctions.h"
#include "MathFunctions.h"

// neutral particle flux from I class ions, outside the cathode, inwards.
double Sfn1_InMinInte(double r, double dr )
{
    double fac  = CrosssecFusion(ParticleEnergy2(r,dr)) * pow(giveAnodeRadius()/r,2) ;
           fac *= ngas * EdgeIonFlux * f(dr);

    return fac;
}

double Sfn1_InMin( double r )
{
    double S;

    double (*Sfn1_InMinIntePtr)(double, double);
    Sfn1_InMinIntePtr = &Sfn1_InMinInte;

    S = NIntegration_2(*Sfn1_InMinIntePtr, r, r, giveAnodeRadius());

    return S;
}

// neutral particle flux from I class ions, outside the cathode, outwards.
double Sfn1_OutMinInte(double r, double dr )
{
    double fac  = CrosssecFusion(ParticleEnergy2(r,dr)) * pow(giveAnodeRadius()/r,2) ;
           fac *= ngas * EdgeIonFlux * f(dr);

    return fac;
}

double Sfn1_OutMin( double r )
{
    double S;

    double (*Sfn1_OutMinIntePtr)(double, double);
    Sfn1_OutMinIntePtr = &Sfn1_OutMinInte;

    S = NIntegration_2(*Sfn1_OutMinIntePtr, r, r, giveAnodeRadius());

    return S;
}


// neutral particle flux from I class ions, inside the cathode, inwards.


// neutral particle flux from I class ions, inside the cathode, outwards.



// neutral particle flux from II class ions, outside the cathode, inwards.


// neutral particle flux from II class ions, outside the cathode, outwards.


// neutral particle flux from II class ions, inside the cathode, inwards.


// neutral particle flux from II class ions, inside the cathode, outwards.


