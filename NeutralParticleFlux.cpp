
#include "math.h"

#include "constants.h"
#include "CrossSections.h"
#include "SurvivalFunctions.h"
#include "PotentialFunctions.h"
#include "MathFunctions.h"

// neutral particle flux from I class ions, outside the cathode, inwards.
/**
  * Equation 40 integrated
*/
double Sfn1_InMinInte(double r, double dr )
{
    return CrosssecFusion(ParticleEnergy2(r,dr)) * f(dr);
}

double Sfn1_InMin( double r )
{
    double S;

    double (*Sfn1_InMinIntePtr)(double, double);
    Sfn1_InMinIntePtr = &Sfn1_InMinInte;

    S = NIntegration_2(*Sfn1_InMinIntePtr, r, r, giveAnodeRadius());

    S *= pow(ngas * giveAnodeRadius()/r,2) * EdgeIonFlux ;

    return S;
}

// neutral particle flux from I class ions, inside the cathode, inwards.
/**
  * Equation 41 integrated
  */
double Sfn1_OutMinInte(double r, double dr )
{
    return CrosssecFusion(ParticleEnergy2(r,dr)) * f(dr);
}

double Sfn1_OutMin( double r )
{
    double S;

    double (*Sfn1_OutMinIntePtr)(double, double);
    Sfn1_OutMinIntePtr = &Sfn1_OutMinInte;

    S = NIntegration_2(*Sfn1_OutMinIntePtr, r, r, giveAnodeRadius());

    S *= pow(ngas *giveAnodeRadius()/r,2) * EdgeIonFlux;

    return S;
}

// neutral particle flux from I class ions, inside the cathode, outwards.
/**
    Equation 47 integrated
*/

double Sfn1_InPlusInte(double r, double dr)
{
    return CrosssecFusion(ParticleEnergy2(r,dr)) * f(dr);
}

// outside cathode inward.
double Sfn1_InPlus(double r)
{
    double S;
    double term1, term2;

    double (*Sfn1_InPlusIntePtr)(double, double);
    Sfn1_InPlusIntePtr = &Sfn1_InPlusInte;

    term1 = NIntegration_2(*Sfn1_InPlusIntePtr, r, 0, giveCathodeRadius());
    term1 *= pow( ngas * giveAnodeRadius()/r,2) * EdgeIonFlux * giveTransparency();

    term2 = ngas * pow(giveAnodeRadius()/r,2) * EdgeIonFlux * giveTransparency() * f(giveAnodeRadius());
    term2 *= ( 1 - exp(ngas * CrosssecCX(-giveVoltage()) * ( r + giveCathodeRadius())));

    S = term1 + term2;

    return S;
}


// neutral particle flux from I class ions, outside the cathode, outwards.
/** \brief equation 51 integrated
 *
 * \param
 * \param
 * \return
 *
 */

double Sfn1_OutPlusInte(double r, double dr)
{
    if ( r < dr )
        return CrosssecFusion(ParticleEnergy2(r,dr)) * (f(dr) + pow(f(0),2) / f(dr));
    else
        return CrosssecFusion(ParticleEnergy2(r,dr)) *  f(dr);
}

// outside cathode inward.
double Sfn1_OutPlus(double r)
{
    double S;
    double term1, term2;

    double (*Sfn1_OutPlusIntePtr)(double, double);
    Sfn1_OutPlusIntePtr = &Sfn1_OutPlusInte;

    term1 = NIntegration_2(*Sfn1_OutPlusIntePtr, r, 0, giveCathodeRadius());
    term1 *= pow( ngas * giveAnodeRadius() * giveTransparency() / r,2) * EdgeIonFlux ;

    term2 = ngas * pow(giveAnodeRadius() * giveTransparency()/r,2) * EdgeIonFlux * f(giveCathodeRadius());
    term2 *= ( 1 - exp(ngas * CrosssecCX(-giveVoltage()) * ( r + giveCathodeRadius())));

    S = term1 + term2;

    return S;
}





// neutral particle flux from II class ions, outside the cathode, inwards.


// neutral particle flux from II class ions, outside the cathode, outwards.


// neutral particle flux from II class ions, inside the cathode, inwards.


// neutral particle flux from II class ions, inside the cathode, outwards.


