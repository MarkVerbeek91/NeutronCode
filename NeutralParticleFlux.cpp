
#include "math.h"

#include "constants.h"
#include "CrossSections.h"
#include "SurvivalFunctions.h"
#include "PotentialFunctions.h"
#include "MathFunctions.h"
#include "NeutralParticleFlux.h"

/*****************************************************************************
`*                     Neutrals from Class I ions                            *
 *****************************************************************************/

/** \brief Flux INSIDE the cathode, INWARDS. From Class I ions.
 * This function integrate equation 40 to get the neutral particle flux at
 * wanted radii r.
 *
 * \param r: the position where the flux in wanted
 * \return S flux.
 *
 */
double Sfn1_InMin( double r )
{
    double S;

    double (*Sfn1_InMinIntePtr)(double, double);
    Sfn1_InMinIntePtr = &Sfn1_InMinInte;

    S = NIntegration_2(*Sfn1_InMinIntePtr, r, r, giveAnodeRadius());

    S *= pow(ngas * giveAnodeRadius()/r,2) * EdgeIonFlux ;

    return S;
}

double Sfn1_InMinInte(double r, double dr )
{
    return CrosssecFusion(ParticleEnergy2(r,dr)) * f(dr);
}

/** \brief Flux OUTSIDE the cathode, INWARDS. From Class I ions.
 * This function integrate equation 41 to get the neutral particle flux at
 * wanted radii r.
 *
 * \param r: the position where the flux in wanted
 * \return S flux.
 *
 */
double Sfn1_OutMin( double r )
{
    double S;

    double (*Sfn1_OutMinIntePtr)(double, double);
    Sfn1_OutMinIntePtr = &Sfn1_OutMinInte;

    S = NIntegration_2(*Sfn1_OutMinIntePtr, r, r, giveAnodeRadius());

    S *= pow(ngas *giveAnodeRadius()/r,2) * EdgeIonFlux;

    return S;
}

double Sfn1_OutMinInte(double r, double dr )
{
    return CrosssecFusion(ParticleEnergy2(r,dr)) * f(dr);
}

/** \brief Flux INSIDE the cathode, OUTWARDS. From Class I ions.
 * This function integrate equation 47 to get the neutral particle flux at
 * wanted radii r.
 *
 * \param r: the position where the flux in wanted
 * \return S flux.
 *
 */
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

double Sfn1_InPlusInte(double r, double dr)
{
    return CrosssecFusion(ParticleEnergy2(r,dr)) * f(dr);
}

/** \brief Flux OUTSIDE the cathode, OUTWARDS. From Class I ions.
 * This function integrate equation 51 to get the neutral particle flux at
 * wanted radii r.
 *
 * \param r: the position where the flux in wanted
 * \return S flux.
 *
 */
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

double Sfn1_OutPlusInte(double r, double dr)
{
    if ( r < dr )
        return CrosssecFusion(ParticleEnergy2(r,dr)) * (f(dr) + pow(f(0),2) / f(dr));
    else
        return CrosssecFusion(ParticleEnergy2(r,dr)) *  f(dr);
}

/*****************************************************************************
`*                      Neutrals from Class II ions                          *
 *****************************************************************************/

// neutral particle flux from II class ions, outside the cathode, outwards.
double Sfn2_InMin( double r );
{
    return
}
double Sfn2_InMinInte(double r, double dr );
{

}

// neutral particle flux from II class ions, outside the cathode, inwards.
double Sfn2_OutMin( double r );
double Sfn2_OutMinInte(double r, double dr );

// neutral particle flux from II class ions, inside the cathode, outwards.
double Sfn2_InPlus(double r);
double Sfn2_InPlusInte(double r, double dr);

// neutral particle flux from II class ions, inside the cathode, inwards.
double Sfn2_OutPlus(double r);
double Sfn2_OutPlusInte(double r, double dr);


