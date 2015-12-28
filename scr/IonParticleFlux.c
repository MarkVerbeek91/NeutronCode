
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

#include "IonParticleFlux.h"

/** \brief The ion particle flux, OUTSIDE the cathode for ions moving INWARDS
 *
 * \param r is the radius to calculate the flux at
 *
 * \return flux
 */
double Sfi_OutMin(double r)
{
    double S;
    double term1, term2;

    double (*Sfi_OutMinIntePtr)(double, double);
    Sfi_OutMinIntePtr = &Sfi_OutMinInte;

    term1 = NIntegration_2(*Sfi_OutMinIntePtr, r, r, giveAnodeRadius());
    term1 *= ngas * EdgeIonFlux;

    term2 = ngas * pow(giveAnodeRadius()/r,2) * EdgeIonFlux * f(r) * CrosssecFusion(ParticleEnergy1(r));

    S = term1 + term2;

    return S;
}

/** \brief The integrand of the function above.
 *
 * \param r:  the radius where in the flux in wanted to be known
 * \param dr: variable of integration
 * \return the flux of ions at r from dr.
 */
double Sfi_OutMinInte(double r, double dr)
{
    return CrosssecFusion(ParticleEnergy2(r,dr)) * pow(dr/r,2) * interpolation(Table.S, dr)
            * g(r,dr)/( 1 - pow(giveTransparency()*g(0,dr),2) );
}

/** \brief The ion particle flux, OUTSIDE the cathode for ions moving OUTWARTS
 *
 * \param r is the radius to calculate the flux at
 *
 * \return flux
 */
double Sfi_OutPlus(double r)
{
    double S;
    double term1, term2;

    double (*Sfi_OutPlusIntePtr)(double, double);
    Sfi_OutPlusIntePtr = &Sfi_OutPlusInte;

    term1 = NIntegration_2(*Sfi_OutPlusIntePtr, r, r, giveAnodeRadius());
    term1 *= ngas * EdgeIonFlux;

    term2 = ngas * pow(giveAnodeRadius()/r,2) * EdgeIonFlux * (pow(f(0),2) / f(r)) * CrosssecFusion(ParticleEnergy1(r));

    S = term1 + term2;

    return S;
}

/** \brief The integrand of the function above.
 *
 * \param r:  the radius where in the flux in wanted to be known
 * \param dr: variable of integration
 * \return the flux of ions at r from dr.
 */
double Sfi_OutPlusInte(double r, double dr)
{
    return CrosssecFusion(ParticleEnergy2(r,dr)) * pow(dr/r,2) * interpolation(Table.S, dr)
           * (pow(giveTransparency() * g(0,dr),2)/( 1 - pow(giveTransparency()*g(0,dr),2) )) * 1 / g(r,dr);
}

/** \brief The ion particle flux, INSIDE the cathode for ions moving INWARTS
 *
 * \param r is the radius to calculate the flux at
 *
 * \return flux
 */
double Sfi_InMin(double r)
{
    double S;
    double term1, term2;

    double (*Sfi_InMinIntePtr)(double, double);
    Sfi_InMinIntePtr = &Sfi_InMinInte;

    term1 = NIntegration_2(*Sfi_InMinIntePtr, r, giveCathodeRadius(), giveAnodeRadius());
    term1 *= ngas * EdgeIonFlux;

    term2  = ngas * pow(giveAnodeRadius()/r,2) * EdgeIonFlux * f(giveCathodeRadius()) * CrosssecFusion(ParticleEnergy1(r));
    term2 *= exp(ngas * CrosssecCX(-giveVoltage()) * ( r - giveCathodeRadius()));

    S = term1 + term2;

    return S;
}

/** \brief The integrand of the function above.
 *
 * \param r:  the radius where in the flux in wanted to be known
 * \param dr: variable of integration
 * \return the flux of ions at r from dr.
 */
double Sfi_InMinInte(double r, double dr)
{
    return CrosssecFusion(ParticleEnergy2(r,dr)) * pow(dr/r,2) * interpolation(Table.S, dr)
           * giveTransparency() * g(r,dr) / ( 1 - pow(giveTransparency()*g(0,dr),2) )
           * exp(ngas * CrosssecCX(ParticleEnergy2(giveCathodeRadius(),dr)) * ( r - giveCathodeRadius()));
}

/** \brief The ion particle flux, INSIDE the cathode for ions moving OUTWARTS
 *
 * \param r is the radius to calculate the flux at
 *
 * \return flux
 */
double Sfi_InPlus(double r)
{
    double S;
    double term1, term2;

    double (*Sfi_InPlusPtr)(double, double);
    Sfi_InPlusPtr = &Sfi_InPlusInte;

    term1 = NIntegration_2(*Sfi_InPlusPtr, r, giveCathodeRadius(), giveAnodeRadius());
    term1 *= ngas * EdgeIonFlux;

    term2  = ngas * pow(giveAnodeRadius()/r,2) * EdgeIonFlux * pow(f(giveCathodeRadius()),2) * CrosssecFusion(ParticleEnergy1(r)) / f(giveCathodeRadius());
    term2 *= exp( -1 * ngas * CrosssecCX(-giveVoltage()) * ( r - giveCathodeRadius()));

    S = term1 + term2;

    return S;
}

/** \brief The integrand of the function above.
 *
 * \param r:  the radius where in the flux in wanted to be known
 * \param dr: variable of integration
 * \return the flux of ions at r from dr.
 */
double Sfi_InPlusInte(double r, double dr)
{
    return CrosssecFusion(ParticleEnergy2(r,dr)) * pow(dr/r,2) * interpolation(Table.S, dr)
           * pow(giveTransparency() * g(0,dr),2) / ( 1 - pow(giveTransparency()*g(0,dr),2) ) * 1 / g(giveCathodeRadius(),dr)
           * exp( -1 * ngas * CrosssecCX(ParticleEnergy2(giveCathodeRadius(),dr)) * ( r - giveCathodeRadius()));
}
