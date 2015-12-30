/**
 * \brief All functions to calculate the spectrum of Ions moving through the 
 * fusor. This are two main functions. One for inward traveling ions and 
 * another for outward traving ions. 
 */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "constants.h"
#include "MathFunctions.h"
#include "PotentialFunctions.h"
#include "SurvivalFunctions.h"
#include "Crosssections.h"

#include "IonSpectrum.h"

/** \brief Flux of ions moving inwards in the fusor.
 *
 * \param r the radius where to calculate the spectrum
 * \param E the energy of ions
 * \return The flux of the ions of on r with E energy.
 *
 */
// equations 28 and 30;
double IonSpectrumInwards(double r, double E)
{
    double flux, dr;
    double term1 = 0, term2 = 0;
	double (*PhiPtr)(double);
    PhiPtr = &Potential_Phi;

	dr = Potential_Phi_Inv(E);
	
    #ifdef DEBUG_PARAMETER
    if ( giveAnodeRadius() < dr || dr < r )
    {
        printf("IonSpectrumInwards function error: AnodeRadius < dr or dr < r\n");
        return NAN;
    }
    #endif

    term1  = 1/giveq();
    term1 *= pow(dr/r,2);
    term1 *= interpolation(Table->S, dr) / abs(differentiat(*PhiPtr, dr));
    term1 /= 1 - pow(giveTransparency() * g(0,dr),2);

    if ( giveCathodeRadius() < r )
    {
        // inwards flux outside the cathode
        term1  *= g(r,dr);

        if ( DELTA(E - ParticleEnergy1(r)) )
            term2  = pow(giveAnodeRadius()/r,2) * EdgeIonFlux * f(r);

        flux = term1 + term2;
    }
    else
    {
        // inwards flux inside the cathode
        term1 *= g(giveCathodeRadius(),dr) * exp(ngas * CrosssecCX(ParticleEnergy2(giveCathodeRadius(),dr)) * (r - giveCathodeRadius()));

        if ( DELTA(E - ParticleEnergy1(r)) )
        {
            term2  = pow(giveAnodeRadius()/r,2) * EdgeIonFlux ;
            term2 *= f(giveCathodeRadius()) * exp(ngas * CrosssecCX(ParticleEnergy1(giveCathodeRadius())) * (r - giveCathodeRadius()));
        }

        flux = giveTransparency() * (term1 + term2);
    }

    return flux;
}

/** \brief Flux of ions moving outwards in the fusor.
 *
 * \param r the radius where to calculate the spectrum
 * \param E the energy of ions
 * \return The flux of the ions of on r with E energy.
 *
 */
// equations 29 and 31
double IonSpectrumOutwards(double r, double E)
{
    double flux;
    double term1 = 0, term2 = 0;
    double dr = Potential_Phi_Inv(E);

    #ifdef DEBUG_PARAMETER
    if ( giveCathodeRadius() < dr || dr < giveAnodeRadius() )
    {
        printf("IonSpectrumOutwards error: dr < b or dr < r\n");
        return NAN;
    }
    #endif

    double (*PhiPtr)(double);
    PhiPtr = &Potential_Phi;

    term1  = 1/giveq();
    term1 *= pow(dr/r,2);
    term1 *= interpolation(Table->S, dr) / abs(differentiat(*PhiPtr, dr));
    term1 *= pow(g(0, dr),2) / ( 1 - pow(giveTransparency() * g(0,dr),2) );

    if ( r < giveCathodeRadius())
    {
        // outwards flux inside the cathode
        term1 /= g(giveCathodeRadius(), dr) * exp(ngas * CrosssecCX(ParticleEnergy2(giveCathodeRadius(),dr)) * ( r - giveCathodeRadius()));

        if ( DELTA(E - ParticleEnergy1(r)) )
        {
            term2  = pow(giveAnodeRadius() * f(0)/r,2) * EdgeIonFlux;
            term2 /= f(giveCathodeRadius()) * exp(ngas * CrosssecCX(ParticleEnergy1(giveCathodeRadius())) * (r - giveCathodeRadius()));
        }

        flux = giveTransparency() * (term1 + term2);
    }
    else
    {
        // outwards flux outside the cathode
        term1 /= g(r, dr);

        if ( DELTA(E - ParticleEnergy1(r)) )
            term2 =  pow(giveAnodeRadius() * f(0)/r,2) * EdgeIonFlux / f(r);

        flux = pow(giveTransparency(),2) * (term1 + term2);
    }

    return flux;
}