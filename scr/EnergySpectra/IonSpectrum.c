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
#include "CrossSections.h"

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

    dr = Potential_Phi_Inv( E/giveq() + Potential_Phi(r) );

    if ( dr < r )
        dr = roundf(dr * 10000) / 10000 + 1e-5;

    term1  = 1/giveq();
    term1 *= pow(dr/r,2);
    term1 *= interpolation(Table->S, dr) / abs(differentiat(*PhiPtr, dr));
    term1 /= 1 - pow(fusor->Tc * g(0,dr),2);

    if ( fusor->a < r )
    {
        // inwards flux outside the cathode
        if ( dr < r || fusor->b < dr )
        {
            printf("IonSpectrumInwards error: dr < r or b < dr:\n \tr = %E\n \tdr = %E\n \tE = %E\n", r, dr, E);

            return NAN;
        }

        term1  *= g(r,dr);

        if ( DELTA(E - ParticleEnergy1(r)) )
            term2  = pow(fusor->b/r,2) * EdgeIonFlux * f(r);

        flux = term1 + term2;
    }
    else
    {
        if ( dr < fusor->a || fusor->b < dr )
        {
            printf("IonSpectrumInwards error: dr < a or b < dr:\n \tr = %E\n \tdr = %E\n \tE = %E\n", r, dr, E);

            return NAN;
        }

        // inwards flux inside the cathode
        term1 *= g(fusor->a,dr) * exp(ngas * CrosssecCX(ParticleEnergy2(fusor->a,dr)) * (r - fusor->a));

        if ( DELTA(E - ParticleEnergy1(r)) )
        {
            term2  = pow(fusor->b/r,2) * EdgeIonFlux ;
            term2 *= f(fusor->a) * exp(ngas * CrosssecCX(ParticleEnergy1(fusor->a)) * (r - fusor->a));
        }

        flux = fusor->Tc * (term1 + term2);
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
    double flux, dr;
    double term1 = 0, term2 = 0;
    double (*PhiPtr)(double);
    PhiPtr = &Potential_Phi;

    dr = Potential_Phi_Inv( E/giveq() + Potential_Phi(r) );

    if ( dr < r )
    	dr = roundf(dr * 10000) / 10000 + 1e-5;

    term1  = 1/giveq();
    term1 *= pow(dr/r,2);
    term1 *= interpolation(Table->S, dr) / abs(differentiat(*PhiPtr, dr));
    term1 *= pow(g(0, dr),2) / ( 1 - pow(fusor->Tc * g(0,dr),2) );

    if ( r < fusor->a)
    {
        // outwards flux inside the cathode
        if ( dr < fusor->a || fusor->b < dr )
        {
            printf("IonSpectrumOutwards error: dr < a or b < dr:\n \tr = %E\n \tdr = %E\n \tE = %E\n", r, dr, E);
            return NAN;
        }

        term1 /= g(fusor->a, dr) * exp(ngas * CrosssecCX(ParticleEnergy2(fusor->a,dr)) * ( r - fusor->a));

        if ( DELTA(E - ParticleEnergy1(r)) )
        {
            term2  = pow(fusor->b * f(0)/r,2) * EdgeIonFlux;
            term2 /= f(fusor->a) * exp(ngas * CrosssecCX(ParticleEnergy1(fusor->a)) * (r - fusor->a));
        }

        flux = fusor->Tc * (term1 + term2);
    }
    else
    {
        // outwards flux outside the cathode
        if ( dr < r || fusor->b < dr )
        {
            printf("IonSpectrumInwards error: dr < a or b < dr:\n \tr = %E\n \tdr = %E\n \tE = %E\n", r, dr, E);
            return NAN;
        }

        term1 /= g(r, dr);

        if ( DELTA(E - ParticleEnergy1(r)) )
            term2 =  pow(fusor->b * f(0)/r,2) * EdgeIonFlux / f(r);

        flux = pow(fusor->Tc,2) * (term1 + term2);
    }

    return flux;
}
