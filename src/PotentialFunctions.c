
/** \brief This function returns the potential on position r. Input is a double r in
 *  meter and return is a double phi in KiloVolt.
 *
 * \param r = radius
 * \return potential in V
 */

#include <math.h>
#include <stdlib.h>
#include "constants.h"

#include "PotentialFunctions.h"

double Potential_Phi(double r)
{
    double phi;

    if ( r <= fusor->a)
        phi = fusor->V0;
    else
        phi = (fusor->a * (fusor->b - r)  * fusor->V0) / (r * ( fusor->b - fusor->a ));

    return phi;
}

double Potential_Phi_Inv(double E)
{
    double r_left, r_right, r;
    double Etmp;

    if ( E > fusor->V0)
        return fusor->b;

    r_left  = fusor->a;
    r_right = fusor->b;

    r = r_left + ( r_right - r_left )/2.0;

    Etmp = Potential_Phi(r);

    while ( abs(Etmp - E) > 0.1 )
    {
        if ( Etmp < E )
        {
            // energy is smaller than needed so the radius should be bigger. From larger distance.
            r_left = r;

            r = r_left + ( r_right - r_left )/2.0;

            Etmp = Potential_Phi(r);
        }
        else
        {
            // energy is bigger so originates from larger distance.
            r_right = r;

            r = r_left + ( r_right - r_left )/2.0;

            Etmp = Potential_Phi(r);
        }
    }


    return r;
}


/** \brief calculates the energy of particle at radius r coming from outer edge
 *
 * \param r = radius where particle is
 * \return energie in eV
 *
 */
double ParticleEnergy1(double r)
{
    double energy;
    energy = -q * Potential_Phi(r) + 0.001;
    return energy;
}

/** \brief Calculates the energy of particle at radius r coming from born radius r1
 *
 * \param r = radius where particle is
 * \param r1= radius where particle is born
 * \return energy in eV.
 *
 */
double ParticleEnergy2(double r, double dr)
{
    double energy;
    energy = q * (Potential_Phi(dr) - Potential_Phi(r)) + 0.1;
    return energy;
}



