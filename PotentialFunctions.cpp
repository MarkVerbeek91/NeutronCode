
/** \brief This function returns the potential on position r. Input is a double r in
 *  meter and return is a double phi in KiloVolt.
 *
 * \param r = radius
 * \return potential in V
 */

#include "constants.h"

#include "PotentialFunctions.h"

double Potential_Phi(double r)
{
    double phi;

    if ( r <= giveCathodeRadius())
        phi = giveVoltage();
    else
        phi = (giveCathodeRadius() * (giveAnodeRadius() - r) * -1 * giveVoltage()) / (r * (giveCathodeRadius() - giveAnodeRadius()));

    return phi;
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
    energy = - Potential_Phi(r) + 0.001;
    return energy;
}

/** \brief Calculates the energy of particle at radius r coming from born radius r1
 *
 * \param r = radius where particle is
 * \param r1= radius where particle is born
 * \return energy in eV.
 *
 */
double ParticleEnergy2(double r, double r1)
{
    double energy;
    energy = (Potential_Phi(r1) - Potential_Phi(r)) + 4;
    return energy;
}


