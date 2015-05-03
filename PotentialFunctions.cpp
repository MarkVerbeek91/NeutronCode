
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
        phi = (giveCathodeRadius() * (giveAnodeRadius() - r)  * giveVoltage()) / (r * ( giveAnodeRadius() - giveCathodeRadius() ));

    return phi;
}

double Potential_Phi_Inv(double E)
{
    double r;

	r  = giveAnodeRadius() * giveCathodeRadius();
	r /= (E/giveVoltage()) * ( giveAnodeRadius() - giveCathodeRadius() ) - giveCathodeRadius();

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
double ParticleEnergy2(double r, double dr)
{
    double energy;
    energy = (Potential_Phi(dr) - Potential_Phi(r)) + 4;
    return energy;
}



