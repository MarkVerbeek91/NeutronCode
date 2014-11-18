/**
    This file contains all cross sections

*/

#include <math.h>

#include "CrossSectionsCont.hpp"
#include "CrossSections.hpp"

/** \brief Calculates the CX crosssection for particle with energy E
 *
 * \param E = energy of particle
 * \param
 * \return crosssection in m2
 *
 */
double CrosssecCX(double E)
{
    double crosssection;
    double energy = E/1000.;
    crosssection = 1e-20 * CS_cx.A1cx;
    crosssection = crosssection * log((CS_cx.A2cx/energy) + CS_cx.A6cx);
    crosssection = crosssection / (1 + CS_cx.A3cx * energy + CS_cx.A4cx * pow(energy,3.5) + CS_cx.A5cx * pow(energy,5.4));

    return crosssection;
}

/** \brief Calculates the Ion crosssection for particle with energy E
 *
 * \param E = energy of particle
 * \param
 * \return crosssection in m2
 *
 */
double CrosssecIon(double E)
{
    double crosssection = 0;

    // because the int E is in eV and the formula in keV, the factor 1/1000 is introduced.
    double energy = E/1000.;
    crosssection = (exp(-CS_Ion.A2Ion/energy) * log(1 + CS_Ion.A3Ion * energy)) / energy;
    crosssection = crosssection + CS_Ion.A4Ion * exp(-CS_Ion.A5Ion * energy) / (exp(CS_Ion.A6Ion) + CS_Ion.A7Ion*exp(CS_Ion.A8Ion));
    crosssection = crosssection * 1e-20 * CS_Ion.A1Ion;
    return crosssection;
}

/** \brief Calculates the Tot crosssection for particle with energy E
 *
 * \param E = energy of particle
 * \param
 * \return crosssection in m2
 *
 */
double CrosssecTot(double energy)
{
    double crosssection = 0;
    crosssection = CrosssecCX(energy) + CrosssecIon(energy);
    return crosssection;
}

/** \brief Calculates the fusion crosssection for particle with energy E
 *
 * \param E = energy of particle
 * \param
 * \return crosssection in m2
 *
 */
double CrosssecFusion(double E)
{
    double crosssection, energy = E/1000.;
    crosssection = 1e-28 * (CS_Fusion.A5Fusion + (CS_Fusion.A2Fusion / (pow(CS_Fusion.A4Fusion - CS_Fusion.A3Fusion * energy,2)+1)));
    crosssection = crosssection /(energy * (exp(CS_Fusion.A1Fusion / sqrt(energy))-1));
    return crosssection;
}
