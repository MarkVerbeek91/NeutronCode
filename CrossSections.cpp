/**
    This file contains all cross sections functions

*/

#include <math.h>
#include <stdio.h>

#include "CrossSectionsCont.h"
#include "CrossSections.h"

void InitCrossSectionConstands(void)
{

    CS_cx.A1cx =   3.245     ;
    CS_cx.A2cx = 235.88      ;
    CS_cx.A3cx =   0.038371  ;
    CS_cx.A4cx =   3.8068e-6 ;
    CS_cx.A5cx =   1.1832e-10;
    CS_cx.A6cx =   2.3713    ;

    CS_Ion.A1Ion = 12.899    ;
    CS_Ion.A2Ion = 61.897    ;
    CS_Ion.A3Ion =  9.27e3   ;
    CS_Ion.A4Ion =  4.9749e-4;
    CS_Ion.A5Ion =  3.989e-2 ;
    CS_Ion.A6Ion = -1.59     ;
    CS_Ion.A7Ion =  3.1834   ;
    CS_Ion.A8Ion = -3.7154   ;

    CS_Fusion.A1Fusion =  47.88     ;
    CS_Fusion.A2Fusion = 482        ;
    CS_Fusion.A3Fusion =   3.08e-4  ;
    CS_Fusion.A4Fusion =   1.177    ;
    CS_Fusion.A5Fusion =   0        ;

}

/** \brief Calculates the CX crosssection for particle with energy E
 *
 * \param E = energy of particle
 * \param
 * \return crosssection in m2
 *
 */
double CrosssecCX(double E)
{
    if ( E < 1e-4 )
    {
        printf("Cross section CX called for zero energy\n");
        E = 1e-5;
    }

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
    if ( E < 1e-4 )
    {
        printf("Cross section Ion called for zero energy\n");
        E = 1e-5;
    }

    double crosssection;

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
    if ( E < 1e-4 )
    {
        printf("Cross section Fusion called for zero energy\n");
        E = 1e-5;
    }

    double crosssection, energy = E/1000.;
    crosssection = 1e-28 * (CS_Fusion.A5Fusion + (CS_Fusion.A2Fusion / (pow(CS_Fusion.A4Fusion - CS_Fusion.A3Fusion * energy,2)+1)));
    crosssection = crosssection /(energy * (exp(CS_Fusion.A1Fusion / sqrt(energy))-1));
    return crosssection;
}
