
/**
    To get the total Neutron production, the neutron production needs to be integrated
    over the volume.
*/

#include <stdio.h>

#include "constants.h"
#include "MathFunctions.h"
#include "FusionReactionRateIons.h"
#include "FusionReactionRateNeutralsClassI.h"
#include "FusionReactionRateNeutralsClassII.h"

#include "NeutronProductionRate.h"

// TODO: update to the new formulas
double Nps(void)
{
    double NPS, tempNPS;
    double (*FunctPtr)(double);

    // Neutrons from Ions
    FunctPtr = &FusionRateIons;

    tempNPS = NIntegration(*FunctPtr, 0.00001, giveAnodeRadius() - 0.000001);
    NPS  = tempNPS;
    fprintf(stdout,"# Neutrons from ...\n");
	fprintf(stdout,"#        Ions          : \t %.2E \n", tempNPS);

    // Neutrons from Neutrals from Class I ions.
    FunctPtr = &FusionRateNeutralsClassI;

    tempNPS = NIntegration(*FunctPtr, 0.00001, giveAnodeRadius() - 0.000001);
    NPS += tempNPS;
	fprintf(stdout,"#        Neutrals 1    : \t %.2E \n", tempNPS);

/*
    // Neutrons from Neutrals from Class II ions.
    FunctPtr = &FusionRateNeutralsClassI;

    tempNPS = NIntegration(*FunctPtr, 0.00001, giveAnodeRadius() - 0.000001);
    NPS += tempNPS;
	fprintf(stdout,"#        Neutrals 2    : \t %.2E \n", tempNPS);

*/
    return NPS;
}

