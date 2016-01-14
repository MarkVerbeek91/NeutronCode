
/**
    To get the total Neutron production, the neutron production needs to be integrated
    over the volume.
*/

#include <stdio.h>

#include "constants.h"
#include "MathFunctions.h"
#include "NeutronProductionIons.h"
#include "NeutronProductionNeutralsClassI.h"
#include "NeutronProductionNeutralsClassII.h"

#include "NeutronProductionRate.h"

// TODO: update to the new formulas
double Nps(void)
{
    double NPS, tempNPS;
    double (*FunctPtr)(double);

    // Neutrons from Ions
    FunctPtr = &NeutronsIonFluxInwards;

    tempNPS = NIntegration(*FunctPtr, 0.00001, giveAnodeRadius() - 0.000001);
    NPS  = tempNPS;
    fprintf(stdout,"# Neutrons from ...\n");
	fprintf(stdout,"# Ions Inwards         : \t %.2E \n", tempNPS);

    FunctPtr = &NeutronsIonFluxOutwards;

    tempNPS = NIntegration(*FunctPtr, 0.00001, giveAnodeRadius() - 0.000001);
    NPS += tempNPS;
	fprintf(stdout,"# Ions Outwards        : \t %.2E \n", tempNPS);

    // Neutrons from Neutrals from Class I ions.
    FunctPtr = &NeutronsNeutralsClassIFluxInwards;

    tempNPS = NIntegration(*FunctPtr, 0.00001, giveAnodeRadius() - 0.000001);
    NPS += tempNPS;
	fprintf(stdout,"# Neutrals 1 Inwards   : \t %.2E \n", tempNPS);

    FunctPtr = &NeutronsNeutralsClassIFluxOutwards;

    tempNPS = NIntegration(*FunctPtr, 0.00001, giveAnodeRadius() - 0.000001);
    NPS += tempNPS;
	fprintf(stdout,"# Neutrals 1 Outwards  : \t %.2E \n", tempNPS);


/*
    // Neutrons from Neutrals from Class II ions.
    FunctPtr = &NeutronsNeutralsClassIIFluxInwards;

    tempNPS = NIntegration(*FunctPtr, 0.00001, giveAnodeRadius() - 0.000001);
    NPS += tempNPS;
	fprintf(stdout,"# Neutrals 2 Inwards   : \t %.2E \n", tempNPS);

    FunctPtr = &NeutronsNeutralsClassIIFluxOutwards;

    tempNPS = NIntegration(*FunctPtr, 0.00001, giveAnodeRadius() - 0.000001);
    NPS += tempNPS;
	fprintf(stdout,"# Neutrals 2 Outwards  : \t %.2E \n", tempNPS);
*/
    return NPS;
}

