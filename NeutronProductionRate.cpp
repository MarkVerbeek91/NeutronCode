
/**
    To get the total Neutron production, the neutron production needs to be integrated
    over the volume.
*/

#include <stdio.h>

#include "constants.h"
#include "MathFunctions.h"
#include "ParticleFlux.h"
#include "NeutronProductionRate.h"


double Nps(void)
{
    double NPS;

    double (*Sfi_InMinPtr)(double);
    Sfi_InMinPtr = &Sfi_InMin;
    double (*Sfi_InPlusPtr)(double);
    Sfi_InPlusPtr = &Sfi_InPlus;
    double (*Sfi_OutMinPtr)(double);
    Sfi_OutMinPtr = &Sfi_OutMin;
    double (*Sfi_OutPlusPtr)(double);
    Sfi_OutPlusPtr = &Sfi_OutPlus;

    printf(" - Ions inwards inside cathode\n");
    NPS  = NIntegration(*Sfi_InMinPtr, 0.01, giveCathodeRadius() - 0.000001);
    printf(" - Ions outwards inside cathode\n");
    NPS += NIntegration(*Sfi_InPlusPtr, 0.01, giveCathodeRadius() - 0.000001);
    printf(" - Ions inward outside cathode\n");
    NPS += NIntegration(*Sfi_OutMinPtr, giveCathodeRadius(), giveAnodeRadius() - 0.000001);
    printf(" - Ions outwards outside cathode\n");
    NPS += NIntegration(*Sfi_OutPlusPtr, giveCathodeRadius(), giveAnodeRadius() - 0.000001);

    return NPS;
}

