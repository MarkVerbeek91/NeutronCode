
/**
    To get the total Neutron production, the neutron production needs to be integrated
    over the volume.
*/

#include <stdio.h>

#include "constants.h"
#include "MathFunctions.h"
#include "IonParticleFlux.h"
#include "NeutralParticleFlux.h"
#include "NeutronProductionRate.h"


double Nps(void)
{
    double NPS, tempNPS;

    double (*Sfi_InMinPtr)(double);
    Sfi_InMinPtr = &Sfi_InMin;
    double (*Sfi_InPlusPtr)(double);
    Sfi_InPlusPtr = &Sfi_InPlus;
    double (*Sfi_OutMinPtr)(double);
    Sfi_OutMinPtr = &Sfi_OutMin;
    double (*Sfi_OutPlusPtr)(double);
    Sfi_OutPlusPtr = &Sfi_OutPlus;

    printf(        " - Ions inwards inside cathode\n");
    tempNPS = NIntegration(*Sfi_InMinPtr, 0.00001, giveCathodeRadius() - 0.000001);
    NPS  = tempNPS;
    printf("\t %E \n - Ions outwards inside cathode\n", tempNPS);
    tempNPS = NIntegration(*Sfi_InPlusPtr, 0.00001, giveCathodeRadius() - 0.000001);
    NPS += tempNPS;
    printf("\t %E \n - Ions inward outside cathode\n", tempNPS);
    tempNPS = NIntegration(*Sfi_OutMinPtr, giveCathodeRadius(), giveAnodeRadius() - 0.000001);
    NPS  += tempNPS;
    printf("\t %E \n - Ions outwards outside cathode\n", tempNPS);
    tempNPS = NIntegration(*Sfi_OutPlusPtr, giveCathodeRadius(), giveAnodeRadius() - 0.000001);
    NPS  += tempNPS;
    printf("\t %E \n", tempNPS);

    double (*Sfn1_InMinPtr)(double);
    Sfn1_InMinPtr = &Sfn1_InMin;

    double (*Sfn1_OutMinPtr)(double);
    Sfn1_OutMinPtr = &Sfn1_OutMin;

    double Sfn1_InPlus(double r);
    double (*Sfn1_InPlusPtr)(double);
    Sfn1_InPlusPtr = &Sfn1_InPlus;

    double (*Sfn1_OutPlusPtr)(double);
    Sfn1_OutPlusPtr = &Sfn1_OutPlus;

    printf(        " - Neutrals inwards inside cathode from Class I\n");
    tempNPS = NIntegration(*Sfn1_InMinPtr, 0.01, giveCathodeRadius() - 0.0001);
    NPS  += tempNPS;
    printf("\t %E \n - Neutrals outwards inside cathode from Class I\n", tempNPS);
    tempNPS = NIntegration(*Sfn1_InPlusPtr, 0.01, giveCathodeRadius() - 0.0001);
    NPS += tempNPS;
    printf("\t %E \n - Neutrals inward outside cathode from Class I\n", tempNPS);
    tempNPS = NIntegration(*Sfn1_OutMinPtr, giveCathodeRadius(), giveAnodeRadius() - 0.0001);
    NPS += tempNPS;
    printf("\t %E \n - Neutrals outwards outside cathode from Class I\n", tempNPS);
    tempNPS = NIntegration(*Sfn1_OutPlusPtr, giveCathodeRadius(), giveAnodeRadius() - 0.0001);
    NPS += tempNPS;
    printf("\t %E \n", tempNPS);

    return NPS;
}

