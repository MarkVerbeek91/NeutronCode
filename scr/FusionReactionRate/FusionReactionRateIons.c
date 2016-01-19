/**
    The function intergrate the particle energy spectra to give a fussion
    rate on every radii. Intergrating these function should give a total
    neutron flux.
*/

#include "constants.h"
#include "MathFunctions.h"
#include "IonSpectrum.h"
#include "CrossSections.h"

#include "FusionReactionRateIons.h"

double FusionRateIons_Inte(double r, double E)
{
    double FRR;

    FRR  = IonSpectrumInwards(r, E);
    FRR += IonSpectrumOutwards(r, E);
    FRR *= CrosssecFusion(E);

    return FRR;
}

double FusionRateIons(double r)
{
    double FRR;
    double (*funcPtr)(double, double);

    funcPtr = &FusionRateIons_Inte;
    FRR  = NIntegration_2(funcPtr, r, 0, -1*fusor->V0);

    FRR *= ngas;
    return FRR;
}
