/**
    The function intergrate the particle energy spectra to give a fussion
    rate on every radii. Intergrating these function should give a total
    neutron flux.
*/

#include "constants.h"
#include "MathFunctions.h"
#include "NeutralsClassIISpectrum.h"
#include "CrossSections.h"

#include "FusionReactionRateNeutralsClassII.h"

double FusionRateNeutralsClassII_Inte(double r, double E)
{
    double FRR;

    FRR  = NeutralsClassIISpectrumInwards(r, E);
    FRR += NeutralsClassIISpectrumOutwards(r, E);
    FRR *= CrosssecFusion(E);

    return FRR;
}

double FusionRateNeutralsClassII(double r)
{
    double FRR;
    double (*funcPtr)(double, double);

    funcPtr = &FusionRateNeutralsClassII_Inte;
    FRR  = NIntegration_2(funcPtr, r, 0, fusor->V0);

    FRR *= ngas;
    return FRR;
}
