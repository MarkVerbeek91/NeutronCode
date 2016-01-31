/**
    The function intergrate the particle energy spectra to give a fussion
    rate on every radii. Intergrating these function should give a total
    neutron flux.
*/

#include "constants.h"
#include "MathFunctions.h"
#include "IonSpectrum.h"
#include "CrossSections.h"
#include "PotentialFunctions.h"

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
    double FRR, Emax;
    double (*funcPtr)(double, double);

    funcPtr = &FusionRateIons_Inte;
    Emax = ParticleEnergy1(r);
    FRR  = NIntegration_2(funcPtr, r, 1e-3, Emax);

    FRR *= ngas;
    return FRR;
}
