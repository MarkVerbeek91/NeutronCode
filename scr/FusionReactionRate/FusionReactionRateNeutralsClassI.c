/**
    The function intergrate the particle energy spectra to give a fussion
    rate on every radii. Intergrating these function should give a total
    neutron flux.
*/

#include "constants.h"
#include "MathFunctions.h"
#include "NeutralsClassISpectrum.h"
#include "CrossSections.h"
#include "PotentialFunctions.h"

#include "FusionReactionRateNeutralsClassI.h"

double FusionRateNeutralsClassI_Inte(double r, double E)
{
    double FRR;

    FRR  = NeutralsClassISpectrumInwards(r, E);
    FRR += NeutralsClassISpectrumOutwards(r, E);
    FRR *= CrosssecFusion(E);

    return FRR;
}

double FusionRateNeutralsClassI(double r)
{
    double FRR, Emax;
    double (*funcPtr)(double, double);

    funcPtr = &FusionRateNeutralsClassI_Inte;
    Emax = ParticleEnergy1(r);
    FRR  = NIntegration_2(funcPtr, r, 1e-3, Emax);

    FRR *= ngas;
    return FRR;
}
