#ifndef CROSSSECTIONS_H_INCLUDED
#define CROSSSECTIONS_H_INCLUDED

double CrosssecCX(double E);
double CrosssecIon(double E);
double CrosssecTot(double energy);
double CrosssecFusion(double E);
void InitCrossSectionConstands(void);

#endif // CROSSSECTIONS_HPP_INCLUDED
