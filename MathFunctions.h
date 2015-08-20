#ifndef MATHFUNCTIONS_H_INCLUDED
#define MATHFUNCTIONS_H_INCLUDED

double NIntegration( double (*funcPtr)(double), double Start, double End);
double NIntegration_2( double (*funcPtr)(double, double), double var, double Start, double End);
double NIntegration_3( double (*funcPtr)(double, double, double), double var1, double var2, double Start, double End);

double interpolation(double *table, double r);

double differentiat( double (*funcPtr)(double), double point);

bool DELTA(double d);

#endif // MATHFUNCTIONS_H_INCLUDED
