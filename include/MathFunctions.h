/**
 * @file MathFunctions.c
 * @author Mark Verbeek
 * @date 30 december 2015
 * @brief This file contains a few pure mathematical functions like interpolation and
 * numerical integration. These numerical integration is in several forms
 * because in this code several different type of functions are needed to be 
 * integrations are done.
 */
#ifndef MATHFUNCTIONS_H_INCLUDED
#define MATHFUNCTIONS_H_INCLUDED

/**
 * @brief Interpolate linear the value in table for location r
 *
 */
double interpolation(double *table, double r);

double NIntegration( double (*funcPtr)(double), double Start, double End);
double NIntegration_2( double (*funcPtr)(double, double), double var, double Start, double End);
double NIntegration_3( double (*funcPtr)(double, double, double), double var1, double var2, double Start, double End);


double differentiat( double (*funcPtr)(double), double point);

bool DELTA(double d);

#endif // MATHFUNCTIONS_H_INCLUDED
