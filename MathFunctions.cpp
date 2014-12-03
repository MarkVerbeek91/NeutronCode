/**
 * This file contains a few pure mathematical functions like interpolation and
 * numerical integration. These numerical integration is in several forms
 * because in this code several different integrations are done.
 */

#include <math.h>

#include "constants.h"

#include "MathFunctions.h"

double interpolation(double r)
{
    double value = -1;
    int low = N_TABLE - 2, upper = N_TABLE - 1;

    for ( int i = 0; i < N_TABLE; i++)
    {
        if ( r < Table.R[i])
        {
            low = i - 1;
            upper = i;
        }
    }

    double slope = ( Table.S[upper] - Table.S[low])/(Table.R[upper] - Table.R[low]);

    value = slope*r + Table.S[low] - Table.R[low] * slope;

    if (value < 0)
        value *= -1;

    return value;
}

double adaptiveSimpsonsAux(double (*f)(double), double a, double b, double epsilon,
                         double S, double fa, double fb, double fc, int bottom) {
  double c = (a + b)/2, h = b - a;
  double d = (a + c)/2, e = (c + b)/2;
  double fd = f(d), fe = f(e);
  double Sleft = (h/12)*(fa + 4*fd + fc);
  double Sright = (h/12)*(fc + 4*fe + fb);
  double S2 = Sleft + Sright;
  if (bottom <= 0 || fabs(S2 - S) <= 15*epsilon)
    return S2 + (S2 - S)/15;
  return adaptiveSimpsonsAux(f, a, c, epsilon/2, Sleft,  fa, fc, fd, bottom-1) +
         adaptiveSimpsonsAux(f, c, b, epsilon/2, Sright, fc, fb, fe, bottom-1);
}

double NIntegration( double (*f)(double), double a, double b)
{
    double epsilon = 1.0 / PRECISION;

    double c = (a + b)/2, h = b - a;
    double fa = f(a), fb = f(b), fc = f(c);
    double S = (h/6)*(fa + 4*fc + fb);

    return adaptiveSimpsonsAux(f, a, b, epsilon, S, fa, fb, fc, MAX_RECURSION_DEPTH);
}

double NIntegration_2( double (*funcPtr)(double, double), double Bar, double Start, double End)
{
    double sum = (*funcPtr)(Bar,Start) + (*funcPtr)(Bar,End), step = (End - Start)/N_PRECISION;

    for (double r=Start; r<End; r += step)
    {
        sum += 2.0 * (*funcPtr)(Bar, r);
    }

    sum = sum * step / 2.0;

    return sum;
}

// special integrator for the Class II ions.
double NIntegration_3( double (*funcPtr)(double, double), double Start, double End)
{
    double sum = (*funcPtr)(Start,End) + (*funcPtr)(End,End), step = (End - Start)/N_PRECISION;

    for (double r=Start; r<End; r += step)
    {
        sum += 2.0 * (*funcPtr)(r, End);
    }

    sum = sum * step / 2.0;

    return sum;
}

// this function gives of the slope of the function at point.
double differentiat( double (*funcPtr)(double), double point)
{
    double value;
    double step = 1.0 / N_PRECISION;

    value = ((*funcPtr)(point + step) + (*funcPtr)(point)) / step;

    return value;
}
