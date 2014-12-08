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

/** This section of code uses the Adaptive Simpsons rule to calculate the
    integral. It uses a recursive way to minimize to error.

*/
double adaptiveSimpsonsAux(double (*f)(double), double a, double b, double epsilon,
                           double S, double fa, double fb, double fc, int bottom)
{
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

double NIntegration( double (*funcPtr)(double), double Start, double End)
{
    double epsilon = 1.0 / PRECISION;

    double Mid = (Start + End)/2, h = End - Start;
    double funcStart = funcPtr(Start), funcEnd= funcPtr(End), funcMid = funcPtr(Mid);
    double S = (h/6)*(funcStart + 4*funcMid + funcEnd);

    return adaptiveSimpsonsAux(funcPtr, Start, End, epsilon, S, funcStart, funcEnd, funcMid, MAX_RECURSION_DEPTH);
}

double adaptiveSimpsonsAux2(double (*f)(double, double), double a, double b, double bar, double epsilon,
                           double S, double fa, double fb, double fc, int bottom)
{
    double c = (a + b)/2, h = b - a;
    double d = (a + c)/2, e = (c + b)/2;
    double fd = f(bar, d), fe = f(bar, e);
    double Sleft = (h/12)*(fa + 4*fd + fc);
    double Sright = (h/12)*(fc + 4*fe + fb);
    double S2 = Sleft + Sright;
    if (bottom <= 0 || fabs(S2 - S) <= 15*epsilon)
        return S2 + (S2 - S)/15;

    return adaptiveSimpsonsAux2(f, a, c, bar, epsilon/2, Sleft,  fa, fc, fd, bottom-1) +
           adaptiveSimpsonsAux2(f, c, b, bar, epsilon/2, Sright, fc, fb, fe, bottom-1);
}

// special integrator for the Class II ions.
double NIntegration_2( double (*funcPtr)(double, double), double bar, double Start, double End)
{
    double epsilon = 1.0 / PRECISION;

    double Mid = (Start + End)/2, h = End - Start;
    double funcStart = funcPtr(bar, Start), funcEnd= funcPtr(bar, End), funcMid = funcPtr(bar, Mid);
    double S = (h/6)*(funcStart + 4*funcMid + funcEnd);

    return adaptiveSimpsonsAux2(funcPtr, Start, End, bar, epsilon, S, funcStart, funcEnd, funcMid, MAX_RECURSION_DEPTH);
}

double adaptiveSimpsonsAux3(double (*f)(double, double), double a, double b, double epsilon,
                           double S, double fa, double fb, double fc, int bottom)
{
    double c = (a + b)/2, h = b - a;
    double d = (a + c)/2, e = (c + b)/2;
    double fd = f(d, b), fe = f(e, b);
    double Sleft = (h/12)*(fa + 4*fd + fc);
    double Sright = (h/12)*(fc + 4*fe + fb);
    double S2 = Sleft + Sright;
    if (bottom <= 0 || fabs(S2 - S) <= 15*epsilon)
        return S2 + (S2 - S)/15;

    return adaptiveSimpsonsAux3(f, a, c, epsilon/2, Sleft,  fa, fc, fd, bottom-1) +
           adaptiveSimpsonsAux3(f, c, b, epsilon/2, Sright, fc, fb, fe, bottom-1);
}

// integrator for the Class II ions surfival function
double NIntegration_3( double (*funcPtr)(double, double), double a, double b)
{
//    double sum = (*funcPtr)(Start,End) + (*funcPtr)(End,End), step = (End - Start)/N_PRECISION;

    double epsilon = 1.0 / PRECISION;

    double c = (a + b)/2, h = b - a;
    double funcStart = funcPtr(a, b), funcEnd= funcPtr(b, b), funcMid = funcPtr(c, b);
    double S = (h/6)*(funcStart + 4*funcMid + funcEnd);

    return adaptiveSimpsonsAux3(funcPtr, a, b, epsilon, S, funcStart, funcEnd, funcMid, MAX_RECURSION_DEPTH);
}

// this function gives of the slope of the function at point.
double differentiat( double (*funcPtr)(double), double point)
{
    double value;
    double step = 1.0 / N_PRECISION;

    value = ((*funcPtr)(point + step) + (*funcPtr)(point)) / step;

    return value;
}
