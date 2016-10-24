/**
 * @file MathFunctions.c
 * @author Mark Verbeek
 * @date 30 december 2015
 * @brief This file contains a few pure mathematical functions like interpolation and
 * numerical integration. These numerical integration is in several forms
 * because in this code several different type of functions are needed to be
 * integrations are done.
 */

#include <math.h>
#include <stdio.h>

#include "constants.h"
#include "MathFunctions.h"

/**
 * @brief Interpolate linear the value in table for location r
 *
 */
double interpolation(double *table, double r)
{
    double value = -1;
    int low = N_TABLE - 2, upper = N_TABLE - 1;

    for ( int i = 0; i < N_TABLE; i++)
    {
        if ( r < Table->R[i])
        {
            low = i - 1;
            upper = i;
        }
    }

    double slope = ( table[upper] - table[low])/(Table->R[upper] - Table->R[low]);

    value = slope*r + table[low] - Table->R[low] * slope;

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

    if (S2 == NAN)
    {
        fprintf(stderr,"NAN detected while integrating, stopping integration");
        return NAN;
    }

    if (bottom <= 0 || fabs(S2 - S) <= 15*epsilon)
        return S2 + (S2 - S)/15;

    return adaptiveSimpsonsAux(f, a, c, epsilon/2, Sleft,  fa, fc, fd, bottom-1) +
           adaptiveSimpsonsAux(f, c, b, epsilon/2, Sright, fc, fb, fe, bottom-1);
}

// an integration function for functions like: X(r) over r
// TODO: make more rebust. no infinte loops posible and such.
double NIntegration( double (*funcPtr)(double), double Start, double End)
{
    double epsilon = 1.0 / PRECISION;

    double Mid = (Start + End)/2, h = End - Start;
    double funcStart = funcPtr(Start), funcEnd= funcPtr(End), funcMid = funcPtr(Mid);
    double S = (h/6)*(funcStart + 4*funcMid + funcEnd);

    return adaptiveSimpsonsAux(funcPtr, Start, End, epsilon, S, funcStart, funcEnd, funcMid, MAX_RECURSION_DEPTH);
}

double adaptiveSimpsonsAux2(double (*f)(double, double), double a, double b, double epsilon, double var,
                           double S, double fa, double fb, double fc, int bottom)
{
    double c = (a + b)/2, h = b - a;
    double d = (a + c)/2, e = (c + b)/2;
    double fd = f(var, d), fe = f(var, e);
    double Sleft = (h/12)*(fa + 4*fd + fc);
    double Sright = (h/12)*(fc + 4*fe + fb);
    double S2 = Sleft + Sright;

    if (S2 == NAN)
    {
        fprintf(stderr,"NAN detected while integrating, stopping integration");
        return NAN;
    }

    if (bottom <= 0 || fabs(S2 - S) <= 15*epsilon)
        return S2 + (S2 - S)/15;

    return adaptiveSimpsonsAux2(f, a, c, epsilon/2, var, Sleft,  fa, fc, fd, bottom-1) +
           adaptiveSimpsonsAux2(f, c, b, epsilon/2, var, Sright, fc, fb, fe, bottom-1);
}

// An integration function for functions like: X(r, dr) over dr
// TODO: make more robust
double NIntegration_2( double (*funcPtr)(double, double), double var, double Start, double End)
{
    double epsilon = 1.0 / PRECISION;

    double Mid = (Start + End)/2, h = End - Start;

    double funcStart = funcPtr(var, Start), funcEnd = funcPtr(var, End), funcMid = funcPtr(var, Mid);
    double S = (h/6)*(funcStart + 4*funcMid + funcEnd);

    return adaptiveSimpsonsAux2(funcPtr, Start, End, epsilon, var, S, funcStart, funcEnd, funcMid, MAX_RECURSION_DEPTH);
 }

double adaptiveSimpsonsAux3(double (*f)(double, double, double), double a, double b, double epsilon, double var1, double var2,
                           double S, double fa, double fb, double fc, int bottom)
{
    double c = (a + b)/2, h = b - a;
    double d = (a + c)/2, e = (c + b)/2;
    double fd = f(var1, var2, d), fe = f(var1, var2, e);
    double Sleft = (h/12)*(fa + 4*fd + fc);
    double Sright = (h/12)*(fc + 4*fe + fb);
    double S2 = Sleft + Sright;

    if (S2 == NAN)
    {
        fprintf(stderr,"NAN detected while integrating, stopping integration");
        return NAN;
    }

    if (bottom <= 0 || fabs(S2 - S) <= 15*epsilon)
        return S2 + (S2 - S)/15;

    return adaptiveSimpsonsAux3(f, a, c, epsilon/2, var1, var2, Sleft,  fa, fc, fd, bottom-1) +
           adaptiveSimpsonsAux3(f, c, b, epsilon/2, var1, var2,  Sright, fc, fb, fe, bottom-1);
}

// An integration function for functions like: X(var1, var2, ddr) over ddr
// TODO: make more robust
double NIntegration_3( double (*funcPtr)(double, double, double), double var1, double var2, double Start, double End)
{
    double epsilon = 1.0 / PRECISION;

    double Mid = (Start + End)/2, h = End - Start;

    double funcStart = funcPtr(var1, var2, Start), funcEnd = funcPtr(var1, var2, End), funcMid = funcPtr(var1, var2, Mid);
    double S = (h/6)*(funcStart + 4*funcMid + funcEnd);

    return adaptiveSimpsonsAux3(funcPtr, Start, End, epsilon, var1, var2, S, funcStart, funcEnd, funcMid, MAX_RECURSION_DEPTH);
 }

// this function gives of the slope of the function at point.
// TODO: make more robust
double differentiat( double (*funcPtr)(double), double point)
{
    double value;
    double step = 1.0 / N_PRECISION;

    value = ((*funcPtr)(point + step) - (*funcPtr)(point)) / step;

    return value;
}

bool DELTA(double d)
{
    if ( d < 0.0000001 )
        return 1;

    return 0;
}

