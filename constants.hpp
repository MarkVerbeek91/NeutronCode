/**
 * Declaration of constants and tables. Here the properties of a fusor are
 * descripted. Some phicical constants defined like charge of a proton and
 * all parameters needed for the cross section calculation.
 *
 * The tables will contain the values that take a long time to calculate.
 *
 */

 #ifndef CONSTANTS_H
 #define CONSTANTS_H

struct Fusor{
    double a = 0.05;
    double b = 0.25;
    double V0 = -40000; // voltage
    double wire_diameter = 0.005;
    double Tc = 0.95;
} fusor;

// the precision of the functions.
#define N_PRECISION 250
#define N_TABLE     250

double q = 1.602e-19;
double pressure = 0.5;  // Pa
double Tgas = 400; // K
double ngas = 9.05401e19; //6.022e23 * pressure / (8.314 * Tgas);
double E0 = 0.0001;          // reducing errors
double Itot = 0.1;
double EdgeIonFlux;

// some data storage is now needed because other wise the calculation becomes hugh
struct Tables{
    double R[N_TABLE];
    double A[N_TABLE];
    double K[N_TABLE][N_TABLE];
    double S_0[N_TABLE];
    double S_1[N_TABLE];
    double S_2[N_TABLE];
    double S_3[N_TABLE];
    double S_4[N_TABLE];
    double S_5[N_TABLE];
} Table;

#endif // CONSTANTS_H
