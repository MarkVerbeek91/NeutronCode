/**
 * Declaration of constants and tables. Here the properties of a fusor are
 * descripted. Some phicical constants defined like charge of a proton and
 * all parameters needed for the cross section calculation.
 *
 * The tables will contain the values that take a long time to calculate.
 *
 */

#ifndef CONSTANTS_H_INCLUDED
#define CONSTANTS_H_INCLUDED

void init(void);
void initBool(void);
double giveCathodeRadius(void);
double giveAnodeRadius(void);
double giveVoltage(void);
double giveTransparency(void);
double giveq(void);
double giveTgas(void);
double giveNgas(void);
double giveItot(void);
double giveEdgeIonFlux(void);

struct Fusor{
    double a;
    double b;
    double V0; // voltage
    double wire_diameter;
    double Tc;
};

// the precision of the functions.
#define N_PRECISION 100
#define PRECISION 100000
#define MAX_RECURSION_DEPTH 5
#define N_TABLE     250
#define Q_ELECTRON 1.602e-19

// declaration of some variables
extern double q;
extern double pressure;  // Pa
extern double Tgas; // K
extern double ngas; //6.022e23 * pressure / (8.314 * Tgas);
extern double E0;          // reducing errors
extern double Itot;
extern double EdgeIonFlux;

// some data storage is now needed because other wise the calculation becomes hugh
struct Tables{
    double R[N_TABLE];
    double A[N_TABLE];
    double K[N_TABLE][N_TABLE];
    double S[N_TABLE];
};

struct PrintBool{
    bool potential;
    bool SIIEE;
    bool Cross_section;
    bool Survival;
    bool Atable;
    bool KernelTable;
    bool Stable;
    bool Spectrum;
    bool NSR;
    bool NPR;
};

extern Fusor fusor;
extern Tables Table;
extern PrintBool printbool;

#endif // CONSTANTS_H_INCLUDED
