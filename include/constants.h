/**
 * Declaration of constants and tables. Here the properties of a fusor are
 * descripted. Some phicical constants defined like charge of a proton and
 * all parameters needed for the cross section calculation.
 *
 * The tables will contain the values that take a long time to calculate.
 *
 */

#include "stdbool.h"

#ifndef CONSTANTS_H_INCLUDED
#define CONSTANTS_H_INCLUDED

void init(void);
void initBool(void);
double giveq(void);
double giveTgas(void);
double giveItot(void);
double giveEdgeIonFlux(void);

struct Fusor{
    double a;				// Cathode radius
    double b;				// Anode radius
    double V0; 				// Voltage on Fusor
//    double wire_diameter;	//
    double Tc;				// Transparency of grid
    double ngas;            // gas density

};

// the precision of the functions.
#define N_PRECISION 100			// number of step in differencation function
#define PRECISION 100000		// max allowed error in integration functions
#define MAX_RECURSION_DEPTH 15	// max depth of recusing integration functions
#define N_TABLE     250			// size of Tables, Kernel is value^2
#define Q_ELECTRON 1.602e-19   	// coulombs

// declaration of some variables
extern double q;			// unit
extern double pressure;  	// Pa
extern double Tgas; 		// K
extern double ngas; 		// #
extern double E0;          	// eV
extern double Itot;			// A
extern double EdgeIonFlux;	// 1/m^2

// some data storage is now needed because other wise the calculation becomes hugh
struct Tables{
    double R[N_TABLE];			// Table filled with radius location of other tables
    double A[N_TABLE];			// ions source rate table
    double K[N_TABLE][N_TABLE];	// kernel for neutral source rate table
    double S[N_TABLE];			// neutral source rate table
};

// Struct for bools of which parameters should be printen to screen and file
// TODO: separate printing to screen and file.
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

extern struct Fusor FusorData, *fusor;
extern struct Tables TableData, *Table;
extern struct PrintBool Printbool, *printbool;
extern struct PrintBool Printbool2, *printbool2;

#endif // CONSTANTS_H_INCLUDED
