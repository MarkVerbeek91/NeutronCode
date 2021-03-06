/**
 * Declaration of constants and tables. Here the properties of a fusor are
 * descripted. Some phicical constants defined like charge of a proton and
 * all parameters needed for the cross section calculation.
 *
 * The tables will contain the values that take a long time to calculate.
 *
 */
#include "stdbool.h"
#include "stdio.h"

#include "constants.h"

void init(void)
{
    fusor->a = 0.05;
    fusor->b = 0.25;
    fusor->V0 = -40000;
//    fusor->wire_diameter = 0.005;
    fusor->Tc = 0.95;
	q = 1;
	pressure = 0.5;
	Itot = 0.1;
	Tgas = 400;
    fusor->ngas = 6.022e23 * pressure / (8.314 * Tgas); //
}

void initBool(void)
{
    printbool->potential     = false;
    printbool->SIIEE         = false;
    printbool->Cross_section = false;
    printbool->Survival      = false;
    printbool->Atable        = false;
    printbool->KernelTable   = false;
    printbool->Stable        = false;
    printbool->Spectrum      = false;
    printbool->NSR           = false;
    printbool->NPR           = false;

	printbool2->potential     = false;
    printbool2->SIIEE         = false;
    printbool2->Cross_section = false;
    printbool2->Survival      = false;
    printbool2->Atable        = false;
    printbool2->KernelTable   = false;
    printbool2->Stable        = false;
    printbool2->Spectrum      = false;
    printbool2->NSR           = false;
    printbool2->NPR           = false;
}

//Fusor fusor;
//Tables Table;
//PrintBool printbool;
struct Fusor FusorData;
struct Fusor *fusor = &FusorData;
struct Tables TableData;
struct Tables *Table =&TableData;
struct PrintBool Printbool;
struct PrintBool *printbool = &Printbool;
struct PrintBool Printbool2;
struct PrintBool *printbool2 = &Printbool2;

// defining of some variables
double q;
double E0 = 0.0001;          // reducing errors

double pressure;  // Pa
double Tgas; // K
double Itot;
double EdgeIonFlux;

