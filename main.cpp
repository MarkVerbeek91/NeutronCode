/**
 *  NeutronCode written by Mark Verbeek, mark.verbeek91(at)gmail(dot)com.
 *
 *  This code calculated the neutron production rate, NPR, of a fusor devise by
 *  using analitical descriptions of the physics in the fusor.
 *
 *  This analithical description is made by Gilbert A. Emmert and
 *  John F. Santarius in the article:
 *      - Atomic and molecular effects on spherically convergent ion flow. I.
 *        Single atomic species, AIP (2010)
 *
 *  \author Mark Verbeek
 *
 */


#include <iostream>
#include <math>
#include <stdio>

#include "constants.hpp"
#include "functions.hpp"
#include "CrossSections.hpp"

#include "MathFunctions.cpp"
#include "PotentialFunctions.cpp"
#include "CrossSections.cpp"
#include "SurvivalFunctions.cpp"
#include "ClassIISourceFunction.cpp"
#include "CathodeCurrents.cpp"

#include "ParticleFlux.cpp"
#include "NeutronProductionRate.cpp"
#include "Auxiliary.cpp"

int main()
{
    // filling the potential array and particle energy
    printf("-- Start of program -- \n");

    // writing the potential to a file for plotting

    if ( false )
    {
        printf("Potential calculation\n");

        double (*Potential_PhiPtr)(double);
        Potential_PhiPtr = &Potential_Phi;

        double (*ParticleEnergy2Ptr)(double,double);
        ParticleEnergy2Ptr = &ParticleEnergy2;

        print_data_ddd(*ParticleEnergy2Ptr, 0.0, fusor.b+0.001, 0.01, 0.0, "Particle2.csv");

        print_data_dd(*Potential_PhiPtr, 0.0, 0.25, 0.001, "Potential.csv");
    }

    // writing the SIIEE to a file for plotting
    if ( false )
    {
        printf("SIIEE calculation\n");

        double (*SIIEEPtr)(double);
        SIIEEPtr = &SIIEE;

        print_data_dd(*SIIEEPtr, 1.0, -fusor.V0, 1, "SIIEE.csv");
    }

    // writing the Cross sections for Charge Exchange, Iononisation and the sum
    // of those to files to a file for plotting
    if ( false )
    {
        printf("Cross section calculation\n");

        double (*CrosssecCXPtr)(double);
        CrosssecCXPtr = &CrosssecCX;

        print_data_dd(*CrosssecCXPtr, 1.0, 500000, 10, "CrosssecCX.csv");

        double (*CrosssecIonPtr)(double);
        CrosssecIonPtr = &CrosssecIon;

        print_data_dd(*CrosssecIonPtr, 1.0, 500000, 10, "CrosssecIon.csv");

        double (*CrosssecTotPtr)(double);
        CrosssecTotPtr = &CrosssecTot;

        print_data_dd(*CrosssecTotPtr, 1.0, 500000, 10, "CrosssecTot.csv");
    }


    // Writing the survival functions to a file for plotting.
    if ( false )
    {
        printf("Survival function calculation\n");

        double (*fPtr)(double);
        fPtr = &f;

        print_data_dd(*fPtr, 0.0, fusor.b+0.001, 0.001, "f.csv");

        double (*gPtr)(double,double);
        gPtr = &g;

        print_data_ddd(*gPtr, 0.0, fusor.b+0.001, 0.001, 0, "g.csv");
    }


    // writing A to a file for plotting
    if ( false )
    {
        printf("doing some thing with A\n");

        double (*APtr)(double);
        APtr = &A;

        print_data_dd(*APtr, fusor.a, fusor.b+0.01, 0.001, "A.csv");
    }

    // building the "Kernel"
    if ( false )
    {
        printf("Kernel\n");
        double (*KPtr)(double,double);
        KPtr = &kernel;

        print_data_ddd(*KPtr, 0.0, fusor.b+0.001, 0.01, 0.0001, "K.csv");
    }

    // source rate for first generation of Class II ions.

    printf("filling tables:\n");
    kernel_to_table();

    printf("Calculating S:\n");
    S();
    printf(" - Done\n");

    // print the S tables to screen
    if ( false )
    {
        printf("Writing tables to files:\n");

        print_table(1, "Atable.csv");
        print_table(2, "Ktable.csv");
        print_table(10, "S0.csv");
        print_table(11, "S1.csv");
        print_table(12, "S2.csv");
        print_table(13, "S3.csv");
        print_table(14, "S4.csv");
        print_table(15, "S5.csv");

    }

    double TotalCurrent = 0, dummie;

    printf("Calculating Currents:\n");
    dummie = I_c1();
    printf("total current: %E\n",dummie);
    TotalCurrent += dummie;
    dummie = I_c2();
    printf("total current: %E\n",dummie);
    TotalCurrent += dummie;
    dummie = I_c3();
    printf("total current: %E\n",dummie);
    TotalCurrent += dummie;
    dummie = I_c4();
    printf("total current: %E\n",dummie);
    TotalCurrent += dummie;

    EdgeIonFlux = Itot / TotalCurrent;

    printf("total current: %E, \n\n EdgeIonFlux: %E\n - Done\n",TotalCurrent, EdgeIonFlux);


    if ( false )
    {
        printf("Printing neutron source rate to file:\n");

        // outside the cathode
        double (*Sfi_OutMinPtr)(double);
        Sfi_OutMinPtr = &Sfi_OutMin;

        printf("Outside cathode, inwards\n");
        print_data_dd(*Sfi_OutMinPtr, fusor.a, fusor.b, 0.001, "Sfi_OutMin.csv");

        double (*Sfi_OutPlusPtr)(double);
        Sfi_OutPlusPtr = &Sfi_OutPlus;

        printf("Outside cathode, outwards\n");
        print_data_dd(*Sfi_OutPlusPtr, fusor.a, fusor.b, 0.001, "Sfi_OutPlus.csv");

        // inside the cathode
        double (*Sfi_InMinPtr)(double);
        Sfi_InMinPtr = &Sfi_InMin;

        printf("In cathode, inwards\n");
        print_data_dd(*Sfi_InMinPtr, 0, fusor.a, 0.001, "Sfi_InMin.csv");

        double (*Sfi_InPlusPtr)(double);
        Sfi_InPlusPtr = &Sfi_InPlus;

        printf("In cathode, outwards\n");
        print_data_dd(*Sfi_InPlusPtr, 0, fusor.a, 0.001, "Sfi_InPlus.csv");

    }

    if ( true )
    {
        printf("Calculating NPR form fast Ions:\n");
        double NPR = Nps();
        printf("NPR: %E \n - Done\n", NPR);

    }

    if ( false )
    {
        printf("Calculation NPR from fast neutrals:\n");


    }

    // program is done


    printf("\n-- Done --");
    return 0;
}
