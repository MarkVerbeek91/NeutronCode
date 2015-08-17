/**
 *  NeutronCode written by Mark Verbeek, mark(dot)verbeek91(at)gmail(dot)com.
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

// standaard libs
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>     /* getenv */

#include "includes.h"

int main()
{
    FILE * input;
    input = fopen("input.ini","r");

    initBool();

    if ( input == NULL)
    {
        printf("# Missing input file, using standard parameters\n\n");
        init();
    }
    else
    {
        printf("# Reading input file:\n\n");
        readfile(input);
    }

    fclose(input);

    // filling the potential array and particle energy
    printf("# -- Start of program -- \n");
    // initialise the fusor parameters
//    init();

    // writing the potential to a file for plotting

    if ( printbool.potential )
    {
        printf("# Potential calculation\n");

        double (*Potential_PhiPtr)(double);
        Potential_PhiPtr = &Potential_Phi;

        double (*ParticleEnergy2Ptr)(double,double);
        ParticleEnergy2Ptr = &ParticleEnergy2;

//        print_data_ddd(*ParticleEnergy2Ptr, 0.0, giveAnodeRadius()+0.001, 0.01, 0.0, "Particle2.csv", 0);

          plot_function_dd(*Potential_PhiPtr, 0.0, 0.25, 0.001, "Potential.csv", "GNU_potential.txt");
//        print_data_dd(*Potential_PhiPtr, 0.0, 0.25, 0.001, "Potential.csv", 0);
    }

    // writing the SIIEE to a file for plotting
    if ( printbool.SIIEE )
    {
        printf("# SIIEE calculation\n");

        double (*SIIEEPtr)(double);
        SIIEEPtr = &SIIEE;

        //print_data_dd(*SIIEEPtr, 1.0, -giveVoltage(), 1, "SIIEE.csv", 1);
        plot_function_dd(*SIIEEPtr, 1.0, -giveVoltage(), 1, "SIIEE.csv", "GNU_SIIEE.txt");
    }

    // writing the Cross sections for Charge Exchange, Iononisation and the sum
    // of those to files to a file for plotting
    if ( printbool.Cross_section )
    {
        printf("# Cross section calculation\n");

        double (*CrosssecCXPtr)(double);
        CrosssecCXPtr = &CrosssecCX;

        plot_function_dd(*CrosssecCXPtr, 1.0, 500000, 10, "CrosssecCX.csv", "GNU_Cross_sections.txt");

        double (*CrosssecIonPtr)(double);
        CrosssecIonPtr = &CrosssecIon;

        plot_function_dd(*CrosssecIonPtr, 1.0, 500000, 10, "CrosssecIon.csv", "GNU_Cross_sections.txt");

        double (*CrosssecTotPtr)(double);
        CrosssecTotPtr = &CrosssecTot;

        plot_function_dd(*CrosssecTotPtr, 1.0, 500000, 10, "CrosssecTot.csv", "GNU_Cross_sections.txt");
    }

    // Writing the survival functions to a file for plotting.
    if ( printbool.Survival )
    {
        printf("# Survival function calculation\n");

        double (*fPtr)(double);
        fPtr = &f;

        plot_function_dd(*fPtr, 0.0, giveAnodeRadius()+0.001, 0.001, "f.csv", "GNU_Survival_funcitons.txt");

        double (*gPtr)(double,double);
        gPtr = &g;

        plot_function_ddd(*gPtr, 0.0, 0.00001, giveAnodeRadius(), 0.001, "g.csv", "GNU_Survival_funcitons.txt");
    }


    // writing A to a file for plotting
    if ( printbool.Atable )
    {
        printf("# Doing some things with A\n");

        double (*APtr)(double);
        APtr = &A;

        plot_function_dd(*APtr, giveCathodeRadius(), giveAnodeRadius(), 0.001, "A.csv", "GNU_Atable.txt");
    }


    // source rate for first generation of Class II ions.

    printf("# filling tables:\n");
    kernel_to_table();

    printf("# Calculating S:\n");
    S();
    printf("#  - Done\n");

    // calculating the energy spectrum of ions

        // building the "Kernel"
    // TODO: replace this function by the table in memory to file. Just like
    // by the S table a few line down. This will save computation time because
    // the kernel does not to be calculated twice.
    if ( printbool.KernelTable )
    {
        printf("# Kernel\n");
        double (*KPtr)(double,double);
        KPtr = &kernel;

        plot_function_ddd(*KPtr, 0.0, 0.0, giveAnodeRadius()+0.001, 0.01, "K.csv", "GNU_Ktable.txt");

 //       plot_table_2D();

    }


    // print the S tables to screen
    if ( printbool.Stable )
    {
        printf("# Writing tables to files:\n");

        plot_table_1D("STable.csv", "GNU_Stable.txt");

    }



    double I1, I2, I3, I4;

    printf("# Calculating Currents:\n");
    I1 = I_c1();
    printf("# 1 current: %E\n",I1);
    I2 = I_c2();
    printf("# 2 current: %E\n",I2);
    I3 = I_c3();
    printf("# 3 current: %E\n",I3);
    I4 = I_c4();
    printf("# 4 current: %E\n",I4);

    EdgeIonFlux = 5e15; //(Itot - I2 - I3) / (I1 + I3);

    printf("# EdgeIonFlux: %E\n# - Done\n", EdgeIonFlux);

    if ( printbool.Spectrum )
    {
        printf("# Printing Energy spectrum to files\n");

        printf("# Energy spectrum of ions going inwards\n");

        double (*IonSpectrumInwardsPtr)(double, double);
        IonSpectrumInwardsPtr = &IonSpectrumInwards;

        double (*IonSpectrumOutwardsPtr)(double, double);
        IonSpectrumOutwardsPtr = &IonSpectrumOutwards;

//        double r = 0.06;

//        print_data_ddd(*IonSpectrumInwardsPtr , 10, -giveVoltage(),10,r,"IonSpectrumInwards.csv");
//        print_data_ddd(*IonSpectrumOutwardsPtr, 10, -giveVoltage(),10,r,"IonSpectrumOutwards.csv");

//        int j = 6;

        int j = 1;
        for ( double r = 0.01;  r < 0.25; r=r+0.01)
        {
            char filename1[100];
            sprintf( filename1, "IonSpectrumInwards%d.csv", j);

            char filename2[100];
            sprintf( filename2, "IonSpectrumOutwards%d.csv", j);

            print_data_ddd(*IonSpectrumInwardsPtr , 10, -giveVoltage(),10,r,filename1, 5);
            print_data_ddd(*IonSpectrumOutwardsPtr, 10, -giveVoltage(),10,r,filename2, 5);

            j++;
        }


    }

    // the neutron source rate
    if ( printbool.NSR )
    {
        printf("# Printing neutron source rate to file:\n");
        double (*FuncPtr)(double);

            // Ions
        // inwards
        FuncPtr = &NeutronsIonFluxInwards;

        printf("# Outside cathode, inwards\n");
        print_data_dd(*FuncPtr, 0.001, giveAnodeRadius(), 0.001, "NSR_1.csv", 10);

        // outwards
        FuncPtr = &NeutronsIonFluxOutwards;

        printf("# Outside cathode, outwards\n");
        print_data_dd(*FuncPtr, 0.001, giveAnodeRadius(), 0.001, "NSR_2.csv", 11);

            // Neutrals Class I
        // inwards
        FuncPtr = &NeutronsNeutralsClassIFluxInwards;

        printf("# In cathode, inwards\n");
        print_data_dd(*FuncPtr, 0.001, giveAnodeRadius(), 0.001, "NSR_3.csv", 12);

        // outwards
        FuncPtr = &NeutronsNeutralsClassIFluxOutwards;

        printf("# In cathode, outwards\n");
        print_data_dd(*FuncPtr, 0.001, giveAnodeRadius(), 0.001, "NSR_4.csv", 13);

            // Neutrals ClassII
        // inwards
        FuncPtr = &NeutronsNeutralsClassIIFluxInwards;

        printf("# In cathode, inwards\n");
  //      print_data_dd(*FuncPtr, 0.001, giveAnodeRadius(), 0.001, "NSR_5.csv", 14);

        // outwards
        FuncPtr = &NeutronsNeutralsClassIIFluxOutwards;

        printf("# In cathode, outwards\n");
  //      print_data_dd(*FuncPtr, 0.001, giveAnodeRadius(), 0.001, "NSR_6.csv", 15);

    }

    if ( printbool.NPR )
    {
        printf("# Calculating NPR:\n");
        double NPR = Nps();
        printf("# NPR: %E \n# - Done\n", NPR);
    }

    // program is done
    printf("# -- Done --");
    return 0;
}
