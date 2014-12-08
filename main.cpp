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

// standaard libs
#include <iostream>
#include <math.h>
#include <stdio.h>

#include "includes.h"

int main()
{
    FILE * input;
    fopen("input.txt","r");

    if ( input == NULL)
    {
        printf("missing input file, using standard parameters\n\n");
        init();
    }
    else
    {
        printf("Reading input file:\n\n");

        readfile(&input);
    }


    // filling the potential array and particle energy
    printf("-- Start of program -- \n");
    // initialise the fusor parameters
    init();

    // writing the potential to a file for plotting

    if ( false )
    {
        printf("Potential calculation\n");

        double (*Potential_PhiPtr)(double);
        Potential_PhiPtr = &Potential_Phi;

        double (*ParticleEnergy2Ptr)(double,double);
        ParticleEnergy2Ptr = &ParticleEnergy2;

        print_data_ddd(*ParticleEnergy2Ptr, 0.0, giveAnodeRadius()+0.001, 0.01, 0.0, "Particle2.csv");

        print_data_dd(*Potential_PhiPtr, 0.0, 0.25, 0.001, "Potential.csv");
    }

    // writing the SIIEE to a file for plotting
    if ( false )
    {
        printf("SIIEE calculation\n");

        double (*SIIEEPtr)(double);
        SIIEEPtr = &SIIEE;

        print_data_dd(*SIIEEPtr, 1.0, -giveVoltage(), 1, "SIIEE.csv");
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
    if ( true )
    {
        printf("Survival function calculation\n");

        double (*fPtr)(double);
        fPtr = &f;

        print_data_dd(*fPtr, 0.0, giveAnodeRadius()+0.001, 0.001, "f.csv");

        double (*gPtr)(double,double);
        gPtr = &g;

        print_data_ddd(*gPtr, 0.0, giveAnodeRadius()+0.001, 0.001, 0, "g.csv");
    }


    // writing A to a file for plotting
    if ( false )
    {
        printf("doing some thing with A\n");

        double (*APtr)(double);
        APtr = &A;

        print_data_dd(*APtr, giveCathodeRadius(), giveAnodeRadius(), 0.001, "A.csv");
    }

    // building the "Kernel"
    if ( false )
    {
        printf("Kernel\n");
        double (*KPtr)(double,double);
        KPtr = &kernel;

        print_data_ddd(*KPtr, 0.0, giveAnodeRadius()+0.001, 0.01, 0.0001, "K.csv");
    }

    // source rate for first generation of Class II ions.

    printf("filling tables:\n");
    kernel_to_table();

    printf("Calculating S:\n");
    S();
    printf(" - Done\n");

    // calculating the energy spectrum of ions


    // print the S tables to screen
    if ( true )
    {
        printf("Writing tables to files:\n");

        print_table(1, "Atable.csv");
        print_table(2, "Ktable.csv");
        print_table(9, "S.csv");
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
        printf("Printing Energy spectrum to files\n");

        printf("Energy spectrum of ions going inwards\n");

        double (*f_minPtr)(double, double);
        f_minPtr = &f_min;

        double (*f_plusPtr)(double, double);
        f_plusPtr = &f_plus;

        int j = 6;

        for ( double i = 0.06; i < 0.25; i=i+0.01, j++)
        {
            char filename1[100];
            sprintf( filename1, "f_min_%d.csv", j);

            char filename2[100];
            sprintf( filename2, "f_plus_%d.csv", j);

            print_data_ddd( *f_minPtr, 10, -giveVoltage(),10,i,filename1);
            print_data_ddd(*f_plusPtr, 10, -giveVoltage(),10,i,filename2);

        }


/*
        for ( int i = 0; i < 10000; i=i+100)
        {
            printf("%d, %E\n",i, f_min(0.06,i));
        }
  */
    }

    if ( true )
    {
        printf("Printing neutron source rate to file:\n");

        // outside the cathode
        double (*Sfi_OutMinPtr)(double);
        Sfi_OutMinPtr = &Sfi_OutMin;

        printf("Outside cathode, inwards\n");
        print_data_dd(*Sfi_OutMinPtr, giveCathodeRadius(), giveAnodeRadius(), 0.001, "Sfi_OutMin.csv");

        double (*Sfi_OutPlusPtr)(double);
        Sfi_OutPlusPtr = &Sfi_OutPlus;

        printf("Outside cathode, outwards\n");
        print_data_dd(*Sfi_OutPlusPtr, giveCathodeRadius(), giveAnodeRadius(), 0.001, "Sfi_OutPlus.csv");

        // inside the cathode
        double (*Sfi_InMinPtr)(double);
        Sfi_InMinPtr = &Sfi_InMin;

        printf("In cathode, inwards\n");
        print_data_dd(*Sfi_InMinPtr, 0, giveCathodeRadius(), 0.001, "Sfi_InMin.csv");

        double (*Sfi_InPlusPtr)(double);
        Sfi_InPlusPtr = &Sfi_InPlus;

        printf("In cathode, outwards\n");
        print_data_dd(*Sfi_InPlusPtr, 0, giveCathodeRadius(), 0.001, "Sfi_InPlus.csv");

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
