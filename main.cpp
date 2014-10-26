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
 */

#include <iostream>
#include <math.h>
#include <stdio.h>

#include "constants.hpp"
#include "functions.hpp"
#include "math_functions.cpp"
#include "CrossSections.cpp"
#include "CathodeCurrents.cpp"

Fusor fusor;

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

/** \brief This function returns the potential on position r. Input is a double r in
 *  meter and return is a double phi in KiloVolt.
 *
 * \param r = radius
 * \return potential in V
 */
double Potential_Phi(double r)
{
    double phi;

    if ( r <= fusor.a)
        phi = fusor.V0;
    else
        phi = (fusor.a * (fusor.b - r) * -1 * fusor.V0) / (r * (fusor.a - fusor.b));

    return phi;
}

/**
    In E in eV
*/
double SIIEE(double energy)
{
    double tmp;
    tmp = 1.5*pow((1.15*(energy/97861)),-0.667)*(1-exp(-1.8*pow(energy/97891,1.2)));
    return tmp;
}

/** \brief calculates the energy of particle at radius r coming from outer edge
 *
 * \param r = radius where particle is
 * \return energie in eV
 *
 */
double ParticleEnergy1(double r)
{
    double energy;
    energy = - Potential_Phi(r) + 0.001;
    return energy;
}

/** \brief Calculates the energy of particle at radius r coming from born radius r1
 *
 * \param r = radius where particle is
 * \param r1= radius where particle is born
 * \return energy in eV.
 *
 */
double ParticleEnergy2(double r, double r1)
{
    double energy;
    energy = (Potential_Phi(r1) - Potential_Phi(r)) + 4;
    return energy;
}

/**
    The next function compute the chance that a particle survise to that radius
*/
double f_inte(double r)
{
    return CrosssecCX(ParticleEnergy1(r));
}

double f(double r)
{
    double (*f_intePtr)(double);
    f_intePtr = &f_inte;

    double sum = NIntegration(*f_intePtr, r, fusor.b);

    sum *= ngas;

    sum = exp(-1 * sum);

    return sum;
}

double Intergrant(double r)
{
    double tmp;

    tmp = ngas * CrosssecCX(ParticleEnergy1(r));

    return tmp;

}

double A(double r)
{
    double tmp = -1;
    tmp = ngas * CrosssecTot(ParticleEnergy1(r)) * gamma(r);

    if (isnan(tmp))
    {
        printf("went to the end\n");
        tmp = 0;
    }
    return tmp;
}

double gamma(double r)
{
    double Gamma = -1;
    Gamma = pow(fusor.b/r,2) * (f(r) + pow(fusor.Tc * f(0),2)/f(r));
    return Gamma;
}

double g_inte(double r, double dr)
{
    return CrosssecCX(ParticleEnergy2(r,dr));
}

double g(double r, double r1)
{
    if ( r > r1)
    {
        printf("error: r >= r1 \n");
        return -2;
    }

    double (*g_intePtr)(double, double);
    g_intePtr = &g_inte;

    double sum = NIntegration_3(*g_intePtr, r, r1);

    sum *= ngas;
    sum = exp(-sum);

    return sum;
}

double kernel(double r, double r1)
{
    double tmp;
    if ( r < r1 )
    {
        tmp = pow(r1/r,2) * ((g(r,r1) + (pow(fusor.Tc*g(0,r1),2)/g(r,r1)))/(1.0-pow(fusor.Tc*g(0,r1),2)));
        tmp = tmp * ngas * CrosssecTot(ParticleEnergy2(r,r1));
    }
    else
        tmp = 0;

    return tmp;
}

void kernel_to_table(void)
{
    double step = (fusor.b - fusor.a)/(N_TABLE-1);

    for ( int i = 0; i < N_TABLE; i++)
    {
        for ( int j = 0; j < N_TABLE; j++)
        {
            Table.K[i][j] = kernel(i*step+fusor.a,j*step+fusor.a);
        }

        if ( i *(N_TABLE / 50) % 20 == 0)
            printf(".");

        Table.A[i] = A(i*step+fusor.a);

        Table.R[i] = i*step+fusor.a;
    }

    Table.A[N_TABLE-1] = 43.9944;

    printf("\n");

    return;
}


/**
    Four function are needed to calculate the neutron production at given radius.

    a function for ingoing ions in the cathode
    a function for outgoing ions in the cathode
    a function for ingoing ions outside the cathode
    a function for outgoing ions outside the cathode
*/

double Sfi_OutMinInte(double r, double dr)
{
    double fac  = CrosssecFusion(ParticleEnergy2(r,dr)) * pow(dr,2) * interpolation(r);
           fac *= g(r,dr)/( 1 - pow(fusor.Tc*g(0,dr),2) );

    return fac;
}

// outside cathode inward.
double Sfi_OutMin(double r)
{
    double S;
    double term1, term2;

    double (*Sfi_OutMinIntePtr)(double, double);
    Sfi_OutMinIntePtr = &Sfi_OutMinInte;

    term1 = NIntegration_2(*Sfi_OutMinIntePtr, r, r, fusor.b);
    term1 *= ngas * EdgeIonFlux;

    term2 = ngas * pow(fusor.b/r,2) * EdgeIonFlux * f(r) * CrosssecFusion(ParticleEnergy1(r));

    S = term1 + term2;

    return S;
}

// outside cathode outward.
double Sfi_OutPlusInte(double r, double dr)
{
    double fac  = CrosssecFusion(ParticleEnergy2(r,dr)) * pow(dr,2) * interpolation(r);
           fac *= (pow(fusor.Tc * g(0,dr),2)/( 1 - pow(fusor.Tc*g(0,dr),2) )) * 1 / g(r,dr);

    return fac;
}

double Sfi_OutPlus(double r)
{
    double S;
    double term1, term2;

    double (*Sfi_OutPlusIntePtr)(double, double);
    Sfi_OutPlusIntePtr = &Sfi_OutPlusInte;

    term1 = NIntegration_2(*Sfi_OutPlusIntePtr, r, r, fusor.b);
    term1 *= ngas * EdgeIonFlux;

    term2 = ngas * pow(fusor.b/r,2) * EdgeIonFlux * (pow(f(0),2) / f(r)) * CrosssecFusion(ParticleEnergy1(r));

    S = term1 + term2;

    return S;
}

double Sfi_InMinInte(double r, double dr)
{
    double fac  = CrosssecFusion(ParticleEnergy2(r,dr)) * pow(dr,2) * interpolation(r);
           fac *= fusor.Tc * g(r,dr) / ( 1 - pow(fusor.Tc*g(0,dr),2) ) ;
           fac *= exp(ngas * CrosssecCX(ParticleEnergy2(fusor.a,dr)) * ( r - fusor.a));

    return fac;
}

double Sfi_InMin(double r)
{
    double S;
    double term1, term2;

    double (*Sfi_InMinIntePtr)(double, double);
    Sfi_InMinIntePtr = &Sfi_InMinInte;

    term1 = NIntegration_2(*Sfi_InMinIntePtr, r, fusor.a, fusor.b);
    term1 *= ngas * EdgeIonFlux;

    term2  = ngas * pow(fusor.b/r,2) * EdgeIonFlux * f(fusor.a) * CrosssecFusion(ParticleEnergy1(r));
    term2 *= exp(ngas * CrosssecCX(-fusor.V0) * ( r - fusor.a));

    S = term1 + term2;

    return S;
}

double Sfi_InPlusInte(double r, double dr)
{
    double fac = CrosssecFusion(ParticleEnergy2(r,dr)) * pow(dr,2) * interpolation(r);
           fac *= pow(fusor.Tc * g(0,dr),2) / ( 1 - pow(fusor.Tc*g(0,dr),2) ) * 1 / g(fusor.a,dr) ;
           fac *= exp( -1 * ngas * CrosssecCX(ParticleEnergy2(fusor.a,dr)) * ( r - fusor.a));

    return fac;
}

double Sfi_InPlus(double r)
{
    double S;
    double term1, term2;

    double (*Sfi_InPlusPtr)(double, double);
    Sfi_InPlusPtr = &Sfi_InPlusInte;

    term1 = NIntegration_2(*Sfi_InPlusPtr, r, fusor.a, fusor.b);
    term1 *= ngas * EdgeIonFlux;

    term2  = ngas * pow(fusor.b/r,2) * EdgeIonFlux * pow(f(fusor.a),2) * CrosssecFusion(ParticleEnergy1(r)) / f(fusor.a);
    term2 *= exp( -1 * ngas * CrosssecCX(-fusor.V0) * ( r - fusor.a));

    S = term1 + term2;

    return S;
}

/**
    To get the total Neutron production, the neutron production needs to be integrated
    over the volume.
*/
double Nps(void)
{
    double NPS;

    double (*Sfi_InMinPtr)(double);
    Sfi_InMinPtr = &Sfi_InMin;
    double (*Sfi_InPlusPtr)(double);
    Sfi_InPlusPtr = &Sfi_InPlus;
    double (*Sfi_OutMinPtr)(double);
    Sfi_OutMinPtr = &Sfi_OutMin;
    double (*Sfi_OutPlusPtr)(double);
    Sfi_OutPlusPtr = &Sfi_OutPlus;

    printf(" - Ions inwards inside cathode\n");
    NPS  = NIntegration(*Sfi_InMinPtr, 0.01, fusor.a - 0.000001);
    printf(" - Ions outwards inside cathode\n");
    NPS += NIntegration(*Sfi_InPlusPtr, 0.01, fusor.a - 0.000001);
    printf(" - Ions inward outside cathode\n");
    NPS += NIntegration(*Sfi_OutMinPtr, fusor.a, fusor.b - 0.000001);
    printf(" - Ions outwards outside cathode\n");
    NPS += NIntegration(*Sfi_OutPlusPtr, fusor.a, fusor.b - 0.000001);

    return NPS;
}




/**
    This function writes a given function to a file.
*/
void print_data_dd(double (*funcPtr)(double), double Start, double End, double step, char name[])
{
    FILE * output;
    output = fopen(name,"w");
    for (double r = Start; r<End; r+=step)
    {
        fprintf (output, "%E,%E\n",r, (*funcPtr)(r));
    }

    printf("%s: writing done\n",name);

    fclose(output);
    return;
}

/**
    This function writes a given function to a file.
*/
void print_data_ddd(double (*funcPtr)(double,double), double Start, double End, double step, double sec, char name[])
{
    FILE * output;
    output = fopen(name,"w");

    for (double r = Start; r<End; r+=step)
    {
        fprintf (output, "%E,%E\n",r, (*funcPtr)(sec, r));
    }

    printf("%s: writing done\n",name);

    fclose(output);
    return;
}

void print_table(int choice, char name[])
{
    FILE * output;
    output = fopen(name,"w");

    switch (choice)
    {
        case 1:
        for ( int i = 0; i < N_TABLE; i++)
        {
            fprintf(output,"%d, %E\n",i,Table.A[i]);
        }
        printf("Printed ATable\n");
        break;
        case 2:
        for ( int i = 0; i < N_TABLE; i++)
        {
            fprintf(output,"%d, %E\n",i,Table.K[N_TABLE-2][i]);
        }
        printf("Printed KTable\n");
        break;
        case 10:
        for ( int i = 0; i < N_TABLE; i++)
        {
            fprintf(output,"%E, %E\n",Table.R[i],Table.S_0[i]);
        }
        printf("Printed STable\n");
        break;
        case 11:
        for ( int i = 0; i < N_TABLE; i++)
        {
            fprintf(output,"%E, %E\n",Table.R[i],Table.S_1[i]);
        }
        printf("Printed S1Table\n");
        break;
        case 12:
        for ( int i = 0; i < N_TABLE; i++)
        {
            fprintf(output,"%E, %E\n",Table.R[i],Table.S_2[i]);
        }
        printf("Printed S2Table\n");
        break;
        case 13:
        for ( int i = 0; i < N_TABLE; i++)
        {
            fprintf(output,"%E, %E\n",Table.R[i],Table.S_3[i]);
        }
        printf("Printed S3Table\n");
        break;
        case 14:
        for ( int i = 0; i < N_TABLE; i++)
        {
            fprintf(output,"%E, %E\n",Table.R[i],Table.S_4[i]);
        }
        printf("Printed S4Table\n");
        break;
        case 15:
        for ( int i = 0; i < N_TABLE; i++)
        {
            fprintf(output,"%E, %E\n",Table.R[i],Table.S_5[i]);
        }
        printf("Printed S5Table\n");
        break;
    }

    fclose(output);

    return;
}
