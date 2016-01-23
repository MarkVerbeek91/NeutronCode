#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "constants.h"
#include "Auxiliary.h"
#include "..\includes.h"

/**
 * This function writes data to files or screen depending on input file.
 */
void output_data(void)
{
	double (*funcPtr)(double);
	double (*funcPtr2)(double, double);

	// writing the potential to a file for plotting
    if ( printbool->potential )
    {
        printf("# Wrinting potential to stdout in gnu format\n");

		funcPtr = &Potential_Phi;
		GNUplot_function_1D(*funcPtr, 0.0, 0.25, 0.001, "GNUplot\\GNU_potential.txt");

        funcPtr2 = &ParticleEnergy2;
        GNUplot_function_2D(*funcPtr2, 0.0, giveAnodeRadius()+0.001, 0.0, giveAnodeRadius()+0.001, 0.01, 0.01, "GNUplot\\GNU_potential.txt");
    }
	if ( printbool2->potential )
	{
		printf("# Writing potential to GNUfile \n");

		funcPtr = &Potential_Phi;
		plot_function_1D(*funcPtr, 0.0, 0.25, 0.001, "output_files\\Potential.csv", "GNUplot\\GNU_potential.txt");

        funcPtr2 = &ParticleEnergy2;
        plot_function_2D(*funcPtr2, 0.0, giveAnodeRadius()+0.001, 0.0, giveAnodeRadius()+0.001, 0.01, 0.01, "output_files\\Particle2.csv", "GNUplot\\GNU_potential.txt");
	}

    // writing the SIIEE to a file for plotting
    if ( printbool->SIIEE )
    {
        printf("# SIIEE calculation\n");

        funcPtr = &SIIEE;
        GNUplot_function_1D(*funcPtr, 1.0, -giveVoltage(), 1, "GNUplot\\GNU_SIIEE.txt");
    }
    if ( printbool2->SIIEE )
    {
        printf("# SIIEE calculation\n");

        funcPtr = &SIIEE;
        plot_function_1D(*funcPtr, 1.0, -giveVoltage(), 1, "output_files\\SIIEE.csv", "GNUplot\\GNU_SIIEE.txt");
    }

    // writing the Cross sections for Charge Exchange, Iononisation and the sum
    // of those to files to a file for plotting
    if ( printbool->Cross_section )
    {
        printf("# Cross section calculation\n");

		funcPtr = &CrosssecCX;
        GNUplot_function_1D(*funcPtr, 1.0, 500000, 10, "GNUplot\\GNU_Cross_sections.txt");

        funcPtr = &CrosssecIon;
        GNUplot_function_1D(*funcPtr, 1.0, 500000, 10, "GNUplot\\GNU_Cross_sections.txt");

        funcPtr = &CrosssecTot;
        GNUplot_function_1D(*funcPtr, 1.0, 500000, 10, "GNUplot\\GNU_Cross_sections.txt");

    }
	if ( printbool2->Cross_section )
    {
        printf("# Cross section calculation\n");

		funcPtr = &CrosssecCX;
        plot_function_1D(*funcPtr, 1.0, 500000, 10, "output_files\\CrosssecCX.txt", "GNUplot\\GNU_Cross_sections.txt");

        funcPtr = &CrosssecIon;
        plot_function_1D(*funcPtr, 1.0, 500000, 10, "output_files\\CrosssecIon.txt", "GNUplot\\GNU_Cross_sections.txt");

        funcPtr = &CrosssecTot;
        plot_function_1D(*funcPtr, 1.0, 500000, 10, "output_files\\CrosssecTot.txt", "GNUplot\\GNU_Cross_sections.txt");
    }

    // Writing the survival functions to a file for plotting.
    if ( printbool->Survival )
    {
        printf("# Survival function calculation\n");

        funcPtr = &f;
        GNUplot_function_1D(*funcPtr, giveCathodeRadius(), giveAnodeRadius()+0.001, 0.001, "GNUplot\\GNU_Survival_functions.txt");

        double (*funcPtr2)(double,double);
        funcPtr2 = &g;
        GNUplot_function_2D(*funcPtr2, 0.0, 0.0, giveCathodeRadius(), giveAnodeRadius(), 0.001, 0.001, "GNUplot\\GNU_Survival_functions.txt");
    }
	if ( printbool2->Survival )
    {
        printf("# Survival function calculation\n");

        funcPtr = &f;
        plot_function_1D(*funcPtr, giveCathodeRadius(), giveAnodeRadius()+0.001, 0.001, "output_files\\f.txt", "GNUplot\\GNU_Survival_functions.txt");

        double (*funcPtr2)(double,double);
        funcPtr2 = &g;
        plot_function_2D(*funcPtr2, 0.0, 0.0, giveCathodeRadius(), giveAnodeRadius(), 0.001, 0.001, "output_files\\g.txt", "GNUplot\\GNU_Survival_functions.txt");
    }

	// writing K to file
    if ( printbool->KernelTable )
    {
        printf("# Writing Kernel to screen \n");
        GNUplot_table_2D(Table->K, "GNUplot\\GNU_Ktable.txt");
    }
	if ( printbool2->KernelTable )
    {
        printf("# Writing Kernel to file \n");
		plot_table_2D(Table->K, "output_files\\KTable.csv", "GNUplot\\GNU_Ktable.txt");
    }

	// writing A to file
    if ( printbool->Atable )
    {
        printf("# Writing A to screen \n");
		GNUplot_table_1D(Table->A, "GNUplot\\GNU_Atable.txt");
    }
    if ( printbool2->Atable )
    {
        printf("# Writing A to file \n");
        plot_table_1D(Table->A, "output_files\\ATable.csv", "GNUplot\\GNU_Atable.txt");
    }

    // writing S to file
    if ( printbool->Stable )
    {
        printf("# Writing S to screen \n");
        GNUplot_table_1D(Table->S, "GNUplot\\GNU_Stable.txt");
    }
    if ( printbool2->Stable )
    {
        printf("# Writing S to file \n");
        plot_table_1D(Table->S, "output_files\\STable.csv", "GNUplot\\GNU_Stable.txt");
    }

	    // write spectrum to files
    // TODO: write to folder, are a lot of files
    //       or find a way of making one file
/*    if ( printbool->Spectrum )
    {
        printf("# Printing Energy spectrum to files\n");

        printf("# Energy spectrum of ions going inwards\n");

        double (*IonSpectrumInwardsPtr)(double, double);
        IonSpectrumInwardsPtr = &IonSpectrumInwards;

        double (*IonSpectrumOutwardsPtr)(double, double);
        IonSpectrumOutwardsPtr = &IonSpectrumOutwards;

//        double r = 0.06;

//        print_data_ddd(*IonSpectrumInwardsPtr , 10, -giveVoltage(),10,r,"IonSpectrumInwards.GNUfile");
//        print_data_ddd(*IonSpectrumOutwardsPtr, 10, -giveVoltage(),10,r,"IonSpectrumOutwards.GNUfile");

//        int j = 6;

        int j = 1;
        for ( double r = 0.01;  r < 0.25; r=r+0.01)
        {
            char filename1[100];
            sprintf( filename1, "IonSpectrumInwards%d.txt", j);

            char filename2[100];
            sprintf( filename2, "IonSpectrumOutwards%d.GNUfile", j);

            //print_data_ddd(*IonSpectrumInwardsPtr , 10, -giveVoltage(),10,r,filename1, 5);
            //print_data_ddd(*IonSpectrumOutwardsPtr, 10, -giveVoltage(),10,r,filename2, 5);

            j++;
        }

    } */

    // Calculate the neutron source rate
    if ( printbool->NSR )
    {
        printf("# Printing neutron source rate to file:\n");

        // Ions
        funcPtr = &FusionRateIons;
        printf("# Neutrons from inwards ions\n");
        GNUplot_function_1D(*funcPtr, 1e-3, giveAnodeRadius(), 1e-3, "GNUplot\\GNU_NSR.txt");

        // Neutrals Class I
        funcPtr = &FusionRateNeutralsClassI;
        printf("# In cathode, inwards\n");
        GNUplot_function_1D(*funcPtr, 1e-3, giveAnodeRadius(), 1e-3, "GNUplot\\GNU_NSR.txt");

        // Neutrals ClassII
        funcPtr = &FusionRateNeutralsClassII;
        printf("# In cathode, inwards\n");
        GNUplot_function_1D(*funcPtr, 1e-3, giveAnodeRadius(), 1e-3, "GNUplot\\GNU_NSR.txt");
    }

	return;
}

/**
    This function writes a given function to a file.
*/
void plot_function_1D(double (*funcPtr)(double), double Start, double End, double step, const char name[], const char input_file_name[])
{
    FILE *  input;
    FILE * output;
    input  = fopen(input_file_name,"r");
    output = fopen(name,"w");

    if ( input == NULL )
    {
        printf("# Could not open file: %s\n", input_file_name);

        fclose(input);
        fclose(output);
        return;
    }
	if ( output == NULL )
	{
        printf("# Could not open file: %s\n", name);

        fclose(input);
        fclose(output);
        return;
    }

    double value;
    char c = '0';

    while ((c = fgetc(input)) != EOF)
    {
        if ( c == '@')
        {
            for (double r = Start; r<End; r+=step)
            {
                value = (*funcPtr)(r);
                fprintf (output, "%E %E\n",r, value );
            }
        }
        else
        {
            fprintf(output, "%c", c);
        }
    }

    printf("# %s: done writing to file\n",name);

    fclose(input);
    fclose(output);

    return;
}

void GNUplot_function_1D(double (*funcPtr)(double), double Start, double End, double step, const char name[])
{
    FILE * output;
    output = fopen(name,"r");

	if ( output == NULL )
	{
        printf("# Could not open file: %s\n", name);
        fclose(output);
        return;
    }

    double value;
    char c = '0';

    while ((c = fgetc(output)) != EOF)
    {
        if ( c == '@')
        {
            for (double r = Start; r<End; r+=step)
            {
                value = (*funcPtr)(r);
                printf("%E \t %E \n",r, value);
            }
        }
        else
            printf("%c",c);
    }

    printf("# \t done writing to screen\n");

    fclose(output);

    return;
}


/**
    This function writes a given function to a file.
*/
void plot_function_2D(double (*funcPtr)(double, double), double Start1, double End1, double Start2, double End2, double step1, double step2, const char name[], const char input_file_name[])
{
    FILE *  input;
    FILE * output;
    input  = fopen(input_file_name,"r");
    output = fopen(name,"w");

    if ( input == NULL )
    {
        printf("# Could not open file: %s\n", input_file_name);

        fclose(input);
        fclose(output);
        return;
    }
	if ( output == NULL )
	{
        printf("# Could not open file: %s\n", name);

        fclose(input);
        fclose(output);
        return;
    }

    double value;
    char c = '0';

    if ( Start1 == End1 )
    {
        while ((c = fgetc(input)) != EOF)
        {
            if ( c == '@')
            {
                for (double r = Start2; r<End2; r+=step2)
                {
                    value = (*funcPtr)(Start1, r);
                    fprintf (output, "%E %E\n",r, value );
                }
            }
            else
            {
                fprintf(output, "%c", c);
            }
        }
    }
    else
    {
        while ((c = fgetc(input)) != EOF)
        {
            for (double r = Start1; r<End1; r+=step1)
            {
                for ( double dr = Start1; dr<End1; dr+=step2)
                {
                    value = (*funcPtr)(r, dr);
                    fprintf(output, "%E ", value );
                }
                fprintf(output, "\n");
            }
        }
    }

    printf("# %s: done writing to file\n",name);

    fclose(input);
    fclose(output);

    return;
}

void GNUplot_function_2D(double (*funcPtr)(double, double), double Start1, double End1, double Start2, double End2, double step1, double step2, const char name[])
{
    FILE * output;
    output = fopen(name,"r");

	if ( output == NULL )
	{
        printf("# Could not open file: %s\n", name);
        fclose(output);
        return;
    }

    double value;
    char c = '0';

    if ( Start1 == End1 )
    {
        while ((c = fgetc(output)) != EOF)
        {
            if ( c == '@')
            {
                for (double r = Start2; r<End2; r+=step2)
                {
                    value = (*funcPtr)(Start1, r);
                    printf("%E \t %E \n",r, value);
                }
            }
            else
            {
                printf("%c",c);
            }
        }
    }
    else
    {
        while ((c = fgetc(output)) != EOF)
        {
            for (double r = Start1; r<End1; r+=step1)
            {
                for ( double dr = Start1; dr<End1; dr+=step2)
                {
                    value = (*funcPtr)(r, dr);
                    printf("%E ", value);
                }
                printf("\n");
            }
        }
    }

    printf("# \t done writing to screen \n");

    fclose(output);

    return;
}

/**
    print data of a table to file
*/
void plot_table_1D(double *table, const char name[], const char input_file_name[])
{
    FILE *  input;
    FILE * output;
    input  = fopen(input_file_name,"r");
    output = fopen(name,"w");

    if ( input == NULL )
    {
        printf("# Could not open file: %s\n", input_file_name);

        fclose(input);
        fclose(output);
        return;
    }
	if ( output == NULL )
	{
        printf("# Could not open file: %s\n", name);

        fclose(input);
        fclose(output);
        return;
    }

    char c = '0';

    while ((c = fgetc(input)) != EOF)
    {
        if ( c == '@')
        {
            for ( int i = 0; i < N_TABLE; i++)
            {
                fprintf(output, "%E %E \n", Table->R[i], table[i]);
            }
        }
        else
        {
            fprintf(output, "%c", c);
        }
    }

    printf("# %s: done writing to file\n",name);

    fclose(input);
    fclose(output);

    return;
}

void GNUplot_table_1D(double *table, const char name[])
{
    FILE * output;
    output = fopen(name,"r");

	if ( output == NULL )
	{
        printf("# Could not open file: %s\n", name);

        fclose(output);
        return;
    }

    char c = '0';

    while ((c = fgetc(output)) != EOF)
    {
        if ( c == '@')
        {
            for ( int i = 0; i < N_TABLE; i++)
            {
                fprintf(output, "%E %E \n", Table->R[i], table[i]);
            }
        }
        else
        {
            printf("%c",c);
        }
    }

    printf("# %s: done writing to file\n",name);

    fclose(output);

    return;
}


/**
    Print data of a 2D table to screen.
*/
void plot_table_2D(double (*table)[N_TABLE], const char name[], const char input_file_name[])
{
    FILE *  input;
    FILE * output;
    input  = fopen(input_file_name,"r");
    output = fopen(name,"w");

    if ( input == NULL )
    {
        printf("# Could not open file: %s\n", name);

        fclose(input);
        fclose(output);
        return;
    }
	if ( output == NULL )
	{
        printf("# Could not open file: %s\n", input_file_name);

        fclose(input);
        fclose(output);
        return;
    }

    char c = '0';

    while ((c = fgetc(input)) != EOF)
    {
        if ( c == '@')
        {
            for ( int j = 0; j < N_TABLE; j++)
            {
                for ( int i = 0; i < N_TABLE; i++)
                {
                    fprintf(output, "%E ", table[i][j]);
                }
                fprintf(output, "\n");
            }
        }
        else
        {
            fprintf(output, "%c", c);
        }
    }

    printf("# %s: done writing to file\n",name);

    fclose(input);
    fclose(output);

	return;
}

void GNUplot_table_2D(double (*table)[N_TABLE], const char name[])
{
    FILE * output;
    output = fopen(name,"r");

	if ( output == NULL )
	{
        printf("# Could not open file: %s\n", name);
        fclose(output);
        return;
    }

    char c = '0';

    while ((c = fgetc(output)) != EOF)
    {
        if ( c == '@')
        {
            for ( int j = 0; j < N_TABLE; j++)
            {
                for ( int i = 0; i < N_TABLE; i++)
                {
                    printf("%E ", table[i][j]);
                }
                printf("\n");
            }
        }
        else
        {
            printf("%c", c);
        }
    }

    printf("# %s: done writing to file\n",name);

    fclose(output);

	return;
}

void arg_parameter_check(char *cvalue)
{
	char parameter_name[6];
	double parameter_value;

	strncpy ( parameter_name, cvalue, 5);
					parameter_name[6] = '\0';

	while (*cvalue && !(isdigit(*cvalue) || ((*cvalue == '-' || *cvalue == '+') && isdigit(*(cvalue + 1)))))
		cvalue++;

	parameter_value = atof (cvalue);

	if (strcmp(parameter_name, "cath=") == 0)
	{
		fusor->a = parameter_value;
	}
	else if (strcmp(parameter_name, "anod=") == 0)
	{
		fusor->b = parameter_value;
	}
	else if (strcmp(parameter_name, "volt=") == 0)
	{
		fusor->V0 = parameter_value;
	}
	else if (strcmp(parameter_name, "pres=") == 0)
	{
		pressure = parameter_value;
	}
	else if (strcmp(parameter_name, "curr=") == 0)
	{
		Itot = parameter_value;
	}
	else if (strcmp(parameter_name, "tran=") == 0)
	{
		fusor->Tc = parameter_value;
	}
	else
	{
		printf ("parameter unknown: %s", parameter_name);
	}

	return;
}

/**
    This function reads date from the input file.
    TODO: make this function so the order of settings in the input file does
    not matter anymore.
*/
void readfile(FILE * input, char* name)
{
    char str [80];

	fprintf(stdout,"# Reading input file: %s \n\n", name);

    fscanf (input, "%s", str); // read word '#fusor-dimentions'

    // cathode
    fscanf (input, "%s", str); // read word 'cathode'
    fscanf (input, "%s", str); // read word '0.05'
    fusor->a = strtof(str,NULL);

    // anode
    fscanf (input, "%s", str); // read word 'anode'
    fscanf (input, "%s", str); // read word '0.25'
    fusor->b = strtof(str,NULL);

    // read comment in input file
    fscanf (input, "%s", str); // read word '#fusor-parameters'

    // Voltage
    fscanf (input, "%s", str); // read word 'voltage'
	fscanf (input, "%s", str); // read number
    fusor->V0 = strtof(str,NULL);

    // Pressure
    fscanf (input, "%s", str); // read word 'pressure'
    fscanf (input, "%s", str); // read number
    pressure = strtof(str,NULL);

    // Current
    fscanf (input, "%s", str); // read word 'Current'
    fscanf (input, "%s", str); // read number
    Itot = strtof(str,NULL);

    // Transparenty
    fscanf (input, "%s", str); // read word 'Transparenty'

    fscanf (input, "%s", str); // read number
    fusor->Tc = strtof(str,NULL);

	q = 1;
	pressure = 0.5;
	Itot = 0.1;
	Tgas = 400;
	ngas = 6.022e23 * pressure / (8.314 * Tgas); //

    /* READ IN BOOLS FROM FILE */
    fscanf (input, "%s", str); // read word '#outfiles'

	printbool->potential = readline_ini_file_bools(input);
	printbool->SIIEE = readline_ini_file_bools(input);
	printbool->Cross_section = readline_ini_file_bools(input);
	printbool->Survival = readline_ini_file_bools(input);
	printbool->Atable = readline_ini_file_bools(input);
	printbool->KernelTable = readline_ini_file_bools(input);
	printbool->Stable = readline_ini_file_bools(input);
	printbool->Spectrum = readline_ini_file_bools(input);
	printbool->NSR = readline_ini_file_bools(input);
	printbool->NPR = readline_ini_file_bools(input);

	fscanf (input, "%s", str); // read word '#print_to_screen'

	printbool2->potential = readline_ini_file_bools(input);
	printbool2->SIIEE = readline_ini_file_bools(input);
	printbool2->Cross_section = readline_ini_file_bools(input);
	printbool2->Survival = readline_ini_file_bools(input);
	printbool2->Atable = readline_ini_file_bools(input);
	printbool2->KernelTable = readline_ini_file_bools(input);
	printbool2->Stable = readline_ini_file_bools(input);
	printbool2->Spectrum = readline_ini_file_bools(input);
	printbool2->NSR = readline_ini_file_bools(input);
	printbool2->NPR = readline_ini_file_bools(input);

}

bool readline_ini_file_bools(FILE* input)
{
    char str[80];
	fscanf (input, "%s", str); // read word 'NPR'
    fscanf (input, "%s", str); // read number

    if ( !strcmp(str, "true;") )
        return true;
    else
        return false;
}

void print_program_parameters(void)
{
	print_comment2scr("Using following parameters");
//    fprintf(stdout,"# - Using following parameters ------------------------- #\n");

	// print fusor parameters
    printf("# Fusor Cathode Radius : \t\t %f \n",fusor->a);
    printf("# Fusor Anode Radius   : \t\t %f \n",fusor->b);
    printf("# Voltage on Fusor     : \t %+11.1f \n",fusor->V0);
    printf("# Pressure in Fusor    : \t\t %f \n",pressure);
    printf("# Current though Fusor : \t\t %f \n",Itot);
    printf("# Transparenty of gird : \t\t %f \n",fusor->Tc);
	printf("# Ngas                 : \t\t %E \n",ngas);

	// print output parameters
	print_comment2scr("Which data is printen to screen");
//	fprintf(stdout,"\n# - Which data is printen to screen -------------------- #\n");

	printf("# Potential            : \t %s \n", printbool->potential ? "true" : "false");
	printf("# SIIEE                : \t %s \n", printbool->SIIEE ? "true" : "false");
	printf("# Cross_sections       : \t %s \n", printbool->Cross_section ? "true" : "false");
	printf("# Survival functions   : \t %s \n", printbool->Survival ? "true" : "false");
	printf("# Source Table         : \t %s \n", printbool->Atable ? "true" : "false");
    printf("# Kernel               : \t %s \n", printbool->KernelTable ? "true" : "false");
	printf("# Neutral Source Table : \t %s \n", printbool->Stable ? "true" : "false");
	printf("# Spectra              : \t %s \n", printbool->Spectrum ? "true" : "false");
    printf("# NSR                  : \t %s \n", printbool->NSR ? "true" : "false");
    printf("# NPR                  : \t %s \n", printbool->NPR ? "true" : "false");

	// print output parameters
	print_comment2scr("Which data is printen to files");

	printf("# Potential            : \t %s \n", printbool2->potential ? "true" : "false");
	printf("# SIIEE                : \t %s \n", printbool2->SIIEE ? "true" : "false");
	printf("# Cross_sections       : \t %s \n", printbool2->Cross_section ? "true" : "false");
	printf("# Survival functions   : \t %s \n", printbool2->Survival ? "true" : "false");
	printf("# Source Table         : \t %s \n", printbool2->Atable ? "true" : "false");
    printf("# Kernel               : \t %s \n", printbool2->KernelTable ? "true" : "false");
	printf("# Neutral Source Table : \t %s \n", printbool2->Stable ? "true" : "false");
	printf("# Spectra              : \t %s \n", printbool2->Spectrum ? "true" : "false");
    printf("# NSR                  : \t %s \n", printbool2->NSR ? "true" : "false");
    printf("# NPR                  : \t %s \n", printbool2->NPR ? "true" : "false");

	return;
}

void print_comment2scr(char* str)
{

	int len = strlen(str);
	fprintf(stdout,"# - %s ",str);

	for (int i = len; i < 71; i++ )
		fprintf(stdout,"-");

	fprintf(stdout," #\n");

	return;
}


