#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "constants.h"
#include "Auxiliary.h"

/**
    This function writes a given function to a file.
*/
void plot_function_1D(double (*funcPtr)(double), double Start, double End, double step, const char name[], const char input_file_name[])
{
    FILE *  input;
    FILE * output;
    input  = fopen(input_file_name,"r");
    output = fopen(name,"w");

    if ( input == NULL || output == NULL)
    {
        printf("could not open file or open new file");

        fclose(input);
        fclose(output);
        return;
    }

    double value;
    char c = '0';

    while ((c = fgetc(input)) != EOF)
    {
       // fread(&c, 1, 1, input);

        if ( c == '@')
        {
            for (double r = Start; r<End; r+=step)
            {
                value = (*funcPtr)(r);
                fprintf (output, "%E %E\n",r, value );
                printf("%E \t %E \n",r, value);
            }
        }
        else
        {
            printf("%c",c);
        }
    }

//    printf("e\n");
    printf("# %s: done writing to file\n",name);

    fclose(input);
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

    if ( input == NULL || output == NULL)
    {
        printf("could not open file or open new file");

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
        while ((c = fgetc(input)) != EOF)
        {
            for (double r = Start1; r<End1; r+=step1)
            {
                for ( double dr = Start1; dr<End1; dr+=step2)
                {
                    value = (*funcPtr)(r, dr);
                    fprintf(output, "%E ", value );
                    printf("%E ", value);
                }
                fprintf(output, "\n");
                printf("\n");
            }
        }
    }

    printf("# %s: done writing to file\n",name);

    fclose(input);
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

    if ( input == NULL || output == NULL)
    {
        printf("could not open file or open new file\n");

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
                printf("%E %E \n", Table->R[i], table[i]);
            }
        }
        else
        {
            printf("%c",c);
        }
    }

    printf("# %s: done writing to file\n",name);

    fclose(input);
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

    if ( input == NULL || output == NULL)
    {
        printf("could not open file or open new file\n");

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
                    printf("%E ", table[i][j]);
                    fprintf(output, "%E ", table[i][j]);
                }
                fprintf(output, "\n");
                printf("\n");
            }
        }
        else
        {
            printf("%c",c);
        }
    }

    printf("# %s: done writing to file\n",name);

    fclose(input);
    fclose(output);

}

/**
    This function reads date from the input file.
    TODO: make this function so the order of settings in the input file does
    not matter anymore.
*/
void readfile(FILE * input)
{
    char str [80];

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

    // potential
    fscanf (input, "%s", str); // read word 'Potential'
    fscanf (input, "%s", str); // read number

    if ( !strcmp(str, "true;") )
        printbool->potential = true;
    else
        printbool->potential = false;

    // SIIEE
    fscanf (input, "%s", str); // read word 'SIIEE'
    fscanf (input, "%s", str); // read number

    if ( !strcmp(str, "true;") )
        printbool->SIIEE = true;
    else
        printbool->SIIEE = false;

    // Cross_section
    fscanf (input, "%s", str); // read word 'Cross_section'
    fscanf (input, "%s", str); // read number

    if ( !strcmp(str, "true;") )
        printbool->Cross_section = true;
    else
        printbool->Cross_section = false;

    // Survival
    fscanf (input, "%s", str); // read word 'Survival'
    fscanf (input, "%s", str); // read number

    if ( !strcmp(str, "true;") )
        printbool->Survival = true;
    else
        printbool->Survival = false;

    // Atable
    fscanf (input, "%s", str); // read word 'Atable'
    fscanf (input, "%s", str); // read number

    if ( !strcmp(str, "true;") )
        printbool->Atable = true;
    else
        printbool->Atable = false;


    // KernelTable
    fscanf (input, "%s", str); // read word 'KernelTable'
    fscanf (input, "%s", str); // read number

    if ( !strcmp(str, "true;") )
        printbool->KernelTable = true;
    else
        printbool->KernelTable = false;

    // Stable
    fscanf (input, "%s", str); // read word 'Stable'
    fscanf (input, "%s", str); // read number

    if ( !strcmp(str, "true;") )
        printbool->Stable = true;
    else
        printbool->Stable = false;

    // Spectrum
    fscanf (input, "%s", str); // read word 'Spectrum'
    fscanf (input, "%s", str); // read number

    if ( !strcmp(str, "true;") )
        printbool->Spectrum = true;
    else
        printbool->Spectrum = false;

    // NSR
    fscanf (input, "%s", str); // read word 'NSR'
    fscanf (input, "%s", str); // read number

    if ( !strcmp(str, "true;") )
        printbool->NSR = true;
    else
        printbool->NSR = false;

    // NPR
    fscanf (input, "%s", str); // read word 'NPR'
    fscanf (input, "%s", str); // read number

    if ( !strcmp(str, "true;") )
        printbool->NPR = true;
    else
        printbool->NPR = false;



}

void print_program_parameters(void)
{
	// print fusor parameters
    printf("# Fusor Cathode Radius : \t %f \n",fusor->a);
    printf("# Fusor Anode Radius   : \t %f \n",fusor->b);
    printf("# Voltage on Fusor     : \t %f \n",fusor->V0);
    printf("# Pressure in Fusor    : \t %f \n",pressure);
    printf("# Current though Fusor : \t %f \n",Itot);
    printf("# Transparenty of gird : \t %f \n",fusor->Tc);
	printf("# Ngas                 : \t %f \n",ngas);
	
	// print output parameters    
	printf("\n# - Which data is printen to files and screen:\n");
	
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
	
	return;
}

