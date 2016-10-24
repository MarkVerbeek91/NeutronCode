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
#include <math.h>
#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>     /* getenv */
#include <unistd.h>
#include <string.h>
#include <ctype.h>

#include "includes.h"

int main( int argc, char **argv )
{
	// notifi user that debug mode is on.
	#ifdef DEBUG_PARAMETER
        fprintf(stdout,"# Debug parameter turned on\n");
    #endif

	// Initalise needed variables.
    initBool();
    InitCrossSectionConstands();

	// open input file, default or supplied by user.
	FILE * input;
    int opt;
	char *cvalue = NULL;

	if( argc >= 2 )
	{
		fprintf(stdout,"# Analising give arguments\n");

		while ((opt = getopt (argc, argv, "f:p:")) != -1)
		{
			switch (opt)
			{
				case 'f': 		// get file
					cvalue = optarg;
					input = fopen(cvalue,"r");
					if ( input == NULL)
					{
						fprintf(stderr,"input file not found, stopping code");
						return 1;
					}
					else
					{
						readfile(input,cvalue);
					}
					fclose(input);

					break;

				case 'p':		// parameter
					arg_parameter_check(optarg);

					break;
				case '?':		// error case
					if (optopt == 'f')
						printf("Wrong or no input file specified \n");
					else if (isprint (optopt))
						printf ("Unknown option `-%c'.\n", optopt);
					else
						printf ("Unknown option character `\\x%x'.\n", optopt);
					return 1;
				default:
					abort ();

			}
		}
	}
	else
	{
		printf("# Using default input file: input.ini\n");
		input = fopen("input.ini","r");
		readfile(input, "input.ini");
		fclose(input);
	}

    fprintf(stdout,"# -------------------------------------------------------------------------- #\n");
    fprintf(stdout,"# ---------------------------- Start of program ---------------------------- #\n");
    fprintf(stdout,"# -------------------------------------------------------------------------- #\n");

	// write parameter that are used to screen
	print_program_parameters();

    // building the Kernel, filling tables A and R.
	print_comment2scr("Filling tables");
    CalculateTables();

    // Calculating: source rate for first generation of Class II ions.
	print_comment2scr("Calculating S");
    S();
    printf("#  - Done\n");

    // calculate the currents in the cathode
    double I1, I2, I3, I4;

    print_comment2scr("Calculating Currents");
    I1 = I_c1();
    printf("# I1 Current           : \t %.2E\n",I1);
    I2 = I_c2();
    printf("# I2 Current           : \t %.2E\n",I2);
    I3 = I_c3();
    printf("# I3 Current           : \t %.2E\n",I3);
    I4 = I_c4();
    printf("# I4 Current           : \t %.2E\n",I4);

    //EdgeIonFlux = (Itot - I2 - I4) / (I1 + I3);
    EdgeIonFlux = Itot / ( I1 + I2 + I3 + I4);

    printf("# EdgeIonFlux: %.2E\n \n", EdgeIonFlux);

    // calculate neutron production rate (NPR)
    if ( printbool->NPR )
    {
        print_comment2scr("Calculating NPR");
        double NPR = Nps();
        printf("# NPR: %.2E \n\n", NPR);
    }

	print_comment2scr("Writing data to screen or files");
	output_data();

    // program is done
    printf("# ---------------------------- Program Finished ---------------------------- #");
    return 0;
}
