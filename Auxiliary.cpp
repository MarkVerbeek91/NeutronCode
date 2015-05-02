#include <stdio.h>
#include <stdlib.h>

#include "constants.h"
#include "Auxiliary.h"

/**
    This function writes a given function to a file.
*/
void print_data_dd(double (*funcPtr)(double), double Start, double End, double step, char name[])
{
    FILE * output;
    output = fopen(name,"w");
    double value;

    // print data to gnu plot
    printf("unset key\n");
    printf("set xlabel \"position\"\n");
    printf("set ylabel \"a.u\"\n");
    printf("set Titel \"%s\"\n",name);
    printf("set xrange [%f:%f]\n",Start,End);
    printf("plot '-' w l\n");

    for (double r = Start; r<End; r+=step)
    {
        value = (*funcPtr)(r);
        fprintf (output, "%E,%E\n",r, value );
        printf("%f \t %f \n",r, value);
    }

    printf("e\n");
    printf("%s: done writing to file\n",name);

    fclose(output);

    return;
}

/**
    This function writes a given function to a file.
*/
void print_data_ddd(double (*funcPtr)(double,double), double var, double Start, double End, double step, char name[])
{
    FILE * output;
    output = fopen(name,"w");
    double value;

    // print data to gnu plot
    printf("unset key\n");
    printf("set xlabel \"position\"\n");
    printf("set ylabel \"a.u\"\n");
    printf("set Titel \"%s\"\n",name);
    printf("set xrange [%f:%f]\n",Start,End);
    printf("plot '-' w l\n");

    for (double r = Start; r<End; r+=step)
    {
        value = (*funcPtr)(var, r);
        fprintf (output, "%E,%E\n",r, value);
        printf("%f \t %f \n",r, value);
    }

    printf("e\n");
    printf("%s: done writing to file\n",name);

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
        case 9:
        for ( int i = 0; i < N_TABLE; i++)
        {
            fprintf(output,"%E, %E\n",Table.R[i],Table.S[i]);
        }
        printf("Printed STable\n");
        break;
        default:
            printf("nothing to do\n");
        break;
    }

    fclose(output);

    return;
}

void readfile(FILE * input)
{
    char str [80];

    fscanf (input, "%s", str); // read word '#fusor-dimentions'

    // cathode
    fscanf (input, "%s", str); // read word 'cathode'
    printf("%s ",str);

    fscanf (input, "%s", str); // read word '0.05'
    printf("%s \n",str);

    fusor.a = strtof(str,NULL);

    // anode
    fscanf (input, "%s", str); // read word 'anode'
    printf("%s ",str);

    fscanf (input, "%s", str); // read word '0.25'

    fusor.b = strtof(str,NULL);
    printf("\t%f \n",fusor.b);

    // read comment in input file
    fscanf (input, "%s", str); // read word '#fusor-dimentions'

    // Voltage
    fscanf (input, "%s", str); // read word 'voltage'
    printf("%s ",str);

    fscanf (input, "%s", str); // read number

    fusor.V0 = strtof(str,NULL);
    printf("\t%f \n",fusor.V0);

    // Pressure
    fscanf (input, "%s", str); // read word 'pressure'
    printf("%s ",str);

    fscanf (input, "%s", str); // read number

    pressure = strtof(str,NULL);
    printf("\t%f \n",pressure);

    // Current
    fscanf (input, "%s", str); // read word 'Current'
    printf("%s ",str);

    fscanf (input, "%s", str); // read number

    Itot = strtof(str,NULL);
    printf("\t%f \n",Itot);

}

