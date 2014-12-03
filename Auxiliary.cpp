#include <stdio.h>

#include "constants.h"
#include "Auxiliary.h"

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

