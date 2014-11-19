/**
    The calculation of the source rate of the source rate for the
    class II ions.
*/


#include "constants.h"
#include "ClassIISourceFunction.h"

void S(void)
{
    int i;
    for ( i = N_TABLE-1; i > -1; i--)
        S_0(i);

    for ( i = 0; i < N_TABLE; i++)
        Table.S_0[i] += Table.A[i];

    for ( i = N_TABLE-1; i > -1; i--)
        S_1(i);

    for ( i = 0; i < N_TABLE; i++)
        Table.S_1[i] += Table.A[i];

    for ( i = N_TABLE-1; i > -1; i--)
        S_2(i);

    for ( i = 0; i < N_TABLE; i++)
        Table.S_2[i] += Table.A[i];

    for ( i = N_TABLE-1; i > -1; i--)
        S_3(i);

    for ( i = 0; i < N_TABLE; i++)
        Table.S_3[i] += Table.A[i];

    for ( i = N_TABLE-1; i > -1; i--)
        S_4(i);

    for ( i = 0; i < N_TABLE; i++)
        Table.S_4[i] += Table.A[i];

    for ( i = N_TABLE-1; i > -1; i--)
        S_5(i);

    for ( i = 0; i < N_TABLE; i++)
        Table.S_5[i] += Table.A[i];

    return;
}

void S_0(int r)
{
    double step = (fusor.b-fusor.a)/(N_TABLE-1);
    double dot = 0;

    if ( r == N_TABLE-1)
    {
        for ( int i = 0; i < N_TABLE; i++)
            dot += Table.A[i] * Table.K[r][i];

        Table.S_0[r] = step * dot;
    }
    else
    {
        for ( int i = 0; i < N_TABLE; i++)
            dot += Table.A[i] * Table.K[r][i];

        Table.S_0[r] = step * dot + Table.S_0[r+1];
    }

    return;
}

void S_1(int r)
{
    double step = (fusor.b-fusor.a)/(N_TABLE-1);
    double dot = 0;

    if ( r == N_TABLE-1)
    {
        for ( int i = 0; i < N_TABLE; i++)
            dot += Table.S_0[i] * Table.K[r][i];

        Table.S_1[r] = step * dot;
    }
    else
    {
        for ( int i = 0; i < N_TABLE; i++)
            dot += Table.S_0[i] * Table.K[r][i];

        Table.S_1[r] = step * dot + Table.S_1[r+1];
    }

    return;
}

void S_2(int r)
{
    double step = (fusor.b-fusor.a)/(N_TABLE-1);
    double dot = 0;

    if ( r == N_TABLE-1)
    {
        for ( int i = 0; i < N_TABLE; i++)
            dot += Table.S_1[i] * Table.K[r][i];

        Table.S_2[r] = step * dot;
    }
    else
    {
        for ( int i = 0; i < N_TABLE; i++)
            dot += Table.S_1[i] * Table.K[r][i];

        Table.S_2[r] = step * dot + Table.S_2[r+1];
    }

    return;
}

void S_3(int r)
{
    double step = (fusor.b-fusor.a)/(N_TABLE-1);
    double dot = 0;

    if ( r == N_TABLE-1)
    {
        for ( int i = 0; i < N_TABLE; i++)
            dot += Table.S_2[i] * Table.K[r][i];

        Table.S_3[r] = step * dot;
    }
    else
    {
        for ( int i = 0; i < N_TABLE; i++)
            dot += Table.S_2[i] * Table.K[r][i];

        Table.S_3[r] = step * dot + Table.S_3[r+1];
    }

    return;
}

void S_4(int r)
{
    double step = (fusor.b-fusor.a)/(N_TABLE-1);
    double dot = 0;

    if ( r == N_TABLE-1)
    {
        for ( int i = 0; i < N_TABLE; i++)
            dot += Table.S_3[i] * Table.K[r][i];

        Table.S_4[r] = step * dot;
    }
    else
    {
        for ( int i = 0; i < N_TABLE; i++)
            dot += Table.S_3[i] * Table.K[r][i];

        Table.S_4[r] = step * dot + Table.S_4[r+1];
    }

    return;
}

void S_5(int r)
{
    double step = (fusor.b-fusor.a)/(N_TABLE-1);
    double dot = 0;

    if ( r == N_TABLE-1)
    {
        for ( int i = 0; i < N_TABLE; i++)
            dot += Table.S_4[i] * Table.K[r][i];

        Table.S_5[r] = step *  dot;
    }
    else
    {
        for ( int i = 0; i < N_TABLE; i++)
            dot += Table.S_4[i] * Table.K[r][i];

        Table.S_5[r] = step * dot + Table.S_5[r+1];
    }

    return;
}


