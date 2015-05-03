/**
    The calculation of the source rate of the source rate for the
    class II ions.
*/


#include "constants.h"
#include "SourceFunction.h"

void S(void)
{
    // a volterra solver
    int i, j;
    double h, sum;

    // case i = 0
    Table.S[0] = Table.A[0];

    for ( i = 1; i < N_TABLE; i++)
    {
        h = giveAnodeRadius() / N_TABLE;
        sum = 0;

        Table.S[i]  = h / ( 1 - (h/2) * Table.K[i][i] );

        for ( j = 1; j < i - 1; j++)
            sum += Table.K[i][j] * Table.S[j];

        Table.S[i] *= ( 0.5 * Table.K[i][0] * Table.S[0] + sum);
        Table.S[i] += Table.A[i];

    }

    return;
}


