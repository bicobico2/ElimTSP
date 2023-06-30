/*
 * This is part of ElimTSP, a set of codes for eliminating edges in the
 * (symmetric) traveling salesman problem. The codes are based on the
 * paper "Local elimination in the traveling salesman problem" by William Cook,
 * Keld Helsgaun, Stefan Hougardy, and Rasmus T. Schroeder.
 */

/*
MIT License

Copyright (c) 2023  William Cook

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/


/****************************************************************************/
/*                                                                          */
/*                 Edge Elimination: Generate k-Swap Code                   */
/*                                                                          */
/*  Generates C code for testing all ways to reconnect a tour after the     */
/*  removal of k edges (for k = 3, 4, or 5).                                */
/*                                                                          */
/*  Written by:  Cook, 2013                                                 */
/*                                                                          */
/****************************************************************************/


#include <stdio.h>

#define SWAP(a,b,t) (((t)=(a)),((a)=(b)),((b)=(t)))
#define KNUM 5

char alph[12] =  {'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L'};
char ralph[12] = {'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l'};

int main (int ac, char **av);
static void perms (int k, int *t, int *count);
static void gentours (int *t, int *count);
static void tours (int k, int *t, int *reverse, int *count);
static void checktour (int *t, int *reverse, int *yesno);

int main (int ac, char **av)
{
    int t[KNUM], i, count = 1;

    for (i = 0; i < KNUM; i++) t[i] = i;
    perms (1, t, &count);

    return 1;
}

static void perms (int k, int *t, int *count)
{
    int i, tmp;

    if (k == KNUM) {
        gentours (t, count);
        return;
    }

    for (i = k; i < KNUM; i++) {
        SWAP (t[i], t[k], tmp);
        perms (k+1, t, count);
        SWAP (t[i], t[k], tmp);
    }
}

static void gentours (int *t, int *count)
{
    int reverse[KNUM], i;

    for (i = 0; i < KNUM; i++) reverse[i] = 0;
    tours (1, t, reverse, count);
}

static void tours (int k, int *t, int *reverse, int *count)
{
    int i, j, p, q, yesno;
    char add[2*KNUM];

    if (k == KNUM) {
        checktour (t, reverse, &yesno);
        if (yesno == 1) {
            for (i = 0; i < KNUM; i++) {
                j = (i+1)%KNUM;
                if (reverse[i]) p = 2*t[i]+1;
                else            p = 2*((t[i]+1)%KNUM);

                if (reverse[j]) q = 2*((t[j]+1)%KNUM);
                else            q = 2*t[j]+1;

                add[2*i]   = ralph[p];
                add[2*i+1] = ralph[q];
            }

            printf ("    /* Case %d: Add", (*count)++);
            for (i = 0; i < KNUM; i++) {
                printf (" %c-%c", add[2*i], add[2*i+1]);
            }
            printf (" for tour ");
            for (i = 0; i < KNUM; i++) {
                if (reverse[i]) printf ("%c", ralph[t[i]]);
                else            printf ("%c", alph[t[i]]); 
            }
            printf (" */\n");

            printf ("    if (M[%c][%c]", add[0], add[1]);
            for (i = 1; i < KNUM; i++) {
                printf (" + M[%c][%c]", add[2*i], add[2*i+1]);
            }
            printf (" < r) {\n");
            for (i = 0; i < KNUM; i++) {
                if (i % 2 == 0) {
                    printf ("       I[%d]=%c; I[%d]=%c;",
                             2*i, add[2*i], 2*i+1, add[2*i+1]);
                } else {
                    printf (" I[%d]=%c; I[%d]=%c;\n",
                             2*i, add[2*i], 2*i+1, add[2*i+1]);
                }
            }
            if (i % 2 == 1) printf ("\n");
            printf ("       return;\n");
            printf ("    }\n\n");
        }
        return;
    }

    tours (k+1, t, reverse, count);
    reverse[k] = 1;
    tours (k+1, t, reverse, count);
    reverse[k] = 0;
}

static void checktour (int *t, int *reverse, int *yesno)
{
    int i;

    *yesno = 1;

    for (i = 0; i < KNUM; i++) {
        if (reverse[i] == 0 && reverse[(i+1)%KNUM] == 0 &&
                                     t[(i+1)%KNUM] == ((t[i]+1) % KNUM)) {
            *yesno = 0; return;
        }

        if (reverse[i] == 1 && reverse[(i+1)%KNUM] == 1 &&
                                     t[(i+1)%KNUM] == ((t[i]-1) % KNUM)) {
            *yesno = 0; return;
        }
    }
}
