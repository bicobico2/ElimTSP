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
/*  Written by:  Cook, 2013                                                 */
/*  Last Update: May 27, 2016; July 28, 2016                                */
/*                                                                          */
/*  Exported Functions                                                      */
/*                                                                          */
/*  int CCelim_read_nonpairs (char *fname, int ncount, int **acounts,       */
/*      int ***apairs)                                                      */
/*    READS a file of non-pairs                                             */
/*     -acounts returns the number of pairs for each node (node is the      */
/*      middle of the pair)                                                 */
/*     -apairs returns a list of the pairs for each node                    */
/*                                                                          */
/*  int CCelim_write_nonpairs (char *fname, int ncount, int *pcounts,       */
/*       int **pairs)                                                       */
/*    WRITE a non-pairs file                                                */
/*                                                                          */
/*  int CCelim_nonpairs_to_pairs (CCelim_graph *G, int *nonpcounts,         */
/*      int **nonpairs, int **pcounts, int ***pairs)                        */
/*    CONVERT list of non-pairs to list of pairs (for graph G)              */
/*                                                                          */
/*  int CCelim_pairlist_to_nonpairs (CCelim_graph *G, int plistcount,       */
/*      int *plist, int *nonpcounts, int **nonpairs, int **outpcounts,      */
/*      int ***outpairs)                                                    */
/*    MAKE a list of non-pairs of G contained in union of nonpairs and      */
/*     and plist (where plist is an array of plistcount triples)            */
/*                                                                          */
/*  int CCelim_build_pairprocesslist (int ncount, int ecount, int *elist,   */
/*      int *nonpcounts, int **nonpairs, int *processcount,                 */
/*      int **processlist)                                                  */
/*    CREATE an array of triples of pairs of edges in elist that are not    */
/*      excluded by the nonpairs list (nonpairs can be NULL)                */
/*                                                                          */
/****************************************************************************/

#include <stdio.h>
#include "ccutil.h"
#include "elim.h"

static int pair_in_list (int count, int *list, int a, int b);

/* A nonpair is a pair of edges (ab, ac) that cannot be together in an    */
/* optimal tour. Node a is the center of the pair; we store pairs in a    */
/* node array; the array node node a will contain (b, c) in positions 2*i */
/* and 2*i+1 for some i.                                                  */

/* The file format for a nonpair list is one line for each node, giving  */
/* the index of the node and the # of pairs where it is the center.  If  */
/* the # is positive, then the pairs are given in the following lines,   */
/* one pair per line. The first line in the file gives the number of     */
/* nodes in the problem and the total number of nonpairs in the list.    */

int CCelim_read_nonpairs (char *fname, int ncount, int **acounts, int ***apairs)
{
    int rval = 0, i, a, b, c, tmp;
    int *counts = (int *) NULL, **pairs = (int **) NULL;
    FILE *in = (FILE *) NULL;

    CC_MALLOC (counts, ncount, int);
    CC_MALLOC (pairs, ncount, int *);
    for (a = 0; a < ncount; a++) pairs[a] = (int *) NULL;

    in = fopen (fname, "r");
    if (!in) {
        fprintf (stderr, "could not open %s for reading\n", fname);
        rval = 1; goto CLEANUP;
    }

    fscanf (in, "%d %d", &i, &tmp);
    if (i != ncount) {
        fprintf (stderr, "nonpair file does not match problem\n");
        rval = 1; goto CLEANUP;
    } 

    for (a = 0; a < ncount; a++) {
        fscanf (in, "%d %d\n", &tmp, &(counts[a]));
        if (tmp != a) {
            fprintf (stderr, "missing node %d in nonpairs file\n", a);
            rval = 1; goto CLEANUP;
        }
        if (counts[a] > 0) {
            CC_MALLOC (pairs[a], 2*counts[a], int);
            for (i = 0; i < counts[a]; i++) {
                fscanf (in, "%d %d\n", &b, &c);
                if (b > c) CC_SWAP (b, c, tmp);
                pairs[a][2*i] = b;  pairs[a][2*i+1] = c;
            }
        }
    }

    *acounts = counts;
    *apairs = pairs;

CLEANUP:
    if (in) fclose (in);
    return rval;
}

int CCelim_write_nonpairs (char *fname, int ncount, int *pcounts, int **pairs)
{
    int rval = 0, i, a, b, c, tmp, pcnt = 0;
    FILE *out = (FILE *) NULL;

    out = fopen (fname, "w");
    if (!out) {
        fprintf (stderr, "could not open %s for writing\n", fname);
        rval = 1; goto CLEANUP;
    }

    for (a = 0; a < ncount; a++) pcnt += pcounts[a];
    fprintf (out, "%d %d\n", ncount, pcnt);

    for (a = 0; a < ncount; a++) {
        fprintf (out, " %d %d\n", a, pcounts[a]);
        for (i = 0; i < pcounts[a]; i++) {
            b = pairs[a][2*i]; c = pairs[a][2*i+1];
            if (b > c) CC_SWAP (b, c, tmp);
            fprintf (out, "%d %d\n", b, c);
        }
    }

CLEANUP:
    if (out) fclose (out);
    return rval;
}

/* During elimination we need the pairs of edges in the graph that may   */
/* possibly be in an optimal tour.  If we have a nonpair list, then we   */
/* convert this to a list of valid pairs: run through all pairs in G     */
/* and keep only those that are not in the nonpair list.                 */

int CCelim_nonpairs_to_pairs (CCelim_graph *G, int *nonpcounts, int **nonpairs,
        int **pcounts, int ***pairs)
{
    int rval = 0, i, j, a, b, c, tmp, ncount = G->ncount, mdeg = 0;
    int *acounts = (int *) NULL, **apairs = (int **) NULL, *atmp = (int *) NULL;
    CCelim_node *pa;

    /* acounts+apairs will be the output lists; the atmp array is long   */
    /* enough to handle all pairs for any node; we will copy its content */
    /* into apairs[a] for each node a                                    */

    CC_MALLOC (acounts, ncount, int);
    CC_MALLOC (apairs, ncount, int *);
    for (a = 0; a < ncount; a++) {
        apairs[a] = (int *) NULL;
        if (G->nodelist[a].deg > mdeg) mdeg = G->nodelist[a].deg;
    }
    CC_MALLOC (atmp, 2*mdeg*mdeg, int);
   
    for (a = 0; a < ncount; a++) {
        pa = &G->nodelist[a];
        /* run through all pairs of neighbors (b,c) of a */
        acounts[a] = 0;
        for (i = 0; i < pa->deg; i++) {
            b = pa->neighbors[i];
            for (j = i+1; j < pa->deg; j++) {
                c = pa->neighbors[j];
                /* if (b,c) not in nonpair list, add it to pair list */
                if (!pair_in_list (nonpcounts[a], nonpairs[a], b, c)) {
                    atmp[2*acounts[a]]   = b;
                    atmp[2*acounts[a]+1] = c;
                    acounts[a]++;
                }
            }
        }
        if (acounts[a]) {
            CC_MALLOC (apairs[a], 2*acounts[a], int);
            for (i = 0; i < acounts[a]; i++) {
                b = atmp[2*i];  c = atmp[2*i+1];
                if (b > c) CC_SWAP (b, c, tmp);
                apairs[a][2*i]   = b;
                apairs[a][2*i+1] = c;
            }
        }
    }

    *pcounts = acounts;
    *pairs = apairs;

CLEANUP:
    CC_IFFREE (atmp, int);
    return rval;
}

/* When we are trying to eliminate further pairs, we create a list of  */
/* triples for the remaining pairs, that is, (a, b, c) where the pair  */
/* is (ab, ac). If a nonpair list is available, then we only keep the  */
/* triples that are non already in the nonpair list.                   */

int CCelim_build_pairprocesslist (int ncount, int ecount, int *elist,
        int *nonpcounts, int **nonpairs, int *processcount,
        int **processlist)
{
    int rval = 0, i, j, *pcounts = (int *) NULL, **pairs = (int **) NULL;
    int pcnt = 0, *plist = (int *) NULL, a, b, c, tmp;
    CCelim_graph G;
    CCelim_node *pa;

    /* build a graph from the edge list to be able to run through neighbors */

    rval = CCelim_build_graph (&G, ncount, ecount, elist, (int *) NULL,
                 (CCelim_distobj *) NULL);
    CCcheck_rval (rval, "CCelim_build_graph failed");

    if (nonpairs) {
        /* convert the nonpairs to pairs and build the triples */
        rval = CCelim_nonpairs_to_pairs (&G, nonpcounts, nonpairs, &pcounts,
                                         &pairs);
        CCcheck_rval (rval, "CCelim_nonpairs_to_pairs failed");

        for (a = 0; a < ncount; a++) pcnt += pcounts[a];
        CC_MALLOC (plist, 3*pcnt, int);
        pcnt = 0;
        for (a = 0; a < ncount; a++) {
            for (i = 0; i < pcounts[a]; i++) {
                plist[3*pcnt] = a;
                b = pairs[a][2*i]; c = pairs[a][2*i+1];
                if (b > c) CC_SWAP (b, c, tmp);
                plist[3*pcnt+1] = b; plist[3*pcnt+2] = c;
                pcnt++;
            }
        }
    } else {
        /* no nonpairs, so create a triple for every pair of neighbors */
        for (a = 0; a < ncount; a++) {
            pcnt += ((G.nodelist[a].deg * (G.nodelist[a].deg - 1)) / 2);
        }
        CC_MALLOC (plist, 3*pcnt, int);
        pcnt = 0;
        for (a = 0; a < ncount; a++) {
            pa = &G.nodelist[a];
            for (i = 0; i < pa->deg; i++) {
                b = pa->neighbors[i];
                for (j = i+1; j < pa->deg; j++) {
                    c = pa->neighbors[j];
                    plist[3*pcnt] = a;
                    /* Note: do not swap b and c, since b is in outer loop */
                    if (b < c) { plist[3*pcnt+1] = b; plist[3*pcnt+2] = c; }
                    else       { plist[3*pcnt+1] = c; plist[3*pcnt+2] = b; }
                    pcnt++;
                }
            }
        }
    }

    *processcount = pcnt;
    *processlist = plist;

CLEANUP:
    CCelim_free_graph (&G);
    CC_IFFREE (pcounts, int);
    if (pairs) {
        for (i = 0; i < ncount; i++) { CC_IFFREE (pairs[i], int); }
        CC_FREE (pairs, int *);
    }
    return rval;
}

/* After a run of pair elimination, we have a list of nonpair triples  */
/* that need to be combined with the existing nonpair list. This func  */
/* goes through all pairs in G and saves those that are either in the  */
/* new triple list or in the old nonpair list.  Note that pairs no     */
/* longer in G (due to the fact that one of the edges has been         */
/* eliminated) will not be saved.                                      */

int CCelim_pairlist_to_nonpairs (CCelim_graph *G, int plistcount, int *plist,
        int *nonpcounts, int **nonpairs, int **outpcounts, int ***outpairs)
{
    int rval = 0, i, *pcounts = (int *) NULL, **pairs = (int **) NULL;
    int *acounts = (int *) NULL, **apairs = (int **) NULL, total = 0;
    int ncount = G->ncount, j, *atmp = (int *) NULL, mdeg = 0, a, b, c, tmp;
    int verbose = 0;
    CCelim_node *pa;

    /* First convert triples in plist to list of pairs for each node.  */
    /* The first node in the triple is the common end of the pair.     */

    CC_MALLOC (pcounts, ncount, int);
    for (a = 0; a < ncount; a++) pcounts[a] = 0;
    for (i = 0; i < plistcount; i++) pcounts[plist[3*i]]++;

    CC_MALLOC (pairs, ncount, int *);
    for (a = 0; a < ncount; a++) pairs[a] = (int *) NULL; 

    for (a = 0; a < ncount; a++) {
        if (pcounts[a] > 0) {
            CC_MALLOC (pairs[a], 2*pcounts[a], int);
            pcounts[a] = 0;
        }
    }

    for (i = 0; i < plistcount; i++) {
        a = plist[3*i]; b = plist[3*i+1]; c = plist[3*i+2];
        pairs[a][2*pcounts[a]]   = b;
        pairs[a][2*pcounts[a]+1] = c;
        pcounts[a]++;
    }

    /* acounts+apairs will be the output lists; the atmp array is long   */
    /* enough to handle all pairs for any node; we will copy its content */
    /* into apairs[a] for each node a                                    */

    CC_MALLOC (acounts, ncount, int);
    CC_MALLOC (apairs, ncount, int *);
    for (a = 0; a < ncount; a++) {
        apairs[a] = (int *) NULL;
        if (G->nodelist[a].deg > mdeg) mdeg = G->nodelist[a].deg;
    }
    CC_MALLOC (atmp, 2*mdeg*mdeg, int);

    for (a = 0; a < ncount; a++) {
        pa = &G->nodelist[a];
        acounts[a] = 0;
        /* run through all pairs of neighbors (b,c) of a */
        for (i = 0; i < pa->deg; i++) {
            b = pa->neighbors[i];
            for (j = i+1; j < pa->deg; j++) {
                c = pa->neighbors[j];
                /* check if (b,c) is in one of the two lists for node a */
                if (pair_in_list (pcounts[a], pairs[a], b, c) || (nonpairs &&
                    pair_in_list (nonpcounts[a], nonpairs[a], b, c))) {
                    atmp[2*acounts[a]] = b;
                    atmp[2*acounts[a]+1] = c;
                    acounts[a]++;
                }
            }
        }
        if (acounts[a] > 0) {
            CC_MALLOC (apairs[a], 2*acounts[a], int);
            for (i = 0; i < acounts[a]; i++) {
                b = atmp[2*i];  c = atmp[2*i+1];
                if (b > c) CC_SWAP (b, c, tmp);
                apairs[a][2*i]   = b;
                apairs[a][2*i+1] = c;
            }
            total += acounts[a];
        }
    }

    *outpcounts = acounts;
    *outpairs = apairs;

    if (verbose) {
        printf ("Total Non-Pairs: %d\n", total); fflush (stdout);
    }

CLEANUP:
    CC_IFFREE (pcounts, int);
    CC_IFFREE (atmp, int);
    if (pairs) {
        for (i = 0; i < ncount; i++) { CC_IFFREE (pairs[i], int); }
        CC_FREE (pairs, int *);
    }
    return rval;
}

static int pair_in_list (int count, int *list, int a, int b) 
{
    int tmp, i;

    if (a > b) CC_SWAP (a, b, tmp);
    for (i = 0; i < count; i++) {
        if (list[2*i] == a && list[2*i+1] == b) return 1;
    }
    return 0;
}
