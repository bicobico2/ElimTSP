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
/*  Last Update: May 12, 2016; December 4, 2015                             */
/*                                                                          */
/*  Exported Functions                                                      */
/*                                                                          */
/*  void CCelim_try_swaps (int *ab, CCelim_path *psys, int nodecount,       */
/*      int *invnames, int *omatch, int **M, int *yesno, int *imatch,       */
/*      int use_tsp)                                                        */
/*    USE exhaustive search for 3, 4, and 5 swaps to immprove path-tour     */
/*     -ab gives nodes of an edge that should be in the swap                */
/*     -use_tsp specifies k = 3, 4, or 5 (if 1 then use up to 5-swaps)      */
/*     -yesno set to 1 if an improving swap is found                        */
/*                                                                          */
/*  void CCelim_tsp_swap (int nodecount, int *nodes, int pathlen,           */
/*      int pathcount, int *omatch, int **M, int *yesno, int *imatch,       */
/*      CCelim_graph *G)                                                    */
/*    USE exact TSP solver to test for a local improvement                  */
/*     -G is used when testing the Bonn idea for a sparse search (it gives  */
/*        a faster search, but misses many swaps; turn it on by setting     */
/*        full = 0 in the routine)                                          */
/*                                                                          */
/****************************************************************************/

#include <stdio.h>
#include "ccutil.h"
#include "elim.h"

static void next_set (int sz, int *Set);
static void get_imatch (int nodecount, int *tour, int pathcount, int *omatch, 
    int *imatch, int *hit);

#define NMAX CCelim_NMAX
#define PMAX CCelim_PMAX

#define HK_NODELIMIT 1000000  /* 1000000 */

void CCelim_try_swaps (int *ab, CCelim_path *psys, int nodecount, int *invnames,
        int *omatch, int **M, int *yesno, int *imatch, int use_tsp)
{
    int Set[5], del[10], add[10], t[2*NMAX], hit[2*NMAX], tour[2*NMAX];
    int et[2*NMAX], i, j, k, u, v, sz, test, pathcount, tcount = 0, rval;

    *yesno = 0;
    pathcount = CCelim_path_system_count (psys);

    if (use_tsp == 1) use_tsp = 5;

    /* hit is used to identify omatch, linking the paths into a tour     */

    for (i = 0; i < nodecount; i++) hit[i] = 0;
    for (i = 0; i < pathcount; i++) {
        hit[omatch[2*i]] = i+1; 
        hit[omatch[2*i+1]] = i+1;
    }

    /* use local names for nodes in path edges; add omatch to get a tour */

    CCelim_paths_to_edges (psys, &tcount, t);
    for (i = 0; i < 2*tcount; i++) t[i] = invnames[t[i]];

    for (i = 0, k = tcount; i < pathcount; i++, k++) {
        t[2*k] = omatch[2*i];
        t[2*k+1] = omatch[2*i+1];
    }
    rval = CCutil_edge_to_cycle (nodecount, t, &test, tour);
    if (rval || !test) {
        printf ("WARNING: paths+omatch is not a tour!\n"); exit (1);
    }

    /* grab the path edges in tour order; put edge ab first in order */

    if (ab) {  
        int a = invnames[ab[0]], b = invnames[ab[1]];
        for (j = 0; j < nodecount; j++) {
            if ((tour[j] == a && tour[(j+1) % nodecount] == b) ||
                (tour[j] == b && tour[(j+1) % nodecount] == a)) {
                break;
            }
        }
        if (j == nodecount) {
            printf ("WARNING: missing the ab edge\n"); exit (1);
        }
    } else {
        j = 0;
    }
    for (k = 0, i = 0; i < nodecount; i++, j++) {
        u = tour[j % nodecount], v = tour[(j+1) % nodecount];
        if (hit[u] == 0 || hit[v] == 0 || hit[u] != hit[v]) {
            et[2*k] = u; et[2*k+1] = v; k++;
        }
    }

    /* try all possible k-swaps that include the first edge (ab) in list */

    for (k = 3; k <= use_tsp && k < tcount && *yesno == 0; k++) {
        sz = k-1;
        for (i = 0; i < sz; i++) Set[i] = i;

        del[0] = et[0]; del[1] = et[1];      /* delete the first edge */
        for (; Set[sz-1] < tcount-1 && *yesno == 0; next_set(sz,Set)) {
            for (i = 1; i <= sz; i++) {      /* and k-1 other edges   */
                del[2*i]   = et[2*(Set[i-1]+1)];
                del[2*i+1] = et[2*(Set[i-1]+1)+1];
            }
            if (k == 3) {
                CCelim_compare_three_swap (del[0], del[1], del[2], del[3],
                 del[4], del[5], M, yesno, add);
            } else if (k == 4) {
                CCelim_compare_four_swap (del[0], del[1], del[2], del[3],
                 del[4], del[5], del[6], del[7], M, yesno, add);
            } else {
                CCelim_compare_five_swap (del[0], del[1], del[2], del[3],
                 del[4], del[5], del[6], del[7], del[8], del[9], M, yesno, add);
            }
        }
        if (*yesno == 1) {
            int diff = 0;
            for (i = 0; i < k; i++) {
                diff += (M[del[2*i]][del[2*i+1]] - M[add[2*i]][add[2*i+1]]);
            }
            if (diff <= 0) {
                printf ("WARNING: k-swap does not improve tour length\n");
                exit (1);
            }
 
            /* replace deleted edges with added edges and get imatch */

            for (i = 0; i < k; i++) {
                for (j = 0; j < tcount; j++) {
                    if ((et[2*j] == del[2*i]   && et[2*j+1] == del[2*i+1]) ||
                        (et[2*j] == del[2*i+1] && et[2*j+1] == del[2*i])) {
                        et[2*j] = add[2*i]; et[2*j+1] = add[2*i+1];
                        break;
                    }
                }
            }

            /* add omatch to get a tour */

            j = tcount;
            for (j = tcount, i = 0; i < pathcount; i++, j++) {
                et[2*j]   = omatch[2*i];
                et[2*j+1] = omatch[2*i+1];
            }
            rval = CCutil_edge_to_cycle (nodecount, et, &test, tour);
            if (rval || !test) {
                printf ("WARNING: k-swapped to non-tour\n"); exit (1);
            }

            get_imatch (nodecount, tour, pathcount, omatch, imatch, hit);
        }
    }
}

void CCelim_tsp_swap (int nodecount, int *nodes, int pathlen, int pathcount,
        int *omatch, int **M, int *yesno, int *imatch, CCelim_graph *G)
{
    int i, j, k, u, v, target, optlen = 0, found, test, rval;
    int  hit[2*NMAX], edgetour[4*NMAX], tour[2*NMAX];
    int *elist = (int *) NULL, *elen = (int *) NULL, ecount = 0;
    double tval, upbound;
    int full = 1;  /* use complete graph */

    /* the nodes array lists the indices of the nodes in the path system     */
    /* omatch links the paths into a circuit; omatch == outside matching     */
    /* imatch will return ends of new paths; imatch == inside matching       */

    *yesno = 0;

    /* yesno set to 1 if TSP tour replaces path sys at cost < pathlen        */

    if (nodecount <= 3) return;

    CC_MALLOC (elist, nodecount * (nodecount-1), int);
    CC_MALLOC (elen,  (nodecount * (nodecount-1))/2, int);

    /* hit marks edges in omatch with the same index; used for elist         */

    for (i = 0; i < nodecount; i++) hit[i] = 0;
    for (i = 0; i < pathcount; i++) {
        hit[omatch[2*i]] = i+1; 
        hit[omatch[2*i+1]] = i+1;
    }

    /* build elist, giving omatch edges cost 0 and adding pathlen to cost of */
    /* all other edges; this forces the optimal tour to use all omatch edges */

    for (i = 0; i < nodecount; i++) {
        for (j = i+1; j < nodecount; j++) {
            if (hit[i] > 0 && hit[j] > 0 && hit[i] == hit[j]) {
                elist[2*ecount] = i;
                elist[2*ecount+1] = j;
                elen[ecount++] = 0;
            } else if (full || CCelim_edge_in_graph (nodes[i], nodes[j], G)) {
                elist[2*ecount] = i;
                elist[2*ecount+1] = j;
                elen[ecount++] = M[i][j] + pathlen;
            }
        }
    }

    target = ((nodecount - pathcount) + 1) * pathlen;

    /* paths+omatch gives a circuit of cost target */
    
    upbound = (double) target; 
    rval = CCheldkarp_small_elist (nodecount, ecount, elist, elen, &upbound,
                            &tval, &found, 1, edgetour, HK_NODELIMIT, 2);
    if (rval) {  /* likely hit node bound */
        printf ("H%d ", nodecount); fflush (stdout); goto CLEANUP;
    }
    optlen = (int) tval;

    if (optlen < target) {
        if (!found) {
            printf ("WARNING: HK did not return tour!\n"); exit (1);
        }

        rval = CCutil_edge_to_cycle (nodecount, edgetour, &test, tour);
        if (rval || !test) {
            printf ("WARNING: HK solution is not a tour!\n"); exit (1);
        }

        for (k = 0, i = 0; i < nodecount; i++) {
            u = tour[i], v = tour[(i+1) % nodecount];
            if (hit[u] == 0 || hit[v] == 0 || hit[u] != hit[v]) {
                k += (M[u][v] + pathlen);
            }
        }
        if (k != optlen) {
            printf ("WARNING: HK tour has wrong length: %d vs %d (%d)\n",
                     k,optlen,target);
            exit (1);
        }

        get_imatch (nodecount, tour, pathcount, omatch, imatch, hit);
        *yesno = 1;
    }

CLEANUP:
    CC_IFFREE (elist, int);
    CC_IFFREE (elen, int);
}

static void get_imatch (int nodecount, int *tour, int pathcount, int *omatch, 
        int *imatch, int *hit)
{
    int i, j, eorder[NMAX];

    /* eorder lists the order in which the tour visits the path ends  */

    for (j = 0, i = 0; i < nodecount; i++) {
        if (hit[tour[i]]) eorder[j++] = tour[i];
    }

    /* test if imatch are the odd edges in eorder or the even edges */

    for (i = 0; i < pathcount; i++) {
        if ((omatch[2*i] == eorder[0] && omatch[2*i+1] == eorder[1]) ||
            (omatch[2*i] == eorder[1] && omatch[2*i+1] == eorder[0])) break;
    }
    if (i == pathcount) {  /* so put the first edge into imatch */
        for (j = 0; j < 2*pathcount; j++)   imatch[j] = eorder[j];
    } else {
        for (j = 0; j < 2*pathcount-1; j++) imatch[j] = eorder[j+1];
        imatch[j] = eorder[0];
    }
}

static void next_set (int sz, int *Set)
{
   int i;
   for (i=0; i < sz-1 && Set[i]+1 == Set[i+1]; i++) Set[i] = i;
   Set[i] = Set[i]+1;
}
