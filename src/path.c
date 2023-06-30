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
/*  Last Update: May 27, 2016                                               */
/*                                                                          */
/*  Exported Functions                                                      */
/*                                                                          */
/*  int CCelim_validate_paths (int ncount, CCelim_path *psys, int *yesno,   */
/*      CCelim_path *pmerge)                                                */
/*    CHECK that psys gives internally disjoint paths and no cycles         */
/*     -pmerge if not NULL will return a system with paths merged           */
/*                                                                          */
/*  int CCelim_test_paths (int *ab, CCelim_path *psystem,                   */
/*      CCelim_distobj *D, CCelim_graph *G, int use_tsp, int *bad)          */
/*    TEST that every tour containing the paths can be improved             */
/*    -bad is set to 1 if some tour orientation cannot be improved          */
/*    -use_tsp should be set to 1 to use exact TSP search, or k = 3, 4, 5   */
/*     for k-opt search                                                     */
/*                                                                          */
/*  int CCelim_path_system_count (CCelim_path *p)                           */
/*   RETURNS number of paths in the system                                  */
/*                                                                          */
/*  void CCelim_copy_path_system (CCelim_path *porig, CCelim_path *pcopy,   */
/*      int reverse)                                                        */
/*   COPIES the path system; if reverse = 1 then reverse each path          */
/*                                                                          */
/*  void CCelim_paths_to_edges (CCelim_path *pathsystem, int *tcount,       */
/*      int *t)                                                             */
/*   RETURNS in array t the edges in the paths in node-node format          */
/*                                                                          */
/*  void CCelim_path_system_nodeset (CCelim_path *psys, int *nodecount,     */
/*      int *nodeset)                                                       */
/*   RETURNS in array nodeset the nodes in the path system                  */ 
/*                                                                          */
/*  void CCelim_path_system_mark (CCelim_path *psys, int *marks,            */
/*      int marker)                                                         */
/*   MARKS the nodes in the path system with value of marker                */
/*                                                                          */
/****************************************************************************/

#include <stdio.h>
#include "ccutil.h"
#include "elim.h"

static void  check_for_circuit (int ncount, int ecount, int *elist, int *circ,
    CCelim_path *pmerge);
static int next_perm (int n, int *p);
static void array_reverse (int x, int y, int *p);
static void compare_path_systems (CCelim_path *pa, CCelim_path *pb, int *yesno);
static void check_match (int ncount, int pcount, int *omatch, int *imatch,
    int *yesno);
static int build_omatch_list (int pathcount, CCelim_path **pathlist,
    int *invnames, int *pocount, int **polist);

#define NMAX CCelim_NMAX
#define PMAX CCelim_PMAX
#define IMAX 100000
#define PMAX_TEST 7  /* max number of paths to test for TSP optimality */

int CCelim_validate_paths (int ncount, CCelim_path *psys, int *yesno,
        CCelim_path *pmerge)
{
    int *hit = (int *) NULL, i, pcount = 0, namer = 0, pelist[2*NMAX];
    int rval = 0, icircuit = 0, names[2*NMAX], a, b, j, k, newnodes[2*NMAX];
    CCelim_path *p, *q;

    *yesno = 1;

    CC_MALLOC (hit, ncount, int);

    for (p = psys; p; p = p->next) {
        for (i = 0; i < p->cnt; i++) hit[p->node[i]] = 0;
        pcount++;
    }

    for (p = psys; p; p = p->next) {
        for (i = 0; i < p->cnt; i++) {
            if (hit[p->node[i]]) {
                *yesno = 0;  /* repeated node in a path */
                goto CLEANUP;
            } else {
                hit[p->node[i]] = 1;
            }
        }
        for (i = 0; i < p->cnt; i++) hit[p->node[i]] = 0;
    }
             
    for (p = psys; p; p = p->next) {
        if (hit[p->node[0]] == 2 || hit[p->node[p->cnt-1]] == 2) {
            *yesno = 0;  /* end node of degree 3 */
            goto CLEANUP;
        }
        hit[p->node[0]]++;
        hit[p->node[p->cnt-1]]++;
    }

    for (p = psys; p; p = p->next) {
        for (i = 1; i < p->cnt-1; i++) {
            if (hit[p->node[i]]) {
                *yesno = 0;  /* interior node appears elsewhere */
                goto CLEANUP;
            } else {
                hit[p->node[i]] = 1;
            }
        }
    }

    CCelim_path_system_mark (psys, hit, -1);

    namer = 0;
    for (p = psys; p; p = p->next) {
        if (hit[p->node[0]] == -1) {
            names[namer] = p->node[0];
            hit[p->node[0]] = namer++;
        }
        if (hit[p->node[p->cnt-1]] == -1) {
            names[namer] = p->node[p->cnt-1];
            hit[p->node[p->cnt-1]] = namer++;
        }
    }

    for (i = 0, p = psys; p; p = p->next, i++) {
        pelist[2*i]   = hit[p->node[0]];
        pelist[2*i+1] = hit[p->node[p->cnt-1]];
    }

    check_for_circuit (namer, pcount, pelist, &icircuit, pmerge);
    if (icircuit == 1) {
        *yesno = 0;  /* paths contain a cirucit */
        goto CLEANUP;
    }

    if (pmerge) {
        /* check_for_circuit has made the merge working with a graph made up */
        /* of edges from the ends of the paths; now expand the graph edges   */
        /* back to the original paths (checking both a to b and b to a)      */

        for (p = pmerge; p; p = p->next) {
            k = 0;
            for (i = 0; i < p->cnt-1; i++) {
                a = names[p->node[i]];
                b = names[p->node[i+1]];
                newnodes[k++] = a;
                for (q = psys; q; q = q->next) {
                     if (q->node[0] == a && q->node[q->cnt-1] == b) {
                         for (j = 1; j < q->cnt-1; j++) {
                             newnodes[k++] = q->node[j];
                         }
                         break;
                     }
                     if (q->node[0] == b && q->node[q->cnt-1] == a) {
                         for (j = q->cnt-2; j > 0; j--) {
                             newnodes[k++] = q->node[j];
                         }
                         break;
                     }
                }
                if (!q) { printf ("Did not find ends of path\n"); exit (1); }
            }
            newnodes[k++] = names[p->node[p->cnt-1]];
            for (i = 0; i < k; i++) { p->node[i] = newnodes[i]; }
            p->cnt = k;
        }
    }

CLEANUP:
   CC_IFFREE (hit, int);
   return rval;
}

int CCelim_test_paths (int *ab, CCelim_path *psystem, CCelim_distobj *D,
        CCelim_graph *G, int use_tsp, int *bad)
{
    int rval = 0, i, j, k, pathcount = 0, valid = 1, *done = (int *) NULL;
    int nodeset[NMAX], nodecount = 0, pathlen = 0, *invnames = (int *) NULL;
    int ocount = 0, *olist = (int *) NULL, *omatch, imatch[2*PMAX];
    int **M = (int **) NULL;
    CCelim_path *p, *plist[PMAX], pmerge[PMAX];

    /* bad set to 1 if paths are valid but some tour cannot be improved */

    *bad = 0;

    rval = CCelim_validate_paths (G->ncount, psystem, &valid, pmerge);
    CCcheck_rval (rval, "CCelim_validate_paths failed");
    if (!valid) goto CLEANUP;

    if (CCelim_path_system_count (pmerge) > PMAX_TEST) {
        /* system has too many paths to check TSP optimality */
        *bad = 1; goto CLEANUP;  /* BICO 221030 */
    }

    /* Note: comparison is a debugging check, can be removed */
    compare_path_systems (psystem, pmerge, &i);
    if (i == 0) { printf ("Merged path is not same as original\n"); exit (1); }

    for (p = pmerge; p; p = p->next) { plist[pathcount++] = p; }

    /* gather info for tsp_swap: nodeset, total path len, distance matrix */

    CCelim_path_system_nodeset (pmerge, &nodecount, nodeset);
    CC_MALLOC (invnames, G->ncount, int);
    for (i = 0; i < nodecount; i++) invnames[nodeset[i]] = i;

    for (p = pmerge; p; p = p->next) {
        for (k = 0; k < p->cnt-1; k++) {
            pathlen += CCelim_dist (p->node[k], p->node[k+1], D);
        }
    }

    CC_MALLOC (M, nodecount, int *);
    for (i = 0; i < nodecount; i++) { M[i] = (int *) NULL; }
    for (i = 0; i < nodecount; i++) { CC_MALLOC (M[i], nodecount, int); }

    for (i = 0; i < nodecount; i++) {
        for (j = i+1; j < nodecount; j++) {
            M[i][j] = M[j][i] = CCelim_dist (nodeset[i], nodeset[j], D);
        }
        M[i][i] = 0;
    }

    rval  = build_omatch_list (pathcount, plist, invnames, &ocount, &olist);
    CCcheck_rval (rval, "build_omatch_list failed");

    CC_MALLOC (done, ocount, int);
    for (i = 0; i < ocount; i++) done[i] = 0;

    for (i = 0; i < ocount; i++) {
        int mhit = 0, good = 0;

        if (done[i]) continue;   /* already covered this matching */

        omatch = olist + (2*i*pathcount);

        if (!mhit) {
            CCelim_try_swaps (ab, pmerge, nodecount, invnames, omatch, M,
                              &good, imatch, use_tsp);
            if (good == 0 && use_tsp == 1) {
                CCelim_tsp_swap (nodecount, nodeset, pathlen, pathcount,
                                 omatch, M, &good, imatch, G);
            }

            if (good == 1) {
                /* check if imatch covers other matchings in our list */
                for (j = i+1; j < ocount; j++) {
                    if (done[j] == 0) {
                        check_match (nodecount, pathcount,
                                     olist + (2*j*pathcount), imatch, &k);
                        if (k) done[j] = 1;
                    }
                }
            } else {
                *bad = 1;
                goto CLEANUP;  /* cound not improve the path system */
            }
        }
    }

CLEANUP:
    CC_IFFREE (invnames, int);
    CC_IFFREE (olist, int);
    CC_IFFREE (done, int);
    if (M) {
        for (i = 0; i < nodecount; i++) CC_IFFREE (M[i], int);
        CC_IFFREE (M, int *);
    }
    return rval;
}

static int build_omatch_list (int pathcount, CCelim_path **pathlist,
    int *invnames, int *pocount, int **polist)
{
    int rval = 0, i, s, perm[NMAX];
    int ocount, *olist = (int *) NULL, *e;
    CCelim_path *p;

    *pocount = 0;
    *polist = (int *) NULL;

    /* number of outside matchings is (pathcount-1)! * 2^(pathcount-1) */

    ocount = (1 << (pathcount-1));
    for (i = pathcount-1; i >= 2; i--) ocount = i * ocount;

    /* each matching has size 2*pathcount */
    
    CC_MALLOC (olist, 2*pathcount*ocount, int);

    /* entries in matchings are indices into the path system nodeset     */
    /* olist is an int array, listing the matchings, one after the other */

    if (pathcount == 1) {
        olist[0] = invnames[pathlist[0]->node[0]];
        olist[1] = invnames[pathlist[0]->node[pathlist[0]->cnt-1]];
        ocount = 1;
        goto CLEANUP;
    }

    ocount = 0;  
    e = olist;   /* e points to next available space for an edge in olist */

    for (i = 0; i < pathcount-1; i++) { perm[i] = i;}
    do {
        for (s = 0; s < (1<<(pathcount-1)); s++) {

            /* the first path always has forward orientation; bits of s   */
            /* are used to run through orientations of the other paths    */

            /* first matching edge starts with the end node of path 0     */
            p = pathlist[0];
            *(e++) = invnames[p->node[p->cnt-1]];

            for (i = 0; i < pathcount-1; i++) {
                p = pathlist[perm[i]+1];  /* p is next path in the tour   */
                if (s & (1<<i)) {         /* bit is set, so reverse p     */
                    /* 2nd node in matching edge i is final node of p     */
                    *(e++) = invnames[p->node[p->cnt-1]]; 
                    /* 1st node in matching edge i+1 is first node of p   */
                    *(e++) = invnames[p->node[0]];
                } else {                  /* bit is not set, so p forward */
                    /* 2nd node in matching edge i is first node of p     */
                    *(e++) = invnames[p->node[0]];
                    /* 1st node in matching edge i+1 is final node of p   */
                    *(e++) = invnames[p->node[p->cnt-1]];
                }
            }

            /* the last matching edge ends with the start node of path 0  */
            p = pathlist[0];
            *(e++) = invnames[p->node[0]];
            ocount++;
        }
    } while (next_perm (pathcount-1, perm) == 0);

CLEANUP:
    if (rval) {
        CC_IFFREE (olist, int);
    } else {
        *pocount = ocount;
        *polist = olist;
    }
    return rval;
}

static void check_match (int ncount, int pcount, int *omatch, int *imatch,
        int *yesno)
{
    int a[NMAX], b[NMAX], hit[NMAX], i, k, v;

    *yesno = 1;

    for (i = 0; i < ncount; i++) hit[i] = 0;
    for (i = 0; i < pcount; i++) {
        a[omatch[2*i]]   = omatch[2*i+1];
        a[omatch[2*i+1]] = omatch[2*i];
        b[imatch[2*i]]   = imatch[2*i+1];
        b[imatch[2*i+1]] = imatch[2*i];
    }

    v = omatch[0], k = 0;
    while (!hit[v]) {
        hit[v] = 1;
        if (k++ % 2) v = a[v];
        else         v = b[v];
    } 
    if (k != 2*pcount || v != omatch[0]) *yesno = 0;
}

static void check_for_circuit (int ncount, int ecount, int *elist, int *circ,
        CCelim_path *pmerge)
{
    int i, n0, n1, deg[2*NMAX], adj[2*NMAX][2], start, prev, next;
    int k = 0, pcount = 0;
    CCelim_path *p = (CCelim_path *) NULL;

    /* graph specified by the edgelist arises from a list of internally */
    /* disjoint paths, with no node of degree greater than 2            */

    *circ = 0;

    for (i = 0; i < ncount; i++) deg[i] = 0;

    for (i = 0; i < ecount; i++) {
        n0 = elist[2*i];
        n1 = elist[2*i+1];
        adj[n0][deg[n0]++] = n1; 
        adj[n1][deg[n1]++] = n0; 
    }

    for (i = 0; i < ncount; i++) { if (deg[i] != 2) break; }
    if (i == ncount) { *circ = 1; return; }

    for (i = 0; i < ncount; i++) {
        if (deg[i] != 1) continue;
        if (pmerge) {
            p = &pmerge[pcount++];
        }
        start = i;
        prev = i;
        next = adj[i][0];
        if (pmerge) {
            p->node[0] = prev;
            p->node[1] = next;
            k = 2;
        }
        deg[start] = 0; /* marks the visit */
        while (next != start && deg[next] == 2) {
            deg[next] = 0;  /* marks the visit */
            if (adj[next][0] == prev) {
                prev = next;
                next = adj[next][1];
            } else {
                prev = next;
                next = adj[next][0];
            }
            if (pmerge) p->node[k++] = next;
        } 
        if (next == start) {
            *circ = 1; return;
        }
        deg[next] = 0;  /* marks the visit */
        if (pmerge) p->cnt = k;
    }

    for (i = 0; i < ncount; i++) {
        if (deg[i] == 2) {
            *circ = 1; return;
        }
    }

    if (pmerge) {
        for (i = 0; i < pcount-1; i++) pmerge[i].next = &pmerge[i+1];
        pmerge[pcount-1].next = (CCelim_path *) NULL;
    }
}

void CCelim_paths_to_edges (CCelim_path *p, int *tcount, int *t)
{
    int k, ecount = 0;

    for (; p; p = p->next) {
        for (k = 0; k < p->cnt-1; k++) {
            t[2*ecount]   = p->node[k];
            t[2*ecount+1] = p->node[k+1];
            ecount++;
        }
    }
    *tcount = ecount;
}

static int next_perm (int n, int *p)
{
    /* Fill in p with the next lexicographicly larger permuation. */

    int k, l, tmp;

    k = n-2;
    while (k > -1 && p[k] >= p[k+1]) k--;

    if (k == -1) return 1;  /* p is max-lex */

    l = n-1;
    while (p[k] >= p[l]) l--;
    CC_SWAP(p[k], p[l], tmp);

    array_reverse (k+1, n-1, p);
    return 0;
}

static void array_reverse (int x, int y, int *p)
{
    int tmp;

    while (x < y) {
        CC_SWAP (p[x], p[y], tmp);
        x++; y--;
    }
}

void CCelim_copy_path_system (CCelim_path *porig, CCelim_path *pcopy,
        int reverse)
{
    int i, k = 0;
    CCelim_path *p;

    /* reverse only reverses individual paths, not the order of the paths */

    for (p = porig; p; p = p->next) {
        for (i = 0; i < p->cnt; i++) {
            if (reverse == 0) {
                pcopy[k].node[i] = p->node[i];
            } else {
                pcopy[k].node[(p->cnt-1)-i] = p->node[i];
            }
        }
        pcopy[k].cnt = p->cnt;
        pcopy[k].next = &(pcopy[k+1]);
        k++;
    }
    pcopy[k-1].next = (CCelim_path *) NULL;
}

static void compare_path_systems (CCelim_path *pa, CCelim_path *pb, int *yesno)
{
    int acount, a[2*NMAX], bcount, b[2*NMAX], i, j, k, tmp;
    CCelim_path *p;

    /* check that the two path systems cover the same set of edges */

    *yesno = 1;

    CCelim_paths_to_edges (pa, &acount, a);
    CCelim_paths_to_edges (pb, &bcount, b);

    if (acount != bcount) {
        printf ("acount = %d, bcount = %d\n", acount, bcount);
        *yesno = 0; goto CLEANUP;
    }

    for (i = 0; i < acount; i++) {
        if (a[2*i] > a[2*i+1]) {
             CC_SWAP(a[2*i], a[2*i+1], tmp);
        }
        if (b[2*i] > b[2*i+1]) {
             CC_SWAP(b[2*i], b[2*i+1], tmp);
        }
    }

    for (i = 0; i < acount; i++) {
        for (j = 0; j < bcount; j++) {
            if (b[2*j] == a[2*i] && b[2*j+1] == a[2*i+1]) break;
        }
        if (j == bcount) {
            printf ("b is missing an edge from a\n");
            *yesno = 0; goto CLEANUP;
        }
    }

CLEANUP:
    if (*yesno == 0) {
        printf ("PATH A: ");
        for (p = pa; p; p = p->next) {
            printf (" [");
            for (k = 0; k < p->cnt; k++) printf ("%d ", p->node[k]);
            printf ("]");
        }
        printf ("\n");
        printf ("PATH B: ");
        for (p = pb; p; p = p->next) {
            printf (" [");
            for (k = 0; k < p->cnt; k++) printf ("%d ", p->node[k]);
            printf ("]");
        }
        printf ("\n");
        fflush (stdout);
    }
}

int CCelim_path_system_count (CCelim_path *p)
{
    int k = 0;
    for (; p; p = p->next) k++;
    return k;
}

void CCelim_path_system_nodeset (CCelim_path *p, int *nodecount,
        int *nodeset)
{
    int k = 0, count = 0;

    for (; p; p = p->next) {
        for (k = 0; k < p->cnt; k++) {
            nodeset[count++] = p->node[k];
        }
    }
    *nodecount = count;
}

void CCelim_path_system_mark (CCelim_path *p, int *marks, int marker)
{
    int k = 0;

    for (; p; p = p->next) {
        for (k = 0; k < p->cnt; k++) { marks[p->node[k]] = marker; }
    }
}
