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
/*  Last Update: May 25, 2016; July 28, 2016                                */
/*                                                                          */
/*  Exported Functions                                                      */
/*                                                                          */
/*  void CCelim_init_graph (CCelim_graph *G)                                */
/*    INITIALIZE the graph to null settings                                 */
/*                                                                          */
/*  void CCelim_free_graph (CCelim_graph *G)                                */
/*    FREE the memory in the graph structure                                */
/*                                                                          */
/*  int CCelim_build_graph (CCelim_graph *G, int ncount, int ecount,        */
/*      int *elist, int *elen, (CCelim_distobj *) D)                        */
/*    BUILD the graph struture (elist is in end-end format)                 */
/*    Note: for the main routines, you must specify either elen or D, to    */
/*    allow the graph to contain valid edge lengths; if lengths are not     */
/*    needed, then set elen and D to NULL).                                 */
/*                                                                          */
/*  int CCelim_getedges_graph (CCelim_graph *G, int *ocount, int **olist)   */
/*    RETURN a list of the edges in the graph (end-end format)              */
/*                                                                          */
/*  int CCelim_edge_in_graph (int a, int b, CCelim_graph *G)                */
/*    CHECK if edge (a,b) is in the graph                                   */
/*                                                                          */
/*  void CCelim_delete_edge (int a, int b, CCelim_graph *G)                 */
/*    REMOVE edge (a,b) from graph                                          */
/*                                                                          */
/*  int CCelim_add_edge (int a, int b, CCelim_graph *G)                     */
/*    ADD edge (a,b) to graph; enough room to add two edges to each node    */
/*                                                                          */
/****************************************************************************/

#include <stdio.h>
#include "ccutil.h"
#include "elim.h"

void CCelim_init_graph (CCelim_graph *G)
{
    if (G) {
        G->ncount = 0;
        G->ecount = 0;
        G->nodelist = (CCelim_node *) NULL;
        G->neighborspace = (int *) NULL;
        G->lenspace = (int *) NULL;
    }
}

void CCelim_free_graph (CCelim_graph *G)
{
    if (G) {
        CC_IFFREE (G->nodelist, CCelim_node);
        CC_IFFREE (G->neighborspace, int);
        CC_IFFREE (G->lenspace, int);
        CCelim_init_graph (G);
    }
}

int CCelim_build_graph (CCelim_graph *G, int ncount, int ecount, int *elist,
        int *elen, CCelim_distobj *D)
{
    int rval = 0, i, t, *p, *q;
    CCelim_node *n0, *n1;

    /* leave enough space to possibly add two edges incident with each     */
    /* node; this allows us to use a graph to store the set of fixed edges */

    CCelim_init_graph (G);
    G->ncount = ncount;
    G->ecount = ecount;
    CC_MALLOC (G->nodelist, ncount, CCelim_node);
    CC_MALLOC (G->neighborspace, 2*ecount + 2*ncount, int);
    CC_MALLOC (G->lenspace, 2*ecount + 2*ncount, int);
    for (i = 0; i < ncount; i++) G->nodelist[i].deg = 0;

    for (i = 0; i < ecount; i++) {
        G->nodelist[elist[2*i]].deg++;
        G->nodelist[elist[2*i+1]].deg++;
    } 

    p = G->neighborspace;
    q = G->lenspace;
    for (i = 0; i < ncount; i++) {
        G->nodelist[i].neighbors = p;
        G->nodelist[i].len = q;
        p += (G->nodelist[i].deg + 2);
        q += (G->nodelist[i].deg + 2);
        G->nodelist[i].maxdeg = G->nodelist[i].deg + 2;  /* for fixed edges */
        G->nodelist[i].deg = 0;
    } 

    for (i = 0; i < ecount; i++) {
        if (elen)   t = elen[i];
        else if (D) t = CCelim_dist (elist[2*i], elist[2*i+1], D);
        else        t = 1;

        n0 = &G->nodelist[elist[2*i]];
        n1 = &G->nodelist[elist[2*i+1]];
        n0->len[n0->deg] = t;
        n1->len[n1->deg] = t;
        n0->neighbors[n0->deg++] = elist[2*i+1];
        n1->neighbors[n1->deg++] = elist[2*i];
    }

CLEANUP:
    return rval;
}

int CCelim_getedges_graph (CCelim_graph *G, int *ocount, int **olist)
{
    int rval = 0, b, i, j, k, ecount = 0, *elist = (int *) NULL;

    *ocount = 0;
    *olist = (int *) NULL;

    for (i = 0, k = 0; i < G->ncount; i++) {
        k += G->nodelist[i].deg;
    }
    ecount = k/2;
    CC_MALLOC (elist, 2*ecount, int);

    for (i = 0, k = 0; i < G->ncount; i++) {
        for (j = 0; j < G->nodelist[i].deg; j++) {
            b = G->nodelist[i].neighbors[j];
            if (i < b) {
                elist[2*k] = i; elist[2*k+1] = b; k++;
            }
        }
    }
    if (k != ecount) {
        fprintf (stderr, "Lost an edge in the graph!\n");
        rval = 1; goto CLEANUP;
    }

    *ocount = ecount;
    *olist = elist;

CLEANUP:
    return rval;
}

int CCelim_edge_in_graph (int a, int b, CCelim_graph *G)
{
    int i;

    for (i = 0; i < G->nodelist[a].deg; i++) {
        if (G->nodelist[a].neighbors[i] == b) return 1;
    }
    return 0; 
}

void CCelim_delete_edge (int a, int b, CCelim_graph *G)
{
    int i, j, p[2];
    CCelim_node *n;

    p[0] = a; p[1] = b;

    for (j = 0; j < 2; j++) {
        n = &G->nodelist[p[j]];
        for (i = 0; i < n->deg; i++) {
            if (n->neighbors[i] == p[1-j]) {
                n->deg--;
                n->neighbors[i] = n->neighbors[n->deg];
                n->len[i] = n->len[n->deg];
                break;
            }
        }
    }
    G->ecount--;
}

int CCelim_add_edge (int a, int b, int len, CCelim_graph *G)
{
    int rval = 0;
    CCelim_node *pa = &G->nodelist[a], *pb = &G->nodelist[b];

    rval = CCelim_edge_in_graph (a, b, G);
    if (rval) {
        printf ("cannot add edge already in graph\n"); goto CLEANUP;
    }
    if (pa->deg >= pa->maxdeg || pb->deg >= pb->maxdeg) {
        fprintf (stderr, "no room to add edge\n");
        rval = 1; goto CLEANUP;
    }

    pa->len[pa->deg] = pb->len[pb->deg] = len;
    pa->neighbors[pa->deg++] = b;
    pb->neighbors[pb->deg++] = a;
    G->ecount++;

CLEANUP:
    return rval;
}
