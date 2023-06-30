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
/*  Last Update: June-July 2016, October 2017                               */
/*                                                                          */
/*  Exported Functions                                                      */
/*                                                                          */
/*  int CCelim_run_elim_edge (int a, int b, CCelim_graph *G,                */
/*      CCelim_graph *F, CCelim_distobj *D, int wtype,                      */
/*      CCrandstate *rstate, int *yesno, double *ezeit, int level,          */
/*      int *pcounts, int **pairs, CCkdtree *kt, CCdatagroup *euclid_dat,   */
/*      CCkdtree *euclid_kt, double timelimit, int longpath, int use_tsp,   */
/*      int max_neighborhood, CCelim_httree *ht)                            */
/*    ELIMINATE edge (a,b); yesno set to 1 if successful                    */
/*     -a and b are ends of the edge                                        */
/*     -G is the sparse graph containing all remaining edges                */
/*     -F is the graph containing all edges fixed to 1 (it can be NULL)     */
/*     -wtype specifies the type of edge used in the first step of the      */
/*      elimination process                                                 */
/*     -yesno is set to 1 if the edge is eliminated                         */
/*     -ezeit is set to the running time for the elimination routine        */
/*     -level controls the depth of the elimination-search tree             */
/*     -pcounts and pairs can specify the pairs of adjacent edges that can  */
/*      be included in an optimal tour (can be NULL)                        */
/*     -kt can specify a kd-tree to permit neighbor searches (can be NULL)  */
/*     -euclid_dat and euclid_kt can specify Euclidean data to build        */
/*      neighborhoods for non-Euclidean instances (like road examples),     */
/*      (they can be NULL)                                                  */
/*     -timelimit is an upperbound on time for the elimination search       */
/*     -longpath can be set to k >= 3 to make the initial witness family    */
/*      paths of length k containing edge (a,b)                             */
/*     -use_tsp can be set 1 to specify that an exact TSP code should be    */
/*      used in the search for improving k-swaps                            */
/*     -max_neighborhood specifies the number of candidates for witness     */
/*      nodes (a neighborhood around the edge (a,b))                        */
/*     -ht if not NULL then the search tree will be saved                   */
/*                                                                          */
/*  int CCelim_run_fast_elim_edge (int a, int b, CCelim_graph *G,           */
/*      CCelim_distobj *D, CCrandstate *rstate, CCkdtree *kt,               */
/*      CCdatagroup *euclid_dat, CCkdtree *euclid_kt, int level,            */
/*      int max_neighborhood, int *yesno, CCelim_httree *ht)                */
/*    ELIMINATE edge (a,b) using simple methods; it should run faster than  */
/*      the full elimination code; it is meant to be used as an initial     */
/*      process on the edge set obtained from LP-based elimination          */
/*      -ht if not NULL, then the search tree will be saved (in this case   */
/*       the fast split_cd code will not be called, to match verify.c)      */
/*                                                                          */
/*  int CCelim_run_fix_edge (int a, int b, CCelim_graph *G,                 */
/*      CCelim_graph *F, CCelim_distobj *D, int wtype,                      */
/*      CCrandstate *rstate, CCkdtree *kt, CCdatagroup *euclid_dat,         */
/*      CCkdtree *euclid_kt, int *yesno, double *ezeit, int level,          */
/*      int *pcounts, int **pairs, double timelimit, int use_tsp,           */
/*      int max_neighborhood, CCelim_httree *ht)                            */
/*    FIX edge (a,b) to 1; yesno set to 1 if successful                     */
/*                                                                          */
/*  int CCelim_run_elim_pair (int a, int b, int c, CCelim_graph *G,         */
/*      CCelim_graph *F, CCelim_distobj *D,                                 */
/*      CCrandstate *rstate, CCkdtree *kt, CCdatagroup *euclid_dat,         */
/*      CCkdtree *euclid_kt, int *yesno, double *ezeit, int level,          */
/*      int *pcounts, int **pairs, double timelimit, int longpath,          */
/*      int use_tsp, int max_neighborhood, Celim_httree *ht)                */
/*    ELIMINATE path a-b-c; yesno set to 1 if successful                    */
/*                                                                          */
/****************************************************************************/

#include <stdio.h>
#include "ccutil.h"
#include "elim.h"

#define MAX_CD  5           /*  5, number of (c,d) witnesses at root of tree */
#define MAX_EXTRA_SEARCH 3  /*  3, number of witnesses at search node */
#define MAX_EXTRA_BRANCH 5  /*  5, number of cases to handle at search node */
#define DEEP_LEVEL 10       /* 10, switch to deep settings */
#define MAX_DEEP_SEARCH 1   /*  1, number of witnesses at deep node */
#define MAX_DEEP_BRANCH 1   /*  1, number of cases to handle at deep node */
#define MAX_PATH_COUNT 5    /*  5, number of paths in system */

#define PMAX CCelim_PMAX
#define NMAX CCelim_NMAX
#define MAX_AB CCelim_MAX_AB

#define CCcheck_doctor if (CCutil_zeit()-startzeit > timelimit) goto CLEANUP;

/* #define CCcheck_doctor */ /* nothing */

static int elim_with_cd (int *ab, int a, int b, int c, int d, CCelim_distobj *D,
    CCelim_graph *G, int *yesno, CCelim_array *ablist, int level,
    CCelim_graph *F, int *pcounts, int **pairs, int *countbad,
    double timelimit, double startzeit, int use_tsp, CCelim_httree *ht,
    CCelim_htnode *ht_parent);
static int elim_via_long_path (int path_len, int *ab, int a, int b,
    CCelim_distobj *D, CCelim_graph *G, int *yesno, CCelim_array *ablist,
    int level, CCelim_graph *F, int *pcounts, int **pairs, double timelimit,
    double startzeit, int use_tsp, CCelim_httree *ht, CCelim_htnode *ht_parent);
static int pair_via_long_path (int path_len, int a, int b, int c,
    CCelim_distobj *D, CCelim_graph *G, int *yesno, CCelim_array *ablist,
    int level, CCelim_graph *F, int *pcounts, int **pairs, double timelimit,
    double startzeit, int use_tsp, CCelim_httree *ht, CCelim_htnode *ht_parent);
static int grow_long_path (int path_len, int *ab, CCelim_distobj *D,
    CCelim_graph *G, int *yesno, CCelim_array *ablist, int level,
    CCelim_graph *F, int *pcounts, int **pairs, CCelim_path *psys,
    double timelimit, double startzeit, int use_tsp, CCelim_httree *ht,
    CCelim_htnode *ht_parent);
static int fix_ab (int a, int b, CCelim_distobj *D, CCelim_graph *G, int *yesno,
    CCelim_array *ablist, int level, CCelim_graph *F, int *pcounts,
    int **pairs, double timelimit, double startzeit, int use_tsp,
    CCelim_httree *ht);
static int try_extra_point (int *ab, CCelim_path *pin, CCelim_distobj *D,
    CCelim_graph *G, CCelim_array *ablist, int *good, int level,
    CCelim_graph *F, int *pcounts, int **pairs, double timelimit,
    double startzeit, int use_tsp, int depth, CCelim_httree *ht,
    CCelim_htnode *ht_parent);
static int try_extra_edge (int *ab, CCelim_path *pin, CCelim_distobj *D,
    CCelim_graph *G, CCelim_array *ablist, int *good, int level,
    CCelim_graph *F, int *pcounts, int **pairs, double timelimit,
    double startzeit, int use_tsp, int depth, CCelim_httree *ht,
    CCelim_htnode *ht_parent);
static int full_path_check (CCelim_path *psys, int *yesno, int *ab,
    CCelim_distobj *D, CCelim_graph *G, CCelim_array *ablist, int level,
    CCelim_graph *F, int *pcounts, int **pairs,
    double timelimit, double startzeit, int use_tsp, int depth,
    CCelim_httree *ht, CCelim_htnode *ht_parent);
static int potential_cd (int a, int b, CCelim_array *ablist, CCelim_graph *G,
    CCelim_distobj *D, int *pcount, int **plist, int wtype, CCelim_graph *F,
    int *pcounts, int **pairs);
static int potential_extra_points (CCelim_array *ablist, CCelim_path *psys,
    CCelim_graph *G, CCelim_distobj *D, CCelim_graph *F, int *pcount,
    int *plist, int *pcounts, int **pairs);
static int rank_extra_points (int *ab, CCelim_path *pin, CCelim_graph *G,
    CCelim_distobj *D, int *pcount, int *plist, int *yesno, CCelim_graph *F,
    int *pcounts, int **pairs, double timelimit, double startzeit, int depth,
    int *winner);
static int compatible_test (int a, int b, int c, int d, CCelim_distobj *D);
static int test_three (int a, int b, int n, int n0, int n1, CCelim_distobj *D);
static int ed2 (int a, int b, CCelim_graph *G);

CC_ELIM_DIST_IN   /* macro builds local copy dist() of CCelim_dist() */

int CCelim_run_elim_edge (int a, int b, CCelim_graph *G, CCelim_graph *F,
        CCelim_distobj *D, int wtype, CCrandstate *rstate, int *yesno,
        double *ezeit, int level, int *pcounts, int **pairs,
        CCkdtree *kt, CCdatagroup *euclid_dat, CCkdtree *euclid_kt,
        double timelimit, int longpath, int use_tsp, int max_neighborhood,
        CCelim_httree *ht)
{
    int rval = 0, j, cdcount = 0, c, d, *cdlist = (int *) NULL, ab[2];
    int countbad, *cdrank = (int *) NULL, *perm = (int *) NULL, k;
    double qzeit = CCutil_zeit(), startzeit = CCutil_zeit();
    CCelim_array ablist;
    CCelim_htnode *root = (CCelim_htnode *) NULL;

    *yesno = 0;
    ablist.arr = (int *) NULL;
    if (F && CCelim_edge_in_graph (a, b, F)) goto CLEANUP;  /* fixed edge */
    ab[0] = a; ab[1] = b;

    if (ht) {
        ht->elimtype = CCelim_HTTREE_EDGE;
        ht->level = level;
        ht->longpath = longpath;
        ht->use_tsp = use_tsp;
        ht->max_neighborhood = max_neighborhood;

        CC_MALLOC (ht->root, 1, CCelim_htnode);
        CCelim_init_htnode (ht->root);
        ht->count = 1;
        ht->depth = 0;
        ht->root->id = 0;
        ht->root->hamilton_type = CCelim_HAMILTON_EDGE;
        ht->root->hamilton_nodes[0] = a;
        ht->root->hamilton_nodes[1] = b;
        root = ht->root;
    }

    /* check edge ab appears in pairs lists for a and b */

    if (pcounts && pairs) {
        int *list = (int *) NULL, acnt, bcnt;
        CC_MALLOC (list, G->nodelist[a].deg + G->nodelist[b].deg, int);
        CCelim_get_partner_pairs (G, a, b, pcounts, pairs, &acnt, list);
        CCelim_get_partner_pairs (G, b, a, pcounts, pairs, &bcnt, list);
        CC_FREE (list, int);

        if (acnt == 0 || bcnt == 0) {
            /* printf ("Can eliminate [%d,%d], not in pairs list\n", a, b); */
            *yesno = 1; goto CLEANUP;
        }
    }

    /* create the neighborhood around edge ab */

    rval = CCelim_build_ablist (a, b, G, D, &ablist, max_neighborhood, kt,
                                rstate, euclid_dat, euclid_kt);  
    CCcheck_rval (rval, "CCelim_build_ablist failed");
    CCcheck_doctor

    if (longpath) {
        if (ht) { 
            ht->root->tutte_type = CCelim_TUTTE_PATH;
            rval = elim_via_long_path (longpath, ab, a, b, D, G, yesno, &ablist,
                level, F, pcounts, pairs, timelimit, startzeit, use_tsp, ht,
                root);
        } else {
            rval = elim_via_long_path (longpath, ab, a, b, D, G, yesno, &ablist,
                level, F, pcounts, pairs, timelimit, startzeit, use_tsp,
                (CCelim_httree *) NULL, (CCelim_htnode *) NULL);
        }
        CCcheck_rval (rval, "elim_via_long_path failed");
        goto CLEANUP;
    }

    rval = potential_cd (a, b, &ablist, G, D, &cdcount, &cdlist, wtype, F,
                         pcounts, pairs);
    CCcheck_rval (rval, "potential_cd failed");
    CCcheck_doctor

    if (cdcount != 0) {
        if (cdcount > 2*MAX_CD) cdcount = 2*MAX_CD; /* Bico 160730 */
        CC_MALLOC (cdrank, cdcount, int);
        CC_MALLOC (perm, cdcount, int);
        for (j = 0; j < cdcount; j++) {
            c = cdlist[2*j]; d = cdlist[2*j+1];
            rval = elim_with_cd (ab, a, b, c, d, D, G, &k, &ablist, level, F,
                  pcounts, pairs, &countbad, timelimit, startzeit, use_tsp,
                  (CCelim_httree *) NULL, (CCelim_htnode *) NULL);
            CCcheck_rval (rval, "elim_with_cd failed");
            CCcheck_doctor
            cdrank[j] = countbad;
        }
        CCutil_int_perm_quicksort (perm, cdrank, cdcount);

        for (j = 0; j < cdcount && j < MAX_CD && !(*yesno); j++) {
            c = cdlist[2*perm[j]]; d = cdlist[2*perm[j]+1];
            if (ht) {
                ht->root->tutte_type = CCelim_TUTTE_CD;
                ht->root->tutte_nodes[0] = c;
                ht->root->tutte_nodes[1] = d;
            }
            rval = elim_with_cd (ab, a, b, c, d, D, G, yesno, &ablist, level,
                  F, pcounts, pairs, (int *) NULL, timelimit, startzeit,
                  use_tsp, ht, root);
            CCcheck_rval (rval, "elim_with_cd failed");
            CCcheck_doctor
        }
    }

    if (*yesno == 0) {
        if (ht) ht->root->tutte_type = CCelim_TUTTE_PATH;
        rval = elim_via_long_path (3, ab, a, b, D, G, yesno, &ablist, level, F,
              pcounts, pairs, timelimit, startzeit, use_tsp, ht, root);
        CCcheck_rval (rval, "elim_via_longpath failed");
        CCcheck_doctor
    }

CLEANUP:
    if (ht && *yesno == 0) {
        CC_IFFREE (ht->root, CCelim_htnode);
        ht->count = 0;
    }
    *ezeit = CCutil_zeit() - qzeit;
    CC_IFFREE (cdlist, int);
    CC_IFFREE (ablist.arr, int);
    CC_IFFREE (cdrank, int);
    CC_IFFREE (perm, int);
    return rval;
}

int CCelim_run_fast_elim_edge (int a, int b, CCelim_graph *G, 
        CCelim_distobj *D, CCrandstate *rstate, CCkdtree *kt,
        CCdatagroup *euclid_dat, CCkdtree *euclid_kt, int level,
        int max_neighborhood, int *yesno, CCelim_httree *ht)
{
    int rval = 0, cdcount = 0, c, d, *cdlist = (int *) NULL, ab[2];
    int *pcounts = (int *) NULL, **pairs = (int **) NULL;
    CCelim_graph *F = (CCelim_graph *) NULL;
    double startzeit = CCutil_zeit(), timelimit = 10000000.0;
    CCelim_array ablist;
    CCelim_htnode *root = (CCelim_htnode *) NULL;
   
    *yesno = 0;
    ablist.arr = (int *) NULL;
    ab[0] = a; ab[1] = b;

    if (ht) {
        ht->elimtype = CCelim_HTTREE_EDGE;
        ht->level = level;
        ht->longpath = 0;
        ht->use_tsp = 5;
        ht->max_neighborhood = max_neighborhood;

        CC_MALLOC (ht->root, 1, CCelim_htnode);
        CCelim_init_htnode (ht->root);
        ht->count = 1;
        ht->depth = 0;
        ht->root->id = 0;
        ht->root->hamilton_type = CCelim_HAMILTON_EDGE;
        ht->root->hamilton_nodes[0] = a;
        ht->root->hamilton_nodes[1] = b;
        root = ht->root;
    }

    rval = CCelim_build_ablist (a, b, G, D, &ablist, max_neighborhood, kt,
                                rstate, euclid_dat, euclid_kt);
    CCcheck_rval (rval, "CCelim_build_ablist failed");

    rval = potential_cd (a, b, &ablist, G, D, &cdcount, &cdlist,
                         CCelim_CD_EDGE, F, pcounts, pairs);
    CCcheck_rval (rval, "potential_cd failed");

    if (cdcount) {
        c = cdlist[0];  d = cdlist[1];
        if (ht) {
            ht->root->tutte_type = CCelim_TUTTE_CD;
            ht->root->tutte_nodes[0] = c;
            ht->root->tutte_nodes[1] = d;
        }
        rval = elim_with_cd (ab, a, b, c, d, D, G, yesno, &ablist, level, F,
              pcounts, pairs, (int *) NULL, timelimit, startzeit, 5,
              ht, root);
        CCcheck_rval (rval, "elim_with_cd failed");
    }

CLEANUP:
    if (ht && *yesno == 0) {
        CC_IFFREE (ht->root, CCelim_htnode);
        ht->count = 0;
    }
    CC_IFFREE (cdlist, int);
    CC_IFFREE (ablist.arr, int);
    return rval;
}

int CCelim_run_fix_edge (int a, int b, CCelim_graph *G, CCelim_graph *F,
        CCelim_distobj *D, CCrandstate *rstate, CCkdtree *kt,
        CCdatagroup *euclid_dat, CCkdtree *euclid_kt, int *yesno,
        double *ezeit, int level, int *pcounts, int **pairs, double timelimit,
        int use_tsp, int max_neighborhood, CCelim_httree *ht)
{
    int rval = 0;
    double qzeit = CCutil_zeit(), startzeit = CCutil_zeit();
    CCelim_array ablist;

    *yesno = 0;
    ablist.arr = (int *) NULL;

    if (ht) {
        ht->elimtype = CCelim_HTTREE_FIX;
        ht->level = level;
        ht->longpath = 0;
        ht->use_tsp = use_tsp;
        ht->max_neighborhood = max_neighborhood;

        CC_MALLOC (ht->root, 1, CCelim_htnode);
        CCelim_init_htnode (ht->root);
        ht->count = 1;
        ht->depth = 0;
        ht->root->id = 0;
        ht->root->hamilton_type = CCelim_HAMILTON_EDGE;
        ht->root->hamilton_nodes[0] = a;
        ht->root->hamilton_nodes[1] = b;
    }

    if (F && CCelim_edge_in_graph (a, b, F)) { 
        *yesno = 1; goto CLEANUP; 
    }

    rval = CCelim_build_ablist (a, b, G, D, &ablist, max_neighborhood, kt,
                                rstate, euclid_dat, euclid_kt);
    CCcheck_rval (rval, "CCelim_build_ablist failed");
    rval = fix_ab (a, b, D, G, yesno, &ablist, level, F, pcounts,
                   pairs, timelimit, startzeit, use_tsp, ht);
    CCcheck_rval (rval, "fix_ab failed");

CLEANUP:
    if (ht && *yesno == 0) {
        CC_IFFREE (ht->root, CCelim_htnode);
        ht->count = 0;
    }
    *ezeit = CCutil_zeit() - qzeit;
    CC_IFFREE (ablist.arr, int);
    return rval;
}

int CCelim_run_elim_pair (int a, int b, int c, CCelim_graph *G,
        CCelim_graph *F, CCelim_distobj *D, CCrandstate *rstate,
        CCkdtree *kt, CCdatagroup *euclid_dat, CCkdtree *euclid_kt,
        int *yesno, double *ezeit, int level, int *pcounts, int **pairs,
        double timelimit, int longpath, int use_tsp, int max_neighborhood,
        CCelim_httree *ht)
{
    int rval = 0, i, j, k, d, good, depth = 0;
    int dcount = 0, *dlist = (int *) NULL, *p;
    CCelim_node *pb = &G->nodelist[b];
    CCelim_path psys[PMAX];
    CCelim_array ablist;
    double qzeit = CCutil_zeit(), startzeit = CCutil_zeit();;
    CCelim_htnode *ht_child = (CCelim_htnode *) NULL;

    *yesno = 0;
    ablist.arr = (int *) NULL;

    if (ht) {
        ht->elimtype = CCelim_HTTREE_PAIR;
        ht->level = level;
        ht->longpath = longpath;
        ht->use_tsp = use_tsp;
        ht->max_neighborhood = max_neighborhood;

        CC_MALLOC (ht->root, 1, CCelim_htnode);
        CCelim_init_htnode (ht->root);
        ht->count = 1;
        ht->depth = 0;
        ht->root->id = 0;
        ht->root->hamilton_type = CCelim_HAMILTON_PATH;
        ht->root->hamilton_nodes[0] = a;
        ht->root->hamilton_nodes[1] = b;
        ht->root->hamilton_nodes[2] = c;
    }

    if (!ht && F && F->nodelist[b].deg > 0) {
        /* Bico 220509 Skip this if we have ht, since we need a tree proof */
        if (!CCelim_edge_in_graph (a, b, F) &&
            !CCelim_edge_in_graph (b, c, F)) {
            *yesno = 1; goto CLEANUP;
        }
    }

    rval = CCelim_build_ablist (a, c, G, D, &ablist, max_neighborhood,
                                kt, rstate, euclid_dat, euclid_kt);  
    CCcheck_rval (rval, "build_ablist failed");

    if (longpath) {
        if (ht) { 
            ht->root->tutte_type = CCelim_TUTTE_PATH;
            rval = pair_via_long_path (longpath, a, b, c, D, G, yesno, &ablist,
                level, F, pcounts, pairs, timelimit, startzeit, use_tsp,
                ht, ht->root);
        } else {
            rval = pair_via_long_path (longpath, a, b, c, D, G, yesno, &ablist,
                level, F, pcounts, pairs, timelimit, startzeit, use_tsp,
                (CCelim_httree *) NULL, (CCelim_htnode *) NULL);
        }
        CCcheck_rval (rval, "pair_via_path failed");
        goto CLEANUP;
    }

    CC_MALLOC (dlist, 2*CCelim_MAX_DEGREE*CCelim_MAX_DEGREE, int);
    psys[0].next = &(psys[1]);
    psys[1].next = (CCelim_path *) NULL;
    psys[0].cnt = 3;
    psys[1].cnt = 3;
    psys[0].node[0] = a; psys[0].node[1] = b; psys[0].node[2] = c;

    /* try as a witness a neighbor d of the middle node b; this fits */
    /* nicely since already b has degree 2 via the path abc          */

    for (i = 0; i < pb->deg; i++) {
        d = pb->neighbors[i];
        if (d == a || d == c || G->nodelist[d].deg > CCelim_MAX_DEGREE)
            continue;
        CCelim_pair_neighbors_three_swap (a, b, c, d, G, D, &dcount, dlist,
                                          pcounts, pairs);
        good = 1;
        if (ht) {
            ht->root->tutte_type = CCelim_TUTTE_POINT;
            ht->root->tutte_nodes[0] = d;
        }
        for (j = 0; j < dcount && good; j++) {
            psys[1].node[0] = dlist[2*j];
            psys[1].node[1] = d;
            psys[1].node[2] = dlist[2*j+1];
            if (ht && ht->count < ht->max_count) {
                CC_MALLOC (ht_child, 1, CCelim_htnode);
                CCelim_init_htnode (ht_child);
                ht_child->id = ht->count++;
                ht_child->hamilton_type = CCelim_HAMILTON_PATH;
                p = ht_child->hamilton_nodes;
                for (k = 0; k < 3; k++) p[k] = psys[1].node[k];
            } else {
                ht_child = (CCelim_htnode *) NULL;
            }
            rval = full_path_check (psys, &good, (int *) NULL, D, G, &ablist,
                  level, F, pcounts, pairs, timelimit, startzeit, use_tsp,
                  depth, ht, ht_child);
            CCcheck_rval (rval, "full_path_check failed");
            CCcheck_doctor
            if (ht && ht_child) {
                if (ht_child->tutte_type != CCelim_TUTTE_NONE) {
                    ht_child->next = ht->root->child;
                    ht->root->child = ht_child;
                    if (depth > ht->depth) ht->depth = depth;
                } else {
                    CC_IFFREE (ht_child, CCelim_htnode);
                    ht->count--;
                }
            }
        }
        if (good == 1) {
            *yesno = 1; goto CLEANUP;
        } else {
            if (ht) {
                CCelim_free_htnode_children (ht->root);
                ht->count = 1;
            }
        }
    }

CLEANUP:
    if (ht && *yesno == 0) {
        CC_IFFREE (ht->root, CCelim_htnode);
        ht->count = 0;
    }
    *ezeit = CCutil_zeit() - qzeit;
    CC_IFFREE (ablist.arr, int);
    CC_IFFREE (dlist, int);
    return rval;
}

/* elim_with_cd attempts to show that edge ab is not in any optimal tour by */
/* building an elimination-search tree starting with ab and all pairs edges */
/* meeting c and all pairs of edges meeting d (so all pairs of pairs)       */

static int elim_with_cd (int *ab, int a, int b, int c, int d, CCelim_distobj *D,
        CCelim_graph *G, int *yesno, CCelim_array *ablist, int level,
        CCelim_graph *F, int *pcounts, int **pairs, int *countbad,
        double timelimit, double startzeit, int use_tsp, CCelim_httree *ht,
        CCelim_htnode *ht_parent)
{
    int rval = 0, ccount, dcount, i, j, k, bad, valid, good, *p;
    int *clist = (int *) NULL, *dlist = (int *) NULL, depth = 2;
    CCelim_path psys[PMAX], pmerge[PMAX];
    CCelim_htnode *ht_child = (CCelim_htnode *) NULL;
    int ht_count = 0;

    *yesno = 0;
    if (ht) ht_count = ht->count;

    CC_MALLOC (clist, 2 * G->nodelist[c].deg * G->nodelist[c].deg, int);
    CC_MALLOC (dlist, 2 * G->nodelist[d].deg * G->nodelist[d].deg, int);
    CCelim_check_neighbors_three_swap (a, b, c, G, D, &ccount, clist, F,
                                       pcounts, pairs);
    CCelim_check_neighbors_three_swap (a, b, d, G, D, &dcount, dlist, F,
                                       pcounts, pairs);
    if ((double) ccount * (double) dcount > 1000000.0) goto CLEANUP;  /* HACK */
    if (ccount == 0 || dcount == 0) {
        *yesno = 1; goto CLEANUP;
    }

    if (level <= -1 && use_tsp != 1 && !countbad && !ht) { 
        /* use specialized code; improves speed in fast elim */
        rval = CCelim_split_cd (a, b, c, d, ccount, clist, dcount, dlist, D,
                                yesno);
        CCcheck_rval (rval, "CCelim_split_cd failed");
        goto CLEANUP;
    }

    psys[0].next = &(psys[1]);
    psys[1].next = &(psys[2]);
    psys[2].next = (CCelim_path *) NULL;
    psys[0].cnt = 2; psys[1].cnt = psys[2].cnt = 3;
    psys[0].node[0] = a; psys[0].node[1] = b;

    if (countbad) {
        int nbad = 0;
        for (i = 0; i < ccount; i++) {
            psys[1].node[0] = clist[2*i];
            psys[1].node[1] = c;
            psys[1].node[2] = clist[2*i+1];
            for (j = 0; j < dcount; j++) {
                psys[2].node[0] = dlist[2*j];
                psys[2].node[1] = d;
                psys[2].node[2] = dlist[2*j+1];
                rval = CCelim_validate_paths (G->ncount, psys, &valid, pmerge);
                CCcheck_rval (rval, "CCelim_validate_paths failed");
                if (valid) {
                    rval = CCelim_test_paths (ab, pmerge, D, G, use_tsp, &bad);
                    CCcheck_rval (rval, "CCelim_test_paths failed");
                    if (bad) nbad++;
                }
            }
        }
        *countbad = nbad;
    } else {
        good = 1;
        for (i = 0; i < ccount && good; i++) {
            psys[1].node[0] = clist[2*i];
            psys[1].node[1] = c;
            psys[1].node[2] = clist[2*i+1];
            for (j = 0; j < dcount && good; j++) {
                psys[2].node[0] = dlist[2*j];
                psys[2].node[1] = d;
                psys[2].node[2] = dlist[2*j+1];
                if (ht && ht_parent && ht->count < ht->max_count) {
                    CC_MALLOC (ht_child, 1, CCelim_htnode);
                    CCelim_init_htnode (ht_child);
                    ht_child->id = ht->count++;
                    ht_child->hamilton_type = CCelim_HAMILTON_PAIR;
                    p = ht_child->hamilton_nodes;
                    for (k = 0; k < 3; k++) p[k]   = psys[1].node[k];
                    for (k = 0; k < 3; k++) p[k+3] = psys[2].node[k];
                } else {
                    ht_child = (CCelim_htnode *) NULL;
                }
                rval = full_path_check (psys, &good, ab, D, G, ablist, level,
                      F, pcounts, pairs, timelimit, startzeit, use_tsp, depth,
                      ht, ht_child);
                CCcheck_rval (rval, "full_path_check failed");
                CCcheck_doctor
                if (ht && ht_parent && ht_child) {
                    if (ht_child->tutte_type != CCelim_TUTTE_NONE) {
                        ht_child->next = ht_parent->child;
                        ht_parent->child = ht_child;
                        if (depth > ht->depth) ht->depth = depth;
                    } else {
                        CC_IFFREE (ht_child, CCelim_htnode);
                        ht->count--;
                    }
                }
            }
        }
        *yesno = good; 
    }

    if (ht && ht_parent && *yesno == 0) {
        CCelim_free_htnode_children (ht_parent);
        ht->count = ht_count;
    }

CLEANUP:
    CC_IFFREE (clist, int);
    CC_IFFREE (dlist, int);
    return rval;
}

/* elim_via_long_path attempts to show that edge ab is not in any tour      */
/* by building an elimination-search tree starting with all paths of length */
/* path_len that contain a-b                                                */

static int elim_via_long_path (int path_len, int *ab, int a, int b,
        CCelim_distobj *D, CCelim_graph *G, int *yesno, CCelim_array *ablist,
        int level, CCelim_graph *F, int *pcounts, int **pairs,
        double timelimit, double startzeit, int use_tsp, CCelim_httree *ht,
        CCelim_htnode *ht_parent)
{
    int rval = 0;
    CCelim_path psys[PMAX];
    int ht_count = 0;

    *yesno = 0;
    if (ht) ht_count = ht->count;

    if (path_len <= 1) goto CLEANUP;

    psys[0].next = (CCelim_path *) NULL;
    psys[0].cnt = 2;
    psys[0].node[0] = a; psys[0].node[1] = b;

    rval = grow_long_path (path_len, ab, D, G, yesno, ablist, level, F,
               pcounts, pairs, psys, timelimit, startzeit, use_tsp, ht,
               ht_parent);
    CCcheck_rval (rval, "grow_long_path failed");

    if (ht && ht_parent && *yesno == 0) {
        CCelim_free_htnode_children (ht_parent);
        ht->count = ht_count;
    }

CLEANUP:
    return rval;
}

/* pair_via_long_path attempts to show that the path a-b-c is incompatible */
/* by building a search three starting with all paths of length path_len   */
/* containing a-b-c                                                        */

static int pair_via_long_path (int path_len, int a, int b, int c,
        CCelim_distobj *D, CCelim_graph *G, int *yesno, CCelim_array *ablist,
        int level, CCelim_graph *F, int *pcounts, int **pairs,
        double timelimit, double startzeit, int use_tsp, CCelim_httree *ht,
        CCelim_htnode *ht_parent)
{
    int rval = 0;
    CCelim_path psys[PMAX];
    int ht_count = 0;

    *yesno = 0;
    if (ht) ht_count = ht->count;

    if (path_len <= 2) goto CLEANUP;

    psys[0].next = (CCelim_path *) NULL;
    psys[0].cnt = 3;
    psys[0].node[0] = a;
    psys[0].node[1] = b;
    psys[0].node[2] = c;

    rval = grow_long_path (path_len, (int *) NULL, D, G, yesno, ablist,
          level, F, pcounts, pairs, psys, timelimit, startzeit, use_tsp,
          ht, ht_parent);
    CCcheck_rval (rval, "grow_long_path failed");

    if (ht && ht_parent && *yesno == 0) {
        CCelim_free_htnode_children (ht_parent);
        ht->count = ht_count;
    }

CLEANUP:
    return rval;
}

/* grow_long_path considers (recursively) all ways to extend psys by */
/* a single edge from one of its end nodes                           */

static int grow_long_path (int path_len, int *ab, CCelim_distobj *D,
        CCelim_graph *G, int *yesno, CCelim_array *ablist, int level,
        CCelim_graph *F, int *pcounts, int **pairs, CCelim_path *psys,
        double timelimit, double startzeit, int use_tsp, CCelim_httree *ht,
        CCelim_htnode *ht_parent)
{
    int rval = 0, good, i, j, c, d, e, ccount = 0, *clist = (int *) NULL;
    int depth = 1;
    CCelim_path preverse[PMAX];
    CCelim_htnode *ht_child = (CCelim_htnode *) NULL;

    *yesno = 0;

    if (psys[0].cnt > path_len) {
        if (ht && ht_parent && ht->count < ht->max_count &&
                  psys[0].cnt < CCelim_HTNODE_MAX_PATH_LEN) {
            CC_MALLOC (ht_child, 1, CCelim_htnode);
            CCelim_init_htnode (ht_child);
            ht_child->id = ht->count++;
            ht_child->hamilton_type = CCelim_HAMILTON_LONG;
            ht_child->hamilton_path_len = psys[0].cnt;
            for (i = 0; i < psys[0].cnt; i++) {
                ht_child->hamilton_nodes[i] = psys[0].node[i];
            }
        } else {
            ht_child = (CCelim_htnode *) NULL;
        }
        rval = full_path_check (psys, yesno, ab, D, G, ablist, level, F,
                      pcounts, pairs, timelimit, startzeit, use_tsp, depth,
                      ht, ht_child);
        CCcheck_rval (rval, "full_path_check failed");
        if (ht && ht_parent && ht_child) {
            if (ht_child->tutte_type != CCelim_TUTTE_NONE) {
                ht_child->next = ht_parent->child;
                ht_parent->child = ht_child;
                if (depth > ht->depth) ht->depth = depth;
            } else {
                CC_IFFREE (ht_child, CCelim_htnode);
                ht->count--;
            }
        }
        goto CLEANUP;
    }

    /* reverse the path to keep a-b in the middle */
    CCelim_copy_path_system (psys, preverse, 1);

    /* grow the path by one edge starting at the last city c */
    c = preverse[0].node[preverse[0].cnt-1];
    d = preverse[0].node[preverse[0].cnt-2];

    /* clist will hold the potential new cities to extend from c */
    CC_MALLOC (clist, G->nodelist[c].deg, int);
    CCelim_get_partner_pairs (G, c, d, pcounts, pairs, &ccount, clist);

    preverse[0].cnt++;
    good = 1;
    for (i = 0; i < ccount && good; i++) {
        e = clist[i];
        for (j = 0; j < preverse[0].cnt-1; j++) {
             if (preverse[0].node[j] == e) break;
        }
        if (j == preverse[0].cnt-1) {
            preverse[0].node[preverse[0].cnt-1] = e;
            rval = grow_long_path (path_len, ab, D, G, &good, ablist,
                     level, F, pcounts, pairs, preverse, timelimit,
                     startzeit, use_tsp, ht, ht_parent);
            CCcheck_rval (rval, "grow_long_path failed");
            CCcheck_doctor
        }
    }
    *yesno = good;

CLEANUP:
    CC_IFFREE (clist, int);
    return rval;
}

/* fix_ab attempts to show that edge ab is in every optimal tour */

static int fix_ab (int a, int b, CCelim_distobj *D, CCelim_graph *G,
        int *yesno, CCelim_array *ablist, int level, CCelim_graph *F,
        int *pcounts, int **pairs, double timelimit, double startzeit,
        int use_tsp, CCelim_httree *ht)
{
    int rval = 0, acount, bcount, i, j, k, good, *ab = (int *) NULL;
    int *alist = (int *) NULL, *blist = (int *) NULL, depth = 2, *p;
    int ht_count = 0;
    CCelim_path psys[PMAX];
    CCelim_htnode *ht_child = (CCelim_htnode *) NULL;

    *yesno = 0;
    if (ht) ht_count = ht->count;

    if (!pairs) {
        CC_MALLOC(alist, 2 * G->nodelist[a].deg * G->nodelist[a].deg,int);
        CC_MALLOC(blist, 2 * G->nodelist[b].deg * G->nodelist[b].deg,int);
        CCelim_all_neighbor_pairs (a, G, F, &acount, alist);
        CCelim_all_neighbor_pairs (b, G, F, &bcount, blist);
    } else {
        alist = pairs[a];  acount = pcounts[a]; 
        blist = pairs[b];  bcount = pcounts[b];
    }

    if (ht) {
        ht->root->tutte_type = CCelim_TUTTE_CD;
        ht->root->tutte_nodes[0] = a;
        ht->root->tutte_nodes[1] = b;
    }

    psys[0].next = &(psys[1]);
    psys[1].next = (CCelim_path *) NULL;
    psys[0].cnt = psys[1].cnt = 3;
    good = 1;
    for (i = 0; i < acount && good; i++) {
        if (alist[2*i] == b || alist[2*i+1] == b) continue;
        psys[0].node[0] = alist[2*i];
        psys[0].node[1] = a;
        psys[0].node[2] = alist[2*i+1];
        for (j = 0; j < bcount && good; j++) {
            if (blist[2*j] == a || blist[2*j+1] == a) continue;
            psys[1].node[0] = blist[2*j];
            psys[1].node[1] = b;
            psys[1].node[2] = blist[2*j+1];
            if (ht && ht->count < ht->max_count) {
                CC_MALLOC (ht_child, 1, CCelim_htnode);
                CCelim_init_htnode (ht_child);
                ht_child->id = ht->count++;
                ht_child->hamilton_type = CCelim_HAMILTON_PAIR;
                p = ht_child->hamilton_nodes;
                for (k = 0; k < 3; k++) p[k]   = psys[0].node[k];
                for (k = 0; k < 3; k++) p[k+3] = psys[1].node[k];
            } else {
                ht_child = (CCelim_htnode *) NULL;
            }
            rval = full_path_check (psys, &good, ab, D, G, ablist, level, F,
                  pcounts, pairs, timelimit, startzeit, use_tsp, depth,
                  ht, ht_child);
            CCcheck_rval (rval, "full_path_check failed");
            CCcheck_doctor

            if (ht_child) {
                if (ht_child->tutte_type != CCelim_TUTTE_NONE) {
                    ht_child->next = ht->root->child;
                    ht->root->child = ht_child;
                    if (depth > ht->depth) ht->depth = depth;
                } else {
                    CC_IFFREE (ht_child, CCelim_htnode);
                    ht->count--;
                }
            }
        }
    }
    *yesno = good; 

    if (ht && *yesno == 0) {
        CCelim_free_htnode_children (ht->root);
        ht->count = ht_count;
    }

CLEANUP:
    if (!pairs) {
        CC_IFFREE (alist, int);
        CC_IFFREE (blist, int);
    }
    return rval;
}

/* try_extra_point grows the elimination-search tree by selecting a node f  */
/* in the neighborhood ablist and adding all pairs of edges meeting f       */

static int try_extra_point (int *ab, CCelim_path *pin, CCelim_distobj *D,
        CCelim_graph *G, CCelim_array *ablist, int *good, int level,
        CCelim_graph *F, int *pcounts, int **pairs, double timelimit,
        double startzeit, int use_tsp, int depth, CCelim_httree *ht,
        CCelim_htnode *ht_parent)
{
    int rval = 0, i, j, f, fcount, *flist = (int *) NULL, okay;
    int pcount, plist[MAX_AB], a = pin->node[0], b = pin->node[1];
    CCelim_path psys[PMAX], pex;
    CCelim_htnode *ht_child = (CCelim_htnode *) NULL;
    int deg = 0, ht_count = 0;

    /* note: could set a and b to be the elements of CCelim_array ab */

    *good = 0;
    if (ht) ht_count = ht->count;

    CCelim_copy_path_system (pin, psys, 0);

    rval = potential_extra_points (ablist, psys, G, D, F, &pcount, plist,
                                   pcounts, pairs);
    CCcheck_rval (rval, "potential_extra_points failed");

    if (pcount > 1 && level > 0) {
        int winner = -1;
        rval = rank_extra_points (ab, psys, G, D, &pcount, plist, good,
           F, pcounts, pairs, timelimit, startzeit, depth, &winner);
        CCcheck_rval (rval, "rank_extra_points failed");
        if (*good == 1) {
            if (ht && ht_parent) {
                ht_parent->tutte_type = CCelim_TUTTE_POINT;
                ht_parent->tutte_nodes[0] = winner;
            }
            goto CLEANUP;
        }

    }
    CCcheck_doctor

    if (pcount == 0) goto CLEANUP;

    for (i = 0, deg = 0; i < pcount; i++) {
        if (G->nodelist[plist[i]].deg > deg) deg = G->nodelist[plist[i]].deg;
    }
    CC_MALLOC (flist, deg*deg, int);

    pex.next = psys;
    pex.cnt = 3;
    for (i = 0; i < pcount; i++) {
        f = plist[i];
        CCelim_check_neighbors_three_swap (a, b, f, G, D, &fcount, flist,
                                           F, pcounts, pairs);
        okay = 1;
        if (ht && ht_parent) {
            ht_parent->tutte_type = CCelim_TUTTE_POINT;
            ht_parent->tutte_nodes[0] = f;
        }
        for (j = 0; j < fcount && okay; j++) {
            pex.node[0] = flist[2*j];
            pex.node[1] = f;
            pex.node[2] = flist[2*j+1];

            if (ht && ht_parent && ht->count < ht->max_count) {
                CC_MALLOC (ht_child, 1, CCelim_htnode);
                CCelim_init_htnode (ht_child);
                ht_child->id = ht->count++;
                ht_child->hamilton_type = CCelim_HAMILTON_PATH;
                ht_child->hamilton_nodes[0] = pex.node[0];
                ht_child->hamilton_nodes[1] = pex.node[1];
                ht_child->hamilton_nodes[2] = pex.node[2];
            } else {
                ht_child = (CCelim_htnode *) NULL;
            }

            rval = full_path_check (&pex, &okay, ab, D, G, ablist, level-1, F,
                  pcounts, pairs, timelimit, startzeit, use_tsp, depth,
                  ht, ht_child);
            CCcheck_rval (rval, "full_path_check failed");
            CCcheck_doctor

            if (ht && ht_parent && ht_child) {
                if (okay == 1 && ht_child->tutte_type != CCelim_TUTTE_NONE) {
                    ht_child->next = ht_parent->child;
                    ht_parent->child = ht_child;
                    if (depth > ht->depth) ht->depth = depth;
                } else {
                    CC_IFFREE (ht_child, CCelim_htnode);
                    ht->count--;
                }
            }
        }
        if (ht && ht_parent && okay == 0) {
            CCelim_free_htnode_children (ht_parent);
            ht_parent->tutte_type = CCelim_TUTTE_NONE;
            ht->count = ht_count;
        }
        if (okay == 1) { *good = 1; goto CLEANUP; }
   }

CLEANUP:
    CC_IFFREE (flist, int);
    return rval;
}

/* try_extra_edge grows the elimination-search tree by selecting an end node */
/* c of the current path system pin and adding all edges meeting c           */

static int try_extra_edge (int *ab, CCelim_path *pin, CCelim_distobj *D,
        CCelim_graph *G, CCelim_array *ablist, int *good, int level,
        CCelim_graph *F, int *pcounts, int **pairs,
        double timelimit, double startzeit, int use_tsp, int depth,
        CCelim_httree *ht, CCelim_htnode *ht_parent)
{
    int rval = 0, i, j, k, cnt, okay, perm[NMAX], wlist[NMAX];
    int wrank[NMAX], nlist[NMAX], c, d, n, try, max_search = 1, max_branch = 1;
    int wcount = 0, bcount = 0, *blist = (int *) NULL, *middle = (int *) NULL;
    CCelim_path psys[PMAX], pex, *p;
    CCelim_htnode *ht_child = (CCelim_htnode *) NULL;
    int ht_count = 0;

    *good = 0;
    if (ht) ht_count = ht->count;

    if (depth >= DEEP_LEVEL) {
        max_search = MAX_DEEP_SEARCH; max_branch = MAX_DEEP_BRANCH;
    } else {
        max_search = MAX_EXTRA_SEARCH; max_branch = MAX_EXTRA_BRANCH;
    }

    CCelim_copy_path_system (pin, psys, 0);

    /* middle marks the interior of the paths */
    CC_MALLOC (middle, G->ncount, int);
    for (p = psys; p; p = p->next) {
        for (k = 0; k < 2; k++) {
            if (k == 0) c = p->node[0];
            else        c = p->node[p->cnt-1];
            for (j = 0; j < G->nodelist[c].deg; j++) {
                middle[G->nodelist[c].neighbors[j]] = 0;
            }
        }
    }
    for (p = psys; p; p = p->next) {
        for (j = 1; j < p->cnt-1; j++) middle[p->node[j]] = 1;
    }

    CC_MALLOC (blist, G->ncount, int);
    for (p = psys; p; p = p->next) {
        for (k = 0; k < 2; k++) {
            if (k == 0) { c = p->node[0];         n = p->node[1]; }
            else        { c = p->node[p->cnt-1];  n = p->node[p->cnt-2]; }
            CCelim_get_partner_pairs (G, c, n, pcounts, pairs, &bcount, blist);
            for (j = 0, cnt = 0; j < bcount; j++) {
                if (middle[blist[j]] == 0) cnt++;
            }
            if (cnt <= 2*max_branch) {  /* Bico: can try 2 here */
                wlist[wcount] = c; nlist[wcount] = n; wrank[wcount] = cnt;
                wcount++;
            }
        }
    }
    if (wcount == 0) goto CLEANUP;
    CCutil_int_perm_quicksort (perm, wrank, wcount);
        
    pex.next = psys;
    pex.cnt = 2;
    for (i = 0, try = 0; i < wcount && try < max_search; i++) {
        int bad, badcount;
        c = wlist[perm[i]];
        n = nlist[perm[i]];  /* n is adj to c in the path */
        CCelim_get_partner_pairs (G, c, n, pcounts, pairs, &bcount, blist);

        badcount = 0;
        for (j = 0; j < bcount && badcount <= max_branch; j++) {
            d = blist[j];
            if (middle[d] == 1) continue;
            pex.node[0] = c; pex.node[1] = d;
            bad = 0;
            rval = CCelim_test_paths (ab, &pex, D, G, 5 /*use_tsp*/, &bad);
            CCcheck_rval (rval, "CCelim_test_paths failed");
            badcount += bad;
            CCcheck_doctor
        }
        if (badcount > max_branch) continue;

        if (ht && ht_parent) {
            ht_parent->tutte_type = CCelim_TUTTE_END;
            ht_parent->tutte_nodes[0] = c;
        }

        okay = 1;
        for (j = 0; j < bcount && okay; j++) {
            d = blist[j];
            if (middle[d] == 1) continue;
            pex.node[0] = c; pex.node[1] = d;

            if (ht && ht_parent && ht->count < ht->max_count) {
                CC_MALLOC (ht_child, 1, CCelim_htnode);
                CCelim_init_htnode (ht_child);
                ht_child->id = ht->count++;
                ht_child->hamilton_type = CCelim_HAMILTON_EDGE;
                ht_child->hamilton_nodes[0] = pex.node[0];
                ht_child->hamilton_nodes[1] = pex.node[1];
            } else {
                ht_child = (CCelim_htnode *) NULL;
            }

            rval = full_path_check (&pex, &okay, ab, D, G, ablist, level-1, F,
                  pcounts, pairs, timelimit, startzeit, use_tsp, depth,
                  ht, ht_child);
            CCcheck_rval (rval, "full_path_check failed");
            CCcheck_doctor

            if (ht && ht_parent && ht_child) {
                if (okay == 1 && ht_child->tutte_type != CCelim_TUTTE_NONE) {
                    ht_child->next = ht_parent->child;
                    ht_parent->child = ht_child;
                    if (depth > ht->depth) ht->depth = depth;
                } else {
                    CC_IFFREE (ht_child, CCelim_htnode);
                    ht->count--;
                }
            }
        }
        if (ht && ht_parent && okay == 0) {
            CCelim_free_htnode_children (ht_parent);
            ht->count = ht_count;
            ht_parent->tutte_type = CCelim_TUTTE_NONE; 
        }
        if (okay == 1) { *good = 1; goto CLEANUP; }
        try++;
    }

CLEANUP:
    CC_IFFREE (blist, int);
    CC_IFFREE (middle, int);
    return rval;
}

/* full_path_check sets yesno to 1 if the path system psys is not in any  */
/* optimal tour. if level is nonegative, then elimation tree is grown     */
/* further in an attempt to prove that psys is not in any optimal tour    */

static int full_path_check (CCelim_path *psys, int *yesno, int *ab,
        CCelim_distobj *D, CCelim_graph *G, CCelim_array *ablist, int level,
        CCelim_graph *F, int *pcounts, int **pairs, double timelimit,
        double startzeit, int use_tsp, int depth, CCelim_httree *ht,
        CCelim_htnode *ht_parent)
{
    int rval = 0, bad = 0, valid, good, pathcnt = 0;
    CCelim_path pmerge[PMAX];

    *yesno = 0;

    if (ht && depth >= ht->max_depth) ht = (CCelim_httree *) NULL;

    rval = CCelim_validate_paths (G->ncount, psys, &valid, pmerge);
    CCcheck_rval (rval, "CCelim_validate_paths failed");
    if (!valid) { *yesno = 1; goto CLEANUP; }

    pathcnt = CCelim_path_system_count (pmerge);
    rval = CCelim_test_paths (ab, pmerge, D, G, use_tsp, &bad);
    CCcheck_rval (rval, "CCelim_test_paths failed");
    if (bad && level >= 0 && pathcnt < MAX_PATH_COUNT) {
        rval = try_extra_point (ab, pmerge, D, G, ablist, &good, level, F,
              pcounts, pairs, timelimit, startzeit, use_tsp, depth+2, ht,
              ht_parent);
        CCcheck_rval (rval, "try_extra_point failed");
        CCcheck_doctor
        if (good) bad = 0;
    }
    if (bad && level >= 0) {
        rval = try_extra_edge (ab, pmerge, D, G, ablist, &good, level, F,
              pcounts, pairs, timelimit, startzeit, use_tsp, depth+1,
              ht, ht_parent);
        CCcheck_rval (rval, "try_extra_edge failed");
        CCcheck_doctor
        if (good) bad = 0;
    }
    *yesno = 1 - bad;

CLEANUP:
    return rval;
}

/* potential_cd builds a lists of pairs of nodes (c,d) such that edge cd is  */
/* incompatible with ab and the nodes c and d are candidates to be the first */
/* two demand nodes in the elimination tree                                  */

static int potential_cd (int a, int b, CCelim_array *ablist, CCelim_graph *G,
        CCelim_distobj *D, int *pcount, int **plist, int wtype, CCelim_graph *F,
        int *pcounts, int **pairs)
{
    int rval = 0, i, j, c, d, in, co, count = 0, count0, count1;
    int *wrank = (int *) NULL, *wlist = (int *) NULL;
    int *mark = (int *) NULL, *perm = (int *) NULL;
    CCelim_node *pa = &G->nodelist[a], *pb = &G->nodelist[b];

    *plist = (int *) NULL;
    *pcount = 0;

    if (ablist->count == 0) goto CLEANUP;

    CC_MALLOC (wlist, 2 * ablist->count * ablist->count, int);
    CC_MALLOC (wrank, ablist->count * ablist->count, int);

    if (wtype == CCelim_CD_TWO) {
        CC_MALLOC (mark, G->ncount, int);
        for (i = 0; i < ablist->count; i++) mark[ablist->arr[i]] = 0;
        for (i = 0; i < pb->deg; i++) mark[pb->neighbors[i]] = 0;
        for (i = 0; i < pa->deg; i++) mark[pa->neighbors[i]] = 1;
        for (i = 0; i < pb->deg; i++) {
            if (mark[pb->neighbors[i]] == 1) mark[pb->neighbors[i]] = 2;
        }
    }

    for (i = 0; i < ablist->count; i++) {
        c = ablist->arr[i];
        if (G->nodelist[c].deg > CCelim_MAX_DEGREE) continue;
        CCelim_check_neighbors_three_swap(a, b, c, G, D, &count0, (int *) NULL,
                                          F, pcounts, pairs);
        for (j = i+1; j < ablist->count; j++) {
            d = ablist->arr[j];
            if (G->nodelist[d].deg > CCelim_MAX_DEGREE) continue;
            in = CCelim_edge_in_graph(c,d,G);
            co = compatible_test (a,b,c,d,D);
            if ((wtype == CCelim_CD_NONEDGE && (!in || !co)) ||
                (wtype == CCelim_CD_EDGE && (in && !co)) ||
                (wtype == CCelim_CD_DIST2 && ((in || ed2(c,d,G)) && !co)) ||
                (wtype == CCelim_CD_TWO && mark[c]==2 && mark[d]==2 && !co)) {
                CCelim_check_neighbors_three_swap(a, b, d, G, D, &count1,
                     (int *) NULL, F, pcounts, pairs);
                wlist[2*count]   = c;
                wlist[2*count+1] = d;
                wrank[count] = count0*count1;
                count++;
            }
        }
    }

    if (count) {
        CC_MALLOC (perm, count, int);
        CCutil_int_perm_quicksort (perm, wrank, count);

        *plist = CC_SAFE_MALLOC (2*count, int);
        CCcheck_NULL (*plist, "out of memory plist");
        for (i = 0; i < count; i++) {
            (*plist)[2*i]   = wlist[2*perm[i]];
            (*plist)[2*i+1] = wlist[2*perm[i]+1];
        }
        *pcount = count;
    } 

CLEANUP:
    CC_IFFREE (mark, int);
    CC_IFFREE (wlist, int);
    CC_IFFREE (wrank, int);
    CC_IFFREE (perm, int);
    return rval;
}

/* ed2 checks if nodes a and b have a common neighbor */

static int ed2 (int a, int b, CCelim_graph *G)
{
    int i, j;
    CCelim_node *n;

    for (i = 0; i < G->nodelist[a].deg; i++) {
        n = &G->nodelist[G->nodelist[a].neighbors[i]];
        for (j = 0; j < n->deg; j++) {
            if (n->neighbors[j] == b) return 1;
        }
    }
    return 0; 
}

/* potential_extra_points builds a list of candidates to use in the func */
/* try_extra_point when growing the elimination tree                     */

static int potential_extra_points (CCelim_array *ablist, CCelim_path *psys,
        CCelim_graph *G, CCelim_distobj *D, CCelim_graph *F, int *pcount,
        int *plist, int *pcounts, int **pairs)
{
    int rval = 0, i, f, fcount, wlist[MAX_AB], wrank[MAX_AB], perm[MAX_AB];
    int count = 0, a = psys->node[0], b = psys->node[1], *marks = (int *) NULL;

    *pcount = 0;

    /* mark the nodes in the path system; only consider unmarked nodes */

    CC_MALLOC (marks, G->ncount, int);
    for (i = 0; i < ablist->count; i++) marks[ablist->arr[i]] = 0;
    CCelim_path_system_mark (psys, marks, 1);

    for (i = 0; i < ablist->count; i++) {
        f = ablist->arr[i];
        if (marks[f] == 0 && G->nodelist[f].deg <= CCelim_MAX_DEGREE) {
            CCelim_check_neighbors_three_swap(a, b, f, G, D, &fcount,
                                              (int *) NULL, F, pcounts, pairs);
            wrank[count] = fcount;
            wlist[count++] = f;
        }
    }

    if (count) {
        CCutil_int_perm_quicksort (perm, wrank, count);
        for (i = 0; i < count; i++) plist[i] = wlist[perm[i]];
        *pcount = count;
    }


CLEANUP:
    CC_IFFREE (marks, int);
    return rval;
}

/* rank_extra_points orders the candidates for try_extra_point by the    */
/* fewest edge-pairs that cannot be immediately shown to be incompatible */
/* with the path system pin                                              */

static int rank_extra_points (int *ab, CCelim_path *pin, CCelim_graph *G,
        CCelim_distobj *D, int *pcount, int *plist, int *yesno,
        CCelim_graph *F, int *pcounts, int **pairs, double timelimit,
        double startzeit, int depth, int *winner)
{
    int rval = 0, i, j, bad, badcount, f, fcount, *flist = (int *) NULL;
    int wlist[MAX_AB], wrank[MAX_AB], perm[MAX_AB], newcount = 0, deg = 0;
    int a = pin->node[0], b = pin->node[1], max_search, max_branch;
    CCelim_path *p, psys[PMAX], pex;

    *yesno = 0;

    if (depth >= DEEP_LEVEL) {
        max_search = MAX_DEEP_SEARCH; max_branch = MAX_DEEP_BRANCH;
    } else {
        max_search = MAX_EXTRA_SEARCH; max_branch = MAX_EXTRA_BRANCH;
    }

    for (i = 0; i < *pcount; i++) {
        if (G->nodelist[plist[i]].deg > deg) deg = G->nodelist[plist[i]].deg;
    }
    CC_MALLOC (flist, deg*deg, int);
    CCelim_copy_path_system (pin, psys, 0);
    for (p = psys; p->next; p = p->next);
    p->next = &pex;
    pex.next = (CCelim_path *) NULL;
    pex.cnt = 3;

    for (i = 0; i < *pcount; i++) {
        f = plist[i];
        CCelim_check_neighbors_three_swap (a, b, f, G, D, &fcount, flist, F,
                                           pcounts, pairs);
        badcount = 0;
        for (j = 0; j < fcount && badcount <= max_branch; j++) {
            pex.node[0] = flist[2*j];
            pex.node[1] = f;
            pex.node[2] = flist[2*j+1];
            bad = 0;
            /* Bico: also possible to run CCelim_test_paths with 1, not 5 */
            rval = CCelim_test_paths (ab, psys, D, G, 5, &bad);
            CCcheck_rval (rval, "CCelim_test__paths failed");
            badcount += bad;
            CCcheck_doctor
        }
        if (badcount == 0) { *yesno = 1; *winner = f; goto CLEANUP; }
        if (badcount <= max_branch) {
            wrank[newcount]   = badcount;
            wlist[newcount++] = f;
        }
    }

    if (newcount) {
        CCutil_int_perm_quicksort (perm, wrank, newcount);
        for (i = 0; i < newcount; i++) plist[i] = wlist[perm[i]];
        *pcount = (newcount <= max_search ? newcount : max_search);
    } else {
        *pcount = 0;
    }

CLEANUP:
    CC_IFFREE (flist, int);
    return rval;
}

/* check_neighbors_three_swap looks for a quick proof that pairs of edges */
/* meeting n are incompatible with ab                                     */

void CCelim_check_neighbors_three_swap (int a, int b, int n, CCelim_graph *G, 
        CCelim_distobj *D, int *tcount, int *tlist, CCelim_graph *F,
        int *pcounts, int **pairs)
{
    /* For each pair of neighbors of n, check for 3-opt move with a-b   */
    /*                                                                  */
    /*    n0      n1                                                    */
    /*     \     /                                                      */
    /*      \   /                                                       */
    /*        n       Consider replacing a-b n-n0 n-1 by a-n n-b n0-n1  */
    /*                                                                  */
    /*    a ------ b                                                    */
    /*                                                                  */

    int i, j, n0, n1, count = 0;
    CCelim_node *pn = &G->nodelist[n], *fn = (CCelim_node *) NULL;

    if (F) fn = &F->nodelist[n];

    if (fn && fn->deg == 2) {   /* have a set of fixed edges */
        if (tlist) { tlist[0] = fn->neighbors[0]; tlist[1] = fn->neighbors[1]; }
        count = 1;
    } else if (fn && fn->deg == 1) {
        n0 = fn->neighbors[0];
        for (j = 0; j < pn->deg; j++) {
            if ((n1 = pn->neighbors[j]) == n0) continue;
            if (compatible_test(a,b,n,n1,D) && test_three(a,b,n,n0,n1,D)) {
                if (tlist) { tlist[2*count] = n0; tlist[2*count+1] = n1; }
                count++;
            }
        }
    } else {
        if (pairs) {
            for (i = 0; i < pcounts[n]; i++) {
                n0 = pairs[n][2*i];
                n1 = pairs[n][2*i+1];
                if (CCelim_edge_in_graph (n, n0, G) &&  
                    CCelim_edge_in_graph (n, n1, G) &&
                    compatible_test(a,b,n,n0,D) &&
                    compatible_test(a,b,n,n1,D) && test_three(a,b,n,n0,n1,D)) {
                    if (tlist) { tlist[2*count] = n0; tlist[2*count+1] = n1; }
                    count++;
                }
            }
        } else {
            int compat[CCelim_MAX_DEGREE], *ndist = pn->len;
            int ab_dist = dist(a, b, D);
            int an_dist = dist(a, n, D);
            int bn_dist = dist(b, n, D);

            for (i = 0; i < pn->deg; i++) {
                n0 = pn->neighbors[i];
                compat[i] = (an_dist+dist (b,n0,D) >= ab_dist+ndist[i] ||
                             bn_dist+dist (a,n0,D) >= ab_dist+ndist[i]);
            }

            for (i = 0; i < pn->deg; i++) {
                if (!compat[i]) continue;
                n0 = pn->neighbors[i];
                for (j = i+1; j < pn->deg; j++) {
                    if (!compat[j]) continue;
                    n1 = pn->neighbors[j];
                    if ((n0 == a && n1 == b) || (n0 == b && n1 == a)) continue;
                    if (ab_dist + ndist[i] + ndist[j] <=
                        an_dist + bn_dist + dist(n0,n1,D)) {
                        if (tlist) { 
                            tlist[2*count] = n0;
                            tlist[2*count+1] = n1;
                        }
                        count++;
                    }
                }
            }
        }
    }
    if (tcount) *tcount = count;
}

/* pair_neighbors_three_swap looks for a quick proof that parirs of edges */
/* meeting n are incompatible with the path a-b-c                         */

void CCelim_pair_neighbors_three_swap (int a, int b, int c, int n,
        CCelim_graph *G, CCelim_distobj *D, int *tcount, int *tlist,
        int *pcounts, int **pairs)
{
    int i, j, n0, n1, count = 0;
    CCelim_node *pn = &G->nodelist[n];

    /* a version of check_neighbors_three_swap for a path a-b-c; in this */
    /* case we check for three swaps with n and both ab and ac.          */

    if (pairs) {
        for (i = 0; i < pcounts[n]; i++) {
            n0 = pairs[n][2*i];
            n1 = pairs[n][2*i+1];
            if (n0 == b || n1 == b || (n0 == a && n1 == c) ||
                                      (n0 == c && n1 == a)) continue;
            if (CCelim_edge_in_graph (n, n0, G) &&  
                CCelim_edge_in_graph (n, n1, G) &&
                compatible_test(a,b,n,n0,D) &&
                compatible_test(a,b,n,n1,D) &&
                compatible_test(b,c,n,n0,D) &&
                compatible_test(b,c,n,n1,D) &&
                test_three(a,b,n,n0,n1,D)   &&
                test_three(b,c,n,n0,n1,D)) {
                if (tlist) { tlist[2*count] = n0; tlist[2*count+1] = n1; }
                count++;
            }
        }
    } else {
        int compat[CCelim_MAX_DEGREE], *ndist = pn->len;
        int ab_dist = dist(a, b, D);
        int bc_dist = dist(b, c, D);
        int an_dist = dist(a, n, D);
        int bn_dist = dist(b, n, D);
        int cn_dist = dist(c, n, D);

        for (i = 0; i < pn->deg; i++) {
            n0 = pn->neighbors[i];
            compat[i] = ((n0 != b) &&
                         (an_dist+dist(b, n0, D) >= ab_dist+ndist[i] ||
                          bn_dist+dist(a, n0, D) >= ab_dist+ndist[i]) && 
                         (bn_dist+dist(c, n0, D) >= bc_dist+ndist[i] ||
                          cn_dist+dist(b, n0, D) >= bc_dist+ndist[i]));
        }

        for (i = 0; i < pn->deg; i++) {
            if (!compat[i]) continue;
            n0 = pn->neighbors[i];
            for (j = i+1; j < pn->deg; j++) {
                if (!compat[j]) continue;
                n1 = pn->neighbors[j];
                if ((n0 == a && n1 == c) || (n0 == c && n1 == a)) continue;
                if (ab_dist + ndist[i] + ndist[j] <=
                    an_dist + bn_dist + dist(n0,n1,D) &&
                    bc_dist + ndist[i] + ndist[j] <=
                    bn_dist + cn_dist + dist(n0,n1,D)) {
                    if (tlist) { 
                        tlist[2*count] = n0;
                        tlist[2*count+1] = n1;
                    }
                    count++;
                }
            }
        }
    }
    if (tcount) *tcount = count;
}

/* compatible_test checks if edges ab and cd cannot be shown to be        */
/* incompatible by a simple 2-swap (that is, one of the two orientations  */
/* is not an improving 2-swap)                                            */

static int compatible_test (int a, int b, int c, int d, CCelim_distobj *D)
{
    int t = dist(a,b,D) + dist(c,d,D);
    return (dist(a,d,D) + dist(b,c,D) >= t || dist(a,c,D) + dist(b,d,D) >= t) ;
}

/* test_three checks that a simple 3-swap is not improving */

static int test_three (int a, int b, int n, int n0, int n1, CCelim_distobj *D) 
{
   return ((n0 != a || n1 != b) && (n0 != b || n1 != a) &&
       dist(a,b,D)   + dist(n,n0,D) + dist(n,n1,D) <=
       dist(n0,n1,D) + dist(a,n,D)  + dist(b,n,D));
}

/* all_neighbor_pairs builds a list of all of the pairs of edges meeting n */
/* that include any possible fixed edges                                   */

void CCelim_all_neighbor_pairs (int n, CCelim_graph *G, CCelim_graph *F,
        int *pcount, int *plist)
{
    int i, j, n0, n1, count = 0;
    CCelim_node *pn = &G->nodelist[n], *fn = (CCelim_node *) NULL;

    if (F) fn = &F->nodelist[n];

    if (fn && fn->deg == 2) {   /* have a set of fixed edges */
        if (plist) { plist[0] = fn->neighbors[0]; plist[1] = fn->neighbors[1]; }
        count = 1;
    } else if (fn && fn->deg == 1) {
        n0 = fn->neighbors[0];
        for (j = 0; j < pn->deg; j++) {
            if ((n1 = pn->neighbors[j]) == n0) continue;
            if (plist) { plist[2*count] = n0; plist[2*count+1] = n1; }
            count++;
        }
    } else {
        for (i = 0; i < pn->deg; i++) {
            n0 = pn->neighbors[i];
            for (j = i+1; j < pn->deg; j++) {
                n1 = pn->neighbors[j];
                if (plist) { plist[2*count] = n0; plist[2*count+1] = n1; }
                count++;
            }
       }
   }
   *pcount = count;
}

/* CCelim_build_ablist builds a list of nodes in a neighborhood of edge ab */

int CCelim_build_ablist (int a, int b, CCelim_graph *G, CCelim_distobj *D,
        CCelim_array *ablist, int max_neighborhood, CCkdtree *kt,
        CCrandstate *rstate, CCdatagroup *euclid_dat, CCkdtree *euclid_kt)  
{
    int rval = 0, i, nwant, *len = (int *) NULL, *perm = (int *) NULL, norm;
    double tx, mx, ty, my, tz = 0.0, mz; 

    /* for Euclidean data, neighborhood is from center of [a,b] line    */
    /* for other data, use a breadth-first-search neighborhood around a */

    nwant = (max_neighborhood <= G->ncount-1 ? max_neighborhood : G->ncount-1);
    ablist->count = 0;
    CC_MALLOC (ablist->arr, nwant, int);
    CC_MALLOC (perm, nwant, int);
    CC_MALLOC (len, nwant, int);

    CCutil_dat_getnorm (D->dat, &norm);
    if (kt && (norm == CC_EUCLIDEAN || norm == CC_EUCLIDEAN_CEIL ||
               norm == CC_EUCLIDEAN_3D)) {
        tx = D->dat->x[a], mx = (D->dat->x[a] + D->dat->x[b]) / 2.0;
        ty = D->dat->y[a], my = (D->dat->y[a] + D->dat->y[b]) / 2.0;
        D->dat->x[a] = mx;
        D->dat->y[a] = my;
        if (D->dat->z) {
            tz = D->dat->z[a], mz = (D->dat->z[a] + D->dat->z[b]) / 2.0;
            D->dat->z[a] = mz;
        }
        CCkdtree_delete (kt, b);
        rval = CCkdtree_node_k_nearest (kt, G->ncount, a, nwant, D->dat,
                                        NULL, ablist->arr, rstate);
        CCkdtree_undelete (kt, b);

        /* sort from closest to furthest */
        for (i = 0; i < nwant; i++) {
            len[i] = CCutil_dat_edgelen (a, ablist->arr[i], D->dat);
        }
        CCutil_int_perm_quicksort (perm, len, nwant);
        for (i = 0; i < nwant; i++) len[i] = ablist->arr[perm[i]];
        for (i = 0; i < nwant; i++) ablist->arr[i] = len[i];
        
        D->dat->x[a] = tx;
        D->dat->y[a] = ty;
        if (D->dat->z) {
            D->dat->z[a] = tz;
        }
        CCcheck_rval (rval, "CCkdtree_node_k_nearest failed");
        ablist->count = nwant;
    } else if (euclid_dat && euclid_kt) {
        tx = euclid_dat->x[a], mx = (euclid_dat->x[a] + euclid_dat->x[b]) / 2.0;
        ty = euclid_dat->y[a], my = (euclid_dat->y[a] + euclid_dat->y[b]) / 2.0;
        euclid_dat->x[a] = mx;
        euclid_dat->y[a] = my;
        CCkdtree_delete (euclid_kt, b);
        rval = CCkdtree_node_k_nearest (euclid_kt, G->ncount, a, nwant, 
                                        euclid_dat, NULL, ablist->arr, rstate);
        CCkdtree_undelete (euclid_kt, b);

        for (i = 0; i < nwant; i++) {
            len[i] = CCutil_dat_edgelen (a, ablist->arr[i], D->dat);
        }
        CCutil_int_perm_quicksort (perm, len, nwant);
        for (i = 0; i < nwant; i++) len[i] = ablist->arr[perm[i]];
        for (i = 0; i < nwant; i++) ablist->arr[i] = len[i];
        
        euclid_dat->x[a] = tx;
        euclid_dat->y[a] = ty;
        CCcheck_rval (rval, "CCkdtree_node_k_nearest failed");
        ablist->count = nwant;
    } else {
        int *mark = (int *) NULL, *Q = (int *) NULL;
        int qhead, qtail, c, n;
        CC_MALLOC (mark, G->ncount, int);
        for (i = 0; i < G->ncount; i++) mark[i] = 0;
        CC_MALLOC (Q, G->ncount, int);
        Q[0] = a;  Q[1] = b;
        qhead = 0; qtail = 2;
        mark[a] = 1; mark[b] = 1;
        while (qhead < qtail && ablist->count < nwant) {
            n = Q[qhead++];
            for (i = 0; i < G->nodelist[n].deg && ablist->count < nwant; i++) {
                c = G->nodelist[n].neighbors[i];
                if (mark[c] == 0) {
                    ablist->arr[ablist->count++] = c;
                    mark[c] = mark[n]+1;
                    /* don't want too far */
                    if (mark[c] < 10) Q[qtail++] = c;
                }
            }
        }
        CC_IFFREE (Q, int);
        CC_IFFREE (mark, int);
    }

CLEANUP:
    CC_IFFREE (len, int);
    CC_IFFREE (perm, int);
    return rval;
}

/* CCelim_alloc_elim allocates memory for data structures used in the */
/* elimination routines                                               */

int CCelim_alloc_elim (int ncount, int ecount, int *elist, int fixcount, 
        int *fixlist, CCdatagroup *dat, CCelim_graph *G, CCelim_graph *F,
        CCelim_distobj *D, CCkdtree **kt, CCdatagroup **euclid_dat,
        CCkdtree **euclid_kt, CCrandstate *rstate)
{
    int i, norm, rval = 0;

    if (G) CCelim_init_graph (G);
    if (F) CCelim_init_graph (F);
    CCelim_init_distobj (D);

    rval = CCelim_build_distobj (D, ncount, dat);
    CCcheck_rval (rval, "CCelim_build_distobj failed");

    if (G) {
        rval = CCelim_build_graph (G, ncount, ecount, elist, (int *) NULL, D);
        CCcheck_rval (rval, "CCelim_build_graph failed");
    }

    if (F) {
        rval = CCelim_build_graph (F, ncount, fixcount, fixlist, (int *) NULL,
                                   D);
        CCcheck_rval (rval, "CCelim_build_graph failed");
    }

    CCutil_dat_getnorm (dat, &norm);
    if (kt && (norm & CC_NORM_BITS) == CC_KD_NORM_TYPE) {
        CC_MALLOC (*kt, 1, CCkdtree);
        rval = CCkdtree_build (*kt, ncount, dat, (double *) NULL, rstate);
        CCcheck_rval (rval, "CCkdtree_build failed");
    }

    if (euclid_dat && norm == CC_ROAD) {
        /* build Euclidean kd-tree for use in elim/build_ablist */
        CC_MALLOC (*euclid_dat, 1, CCdatagroup);
        CCutil_init_datagroup (*euclid_dat);

        CCutil_dat_setnorm (*euclid_dat, CC_EUCLIDEAN);
        CC_MALLOC ((*euclid_dat)->x, ncount, double);
        CC_MALLOC ((*euclid_dat)->y, ncount, double);
        for (i = 0; i < ncount; i++) {
            (*euclid_dat)->x[i] = 1000.0 * dat->x[i];
            (*euclid_dat)->y[i] = 1000.0 * dat->y[i];
        }
        if (euclid_kt) {
            CC_MALLOC (*euclid_kt, 1, CCkdtree);
            rval = CCkdtree_build (*euclid_kt, ncount, *euclid_dat,
                                  (double *) NULL, rstate);
            CCcheck_rval (rval, "CCkdtree_build failed");
        }
    }

CLEANUP:
    return rval;
}

/* CCelim_free_elim frees the memory allocated by CCelim_alloc_elim */

void CCelim_free_elim (CCelim_graph *G, CCelim_graph *F, CCelim_distobj *D,
        CCkdtree **kt, CCdatagroup **euclid_dat, CCkdtree **euclid_kt)
{
    CCelim_free_distobj (D);
    CCelim_free_graph (G);
    if (F) CCelim_free_graph (F);
    if (kt && *kt) {
        CCkdtree_free (*kt);
        CC_FREE (*kt, CCkdtree);
    }
    if (euclid_dat && *euclid_dat) {
        CCutil_freedatagroup (*euclid_dat);
        CC_FREE (*euclid_kt, CCkdtree);
    }
    if (euclid_kt && *euclid_kt) {
        CCkdtree_free (*euclid_kt);
        CC_FREE (*euclid_kt, CCkdtree);
    }
}

/* get_partner_pairs builds a list of nodes b such that an and ab are in */
/* the pairs list specified by pcounts and pairs                         */

void CCelim_get_partner_pairs (CCelim_graph *G, int a, int n, int *pcounts,
        int **pairs, int *bcount, int *blist)
{
    int i, b, cnt = 0;

    if (pcounts && pairs) {
        for (i = 0; i < pcounts[a]; i++) {
            b = -1;
            if      (pairs[a][2*i]   == n) b = pairs[a][2*i+1];
            else if (pairs[a][2*i+1] == n) b = pairs[a][2*i];
            if (b != -1) {
                if (CCelim_edge_in_graph (a, b, G)) {
                    blist[cnt++] = b;
                }
            }
        }
    } else {
        for (i = 0; i < G->nodelist[a].deg; i++) {
            if (G->nodelist[a].neighbors[i] != n) {
                blist[cnt++] = G->nodelist[a].neighbors[i];
            }
        }
    }
    *bcount = cnt;
}

void CCelim_init_distobj (CCelim_distobj *D)
{
    D->dat = (CCdatagroup *) NULL;
    D->cacheind  = (int *) NULL;
    D->cacheval  = (int *) NULL;
}

void CCelim_free_distobj (CCelim_distobj *D)
{
    if (D) {
         D->dat = (CCdatagroup *) NULL;
         CC_IFFREE (D->cacheind, int);
         CC_IFFREE (D->cacheval, int);
    }
}

int CCelim_build_distobj (CCelim_distobj *D, int ncount, CCdatagroup *dat)
{
    int rval = 0;
    int i, cacheM;

    CCelim_init_distobj (D);
    D->dat = dat;

    i = 0;
    while ((1 << i) < ncount)
        i++;
    cacheM = (1 << i);

    CC_MALLOC (D->cacheind, cacheM, int);
    CC_MALLOC (D->cacheval, cacheM, int);
    for (i = 0; i < cacheM; i++) {
        D->cacheind[i] = -1;
    }

CLEANUP:
    if (rval) CCelim_free_distobj (D);
    return rval; 
}

/* NOTE: Code in CCelim_dist is duplicated in macro CC_ELIM_DIST */

int CCelim_dist (int i, int j, CCelim_distobj *D)  /* Bentley's kdtree paper */
{
    int ind;

    if (i > j) { CC_SWAP (i, j, ind); }

    ind = i ^ j;

    if (D->cacheind[ind] != i) {
        D->cacheind[ind] = i;
        D->cacheval[ind] = CCutil_dat_edgelen (i, j, D->dat);
    }
    return D->cacheval[ind];
}
