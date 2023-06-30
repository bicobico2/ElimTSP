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
/*        Edge Elimination: Verify results using Hamilton-Tutte trees       */
/*                                                                          */
/*  The verification will use the H-T trees to choose the Tutte moves in    */
/*  the elimination search.                                                 */
/*                                                                          */
/*  Written by:  Cook, 2018                                                 */
/*                                                                          */
/****************************************************************************/

#include <stdio.h>
#include "ccutil.h"
#include "elim.h"

static char *edgefname = (char *) NULL;
static char *tspfname = (char *) NULL;
static char *fixedfname = (char *) NULL;
static char *nonpairsfname = (char *) NULL;
static char *verifyfname = (char *) NULL;
static char *outedgefname = (char *) NULL;
static char *outfixfname = (char *) NULL;
static char *outpairfname = (char *) NULL;
static int beverbose = 0;
static int binary_in = 0;
static int show_trace_miss = 0;
static int seed = 0;
static double elim_time_limit = 10000000.0;

int main (int ac, char **av);
static int run_verify (char *fname, int ncount, int ecount, int *elist,
    int fixcount, int *fixlist, int *nonpcounts, int **nonpairs, int loud,
    CCdatagroup *dat, CCrandstate *rstate);
static int verify_edge (int a, int b, CCelim_graph *G, CCelim_graph *F,
    CCelim_distobj *D, int *yesno, double *ezeit, int *pcounts, int **pairs,
    int longpath, int use_tsp, CCelim_httree *vt,
    int *htcounter, int *hamcounter);
static int verify_fixed_edge (int a, int b, CCelim_graph *G, CCelim_graph *F,
    CCelim_distobj *D, int *yesno, double *ezeit, int *pcounts, int **pairs,
    int use_tsp, CCelim_httree *vt, int *htcounter, int *hamcounter);
static int verify_pair (int a, int b, int c, CCelim_graph *G, CCelim_graph *F,
    CCelim_distobj *D, int *yesno, double *ezeit, int *pcounts, int **pairs,
    int longpath, int use_tsp, CCelim_httree *vt,
    int *htcounter, int *hamcounter);
static int verify_with_cd (int a, int b, int c, int d,
    CCelim_distobj *D, CCelim_graph *G, int *yesno,
    CCelim_graph *F, int *pcounts, int **pairs, double timelimit,
    int use_tsp, CCelim_htnode *vt, int *htcounter, int *hamcounter);
static int verify_via_long_path (int path_len, int a, int b,
    CCelim_distobj *D, CCelim_graph *G, int *yesno,
    CCelim_graph *F, int *pcounts, int **pairs, double timelimit,
    int use_tsp, CCelim_htnode *vt, int *htcounter, int *hamcounter);
static int verify_pair_via_long_path (int path_len, int a, int b, int c,
    CCelim_distobj *D, CCelim_graph *G, int *yesno,
    CCelim_graph *F, int *pcounts, int **pairs, double timelimit, int use_tsp,
    CCelim_htnode *vt, int *htcounter, int *hamcounter);
static int verify_grow_long_path (int path_len, int *ab, CCelim_distobj *D,
    CCelim_graph *G, int *yesno, CCelim_graph *F,
    int *pcounts, int **pairs, CCelim_path *psys, double timelimit,
    int use_tsp, CCelim_htnode *vt, int *htcounter, int *hamcounter);
static int verify_full_path_check (CCelim_path *psys, int *yesno, int *ab,
    CCelim_distobj *D, CCelim_graph *G, CCelim_graph *F,
    int *pcounts, int **pairs, double timelimit, int use_tsp, int depth,
    CCelim_htnode *vt, int *htcounter, int *hamcounter);
static int verify_try_extra_point (int *ab, CCelim_path *pin, CCelim_distobj *D,
    CCelim_graph *G, int *good, CCelim_graph *F,
    int *pcounts, int **pairs, double timelimit, int use_tsp, int depth,
    CCelim_htnode *vt, int *htcounter, int *hamcounter);
static int verify_try_extra_edge (int *ab, CCelim_path *pin, CCelim_distobj *D,
    CCelim_graph *G, int *good, CCelim_graph *F,
    int *pcounts, int **pairs, double timelimit, int use_tsp, int depth,
    CCelim_htnode *vt, int *htcounter, int *hamcounter);
static void free_pair_info (int ncount, int **pcounts, int ***pairs);
static int compatible_test (int a, int b, int c, int d, CCelim_distobj *D);
static int get_file_counts (char *fname, int binary_in, int *count_edge,
    int *count_fix, int *count_pair);
static int parseargs (int ac, char **av);
static void usage (char *f);

int main (int ac, char **av)
{
    int rval = 0, i;
    int ncount = 0, ecount = 0, *elist = (int *) NULL, *elen = (int *) NULL;
    int fixcount = 0, *fixlist = (int *) NULL, *fixlen = (int *) NULL;
    int *nonpaircounts = (int *) NULL, **nonpairs = (int **) NULL;
    CCdatagroup dat;
    CCrandstate rstate;

    CCutil_init_datagroup (&dat);

    rval = parseargs (ac, av);
    if (rval) goto CLEANUP;

    rval = CCutil_print_command (ac, av);
    CCcheck_rval (rval, "CCutil_print_command failed");
    CCutil_printlabel ();
    CCutil_sprand (seed, &rstate);

    if (!tspfname) {
        fprintf (stderr, "must specify TSPLIB file with -T\n");
        rval = 1; goto CLEANUP;
    }

    if (!verifyfname) {
        fprintf (stderr, "must specify Hamilton-Tutte file with -V\n");
        rval = 1; goto CLEANUP;
    }

    rval = CCutil_gettsplib (tspfname, &ncount, &dat);
    CCcheck_rval (rval, "CCutil_gettsplib failed");

    rval = CCutil_getedgelist (ncount, edgefname, &ecount, &elist, &elen, 0);
    CCcheck_rval (rval, "CCutil_getedgelist failed");
    printf ("Edges:      %d\n", ecount); fflush (stdout);
    for (i = 0; i < ecount; i++) {
        if (elen[i] != CCutil_dat_edgelen (elist[2*i], elist[2*i+1], &dat)) {
            fprintf (stderr, "edge lengths do not match TSPLIB file\n");
            rval = 1; goto CLEANUP;
        }
    }

    if (fixedfname) {
        rval = CCutil_getedgelist (ncount, fixedfname, &fixcount,
                                   &fixlist, &fixlen, 0);
        CCcheck_rval (rval, "CCutil_getedgelist failed");
        printf ("Fixed:      %d\n", fixcount); fflush (stdout);

        for (i = 0; i < fixcount; i++) {
            if (fixlen[i] != CCutil_dat_edgelen (fixlist[2*i],
                                                 fixlist[2*i+1], &dat)) {
                fprintf (stderr, "fixed lengths do not match TSPLIB file\n");
                rval = 1; goto CLEANUP;
            }
        }
        CC_IFFREE (fixlen, int);
    }

    if (nonpairsfname) {
        rval = CCelim_read_nonpairs (nonpairsfname, ncount, &nonpaircounts,
                                     &nonpairs);
        CCcheck_rval (rval, "CCelim_read_nonpairs failed");
    }

    rval = run_verify (verifyfname, ncount, ecount, elist, fixcount, fixlist,
                       nonpaircounts, nonpairs, beverbose, &dat, &rstate);
    CCcheck_rval (rval, "run_verify failed");

CLEANUP:
    CCutil_freedatagroup (&dat);
    CC_IFFREE (elist, int);
    CC_IFFREE (elen, int);
    CC_IFFREE (fixlist, int);
    CC_IFFREE (fixlen, int);
    free_pair_info (ncount, &nonpaircounts, &nonpairs);
    return rval;
}

static int run_verify (char *fname, int ncount, int ecount, int *elist,
        int fixcount, int *fixlist, int *nonpcounts, int **nonpairs, int loud,
        CCdatagroup *dat, CCrandstate *rstate)
{
    int rval = 0, count = 0, a, b, c = -1, i, yesno = 0, verified_count = 0;
    int *pcounts = (int *) NULL, **pairs = (int **) NULL, htcounter = 0;
    int count_edge = 0, count_fix = 0, count_pair = 0;
    int hamcounter = 0, bad_htcounter = 0, donepaircount = 0;
    int nremain = 0, *remainlist = (int *) NULL, *donepairlist = (int *) NULL;
    int newfix = 0, *newfixlist = (int *) NULL;
    int *outpcounts = (int *) NULL, **outpairs = (int **) NULL;
    double ezeit = 0.0, szeit = CCutil_zeit ();
    FILE *in = (FILE *) NULL;
    CC_SFILE *sin = (CC_SFILE *) NULL;
    CCelim_httree *ht = (CCelim_httree *) NULL;
    CCelim_graph G, fixG, *F = &fixG, currentF;
    CCelim_distobj D;

    CCelim_init_graph (&currentF);  /* stores new fixed edges */
    rval = CCelim_alloc_elim (ncount, ecount, elist, fixcount, fixlist, dat,
               &G, F, &D, (CCkdtree **) NULL, (CCdatagroup **) NULL,
               (CCkdtree **) NULL, rstate);
    CCcheck_rval (rval, "CCelim_alloc_elim failed");
    rval = CCelim_build_graph (&currentF, ncount, fixcount, fixlist,
                               (int *) NULL, &D);
    CCcheck_rval (rval, "CCelim_build_graph failed\n");

    if (nonpairs) {
        rval = CCelim_nonpairs_to_pairs (&G, nonpcounts, nonpairs, &pcounts,
                                         &pairs);
        CCcheck_rval (rval, "CCelim_nonpairs_to_pairs failed");
        a = b = 0;
        for (i = 0; i < ncount; i++) {
            a += ((G.nodelist[i].deg * (G.nodelist[i].deg - 1)) / 2);
            b += pcounts[i];
        }
        printf ("Pairs:      %d (%.1f%%)\n", b, ((double)b/(double)a)*100.0);
        fflush (stdout);
    }

    rval = get_file_counts (fname, binary_in, &count_edge, &count_fix,
                            &count_pair);
    CCcheck_rval (rval, "get_file_counts failed\n");
    if (count_edge) printf ("Verify Edges: %d\n", count_edge);
    if (count_fix)  printf ("Verify Fixed: %d\n", count_fix);
    if (count_pair) printf ("Verify Pairs: %d\n", count_pair);
    fflush (stdout);

    if (count_pair) {
        CC_MALLOC (donepairlist, 3*count_pair, int);
    }

    if (binary_in == 0) {
        if ((in = fopen (fname, "r")) == (FILE *) NULL) {
            fprintf (stderr, "could not open %s for reading\n", fname);
            rval = 1; goto CLEANUP;
        }
    } else {
        sin = CCutil_sopen (fname, "r");
        CCcheck_NULL (sin, "could not open file for binary read");
    }

    do {
        rval = CCelim_read_httree (&ht, in, sin);
        CCcheck_rval (rval, "CCelim_read_httree failed");
        if (!ht) break;

        if (loud) {
            switch (ht->elimtype) {
            case CCelim_HTTREE_EDGE: printf ("EDGE "); break;
            case CCelim_HTTREE_FIX:  printf ("FIX "); break;
            case CCelim_HTTREE_PAIR: printf ("PAIR "); break;
            default:
                fprintf (stderr, "unknown elim type %d\n", ht->elimtype);
                rval = 1; goto CLEANUP;
            }
        }

        a = ht->root->hamilton_nodes[0];
        b = ht->root->hamilton_nodes[1];
        if (ht->elimtype == CCelim_HTTREE_PAIR) {
            c = ht->root->hamilton_nodes[2];
        }
        if (loud) {
            if (ht->elimtype == CCelim_HTTREE_PAIR) {
                printf ("[%d,%d,%d]", a, b, c);
            } else {
                printf ("[%d,%d]", a, b);
            }
            printf (", count = %d, depth = %d\n", ht->count, ht->depth);
            fflush (stdout);
        }
        htcounter = 1;
        hamcounter = 0;

        if (ht->elimtype == CCelim_HTTREE_EDGE) {
            rval = verify_edge (a, b, &G, F, &D, &yesno, &ezeit,
                pcounts, pairs, ht->longpath, ht->use_tsp, ht,
                &htcounter, &hamcounter);
            CCcheck_rval (rval, "verify_edge failed");
            if (yesno) CCelim_delete_edge (a, b, &G);
        } else if (ht->elimtype == CCelim_HTTREE_FIX) {
            rval = verify_fixed_edge (a, b, &G, F, &D, &yesno, &ezeit,
                pcounts, pairs, ht->use_tsp, ht, &htcounter, &hamcounter);
            CCcheck_rval (rval, "verify_fixed_edge failed");
            if (yesno) CCelim_add_edge (a, b, CCelim_dist(a, b, &D), &currentF);
        } else if (ht->elimtype == CCelim_HTTREE_PAIR) {
            rval = verify_pair (a, b, c, &G, F, &D, &yesno, &ezeit,
                pcounts, pairs, ht->longpath, ht->use_tsp, ht,
                &htcounter, &hamcounter);
                CCcheck_rval (rval, "verify_pair failed");
        } else {
            fprintf (stderr, "unknown elim type %d\n", ht->elimtype);
            rval = 1; goto CLEANUP;
        }

        if (yesno) {
            verified_count++;
            if (ht->elimtype == CCelim_HTTREE_PAIR) {
                donepairlist[3*donepaircount]   = b;
                donepairlist[3*donepaircount+1] = a;
                donepairlist[3*donepaircount+2] = c;
                donepaircount++;
            }
        }

        if (loud) {
            if (yesno) printf ("Verified");
            else       printf ("Not Verified");
            printf (", %d HT-nodes, %d H-leafs, %.2f seconds\n",
                         htcounter, hamcounter, ezeit);
        } else {
            printf ("%d", yesno); fflush (stdout);
            if (count % 1000 == 999) {
                printf ("\nVerifed %d/%d (%.2f seconds)\n", verified_count,
                                           count+1, CCutil_zeit() - szeit);
                fflush (stdout);
            }
            if (yesno == 0) {
                if (ht->elimtype == CCelim_HTTREE_EDGE) {
                    printf ("\nDid not verify: EDGE %d %d\n", a, b);
                } else if (ht->elimtype == CCelim_HTTREE_FIX) {
                    printf ("\nDid not verify: FIX %d %d\n", a, b);
                } else {
                    printf ("\nDid not verify: PAIR %d %d %d\n", a, b, c);
                }
                fflush (stdout);
            }
        }
        if (show_trace_miss && ht->count != htcounter) {
            printf ("\nNOTE: HT nodes traced does not match ht->count\n");
            fflush (stdout);
            bad_htcounter++;
        }

        CCelim_free_httree (ht);
        CC_IFFREE (ht, CCelim_httree);
        count++;
    } while (1);

    printf ("Number of HT trees: %d\n", count);
    printf ("Verified eliminations: %d\n", verified_count);
    printf ("Running Time: %.2f seconds\n", CCutil_zeit() - szeit);
    if (show_trace_miss) {
        if (bad_htcounter) {
            printf ("WARNING: %d did not hit the Tutte-node count\n",
                     bad_htcounter);
        } else {
            printf ("All runs matched the Tutte-node counts\n");
        }
    }
    fflush (stdout);

    if (outedgefname) {
        printf ("Writing remaining edges to file %s\n", outedgefname);
        fflush (stdout);
        rval = CCelim_getedges_graph (&G, &nremain, &remainlist);
        CCcheck_rval (rval, "CCelim_getedges_graph failed");
        printf ("Remaining:  %d\n", nremain);
        rval = CCutil_writeedges (ncount, outedgefname, nremain, remainlist,
                                  dat, 0);
        CCcheck_rval (rval, "CCtuil_writeedges failed");
    }

    if (outfixfname) {
        printf ("Writing fixed edges to file %s\n", outfixfname);
        fflush (stdout);
        rval = CCelim_getedges_graph (&currentF, &newfix, &newfixlist);
        CCcheck_rval (rval, "CCelim_getedges_graph failed");
        printf ("Fixed:  %d\n", newfix);
        rval = CCutil_writeedges (ncount, outfixfname, newfix, newfixlist,
                                  dat, 0);
        CCcheck_rval (rval, "CCtuil_writeedges failed");
    }

    if (outpairfname && donepaircount > 0) {
        printf ("Writing non-pairs to file %s\n", outpairfname);
        fflush (stdout);
        rval = CCelim_pairlist_to_nonpairs (&G, donepaircount, donepairlist,
                            nonpcounts, nonpairs, &outpcounts, &outpairs);
        CCcheck_rval (rval, "pairlist_to_nonpairs failed");
        rval = CCelim_write_nonpairs (outpairfname, ncount, outpcounts,
                                      outpairs);
        CCcheck_rval (rval, "CCelim_write_nonpairs failed");
    }

CLEANUP:
    CCelim_free_elim (&G, F, &D, (CCkdtree **) NULL, (CCdatagroup **) NULL,
                     (CCkdtree **) NULL);
    CCelim_free_graph (&currentF);
    free_pair_info (ncount, &pcounts, &pairs);
    free_pair_info (ncount, &outpcounts, &outpairs);
    if (in) fclose (in);
    if (sin) CCutil_sclose (sin);
    CCelim_free_httree (ht);
    CC_IFFREE (ht, CCelim_httree);
    CC_IFFREE (remainlist, int);
    CC_IFFREE (newfixlist, int);
    CC_IFFREE (donepairlist, int);
    return rval;
}

static int verify_edge (int a, int b, CCelim_graph *G, CCelim_graph *F,
        CCelim_distobj *D, int *yesno, double *ezeit, int *pcounts, int **pairs,
        int longpath, int use_tsp, CCelim_httree *vt,
        int *htcounter, int *hamcounter)
{
    int rval = 0;
    double qzeit = CCutil_zeit();

    *yesno = 0;
    if (F && CCelim_edge_in_graph (a, b, F)) goto CLEANUP;  /* fixed edge */

    /* check if edge ab appears in pairs lists for a and b */
    if (pcounts && pairs) {
        int *list = (int *) NULL, acnt, bcnt;
        CC_MALLOC (list, G->nodelist[a].deg + G->nodelist[b].deg, int);
        CCelim_get_partner_pairs (G, a, b, pcounts, pairs, &acnt, list);
        CCelim_get_partner_pairs (G, b, a, pcounts, pairs, &bcnt, list);
        CC_FREE (list, int);

        if (acnt == 0 || bcnt == 0) {
            *yesno = 1; goto CLEANUP;
        }
    }

    if (vt->root->tutte_type == CCelim_TUTTE_PATH) {
        int k = (longpath < 1 ? 3 : longpath);  /* BICO: was < 3, close */
        rval = verify_via_long_path (k, a, b, D, G, yesno,
                                     F, pcounts, pairs, elim_time_limit,
                                     use_tsp, vt->root, htcounter, hamcounter);
        CCcheck_rval (rval, "verify_via_long_path failed");
    } else if (vt->root->tutte_type == CCelim_TUTTE_CD) {
        int c = vt->root->tutte_nodes[0], d = vt->root->tutte_nodes[1];
        rval = verify_with_cd (a, b, c, d, D, G, yesno, F, pcounts,
              pairs, elim_time_limit, use_tsp, vt->root, htcounter, hamcounter);
        CCcheck_rval (rval, "verify_with_cd failed");
    } else {
        fprintf (stderr, "Tutte type for elim is not path or cd\n");
        rval = 1; goto CLEANUP;
    }

CLEANUP:
    *ezeit = CCutil_zeit() - qzeit;
    return rval;
}

static int verify_fixed_edge (int a, int b, CCelim_graph *G, CCelim_graph *F,
        CCelim_distobj *D, int *yesno, double *ezeit, int *pcounts,
        int **pairs, int use_tsp, CCelim_httree *vt, int *htcounter,
        int *hamcounter)
{
    int rval = 0, ham_nodes[6], i, j, k, good, depth = 2;
    int *alist = (int *) NULL, *blist = (int *) NULL, acount, bcount;
    double qzeit = CCutil_zeit();
    CCelim_path psys[CCelim_PMAX];
    CCelim_htnode *vtchild = (CCelim_htnode *) NULL;
   
    *yesno = 0;

    if (F && CCelim_edge_in_graph (a, b, F)) {
        printf ("\nEdge already fixed\n"); fflush (stdout);
        *yesno = 1; goto CLEANUP;
    } 

    if (!CCelim_edge_in_graph (a, b, G)) {
        printf ("\nTrying to fix edge that is not in graph\n"); fflush (stdout);
        goto CLEANUP;
    }

    if (vt->root->tutte_type != CCelim_TUTTE_CD) {
        fprintf (stderr, "Root node is not of Tutte Type cd\n");
        rval = 1; goto CLEANUP;
    }

    if ((a != vt->root->tutte_nodes[0] || b!=vt->root->tutte_nodes[1]) &&
        (a != vt->root->tutte_nodes[1] || b!=vt->root->tutte_nodes[0])) {
        fprintf (stderr, "Root Tutte nodes are not ends of the edge\n");
        rval = 1; goto CLEANUP;
    }

    if (!pairs) {
        CC_MALLOC(alist, 2 * G->nodelist[a].deg * G->nodelist[a].deg,int);
        CC_MALLOC(blist, 2 * G->nodelist[b].deg * G->nodelist[b].deg,int);
        CCelim_all_neighbor_pairs (a, G, F, &acount, alist);
        CCelim_all_neighbor_pairs (b, G, F, &bcount, blist);
    } else {
        alist = pairs[a];  acount = pcounts[a];
        blist = pairs[b];  bcount = pcounts[b];
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
        for (k = 0; k < 3; k++) ham_nodes[k] = psys[0].node[k];
        for (j = 0; j < bcount && good; j++) {
            if (blist[2*j] == a || blist[2*j+1] == a) continue;
            psys[1].node[0] = blist[2*j];
            psys[1].node[1] = b;
            psys[1].node[2] = blist[2*j+1];
            for (k = 0; k < 3; k++) ham_nodes[k+3] = psys[1].node[k];

            CCelim_hamilton_find (vt->root, 3, 6, ham_nodes, &vtchild);
            if (vtchild) (*htcounter)++;
            else         (*hamcounter)++;

            rval = verify_full_path_check (psys, &good, (int *) NULL, D, G,
                  F, pcounts, pairs, elim_time_limit, use_tsp,
                  depth, vtchild, htcounter, hamcounter);
            CCcheck_rval (rval, "full_path_check failed");
        }
    }
    *yesno = good;


CLEANUP:
    *ezeit = CCutil_zeit() - qzeit;
    if (!pairs) {
        CC_IFFREE (alist, int);
        CC_IFFREE (blist, int);
    }
    return rval;
}

static int verify_pair (int a, int b, int c, CCelim_graph *G, CCelim_graph *F,
        CCelim_distobj *D, int *yesno, double *ezeit, int *pcounts,
        int **pairs, int longpath, int use_tsp, CCelim_httree *vt,
        int *htcounter, int *hamcounter)
{
    int rval = 0, j, depth = 1, dcount = 0, *dlist = (int *) NULL;
    double qzeit = CCutil_zeit();
    CCelim_htnode *vtchild = (CCelim_htnode *) NULL;

    *yesno = 0;

    if (F && F->nodelist[b].deg > 0) {
        if (!CCelim_edge_in_graph (a, b, F) &&
            !CCelim_edge_in_graph (b, c, F)) {
            *yesno = 1; goto CLEANUP;
        }
    }

    if (vt->root->tutte_type == CCelim_TUTTE_PATH) {
        int k = (longpath < 3 ? 3 : longpath);
        rval = verify_pair_via_long_path (k, a, b, c, D, G, yesno,
                                     F, pcounts, pairs, elim_time_limit,
                                     use_tsp, vt->root, htcounter, hamcounter);
        CCcheck_rval (rval, "verify_pair_via_long_path failed");
    } else if (vt->root->tutte_type == CCelim_TUTTE_POINT) {
        int d = vt->root->tutte_nodes[0];
        CCelim_path psys[CCelim_PMAX];
        int good = 1;

        CC_MALLOC (dlist, 2*CCelim_MAX_DEGREE*CCelim_MAX_DEGREE, int);
        CCelim_pair_neighbors_three_swap (a, b, c, d, G, D, &dcount, dlist,
                                          pcounts, pairs);

        psys[0].next = &(psys[1]);
        psys[1].next = (CCelim_path *) NULL;
        psys[0].cnt = 3;
        psys[1].cnt = 3;
        psys[0].node[0] = a; psys[0].node[1] = b; psys[0].node[2] = c;
        for (j = 0; j < dcount && good; j++) {
            psys[1].node[0] = dlist[2*j];
            psys[1].node[1] = d;
            psys[1].node[2] = dlist[2*j+1];

            CCelim_hamilton_find (vt->root, CCelim_HAMILTON_PATH, 3,
                                  psys[1].node, &vtchild);
            if (vtchild) (*htcounter)++;
            else         (*hamcounter)++;

            rval = verify_full_path_check (psys, &good, (int *) NULL, D, G,
                  F, pcounts, pairs, elim_time_limit, use_tsp,
                  depth, vtchild, htcounter, hamcounter);
            CCcheck_rval (rval, "verify_full_path_check failed");
        }
        *yesno = good;
    } else {
        fprintf (stderr,"pair Tutte type not path or point\n");
        rval = 1; goto CLEANUP;
    }

CLEANUP:
    *ezeit = CCutil_zeit() - qzeit;
    CC_IFFREE (dlist, int);
    return rval;
}

/* verify_with_cd attempts to show that edge ab is not in any optimal tour  */
/* by building an elimination-search tree starting with ab and all pairs of */
/* edges meeting c and all pairs of edges meeting d, so all pairs of pairs  */

static int verify_with_cd (int a, int b, int c, int d, CCelim_distobj *D,
        CCelim_graph *G, int *yesno, CCelim_graph *F,
        int *pcounts, int **pairs, double timelimit, int use_tsp,
        CCelim_htnode *vt, int *htcounter, int *hamcounter)
{
    int rval = 0, ccount, dcount, i, j, k, good, ham_nodes[6], ab[2];
    int *clist = (int *) NULL, *dlist = (int *) NULL, depth = 2;
    CCelim_path psys[CCelim_PMAX];
    CCelim_htnode *vtchild = (CCelim_htnode *) NULL;

    *yesno = 0;

    /* first verify that edge cd is not compatiable with edge ab */

    if (CCelim_edge_in_graph (c, d, G) && compatible_test (a, b, c, d, D)) {
        fprintf (stderr, "ERROR: cd is compatible with ab\n");
        rval = 1; goto CLEANUP;
    }

    ab[0] = a; ab[1] = b;
    CC_MALLOC (clist, 2 * G->nodelist[c].deg * G->nodelist[c].deg, int);
    CC_MALLOC (dlist, 2 * G->nodelist[d].deg * G->nodelist[d].deg, int);
    CCelim_check_neighbors_three_swap (a, b, c, G, D, &ccount, clist, F,
                                       pcounts, pairs);
    CCelim_check_neighbors_three_swap (a, b, d, G, D, &dcount, dlist, F,
                                       pcounts, pairs);
    if (ccount == 0 || dcount == 0) {
        *yesno = 1; goto CLEANUP;
    }

    psys[0].next = &(psys[1]);
    psys[1].next = &(psys[2]);
    psys[2].next = (CCelim_path *) NULL;
    psys[0].cnt = 2; psys[1].cnt = psys[2].cnt = 3;
    psys[0].node[0] = a; psys[0].node[1] = b;

    good = 1;
    for (i = 0; i < ccount && good; i++) {
        psys[1].node[0] = clist[2*i];
        psys[1].node[1] = c;
        psys[1].node[2] = clist[2*i+1];
        for (k = 0; k < 3; k++) ham_nodes[k] = psys[1].node[k];
        for (j = 0; j < dcount && good; j++) {
            psys[2].node[0] = dlist[2*j];
            psys[2].node[1] = d;
            psys[2].node[2] = dlist[2*j+1];
            for (k = 0; k < 3; k++) ham_nodes[k+3] = psys[2].node[k];

            CCelim_hamilton_find (vt, 3, 6, ham_nodes, &vtchild);
            if (vtchild) (*htcounter)++;
            else         (*hamcounter)++;

            rval = verify_full_path_check (psys, &good, ab, D, G,
                  F, pcounts, pairs, timelimit, use_tsp, depth, vtchild,
                  htcounter, hamcounter);
            CCcheck_rval (rval, "verify_full_path_check failed");
        }
    }
    *yesno = good;

CLEANUP:
    CC_IFFREE (clist, int);
    CC_IFFREE (dlist, int);
    return rval;
}

static int verify_via_long_path (int path_len, int a, int b,
        CCelim_distobj *D, CCelim_graph *G, int *yesno,
        CCelim_graph *F, int *pcounts, int **pairs, double timelimit,
        int use_tsp, CCelim_htnode *vt, int *htcounter, int *hamcounter)
{
    int rval = 0, ab[2];
    CCelim_path psys[CCelim_PMAX];

    *yesno = 0;

    ab[0] = a; ab[1] = b;
    psys[0].next = (CCelim_path *) NULL;
    psys[0].cnt = 2;
    psys[0].node[0] = a; psys[0].node[1] = b;

    rval = verify_grow_long_path (path_len, ab, D, G, yesno, F,
                                  pcounts, pairs, psys, timelimit, use_tsp, vt,
                                  htcounter, hamcounter);
    CCcheck_rval (rval, "verify_grow_long_path failed");

CLEANUP:
    return rval;
}

static int verify_pair_via_long_path (int path_len, int a, int b, int c,
        CCelim_distobj *D, CCelim_graph *G, int *yesno,
        CCelim_graph *F, int *pcounts, int **pairs, double timelimit,
        int use_tsp, CCelim_htnode *vt, int *htcounter, int *hamcounter)
{
    int rval = 0;
    CCelim_path psys[CCelim_PMAX];

    *yesno = 0;

    if (path_len <= 2) goto CLEANUP;

    psys[0].next = (CCelim_path *) NULL;
    psys[0].cnt = 3;
    psys[0].node[0] = a;
    psys[0].node[1] = b;
    psys[0].node[2] = c;

    rval = verify_grow_long_path (path_len, (int *) NULL, D, G, yesno,
                                  F, pcounts, pairs, psys, timelimit, use_tsp,
                                  vt, htcounter, hamcounter);
    CCcheck_rval (rval, "verify_grow_long_path failed");

CLEANUP:
    return rval;
}

static int verify_grow_long_path (int path_len, int *ab, CCelim_distobj *D,
        CCelim_graph *G, int *yesno, CCelim_graph *F,
        int *pcounts, int **pairs, CCelim_path *psys, double timelimit,
        int use_tsp, CCelim_htnode *vt, int *htcounter, int *hamcounter)
{
    int rval = 0, good, i, j, c, d, e, ccount = 0, *clist = (int *) NULL;
    int depth = 1;
    CCelim_path preverse[CCelim_PMAX];
    CCelim_htnode *vtchild = (CCelim_htnode *) NULL;

    *yesno = 0;

    if (!vt) {
        printf ("WARNING: verify_grow_long_path without vt node\n");
        fflush (stdout); goto CLEANUP;
    }

    if (psys[0].cnt > path_len) {
        CCelim_hamilton_find (vt, CCelim_HAMILTON_LONG, psys[0].cnt,
                              psys[0].node, &vtchild);
        if (vtchild) (*htcounter)++;
        else         (*hamcounter)++;
        rval = verify_full_path_check (psys, yesno, ab, D, G, F,
                                       pcounts, pairs, timelimit, use_tsp,
                                       depth, vtchild, htcounter, hamcounter);
        CCcheck_rval (rval, "verify_full_path_check failed");
        goto CLEANUP;
    }

    /* reverse the path to keep a-b in the middle */
    CCelim_copy_path_system (psys, preverse, 1);

    /* grow the path by one edge starting at the last city c */
    c = preverse[0].node[preverse[0].cnt-1];
    d = preverse[0].node[preverse[0].cnt-2];

    /* clist will hold the potential new points to extend from c */
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
            rval = verify_grow_long_path (path_len, ab, D, G, &good,
                  F, pcounts, pairs, preverse, timelimit, use_tsp, vt,
                  htcounter, hamcounter);
            CCcheck_rval (rval, "verify_grow_long_path failed");
        }
    }
    *yesno = good;

CLEANUP:
    CC_IFFREE (clist, int);
    return rval;
}

static int verify_full_path_check (CCelim_path *psys, int *yesno, int *ab,
        CCelim_distobj *D, CCelim_graph *G,
        CCelim_graph *F, int *pcounts, int **pairs, double timelimit,
        int use_tsp, int depth, CCelim_htnode *vt, int *htcounter,
        int *hamcounter)
{
    int rval = 0, bad = 0, valid;
    CCelim_path pmerge[CCelim_PMAX];

    *yesno = 0;

    rval = CCelim_validate_paths (G->ncount, psys, &valid, pmerge);
    CCcheck_rval (rval, "CCelim_validate_paths failed");
    if (!valid) { *yesno = 1; goto CLEANUP; }

    rval = CCelim_test_paths (ab, pmerge, D, G, use_tsp, &bad);
    CCcheck_rval (rval, "CCelim_test_paths failed");
    if (!bad) { *yesno = 1; goto CLEANUP; }

    if (!vt) {
        printf ("WARNING: Did not verify with tree nodes only\n");
        {
            CCelim_path *p;
            int i = 0;
            for (p = pmerge; p; p = p->next) {
                int k;
                printf ("Path %d: ", i++);
                for (k = 0; k < p->cnt; k++) printf ("%d ", p->node[k]);
                printf ("\n");
                fflush (stdout);
            }
        }
        fflush (stdout); goto CLEANUP;
    }

    if (vt->tutte_type == CCelim_TUTTE_END) {
        rval = verify_try_extra_edge (ab, pmerge, D, G, yesno,
                                      F, pcounts, pairs, timelimit, use_tsp,
                                      depth+1, vt, htcounter, hamcounter);
        CCcheck_rval (rval, "try_extra_edge failed");
    } else if (vt->tutte_type == CCelim_TUTTE_POINT) {
        rval = verify_try_extra_point (ab, pmerge, D, G, yesno,
                                       F, pcounts, pairs, timelimit, use_tsp,
                                       depth+2, vt, htcounter, hamcounter);
        CCcheck_rval (rval, "verify_try_extra_point failed");
    } else {
        fprintf (stderr, "Tutte move is not an edge or point\n");
        rval = 1; goto CLEANUP;
    }

CLEANUP:
    return rval;
}

static int verify_try_extra_point (int *ab, CCelim_path *pin, CCelim_distobj *D,
        CCelim_graph *G, int *good, CCelim_graph *F,
        int *pcounts, int **pairs, double timelimit, int use_tsp, int depth,
        CCelim_htnode *vt, int *htcounter, int *hamcounter)
{
    int rval = 0, j, f, fcount, *flist = (int *) NULL, okay;
    int a = pin->node[0], b = pin->node[1];
    CCelim_path psys[CCelim_PMAX], pex;
    CCelim_htnode *vtchild  = (CCelim_htnode *) NULL;

    /* note: could set a and be to be the elements of CCelim_array ab */

    *good = 0;

    CC_MALLOC (flist, 2*CCelim_MAX_DEGREE*CCelim_MAX_DEGREE, int);
    CCelim_copy_path_system (pin, psys, 0);

    if (!vt) {
        printf ("WARNING: verify_try_extra_point without vt node\n");
        fflush (stdout); goto CLEANUP;
    }

    if (vt->tutte_type != CCelim_TUTTE_POINT) {
        fprintf (stderr, "Tutte node not of single point\n");
        rval = 1; goto CLEANUP;
    }

    f = vt->tutte_nodes[0];
    CCelim_check_neighbors_three_swap (a, b, f, G, D, &fcount, flist, F,
                                       pcounts, pairs);

    pex.next = psys;
    pex.cnt = 3;
    okay = 1;

    for (j = 0; j < fcount && okay; j++) {
        pex.node[0] = flist[2*j];
        pex.node[1] = f;
        pex.node[2] = flist[2*j+1];

        CCelim_hamilton_find (vt, 2, 3, pex.node, &vtchild);
        if (vtchild) (*htcounter)++;
        else         (*hamcounter)++;
        rval = verify_full_path_check (&pex, &okay, ab, D, G,
                                       F, pcounts, pairs, timelimit, use_tsp,
                                       depth, vtchild, htcounter, hamcounter);
        CCcheck_rval (rval, "verify_full_path_check failed");
    }
    if (okay == 1) *good = 1;

CLEANUP:
    CC_IFFREE (flist, int);
    return rval;
}

static int verify_try_extra_edge (int *ab, CCelim_path *pin, CCelim_distobj *D,
        CCelim_graph *G, int *good, CCelim_graph *F,
        int *pcounts, int **pairs, double timelimit, int use_tsp, int depth,
        CCelim_htnode *vt, int *htcounter, int *hamcounter)
{
    int rval = 0, j, k, okay, c, d, n;
    int bcount = 0, *blist = (int *) NULL, *middle = (int *) NULL;
    CCelim_path psys[CCelim_PMAX], pex, *p;
    CCelim_htnode *vtchild  = (CCelim_htnode *) NULL;

    *good = 0;

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

    pex.next = psys;
    pex.cnt = 2;

    CC_MALLOC (blist, G->ncount, int);

    if (!vt) {
        printf ("WARNING: verify_try_extra_point without vt node\n");
        fflush (stdout); goto CLEANUP;
    }

    if (vt->tutte_type != CCelim_TUTTE_END) {
        fprintf (stderr, "Tutte node not of type end\n");
        rval = 1; goto CLEANUP;
    }

    c = vt->tutte_nodes[0];
    n = -1;
    for (p = psys; p; p = p->next) {
        if (c == p->node[0])        { n = p->node[1]; break; }
        if (c == p->node[p->cnt-1]) { n = p->node[p->cnt-2]; break; }
    }
    if (n == -1) {
        fprintf (stderr, "Verify: Did not find node next to c\n");
        rval = 1; goto CLEANUP;
    }
    CCelim_get_partner_pairs (G, c, n, pcounts, pairs, &bcount, blist);

    okay = 1;
    for (j = 0; j < bcount && okay; j++) {
        d = blist[j];
        if (middle[d] == 1) continue;
        pex.node[0] = c; pex.node[1] = d;

        CCelim_hamilton_find (vt, 1, 2, pex.node, &vtchild);
        if (vtchild) (*htcounter)++;
        else         (*hamcounter)++;
        rval = verify_full_path_check (&pex, &okay, ab, D, G,
                                       F, pcounts, pairs, timelimit, use_tsp,
                                       depth, vtchild, htcounter, hamcounter);
        CCcheck_rval (rval, "full_path_check failed");
    }
    if (okay == 1) *good = 1;

CLEANUP:
    CC_IFFREE (blist, int);
    CC_IFFREE (middle, int);
    return rval;
}

static void free_pair_info (int ncount, int **pcounts, int ***pairs)
{
    if (*pairs) {
        int i;
        for (i = 0; i < ncount; i++) CC_IFFREE ((*pairs)[i], int);
        CC_FREE ((*pairs), int *);
    }
    CC_IFFREE ((*pcounts), int);
}

/* compatible_test checks if edges ab and cd cannot be shown to be        */
/* incompatible by a simple 2-swap (that is, one of the two orientations  */
/* is not an improving 2-swap)                                            */

static int compatible_test (int a, int b, int c, int d, CCelim_distobj *D)
{
    int t = CCelim_dist(a,b,D) + CCelim_dist(c,d,D);
    return (CCelim_dist(a,d,D) + CCelim_dist(b,c,D) >= t ||
            CCelim_dist(a,c,D) + CCelim_dist(b,d,D) >= t) ;
}

static int get_file_counts (char *fname, int binary_in, int *count_edge,
        int *count_fix, int *count_pair)
{
    int rval = 0, counte = 0, countf = 0, countp = 0;
    CCelim_httree *ht = (CCelim_httree *) NULL;
    FILE *in = (FILE *) NULL;
    CC_SFILE *sin = (CC_SFILE *) NULL;

    if (binary_in == 0) {
        if ((in = fopen (fname, "r")) == (FILE *) NULL) {
            fprintf (stderr, "could not open %s for reading\n", fname);
            rval = 1; goto CLEANUP;
        }
    } else {
        sin = CCutil_sopen (fname, "r");
        CCcheck_NULL (sin, "could not open file for binary read");
    }

    do {
        rval = CCelim_read_httree (&ht, in, sin);
        CCcheck_rval (rval, "CCelim_read_httree failed");
        if (!ht) break;

        switch (ht->elimtype) {
        case CCelim_HTTREE_EDGE: counte++; break;
        case CCelim_HTTREE_FIX:  countf++; break;
        case CCelim_HTTREE_PAIR: countp++; break;
        default:
            fprintf (stderr, "unknown elim type %d\n", ht->elimtype);
            rval = 1; goto CLEANUP;
        }

        CCelim_free_httree (ht);
        CC_IFFREE (ht, CCelim_httree);
    } while (1);

    *count_edge = counte;
    *count_fix  = countf;
    *count_pair = countp;

CLEANUP:
    if (in) fclose (in);
    if (sin) CCutil_sclose (sin);
    return rval;
}

static int parseargs (int ac, char **av)
{
    int c, boptind = 1;
    char *boptarg = (char *) NULL;

    while ((c = CCutil_bix_getopt (ac, av, "be:f:F:mp:P:T:vV:", &boptind, &boptarg)) != EOF) {
        switch (c) {
        case 'b': binary_in = 1; break;
        case 'e': outedgefname = boptarg; break;
        case 'f': outfixfname = boptarg; break;
        case 'F': fixedfname = boptarg; break;
        case 'm': show_trace_miss = 1; break;
        case 'p': outpairfname = boptarg; break;
        case 'P': nonpairsfname = boptarg; break;
        case 'T': tspfname = boptarg; break;
        case 'V': verifyfname = boptarg; break;
        case 'v': beverbose = 1; break;
        case CC_BIX_GETOPT_UNKNOWN:
        case '?':
        default:
            usage (av[0]);
            return 1;
        }
    }

    if (boptind < ac) { edgefname = av[boptind++]; }
    else              { usage (av[0]); return 1; }

    return 0;
}

static void usage (char *f)
{
    fprintf (stderr, "Usage: %s [-see below-] -V ht_file -T tsp_file edge_file\n", f);
    fprintf (stderr, "   -b    ht_file is a ht_binary_file\n");
    fprintf (stderr, "   -e f  dump remaining edges to file f\n");
    fprintf (stderr, "   -f f  dump fixed edges to file f\n");
    fprintf (stderr, "   -F f  input fixed edges in edge-file f\n");
    fprintf (stderr, "   -m    show trace misses (not accurate if ht_file is from a parallel run)\n");
    fprintf (stderr, "   -p f  dump non-pairs to files f\n");
    fprintf (stderr, "   -P f  input eliminated pairs in file f\n");
    fprintf (stderr, "   -v    lots of output\n");
}
