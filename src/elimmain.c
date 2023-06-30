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
/*                             Edge Elimination                             */
/*                                                                          */
/*  Local-search technique to fix edges to 0 or 1 in an instance of the     */
/*  TSP. This file contains calling programs for the routines in elim.c.    */
/*  The parallel codes use sockets for communication, employing routines    */
/*  from the Concorde TSP code.                                             */
/*                                                                          */
/*  Written by:  Cook, 2013                                                 */
/*  Updates: June-July 2016, October-November 2017, July 2018,              */
/*           August 2022 - June 2023                                        */
/*                                                                          */
/*  Adopts ideas from the paper "Edge elimination in TSP instances" by      */
/*  S. Hougardy and R. T. Schroeder (2014).                                 */
/*                                                                          */
/****************************************************************************/

#include <stdio.h>
#include "ccutil.h"
#include "elim.h"

#define ELIM_PORT ((unsigned short) 24872)
#define ELIM_HELLO (1)
#define ELIM_WORK  (2)
#define ELIM_DONE  (3)

#define WORK_ELIM (1)
#define WORK_FIX  (2)
#define WORK_PAIR (3)

#define HT_MAX_COUNT 1000000
#define HT_MAX_DEPTH 100

#define DEFAULT_TIME_LIMIT 60.0  /* 10000000.0 */

static char *edgefname = (char *) NULL;
static char *tspfname = (char *) NULL;
static char *outfname = (char *) NULL;
static char *infname = (char *) NULL;
static char *elimfname = (char *) NULL;
static char *fixedfname = (char *) NULL;
static char *nonpairsfname = (char *) NULL;
static char *tourfname = (char *) NULL;
static char *htfname = (char *) NULL;
static int point_levels = 2;
static int beverbose = 0;
static int use_longpath = 0;
static int use_tsp_swapping = 1;
static int use_max_neighborhood = 10;
static int use_fast_elim = 0;
static int use_single_fast = 0;
static int use_level_loop = 0;
static int use_full_loop = 0;
static int witness_type = CCelim_CD_EDGE;
static int only_crossing_edges = 0;
static double elim_time_limit = DEFAULT_TIME_LIMIT;
static int be_nethost = 0;
static char *grunthostname = (char *) NULL;
static int worktype = WORK_ELIM;
static int seed = 0;
static int use_norm = CC_EUCLIDEAN;
static int tsplib_in = 1;
static int save_hamilton_tutte_tree = 0;

int main (int ac, char **av);
static int run_full_loop (int ncount, CCdatagroup *dat, int ecount, int *elist,
    CCrandstate *rstate, int loud, char *probname, int *incyc, int net_boss);
static int run_elim (int gotype, int ncount, int ecount, int *elist,
    int processcount, int *processlist, int fixcount, int *fixlist,
    CCdatagroup *dat, CCrandstate *rstate, int level_count, int wtype,
    int *remaincount, int **remainlist, int *pelimcount, int *elimlist, 
    int loud, int *nonpaircounts, int **nonpairs, int **outpcounts,
    int ***outpairs, double timelimit, int longpath, int use_tsp,
    int max_neighborhood, int call_fast, int save_tree);
static int run_elim_loop (int gotype, int ncount, int ecount, int *elist,
    int processcount, int *processlist, int fixcount, int *fixlist,
    CCdatagroup *dat, CCrandstate *rstate,
    int *remaincount, int **remainlist, int *pelimcount, int *elimlist,
    int loud, int *nonpcounts, int **nonpairs, int **outpcounts,
    int ***outpairs, double timelimit, int longpath, int use_tsp,
    int call_fast, int single_fast);
static int elim_boss (int gotype, int ncount, int ecount, int *elist,
    int processcount, int *processlist, int fixcount, int *fixlist,
    int *nonpcounts, int **nonpairs, CCdatagroup *dat,
    int *remaincount, int **remainlist, int *pelimcount, int *elimlist,
    int **outpcounts, int ***outpairs, int level_count, int wtype,
    int longpath, int use_tsp, double timelimit, int max_neighborhood,
    int use_loop, int call_fast, int single_fast, int loud, 
    CC_SPORT *lport, int *round_id, double *boss_wall_time, int save_tree,
    int *incyc);
static int elim_grunt (char *theboss, CCrandstate *rstate, int loud,
    double local_timelimit);
static int processlist_cycle (int ncount, int *incyc, int worktype,
    int processcount, int *processlist, int *new_processcount,
    int **new_processlist);
static int read_pair_list (char *fname, int ncount, int ecount, int *elist,
    int *pcount, int **plist, int *nonpcounts, int **nonpairs);

static int open_connection (char *theboss, CC_SFILE **s);
static void free_pair_info (int ncount, int **pcounts, int ***pairs);
static int count_crossings (int ncount, CCdatagroup *dat, CCrandstate *rstate,
    int ecount, int *elist, int *ccount, int **clist);
static int segments_intersect (double p0_x, double p0_y, double p1_x,
    double p1_y, double p2_x, double p2_y, double p3_x, double p3_y);
static int parseargs (int ac, char **av);
static void usage (char *f);

int main (int ac, char **av)
{
    int rval = 0, i, k, nremain = 0, *remain = (int *) NULL;
    int ncount = 0, ecount, *elist = (int *) NULL, *elen = (int *) NULL;
    int incount = 0, *inlist = (int *) NULL, *inlen = (int *) NULL;
    int fixincount = 0, *fixinlist = (int *) NULL, *fixinlen = (int *) NULL;
    int nelim = 0, *elimlist = (int *) NULL, *newlist = (int *) NULL;
    int processcount, *processlist = (int *) NULL, *perm = (int *) NULL;
    int *nonpaircounts = (int *) NULL, **nonpairs = (int **) NULL;
    int *outpcounts = (int *) NULL, **outpairs = (int **) NULL;
    int *pairprocesslist = (int *) NULL, temp, *incyc = (int *) NULL;
    int tprocesscount = 0, *tprocesslist = (int *) NULL, round_id = 0;
    char *probname = (char *) NULL;
    CCdatagroup dat;
    CCrandstate rstate;
    double szeit;
    CC_SPORT *lport = (CC_SPORT *) NULL;

    CCutil_init_datagroup (&dat);

    rval = parseargs (ac, av);
    if (rval) goto CLEANUP;

    rval = CCutil_print_command (ac, av);
    CCcheck_rval (rval, "CCutil_print_command failed");
    CCutil_printlabel ();
    CCutil_sprand (seed, &rstate);

    if (grunthostname) {  /* parallel code, process is a worker */
        rval = elim_grunt (grunthostname, &rstate, beverbose, elim_time_limit);
        CCcheck_rval (rval, "elim_grunt failed");
        goto CLEANUP;
    }

    if (tsplib_in == 1) {
        rval = CCutil_gettsplib (tspfname, &ncount, &dat);
        CCcheck_rval (rval, "CCutil_gettsplib failed");
    } else {
        rval = CCutil_getdata (tspfname, 0, use_norm, &ncount, &dat, 1000, 0,
                               &rstate);
        CCcheck_rval (rval, "CCutil_getdata failed");
    }

    probname = CCutil_problabel (tspfname);
    printf ("Prob File Name: %s\n", probname);

    if (save_hamilton_tutte_tree && (use_full_loop || use_level_loop)) {
        fprintf (stderr, "cannot save H-T tree with loops\n");
        rval = 1; goto CLEANUP;
    }

    if (save_hamilton_tutte_tree &&
        use_longpath > CCelim_HTNODE_MAX_PATH_LEN) {
        fprintf (stderr, "cannot save H-T tree for paths longer that %d\n",
                          CCelim_HTNODE_MAX_PATH_LEN);
        rval = 1; goto CLEANUP;
    }

    rval = CCutil_getedgelist (ncount, edgefname, &ecount, &elist, &elen, 0);
    CCcheck_rval (rval, "CCutil_getedgelist failed");
    printf ("Edges:      %d\n", ecount); fflush (stdout);
    for (i = 0; i < ecount; i++) {
        if (elen[i] != CCutil_dat_edgelen (elist[2*i], elist[2*i+1], &dat)) {
            fprintf (stderr, "edge lengths do not match TSPLIB file\n");
            rval = 1; goto CLEANUP;
        }
    }

    if (tourfname) {
        CC_MALLOC (incyc, ncount, int);
        rval = CCutil_getcycle_tsplib (ncount, tourfname, incyc);
        CCcheck_rval (rval, "CCutil_getcycle_tsplib failed");
    }

    if (use_full_loop) {   /* run fast-elim, pairs, elim, and fixed */
        szeit = CCutil_zeit ();
        rval = run_full_loop (ncount, &dat, ecount, elist, &rstate,
                              beverbose, probname, incyc, be_nethost);
        CCcheck_rval (rval, "run_full_loop failed");
        goto CLEANUP;
    }

    if (nonpairsfname) {
        rval = CCelim_read_nonpairs (nonpairsfname, ncount, &nonpaircounts,
                                     &nonpairs);
        CCcheck_rval (rval, "CCelim_read_nonpairs failed");
    }

    /* build the list of edges (or pairs) to process. This is either the  */
    /* full edge set, an edge set specified by -I, or edges that cross in */
    /* the embedding determined by the x-y coordinates for geom instances */

    if (worktype == WORK_PAIR) {
        if (infname) {
            rval = read_pair_list (infname, ncount, ecount, elist,
                   &processcount, &pairprocesslist, nonpaircounts, nonpairs);
            CCcheck_rval (rval, "read_pair_list failed");
            processlist = pairprocesslist;
        } else {
            rval = CCelim_build_pairprocesslist (ncount, ecount, elist,
                      nonpaircounts, nonpairs, &processcount, &pairprocesslist);
            CCcheck_rval (rval, "CCelim_build_pairprocesslist failed");
            processlist = pairprocesslist;
        }
    } else if (infname) {
        rval = CCutil_getedgelist (ncount,infname,&incount,&inlist,&inlen,0);
        CCcheck_rval (rval, "CCutil_getedgelist failed");
        for (i = 0; i < incount; i++) {
            if (inlen[i] != CCutil_dat_edgelen (inlist[2*i], inlist[2*i+1],
                                                &dat)) {
                fprintf (stderr, "process lengths do not match TSPLIB file\n");
                rval = 1; goto CLEANUP;
            }
        }
        processcount = incount;
        processlist = inlist;
    } else if (only_crossing_edges) {
        rval = count_crossings (ncount, &dat, &rstate, ecount, elist, &incount,
                                &inlist);
        CCcheck_rval (rval, "count_crossings failed");
        printf ("Crossing:   %d\n", incount); fflush (stdout);
        processcount = incount;
        processlist = inlist;
    } else {
        processcount = ecount;
        processlist = elist;
    }

    /* if a tour is specified with -t, then only try to eliminate edges or  */
    /* pairs not in the tour (or only fix edges that are in the tour).      */

    if (incyc) {
        rval = processlist_cycle (ncount, incyc, worktype, processcount,
                                  processlist, &tprocesscount, &tprocesslist);
        CCcheck_rval (rval, "processlist_cycle failed");
        processcount = tprocesscount;
        processlist = tprocesslist;
        if (worktype == WORK_PAIR) {
            printf ("Using tour, processlist now has %d pairs\n", processcount);
        } else {
            printf ("Using tour, processlist now has %d edges\n", processcount);
        }
        fflush (stdout);
    }

    /* if a seed is specified with -s, then randomize the process order.   */

    if (seed) {  /* randomize the process order */
        if (worktype == WORK_PAIR) {
            CC_MALLOC (newlist, 3*processcount, int);
        } else {
            CC_MALLOC (newlist, 2*processcount, int);
        }
        CC_MALLOC (perm, processcount, int);
        for (i = 0; i < processcount; i++) perm[i] = i;
        for (i = processcount; i > 1; i--) {
            k = CCutil_lprand (&rstate) % i;
            CC_SWAP (perm[i - 1], perm[k], temp);
        }
        for (i = 0; i < processcount; i++) {
            k = perm[i];
            if (worktype == WORK_PAIR) {
                newlist[3*i]   = processlist[3*k];
                newlist[3*i+1] = processlist[3*k+1];
                newlist[3*i+2] = processlist[3*k+2];
            } else {
                newlist[2*i]   = processlist[2*k];
                newlist[2*i+1] = processlist[2*k+1];
            }
        }
        processlist = newlist;
    }

    /* read in edges already fixed to 1 as specifed by the file given in -F */

    if (fixedfname) {
        rval = CCutil_getedgelist (ncount, fixedfname, &fixincount,
                                   &fixinlist, &fixinlen, 0);
        CCcheck_rval (rval, "CCutil_getedgelist failed");
        printf ("Fixed:      %d\n", fixincount); fflush (stdout);

        for (i = 0; i < fixincount; i++) {
            if (fixinlen[i] != CCutil_dat_edgelen (fixinlist[2*i],
                                                   fixinlist[2*i+1], &dat)) {
                fprintf (stderr, "fixed lengths do not match TSPLIB file\n");
                rval = 1; goto CLEANUP;
            }
        }
    }

    CC_MALLOC (elimlist, 2*ecount, int);  /* For eliminated/fixed edges */

    if (be_nethost == 1) {  /* process is the boss for parallel run */
        lport = CCutil_snet_listen (ELIM_PORT);
        if (lport == (CC_SPORT *) NULL) {
            fprintf (stderr, "CCutil_snet_listen failed\n");
            rval = 1; goto CLEANUP;
        }
        rval = elim_boss (worktype, ncount, ecount, elist,
            processcount, processlist, fixincount, fixinlist,
            nonpaircounts, nonpairs, &dat, &nremain, &remain,
            &nelim, elimlist, &outpcounts, &outpairs,
            point_levels-3, witness_type, use_longpath, use_tsp_swapping,
            elim_time_limit, use_max_neighborhood, use_level_loop,
            use_fast_elim, use_single_fast, beverbose, lport,
            &round_id, (double *) NULL, save_hamilton_tutte_tree, incyc);
        CCcheck_rval (rval, "elim_boss failed");
        CCutil_snet_unlisten (lport);
    } else if (use_level_loop) {  /* make repeated runs, increasing params */
        szeit = CCutil_zeit ();
        rval = run_elim_loop (worktype, ncount, ecount, elist, processcount,
            processlist, fixincount, fixinlist, &dat, &rstate,
            &nremain, &remain, &nelim, elimlist, beverbose,
            nonpaircounts, nonpairs, &outpcounts, &outpairs, elim_time_limit,
            use_longpath, use_tsp_swapping, use_fast_elim, use_single_fast);
        CCcheck_rval (rval, "run_elim_loop failed");
        printf ("Time: %.2f seconds\n", CCutil_zeit () - szeit);
        fflush (stdout);
    } else {  /* make the single run with command-line parms */
        szeit = CCutil_zeit ();
        rval = run_elim (worktype, ncount, ecount, elist, processcount,
            processlist, fixincount, fixinlist, &dat, &rstate,
            point_levels-3, witness_type, &nremain, &remain,
            &nelim, elimlist, beverbose, nonpaircounts, nonpairs, &outpcounts,
            &outpairs, elim_time_limit, use_longpath, use_tsp_swapping,
            use_max_neighborhood, use_fast_elim, save_hamilton_tutte_tree);
        CCcheck_rval (rval, "run_elim failed");
        printf ("Time: %.2f seconds\n", CCutil_zeit () - szeit);
        fflush (stdout);
    }

    /* if a file is specifed by -o, write the remaining nonpairs or edges. */
    /* note: if fixing edges to 1, the fixed list is the remain array      */

    if (outfname) {
        if (worktype == WORK_PAIR) {
            rval = CCelim_write_nonpairs (outfname,ncount,outpcounts,outpairs);
            CCcheck_rval (rval, "CCelim_write_nonpairs failed");
        } else {
            rval = CCutil_writeedges (ncount,outfname,nremain,remain,&dat,0);
            CCcheck_rval (rval, "CCtuil_writeedges failed");
        }
    }

    /* always save the non-eliminated and fixed edges to default files */

    if (worktype == WORK_ELIM) {
        char buf[1024];
        sprintf (buf, "%s.remain", probname);
        printf ("Writing remaining edges to %s\n", buf); fflush (stdout);
        rval = CCutil_writeedges (ncount, buf, nremain, remain, &dat, 0);
        CCcheck_rval (rval, "CCtuil_writeedges failed");
    } else if (worktype == WORK_FIX) {
        char buf[1024];
        sprintf (buf, "%s.fix", probname);
        printf ("Writing fixed edges to %s\n", buf); fflush (stdout);
        rval = CCutil_writeedges (ncount, buf, nremain, remain, &dat, 0);
        CCcheck_rval (rval, "CCtuil_writeedges failed");
    }

    /* if a file is specifed by -e, write the eliminated edges */

    if (elimfname) {
        if (worktype == WORK_PAIR) {
            printf ("Warning: no function to output newly eliminated pairs\n");
        } else {
            rval = CCutil_writeedges (ncount,elimfname,nelim,elimlist,&dat,0);
            CCcheck_rval (rval, "CCtuil_writeedges failed");
        }
    }

CLEANUP:
    CCutil_freedatagroup (&dat);
    CC_IFFREE (probname, char);
    CC_IFFREE (elist, int);
    CC_IFFREE (elen, int);
    CC_IFFREE (inlist, int);
    CC_IFFREE (inlen, int);
    CC_IFFREE (fixinlist, int);
    CC_IFFREE (fixinlen, int);
    CC_IFFREE (remain, int);
    CC_IFFREE (elimlist, int);
    CC_IFFREE (pairprocesslist, int);
    CC_IFFREE (newlist, int);
    CC_IFFREE (tprocesslist, int);
    CC_IFFREE (perm, int);
    CC_IFFREE (incyc, int);
    free_pair_info (ncount, &nonpaircounts, &nonpairs);
    free_pair_info (ncount, &outpcounts, &outpairs);
    return rval;
}

/* run_full_loop will first process edges with the fast_elim code in a    */
/* loop that gradually increases the parameter settings, then eliminates  */
/* pairs of edges (that is, determines pairs (a, b) and (a, c) such that  */
/* the path b-a-c cannot be in an optimal tour), then calls the default   */
/* elim code (again in a loop), followed by code to fix edges to one. to  */
/* adjust parameter settings, you can modify code in the function         */
/* run_elim_loop                                                          */

static int run_full_loop (int ncount, CCdatagroup *dat, int ecount, int *elist,
        CCrandstate *rstate, int loud, char *probname, int *incyc,
        int net_boss)
{
    int rval = 0, use_tsp = 1, longpath = 0, round_id = 0;
    int fast_nremain = 0, *fast_remain = (int *) NULL;
    int loop_nremain = 0, *loop_remain = (int *) NULL;
    int nfixed = 0, *fixedlist = (int *) NULL;
    int *nonpcounts = (int *) NULL, **nonpairs = (int **) NULL;
    int *elimlist = (int *) NULL;
    int processcount = ecount, *processlist = elist;
    int tprocesscount, *tprocesslist = (int *) NULL;
    int *pairprocesslist = (int *) NULL;
    int fast_nelim = 0, loop_nelim = 0, pair_nelim = 0;
    char buf[1024];
    double szeit, fast_time, pair_time, loop_time, fix_time, start_szeit;
    CC_SPORT *lport = (CC_SPORT *) NULL;

    if (net_boss == 0) {
        printf ("Run fast loop, pairs loop, elim loop, and fix loop\n");
        fflush (stdout);
        start_szeit = CCutil_zeit ();
    } else {
        printf ("Parallel fast loop, pairs loop, elim loop, and fix loop\n");
        fflush (stdout);
        start_szeit = CCutil_real_zeit ();
        lport = CCutil_snet_listen (ELIM_PORT);
        if (lport == (CC_SPORT *) NULL) {
            fprintf (stderr, "CCutil_snet_listen failed\n");
            rval = 1; goto CLEANUP;
        }
    }

    /* Fast Loop */

    CC_MALLOC (elimlist, 2*ecount, int);  /* For eliminated/fixed edges */

    printf ("\nFast edge-elimination loop ...\n");
    printf ("Total edges to process: %d\n", processcount); fflush (stdout);
    if (incyc) {
        rval = processlist_cycle (ncount, incyc, WORK_ELIM, processcount,
                                  processlist, &tprocesscount, &tprocesslist);
        CCcheck_rval (rval, "processlist_cycle failed");
        processcount = tprocesscount;
        processlist = tprocesslist;
        printf ("Using tour, processlist now has %d edges\n", processcount);
        fflush (stdout);
    }

    if (net_boss == 0) {
        szeit = CCutil_zeit ();
        rval = run_elim_loop (WORK_ELIM, ncount, ecount, elist, processcount,
            processlist, 0, (int *) NULL, dat, rstate,
            &fast_nremain, &fast_remain, &fast_nelim, elimlist, loud,
            (int *) NULL, (int **) NULL, (int **) NULL, (int ***) NULL,
            elim_time_limit, longpath, use_tsp, 1, 0);
        CCcheck_rval (rval, "run_elim_loop failed");
        fast_time = CCutil_zeit () - szeit;
        printf ("Time: %.2f seconds\n", fast_time); fflush (stdout);
    } else {
        rval = elim_boss (WORK_ELIM, ncount, ecount, elist,
            processcount, processlist, 0, (int *) NULL, 
            (int *) NULL, (int **) NULL, dat, &fast_nremain, &fast_remain,
            &fast_nelim, elimlist, (int **) NULL, (int ***) NULL,
            0, 0, longpath, use_tsp,
            elim_time_limit, 0, 1, 1 /* call_fast */,
            0, loud, lport, &round_id, &fast_time, 0, incyc);
        CCcheck_rval (rval, "elim_boss failed");
    }

    sprintf (buf, "%s.fast.remain", probname);
    printf ("Writing edges to %s\n", buf); fflush (stdout);
    rval = CCutil_writeedges (ncount, buf, fast_nremain, fast_remain, dat, 0);
    CCcheck_rval (rval, "CCutil_writeedges failed");
    CC_IFFREE (tprocesslist, int);

    printf ("\nPair-elimination loop ...\n"); fflush (stdout);

    ecount = fast_nremain;
    elist = fast_remain;
    processcount = 0;
    rval = CCelim_build_pairprocesslist (ncount, ecount, elist,
                 (int *) NULL, (int **) NULL, &processcount, &pairprocesslist);
    CCcheck_rval (rval, "CCelim_build_pairprocesslist failed");
    processlist = pairprocesslist;

    printf ("Total pairs to process: %d\n", processcount); fflush (stdout);

    if (net_boss == 0) {
        szeit = CCutil_zeit ();
        rval = run_elim_loop (WORK_PAIR, ncount, ecount, elist, processcount,
            processlist, 0, (int *) NULL, dat, rstate,
            (int *) NULL, (int **) NULL, &pair_nelim, (int *) NULL, loud,
            (int *) NULL, (int **) NULL, &nonpcounts, &nonpairs,
            elim_time_limit, longpath, use_tsp, 0, 0);
        CCcheck_rval (rval, "run_elim_loop failed");
        pair_time = CCutil_zeit () - szeit;
        printf ("Time: %.2f seconds\n", pair_time);
        fflush (stdout);
    } else {
        rval = elim_boss (WORK_PAIR, ncount, ecount, elist,
            processcount, processlist, 0, (int *) NULL, 
            (int *) NULL, (int **) NULL, dat, (int *) NULL, (int **) NULL,
            &pair_nelim, (int *) NULL, &nonpcounts, &nonpairs,
            0, 0, longpath, use_tsp,
            elim_time_limit, 0, 1, 0,
            0, loud, lport, &round_id, &pair_time, 0, incyc);
        CCcheck_rval (rval, "elim_boss failed");
    }

    sprintf (buf, "%s.loop.nonpairs", probname);
    printf ("Writing pairs to %s\n", buf); fflush (stdout);
    rval = CCelim_write_nonpairs (buf, ncount, nonpcounts, nonpairs);
    CCcheck_rval (rval, "CCelim_write_nonpairs failed");

    processcount = ecount;
    processlist = elist;

    printf ("\nEdge-elimination loop ...\n"); fflush (stdout);
    printf ("Total edges to process: %d\n", processcount); fflush (stdout);
    if (incyc) {
        rval = processlist_cycle (ncount, incyc, WORK_ELIM, processcount,
                                  processlist, &tprocesscount, &tprocesslist);
        CCcheck_rval (rval, "processlist_cycle failed");
        processcount = tprocesscount;
        processlist = tprocesslist;
        printf ("Using tour, processlist now has %d edges\n", processcount);
        fflush (stdout);
    }

    if (net_boss == 0) {
        szeit = CCutil_zeit ();
        rval = run_elim_loop (WORK_ELIM, ncount, ecount, elist, processcount,
            processlist, 0, (int *) NULL, dat, rstate,
            &loop_nremain, &loop_remain, &loop_nelim, elimlist, loud,
            nonpcounts, nonpairs, (int **) NULL, (int ***) NULL,
            elim_time_limit, longpath, use_tsp, 0, 0);
        CCcheck_rval (rval, "run_elim_loop failed");
        loop_time = CCutil_zeit () - szeit;
        printf ("Time: %.2f seconds\n", loop_time);
        fflush (stdout);
    } else {
        rval = elim_boss (WORK_ELIM, ncount, ecount, elist,
            processcount, processlist, 0, (int *) NULL, 
            nonpcounts, nonpairs, dat, &loop_nremain, &loop_remain,
            &loop_nelim, elimlist, (int **) NULL, (int ***) NULL,
            0, 0, longpath, use_tsp,
            elim_time_limit, 0, 1, 0,
            0, loud, lport, &round_id, &loop_time, 0, incyc);
        CCcheck_rval (rval, "elim_boss failed");
    }

    sprintf (buf, "%s.loop.remain", probname);
    printf ("Writing edges to %s\n", buf); fflush (stdout);
    rval = CCutil_writeedges (ncount, buf, loop_nremain, loop_remain, dat, 0);
    CCcheck_rval (rval, "CCutil_writeedges failed");
    CC_IFFREE (tprocesslist, int);

    ecount = loop_nremain;
    elist = loop_remain;
    processcount = ecount;
    processlist = elist;
    nfixed = 0;
    printf ("\nFix-edge loop ...\n"); fflush (stdout);
    printf ("Total edges to process: %d\n", processcount); fflush (stdout);
    if (incyc) {
        rval = processlist_cycle (ncount, incyc, WORK_FIX, processcount,
                                  processlist, &tprocesscount, &tprocesslist);
        CCcheck_rval (rval, "processlist_cycle failed");
        processcount = tprocesscount;
        processlist = tprocesslist;
        printf ("Using tour, processlist now has %d edges\n", processcount);
        fflush (stdout);
    }

    if (net_boss == 0) {
        szeit = CCutil_zeit ();
        rval = run_elim_loop (WORK_FIX, ncount, ecount, elist, processcount,
            processlist, 0, (int *) NULL, dat, rstate,
            &nfixed, &fixedlist, (int *) NULL, elimlist, loud,
            nonpcounts, nonpairs, (int **) NULL, (int ***) NULL,
            elim_time_limit, longpath, use_tsp, 0, 0);
        CCcheck_rval (rval, "run_elim_loop failed");
        fix_time = CCutil_zeit () - szeit;
        printf ("Time: %.2f seconds\n", fix_time);
        fflush (stdout);
    } else {
        rval = elim_boss (WORK_FIX, ncount, ecount, elist,
            processcount, processlist, 0, (int *) NULL, 
            nonpcounts, nonpairs, dat, &nfixed, &fixedlist,
            (int *) NULL, elimlist, (int **) NULL, (int ***) NULL,
            0, 0, longpath, use_tsp,
            elim_time_limit, 0, 1, 0,
            0, loud, lport, &round_id, &fix_time, 0, incyc);
        CCcheck_rval (rval, "elim_boss failed");
    }

    sprintf (buf, "%s.loop.fix", probname);
    printf ("Writing fixed edges to %s\n", buf); fflush (stdout);
    rval = CCutil_writeedges (ncount, buf, nfixed, fixedlist, dat, 0);
    CCcheck_rval (rval, "CCutil_writeedges failed");
    CC_IFFREE (tprocesslist, int);

    printf ("\nSummary\n");
    printf ("Fast Loop: %d eliminated, %d remaining, %.2f seconds\n",
                fast_nelim, fast_nremain, fast_time);
    printf ("Pair Loop: %d nonpairs, %.2f seconds\n", pair_nelim, pair_time);
    printf ("Edge Loop: %d eliminated, %d remaining, %.2f seconds\n",
                loop_nelim, loop_nremain, loop_time);
    printf ("Fix Loop:  %d fixed, %.2f seconds\n", nfixed, fix_time);
    fflush (stdout);

    if (net_boss == 0) {
        printf ("Total Time: %.2f seconds\n", CCutil_zeit () - start_szeit);
    } else {
        printf ("Total Wall Clock: %.0f seconds\n",
                                         CCutil_real_zeit () - start_szeit);
    }
    fflush (stdout);

CLEANUP:
    CC_IFFREE (fast_remain, int);
    CC_IFFREE (loop_remain, int);
    CC_IFFREE (elimlist, int);
    CC_IFFREE (fixedlist, int);
    CC_IFFREE (tprocesslist, int);
    CC_IFFREE (pairprocesslist, int);
    free_pair_info (ncount, &nonpcounts, &nonpairs);
    if (lport) CCutil_snet_unlisten (lport);
    return rval;
}

/* run_elim is a single pass of the elimination code, using the parameters */
/* specified on the command line                                           */

static int run_elim (int gotype, int ncount, int ecount, int *elist,
        int processcount, int *processlist, int fixcount, int *fixlist,
        CCdatagroup *dat, CCrandstate *rstate, int level_count, int wtype,
        int *remaincount, int **remainlist, int *pelimcount, int *elimlist,
        int loud, int *nonpcounts, int **nonpairs, int **outpcounts,
        int ***outpairs, double timelimit, int longpath, int use_tsp,
        int max_neighborhood, int call_fast, int save_tree)
{
    int rval = 0, i, a, b, c = 0, yesno, elimcount = 0;
    int *donepairlist = (int *) NULL;
    int *pcounts = (int *) NULL, **pairs = (int **) NULL, processed = 0;
    CCelim_graph G, fixG, *F = &fixG;
    CCelim_distobj D;
    CCdatagroup *euclid_dat = (CCdatagroup *) NULL;
    CCkdtree *kt = (CCkdtree *) NULL, *euclid_kt = (CCkdtree *) NULL;
    CCelim_httree *ht = (CCelim_httree *) NULL;
    double ezeit = 0.0, szeit, roundzeit;
    char s = '?';

    if (gotype == WORK_PAIR) { CC_MALLOC (donepairlist, 3*processcount, int); }
    switch (gotype) {
        case WORK_PAIR: s = 'P'; break;
        case WORK_FIX:  s = 'F'; break;
        case WORK_ELIM: s = 'E';
    }

    rval = CCelim_alloc_elim (ncount, ecount, elist, fixcount, fixlist, dat,
               &G, F, &D, &kt, &euclid_dat, &euclid_kt, rstate);
    CCcheck_rval (rval, "CCelim_alloc_elim failed");

    if (save_tree) {
        CC_MALLOC (ht, 1, CCelim_httree);
        CCelim_init_httree (ht);
        ht->max_count = HT_MAX_COUNT;
        ht->max_depth = HT_MAX_DEPTH;
    }

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

    szeit = roundzeit =  CCutil_zeit ();

    printf ("Process:    %d\n", processcount); fflush (stdout);
    for (i = 0; i < processcount; i++) {
        if (gotype == WORK_PAIR) {
            a=processlist[3*i]; b=processlist[3*i+1]; c=processlist[3*i+2];
            rval = CCelim_run_elim_pair (b, a, c, &G, F, &D,
                   rstate, kt, euclid_dat, euclid_kt, &yesno, &ezeit,
                   level_count, pcounts, pairs, timelimit, longpath, use_tsp,
                   max_neighborhood, ht);
            CCcheck_rval (rval, "CCelim_run_elim_pair failed");
        } else {
            a = processlist[2*i]; b = processlist[2*i+1];
            if (gotype == WORK_FIX) {
                if (CCelim_edge_in_graph (a, b, F)) continue;
                rval = CCelim_run_fix_edge (a, b, &G, F, &D, rstate, kt,
                   euclid_dat, euclid_kt, &yesno, &ezeit, level_count, pcounts,
                   pairs, timelimit, use_tsp, max_neighborhood, ht);
                CCcheck_rval (rval, "CCelim_run_fix_edge failed");

            } else {
                if (call_fast) {
                    rval = CCelim_run_fast_elim_edge (a, b, &G, &D, rstate, kt,
                       euclid_dat, euclid_kt, level_count, max_neighborhood,
                       &yesno, ht);
                    CCcheck_rval (rval, "CCelim_run_fast_elim_edge failed");
                } else {
                    rval = CCelim_run_elim_edge (a, b, &G, F, &D, wtype,
                       rstate, &yesno, &ezeit, level_count, pcounts, pairs,
                       kt, euclid_dat, euclid_kt, timelimit, longpath, use_tsp,
                       max_neighborhood, ht);
                    CCcheck_rval (rval, "CCelim_run_elim_edge failed");
                }
            }
        }
        if (yesno == 1) { 
            if (ht) {
/*
                printf ("Hamilton-Tutte Tree: %d %d\n", ht->count, ht->depth);
                fflush (stdout);
                rval = CCelim_print_httree (ht, stdout, 0);
                CCcheck_rval (rval, "CCelim_print_httree failed");
*/
                if (elimcount == 0) {
                    /* create new file, overwriting if it exists */
                    rval = CCelim_write_httree (ht, htfname, 0, 0);
                } else {
                    /* append to existing file */
                    rval = CCelim_write_httree (ht, htfname, 1, 0);
                }
                CCcheck_rval (rval, "CCelim_write_httree failed");

                CCelim_free_htnode_children (ht->root);
                CC_IFFREE (ht->root, CCelim_htnode);
            }

            if (loud == 1) {
                printf ("%c", s); fflush (stdout);
            } else if (loud > 1) {
                if (gotype == WORK_PAIR) {
                    printf ("Z %d %d %d\n", a, b, c);
                } else {
                    printf ("Z %d %d %d\n", a, b, CCelim_dist (a, b, &D));
                }
                fflush (stdout);
            }
            if (gotype == WORK_PAIR) {
                donepairlist[3*elimcount]   = a;
                donepairlist[3*elimcount+1] = b;
                donepairlist[3*elimcount+2] = c;
            } else {
                elimlist[2*elimcount]   = a;
                elimlist[2*elimcount+1] = b;
                if (gotype == WORK_ELIM) {
                    CCelim_delete_edge (a, b, &G);
                } else {
                    CCelim_add_edge (a, b, CCelim_dist (a, b, &D), F);
                }
            }
            elimcount++;
        } else {
            if (ht) {
                CCelim_free_htnode_children (ht->root);
                CC_IFFREE (ht->root, CCelim_htnode);
            }
        }

        if (loud > 1) { printf ("%.2f\n", ezeit); fflush (stdout); }
        if (processed++ % 1000 == 999) {
            if (loud == 1) printf ("\n");
            printf ("%d %d %.2f (%.2f seconds, %.2f total)\n",
                  elimcount, processed, (double) elimcount / (double) processed,
                  CCutil_zeit () - roundzeit, CCutil_zeit () - szeit);
            fflush (stdout);
            roundzeit = CCutil_zeit ();
        }
    }

    if (gotype == WORK_PAIR) {
        printf ("\nEliminated Pairs: %d\n", elimcount);
        rval = CCelim_pairlist_to_nonpairs (&G, elimcount, donepairlist,
                            nonpcounts, nonpairs, outpcounts, outpairs);
        CCcheck_rval (rval, "pairlist_to_nonpairs failed");
    } else if (gotype == WORK_FIX) {
        printf ("\nFixed: %d\n", elimcount);
        rval = CCelim_getedges_graph (F, remaincount, remainlist);
        CCcheck_rval (rval, "CCelim_getedges_graph failed");
        printf ("Total: %d\n", *remaincount);
    } else {
        printf ("\nEliminated: %d\n", elimcount);
        rval = CCelim_getedges_graph (&G, remaincount, remainlist);
        CCcheck_rval (rval, "CCelim_getedges_graph failed");
        printf ("Remaining:  %d\n", *remaincount);
    }
    fflush (stdout);
    if (pelimcount) *pelimcount = elimcount;

CLEANUP:
    CCelim_free_elim (&G, F, &D, &kt, &euclid_dat, &euclid_kt);
    free_pair_info (ncount, &pcounts, &pairs);
    CC_IFFREE (donepairlist, int);
    CC_IFFREE (ht, CCelim_httree);
    return rval;
}

/* run_elim_loop calls repeatedly the elimination code, with parameters     */
/* gradually increased (to more effective, but more time-consuming values). */
/* the settings of the parameters are increased each time the previous      */
/* round of elimination failed to eliminate a fixed fraction (specifed by   */
/* the variable target) of the existing edges.                              */

/* the levels of the parameters are hard-coded into the function. to tune   */
/* the code for your problems, you can adjust the values of wtype, level,   */
/* and max_neighborhood: below for "SETTINGS FOR LOOPS" comment             */

static int run_elim_loop (int gotype, int ncount, int ecount, int *elist,
        int processcount, int *processlist, int fixcount, int *fixlist,
        CCdatagroup *dat, CCrandstate *rstate,
        int *remaincount, int **remainlist, int *pelimcount, int *elimlist,
        int loud, int *nonpcounts, int **nonpairs, int **outpcounts,
        int ***outpairs, double timelimit, int longpath, int use_tsp,
        int call_fast, int single_fast)
{
    int rval = 0, i, a, b, c, yesno, elimcount = 0;
    int processed = 0, level = -1, maxlevel = 1, setlevel = 1, round = 1;
    int round_elimcount, round_total, *donepairlist = (int *) NULL;
    int *hit = (int *) NULL, *pcounts = (int *) NULL, **pairs = (int **) NULL;
    int *round_nonpcounts = (int *) NULL, **round_nonpairs = (int **) NULL;
    int *mynonpcounts = (int *) NULL, **mynonpairs = (int **) NULL;
    int wtype = CCelim_CD_EDGE, max_neighborhood = 5;
    CCelim_graph G, F;
    CCelim_distobj D;
    CCdatagroup *euclid_dat = (CCdatagroup *) NULL;
    CCkdtree *kt = (CCkdtree *) NULL, *euclid_kt = (CCkdtree *) NULL;
    double ezeit, szeit, roundzeit, zeit1000;
    double target = (gotype == WORK_PAIR ? 0.25 : 0.05);
    char t = 'E';

    rval = CCelim_alloc_elim (ncount, ecount, elist, fixcount, fixlist, dat,
               &G, &F, &D, &kt, &euclid_dat, &euclid_kt, rstate);
    CCcheck_rval (rval, "CCelim_alloc_elim failed");

    /* in the loop, hit[i] will mark if edge i has already been eliminated */

    CC_MALLOC (hit, processcount, int);
    for (i = 0; i < processcount; i++) hit[i] = 0;
    if (gotype == WORK_PAIR) CC_MALLOC (donepairlist, 3*processcount, int);

    switch (gotype) {
        case WORK_PAIR: t = 'P'; break;
        case WORK_FIX:  t = 'F'; break;
        case WORK_ELIM: t = 'E';
    }

    switch (gotype) {
        case WORK_PAIR: maxlevel = 3; break;
        case WORK_FIX:  maxlevel = 1; break;
        case WORK_ELIM: 
            if (single_fast) maxlevel = single_fast;
            else             maxlevel = (call_fast ? 8 : 3);
            break;
    }

    szeit = roundzeit = zeit1000 = CCutil_zeit ();

    while (setlevel <= maxlevel && processcount - elimcount > 0) {
        if (round_nonpairs) {
            mynonpcounts = round_nonpcounts; mynonpairs = round_nonpairs;
        } else {
            mynonpcounts = nonpcounts; mynonpairs = nonpairs;
        }

        if (mynonpairs) {
            rval = CCelim_nonpairs_to_pairs (&G, mynonpcounts, mynonpairs, 
                                             &pcounts, &pairs);
            CCcheck_rval (rval, "CCelim_nonpairs_to_pairs failed");
            a = b = 0;
            for (i = 0; i < ncount; i++) {
                a += ((G.nodelist[i].deg * (G.nodelist[i].deg - 1)) / 2);
                b += pcounts[i];
            }
            printf ("Pairs:      %d (%.1f%%)\n", b,((double)b/(double)a)*100.0);
            fflush (stdout);
        }
        free_pair_info (ncount, &round_nonpcounts, &round_nonpairs);

        /* SETTINGS FOR LOOPS */

        if (gotype == WORK_FIX) {
            level = 3; max_neighborhood = 25;
        } else if (gotype == WORK_PAIR) {
            switch (setlevel) {
                case 1: level = -1; max_neighborhood = 5;  break;
                case 2: level = 0;  max_neighborhood = 10; break;
                case 3: level = 1;  max_neighborhood = 25; break;
            }
        } else if (call_fast) {
            switch (setlevel) {
                case 1: level = -1; max_neighborhood = 5;  break;
                case 2: level = -1; max_neighborhood = 10; break;
                case 3: level = 0;  max_neighborhood = 5;  break;
                case 4: level = 0;  max_neighborhood = 10; break;
                case 5: level = 1;  max_neighborhood = 10; break;
                case 6: level = 1;  max_neighborhood = 15; break;
                case 7: level = 2;  max_neighborhood = 25; break;
                case 8: level = 3;  max_neighborhood = 25; break;
            }
        } else {
            switch (setlevel) {
                case 1: level=2; max_neighborhood=25; wtype=CCelim_CD_NONEDGE;
                        break;
                case 2: level=3; max_neighborhood=25; wtype=CCelim_CD_EDGE;
                        break;
                case 3: level=3; max_neighborhood=50; wtype=CCelim_CD_EDGE; 
                        break;
            }
        }

        printf ("STARTING ROUND %d, LEVEL %d (z = %d, n = %d, w = %d)\n",
                round, setlevel, level+3, max_neighborhood, wtype);
        fflush (stdout);

        round_elimcount = 0;
        round_total = processcount - elimcount;
        processed = 0;

        for (i = 0; i < processcount; i++) {
            if (hit[i]) continue;
            if (gotype == WORK_PAIR) {
                a = processlist[3*i];
                b = processlist[3*i+1];
                c = processlist[3*i+2];
                rval = CCelim_run_elim_pair (b, a, c, &G, &F, &D,
                   rstate, kt, euclid_dat, euclid_kt, &yesno, &ezeit, level,
                   pcounts, pairs, timelimit, longpath, use_tsp,
                   max_neighborhood, (CCelim_httree *) NULL);
                CCcheck_rval (rval, "CCelim_run_elim_pair failed");
                if (yesno == 1) {
                    donepairlist[3*elimcount]   = a;
                    donepairlist[3*elimcount+1] = b;
                    donepairlist[3*elimcount+2] = c;
                }
            } else if (gotype == WORK_FIX) {
                a = processlist[2*i];
                b = processlist[2*i+1];
                if (!CCelim_edge_in_graph (a, b, &G) ||
                     CCelim_edge_in_graph (a, b, &F)) continue;
                rval = CCelim_run_fix_edge (a, b, &G, &F, &D, rstate, kt,
                   euclid_dat, euclid_kt, &yesno, &ezeit, level, pcounts,
                   pairs, timelimit, use_tsp, max_neighborhood,
                   (CCelim_httree *) NULL);
                CCcheck_rval (rval, "CCelim_run_fix_edge failed");
                if (yesno == 1) {
                    elimlist[2*elimcount] = a;
                    elimlist[2*elimcount+1] = b;
                    CCelim_add_edge (a, b, CCelim_dist (a, b, &D), &F);
                }
            } else {
                a = processlist[2*i];
                b = processlist[2*i+1];
                if (!CCelim_edge_in_graph (a, b, &G)) continue;
                if (call_fast) {
                    rval = CCelim_run_fast_elim_edge (a, b, &G, &D, rstate, kt,
                       euclid_dat, euclid_kt, level, max_neighborhood, &yesno,
                       (CCelim_httree *) NULL);
                    CCcheck_rval (rval, "CCelim_run_fast_elim_edge failed");
                } else {
                    rval = CCelim_run_elim_edge (a, b, &G, &F, &D, wtype,
                       rstate, &yesno, &ezeit, level, pcounts, pairs,
                       kt, euclid_dat, euclid_kt, timelimit, longpath, use_tsp,
                       max_neighborhood, (CCelim_httree *) NULL);
                    CCcheck_rval (rval, "CCelim_run_elim_edge failed");
                }
                if (yesno == 1) {
                    elimlist[2*elimcount] = a;
                    elimlist[2*elimcount+1] = b;
                    CCelim_delete_edge (a, b, &G);
                }
            }

            if (yesno == 1) { 
                hit[i] = 1;
                if (loud) { printf ("%c", t); fflush (stdout); }
                elimcount++;
                round_elimcount++;
            }

            if (processed++ % 1000 == 999) {
                if (loud) printf ("\n");
                printf ("%d %d/%d %.2f %d total (%.2f seconds, %.2f total)\n",
                  round_elimcount, processed, round_total,
                  (double) round_elimcount / (double) processed, elimcount,
                  CCutil_zeit () - zeit1000, CCutil_zeit () - szeit);
                fflush (stdout);
                zeit1000 = CCutil_zeit ();
            }
        }
        if (loud) printf ("\n");
        printf ("Round %d, level %d: %d new, %d total, %.2f seconds\n", round,
                    setlevel, round_elimcount, elimcount,
                    CCutil_zeit () - roundzeit);
        fflush (stdout);

        if (gotype == WORK_PAIR) {
            rval = CCelim_pairlist_to_nonpairs (&G, elimcount, donepairlist,
                      nonpcounts, nonpairs, &round_nonpcounts, &round_nonpairs);
            CCcheck_rval (rval, "pairlist_to_nonpairs failed");
        }
        free_pair_info (ncount, &pcounts, &pairs);

        roundzeit = CCutil_zeit ();
        round++;
        if ((double) round_elimcount  / (double) round_total < target) {
            setlevel++;
        }
    }

    if (gotype == WORK_PAIR) {
        printf ("\nEliminated Pairs: %d\n", elimcount);
        rval = CCelim_pairlist_to_nonpairs (&G, elimcount, donepairlist,
                        nonpcounts, nonpairs, outpcounts, outpairs);
        CCcheck_rval (rval, "pairlist_to_nonpairs failed");
        {
            /* this call is just to output remaining count */
            rval = CCelim_nonpairs_to_pairs (&G, *outpcounts, *outpairs, 
                                             &pcounts, &pairs);
            CCcheck_rval (rval, "CCelim_nonpairs_to_pairs failed");
            a = b = 0;
            for (i = 0; i < ncount; i++) {
                a += ((G.nodelist[i].deg * (G.nodelist[i].deg - 1)) / 2);
                b += pcounts[i];
            }
            printf ("Remaining Pairs:  %d (%.1f%%)\n",
                                       b,((double)b/(double)a)*100.0);
            fflush (stdout);
            free_pair_info (ncount, &pcounts, &pairs);
        }
    } else if (gotype == WORK_FIX) {
        printf ("\nFixed: %d\n", elimcount);
        rval = CCelim_getedges_graph (&F, remaincount, remainlist);
        CCcheck_rval (rval, "CCelim_getedges_graph failed");
        printf ("Total: %d\n", *remaincount);
    } else {
        printf ("\nEliminated: %d\n", elimcount);
        rval = CCelim_getedges_graph (&G, remaincount, remainlist);
        CCcheck_rval (rval, "CCelim_getedges_graph failed");
        printf ("Remaining:  %d\n", *remaincount);
    }
    fflush (stdout);
    if (pelimcount) *pelimcount = elimcount;

CLEANUP:
    CCelim_free_elim (&G, &F, &D, &kt, &euclid_dat, &euclid_kt);
    CC_IFFREE (donepairlist, int);
    free_pair_info (ncount, &pcounts, &pairs);
    free_pair_info (ncount, &round_nonpcounts, &round_nonpairs);
    return rval;
}


#define MAX_GRUNT_GROUP 10000  /* maximun for the grunt_group variable */

/* elim_boss executes the boss code for a parallel run, passing input graph */
/* and settings to workers, and handing out tasks. if use_loop is set to 1, */
/* then the boss will use repeated loops, as in the run_elim_loop function  */

static int elim_boss (int gotype, int ncount, int ecount, int *elist,
        int processcount, int *processlist, int fixcount, int *fixlist,
        int *nonpcounts, int **nonpairs, CCdatagroup *dat,
        int *remaincount, int **remainlist,
        int *pelimcount, int *elimlist, int **outpcounts, int ***outpairs,
        int level_count, int wtype, int longpath, int use_tsp,
        double timelimit, int max_neighborhood, int use_loop, int call_fast,
        int single_fast, int loud,
        CC_SPORT *lport, int *round_id, double *boss_wall_time, int save_tree,
        int *incyc)
{
    int rval = 0, ngrunt = 0, i, j, id, yesno, next, counter = 0, done = 0;
    int a, b, c, receivecount = 0, grunt_group, level = -1;
    int elimcount = 0, extra = 0, setlevel = 1, round_total, maxlevel = 1;
    int *hits = (int *) NULL, *perm = (int *) NULL, *fini = (int *) NULL;
    int sendcount = 0, sendlist[MAX_GRUNT_GROUP], start_extra = 0;
    int gruntround_id = -1, round = 1, round_ecount = 0, round_elimcount = 0;
    int *round_elist = (int *) NULL, *round_elen = (int *) NULL;
    int *round_nonpcounts = (int *) NULL, **round_nonpairs = (int **) NULL;
    int *mynonpcounts = (int *) NULL, **mynonpairs = (int **) NULL;
    int *donepairlist = (int *) NULL, *newplist = (int *) NULL, newpcount = 0;
    int grunt_gotype, got_tree = 0;
    double target = (gotype == WORK_PAIR ? 0.25 : 0.05);
    double wallzeit, round_wallzeit;
    char request, t = 'E';
    CCelim_graph G, F;
    CCelim_distobj D;
    CC_SFILE *s = (CC_SFILE *) NULL;
    CCelim_httree *ht = (CCelim_httree *) NULL;

    if (boss_wall_time) *boss_wall_time = 0.0;
    wallzeit = CCutil_real_zeit ();

    CCelim_init_graph (&G);
    CCelim_init_graph (&F);
    CCelim_init_distobj (&D);
    rval = CCelim_build_distobj (&D, ncount, dat);
    CCcheck_rval (rval, "CCelim_build_distobj failed");
    rval = CCelim_build_graph (&G, ncount, ecount, elist, (int *) NULL, &D);
    CCcheck_rval (rval, "CCelim_build_graph failed");
    rval = CCelim_build_graph (&F, ncount, fixcount, fixlist, (int *) NULL, &D);
    CCcheck_rval (rval, "CCelim_build_graph failed");

    if (gotype == WORK_FIX && fixcount > 0) {
        CC_MALLOC (newplist, 2*processcount, int);
        for (i = 0; i < processcount; i++) {
            a = processlist[2*i]; b = processlist[2*i+1];
            if (!CCelim_edge_in_graph (a, b, &F)) {
                newplist[2*newpcount] = a; newplist[2*newpcount+1] = b;
                newpcount++;
            }
        }
        if (newpcount < processcount) {
            printf ("Reduced:    %d\n", newpcount);
            processlist = newplist;
            processcount = newpcount;
        }
    }
 
    CC_MALLOC (perm, ncount, int);
    for (i = 0; i < ncount; i++) perm[i] = 0;
    CC_MALLOC (hits, processcount, int);
    for (i = 0; i < processcount; i++) hits[i] = 0;
    CC_MALLOC (fini, processcount, int);
    if (gotype == WORK_PAIR) CC_MALLOC (donepairlist, 3*processcount, int);

    switch (gotype) {
        case WORK_PAIR: t = 'P'; break;
        case WORK_FIX:  t = 'F'; break;
        case WORK_ELIM: t = 'E';
    }

    switch (gotype) {
        case WORK_PAIR: maxlevel = 3; break;
        case WORK_FIX:  maxlevel = 1; break;
        case WORK_ELIM:
            if (single_fast) maxlevel = single_fast;
            else             maxlevel = (call_fast ? 8 : 5);
            break;
    }

    printf ("BEGINNING ELIMINATION NET PROCESSING: %d %s\n\n", processcount,
           (gotype == WORK_PAIR ? "pairs" : "edges"));
    fflush (stdout);

    if (use_loop == 0) setlevel = maxlevel = 1;

    while (setlevel <= maxlevel && processcount - elimcount > 0) {
        if (use_loop != 0) {
            CC_IFFREE (round_elist, int);
            CC_IFFREE (round_elen, int);
        }
        /* grunt_group = (level_count <= 2 ? 1000 : 100); */
        grunt_group = (level_count <= 2 ? 100 : 10);
        /* grunt_group = 1; */

        rval = CCelim_getedges_graph (&G, &round_ecount, &round_elist);
        CCcheck_rval (rval, "CCelim_getedges_graph failed");
        CC_MALLOC (round_elen, round_ecount, int);
        for (i = 0; i < round_ecount; i++) {
            round_elen[i] = CCutil_dat_edgelen (round_elist[2*i],
                                                round_elist[2*i+1], dat);
        }

        if (gotype == WORK_PAIR && elimcount > 0) {
            free_pair_info (ncount, &round_nonpcounts, &round_nonpairs);
            rval = CCelim_pairlist_to_nonpairs (&G, elimcount, donepairlist,
                      nonpcounts, nonpairs, &round_nonpcounts, &round_nonpairs);
            CCcheck_rval (rval, "pairlist_to_nonpairs failed");
            mynonpcounts = round_nonpcounts; mynonpairs = round_nonpairs;
        } else {
            mynonpcounts = nonpcounts; mynonpairs = nonpairs;
        }

        round_total = processcount - elimcount;
        printf ("Remaining %s to process: %d\n",
                        (gotype==WORK_PAIR ? "pairs" : "edges"), round_total);
        fflush (stdout);

        extra = counter = done = round_elimcount = 0;
        round_wallzeit = CCutil_real_zeit ();
        for (i = 0; i < processcount; i++) fini[i] = 0;

        if (use_loop == 0) {
            level = level_count;
        } else if (gotype == WORK_FIX) {
            level = 3; max_neighborhood = 25;
        } else if (gotype == WORK_PAIR) {
            switch (setlevel) {
                case 1: level = -1; max_neighborhood = 5;  break;
                case 2: level = 0;  max_neighborhood = 10; break;
                case 3: level = 1;  max_neighborhood = 25; break;
            }
        } else if (call_fast) {
            switch (setlevel) {
                case 1: level = -1; max_neighborhood = 5;  break;
                case 2: level = -1; max_neighborhood = 10; break;
                case 3: level = 0;  max_neighborhood = 5;  break;
                case 4: level = 0;  max_neighborhood = 10; break;
                case 5: level = 1;  max_neighborhood = 10; break;
                case 6: level = 1;  max_neighborhood = 15; break;
                case 7: level = 2;  max_neighborhood = 25; break;
                case 8: level = 3;  max_neighborhood = 25; break;
            }
        } else {
            switch (setlevel) {
                case 1: level=1; max_neighborhood=25; wtype=CCelim_CD_NONEDGE;
                        break;
                case 2: level=1; max_neighborhood=50; wtype=CCelim_CD_NONEDGE;
                        break;
                case 3: level=2; max_neighborhood=25; wtype=CCelim_CD_NONEDGE;
                        break;
                case 4: level=3; max_neighborhood=25; wtype=CCelim_CD_EDGE;
                        break;
                case 5: level=3; max_neighborhood=50; wtype=CCelim_CD_EDGE;
                        break;
            }
        }

        if (use_loop != 0) {
            printf ("STARTING ROUND %d, LEVEL %d (z = %d, n = %d, w = %d)\n",
                round, setlevel, level+3, max_neighborhood, wtype);
            fflush (stdout);
        }

        while (done < round_total) {
            s = CCutil_snet_receive (lport);
            if (!s) {
                fprintf (stderr, "CCutil_snet_receive failed, ignoring\n");
                continue;
            }

            rval = CCutil_sread_char (s, &request);
            CCcheck_rval (rval, "CCutil_sread_char failed (request)");

            switch (request) {

            case ELIM_HELLO:
            printf ("Welcome grunt %d\n", ngrunt++); fflush (stdout);
            rval = CCutil_writemaster (s, ncount, dat, perm);
            CCcheck_rval (rval, "CCutil_writemaster failed");
            if (incyc) {
                rval = CCutil_swrite_int (s, 1);
                CCcheck_rval (rval, "CCutil_swrite_int failed (incyc)");
                for (i = 0; i < ncount; i++) {
                    rval = CCutil_swrite_int (s, incyc[i]);
                    CCcheck_rval (rval, "CCutil_swrite_int failed (incyc)");
                }
            } else {
                rval = CCutil_swrite_int (s, 0);
                CCcheck_rval (rval, "CCutil_swrite_int failed (incyc)");
            }
            break;

            case ELIM_WORK:
            rval = CCutil_sread_int (s, &gruntround_id);
            CCcheck_rval (rval, "CCutil_sread_int failed (gruntround_id)");
            if (gruntround_id == *round_id) {
                /* grunt is on the current ID */
                rval = CCutil_swrite_int (s, 0);
                CCcheck_rval (rval, "CCutil_swrite_int failed (newgraph)");
            } else {
                /* grunt is on an old ID, send current graph and ID */
                rval = CCutil_swrite_int (s, 1);
                CCcheck_rval (rval, "CCutil_swrite_int failed (newgraph)");
                rval = CCutil_swrite_int (s, *round_id);
                CCcheck_rval (rval, "CCutil_swrite_int failed (round_id)");
                rval = CCutil_swrite_int (s, gotype);
                CCcheck_rval (rval, "CCutil_swrite_int failed (gotype)");
                rval = CCutil_swrite_int (s, ncount);
                CCcheck_rval (rval, "CCutil_swrite_int failed (ncount)");
                rval = CCutil_swrite_int (s, round_ecount);
                CCcheck_rval (rval, "CCutil_swrite_int failed (ecount)");
                for (i = 0; i < round_ecount; i++) {
                    rval = CCutil_swrite_int (s, round_elist[2*i]);
                    CCcheck_rval (rval, "CCutil_swrite_int failed (elist)");
                    rval = CCutil_swrite_int (s, round_elist[2*i+1]);
                    CCcheck_rval (rval, "CCutil_swrite_int failed (elist)");
                }
                for (i = 0; i < round_ecount; i++) {
                    rval = CCutil_swrite_int (s, round_elen[i]);
                    CCcheck_rval (rval, "CCutil_swrite_int failed (elen)");
                }
                rval = CCutil_swrite_int (s, fixcount);
                CCcheck_rval (rval, "CCutil_swrite_int failed (fixcount)");
                for (i = 0; i < fixcount; i++) {
                    rval = CCutil_swrite_int (s, fixlist[2*i]);
                    CCcheck_rval (rval, "swrite_int failed (fixlist[0])");
                    rval = CCutil_swrite_int (s, fixlist[2*i+1]);
                    CCcheck_rval (rval, "swrite_int failed (fixlist[1])");
                }

                rval = CCutil_swrite_double (s, timelimit);
                CCcheck_rval (rval, "CCutil_swrite_double failed (timelimit)");
                rval = CCutil_swrite_int (s, level);
                CCcheck_rval (rval, "CCutil_swrite_int failed (level)");
                rval = CCutil_swrite_int (s, wtype);
                CCcheck_rval (rval, "CCutil_swrite_int failed (wtype)");
                rval = CCutil_swrite_int (s, longpath);
                CCcheck_rval (rval, "CCutil_swrite_int failed (longpath)");
                rval = CCutil_swrite_int (s, use_tsp);
                CCcheck_rval (rval, "CCutil_swrite_int failed (use_tsp)");
                rval = CCutil_swrite_int (s, max_neighborhood);
                CCcheck_rval (rval, "CCutil_swrite_int failed (max)");
                rval = CCutil_swrite_int (s, call_fast);
                CCcheck_rval (rval, "CCutil_swrite_int failed (call_fast)");
                rval = CCutil_swrite_int (s, save_tree);
                CCcheck_rval (rval, "CCutil_swrite_int failed (save_tree)");

                if (mynonpairs) {
                    rval = CCutil_swrite_int (s, 1);
                    CCcheck_rval (rval, "CCutil_swrite_int failed (1)");
                    for (i = 0; i < ncount; i++) {
                        rval = CCutil_swrite_int (s, mynonpcounts[i]);
                        CCcheck_rval (rval, "CCutil_swrite_int failed (non)");
                        for (j = 0; j < mynonpcounts[i]; j++) {
                            rval=CCutil_swrite_int(s,mynonpairs[i][2*j]);
                            CCcheck_rval (rval, "swrite_int failed (pairs)");
                            rval=CCutil_swrite_int(s,mynonpairs[i][2*j+1]);
                            CCcheck_rval (rval, "swrite_int failed (pairs)");
                        }
                    }
                } else {
                    rval = CCutil_swrite_int (s, 0);
                    CCcheck_rval (rval, "CCutil_swrite_int failed (0)");
                }
            }

            sendcount = 0;
            if (counter < processcount) {
                do {
                    if (hits[counter] == 0 && fini[counter] == 0) {
                        sendlist[sendcount++] = counter;
                    }
                    counter++;
                } while (sendcount < grunt_group && counter < processcount);
            }

            if (sendcount == 0) {  /* send any outstanding old work */
                start_extra = extra;
                do {
                    while (hits[extra] == 1 || fini[extra] == 1) {
                       if (++extra >= processcount) extra = 0;
                    }
                    sendlist[sendcount++] = extra;
                    if (++extra >= processcount) extra = 0;
                } while (sendcount < grunt_group && extra != start_extra);
            }

            rval = CCutil_swrite_int (s, sendcount);
            CCcheck_rval (rval, "CCutil_swrite_int failed");
            for (i = 0; i < sendcount; i++) {
                next = sendlist[i];
                rval = CCutil_swrite_int (s, next);
                CCcheck_rval (rval, "CCutil_swrite_int failed (next)");
                if (gotype == WORK_PAIR) {
                    rval = CCutil_swrite_int (s, processlist[3*next]);
                    CCcheck_rval (rval, "CCutil_swrite_int failed (a)");
                    rval = CCutil_swrite_int (s, processlist[3*next+1]);
                    CCcheck_rval (rval, "CCutil_swrite_int failed (b)");
                    rval = CCutil_swrite_int (s, processlist[3*next+2]);
                    CCcheck_rval (rval, "CCutil_swrite_int failed (c)");
                } else {
                    rval = CCutil_swrite_int (s, processlist[2*next]);
                    CCcheck_rval (rval, "CCutil_swrite_int failed (a)");
                    rval = CCutil_swrite_int (s, processlist[2*next+1]);
                    CCcheck_rval (rval, "CCutil_swrite_int failed (b)");
                }
            }
            break;

            case ELIM_DONE:
            rval = CCutil_sread_int (s, &gruntround_id);
            CCcheck_rval (rval, "CCutil_sread_int failed (gruntround_id)");
            rval = CCutil_sread_int (s, &grunt_gotype);
            CCcheck_rval (rval, "CCutil_sread_int failed (gruntround_id)");
            rval = CCutil_sread_int (s, &receivecount);
            CCcheck_rval (rval, "CCutil_sread_int failed (receivecount)");
            for (i = 0; i < receivecount; i++) {
                rval = CCutil_sread_int (s, &id);
                CCcheck_rval (rval, "CCutil_sread_int failed (id)");
                rval = CCutil_sread_int (s, &a);
                CCcheck_rval (rval, "CCutil_sread_int failed (a)");
                rval = CCutil_sread_int (s, &b);
                CCcheck_rval (rval, "CCutil_sread_int failed (b)");
                if (grunt_gotype == WORK_PAIR) {
                    rval = CCutil_sread_int (s, &c);
                    CCcheck_rval (rval, "CCutil_sread_int failed (c)");
                }
                rval = CCutil_sread_int (s, &yesno);
                CCcheck_rval (rval, "CCutil_sread_int failed (yesno)");
                if (yesno && save_tree) {
                    rval = CCutil_sread_int (s, &got_tree);
                    CCcheck_rval (rval, "CCutil_sread_int failed (got_tree)");
                    if (got_tree) {
                        rval = CCelim_read_httree (&ht, (FILE *) NULL, s);
                        CCcheck_rval (rval, "CCelim_read_httree failed");
                        if (gruntround_id == *round_id && fini[id] == 0) {
                            if (elimcount == 0) {
                                /* create new file, overwriting if it exists */
                                rval = CCelim_write_httree (ht, htfname, 0, 0);
                            } else {
                                /* append to existing file */
                                rval = CCelim_write_httree (ht, htfname, 1, 0);
                            }
                            CCcheck_rval (rval, "CCelim_write_httree failed");
                        }
                        CCelim_free_httree (ht);
                        CC_IFFREE (ht, CCelim_httree);
                    } else { 
                        printf ("WARNING: did not receive the HT tree\n");
                        fflush (stdout);
                    }
                }
                if (gruntround_id == *round_id && fini[id] == 0) {
                    if (grunt_gotype != gotype) {
                        printf ("gotype out of sync!\n"); exit (1);
                    }
                    if (hits[id] == 1) {
                        printf ("Already finished!\n"); exit (1);
                    }
                    if (yesno) {
                        if (yesno == 1) {
                            if (loud == 1) {
                                printf ("%c", t);
                            } else if (loud > 1) {
                                if (gotype == WORK_PAIR) {
                                    printf ("Z %d %d %d\n", a, b, c);
                                } else {
                                    printf ("Z %d %d %d\n", a, b,
                                              CCelim_dist (a, b, &D));
                                }
                            }
                            fflush (stdout);
                        }
                        if (gotype == WORK_PAIR) {
                            if (processlist[3*id]    != a ||
                                processlist[3*id+1]  != b || 
                                processlist[3*id+2]  != c) {
                                fprintf (stderr, "ERROR: work mismatch\n");
                                rval = 1; goto CLEANUP;
                            }
                            donepairlist[3*elimcount] = a;
                            donepairlist[3*elimcount+1] = b;
                            donepairlist[3*elimcount+2] = c;
                        } else {
                            if (processlist[2*id]    != a ||
                                processlist[2*id+1]  != b) {
                                fprintf (stderr, "ERROR: work mismatch\n");
                                rval = 1; goto CLEANUP;
                            }
                            elimlist[2*elimcount] = a;
                            elimlist[2*elimcount+1] = b;
                            if (gotype == WORK_ELIM) {
                                CCelim_delete_edge (a, b, &G);
                            } else {                   
                                CCelim_add_edge (a, b,
                                    CCelim_dist (a, b, &D), &F);
                            }
                        }
                        round_elimcount++;
                        elimcount++;
                        hits[id] = 1;
                    }
                    fini[id] = 1;
                    done++;
                    if (done % 10000 == 9999) { 
                        if (loud) printf ("\n");
                        printf ("%d %d/%d %.2f %d total, %.0f seconds\n",
                            round_elimcount, done, round_total,
                            (double) round_elimcount / (double) done,
                            elimcount, CCutil_real_zeit() - wallzeit);
                        fflush (stdout);
                    }
                }
            }
            break;

            default:
                fprintf (stderr, "Unknown elimtask type %c\n", request);
                rval = 1; goto CLEANUP;
            } /* end switch */

            rval = CCutil_sclose (s);
            CCcheck_rval (rval, "CCutil_sclose failed");
            s = (CC_SFILE *) NULL;
        }
    
        if (loud) printf ("\n");
        if (use_loop == 0) break;

        printf ("Completed Round %d, Wall Time %.0f seconds (%.0f total)\n",
                   round, CCutil_real_zeit () - round_wallzeit,
                   CCutil_real_zeit () - wallzeit);
        printf ("Round %d Elimation: %d (%d total)\n", round, round_elimcount,
                   elimcount);
        fflush (stdout);

        round++;
        (*round_id)++;
        if ((double) round_elimcount  / (double) round_total < target ||
                                                   gotype == WORK_FIX) {
            setlevel++;
        }

        if (CCutil_real_zeit () - wallzeit > 86400.0) {
            printf ("Total wall time greater than 24 hours.  Breaking loop\n");
            break;
        }
    }

    if (use_loop != 0) printf ("\n");

    if (gotype == WORK_PAIR) {
        printf ("Eliminated Pairs: %d\n", elimcount); fflush (stdout);

        rval = CCelim_pairlist_to_nonpairs (&G, elimcount, donepairlist,
             nonpcounts, nonpairs, outpcounts, outpairs);
        CCcheck_rval (rval, "pairlist_to_nonpairs failed");
        if (pelimcount) *pelimcount = elimcount;
    } else if (gotype == WORK_FIX) {
        printf ("\nFixed: %d\n", elimcount); fflush (stdout);
        rval = CCelim_getedges_graph (&F, remaincount, remainlist);
        CCcheck_rval (rval, "CCelim_getedges_graph failed");
        printf ("Total: %d\n", *remaincount);
    } else {
        printf ("Eliminated: %d\n", elimcount); fflush (stdout);
        rval = CCelim_getedges_graph (&G, remaincount, remainlist);
        CCcheck_rval (rval, "CCelim_getedges_graph failed");
        printf ("Remaining: %d\n", *remaincount);
        if (pelimcount) *pelimcount = elimcount;
    }

    if (boss_wall_time) *boss_wall_time = CCutil_real_zeit () - wallzeit;
    printf ("Wall clock: %.0f seconds\n", CCutil_real_zeit () - wallzeit);
    fflush (stdout);

CLEANUP:
    CC_IFFREE (perm, int);
    CC_IFFREE (round_elist, int);
    CC_IFFREE (round_elen, int);
    CC_IFFREE (hits, int);
    CC_IFFREE (fini, int);
    CC_IFFREE (newplist, int);
    if (s != (CC_SFILE *) NULL) CCutil_sclose (s);
    CCelim_free_graph (&G);
    CCelim_free_graph (&F);
    CCelim_free_distobj (&D);
    CC_IFFREE (donepairlist, int);
    free_pair_info (ncount, &round_nonpcounts, &round_nonpairs);
    return rval;
}

/* elim_grunt is the worker code for a parallel run */

static int elim_grunt (char *theboss, CCrandstate *rstate, int loud,
        double local_timelimit)
{
    int rval = 0, i, j, a, b, c, id, ncount = 0, ecount = 0, level_count = -1;
    int *elist = (int *) NULL, *elen = (int *) NULL, *perm = (int *) NULL;
    int workcount = 0, worklist[MAX_GRUNT_GROUP*4], yesno_list[MAX_GRUNT_GROUP];
    int yesno, max_neighborhood = 5, newgraph, gotype;
    int havepairs = 0, *nonpcounts = (int *) NULL, **nonpairs = (int **) NULL;
    int *pcounts = (int *) NULL, **pairs = (int **) NULL, ht_list_count = 0;
    int fixcount = 0, *fixlist = (int *) NULL, save_tree = 0;
    int pcnt, use_tsp = 1, round_id = -1, wtype, longpath, call_fast;
    int *incyc = (int *) NULL, havecyc = 0;
    CC_SFILE *s = (CC_SFILE *) NULL;
    CCdatagroup dat, *euclid_dat = (CCdatagroup *) NULL;
    CCkdtree *kt = (CCkdtree *) NULL, *euclid_kt = (CCkdtree *) NULL;
    CCelim_graph G, fixG, *F = (CCelim_graph *) NULL;
    CCelim_distobj D;
    double ezeit = 0.0, timelimit = 100000.0;
    CCelim_httree *ht_list = (CCelim_httree *) NULL;
    CCelim_httree *ht = (CCelim_httree *) NULL;

    CCelim_init_graph (&G);
    CCelim_init_graph (&fixG);
    CCutil_init_datagroup (&dat);
    CCelim_init_distobj (&D);

    if (!theboss) {
        fprintf (stderr, "grunt does not know the boss\n");
        rval = 1; goto CLEANUP;
    }

    printf ("Grunt for boss %s\n", theboss); fflush (stdout);

    /* open a connection to boss, say hello, and receive the datagroup */

    rval = open_connection (theboss, &s);
    CCcheck_rval (rval, "open_connection failed");
    rval = CCutil_swrite_char (s, ELIM_HELLO);
    CCcheck_rval (rval, "CCutil_swrite_char failed (hello)");

    rval = CCutil_readmaster (s, &ncount, &dat, &perm);
    CCcheck_rval (rval, "CCutil_readmaster failed");
    printf ("have the datagroup (%d)\n", ncount); fflush (stdout);
    rval = CCutil_sread_int (s, &havecyc);
    CCcheck_rval (rval, "CCutil_sread_int failed (havecyc)");
    if (havecyc) {
        CC_MALLOC (incyc, ncount, int);
        for (i = 0; i < ncount; i++) {
            rval = CCutil_sread_int (s, &(incyc[i]));
            CCcheck_rval (rval, "incyc[] read failed");
        }
        printf ("have the input tour\n"); fflush (stdout);
    }

    CCutil_sclose (s); s = (CC_SFILE *) NULL;

    /* standard alloc, but not graphs G and F (they will be sent later) */

    rval = CCelim_alloc_elim (ncount, ecount, elist, 0, (int *) NULL, &dat,
                  (CCelim_graph *) NULL, (CCelim_graph *) NULL, &D,
                  &kt, &euclid_dat, &euclid_kt, rstate);
    CCcheck_rval (rval, "CCelim_alloc_elim failed");

    while (1) {
        /* ask the boss for a list of edges or 3-paths to eliminate */
        rval = open_connection (theboss, &s);
        CCcheck_rval (rval, "open_connection failed");
        rval = CCutil_swrite_char (s, ELIM_WORK);
        CCcheck_rval (rval, "CCutil_swrite_char failed (work)");
        rval = CCutil_swrite_int (s, round_id);
        CCcheck_rval (rval, "CCutil_swrite_int failed (round_id)");

        rval = CCutil_sread_int (s, &newgraph);
        CCcheck_rval (rval, "CCutil_sread_int failed (newgraph)");
        if (newgraph) {
            /* boss is on a different round, so get the current graph */
            CCelim_free_graph (&G);
            CCelim_free_graph (&fixG);
            CC_IFFREE (elist, int);
            CC_IFFREE (fixlist, int);
            free_pair_info (ncount, &nonpcounts, &nonpairs);

            rval = CCutil_sread_int (s, &round_id);
            CCcheck_rval (rval, "CCutil_sread_int failed (round_id)");
            printf ("\nRound ID: %d\n", round_id); fflush (stdout);

            rval = CCutil_sread_int (s, &gotype);
            CCcheck_rval (rval, "CCutil_swrite_read failed (gotype)");
            switch (gotype) {
                case WORK_PAIR: printf ("Processing pair elimination\n"); break;
                case WORK_FIX:  printf ("Processing edge fixing\n"); break;
                case WORK_ELIM: printf ("Processing edge elimination\n");
            }

            rval = CCutil_sread_int (s, &i);
            CCcheck_rval (rval, "CCutil_sread_int failed (i)");
            if (i != ncount) {
                fprintf (stderr, "new graph does not match datagroup\n");
                rval = 1; goto CLEANUP;
            }
            rval = CCutil_sread_int (s, &ecount);
            CCcheck_rval (rval, "CCutil_sread_int failed (ecount)");
            printf ("ncount = %d, ecount = %d\n", ncount, ecount);
            fflush (stdout);

            CC_MALLOC (elist, 2*ecount, int);
            for (i = 0; i < ecount; i++) {
                rval = CCutil_sread_int (s, &(elist[2*i]));
                CCcheck_rval (rval, "elist[0] read failed");
                rval = CCutil_sread_int (s, &(elist[2*i+1]));
                CCcheck_rval (rval, "elist[1] read failed");
            }

            CC_MALLOC (elen, ecount, int);
            for (i = 0; i < ecount; i++) {
                rval = CCutil_sread_int (s, &(elen[i]));
                CCcheck_rval (rval, "elen read failed");
            }
            for (i = 0; i < ecount; i++) {
                if (elen[i] != CCutil_dat_edgelen (elist[2*i], elist[2*i+1],
                                                   &dat)) {
                    fprintf (stderr, "edge lengths do not match datagroup\n");
                    rval = 1; goto CLEANUP;
                }
            }
            rval = CCelim_build_graph (&G, ncount, ecount, elist, elen,
                                       (CCelim_distobj *) NULL);
            CCcheck_rval (rval, "CCelim_build_graph failed");
            CC_IFFREE (elen, int);
            printf ("have the edge set (%d), elen checks out\n", ecount);
            fflush (stdout);

            rval = CCutil_sread_int (s, &fixcount);
            CCcheck_rval (rval, "CCutil_sread_int failed (fixcount)");
            if (fixcount > 0) {
                CC_MALLOC (fixlist, 2*fixcount, int);
                for (i = 0; i < fixcount; i++) {
                    rval = CCutil_sread_int (s, &(fixlist[2*i]));
                    CCcheck_rval (rval, "CCutil_sread_int failed (fixlist[0])");
                    rval = CCutil_sread_int (s, &(fixlist[2*i+1]));
                    CCcheck_rval (rval, "CCutil_sread_int failed (fixlist[1])");
                }
                F = &fixG;
                rval = CCelim_build_graph (F, ncount, fixcount, fixlist,
                          (int *) NULL, &D);
                CCcheck_rval (rval, "CCelim_build_graph failed");
                printf ("have the fixed edges (%d)\n", fixcount);
                fflush (stdout);
            }

            rval = CCutil_sread_double (s, &timelimit);
            CCcheck_rval (rval, "CCutil_sread_double failed (timelimit)");
            if (local_timelimit != DEFAULT_TIME_LIMIT) {
                timelimit = local_timelimit;
            }
            rval = CCutil_sread_int (s, &level_count);
            CCcheck_rval (rval, "CCutil_sread_int failed (level_count)");
            rval = CCutil_sread_int (s, &wtype);
            CCcheck_rval (rval, "CCutil_sread_int failed (wtype)");
            rval = CCutil_sread_int (s, &longpath);
            CCcheck_rval (rval, "CCutil_sread_int failed (longpath)");
            rval = CCutil_sread_int (s, &use_tsp);
            CCcheck_rval (rval, "CCutil_sread_int failed (use_tsp)");
            rval = CCutil_sread_int (s, &max_neighborhood);
            CCcheck_rval (rval, "CCutil_sread_int failed (max_neighborhood)");
            rval = CCutil_sread_int (s, &call_fast);
            CCcheck_rval (rval, "CCutil_sread_int failed (call_fast)");
            rval = CCutil_sread_int (s, &save_tree);
            CCcheck_rval (rval, "CCutil_sread_int failed (save_tree)");

            printf ("timelimit = %f\n", timelimit);
            printf ("level_count = %d\n", level_count);
            printf ("wtype = %d\n", wtype);
            printf ("longpath = %d\n", longpath);
            printf ("use_tsp = %d\n", use_tsp);
            printf ("max_neighborhood = %d\n", max_neighborhood);
            printf ("call_fast = %d\n", call_fast);
            printf ("save_tree = %d\n", save_tree);
            fflush (stdout);

            rval = CCutil_sread_int (s, &havepairs);
            CCcheck_rval (rval, "CCutil_sread_int failed (havepairs)");
            printf ("havepairs = %d\n", havepairs);
            if (havepairs) {
                CC_MALLOC (nonpcounts, ncount, int);
                CC_MALLOC (nonpairs, ncount, int *);
                for (i = 0; i < ncount; i++) nonpairs[i] = (int *) NULL;
                pcnt = 0;
                for (i = 0; i < ncount; i++) {
                    rval = CCutil_sread_int (s, &(nonpcounts[i]));
                    CCcheck_rval (rval, "CCutil_sread_int failed (nonpcounts)");
                    if (nonpcounts[i]) {
                        pcnt += nonpcounts[i];
                        CC_MALLOC (nonpairs[i], 2*nonpcounts[i], int);
                        for (j = 0; j < nonpcounts[i]; j++) {
                            rval = CCutil_sread_int (s, &(nonpairs[i][2*j]));
                            CCcheck_rval (rval, "sread_int failed (nonpairs)");
                            rval = CCutil_sread_int (s, &(nonpairs[i][2*j+1]));
                            CCcheck_rval (rval, "sread_int failed (nonpairs)");
                        }
                    }
                }
                printf ("have the pairs (%d)\n", pcnt); fflush (stdout);
            }
        }

        rval = CCutil_sread_int (s, &workcount);
        CCcheck_rval (rval, "CCutil_sread_int failed (workcount)");
        printf ("\nProcess %d %s\n", workcount,
                         (gotype == WORK_PAIR ? "pairs" : "edges"));
        fflush (stdout);
        for (i = 0; i < workcount; i++) {
            rval = CCutil_sread_int (s, &id);
            CCcheck_rval (rval, "CCutil_sread_int failed (id)");
            rval = CCutil_sread_int (s, &a);
            CCcheck_rval (rval, "CCutil_sread_int failed (a)");
            rval = CCutil_sread_int (s, &b);
            CCcheck_rval (rval, "CCutil_sread_int failed (b)");
            if (gotype == WORK_PAIR) {
                rval = CCutil_sread_int (s, &c);
                CCcheck_rval (rval, "CCutil_sread_int failed (c)");
                worklist[4*i] = id;
                worklist[4*i+1] = a;
                worklist[4*i+2] = b;
                worklist[4*i+3] = c;
            } else {
                worklist[3*i] = id;
                worklist[3*i+1] = a;
                worklist[3*i+2] = b;
            }
        }
        CCutil_sclose (s); s = (CC_SFILE *) NULL;

        if (workcount == 0) break;

        if (nonpairs) {
            free_pair_info (ncount, &pcounts, &pairs);
            rval = CCelim_nonpairs_to_pairs (&G, nonpcounts, nonpairs, &pcounts,
                                             &pairs);
            CCcheck_rval (rval, "CCelim_nonpairs_to_pairs failed");
        }

        if (save_tree) {
            CC_MALLOC (ht_list, workcount, CCelim_httree);
            for (i = 0; i < workcount; i++) {
                CCelim_init_httree (&(ht_list[i]));
                ht_list[i].max_count = HT_MAX_COUNT;
                ht_list[i].max_depth = HT_MAX_DEPTH;
            }
            ht_list_count = workcount;
        }

        if (gotype == WORK_ELIM || gotype == WORK_FIX) {
            for (i = 0; i < workcount; i++) {
                a = worklist[3*i+1];
                b = worklist[3*i+2];
                if (!CCelim_edge_in_graph (a, b, &G)) {
                    printf ("Working edge (%d, %d) not in graph!\n", a, b);
                    exit (1);
                }

                if (save_tree) ht = &(ht_list[i]);
                else           ht = (CCelim_httree *) NULL;

                if (gotype == WORK_FIX) {
                    rval = CCelim_run_fix_edge (a, b, &G, F, &D, rstate, kt,
                      euclid_dat, euclid_kt, &yesno, &ezeit, level_count,
                      pcounts, pairs, timelimit, use_tsp, max_neighborhood,
                      ht);
                    CCcheck_rval (rval, "CCelim_run_fix_edge failed");
                } else if (call_fast == 1) {
                    rval = CCelim_run_fast_elim_edge (a, b, &G, &D, rstate, kt,
                       euclid_dat, euclid_kt, level_count, max_neighborhood,
                       &yesno, ht);
                    CCcheck_rval (rval, "CCelim_run_fast_elim_edge failed");
                } else {
                    rval = CCelim_run_elim_edge (a, b, &G, F, &D, wtype,
                       rstate, &yesno, &ezeit, level_count, pcounts, pairs,
                       kt, euclid_dat, euclid_kt, timelimit, longpath, use_tsp,
                       max_neighborhood, ht);
                    CCcheck_rval (rval, "CCelim_run_elim_edge failed");
                }
                printf ("%d", yesno); fflush (stdout);
                if (loud) {
                    printf (" (%d, %d) %.2lf\n", a, b, ezeit); fflush (stdout);
                }
                if (yesno && ht_list) {
/* BICO 220927 
                    if (!ht_list[i].root || (ht_list[i].root &&
                         ht_list[i].root->tutte_type == CCelim_TUTTE_NONE)) {
*/
                    if (!ht_list[i].root) {
                        fprintf (stderr, "missing HT tree info\n");
                        rval = 1; goto CLEANUP;
                    }
                }
                yesno_list[i] = yesno;
                if (yesno == 1 && gotype == WORK_ELIM) {
                    CCelim_delete_edge (a, b, &G);
                }
            }
        } else {
            for (i = 0; i < workcount; i++) {
                a = worklist[4*i+1];
                b = worklist[4*i+2];
                c = worklist[4*i+3];

                if (save_tree) ht = &(ht_list[i]);
                else           ht = (CCelim_httree *) NULL;

                /* BICO Why is longpath not being used? */
                rval = CCelim_run_elim_pair (b, a, c, &G, F, &D,
                   rstate, kt, euclid_dat, euclid_kt, &yesno, &ezeit,
                   level_count, pcounts, pairs, timelimit, longpath, use_tsp,
                   max_neighborhood, ht);
                CCcheck_rval (rval, "CCelim_elim_pair failed");
                printf ("%d", yesno); fflush (stdout);
                if (loud) {
                    printf (" (%d, %d, %d) %.2lf\n", a, b, c, ezeit);
                    fflush (stdout);
                }
                yesno_list[i] = yesno;
            }
        }

        /* send the results to the boss, together with our ID */

        rval = open_connection (theboss, &s);
        CCcheck_rval (rval, "open_connection failed");

        rval = CCutil_swrite_char (s, ELIM_DONE);
        CCcheck_rval (rval, "CCutil_swrite_char failed (work)");
        rval = CCutil_swrite_int (s, round_id);
        CCcheck_rval (rval, "CCutil_swrite_int failed (round_id)");
        rval = CCutil_swrite_int (s, gotype);
        CCcheck_rval (rval, "CCutil_swrite_int failed (gotype)");

        rval = CCutil_swrite_int (s, workcount);
        CCcheck_rval (rval, "CCutil_swrite_int failed (workcount)");
        for (i = 0; i < workcount; i++) {
            if (gotype == WORK_ELIM || gotype == WORK_FIX) {
                id = worklist[3*i];
                a = worklist[3*i+1];
                b = worklist[3*i+2];
                yesno = yesno_list[i];
            } else {
                id = worklist[4*i];
                a = worklist[4*i+1];
                b = worklist[4*i+2];
                c = worklist[4*i+3];
                yesno = yesno_list[i];
            }
            rval = CCutil_swrite_int (s, id);
            CCcheck_rval (rval, "CCutil_swrite_int failed (id)");
            rval = CCutil_swrite_int (s, a);
            CCcheck_rval (rval, "CCutil_swrite_int failed (a)");
            rval = CCutil_swrite_int (s, b);
            CCcheck_rval (rval, "CCutil_swrite_int failed (b)");
            if (gotype == WORK_PAIR) {
                rval = CCutil_swrite_int (s, c);
                CCcheck_rval (rval, "CCutil_swrite_int failed (c)");
            }
            rval = CCutil_swrite_int (s, yesno);
            CCcheck_rval (rval, "CCutil_swrite_int failed (yesno)");
            if (yesno && save_tree) {
/* BICO 220927
                if (ht_list && ht_list[i].root &&
                    ht_list[i].root->tutte_type != CCelim_TUTTE_NONE) {
*/
                if (ht_list && ht_list[i].root) {
                    rval = CCutil_swrite_int (s, 1);
                    CCcheck_rval (rval, "CCutil_swrite_int failed (1)");
                    rval = CCelim_swrite_httree (&(ht_list[i]), s);
                    CCcheck_rval (rval, "CCelim_swrite_httree failed");
                } else {
                    rval = CCutil_swrite_int (s, 0);
                    CCcheck_rval (rval, "CCutil_swrite_int failed (0)");
                }
            }
        }

        CCutil_sclose (s);
        s = (CC_SFILE *) NULL;

        if (ht_list) {
            for (i = 0; i < ht_list_count; i++) {
                CCelim_free_htnode_children (ht_list[i].root);
                CC_IFFREE (ht_list[i].root, CCelim_htnode);
            }
            CC_IFFREE (ht_list, CCelim_httree);
        }
    }

    printf ("No work available. Shutting down.\n"); fflush (stdout);

CLEANUP:
    CC_IFFREE (elist, int);
    CC_IFFREE (elen, int);
    CC_IFFREE (perm, int);
    CC_IFFREE (incyc, int);
    if (s != (CC_SFILE *) NULL) CCutil_sclose (s);
    CCutil_freedatagroup (&dat);
    CCelim_free_elim (&G, (CCelim_graph *) NULL, &D, &kt, &euclid_dat,
                      &euclid_kt);
    free_pair_info (ncount, &nonpcounts, &nonpairs);
    free_pair_info (ncount, &pcounts, &pairs);
    if (ht_list) {
        for (i = 0; i < ht_list_count; i++) {
            CCelim_free_htnode_children (ht_list[i].root);
            CC_IFFREE (ht_list[i].root, CCelim_htnode);
        }
        CC_IFFREE (ht_list, CCelim_httree);
    }

    return rval;
}

static int open_connection (char *theboss, CC_SFILE **s)
{
    *s = CCutil_snet_open (theboss, ELIM_PORT);
    if (!(*s)) {
        fprintf (stderr, "CCutil_snet_open failed\n");
        return 1;
    } else {
        return 0;
    }
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

/* processlist_cycle updates the list of edges to process: when eliminating */
/* edges, only process edges not in the tour; when fixing edges, only       */
/* process edges in the tour                                                */

static int processlist_cycle (int ncount, int *incyc, int worktype,
        int processcount, int *processlist, int *new_processcount,
        int **new_processlist)
{
    CCutil_edgehash h;
    int rval = 0, i, a, b, c, tmp, test;
    int tprocesscount = 0, *tprocesslist = (int *) NULL;

    /* load the tour edges into a hash table */

    rval = CCutil_edgehash_init (&h, 10*ncount);
    CCcheck_rval (rval, "CCutil_edgehash_init failed");
    for (i = 0; i < ncount; i++) {
        rval = CCutil_edgehash_set (&h, incyc[i], incyc[(i+1)%ncount], 1);
        CCcheck_rval (rval, "CCutil_edgehash_set failed");
    }

    if (worktype == WORK_FIX) {
        CC_MALLOC (tprocesslist, 2*ncount, int);
    } else if (worktype == WORK_ELIM) {
        CC_MALLOC (tprocesslist, 2*processcount, int);
    } else {
        CC_MALLOC (tprocesslist, 3*processcount, int);
    }

    /* edgehash_find returns -1 if edge is not in the hash table */

    if (worktype != WORK_PAIR) {
        for (i = 0; i < processcount; i++) {
            a = processlist[2*i];  b = processlist[2*i+1];
            test = CCutil_edgehash_find (&h, a, b, &tmp);
            if ((worktype == WORK_ELIM && test == -1) ||
                (worktype == WORK_FIX  && test != -1)) {
                tprocesslist[2*tprocesscount] = a;
                tprocesslist[2*tprocesscount+1] = b;
                tprocesscount++;
            }
        }
    } else {
        for (i = 0; i < processcount; i++) {
            a = processlist[3*i];  b = processlist[3*i+1];
            c = processlist[3*i+2];
            if (CCutil_edgehash_find (&h, a, b, &tmp) == -1 ||
                CCutil_edgehash_find (&h, a, c, &tmp) == -1) {
                tprocesslist[3*tprocesscount] = a;
                tprocesslist[3*tprocesscount+1] = b;
                tprocesslist[3*tprocesscount+2] = c;
                tprocesscount++;
            }
        }
    }

    CCutil_edgehash_free (&h);
    *new_processcount = tprocesscount;
    *new_processlist = tprocesslist;

CLEANUP:
    if (rval) { CC_IFFREE (tprocesslist, int); }
    return rval;
}

static int read_pair_list (char *fname, int ncount, int ecount, int *elist,
        int *pcount, int **plist, int *nonpcounts, int **nonpairs)
{
    int rval = 0, count, i, j, icount = 0, a, b, c, tmp, *list = (int *) NULL;
    FILE *in = (FILE *) NULL;
    CCelim_graph G;

    rval = CCelim_build_graph (&G, ncount, ecount, elist, (int *) NULL,
                 (CCelim_distobj *) NULL);
    CCcheck_rval (rval, "CCelim_build_graph failed");

    in = fopen (fname, "r");
    if (!in) {
        fprintf (stderr, "could not open %s for reading\n", fname);
        rval = 1; goto CLEANUP;
    }

    fscanf (in, "%d %d", &i, &count);
    if (i != ncount) {
        fprintf (stderr, "nonpair list file does not match problem\n");
        rval = 1; goto CLEANUP;
    }

    CC_MALLOC (list, 3*count, int)

    for (i = 0; i < count; i++) {
        fscanf (in, "%d %d %d", &a, &b, &c);
        if (b > c) { CC_SWAP (b, c, tmp); }
        if (CCelim_edge_in_graph (a,b,&G) && CCelim_edge_in_graph (a,c,&G)) {
            if (nonpairs) {
                for (j = 0; j < nonpcounts[a]; j++) {
                    if (b == nonpairs[a][2*j] && c == nonpairs[a][2*j+1]) break;
                }
                if (j == nonpcounts[a]) {
                   list[3*icount] = a;
                   list[3*icount+1] = b;
                   list[3*icount+2] = c;
                   icount++;
                } else {
                   printf ("Note: pair %d %d %d in nonpair list\n", a, b, c);
                   fflush (stdout);
                }
            } else {
                list[3*icount] = a;
                list[3*icount+1] = b;
                list[3*icount+2] = c;
                icount++;
            }
        } else {
           printf ("Note: pair %d %d %d not in graph\n", a, b, c);
           fflush (stdout);
        }
    }

    *pcount = icount;
    *plist = list;

    printf ("Kept %d/%d pairs\n", icount, count); fflush (stdout);

CLEANUP:
    CCelim_free_graph (&G);
    if (in) fclose (in);
    return rval;
}

#define KNEIGH 50

/* count_crossings is used when -x is specifed to eliminate crossing edges */

static int count_crossings (int ncount, CCdatagroup *dat, CCrandstate *rstate,
        int ecount, int *elist, int *ccount, int **clist)
{
    int rval = 0, neigh[KNEIGH], a, b, c, d, i, j, k, tmp, hit, havetree = 0;
    int norm;
    CCelim_graph G;
    CCkdtree kt;

    *ccount = 0;
    CCelim_init_graph (&G);

    CCutil_dat_getnorm (dat, &norm);
    if ((norm & CC_NORM_BITS) != CC_KD_NORM_TYPE) {
        printf ("To count crossings the norm must be of KD-tree type\n");
        rval = 1; goto CLEANUP;
    }

    rval = CCelim_build_graph (&G, ncount, ecount, elist, (int *) NULL,
                               (CCelim_distobj *) NULL);
    CCcheck_rval (rval, "CCelim_build_graph failed");

    rval = CCkdtree_build (&kt, ncount, dat, (double *) NULL, rstate);
    CCcheck_rval (rval, "CCkdtree_build failed");
    havetree = 1;
    CC_MALLOC (*clist, 2*ecount, int);

    for (i = 0; i < ecount; i++) {
        a = elist[2*i];  b = elist[2*i+1];
        if (a > b) CC_SWAP (a, b, tmp);
        rval = CCkdtree_node_k_nearest (&kt, ncount, a, KNEIGH, dat,
                          (double *) NULL, neigh, rstate);
        CCcheck_rval (rval, "CCkdtree_node_k_nearest failed");
        hit = 0;
        for (j = 0; j < KNEIGH && !hit; j++) {
            c = neigh[j];
            if (c == a || c == b) continue;
            for (k = j+1; k < KNEIGH && !hit; k++) {
                d = neigh[k];
                if (d == a || d == b) continue;
                if (CCelim_edge_in_graph (c, d, &G) &&
                    segments_intersect (dat->x[a], dat->y[a], dat->x[b],
                                        dat->y[b], dat->x[c], dat->y[c],
                                        dat->x[d], dat->y[d])) {
                    (*clist)[2*(*ccount)]   = a;
                    (*clist)[2*(*ccount)+1] = b;
                    (*ccount)++;
                    hit = 1;
                }
            }
        }
    }

CLEANUP:
    if (havetree) CCkdtree_free (&kt);
    CCelim_free_graph (&G);
    return rval;
}

static int segments_intersect (double p0_x, double p0_y, double p1_x,
        double p1_y, double p2_x, double p2_y, double p3_x, double p3_y)
{
    double s1_x, s1_y, s2_x, s2_y, s, t;

    s1_x = p1_x - p0_x;     s1_y = p1_y - p0_y;
    s2_x = p3_x - p2_x;     s2_y = p3_y - p2_y;

    s = (-s1_y * (p0_x - p2_x) + s1_x * (p0_y - p2_y)) /
        (-s2_x * s1_y + s1_x * s2_y);

    t = ( s2_x * (p0_y - p2_y) - s2_y * (p0_x - p2_x)) /
        (-s2_x * s1_y + s1_x * s2_y);

    if (s >= 0.0 && s <= 1.0 && t >= 0.0 && t <= 1.0)
        return 1;
    else
        return 0;
}

static int parseargs (int ac, char **av)
{
    int c, inorm, boptind = 1;
    char *boptarg = (char *) NULL;

    while ((c = CCutil_bix_getopt (ac, av, "Aabcde:EfF:g:hH:I:j:k:l:M:N:n:o:pP:s:t:T:vw:xz:", &boptind, &boptarg)) != EOF) {
        switch (c) {
        case 'A': use_full_loop = 1; break;
        case 'a': use_single_fast=1; use_fast_elim=1; use_level_loop=1; break;
        case 'b': use_single_fast=3; use_fast_elim=1; use_level_loop=1; break;
        case 'c': use_fast_elim = 1; break;
        case 'd': use_level_loop = 1; break;
        case 'e': elimfname = boptarg; break;
        case 'E': beverbose++; break;
        case 'f': worktype = WORK_FIX; break;
        case 'F': fixedfname = boptarg; break;
        case 'g': grunthostname = boptarg; break;
        case 'h': be_nethost = 1; break;
        case 'H': save_hamilton_tutte_tree = 1; htfname = boptarg;  break;
        case 'I': infname = boptarg; break;
        case 'j': use_longpath = atoi (boptarg);
                  if (use_longpath < 4) use_longpath = 4;
                  break;
        case 'k': use_tsp_swapping = atoi (boptarg); break;
        case 'l': elim_time_limit = atof (boptarg); break;
        case 'n': use_max_neighborhood = atoi (boptarg);
                  if (use_max_neighborhood < 5)  use_max_neighborhood = 5;
                  if (use_max_neighborhood > CCelim_MAX_AB) 
                      use_max_neighborhood = CCelim_MAX_AB;
                  break;
        case 'N':
            inorm = atoi (boptarg);
            switch (inorm) {
            case 0: use_norm = CC_MAXNORM;  break;
            case 1: use_norm = CC_MANNORM;  break;
            case 2: use_norm = CC_EUCLIDEAN;  break;
            case 3: use_norm = CC_EUCLIDEAN_3D;  break;
            case 4: use_norm = CC_USER;  break;
            case 5: use_norm = CC_ATT;  break;
            case 6: use_norm = CC_GEOGRAPHIC;  break;
            case 7: use_norm = CC_MATRIXNORM;  break;
            case 8: use_norm = CC_DSJRANDNORM;  break;
            case 9: use_norm = CC_CRYSTAL;  break;
            case 10: use_norm = CC_SPARSE;  break;
            case 11: use_norm = CC_RHMAP1;  break;
            case 12: use_norm = CC_RHMAP2;  break;
            case 13: use_norm = CC_RHMAP3;  break;
            case 14: use_norm = CC_RHMAP4;  break;
            case 15: use_norm = CC_RHMAP5;  break;
            case 16: use_norm = CC_EUCTOROIDAL;  break;
            case 17: use_norm = CC_GEOM;  break;
            case 18: use_norm = CC_EUCLIDEAN_CEIL;  break;
            case 20: use_norm = CC_ROAD;  break;
            default:
                usage (av[0]);
                return 1;
            }
            tsplib_in = 0;
            break;
        case 'o': outfname = boptarg; break;
        case 'p': worktype = WORK_PAIR; break;
        case 'P': nonpairsfname = boptarg; break;
        case 's': seed = atoi (boptarg); break;
        case 't': tourfname = boptarg; break;
        case 'T': tspfname = boptarg; break;
        case 'v': beverbose = 2; break;
        case 'w': witness_type = atoi (boptarg); break;
        case 'x': only_crossing_edges = 1; break;
        case 'z': point_levels = atoi (boptarg);
                  if      (point_levels < 2) point_levels = 2;
                  /* else if (point_levels > 9) point_levels = 9; */
                  break;
        case CC_BIX_GETOPT_UNKNOWN:
        case '?':
        default:
            usage (av[0]);
            return 1;
        }
    }

    if (boptind < ac) { edgefname = av[boptind++]; }

    if (grunthostname == (char *) NULL &&
       (tspfname == (char *) NULL || edgefname == (char *) NULL)) {
        usage (av[0]);
        return 1;
    }

    return 0;
}

static void usage (char *f)
{
    fprintf (stderr, "Usage: %s [-see below-] -T tsp_file edge_file\n", f);
    fprintf (stderr, "   -A    full loop (fast+pairs+elim+fix: can be slow)\n");
    fprintf (stderr, "   -a    first level of -cd elim loop (usually fast)\n");
    fprintf (stderr, "   -b    first three levels of -cd elim loop\n");
    fprintf (stderr, "   -c    call fast elim code\n");
    fprintf (stderr, "   -d    loop through pre-set levels (can be slow if not called with -b)\n");
    fprintf (stderr, "   -e f  dump eliminated (fixed) edges to file f\n");
    fprintf (stderr, "   -E    print an E for each eliminated edge\n");
    fprintf (stderr, "   -f    fix edges to 1 (not eliminate)\n");
    fprintf (stderr, "   -F f  input fixed edges in edge-file f\n");
    fprintf (stderr, "   -g h  be a grunt for host h\n");
    fprintf (stderr, "   -h    be a boss for a parallel elimination run\n");
    fprintf (stderr, "   -H f  save Hamilton-Tutte trees for eliminated edges infile f\n");
    fprintf (stderr, "   -I f  only process edges in edge-file f\n");
    fprintf (stderr, "   -j #  use long path for elimination, set # >= 4\n");
    fprintf (stderr, "   -k #  use up to #-swaps (default is TSP swapping)\n");
    fprintf (stderr, "   -l d  double, time limit in seconds for each edge\n");
    fprintf (stderr, "   -n #  size of search area (default 10)\n");
    fprintf (stderr, "   -o f  dump the remaining edges to file f\n");
    fprintf (stderr, "   -P f  input eliminated pairs in file f\n");
    fprintf (stderr, "   -p    try to eliminate adjacent pairs\n");
    fprintf (stderr, "   -s #  nonzero random seed (randomize process order)\n");
    fprintf (stderr, "   -t f  exclude (include for -f) edges in tour-file f\n");
    fprintf (stderr, "   -T f  TSPLIB file (or dat file with -N)\n");
    fprintf (stderr, "   -v    print run times and eliminated edges\n");
    fprintf (stderr, "   -w #  specify witness type (0=nonedge, 1=edge, 2=dist2, 3=quick), default 1\n");
    fprintf (stderr, "   -x    only try crossing edges\n");
    fprintf (stderr, "   -z #  points to analyse, range 2 to 9, default 2\n");
    fprintf (stderr, "   -N #  norm (must specify if -T gives a dat file)\n");
    fprintf (stderr, "         0=MAX, 1=L1, 2=L2, 3=3D, 4=USER, 5=ATT, 6=GEO, 7=MATRIX,\n");
    fprintf (stderr, "         8=DSJRAND, 9=CRYSTAL, 10=SPARSE, 11-15=RH-norm 1-5, 16=TOROIDAL\n");
    fprintf (stderr, "         17=GEOM, 18=JOHNSON, 20=ROAD\n");
}
