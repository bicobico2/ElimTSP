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


#ifndef __ELIM_H
#define __ELIM_H

#define CCelim_PMAX 20    /* 100 */
#define CCelim_NMAX 200   /* 100 */
#define CCelim_MAX_AB 100 /* 100 */
#define CCelim_MAX_DEGREE 5000 /* 50, skip nodes of high degree */

#define CCelim_CD_NONEDGE (0)
#define CCelim_CD_EDGE    (1)
#define CCelim_CD_DIST2   (2)
#define CCelim_CD_TWO     (3)

typedef struct CCelim_node {
    int *neighbors;
    int *len;    /* lengths of edges to neighbors */
    int deg;
    int maxdeg;  /* max space to add edges */
} CCelim_node;

typedef struct CCelim_graph {
    int ncount;
    int ecount;
    CCelim_node *nodelist;
    int *neighborspace;
    int *lenspace;
} CCelim_graph;

typedef struct CCelim_distobj {
    CCdatagroup *dat;
    int       *cacheval;
    int       *cacheind;
} CCelim_distobj;

typedef struct CCelim_array {
    int count;
    int *arr;
} CCelim_array;

#define CC_ELIM_DIST_IN                                               \
static int dist (int i, int j, CCelim_distobj *D);                    \
static int dist (int i, int j, CCelim_distobj *D)                     \
{                                                                     \
    int ind;                                                          \
                                                                      \
    if (i > j) { CC_SWAP (i, j, ind); }                               \
                                                                      \
    ind = i ^ j;                                                      \
                                                                      \
    if (D->cacheind[ind] != i) {                                      \
        D->cacheind[ind] = i;                                         \
        D->cacheval[ind] = CCutil_dat_edgelen (i, j, D->dat);         \
    }                                                                 \
    return D->cacheval[ind];                                          \
}                                                                     \

typedef struct CCelim_path {
    int          cnt;
    int          node[CCelim_NMAX];
    struct CCelim_path *next;
} CCelim_path;

#define CCelim_HTNODE_MAX_PATH_LEN 48  /* 20 */
#define CCelim_HTTREE_EDGE 1
#define CCelim_HTTREE_FIX  2
#define CCelim_HTTREE_PAIR 3

#define CCelim_HAMILTON_NONE 0
#define CCelim_HAMILTON_EDGE 1
#define CCelim_HAMILTON_PATH 2
#define CCelim_HAMILTON_PAIR 3
#define CCelim_HAMILTON_LONG 4

#define CCelim_TUTTE_NONE -1
#define CCelim_TUTTE_PATH  0
#define CCelim_TUTTE_END   1
#define CCelim_TUTTE_POINT 2
#define CCelim_TUTTE_CD    3

typedef struct CCelim_htnode {
    int hamilton_type;           /* 1:edg 2:3path 3:pair-3path 4:long path */
    int hamilton_path_len;       /* length of long path                    */
    int hamilton_nodes[CCelim_HTNODE_MAX_PATH_LEN];   
                                 /* node indices for hamilton move         */
    int tutte_type;              /* -1:none 0:path 1:end 2:point 3:cd-pair */
    int tutte_nodes[2];          /* node indices for tutte move            */
    int id;                      /* used when saving elim tree to a file   */
    struct CCelim_htnode *child; /* a list of the next hamilton moves      */
    struct CCelim_htnode *next;  /* used to link the children              */
} CCelim_htnode;

typedef struct CCelim_httree {
    CCelim_htnode *root;         /* root of the Hamilton-Tutte tree        */
    int elimtype;                /* specify edge, fix, or pair             */
    int count;                   /* number of nodes in the tree            */
    int depth;                   /* depth (max over all nodes) of the tree */
    int max_count;               /* upper bound on # nodes in the tree     */
    int max_depth;               /* upper bound on depth of the tree       */
    int level;                   /* parameter for CCelim_run_elim edge     */
    int longpath;                /* parameter for CCelim_run_elim edge     */
    int use_tsp;                 /* parameter for CCelim_run_elim edge     */
    int max_neighborhood;        /* parameter for CCelim_run_elim edge     */
} CCelim_httree;

int CCelim_alloc_elim (int ncount, int ecount, int *elist, int fixcount,
    int *fixlist, CCdatagroup *dat, CCelim_graph *G, CCelim_graph *F,
    CCelim_distobj *D, CCkdtree **kt, CCdatagroup **euclid_dat,
    CCkdtree **euclid_kt, CCrandstate *rstate);
void CCelim_free_elim (CCelim_graph *G, CCelim_graph *F, CCelim_distobj *D,
    CCkdtree **kt, CCdatagroup **euclid_dat, CCkdtree **euclid_kt);

int CCelim_run_elim_edge (int a, int b, CCelim_graph *G, CCelim_graph *F,
    CCelim_distobj *D, int wtype, CCrandstate *rstate, int *yesno,
    double *ezeit, int level_count, int *pcounts, int **pairs, CCkdtree *kt,
    CCdatagroup *euclid_dat, CCkdtree *euclid_kt, double timelimit,
    int longpath, int use_tsp, int max_neighborhood, CCelim_httree *ht);
int CCelim_run_fix_edge (int a, int b, CCelim_graph *G, CCelim_graph *F,
    CCelim_distobj *D, CCrandstate *rstate, CCkdtree *kt,
    CCdatagroup *euclid_dat, CCkdtree *euclid_kt, int *yesno, double *ezeit,
    int level, int *pcount, int **pairs, double timelimit, int use_tsp,
    int max_neighborhood, CCelim_httree *ht);
int CCelim_run_elim_pair (int a, int b, int c, CCelim_graph *G,
    CCelim_graph *F, CCelim_distobj *D, CCrandstate *rstate, CCkdtree *kt,
    CCdatagroup *euclid_dat, CCkdtree *euclid_kt, int *yesno, double *ezeit,
    int level, int *pcounts, int **pairs, double timelimit, int longpath,
    int use_tsp, int max_neighborhood, CCelim_httree *ht);
int CCelim_run_fast_elim_edge (int a, int b, CCelim_graph *G,
    CCelim_distobj *D, CCrandstate *rstate, CCkdtree *kt,
    CCdatagroup *euclid_dat, CCkdtree *euclid_kt, int level,
    int max_neighborhood, int *yesno, CCelim_httree *ht);
void CCelim_all_neighbor_pairs (int n, CCelim_graph *G, CCelim_graph *F,
    int *pcount, int *plist);
int CCelim_build_ablist (int a, int b, CCelim_graph *G, CCelim_distobj *D,
    CCelim_array *ablist, int max_neighborhood, CCkdtree *kt,
    CCrandstate *rstate, CCdatagroup *euclid_dat, CCkdtree *euclid_kt);
void CCelim_get_partner_pairs (CCelim_graph *G, int a, int n, int *pcounts,
    int **pairs, int *bcount, int *blist);
void CCelim_check_neighbors_three_swap (int a, int b, int n, CCelim_graph *G,
    CCelim_distobj *D, int *tcount, int *tlist, CCelim_graph *F,
    int *pcounts, int **pairs);
void CCelim_pair_neighbors_three_swap (int a, int b, int c, int n,
    CCelim_graph *G, CCelim_distobj *D, int *tcount, int *tlist,
    int *pcounts, int **pairs);

void CCelim_init_graph (CCelim_graph *G);
void CCelim_free_graph (CCelim_graph *G);
void CCelim_delete_edge (int a, int b, CCelim_graph *G);
int CCelim_add_edge (int a, int b, int len, CCelim_graph *G);

int CCelim_build_graph (CCelim_graph *G, int ncount, int ecount, int *elist,
    int *elen, CCelim_distobj *D);
int CCelim_getedges_graph (CCelim_graph *G, int *ocount, int **olist);
int CCelim_edge_in_graph (int a, int b, CCelim_graph *G);
void CCelim_get_partner_pairs (CCelim_graph *G, int a, int n, int *pcounts,
    int **pairs, int *bcount, int *blist);

void CCelim_init_distobj (CCelim_distobj *D);
void CCelim_free_distobj (CCelim_distobj *D);
int CCelim_build_distobj (CCelim_distobj *D, int ncount, CCdatagroup *dat);

int CCelim_split_cd (int a, int b, int c, int d, int ccount, int *clist,
    int dcount, int *dlist, CCelim_distobj *D, int *good);

int CCelim_dist (int i, int j, CCelim_distobj *D);
int CCelim_tour_solver (int n, int **M, int *optlen, int *otour, int target);

void CCelim_compare_three_swap (int a, int b, int c, int d, int e, int f,
    int **M, int *good, int *I);
void CCelim_compare_four_swap (int a, int b, int c, int d, int e, int f,
    int g, int h, int **M, int *good, int *I);
void CCelim_compare_five_swap (int a, int b, int c, int d, int e, int f, int g,
    int h, int i, int j, int **M, int *good, int *I);

void CCelim_tsp_swap (int nodecount, int *nodes, int pathlen, int pathcount,
    int *omatch, int **M, int *yesno, int *imatch, CCelim_graph *G);
void CCelim_lk_swap (CCelim_path *psys, int nodecount, int *invnames,
    int pathlen, int pathcount, int *omatch, int **M, int *yesno,
    int *imatch);
void CCelim_try_swaps (int *ab, CCelim_path *psys, int nct, int *invnames,
    int *omatch, int **M, int *yesno, int *imatch, int use_tsp);

int CCelim_validate_paths (int ncount, CCelim_path *psys, int *yesno,
    CCelim_path *pmerge);
int CCelim_test_paths (int *ab, CCelim_path *psystem, CCelim_distobj *D,
    CCelim_graph *G, int use_tsp, int *bad);
void CCelim_copy_path_system (CCelim_path *porig, CCelim_path *pcopy,
    int reverse);
void CCelim_paths_to_edges (CCelim_path *pathsystem, int *tcount, int *t);
int CCelim_path_system_count (CCelim_path *p);
void CCelim_path_system_nodeset (CCelim_path *psys, int *nodecount,
    int *nodeset);
void CCelim_path_system_mark (CCelim_path *psys, int *marks, int marker);

int CCelim_read_nonpairs (char *fname, int ncount, int **acounts,
    int ***apairs);
int CCelim_write_nonpairs (char *fname, int ncount, int *pcounts, int **pairs);
int CCelim_nonpairs_to_pairs (CCelim_graph *G, int *nonpcounts, int **nonpairs,
    int **pcounts, int ***pairs);
int CCelim_build_pairprocesslist (int ncount, int ecount, int *elist,
    int *nonpcounts, int **nonpairs, int *processcount, int **processlist);
int CCelim_pairlist_to_nonpairs (CCelim_graph *G, int count, int *plist,
    int *nonpcounts, int **nonpairs, int **outpcounts, int ***outpairs);

void CCelim_hamilton_find (CCelim_htnode *ht, int ham_type, int ham_path_len,
    int ham_nodes[6], CCelim_htnode **child);
void CCelim_init_htnode (CCelim_htnode *n);
void CCelim_free_htnode_children (CCelim_htnode *n);
void CCelim_init_httree (CCelim_httree *ht);
void CCelim_free_httree (CCelim_httree *ht);
int CCelim_write_httree (CCelim_httree *ht, char *fname, int append,
    int only_edges);
int CCelim_swrite_httree (CCelim_httree *ht, CC_SFILE *s);
int CCelim_print_httree (CCelim_httree *ht, FILE *out, int only_edges);
int CCelim_read_httree (CCelim_httree **ht, FILE *in, CC_SFILE *s);
int CCelim_copy_httree (CCelim_httree *ht, CCelim_httree **pout);

#endif  /* __ELIM_H */

