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
/*  Written by:  Cook, 2018                                                 */
/*  Updates: February 2018, July 2018, August 2022                          */
/*                                                                          */
/*  Exported Functions                                                      */
/*                                                                          */
/*  void CCelim_hamilton_find (CCelim_htnode *ht, int ham_type,             */
/*      int ham_path_len, int ham_nodes[CCelim_HTNODE_MAX_PATH_LEN],        */
/*      CCelim_htnode **child)                                              */
/*    SEARCH for the Hamilton move in the child list of htnode ht           */
/*    -child returns node for the Hamilton move (or NULL if not in list)    */
/*                                                                          */
/*  int CCelim_write_httree (CCelim_httree **ht, char *fname,               */
/*      int append, int only_parents)                                       */
/*    WRITE the Hamilton-Tutte tree to file fname                           */
/*    -append (0 == new file, 1 = append)                                   */
/*    -only_parents (0 = full tree, 1 = just tree parents, 2 = edge file)   */
/*                                                                          */
/*  int CCelim_print_httree (CCelim_httree *ht, FILE *out,                  */
/*       int only_parents)                                                  */
/*    WRITE the H-T tree to the open file out (use stdout for printing)     */
/*                                                                          */
/*  int CCelim_swrite_httree (CCelim_httree *ht, CC_SFILE *s)               */
/*    WRITE the Hamilton-Tutte (full) tree to open binary file s            */
/*                                                                          */
/*  int CCelim_read_httree (CCelim_httree **pht, FILE *in, CC_SFILE *s)     */
/*    READ the next H-T tree in from the open file in or the open binary    */
/*      file s                                                              */
/*    -ht returns the next tree (or NULL if end of file)                    */
/*                                                                          */
/*  int CCelim_copy_httree (CCelim_httree *ht, CCelim_httree **pout)        */
/*    COPY the H-T tree from ht to pout                                     */
/*                                                                          */
/*  void CCelim_init_htnode (CCelim_htnode *n)                              */
/*    INITALIZE the Hamilton-Tutte node                                     */
/*                                                                          */
/*  void CCelim_free_htnode_children (CCelim_htnode *n)                     */
/*    FREE all nodes in Hamilton-Tutte subtree rooted at n                  */
/*                                                                          */
/*  void CCelim_init_httree (CCelim_httree *ht)                             */
/*    INITALIZE the Hamilton-Tutte tree structure                           */
/*                                                                          */
/*  void CCelim_free_httree (CCelim_httree *ht)                             */
/*    FREE the Hamiltion-Tutte tree                                         */
/*                                                                          */
/****************************************************************************/

#include <stdio.h>
#include "ccutil.h"
#include "elim.h"

static int print_htnode (CCelim_htnode *root, FILE *out, int parent_id,
    int only_print_edge);
static int swrite_htnode (CCelim_htnode *n, CC_SFILE *s, int parent_id);
static int copy_htnode (CCelim_htnode *n, CCelim_htnode **pcopy);

void CCelim_hamilton_find (CCelim_htnode *ht, int ham_type, int ham_path_len,
        int ham_nodes[CCelim_HTNODE_MAX_PATH_LEN], CCelim_htnode **child)
{
    CCelim_htnode *n;
    int *a = ham_nodes, *b = (int *) NULL, *c, *d, i;

    *child = (CCelim_htnode *) NULL;
    if (ham_type == 3) b = ham_nodes + 3;  /* a and b are the 3-paths */

    for (n = ht->child; n; n = n->next) {
        if (n->hamilton_type == ham_type) {
            if (ham_type == CCelim_HAMILTON_EDGE) {
                c = n->hamilton_nodes;
                if ((a[0] == c[0] && a[1] == c[1]) ||
                    (a[0] == c[1] && a[1] == c[0])) {
                    *child = n; goto CLEANUP;
                }
            } else if (ham_type == CCelim_HAMILTON_PATH) {
                c = n->hamilton_nodes;
                if (a[1] == c[1]) {  /* middle nodes of 3-paths agree */
                    if ((a[0] == c[0] && a[2] == c[2]) ||
                        (a[0] == c[2] && a[2] == c[0])) {
                        *child = n; goto CLEANUP;
                    }
                }
            } else if (ham_type == CCelim_HAMILTON_PAIR) {
                for (i = 0; i < 2; i++) {
                    if (i == 0) { c = n->hamilton_nodes; d = c+3; }
                    else        { d = n->hamilton_nodes; c = d+3; }

                    if (a[1] == c[1] && b[1] == d[1]) {
                        if (((a[0] == c[0] && a[2] == c[2]) ||
                             (a[0] == c[2] && a[2] == c[0])) &&
                            ((b[0] == d[0] && b[2] == d[2]) ||
                             (b[0] == d[2] && b[2] == d[0]))) {
                            *child = n; goto CLEANUP;
                        }
                    }
                }
            } else if (ham_type == CCelim_HAMILTON_LONG) {
                if (n->hamilton_path_len == ham_path_len) {
                    c = n->hamilton_nodes;
                    for (i = 0; i < ham_path_len; i++) {
                        if (a[i] != c[i]) break;
                    }
                    if (i == ham_path_len) { *child = n; goto CLEANUP; }
                    for (i = 0; i < ham_path_len; i++) {
                        if (a[i] != c[ham_path_len - i - 1]) break;
                    }
                    if (i == ham_path_len) { *child = n; goto CLEANUP; }
                }
            }
        }
    }

CLEANUP:
    return;
}

void CCelim_init_htnode (CCelim_htnode *n)
{
    if (n) {
        n->hamilton_type = CCelim_HAMILTON_NONE;
        n->hamilton_path_len = 0;
        n->tutte_type = CCelim_TUTTE_NONE;
        n->id = -1;
        n->child = (CCelim_htnode *) NULL;
        n->next = (CCelim_htnode *) NULL;
    }
}

void CCelim_free_htnode_children (CCelim_htnode *n)
{
    CCelim_htnode *c, *cnext = (CCelim_htnode *) NULL;

    if (n) {
        for (c = n->child; c; c = cnext) {
            cnext = c->next;
            CCelim_free_htnode_children (c);
            CC_FREE (c, CCelim_htnode);
        } 
        n->child = (CCelim_htnode *) NULL;
    }
}

void CCelim_init_httree (CCelim_httree *ht)
{
    if (ht) {
        ht->root = (CCelim_htnode *) NULL;
        ht->elimtype = CCelim_HTTREE_EDGE;
        ht->count = 0;
        ht->depth = 0;
        ht->max_count = 0;
        ht->max_depth = 0;
        ht->level = 0;
        ht->longpath = 0;
        ht->use_tsp = 0;
        ht->max_neighborhood = 0;
    }
}

void CCelim_free_httree (CCelim_httree *ht)
{
    if (ht) {
        CCelim_free_htnode_children (ht->root);
        CC_IFFREE (ht->root, CCelim_htnode);
        CCelim_init_httree (ht);
    }
}

int CCelim_write_httree (CCelim_httree *ht, char *fname, int append,
        int only_parents)
{
    int rval = 0;
    FILE *out = (FILE *) NULL;

    if (append) out = fopen (fname, "a");
    else        out = fopen (fname, "w");
    if (!out) {
        fprintf (stderr, "could not open %s for writing\n", fname);
        rval = 1; goto CLEANUP;
    }

    rval = CCelim_print_httree (ht, out, only_parents);
    CCcheck_rval (rval, "CCelim_print_httree failed");

CLEANUP:
    if (out) fclose (out);
    return rval;
}

int CCelim_swrite_httree (CCelim_httree *ht, CC_SFILE *s)
{
    int rval = 0, i;

    if (!ht || !ht->root ||
               (ht->root->hamilton_type != CCelim_HAMILTON_EDGE &&
                ht->root->hamilton_type != CCelim_HAMILTON_PATH)) {
        fprintf (stderr, "Hamilton-Tutte tree not initialized\n");
        rval = 1; goto CLEANUP;
    }

    rval = CCutil_swrite_int (s, ht->elimtype);
    CCcheck_rval (rval, "CCutil_swrite_int failed (elimtype)");
   
    for (i = 0; i < 2; i++) {
        rval = CCutil_swrite_int (s, ht->root->hamilton_nodes[i]);
        CCcheck_rval (rval, "CCutil_swrite_int_failed (nodes[i])");
    }
    if (ht->elimtype == CCelim_HTTREE_PAIR) {
        rval = CCutil_swrite_int (s, ht->root->hamilton_nodes[2]);
        CCcheck_rval (rval, "CCutil_swrite_int_failed (nodes[2])");
    }

    rval = CCutil_swrite_int (s, ht->count);
    CCcheck_rval (rval, "CCutil_swrite_int_failed (ht->count)");
    rval = CCutil_swrite_int (s, ht->depth);
    CCcheck_rval (rval, "CCutil_swrite_int_failed (ht->depth)");
    rval = CCutil_swrite_int (s, ht->max_count);
    CCcheck_rval (rval, "CCutil_swrite_int_failed (ht->max_count)");
    rval = CCutil_swrite_int (s, ht->max_depth);
    CCcheck_rval (rval, "CCutil_swrite_int_failed (ht->max_depth)");
    rval = CCutil_swrite_int (s, ht->level);
    CCcheck_rval (rval, "CCutil_swrite_int_failed (ht->level)");
    rval = CCutil_swrite_int (s, ht->longpath);
    CCcheck_rval (rval, "CCutil_swrite_int_failed (ht->longpath)");
    rval = CCutil_swrite_int (s, ht->use_tsp);
    CCcheck_rval (rval, "CCutil_swrite_int_failed (ht->use_tsp)");
    rval = CCutil_swrite_int (s, ht->max_neighborhood);
    CCcheck_rval (rval, "CCutil_swrite_int_failed (ht->max_neighborhood)");

    rval = swrite_htnode (ht->root, s, -1);
    CCcheck_rval (rval, "swrite_htnode failed");

CLEANUP:
    return rval;
}

int CCelim_copy_httree (CCelim_httree *ht, CCelim_httree **pout)
{
    int rval = 0;
    CCelim_httree *out = (CCelim_httree *) NULL;
    
    CC_MALLOC (out, 1, CCelim_httree);
    CCelim_init_httree (out);

    if (!ht || ht->count == 0 || !ht->root ||
               (ht->root->hamilton_type != CCelim_HAMILTON_EDGE &&
                ht->root->hamilton_type != CCelim_HAMILTON_PATH)) {
        fprintf (stderr, "Hamilton-Tutte tree not initialized\n");
        rval = 1; goto CLEANUP;
    }

    out->elimtype = ht->elimtype;
    out->count = ht->count;
    out->depth = ht->depth;
    out->max_count = ht->max_count;
    out->max_depth = ht->max_depth;
    out->level = ht->level;
    out->longpath = ht->longpath;
    out->use_tsp = ht->use_tsp;
    out->max_neighborhood = ht->max_neighborhood;

    rval = copy_htnode (ht->root, &out->root);
    CCcheck_rval (rval, "copy_htnode failed");
    if (!out->root) {
        fprintf (stderr, "No root node!\n"); exit (1);
    }
    *pout = out;

CLEANUP:
    return rval;
}

static int copy_htnode (CCelim_htnode *n, CCelim_htnode **pcopy)
{
    int rval = 0, i, hcount = 0;
    CCelim_htnode *copy, *c, *d, *last = (CCelim_htnode *) NULL;

    CC_MALLOC (copy, 1, CCelim_htnode);
    CCelim_init_htnode (copy);

    copy->hamilton_type = n->hamilton_type;
    copy->hamilton_path_len = n->hamilton_path_len;
    switch (n->hamilton_type) {
    case CCelim_HAMILTON_EDGE: hcount = 2; break;
    case CCelim_HAMILTON_PATH: hcount = 3; break;
    case CCelim_HAMILTON_PAIR: hcount = 6; break;
    case CCelim_HAMILTON_LONG: hcount = n->hamilton_path_len; break;
    }
    for (i = 0; i < hcount; i++) {
        copy->hamilton_nodes[i] = n->hamilton_nodes[i];
    }
    copy->tutte_type = n->tutte_type;
    if (n->tutte_type == CCelim_TUTTE_END ||
        n->tutte_type == CCelim_TUTTE_POINT) {
        copy->tutte_nodes[0] = n->tutte_nodes[0];
    } else if (n->tutte_type == CCelim_TUTTE_CD) {
        copy->tutte_nodes[0] = n->tutte_nodes[0];
        copy->tutte_nodes[1] = n->tutte_nodes[1];
    }
    copy->id = n->id;
    
    for (d = n->child; d; d = d->next) {
        rval = copy_htnode (d, &c);
        CCcheck_rval (rval, "copy_htnode failed");
        if (!last) {
            copy->child = c; last = c;
        } else {
            last->next = c; last = c;
        }
    }
    *pcopy = copy;

CLEANUP:
    return rval;
}

int CCelim_print_httree (CCelim_httree *ht, FILE *out, int only_parents)
{
    int rval = 0;

    if (!ht || !ht->root ||
               (ht->root->hamilton_type != CCelim_HAMILTON_EDGE &&
                ht->root->hamilton_type != CCelim_HAMILTON_PATH)) {
        fprintf (stderr, "Hamilton-Tutte tree not initialized\n");
        rval = 1; goto CLEANUP;
    }

    if (only_parents) {
        if (only_parents == 1) fprintf (out, "%d\n", ht->count);
        else                   fprintf (out, "%d %d\n", ht->count, ht->count-1);
    } else {
        switch (ht->elimtype) {
        case CCelim_HTTREE_EDGE: fprintf (out, "EDGE "); break;
        case CCelim_HTTREE_FIX:  fprintf (out, "FIX ");  break;
        case CCelim_HTTREE_PAIR: fprintf (out, "PAIR "); break;
        default:
            fprintf (stderr, "Unknown elimtype for tree\n");
            rval = 1; goto CLEANUP;
        }

        if (ht->elimtype == CCelim_HTTREE_PAIR) {
            fprintf (out, "%d %d %d\n", ht->root->hamilton_nodes[0],
                                        ht->root->hamilton_nodes[1],
                                        ht->root->hamilton_nodes[2]);
        } else {
            fprintf (out, "%d %d\n", ht->root->hamilton_nodes[0],
                                     ht->root->hamilton_nodes[1]);
        }
        fprintf (out, "%d %d %d %d %d %d %d %d\n", ht->count, ht->depth,
                   ht->max_count, ht->max_depth, ht->level, ht->longpath,
                   ht->use_tsp, ht->max_neighborhood);
    }

    rval = print_htnode (ht->root, out, -1, only_parents);
    CCcheck_rval (rval, "print_htnode failed");

    fflush (out);

CLEANUP:
    return rval;
}

static int print_htnode (CCelim_htnode *n, FILE *out, int parent_id,
        int only_print_edge)
{
    int rval = 0, i, hcount = 0;
    CCelim_htnode *c;

    if (n) {
        if (only_print_edge) {
            if (only_print_edge == 1) {
                fprintf (out, "%d %d\n", parent_id, n->id);
            } else if (parent_id != -1) {
                fprintf (out, "%d %d 1\n", parent_id, n->id);
            }
        } else {
            fprintf (out, "%d %d", parent_id, n->id);
            fprintf (out, " %d", n->hamilton_type);
            if (n->hamilton_type == CCelim_HAMILTON_LONG) {
                fprintf (out, " %d", n->hamilton_path_len);
            }
            switch (n->hamilton_type) {
            case CCelim_HAMILTON_EDGE: hcount = 2; break;
            case CCelim_HAMILTON_PATH: hcount = 3; break;
            case CCelim_HAMILTON_PAIR: hcount = 6; break;
            case CCelim_HAMILTON_LONG: hcount = n->hamilton_path_len; break;
            }
            for (i = 0; i < hcount; i++) {
                fprintf (out, " %d", n->hamilton_nodes[i]);
            }  
            fprintf (out, " %d", n->tutte_type);
            if (n->tutte_type == CCelim_TUTTE_END ||
                n->tutte_type == CCelim_TUTTE_POINT) {
                fprintf (out, " %d", n->tutte_nodes[0]);
            } else if (n->tutte_type == CCelim_TUTTE_CD) {
                fprintf (out, " %d %d", n->tutte_nodes[0], n->tutte_nodes[1]);
            }
            fprintf (out, "\n");
        }

        for (c = n->child; c; c = c->next) {
            rval = print_htnode (c, out, n->id, only_print_edge);
            CCcheck_rval (rval, "print_htnode failed");
        }
    }

CLEANUP:
    return rval;
}

static int swrite_htnode (CCelim_htnode *n, CC_SFILE *s, int parent_id)
{
    int rval = 0, i, hcount = 0;
    CCelim_htnode *c;

    if (!n) {
        fprintf (stderr, "trying to swrite an empty htnode\n");
        rval = 1; goto CLEANUP;
    }

    rval = CCutil_swrite_int (s, parent_id);
    CCcheck_rval (rval, "CCutil_swrite_int failed (parent_id)");
    rval = CCutil_swrite_int (s, n->id);
    CCcheck_rval (rval, "CCutil_swrite_int failed (n->id)");
    rval = CCutil_swrite_int (s, n->hamilton_type);
    CCcheck_rval (rval, "CCutil_swrite_int failed (hamilton_type)");
    if (n->hamilton_type == CCelim_HAMILTON_LONG) {
        rval = CCutil_swrite_int (s, n->hamilton_path_len);
        CCcheck_rval (rval, "CCutil_swrite_int failed (hamilton_path_len)");
    }

    switch (n->hamilton_type) {
    case CCelim_HAMILTON_EDGE: hcount = 2; break;
    case CCelim_HAMILTON_PATH: hcount = 3; break;
    case CCelim_HAMILTON_PAIR: hcount = 6; break;
    case CCelim_HAMILTON_LONG: hcount = n->hamilton_path_len; break;
    }

    for (i = 0; i < hcount; i++) {
        rval  = CCutil_swrite_int (s, n->hamilton_nodes[i]);
        CCcheck_rval (rval, "CCutil_swrite_int failed (hamilton_nodes[i])");
    }  

    rval  = CCutil_swrite_int (s, n->tutte_type);
    CCcheck_rval (rval, "CCutil_swrite_int failed (tutte_type)");

    if (n->tutte_type==CCelim_TUTTE_END || n->tutte_type==CCelim_TUTTE_POINT) {
        rval  = CCutil_swrite_int (s, n->tutte_nodes[0]);
        CCcheck_rval (rval, "CCutil_swrite_int failed (tutte_nodes[0])");
    } else if (n->tutte_type==CCelim_TUTTE_CD) {
        for (i = 0; i < 2; i++) {
            rval  = CCutil_swrite_int (s, n->tutte_nodes[i]);
            CCcheck_rval (rval, "CCutil_swrite_int failed (tutte_nodes[i])");
        }
    }

    for (c = n->child; c; c = c->next) {
        rval = swrite_htnode (c, s, n->id);
        CCcheck_rval (rval, "swrite_htnode failed");
    }

CLEANUP:
    return rval;
}

int CCelim_read_httree (CCelim_httree **pht, FILE *in, CC_SFILE *s)
{
    int rval = 0, i, k, parent_id, hcount, elimtype, topnodes[3];
    CCelim_htnode **htnodes = (CCelim_htnode **) NULL, *n, *p;
    CCelim_httree *ht = (CCelim_httree *) NULL;
    char buf[1024], buftype[1024], *q;
   
    *pht = (CCelim_httree *) NULL;

    if ((in && s) || (!in && !s)) {
        fprintf (stderr, "must be called with exactly one open file\n");
        rval = 1; goto CLEANUP;
    }

    if (in) {
        if (!fgets(buf, sizeof(buf), in)) {
            printf ("\n"); fflush (stdout);
            goto CLEANUP;
        }
    } else {
        if (CCutil_sread_int (s, &elimtype)) {
            printf ("\n"); fflush (stdout);
            goto CLEANUP;
        }
    }

    CC_MALLOC (ht, 1, CCelim_httree);
    CCelim_init_httree (ht);

    if (in) {
        q = buf; i = 0;
        while (*q != ' ') {
            if (*q == '\0') {
                fprintf (stderr, "elimtype line has error %s\n", buf);
                rval = 1; goto CLEANUP;
            }
            buftype[i++] = *(q++);
        }
        buftype[i] = '\0';

        if      (strcmp(buftype, "EDGE") == 0) elimtype = CCelim_HTTREE_EDGE;
        else if (strcmp(buftype, "FIX") == 0)  elimtype = CCelim_HTTREE_FIX;
        else if (strcmp(buftype, "PAIR") == 0) elimtype = CCelim_HTTREE_PAIR;
        else {
            fprintf (stderr, "unknown elimtype: %s\n", buftype);
            rval = 1; goto CLEANUP;
        }
    }
    ht->elimtype = elimtype;

    if (in) {
        if (elimtype == CCelim_HTTREE_EDGE || elimtype == CCelim_HTTREE_FIX) {
            if (sscanf (q, "%d %d\n", &topnodes[0], &topnodes[1]) != 2) {
                fprintf (stderr, "missing data in elimtype line %s\n", buf);
                rval = 1; goto CLEANUP;
            }
        } else {
            if (sscanf (q, "%d %d %d\n", &topnodes[0], &topnodes[1],
                                                       &topnodes[2]) != 3) {
                fprintf (stderr, "missing data in elimtype line %s\n", buf);
                rval = 1; goto CLEANUP;
            }
        }
    } else {
        for (i = 0; i < 2; i++) {
            rval = CCutil_sread_int (s, &topnodes[i]);
            CCcheck_rval (rval, "CCutil_sread_failed (topnodes[i])");
        }
        if (ht->elimtype == CCelim_HTTREE_PAIR) {
            rval = CCutil_sread_int (s, &topnodes[2]);
            CCcheck_rval (rval, "CCutil_sread_int_failed (topnodes[2])");
        }
    }

    /* note: topnodes are not needed, they are also in the root htnode */

    if (in) {
        fscanf (in, "%d %d %d %d %d %d %d %d\n", &ht->count, &ht->depth,
                   &ht->max_count, &ht->max_depth, &ht->level, &ht->longpath,
                   &ht->use_tsp, &ht->max_neighborhood);
    } else {
        rval = CCutil_sread_int (s, &ht->count);
        CCcheck_rval (rval, "CCutil_sread_int_failed (ht->count)");
        rval = CCutil_sread_int (s, &ht->depth);
        CCcheck_rval (rval, "CCutil_sread_int_failed (ht->depth)");
        rval = CCutil_sread_int (s, &ht->max_count);
        CCcheck_rval (rval, "CCutil_sread_int_failed (ht->max_count)");
        rval = CCutil_sread_int (s, &ht->max_depth);
        CCcheck_rval (rval, "CCutil_sread_int_failed (ht->max_depth)");
        rval = CCutil_sread_int (s, &ht->level);
        CCcheck_rval (rval, "CCutil_sread_int_failed (ht->level)");
        rval = CCutil_sread_int (s, &ht->longpath);
        CCcheck_rval (rval, "CCutil_sread_int_failed (ht->longpath)");
        rval = CCutil_sread_int (s, &ht->use_tsp);
        CCcheck_rval (rval, "CCutil_sread_int_failed (ht->use_tsp)");
        rval = CCutil_sread_int (s, &ht->max_neighborhood);
        CCcheck_rval (rval, "CCutil_sread_int_failed (ht->max_neighborhood)");
    }

    if (ht->count == 0) {
        fprintf (stderr, "Hamilton-Tutte tree has no nodes\n");
        rval = 1; goto CLEANUP;
    }

    CC_MALLOC (htnodes, ht->count, CCelim_htnode *);
    for (i = 0; i < ht->count; i++) htnodes[i] = (CCelim_htnode *) NULL;

    for (k = 0; k < ht->count; k++) {
        CC_MALLOC (n, 1, CCelim_htnode);
        CCelim_init_htnode (n);

        if (in) {
            fscanf (in, "%d %d", &parent_id, &n->id);
        } else {
            rval = CCutil_sread_int (s, &parent_id);
            CCcheck_rval (rval, "CCutil_sread_int failed (parent_id)");
            rval = CCutil_sread_int (s, &n->id);
            CCcheck_rval (rval, "CCutil_sread_int failed (n->id)");
        }

        if (n->id < 0 || n->id >= ht->count) {
            fprintf (stderr, "Node has bad id: %d\n", n->id);
            rval = 1; goto CLEANUP;
        }
        if (parent_id < -1 || parent_id >= ht->count) {
            fprintf (stderr, "Node has bad parent id: %d\n", n->id);
            rval = 1; goto CLEANUP;
        }
        if (parent_id == -1 && ht->root) {
            fprintf (stderr, "Tree has two root nodes: %d\n", n->id);
            rval = 1; goto CLEANUP;
        }
        if (htnodes[n->id]) {
            fprintf (stderr, "Node reused: %d\n", n->id);
            rval = 1; goto CLEANUP;
        }

        htnodes[n->id] = n;

        if (parent_id != -1) {
            if (!htnodes[parent_id]) {
                fprintf (stderr, "Node missing parent (check edge order)\n");
                rval = 1; goto CLEANUP;
            }
            p = htnodes[parent_id];
            n->next = p->child;
            p->child = n;
        } else {
            ht->root = n;
        }
        
        if (in) {
            fscanf (in, " %d", &n->hamilton_type);
        } else {
            rval = CCutil_sread_int (s, &n->hamilton_type);
            CCcheck_rval (rval, "CCutil_sread_int failed (hamilton_type)");
        }

        if (n->hamilton_type == CCelim_HAMILTON_LONG) {
            if (in) {
                fscanf (in, " %d", &n->hamilton_path_len);
            } else {
                rval = CCutil_sread_int (s, &n->hamilton_path_len);
                CCcheck_rval (rval, "CCutil_sread_int failed (path_len)");
            }
            if (n->hamilton_path_len > CCelim_HTNODE_MAX_PATH_LEN) {
                fprintf (stderr, "Node path len exceeds maximum\n");
                rval = 1; goto CLEANUP;
            }
        }
        switch (n->hamilton_type) {
        case CCelim_HAMILTON_EDGE: hcount = 2; break;
        case CCelim_HAMILTON_PATH: hcount = 3; break;
        case CCelim_HAMILTON_PAIR: hcount = 6; break;
        case CCelim_HAMILTON_LONG: hcount = n->hamilton_path_len; break;
        default:
            fprintf (stderr, "Node has bad hamilton_type\n");
            rval = 1; goto CLEANUP;
        }
        for (i = 0; i < hcount; i++) {
            if (in) {
                fscanf (in, " %d", &n->hamilton_nodes[i]);
            } else {
                rval = CCutil_sread_int (s, &n->hamilton_nodes[i]);
                CCcheck_rval (rval, "CCutil_sread_int failed (nodes[i])");
            }
        }  
        if (in) {
             fscanf (in, " %d", &n->tutte_type);
             if (n->tutte_type == CCelim_TUTTE_PATH || 
                 n->tutte_type == CCelim_TUTTE_NONE) {
                 fscanf (in, "\n"); /* to move to new line */
             } else if (n->tutte_type == CCelim_TUTTE_END ||
                        n->tutte_type == CCelim_TUTTE_POINT) {
                 fscanf (in, " %d\n", &n->tutte_nodes[0]);
             } else if (n->tutte_type == CCelim_TUTTE_CD) {
                 fscanf (in, " %d %d\n", &n->tutte_nodes[0],
                                         &n->tutte_nodes[1]);
             } else if (n->tutte_type != CCelim_TUTTE_NONE) {
                 fprintf (stderr, "Node has bad tutte_type\n");
                 rval = 1; goto CLEANUP;
             }
        } else {
            rval = CCutil_sread_int (s, &n->tutte_type);
            CCcheck_rval (rval, "CCutil_sread_int failed (tutte_type)");
            if (n->tutte_type==CCelim_TUTTE_END ||
                n->tutte_type==CCelim_TUTTE_POINT) {
                rval  = CCutil_sread_int (s, &n->tutte_nodes[0]);
                CCcheck_rval (rval, "CCutil_sread_int failed (nodes[0])");
            } else if (n->tutte_type==CCelim_TUTTE_CD) {
                for (i = 0; i < 2; i++) {
                    rval  = CCutil_sread_int (s, &n->tutte_nodes[i]);
                    CCcheck_rval (rval, "CCutil_sread_int failed (nodes[i])");
                }
            }
        }
    }

    *pht = ht;

CLEANUP:
    if (rval && htnodes) {
        for (i = 0; i < ht->count; i++) {
            CC_IFFREE (htnodes[i], CCelim_htnode);
        }
    }
    CC_IFFREE (htnodes, CCelim_htnode *);
    return rval;
}
