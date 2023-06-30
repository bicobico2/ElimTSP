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
/*  Written by:  May/June, 2016                                             */
/*  Last Update: June 30, 2016                                              */
/*                                                                          */
/*  Exported Function                                                       */
/*                                                                          */
/*  int CCelim_split_cd (int a, int b, int c, int d, int ccount,            */
/*      int *clist, int dcount, int *dlist, CCelim_distobj *D, int *good)   */
/*    ELIMINATE edge (a,b); using witness c and d; clist/dlist contain all  */
/*    all valid pairs of neighbors of c/d; good set to 1 if successful      */
/*                                                                          */
/*  NOTE: Uses speacilized code that should improve speed in fast elim. The */
/*    function elim_with_cd can be used instead, calling the standard test  */
/*    CCelim_test_paths().                                                  */
/*                                                                          */
/****************************************************************************/

#include <stdio.h>
#include "ccutil.h"
#include "elim.h"

static void split_list (int a, int b, int c, int ccount, int *clist,
    int *clist_2_count, int *clist_2, int *clist_a_count, int *clist_a,
    int *clist_b_count, int *clist_b);
static void check_2_3_3 (int *nodeset, int **M, int *I, CCelim_distobj *D,
    int *yesno);
static void check_2_5 (int *nodeset, int **M, int *I, CCelim_distobj *D,
    int *yesno);
static void run_3_4 (int count3, int *list3, int count4, int *list4, int **M,
    int *I, int *nodeset, CCelim_distobj *D, int *yesno);
static void run_4_4 (int count4c, int *list4c, int count4d, int *list4d,
    int **M, int *I, int *nodeset, CCelim_distobj *D, int *yesno);
static void build_M (int count, int *nodeset, int **M, CCelim_distobj *D);
static void swapper5 (int n0, int n1, int n2, int n3, int n4, int n5, int n6,
    int n7, int n8, int n9, int **M, int *I, int *yesno);
static void next_set (int sz, int *Set);

/***************************************************************************/
/*                                                                         */
/*  To eliminate ab, a witness pair (c,d) is chosen such that c-d is not   */
/*  compatiable with a-b. The function check_neighbor_three_swap returns   */
/*  lists of triples (c0, c, c1) that can possibly be in an optimal tour   */
/*  with a-b. Edges c-c0 and c-c1 are compatiable with a-b and the global  */
/*  three-swap is not improving.  Similarly for (d0, d, d1).               */
/*                                                                         */
/*  We test here if, for each pair of triples (c0, c, c1) and (d0, d, d1), */
/*  every tour orientation of the three paths a-b, c0-c-c1, and d0-d-d1    */
/*  yields an improving 3, 4, or 5-swap. The aim is to take advantage of   */
/*  the small structure to obtain a faster run time than calling the       */
/*  general function CCelim_test_paths ().                                 */
/*                                                                         */
/*  In these comments, the paths a-b, c0-c-c1, d0-d-d1 are recorded as     */
/*  0-1, 2-3-4, 5-6-7.                                                     */
/*                                                                         */
/*  The first step is to test if any of the paths are linked, that is,     */
/*  the end nodes are not distinct.                                        */
/*                                                                         */
/*  Case A: All nodes are distinct.  We have three paths, yielding 8 tour  */
/*          orientations.                                                  */
/*   TA1: 0-1  2-3-4  5-6-7                                                */
/*   TA2: 0-1  2-3-4  7-6-5                                                */
/*   TA3: 0-1  4-3-2  5-6-7                                                */
/*   TA4: 0-1  4-3-2  7-6-5                                                */
/*   TA5: 0-1  5-6-7  2-3-4                                                */
/*   TA6: 0-1  5-6-7  4-3-2                                                */
/*   TA7: 0-1  7-6-5  2-3-4                                                */
/*   TA8: 0-1  7-6-5  4-3-2                                                */
/*                                                                         */
/*  Case B: We have two paths, yielding 2 tour orientations. There are     */
/*          12 ways this can happen                                        */
/*   Long path includes 0-1                                                */
/*   B1:  0 == 2   B1a:  1-[0,2]-3-4  5-6-7    B1b:  1-[0,2]-3-4  7-6-5    */
/*   B2:  0 == 4   B2a:  1-[0,4]-3-2  5-6-7    B2b:  1-[0,4]-3-2  7-6-5    */
/*   B3:  0 == 5   B3a:  1-[0,5]-6-7  2-3-4    B3b:  1-[0,5]-6-7  4-3-2    */
/*   B4:  0 == 7   B4a:  1-[0,7]-6-5  2-3-4    B4b:  1-[0,7]-6-5  4-3-2    */
/*   B5:  1 == 2   B5a:  0-[1,2]-3-4  5-6-7    B5b:  0-[1,2]-3-4  7-6-5    */
/*   B6:  1 == 4   B6a:  0-[1,4]-2-3  5-6-7    B6b:  0-[1,4]-2-3  7-6-5    */
/*   B7:  1 == 5   B7a:  0-[1,5]-6-7  2-3-4    B7b:  0-[1,5]-6-7  4-3-2    */
/*   B8:  1 == 7   B8a:  0-[1,7]-6-5  2-3-4    B8b   0-[1,7]-6-5  4-3-2    */
/*   Long path does not include 0-1                                        */
/*   B9:  2 == 5   B9a:  0-1  4-3-[2,5]-6-7    B9b:  1-0  4-3-[2,5]-6-7    */
/*   B10: 2 == 7   B10a: 0-1  4-3-[2,7]-6-5    B10b: 1-0  4-3-[2,7]-6-5    */
/*   B11: 4 == 5   B11a: 0-1  2-3-[4,5]-6-7    B11b: 1-0  2-3-[4,5]-6-7    */
/*   B12: 4 == 7   B12a: 0-1  2-3-[4,7]-6-5    B12b: 1-0  2-3-[4,7]-6-5    */
/*                                                                         */
/*   Case C: We have one path, yielding 1 tour orientation. There are 24   */
/*           way this can happen.                                          */
/*     NOTE: We must check that 2-3-4 and 5-6-7 do not form a circuit.     */
/*           That is, 2 != 5 || 2 != 7 || 4 != 5 || 4 != 7.                */
/*     Orientation: Fix direction to be 0-1                                */
/*   Path starts with 0-1                                                  */
/*   C1:  1 == 2 && 4 == 5   TC1:  0-[1,2]-3-[4,5]-6-7                     */
/*   C2:  1 == 2 && 4 == 7   TC2:  0-[1,2]-3-[4,7]-6-5                     */
/*   C3:  1 == 4 && 2 == 5   TC3:  0-[1,4]-3-[2,5]-6-7                     */
/*   C4:  1 == 4 && 2 == 7   TC4:  0-[1,4]-3-[2,7]-6-5                     */
/*   C5:  1 == 5 && 7 == 2   TC5:  0-[1,5]-6-[7,2]-3-4                     */
/*   C6:  1 == 5 && 7 == 4   TC6:  0-[1,5]-6-[7,4]-3-2                     */
/*   C7:  1 == 7 && 5 == 2   TC7:  0-[1,7]-6-[5,2]-3-4                     */
/*   C8:  1 == 7 && 5 == 4   TC8:  0-[1,7]-6-[5,4]-3-2                     */
/*   Path ends with 0-1                                                    */
/*   C9:  4 == 5 && 7 == 0   TC9:  2-3-[4,5]-6-[7,0]-1                     */
/*   C10: 4 == 7 && 5 == 0   TC10: 2-3-[4,7]-6-[5,0]-1                     */
/*   C11: 2 == 5 && 7 == 0   TC11: 4-3-[2,5]-6-[7,0]-1                     */
/*   C12: 2 == 7 && 5 == 0   TC12: 4-3-[2,7]-6-[5,0]-1                     */
/*   C13: 7 == 2 && 4 == 0   TC13: 5-6-[7,2]-3-[4,0]-1                     */
/*   C14: 7 == 4 && 2 == 0   TC14: 5-6-[7,4]-3-[2,0]-1                     */
/*   C15: 5 == 2 && 4 == 0   TC15: 7-6-[5,2]-3-[4,0]-1                     */
/*   C16: 5 == 4 && 2 == 0   TC16: 7-6-[5,4]-3-[2,0]-1                     */
/*   Path has 0-1 in the middle                                            */
/*   C17: 4 == 0 && 1 == 5   TC17: 2-3-[4,0]-[1,5]-6-7                     */
/*   C18: 4 == 0 && 1 == 7   TC18: 2-3-[4,0]-[1,7]-6-5                     */
/*   C19: 2 == 0 && 1 == 5   TC19: 4-3-[2,0]-[1,5]-6-7                     */
/*   C20: 2 == 0 && 1 == 7   TC20: 4-3-[2,0]-[1,7]-6-5                     */
/*   C21: 7 == 0 && 1 == 2   TC21: 5-6-[7,0]-[1,2]-3-4                     */
/*   C22: 7 == 0 && 1 == 4   TC22: 5-6-[7,0]-[1,4]-3-2                     */
/*   C23: 5 == 0 && 1 == 2   TC23: 7-6-[5,0]-[1,2]-3-4                     */
/*   C24: 5 == 0 && 1 == 4   TC24: 7-6-[5,0]-[1,4]-3-2                     */
/*                                                                         */
/* To cut down on the number of cases, we preprocess the c+d triples to    */
/* to determine their intersection with edge a-b. This will split clist    */
/* into three lists: clist_2 (two paths), clist_a (single path with a as   */
/* an end node, and clist_b (single path with b as an end node). Similarly */
/* we have dlist_2, dlist_a, and dlist_b                                   */
/*                                                                         */
/* To evaluate ab and the c+d pair of triples, we use three loops          */
/*   Loop 1: clist_2 and dlist_2 => three paths or ab and c+d path         */
/*   Loop 2: clist_2 and dlist_a => single path or two paths               */
/*   Loop 3: clist_2 and dlist_b => single path or two paths               */
/*   Loop 4: clist_a and dlist_a => single path,                           */
/*   Loop 5: clist_b and dlist_b => single path                            */
/*                                                                         */
/***************************************************************************/

CC_ELIM_DIST_IN   /* macro builds local copy dist() of CCelim_dist() */

int CCelim_split_cd (int a, int b, int c, int d, int ccount, int *clist,
        int dcount, int *dlist, CCelim_distobj *D, int *good)
{
    int rval = 0, i, j, c0, c1, d0, d1, yesno = 1, I[10];
    int clist_2_count = 0, clist_a_count = 0, clist_b_count = 0;
    int dlist_2_count = 0, dlist_a_count = 0, dlist_b_count = 0;
    int *clist_2 = (int *) NULL, *clist_a = (int *) NULL;
    int *clist_b = (int *) NULL, *dlist_2 = (int *) NULL;
    int *dlist_a = (int *) NULL, *dlist_b = (int *) NULL;
    int nodeset[8], **M = (int **) NULL;

    *good = 0;

    if (ccount == 0 || dcount == 0) { *good = 1; goto CLEANUP; }

    CC_MALLOC (clist_2, 3*ccount, int);
    CC_MALLOC (clist_a, 4*ccount, int); CC_MALLOC (clist_b, 4*ccount, int);
    CC_MALLOC (dlist_2, 3*dcount, int);
    CC_MALLOC (dlist_a, 4*dcount, int); CC_MALLOC (dlist_b, 4*dcount, int);

    split_list (a, b, c, ccount, clist, &clist_2_count, clist_2, &clist_a_count,
                clist_a, &clist_b_count, clist_b);
    split_list (a, b, d, dcount, dlist, &dlist_2_count, dlist_2, &dlist_a_count,
                dlist_a, &dlist_b_count, dlist_b);

    CC_MALLOC (M, 8, int *);
    for (i = 0; i < 8; i++) { M[i] = (int *) NULL; }
    for (i = 0; i < 8; i++) { CC_MALLOC (M[i], 8, int); }

    /* LOOP 1 */

    for (i = 0; i < clist_2_count; i++) {
        nodeset[0] = a; nodeset[1] = b;
        c0 = clist_2[3*i]; c1 = clist_2[3*i+2];
        for (j = 0; j < dlist_2_count; j++) {
            d0 = dlist_2[3*j]; d1 = dlist_2[3*j+2];
            if ((c0 == d0 && c1 == d1) || (c0 == d1 && c1 == d0)) 
                continue;          /* circuit */
            if (c0 == d0) {        /* a-b, c1-c-c0-d-d1 */
                nodeset[2] = c1; nodeset[3] = c;  nodeset[4] = c0;
                nodeset[5] = d;  nodeset[6] = d1;
                check_2_5 (nodeset, M, I, D, &yesno);
            } else if (c0 == d1) { /* a-b, c1-c-c0-d-d0 */
                nodeset[2] = c1; nodeset[3] = c;  nodeset[4] = c0;
                nodeset[5] = d;  nodeset[6] = d0;
                check_2_5 (nodeset, M, I, D, &yesno);
            } else if (c1 == d0) { /* a-b, c0-c-c1-d-d1 */
                nodeset[2] = c0; nodeset[3] = c;  nodeset[4] = c1;
                nodeset[5] = d;  nodeset[6] = d1;
                check_2_5 (nodeset, M, I, D, &yesno);
            } else if (c1 == d1) { /* a-b, c0-c-c1-d-d0 */
                nodeset[2] = c0; nodeset[3] = c;  nodeset[4] = c1;
                nodeset[5] = d;  nodeset[6] = d0;
                check_2_5 (nodeset, M, I, D, &yesno);
            } else {               /* a-b, c0-c-c1, d0-d-d1 */
                nodeset[2] = c0; nodeset[3] = c;  nodeset[4] = c1;
                nodeset[5] = d0; nodeset[6] = d;  nodeset[7] = d1;
                check_2_3_3 (nodeset, M, I, D, &yesno);
            }
            if (yesno == 0) goto CLEANUP;
        }
    }

    /* LOOP 2 */

    run_3_4 (clist_2_count, clist_2, dlist_a_count, dlist_a,
             M, I, nodeset, D, &yesno);
    if (yesno == 0) goto CLEANUP;

    run_3_4 (clist_2_count, clist_2, dlist_b_count, dlist_b,
             M, I, nodeset, D, &yesno);
    if (yesno == 0) goto CLEANUP;

    /* LOOP 3 */

    run_3_4 (dlist_2_count, dlist_2, clist_a_count, clist_a,
             M, I, nodeset, D, &yesno);
    if (yesno == 0) goto CLEANUP;

    run_3_4 (dlist_2_count, dlist_2, clist_b_count, clist_b,
             M, I, nodeset, D, &yesno);
    if (yesno == 0) goto CLEANUP;

    /* LOOP 4 */

    run_4_4 (clist_a_count, clist_a, dlist_b_count, dlist_b,
             M, I, nodeset, D, &yesno);
    if (yesno == 0) goto CLEANUP;

    /* LOOP 5 */

    run_4_4 (clist_b_count, clist_b, dlist_a_count, dlist_a,
             M, I, nodeset, D, &yesno);
    if (yesno == 0) goto CLEANUP;

    *good = 1;

CLEANUP:
    CC_IFFREE (clist_2, int);
    CC_IFFREE (clist_a, int);
    CC_IFFREE (clist_b, int);
    CC_IFFREE (dlist_2, int);
    CC_IFFREE (dlist_a, int);
    CC_IFFREE (dlist_b, int);
    if (M) {
        for (i = 0; i < 8; i++) CC_IFFREE (M[i], int);
        CC_IFFREE (M, int *);
    }
    return rval;
}

static void run_3_4 (int count3, int *list3, int count4, int *list4, int **M,
        int *I, int *nodeset, CCelim_distobj *D, int *yesno)
{
    /* list3 has c0-c-c1, list4 has d0-d-d1-n (where n = a or b) */
    /* we know that c0, c, and c1 do not meet a-b                */
    int i, j, c0, c, c1;

    *yesno = 1;

    for (i = 0; i < count3; i++) {
        c0 = list3[3*i];  c = list3[3*i+1], c1 = list3[3*i+2];
        for (j = 0; j < count4; j++) {
            nodeset[0] = list4[4*j+3]; 
            nodeset[1] = list4[4*j+2]; 
            nodeset[2] = list4[4*j+1]; 
            nodeset[3] = list4[4*j]; 
            if (nodeset[3] == c0) {           /* One path */
                nodeset[4] = c;
                nodeset[5] = c1;
                build_M (6, nodeset, M, D);
                swapper5 (0, 1, 1, 2, 2, 3, 3, 4, 4, 5, M, I, yesno);
                if (*yesno == 0) return;
            } else if (nodeset[3] == c1) {    /* One path */
                nodeset[4] = c;
                nodeset[5] = c0;
                build_M (6, nodeset, M, D);
                swapper5 (0, 1, 1, 2, 2, 3, 3, 4, 4, 5, M, I, yesno);
                if (*yesno == 0) return;
            } else {                          /* Two paths */
                nodeset[4] = c0;
                nodeset[5] = c;
                nodeset[6] = c1;
                build_M (7, nodeset, M, D);
                swapper5 (0, 1, 1, 2, 2, 3, 4, 5, 5, 6, M, I, yesno);
                if (*yesno == 0) return;
                swapper5 (0, 1, 1, 2, 2, 3, 6, 5, 5, 4, M, I, yesno);
                if (*yesno == 0) return;
            }
        }
    }
}

static void run_4_4 (int count4c, int *list4c, int count4d, int *list4d,
        int **M, int *I, int *nodeset, CCelim_distobj *D, int *yesno)
{
    /* path 1 ends with ab and path 2 ends with ba (or vice versa) */
    int i, j;

    *yesno = 1;

    for (i = 0; i < count4c; i++) {
        nodeset[0] = list4c[4*i];
        nodeset[1] = list4c[4*i+1];
        nodeset[2] = list4c[4*i+2];
        nodeset[3] = list4c[4*i+3];
        for (j = 0; j < count4d; j++) {
            nodeset[4] = list4d[4*j+1];
            nodeset[5] = list4d[4*j];
            if (nodeset[0] == nodeset[5]) continue;  /* circuit */
            build_M (6, nodeset, M, D);
            swapper5 (0, 1, 1, 2, 2, 3, 3, 4, 4, 5, M, I, yesno);
            if (*yesno == 0) return;
        }
    }
}

static void check_2_3_3 (int *nodeset, int **M, int *I, CCelim_distobj *D,
        int *yesno)
{
    /* 2-node path, 3-node path, 3-node path */
    build_M (8, nodeset, M, D);
    swapper5 (0, 1, 2, 3, 3, 4, 5, 6, 6, 7, M, I, yesno);
    if (*yesno == 0) return;
    swapper5 (0, 1, 2, 3, 3, 4, 7, 6, 6, 5, M, I, yesno);
    if (*yesno == 0) return;
    swapper5 (0, 1, 4, 3, 3, 2, 5, 6, 6, 7, M, I, yesno);
    if (*yesno == 0) return;
    swapper5 (0, 1, 4, 3, 3, 2, 7, 6, 6, 5, M, I, yesno);
    if (*yesno == 0) return;
    swapper5 (0, 1, 5, 6, 6, 7, 2, 3, 3, 4, M, I, yesno);
    if (*yesno == 0) return;
    swapper5 (0, 1, 5, 6, 6, 7, 4, 3, 3, 2, M, I, yesno);
    if (*yesno == 0) return;
    swapper5 (0, 1, 7, 6, 6, 5, 2, 3, 3, 4, M, I, yesno);
    if (*yesno == 0) return;
    swapper5 (0, 1, 7, 6, 6, 5, 4, 3, 3, 2, M, I, yesno);
}

static void check_2_5 (int *nodeset, int **M, int *I, CCelim_distobj *D,
        int *yesno)
{
    /* 2-node path and 5-node path */
    build_M (7, nodeset, M, D);
    swapper5 (0, 1, 2, 3, 3, 4, 4, 5, 5, 6, M, I, yesno);
    if (*yesno == 0) return;
    swapper5 (1, 0, 2, 3, 3, 4, 4, 5, 5, 6, M, I, yesno);
}

static void split_list (int a, int b, int c, int ccount, int *clist,
       int *clist_2_count, int *clist_2, int *clist_a_count, int *clist_a,
       int *clist_b_count, int *clist_b)
{
    int i, c0, c1, x_count = 0, a_count = 0, b_count = 0;

    for (i = 0; i < ccount; i++) {
        c0 = clist[2*i], c1 = clist[2*i+1];
        if (c0 == a) {
            clist_b[4*b_count]   = c1; clist_b[4*b_count+1] = c;
            clist_b[4*b_count+2] = c0; clist_b[4*b_count+3] = b;
            b_count++;
        } else if (c1 == a) {
            clist_b[4*b_count]   = c0; clist_b[4*b_count+1] = c;
            clist_b[4*b_count+2] = c1; clist_b[4*b_count+3] = b;
            b_count++;
        } else if (c0 == b) {
            clist_a[4*a_count]   = c1; clist_a[4*a_count+1] = c;
            clist_a[4*a_count+2] = c0; clist_a[4*a_count+3] = a;
            a_count++;
        } else if (c1 == b) {
            clist_a[4*a_count]   = c0; clist_a[4*a_count+1] = c;
            clist_a[4*a_count+2] = c1; clist_a[4*a_count+3] = a;
            a_count++;
        } else {
            clist_2[3*x_count]   = c0; clist_2[3*x_count+1] = c;
            clist_2[3*x_count+2] = c1;
            x_count++;
        }
    }
    *clist_2_count = x_count;
    *clist_a_count = a_count;
    *clist_b_count = b_count;
}

static void build_M (int count, int *nodeset, int **M, CCelim_distobj *D)
{
    int i, j;

    for (i = 0; i < count; i++) {
        for (j = i+1; j < count; j++) {
            M[i][j] = M[j][i] = dist(nodeset[i], nodeset[j], D);
        }
        M[i][i] = 0;
    }
}

static void swapper5 (int n0, int n1, int n2, int n3, int n4, int n5, int n6,
        int n7, int n8, int n9, int **M, int *I, int *yesno)
{
    int i, k, Set[5], h[10], del[10];

    h[0] = n0; h[1] = n1; h[2] = n2; h[3] = n3; h[4] = n4;
    h[5] = n5; h[6] = n6; h[7] = n7; h[8] = n8; h[9] = n9;

    *yesno = 0;

    /* handle 3-swaps first */

    k = 3;
    for (i = 0; i < k; i++) Set[i] = i;
    for (; Set[k-1] < 5 && *yesno == 0; next_set (k, Set)) {
        int a, b, c, d, e, f, r;
        a = h[2*Set[0]]; b = h[2*Set[0]+1];
        c = h[2*Set[1]]; d = h[2*Set[1]+1];
        e = h[2*Set[2]]; f = h[2*Set[2]+1];
        r =  M[a][b] + M[c][d] + M[e][f];
        if (M[c][e] + M[d][a] + M[f][b] < r ||
            M[c][f] + M[a][d] + M[e][b] < r ||
            M[c][f] + M[a][e] + M[d][b] < r ||
            M[c][a] + M[f][d] + M[e][b] < r) {
            *yesno = 1;
        }
    }

    for (k = 4; k <= 5 && *yesno == 0; k++) {
        for (i = 0; i < k; i++) Set[i] = i;
        for (; Set[k-1] < 5 && *yesno == 0; next_set (k, Set)) {
            for (i = 0; i < k; i++) {
                del[2*i]   = h[2*Set[i]];
                del[2*i+1] = h[2*Set[i]+1];
            }
            if (k == 4) {
                CCelim_compare_four_swap (del[0], del[1], del[2], del[3],
                 del[4], del[5], del[6], del[7], M, yesno, I);
            } else {
                CCelim_compare_five_swap (del[0], del[1], del[2], del[3],
                 del[4], del[5], del[6], del[7], del[8], del[9], M, yesno, I);
            }
        }
    }
}

static void next_set (int sz, int *Set)
{
   int i;
   for (i=0; i < sz-1 && Set[i]+1 == Set[i+1]; i++) Set[i] = i;
   Set[i] = Set[i]+1;
}
