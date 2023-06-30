/*
MIT License

Copyright (c) 1995--2023 David Applegate, Robert Bixby, Vasek Chvatal, 
and William Cook                                         

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
/*  Routines from the Concorde library that are used in edge elimination.   */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/*  This file is part of CONCORDE                                           */
/*                                                                          */
/*  Use at your own risk.  We make no guarantees about the                  */
/*  correctness or usefulness of this code.                                 */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/*                   MEMORY ALLOCATION MACROS                               */
/*                                                                          */
/*                           TSP CODE                                       */
/*                                                                          */
/*                                                                          */
/*  Written by:  Applegate, Bixby, Chvatal, and Cook                        */
/*  Date: February 2, 1995 (cofeb16)                                        */
/*                                                                          */
/*    EXPORTED FUNCTIONS:                                                   */
/*                                                                          */
/*  void *CCutil_allocrus (size_t size)                                     */
/*    RETURNS a pointer to an allocated block of "size" memory.             */
/*                                                                          */
/*  void CCutil_freerus (void *ptr)                                         */
/*    FREES ptr.                                                            */
/*                                                                          */
/*  void *CCutil_reallocrus (void *ptr, size_t size)                        */
/*    REALLOCS ptr to size bytes.                                           */
/*                                                                          */
/*  int CCutil_reallocrus_scale (void **pptr, int *pnnum, int count,        */
/*      double scale, size_t size)                                          */
/*    void **pptr (a reference to the pointer to the allocated space)       */
/*    int *pnnum (a reference to the number of objects in the               */
/*                allocated space)                                          */
/*    int count (a minimum value for the new nnum)                          */
/*    double scale (a scale factor to apply to nnum)                        */
/*    int size (the size of objects to be realloced)                        */
/*    RETURNS 0 if *pptr was successfully changed to point to at            */
/*            least max(*pnnum*scale, *pnnum+1000, count) objects.          */
/*            *pnnum is changed to the new object count.                    */
/*            Otherwise, prints an error message, leaves *pptr and          */
/*            *pnnum alone, and returns nonzero.                            */
/*                                                                          */
/*  int CCutil_reallocrus_count (void **pptr, int count,                    */
/*      size_t size)                                                        */
/*    void **pptr (a reference to the pointer to the allocated space)       */
/*    int count (number of objects to be realloced)                         */
/*    int size (the size of the objects to be realloced)                    */
/*    RETURNS 0 is successful, and 1 if the realloc failed.                 */
/*                                                                          */
/*  CCbigchunkptr *CCutil_bigchunkalloc (void)                              */
/*         RETURNS a CCbigchunkptr with the "this_one" field loaded with a  */
/*                 a pointer to a bigchunk of memory.                       */
/*    NOTES:                                                                */
/*       The idea is to use bigchunks (the size of a bigchunk is defined    */
/*       by CC_BIGCHUNK in util.h) to supply local routines with memory     */
/*       for ptrs, so the memory can be shared with other                   */
/*       local routines.                                                    */
/*                                                                          */
/*  CCutil_bigchunkfree (CCbigchunkptr *bp)                                 */
/*    ACTION: Frees a CCbigchunkptr.                                        */
/*                                                                          */
/*  void CCptrworld_init (CCptrworld *world)                                */
/*     initialize a CCptrworld with 1 reference                             */
/*                                                                          */
/*  void CCptrworld_add (CCptrworld *world)                                 */
/*     add a reference to a CCptrworld                                      */
/*                                                                          */
/*  void CCptrworld_delete (CCptrworld *world)                              */
/*     delete a reference to a ptrworld, and free if no more references     */
/*                                                                          */
/****************************************************************************/

#include "ccutil.h"

typedef struct CCbigchunk {
    char space[CC_BIGCHUNK];
    CCbigchunkptr ptr;
} CCbigchunk;
    
void *CCutil_allocrus (size_t size)
{
    void *mem = (void *) NULL;

    if (size == 0) {
        fprintf (stderr, "Warning: 0 bytes allocated\n");
    }

    mem = (void *) malloc (size);
    if (mem == (void *) NULL) {
        fprintf (stderr, "Out of memory. Asked for %d bytes\n", (int) size);
    }
    return mem;
}

void CCutil_freerus (void *p)
{
    if (!p) {
        fprintf (stderr, "Warning: null pointer freed\n");
        return;
    }

    free (p);
}

void *CCutil_reallocrus (void *ptr, size_t size)
{
    void *newptr;

    if (!ptr) {
        return CCutil_allocrus (size);
    } else {
        newptr = (void *) realloc (ptr, size);
        if (!newptr) {
            fprintf (stderr, "Out of memory.  Tried to grow to %d bytes\n",
                     (int) size);
        }
        return newptr;
    }
}

int CCutil_reallocrus_scale (void **pptr, int *pnnum, int count, double scale,
        size_t size)
{
    int newsize = (int) (((double) *pnnum) * scale);
    void *p;

    if (newsize < *pnnum+1000) newsize = *pnnum+1000;
    if (newsize < count) newsize = count;
    p = CCutil_reallocrus (*pptr, newsize * size);
    if (!p) {
        return 1;
    } else {
        *pptr = p;
        *pnnum = newsize;
        return 0;
    }
}

int CCutil_reallocrus_count (void **pptr, int count, size_t size)
{
    void *p = CCutil_reallocrus (*pptr, count * size);

    if (!p) {
        return 1;
    } else {
        *pptr = p;
        return 0;
    }
}


CCbigchunkptr *CCutil_bigchunkalloc (void)
{
    CCbigchunk *p = CC_SAFE_MALLOC (1, CCbigchunk);

    if (p == (CCbigchunk *) NULL) {
        fprintf (stderr, "Out of memory in CCutil_bigchunkalloc\n");
        return (CCbigchunkptr *) NULL;
    }
    p->ptr.this_chunk = p;
    p->ptr.this_one = (void *) p->space;
    return &(p->ptr);
}

void CCutil_bigchunkfree (CCbigchunkptr *bp)
{
    /* This copy is necessary since CC_FREE zeros its first argument */
    CCbigchunk *p = bp->this_chunk;
    
    CC_FREE (p, CCbigchunk);
}

void CCptrworld_init (CCptrworld *world)
{
    world->refcount = 1;
    world->freelist = (void *) NULL;
    world->chunklist = (CCbigchunkptr *) NULL;
}

void CCptrworld_add (CCptrworld *world)
{
    world->refcount++;
}

void CCptrworld_delete (CCptrworld *world)
{
    world->refcount--;
    if (world->refcount <= 0) {
        CCbigchunkptr *bp, *bpnext;

        for (bp = world->chunklist ; bp; bp = bpnext) {
            bpnext = bp->next;
            CCutil_bigchunkfree (bp);
        }
        world->chunklist = (CCbigchunkptr *) NULL;
        world->freelist = (void *) NULL;
        world->refcount = 0;
    }
}

/****************************************************************************/
/*                                                                          */
/*                        PORTABLE GETOPT                                   */
/*                                                                          */
/*                           TSP CODE                                       */
/*                                                                          */
/*                                                                          */
/*  Written by:  Applegate, Bixby, Chvatal, and Cook                        */
/*  Date: 1993 (?) (fmfeb02)                                                */
/*  Modified: 15 February 1995 (bico)  - added warning                      */
/*                                                                          */
/*    EXPORTED FUNCTIONS:                                                   */
/*                                                                          */
/*  int CCutil_bix_getopt (int argc, char **argv, const char *def,          */
/*      int *p_optind, char **p_optarg)                                     */
/*     parse an argument list                                               */
/*                                                                          */
/****************************************************************************/

#include "ccutil.h"

int CCutil_bix_getopt (int ac, char **av, const char *def, int *p_optind,
        char **p_optarg)
{
    int c;
    char *sp = av[*p_optind];
    char bwarn[2];

    if (*p_optind < 1 || *p_optind >= ac) {
        *p_optind = ac;
        return (EOF);
    }
    if ((int) *sp != (int) '-')
        return (EOF);
    if ((int) *(sp + 1) == (int) '-') {
        (*p_optind)++;
        return (EOF);
    }
    (av[*p_optind])++;
    sp++;
    while ((int) *sp != (int) *def && (int) *def != (int) '\0')
            def++;
    if ((int) *def == (int) '\0') {
        *p_optind = ac;
        bwarn[0] = *sp;                          /* Bico: February 8, 1995 */
        bwarn[1] = '\0';
        printf ("Illegal option: -%s\n", bwarn);
        return CC_BIX_GETOPT_UNKNOWN;
    }
    if ((int) *(def + 1) != (int) ':') {
        c = *sp;
        if ((int) *(sp + 1) != (int) '\0')
            *sp = '-';
        else
            (*p_optind)++;
        return (c);
    } else {
        if ((int) *(sp + 1) != (int) '\0') {
            *p_optarg = sp + 1;
            c = *sp;
            (*p_optind)++;
            return (c);
        } else if (*p_optind >= ac - 1) {
            *p_optind = ac;
            return (EOF);
        } else {
            *p_optarg = av[*p_optind + 1];
            c = *sp;
            *p_optind += 2;
            return (c);
        }
    }
}

/****************************************************************************/
/*                                                                          */
/*                           DHEAP ROUTINES                                 */
/*                                                                          */
/*                                                                          */
/*                              TSP CODE                                    */
/*  Written by:  Applegate, Bixby, Chvatal, and Cook                        */
/*  Date: February 9, 1995                                                  */
/*  Reference: R.E. Tarjan, Data Structures and Network Algorithms          */
/*                                                                          */
/*    EXPORTED FUNCTIONS:                                                   */
/*                                                                          */
/*  int CCutil_dheap_init (CCdheap *h, int k)                               */
/*        -h should point to a CCdheap struct.                              */
/*        -k the max number of elements in the dheap.                       */
/*                                                                          */
/*  void CCutil_dheap_free (CCdheap *h)                                     */
/*    -frees the spaces allocated by CCutil_dheap_init                      */
/*                                                                          */
/*  int CCutil_dheap_resize (CCdheap *h, int newsize)                       */
/*    -REALLOCs h so it can contain newsize elements.                       */
/*    -returns -1 if it can't resize the heap.                              */
/*                                                                          */
/*  int CCutil_dheap_findmin (CCdheap *h)                                   */
/*    -returns the index of the element with min value h->key[i]            */
/*    -returns -1 if no elements in heap.                                   */
/*                                                                          */
/*  int CCutil_dheap_insert (CCdheap *h, int i)                             */
/*    -inserts the element with index i (so its key should be loaded        */
/*     beforehand in h->key[i]).                                            */
/*                                                                          */
/*  void CCutil_dheap_delete (CCdheap *h, int i)                            */
/*    -deletes the element with index i.                                    */
/*                                                                          */
/*  int CCutil_dheap_deletemin (CCdheap *h)                                 */
/*    -returns the min element in the heap, and deletes the min element     */
/*    -returns -1 if no elements in heap.                                   */
/*                                                                          */
/*  void CCutil_dheap_changekey (CCdheap *h, int i, double newkey)          */
/*    -changes the key of the element with index i to newkey.               */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/*  NOTES:                                                                  */
/*      A k-element heap will malloc 16k bytes of memory. If memory is      */
/*  tight, using integer keys (instead of doubles), brings it down to       */
/*  12k bytes, and if arbitrary deletions are not required, with a little   */
/*  rewriting, the h->loc field can be eliminated, bring the space down     */
/*  to 8k bytes.                                                            */
/*      These routines work with indices into the h->key array, so in       */
/*  some cases, you will need to maintain a separate names array to know    */
/*  what element belongs to index i. For an example, see the k_nearest      */
/*  code in kdnear.c.                                                       */
/*                                                                          */
/****************************************************************************/

#include "ccutil.h"

#define HEAP_D 3
#define HEAP_UP(x) (((x)-1)/HEAP_D)
#define HEAP_DOWN(x) (((x)*HEAP_D)+1)


static void
    dheap_siftup (CCdheap *h, int i, int x),
    dheap_siftdown (CCdheap *h, int i, int x);

static int
    dheap_minchild (int x, CCdheap *h);



int CCutil_dheap_init (CCdheap *h, int k)
{
    h->loc = (int *) NULL;
    h->key = (double *) NULL;
    h->entry = CC_SAFE_MALLOC (k, int);
    if (!h->entry)
        return 1;
    h->loc = CC_SAFE_MALLOC (k, int);
    if (!h->loc) {
        CC_FREE (h->entry, int);
        return 1;
    }
    h->key = CC_SAFE_MALLOC (k, double);
    if (!h->key) {
        CC_FREE (h->entry, int);
        CC_FREE (h->loc, int);
        return 1;
    }
    h->total_space = k;
    h->size = 0;
    return 0;
}

void CCutil_dheap_free (CCdheap *h)
{
    CC_IFFREE (h->entry, int);
    CC_IFFREE (h->loc, int);
    CC_IFFREE (h->key, double);
}

int CCutil_dheap_resize (CCdheap *h, int newsize)
{
    void *tmp_ptr;
    if (newsize < h->size || newsize < h->total_space) return 0;

    tmp_ptr = (void *) h->key;
    if (CCutil_reallocrus_count (&tmp_ptr, newsize, sizeof (double))) {
        return -1;
    }
    h->key = (double *) tmp_ptr;

    tmp_ptr = (void *) h->entry;
    if (CCutil_reallocrus_count (&tmp_ptr,  newsize, sizeof (int))) {
        return -1;
    }
    h->entry = (int *) tmp_ptr;

    tmp_ptr = (void *) h->loc;
    if (CCutil_reallocrus_count (&tmp_ptr, newsize, sizeof (int))) {
        return -1;
    }
    h->loc = (int *) tmp_ptr;

    h->total_space = newsize;

    return 0;
}

int CCutil_dheap_findmin (CCdheap *h)
{
    if (h->size == 0)
        return -1;
    else
        return h->entry[0];
}

int CCutil_dheap_insert (CCdheap *h, int i)
{
    if (h->size >= h->total_space) {
        fprintf (stderr, "Error - heap already full\n");
        return 1;
    }
    h->size++;
    dheap_siftup (h, i, h->size - 1);
    return 0;
}

void CCutil_dheap_delete (CCdheap *h, int i)
{
    int j;

    h->size--;
    j = h->entry[h->size];
    h->entry[h->size] = -1;

    if (j != i) {
        if (h->key[j] <= h->key[i]) {
            dheap_siftup (h, j, h->loc[i]);
        } else {
            dheap_siftdown (h, j, h->loc[i]);
        }
    }
}

int  CCutil_dheap_deletemin (CCdheap *h)
{
    int i;

    if (h->size == 0)
        return -1;
    else {
        i = h->entry[0];
        CCutil_dheap_delete (h, i);
        return i;
    }
}

void CCutil_dheap_changekey (CCdheap *h, int i, double newkey)
{
    if (newkey < h->key[i]) {
        h->key[i] = newkey;
        dheap_siftup (h, i, h->loc[i]);
    } else if (newkey > h->key[i]) {
        h->key[i] = newkey;
        dheap_siftdown (h, i, h->loc[i]);
    }
}

static void dheap_siftup (CCdheap *h, int i, int x)
{
    int p;

    p = HEAP_UP (x);
    while (x && h->key[h->entry[p]] > h->key[i]) {
        h->entry[x] = h->entry[p];
        h->loc[h->entry[p]] = x;
        x = p;
        p = HEAP_UP (p);
    }
    h->entry[x] = i;
    h->loc[i] = x;
}

static void dheap_siftdown (CCdheap *h, int i, int x)
{
    int c;

    c = dheap_minchild (x, h);

    while (c >= 0 && h->key[h->entry[c]] < h->key[i]) {
        h->entry[x] = h->entry[c];
        h->loc[h->entry[c]] = x;
        x = c;
        c = dheap_minchild (c, h);
    }
    h->entry[x] = i;
    h->loc[i] = x;
}

static int dheap_minchild (int x, CCdheap *h)
{
    int c = HEAP_DOWN (x);
    int cend;
    double minval;
    int minloc;

    if (c >= h->size)
        return -1;
    minval = h->key[h->entry[c]];
    minloc = c;
    cend = c + HEAP_D;
    if (h->size < cend)
        cend = h->size;
    for (c++; c < cend; c++) {
        if (h->key[h->entry[c]] < minval) {
            minval = h->key[h->entry[c]];
            minloc = c;
        }
    }
    return minloc;
}

/****************************************************************************/
/*                                                                          */
/*            ROUTINES TO MAP NODE PAIRS TO EDGES                           */
/*                                                                          */
/*                           TSP CODE                                       */
/*                                                                          */
/*                                                                          */
/*  Written by:  Applegate, Bixby, Chvatal, and Cook                        */
/*  Date: September 27, 1995                                                */
/*                                                                          */
/*    EXPORTED FUNCTIONS:                                                   */
/*                                                                          */
/*  int CCutil_edgehash_init (CCutil_edgehash *h, int size)                 */
/*    MISSING                                                               */
/*                                                                          */
/*  int CCutil_edgehash_add (CCutil_edgehash *h, int end1, int end2,        */
/*      int val)                                                            */
/*    MISSING                                                               */
/*                                                                          */
/*  int CCutil_edgehash_set (CCutil_edgehash *h, int end1, int end2,        */
/*      int val)                                                            */
/*    MISSING                                                               */
/*                                                                          */
/*  int CCutil_edgehash_del (CCutil_edgehash *h, int end1, int end2)        */
/*    MISSING                                                               */
/*                                                                          */
/*  int CCutil_edgehash_find (CCutil_edgehash *h, int end1, int end2,       */
/*      int *val)                                                           */
/*    MISSING                                                               */
/*                                                                          */
/*  int CCutil_edgehash_getall (CCutil_edgehash *h, int *ecount,            */
/*      int **elist, int **elen);                                           */
/*    MISSING                                                               */
/*                                                                          */
/*  void CCutil_edgehash_delall (CCutil_edgehash *h)                        */
/*    MISSING                                                               */
/*                                                                          */
/*  void CCutil_edgehash_free (CCutil_edgehash *h)                          */
/*    MISSING                                                               */
/*                                                                          */
/****************************************************************************/

#include "ccutil.h"

CC_PTRWORLD_ROUTINES (CCutil_edgeinf, edgeinf_alloc, edgeinf_bulkalloc,
        edgeinf_free)
CC_PTRWORLD_LISTFREE_ROUTINE (CCutil_edgeinf, edgeinf_listfree, edgeinf_free)

int CCutil_edgehash_init (CCutil_edgehash *h, int size)
{
    unsigned int i;

    h->size = CCutil_nextprime ((unsigned int) size);
    h->mult = (int) sqrt ((double) h->size);
    h->table = CC_SAFE_MALLOC ((int) h->size, CCutil_edgeinf *);
    CCptrworld_init (&h->edgeinf_world);
    
    if (!h->table) {
        h->size = 0;
        return 1;
    }
    for (i=0; i<h->size; i++) {
        h->table[i] = (CCutil_edgeinf *) NULL;
    }
    return 0;
}

int CCutil_edgehash_add (CCutil_edgehash *h, int end1, int end2, int val)
{
    int t;
    unsigned int loc;
    CCutil_edgeinf *e;

    if (h->size == 0) return 1;
    e = edgeinf_alloc(&h->edgeinf_world);
    if (!e) return 1;

    if (end1 > end2) CC_SWAP (end1, end2, t);
    loc = (end1 * h->mult + end2) % h->size;
    e->ends[0] = end1;
    e->ends[1] = end2;
    e->val = val;
    e->next = h->table[loc];
    h->table[loc] = e;
    return 0;
}

int CCutil_edgehash_set (CCutil_edgehash *h, int end1, int end2, int val)
{
    int t;
    unsigned int loc;
    CCutil_edgeinf *e;

    if (h->size == 0) return 1;
    if (end1 > end2) CC_SWAP (end1, end2, t);
    loc = (end1 * h->mult + end2) % h->size;
    for (e = h->table[loc]; e; e = e->next) {
        if (e->ends[0] == end1 && e->ends[1] == end2) {
            e->val = val;
            return 0;
        }
    }
    e = edgeinf_alloc(&h->edgeinf_world);
    if (!e) return 1;

    e->ends[0] = end1;
    e->ends[1] = end2;
    e->val = val;
    e->next = h->table[loc];
    h->table[loc] = e;
    return 0;
}

int CCutil_edgehash_del (CCutil_edgehash *h, int end1, int end2)
{
    int t;
    CCutil_edgeinf **prev;
    CCutil_edgeinf *p;

    if (end1 > end2) CC_SWAP (end1, end2, t);
    if (h->size == 0) return 1;

    prev  = &h->table[(end1 * h->mult + end2) % h->size];
    p = *prev;
    while (p) {
        if (p->ends[0] == end1 && p->ends[1] == end2) {
            *prev = p->next;
            edgeinf_free (&h->edgeinf_world, p);
            return 0;
        }
        prev = &p->next;
        p = *prev;
    }
    return 1;
}

int CCutil_edgehash_getall (CCutil_edgehash *h, int *ecount, int **elist,
        int **elen)
{
    unsigned int i;
    int k = 0;
    CCutil_edgeinf *p;

    for (i = 0; i < h->size; i++) {
        for (p = h->table[i]; p; p = p->next) k++;
    } 
   
    if (k > 0) {
        *elist = CC_SAFE_MALLOC (2*k, int);
        *elen  = CC_SAFE_MALLOC (k, int);
        if (!(*elist) || !(*elen)) {
            fprintf (stderr, "out of memory in CCutil_edgehash_getall\n");
            CC_IFFREE (*elist, int);
            CC_IFFREE (*elen, int);
            return 1;
        }
        *ecount = k;
        k = 0;
        for (i = 0; i < h->size; i++) {
            for (p = h->table[i]; p; p = p->next) {
                (*elist)[2*k]   = p->ends[0];
                (*elist)[2*k+1] = p->ends[1];
                (*elen)[k++]    = p->val;
            }
        } 
    } else {
        *elist = (int *) NULL;
        *elen = (int *) NULL;
        *ecount = 0;
    }
    
    return 0;
}

void CCutil_edgehash_delall (CCutil_edgehash *h)
{
    unsigned int i;

    for (i=0; i<h->size; i++) {
        if (h->table[i]) {
            edgeinf_listfree (&h->edgeinf_world, h->table[i]);
            h->table[i] = (CCutil_edgeinf *) NULL;
        }
    }
}

int CCutil_edgehash_find (CCutil_edgehash *h, int end1, int end2, int *val)
{
    int t;
    CCutil_edgeinf *p;

    *val = 0;
    if (h->size == 0) return -1;
    if (end1 > end2) CC_SWAP (end1, end2, t);

    for (p = h->table[(end1 * h->mult + end2) % h->size]; p; p = p->next) {
        if (p->ends[0] == end1 && p->ends[1] == end2) {
            *val = p->val;
            return 0;
        }
    }
    return -1;
}

void CCutil_edgehash_free (CCutil_edgehash *h)
{
    CCutil_edgehash_delall (h);
    CC_FREE (h->table, CCutil_edgeinf *);
    CCptrworld_delete (&h->edgeinf_world);
    h->size = 0;
}

/****************************************************************************/
/*                                                                          */
/*                    EDGE LIST UTILITY ROUTINES                            */
/*                                                                          */
/*                              TSP CODE                                    */
/*                                                                          */
/*  Written by: Applegate, Bixby, Chvatal, and Cook                         */
/*  Date: February 8, 1995                                                  */
/*                                                                          */
/*    EXPORTED FUNCTIONS:                                                   */
/*                                                                          */
/*  int CCutil_edge_to_cycle (int ncount, int *elist, int *yesno,           */
/*      int *cyc)                                                           */
/*    CONVERTS an edgelist to a cycle.                                      */
/*     -ncount is the number of nodes.                                      */
/*     -elist is an edgelist in end1 end2 format.                           */
/*     -yesno returns 1 if elist describes a tour and 0 otherwise.          */
/*     -cyc returns the cycle in permutation format if it is not NULL       */
/*      (if cyc is not NULL, then it should point to an array of            */
/*      length at least ncount).                                            */
/*     Returns a nonzero value if there was an error.                       */
/*                                                                          */
/****************************************************************************/

#include "ccutil.h"

int CCutil_edge_to_cycle (int ncount, int *elist, int *yesno, int *cyc)
{
    int *Lside, *Rside;
    int i, k, end1, end2, prev, this, next, start, okfirst, first = 0;
    int rval = 0;

    *yesno = 0;

    Lside = CC_SAFE_MALLOC (ncount, int);
    Rside = CC_SAFE_MALLOC (ncount, int);
    if (!Lside || !Rside) {
        fprintf (stderr, "out of memory in CCutil_edge_to_cycle\n");
        rval = 1; goto CLEANUP;
    }

    for (i = 0; i < ncount; i++) {
        Lside[i] = Rside[i] = -1;
    }

    for (i = 0, k = 0; i < ncount; i++) {
        end1 = elist[k++];
        end2 = elist[k++];
        if (Lside[end1] == -1)
            Lside[end1] = end2;
        else
            Rside[end1] = end2;
        if (Lside[end2] == -1)
            Lside[end2] = end1;
        else
            Rside[end2] = end1;
    }

    for (i = 0, k = 0; i < ncount; i++) {
        end1 = elist[k++];
        end2 = elist[k++];
        if (Lside[end1] == -1 || Rside[end1] == -1 ||
            Lside[end2] == -1 || Rside[end2] == -1) {
            *yesno = 0;  goto CLEANUP;
        }
    }
    start = elist[0];
    prev = -2;
    this = start;
    k = 0;
    okfirst = 0;
    do {
        if (this == first)
           okfirst = 1;
        if (Lside[this] != prev)
            next = Lside[this];
        else
            next = Rside[this];
        prev = this;
        this = next;
        k++;
    } while (next != start && k < ncount);

    if (k != ncount || !okfirst) {
        *yesno = 0;  goto CLEANUP;
    }

    *yesno = 1;

    if (cyc) {
        start = first;
        prev = -2;
        this = start;
        k = 0;
        do {
            cyc[k++] = this;
            if (Lside[this] != prev)
                next = Lside[this];
            else
                next = Rside[this];
            prev = this;
            this = next;
        } while (next != start && k < ncount);
    }


CLEANUP:

    CC_IFFREE (Lside, int);
    CC_IFFREE (Rside, int);

    return rval;
}

/****************************************************************************/
/*                                                                          */
/*                   CODE TO READ ASCI INTS FAST                            */
/*                                                                          */
/*                           TSP CODE                                       */
/*                                                                          */
/*                                                                          */
/*  Written by:  Dave                                                       */
/*  Date: Septmember 1994 (Bonn)  (cofeb16)                                 */
/*                                                                          */
/*                                                                          */
/*    EXPORTED FUNCTIONS:                                                   */
/*                                                                          */
/*  int CCutil_readint (FILE *f)                                            */
/*      - Returns the next int in the file f.                               */
/*    This is much faster that scanf. It is useful for big files and        */
/*    and for profiling.                                                    */
/*                                                                          */
/****************************************************************************/

#include "ccutil.h"

int CCutil_readint (FILE *f)
{
    int v = 0;
    int c;

    while (( c = getc(f)) != EOF && !((c >= '0' && c <= '9') || c == '-'));
    if (c == '-') {
        v = 0;
        while ((c = getc(f)) != EOF && c >= '0' && c <= '9') {
            v = v * 10 + c - '0';
        }
        return -v;
    } else {
        v = c - '0';
        while ((c = getc(f)) != EOF && c >= '0' && c <= '9') {
            v = v * 10 + c - '0';
        }
        return v;
    }
}


/****************************************************************************/
/*                                                                          */
/*  This file is part of CONCORDE                                           */
/*                                                                          */
/*  (c) Copyright 1995--1999 by David Applegate, Robert Bixby,              */
/*  Vasek Chvatal, and William Cook                                         */
/*                                                                          */
/*  Permission is granted for academic research use.  For other uses,       */
/*  contact the authors for licensing options.                              */
/*                                                                          */
/*  Use at your own risk.  We make no guarantees about the                  */
/*  correctness or usefulness of this code.                                 */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/*                   ROUTINE FOR BUILDING KDTREES                           */
/*                                                                          */
/*  (Based on Jon Bentley's paper "K-d trees for semidynamic point sets")   */
/*                                                                          */
/*                           TSP CODE                                       */
/*                                                                          */
/*                                                                          */
/*  Written by:  Applegate, Bixby, Chvatal, and Cook                        */
/*  Date: February 24, 1995 (cofeb24)                                       */
/*        Modfied for 3D data, June 21, 2018                                */
/*                                                                          */
/*                                                                          */
/*    EXPORTED FUNCTIONS:                                                   */
/*                                                                          */
/*  int CCkdtree_build (CCkdtree *kt, int ncount, CCdatagroup *dat,         */
/*      double *wcoord, CCrandstate *rstate)                                */
/*    -When called, intree should point to a CCkdtree struct that the       */
/*     funtion will load with the tree it builds. The wcoord array          */
/*     is used for node weights (like in Held-Karp), it can be NULL.        */
/*     The node weights must be nonegative (for cutoffs).                   */
/*                                                                          */
/*  void CCkdtree_free (CCkdtree *kt)                                       */
/*    -Frees the space (including the ptrs) used by kt.                     */
/*                                                                          */
/*  void CCkdtree_delete (CCkdtree *kt, int k)                              */
/*    -Deletes the point k from the CCkdtree kt.                            */
/*                                                                          */
/*  void CCkdtree_undelete (CCkdtree *kt, int k)                            */
/*    -Puts the previously deleted point k back into kt.                    */
/*                                                                          */
/*  void CCkdtree_delete_all (CCkdtree *kt, int ncount)                     */
/*    -Deletes all points in kt.                                            */
/*                                                                          */
/*  void CCkdtree_undelete_all (CCkdtree *kt, int ncount)                   */
/*         -Puts all deleted points back in kt. Used to cleanup trees.      */
/*                                                                          */
/*    NOTES:                                                                */
/*       On a 32 bit machine, a CCkdtree on n nodes needs about 52n         */
/*     bytes of memory. CCkdtree_build will return 1 if an error            */
/*     occurs (most likely running out of memory).                          */
/*       CCutil_sprand () should be called before calling                   */
/*     CCkdtree_build ().                                                   */
/*                                                                          */
/****************************************************************************/

#include "ccutil.h"

#define CUTOFF 5
#define BNDS_DEPTH 5   /* When bnds info is recorded */
#define BIGDOUBLE (1e30)


static void
    kdtree_free_work (CCkdnode *p, CCptrworld *kdnode_world,
        CCptrworld *kdbnds_world),
    kdtree_free_world (CCptrworld *kdnode_world, CCptrworld *kdbnds_world);
static unsigned char
    findmaxspread (int l, int u, CCkdtree *thetree, double *datx,
           double *daty, double *datz, double *datw);
static CCkdnode
    *build (int l, int u, int *depth, double *current_bnds_x,
           double *current_bnds_y, double *current_bnds_z, CCkdtree *thetree,
           double *datx, double *daty, double *datz, double *datw,
           CCrandstate *rstate);


CC_PTRWORLD_ROUTINES (CCkdnode, kdnodealloc, kdnode_bulk_alloc, kdnodefree)
CC_PTRWORLD_LEAKS_ROUTINE (CCkdnode, kdnode_check_leaks, empty, char)

CC_PTRWORLD_ROUTINES (CCkdbnds, kdbndsalloc, kdbnds_bulk_alloc, kdbndsfree)
CC_PTRWORLD_LEAKS_ROUTINE (CCkdbnds, kdbnds_check_leaks, x[0], double)

int CCkdtree_build (CCkdtree *intree, int ncount, CCdatagroup *dat,
        double *wcoord, CCrandstate *rstate)
{
    int rval = 0, i, depth;
    double current_bnds_x[2], current_bnds_y[2], current_bnds_z[2];
    CCkdtree *thetree = (CCkdtree *) NULL;

    CCptrworld_init (&intree->kdnode_world);
    CCptrworld_init (&intree->kdbnds_world);
    
    if (wcoord != (double *) NULL) {
        for (i = 0; i < ncount; i++) {
            if (wcoord[i] < -0.00000001) {
                fprintf (stderr, "Cannot build with negative node weights\n");
                rval = 1; goto CLEANUP;
            }
        }
    }

    thetree = intree;
    thetree->perm = (int *) NULL;
    thetree->bucketptr = (CCkdnode **) NULL;

    CC_MALLOC (thetree->perm, ncount, int);
    for (i = 0; i < ncount; i++) thetree->perm[i] = i;
    CC_MALLOC (thetree->bucketptr, ncount, CCkdnode *);

    depth = 0;
    current_bnds_x[0] = -BIGDOUBLE;
    current_bnds_x[1] =  BIGDOUBLE;
    current_bnds_y[0] = -BIGDOUBLE;
    current_bnds_y[1] =  BIGDOUBLE;
    current_bnds_z[0] = -BIGDOUBLE;
    current_bnds_z[1] =  BIGDOUBLE;

    thetree->root = build (0, ncount - 1, &depth, current_bnds_x,
         current_bnds_y, current_bnds_z, thetree, dat->x, dat->y,
         dat->z, wcoord, rstate);
    if (!thetree->root) {
        rval = 1; goto CLEANUP;
    } else {
        thetree->root->father = (CCkdnode *) NULL;
    }

CLEANUP:
    if (rval && thetree) {
        CC_FREE (thetree->perm, int);
        CC_FREE (thetree->bucketptr, CCkdnode *);
    }
    return rval;
}

void CCkdtree_free (CCkdtree *kt)
{
    if (kt->perm)
        CC_FREE (kt->perm, int);
    if (kt->bucketptr)
        CC_FREE (kt->bucketptr, CCkdnode *);
    kdtree_free_work (kt->root, &kt->kdnode_world, &kt->kdbnds_world);
    kt->root = (CCkdnode *) NULL;
    
    kdtree_free_world (&kt->kdnode_world, &kt->kdbnds_world);
}

static void kdtree_free_world (CCptrworld *kdnode_world,
        CCptrworld *kdbnds_world)
{
    int total, onlist;

    if (kdnode_check_leaks (kdnode_world, &total, &onlist)) {
        fprintf (stderr, "WARNING: %d outstanding kdnodes\n",
                 total - onlist);
    }
    if (kdbnds_check_leaks (kdbnds_world, &total, &onlist)) {
        fprintf (stderr, "WARNING: %d outstanding kdbnds\n", total - onlist);
    }
    CCptrworld_delete (kdnode_world);
    CCptrworld_delete (kdbnds_world);
}

static void kdtree_free_work (CCkdnode *p, CCptrworld *kdnode_world,
        CCptrworld *kdbnds_world)
{
    if (p->bucket) {
        if (p->bnds)
            kdbndsfree (kdbnds_world, p->bnds);
        kdnodefree (kdnode_world, p);
    } else {
        kdtree_free_work (p->loson, kdnode_world, kdbnds_world);
        kdtree_free_work (p->hison, kdnode_world, kdbnds_world);
        if (p->bnds)
            kdbndsfree (kdbnds_world, p->bnds);
        kdnodefree (kdnode_world, p);
    }
}

static CCkdnode *build (int l, int u, int *depth, double *current_bnds_x,
        double *current_bnds_y, double *current_bnds_z, CCkdtree *thetree,
        double *datx, double *daty, double *datz, double *datw,
        CCrandstate *rstate)
{
    CCkdnode *p;
    int i, m;
    double savebnd;

    (*depth)++;
    p = kdnodealloc (&thetree->kdnode_world);
    if (!p) {
        (*depth)--;
        return (CCkdnode *) NULL;
    }
    p->empty = 0;

    if (u - l + 1 < CUTOFF) {
        p->bucket = 1;
        p->lopt = l;
        p->hipt = u;
        for (i = l; i <= u; i++)
            thetree->bucketptr[thetree->perm[i]] = p;
        p->bnds = (CCkdbnds *) NULL;
    } else {
        p->bucket = 0;
        if (!((*depth) % BNDS_DEPTH)) {
            p->bnds = kdbndsalloc (&thetree->kdbnds_world);
            if (!p->bnds) {
                (*depth)--;
                kdnodefree (&thetree->kdbnds_world, p);
                return (CCkdnode *) NULL;
            }
            p->bnds->x[0] = current_bnds_x[0];
            p->bnds->x[1] = current_bnds_x[1];
            p->bnds->y[0] = current_bnds_y[0];
            p->bnds->y[1] = current_bnds_y[1];
        } else {
            p->bnds = (CCkdbnds *) NULL;
        }

        p->cutdim = findmaxspread (l, u, thetree, datx, daty, datz, datw);
        m = (l + u) / 2;
        switch (p->cutdim) {
        case 0:
            CCutil_rselect (thetree->perm, l, u, m, datx, rstate);
            p->cutval = datx[thetree->perm[m]];

            savebnd = current_bnds_x[1];
            current_bnds_x[1] = p->cutval;
            p->loson = build (l, m, depth, current_bnds_x, current_bnds_y,
                      current_bnds_z, thetree, datx, daty, datz, datw, rstate);
            if (!p->loson) {
                (*depth)--;
                kdnodefree (&thetree->kdnode_world, p);
                return (CCkdnode *) NULL;
            }
            current_bnds_x[1] = savebnd;

            savebnd = current_bnds_x[0];
            current_bnds_x[0] = p->cutval;
            p->hison = build (m + 1, u, depth, current_bnds_x, current_bnds_y,
                      current_bnds_z, thetree, datx, daty, datz, datw, rstate);
            if (!p->hison) {
                (*depth)--;
                kdnodefree (&thetree->kdnode_world, p);
                return (CCkdnode *) NULL;
            }
            current_bnds_x[0] = savebnd;

            break;
        case 1:
            CCutil_rselect (thetree->perm, l, u, m, daty, rstate);
            p->cutval = daty[thetree->perm[m]];

            savebnd = current_bnds_y[1];
            current_bnds_y[1] = p->cutval;
            p->loson = build (l, m, depth, current_bnds_x, current_bnds_y,
                      current_bnds_z, thetree, datx, daty, datz, datw, rstate);
            if (!p->loson) {
                (*depth)--;
                kdnodefree (&thetree->kdnode_world, p);
                return (CCkdnode *) NULL;
            }
            current_bnds_y[1] = savebnd;

            savebnd = current_bnds_y[0];
            current_bnds_y[0] = p->cutval;
            p->hison = build (m + 1, u, depth, current_bnds_x, current_bnds_y,
                      current_bnds_z, thetree, datx, daty, datz, datw, rstate);
            if (!p->hison) {
                (*depth)--;
                kdnodefree (&thetree->kdnode_world, p);
                return (CCkdnode *) NULL;
            }
            current_bnds_y[0] = savebnd;

            break;
        case 2:  
            CCutil_rselect (thetree->perm, l, u, m, datz, rstate);
            p->cutval = datz[thetree->perm[m]];

            savebnd = current_bnds_z[1];
            current_bnds_z[1] = p->cutval;
            p->loson = build (l, m, depth, current_bnds_x, current_bnds_y,
                      current_bnds_z, thetree, datx, daty, datz, datw, rstate);
            if (!p->loson) {
                (*depth)--;
                kdnodefree (&thetree->kdnode_world, p);
                return (CCkdnode *) NULL;
            }
            current_bnds_z[1] = savebnd;

            savebnd = current_bnds_z[0];
            current_bnds_z[0] = p->cutval;
            p->hison = build (m + 1, u, depth, current_bnds_x, current_bnds_y,
                      current_bnds_z, thetree, datx, daty, datz, datw, rstate);
            if (!p->hison) {
                (*depth)--;
                kdnodefree (&thetree->kdnode_world, p);
                return (CCkdnode *) NULL;
            }
            break;
        case 3:
            CCutil_rselect (thetree->perm, l, u, m, datw, rstate);
            p->cutval = datw[thetree->perm[m]];

            p->loson = build (l, m, depth, current_bnds_x, current_bnds_y,
                      current_bnds_z, thetree, datx, daty, datz, datw, rstate);
            if (!p->loson) {
                (*depth)--;
                kdnodefree (&thetree->kdnode_world, p);
                return (CCkdnode *) NULL;
            }

            p->hison = build (m + 1, u, depth, current_bnds_x, current_bnds_y,
                      current_bnds_z, thetree, datx, daty, datz, datw, rstate);
            if (!p->hison) {
                (*depth)--;
                kdnodefree (&thetree->kdnode_world, p);
                return (CCkdnode *) NULL;
            }

            break;
        }
        p->loson->father = p;
        p->hison->father = p;
    }
    (*depth)--;
    return p;
}

static unsigned char findmaxspread (int l, int u, CCkdtree *thetree,
        double *datx, double *daty, double *datz, double *datw)
{
    int i;
    double xmax, xmin, xval, xspread;
    double ymax, ymin, yval, yspread;
    double zmax, zmin, zval, zspread;
    double wmax = 0.0, wmin = 0.0, wval, wspread;

    xmin = datx[thetree->perm[l]];
    xmax = xmin;
    ymin = daty[thetree->perm[l]];
    ymax = ymin;
    if (datz) {
        zmin = datz[thetree->perm[l]];
        zmax = zmin;
    }
    if (datw != (double *) NULL) {
        wmin = datw[thetree->perm[l]];
        wmax = wmin;
    }

    for (i = l + 1; i <= u; i++) {
        xval = datx[thetree->perm[i]];
        if (xval < xmin)      xmin = xval;
        else if (xval > xmax) xmax = xval;

        yval = daty[thetree->perm[i]];
        if (yval < ymin)      ymin = yval;
        else if (yval > ymax) ymax = yval;

        if (datz) {
            zval = daty[thetree->perm[i]];
            if (zval < zmin)      zmin = zval;
            else if (zval > zmax) zmax = zval;
        }

        if (datw) {
            wval = datw[thetree->perm[i]];
            if (wval < wmin)      wmin = wval;
            else if (wval > wmax) wmax = wval;
        }
    }

    xspread = xmax - xmin;
    yspread = ymax - ymin;
    if (datz) zspread = zmax - zmin;
    if (datw) wspread = (wmax - wmin);

    if (datw) {
        if (datz) {
            if (xspread >= yspread && xspread >= zspread &&
                xspread >= wspread) return (unsigned char) 0;
            else if (yspread >= xspread && yspread >= zspread &&
                yspread >= wspread) return (unsigned char) 1;
            else if (zspread >= xspread && zspread >= yspread &&
                zspread >= wspread) return (unsigned char) 2;
            else                    return (unsigned char) 3;

        } else {
            if (xspread >= yspread && xspread >= wspread)
                return (unsigned char) 0;
            else if (yspread >= xspread && yspread >= wspread)
                return (unsigned char) 1;
            else 
                return (unsigned char) 3;
        }
    } else {
        if (datz) {
            if (xspread >= yspread && xspread >= zspread) 
                return (unsigned char) 0;
            else if (yspread >= xspread && yspread >= zspread)
                return (unsigned char) 1;
            else 
                return (unsigned char) 2;
        } else {
            if (xspread >= yspread) return (unsigned char) 0;
            else                    return (unsigned char) 1;
        }
    }
}

void CCkdtree_delete (CCkdtree *kt, int k)
{
    int j, temp;
    CCkdnode *p;

    p = kt->bucketptr[k];
    j = p->lopt;
    while (kt->perm[j] != k)
        j++;
    CC_SWAP (kt->perm[j], kt->perm[p->hipt], temp);
    (p->hipt)--;
    if (p->lopt > p->hipt) {
        p->empty = 1;
        while ((p = p->father) != (CCkdnode *) NULL &&
                p->loson->empty && p->hison->empty)
            p->empty = 1;
    }
}

void CCkdtree_delete_all (CCkdtree *kt, int ncount)
{
    int k;

    for (k = 0; k < ncount; k++) CCkdtree_delete (kt, k);
}


void CCkdtree_undelete (CCkdtree *kt, int k)
{
    int j, temp;
    CCkdnode *p;

    p = kt->bucketptr[k];
    j = p->lopt;
    while (kt->perm[j] != k)
        j++;
    if (j > p->hipt) {
        (p->hipt)++;
        CC_SWAP (kt->perm[j], kt->perm[p->hipt], temp);
        if (p->empty) {
            p->empty = 0;
            while ((p = p->father) != (CCkdnode *) NULL && p->empty)
                p->empty = 0;
        }
    }
}

void CCkdtree_undelete_all (CCkdtree *kt, int ncount)
{
    int k;

    for (k = 0; k < ncount; k++) CCkdtree_undelete (kt, k);
}

/****************************************************************************/
/*                                                                          */
/*  This file is part of CONCORDE                                           */
/*                                                                          */
/*  (c) Copyright 1995--2018 by David Applegate, Robert Bixby,              */
/*  Vasek Chvatal, and William Cook                                         */
/*                                                                          */
/*  Permission is granted for academic research use.  For other uses,       */
/*  contact the authors for licensing options.                              */
/*                                                                          */
/*  Use at your own risk.  We make no guarantees about the                  */
/*  correctness or usefulness of this code.                                 */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/*                 ROUTINES FOR FINDING NEAREST NEIGHBORS                   */
/*                                                                          */
/*  (Based on Jon Bentley's paper "K-d trees for semidynamic point sets")   */
/*                                                                          */
/*                           TSP CODE                                       */
/*                                                                          */
/*  Written by:  Applegate, Bixby, Chvatal, and Cook                        */
/*  Date: February 24, 1995 (cofeb24)                                       */
/*  Changes: August 6, 2018 -  added wcoord to fixed radius search          */
/*           June 21, 2018 - support for 3D data                            */
/*                                                                          */
/*                                                                          */
/*  EXPORTED FUNCTIONS                                                      */
/*                                                                          */
/*  int CCkdtree_k_nearest (CCkdtree *kt, int ncount, int k,                */
/*      CCdatagroup *dat, double *wcoord, int wantlist, int *ocount,        */
/*      int **olist, int silent, CCrandstate *rstate)                       */
/*    RETURNS the k-nearest neighbor graph.                                 */
/*      -kt can be NULL, otherwise it should point to a CCkdtree built      */
/*       by a call to kdbuild ()                                            */
/*      -ncount is the number of points.                                    */
/*      -k is the number of nearest neighbors wanted.                       */
/*      -wcoord is an array of node weights (like Held-Karp), it can        */
/*       be NULL. The weights should be nonnegative.                        */
/*      -wantlist is 1 if you want the function to return the edges.        */
/*      -ocount returns the number of edges (if wantlist is 1) and          */
/*       olist returns the edgelist is end1 end2 format.                    */
/*      -silent will turn off print messages if set to nonzero value.       */
/*                                                                          */
/*  int CCkdtree_quadrant_k_nearest (CCkdtree *kt, int ncount, int k,       */
/*      CCdatagroup *dat, double *wcoord,                                   */
/*      int wantlist, int *ocount, int **olist, int silent,                 */
/*      CCrandstate *rstate)                                                */
/*    RETURNS the quadrant k-nearest neighbor graph.                        */
/*      -see CCkdtree_k_nearest.                                            */
/*                                                                          */
/*  int CCkdtree_node_k_nearest (CCkdtree *kt, int ncount, int n, int k,    */
/*      CCdatagroup *dat, double *wcoord, int *list, CCrandstate *rstate)   */
/*    RETURNS the k nearest points to point n.                              */
/*      -The k points are return in list (and list must be allocated by     */
/*       calling routine.                                                   */
/*      -kt is a pointer to a CCkdtree previously built by                  */
/*       CCkdtree_build.                                                    */
/*                                                                          */
/*  int CCkdtree_node_quadrant_k_nearest (CCkdtree *kt, int ncount,         */
/*      int n, int k, CCdatagroup *dat, double *wcoord, int *list,          */
/*      CCrandstate *rstate)                                                */
/*    RETURNS the quadrant k nearest point to point n.                      */
/*      -see CCkdtree_node_k_nearest.                                       */
/*                                                                          */
/*  int CCkdtree_node_nearest (ktree *kt, int n, CCdatagroup *dat,          */
/*      double *wcoord)                                                     */
/*    RETURNS the nearest point to point n.                                 */
/*      -kt CANNOT be NULL.                                                 */
/*      -The point is returned as the function value. kt is a pointer       */
/*       to a CCkdtree (previously buildt by a call to CCkdtree_build)      */
/*                                                                          */
/*  int CCkdtree_fixed_radius_nearest (CCkdtree *kt, CCdatagroup *dat,      */
/*      double *wcoord, int n, double rad,                                  */
/*      int (*doit_fn) (int, int, void *), void *pass_param)                */
/*    ACTION: Calls the function doit_fn (n, a, void *), where a ranges     */
/*            over all points within distance rad of the point n. The       */
/*            void * field can be used to bundle a group of parmeters       */
/*            into pass_param that will be passed to doit_fn.               */
/*      -kt CANNOT be NULL.                                                 */
/*      -doit_fn can also call CCkdtree_fixed_radius_nearest (no globals    */
/*       are set by the function calls)                                     */
/*      -pass_param can be NULL or used to point to a structure with        */
/*       with parameters for doit_fn.                                       */
/*                                                                          */
/*  int CCkdtree_nearest_neighbor_tour (CCkdtree *kt, int ncount,           */
/*      int start, CCdatagroup *dat, int *outcycle, double *val,            */
/*      CCrandstate *rstate)                                                */
/*    -kt can be NULL.                                                      */
/*    -Node weights are not used.                                           */
/*    -start is the starting node for the tour.                             */
/*    -if outcycle is not NULL, then it should point to a array of          */
/*     length at least ncount (allocated by the calling routine). The       */
/*     cycle will be returned in the array in node node node format.        */
/*    -the length of the tour is return in val.                             */
/*                                                                          */
/*  int CCkdtree_nearest_neighbor_2match (CCkdtree *kt, int ncount,         */
/*      int start, CCdatagroup *dat, int *outmatch, double *val,            */
/*      CCrandstate *rstate)                                                */
/*    -Like CCkdtree_nearest_neighbor_tour. If outmatch is not NULL         */
/*     then it should point to an array of length at least 2*ncount.        */
/*                                                                          */
/*    NOTES:                                                                */
/*       If memory is tight, use CCkdtree_node_k_nearest to get the         */
/*    edges one node at a time. (CCkdtree_k_nearest () builds a hash        */
/*    table to avoid duplicate edges, and it will use 8 * nedges            */
/*    bytes.)                                                               */
/*       CCkdtree_node_nearest returns the nearest point as the             */
/*    function value; CCkdtree_fixed_radius_nearest returns 1 if            */
/*    doit_fn returns a nonzero value, otherwise it returns 0; all          */
/*    other routines return 0 if successful and 1 otherwise.                */
/*                                                                          */
/****************************************************************************/


#include "ccutil.h"

#define BIGDOUBLE (1e30)
#define NEAR_HEAP_CUTOFF 100  /* When to switch from list to heap       */
                              /* On an RS6000, the heap started winning */
                              /* at 100 (by 200 it was 3 times faster)  */

typedef struct shortedge {
    int end;
    double length;
} shortedge;

typedef struct intptr {
    int this;
    struct intptr *next;
} intptr;

#define Fedgelen(n1, n2)                                                     \
    (datw != (double *) NULL ?                                               \
      CCutil_dat_edgelen ((n1), (n2), dat)                                   \
            + datw[(n1)] + datw[(n2)] :                                      \
      CCutil_dat_edgelen ((n1), (n2), dat))

#define dtrunc(x) (((x)>0.0)?floor(x):ceil(x))

static void
    node_k_nearest_work (CCkdtree *thetree, CCdatagroup *dat, double *datw,
        CCkdnode *p, CCdheap *near_heap, int *heap_names, int *heap_count,
        int target, int num, shortedge *nearlist, double *worst_on_list,
        CCkdbnds *box),
    node_nearest_work (CCkdtree *thetree, CCdatagroup *dat, double *datw,
         CCkdnode *p, int target, double *ndist, int *nnode);
static int
    run_kdtree_k_nearest (CCkdtree *kt, int ncount, int k, CCdatagroup *dat,
        double *wcoord, int wantlist, int *ocount, int **olist, int doquad,
        int silent, CCrandstate *rstate),
    put_in_table (int i, int j, int *added, intptr **table,
        CCptrworld *intptr_world),
    q_run_it (CCkdtree *thetree, CCdatagroup *dat, double *datw, int *llist,
         int *lcount, int *list, int target, int num, CCkdbnds *box),
    run_kdtree_node_k_nearest (CCkdtree *thetree, CCdatagroup *dat,
         double *datw, int *list, int target, int num, CCkdbnds *box),
    ball_in_bounds (CCdatagroup *dat, CCkdbnds *bnds, int n, double dist),
    fixed_radius_nearest_work (CCkdtree *thetree, CCkdnode *p,
         int (*doit_fn) (int, int, void *),
         int target, double dist, CCdatagroup *dat, double *wcoord,
         double xtarget, double ytarget, double ztarget, void *pass_param);


CC_PTRWORLD_LIST_ROUTINES (intptr, int, intptralloc, intptr_bulk_alloc,
        intptrfree, intptr_listadd, intptr_listfree)
CC_PTRWORLD_LEAKS_ROUTINE (intptr, intptr_check_leaks, this, int)


int CCkdtree_k_nearest (CCkdtree *kt, int ncount, int k, CCdatagroup *dat,
        double *wcoord, int wantlist, int *ocount, int **olist,
        int silent, CCrandstate *rstate)
{
    return run_kdtree_k_nearest (kt, ncount, k, dat, wcoord,
                   wantlist, ocount, olist, 0, silent, rstate);
}

int CCkdtree_quadrant_k_nearest (CCkdtree *kt, int ncount, int k,
        CCdatagroup *dat, double *wcoord, int wantlist, int *ocount,
        int **olist, int silent, CCrandstate *rstate)
{
    return run_kdtree_k_nearest (kt, ncount, k, dat, wcoord,
                   wantlist, ocount, olist, 1, silent, rstate);
}

static int run_kdtree_k_nearest (CCkdtree *kt, int ncount, int k,
        CCdatagroup *dat, double *wcoord, int wantlist, int *ocount,
        int **olist, int doquad, int silent, CCrandstate *rstate)
{
    int rval = 0, i, n, total, onlist, added, ntotal = 0;
    int newtree = 0, goal = k, *list = (int *) NULL;
    intptr *ip, *ipnext, **table = (intptr **) NULL;
    CCkdtree localkt, *mykt;
    CCptrworld intptr_world;

    CCptrworld_init (&intptr_world);
    
    if (wcoord != (double *) NULL) {
        for (i = 0; i < ncount; i++) {
            if (wcoord[i] < -0.00000001) {
                fprintf (stderr, "Cannot run CCkdtree with neg w-weights\n");
                return 1;
            }
        }
    }

    if (wantlist) {
        *ocount = 0;
        *olist = (int *) NULL;
    }

    if (doquad) {
        if (dat->z) goal *= 8;
        else        goal *= 4;
    }

    if (kt == (CCkdtree *) NULL) {
        rval = CCkdtree_build (&localkt, ncount, dat, wcoord, rstate);
        CCcheck_rval (rval, "CCkdtree_build failed");
        mykt = &localkt;
        newtree = 1;
    } else {
        mykt = kt;
    }

    CC_MALLOC (table, ncount, intptr *);
    for (i = 0; i < ncount; i++) table[i] = (intptr *) NULL;
    CC_MALLOC (list, goal, int);

    for (n = 0; n < ncount; n++) {
        if (doquad) {
            rval = CCkdtree_node_quadrant_k_nearest (mykt, ncount, n, k, dat,
                                                     wcoord, list, rstate);
            CCcheck_rval (rval, "CCkdtree_node_quadrant_k_nearest failed");
        } else {
            rval = CCkdtree_node_k_nearest (mykt, ncount, n, k, dat, wcoord,
                                            list, rstate);
            CCcheck_rval (rval, "CCkdtree_node_k_nearest failed");
        }
        for (i = 0; i < goal; i++) {
            if (list[i] != -1) {
                if (put_in_table (n, list[i], &added, table, &intptr_world))  {
                    rval = 1; goto CLEANUP;
                } else {
                    ntotal += added;
                }
            }
        }
        if (!silent && n % 100000 == 99999) { printf ("."); fflush (stdout); }
    }
  
    if (!silent) { printf (" %d edges\n", ntotal); fflush (stdout); }

    if (wantlist) {
        int j = 0;
        CC_MALLOC (*olist, 2 * ntotal, int);
        *ocount = ntotal;
        for (i = 0; i < ncount; i++) {
            for (ip = table[i]; ip; ip = ipnext) {
                ipnext =  ip->next;
                (*olist)[j++] = i;
                (*olist)[j++] = ip->this;
                intptrfree (&intptr_world, ip);
            }
            table[i] = (intptr *) NULL;
        }
    } else {
        for (i = 0; i < ncount; i++) {
            intptr_listfree (&intptr_world, table[i]);
            table[i] = (intptr *) NULL;
        }
    }
    if (intptr_check_leaks (&intptr_world, &total, &onlist)) {
        fprintf (stderr, "WARNING: %d outstanding intptrs in kdnear\n",
                 total - onlist);
    }

CLEANUP:
    CCptrworld_delete (&intptr_world);
    CC_IFFREE(table, intptr *);
    CC_IFFREE (list, int);
    if (newtree) CCkdtree_free (&localkt);
    return rval;
}

static int put_in_table (int i, int j, int *added, intptr **table,
        CCptrworld *intptr_world)
{
    intptr *ip;

    if (j < i) {
        int temp;
        CC_SWAP(i, j, temp);
    }

    for (ip = table[i]; ip; ip = ip->next) {
        if (ip->this == j) {
            *added = 0;
            return 0;
        }
    }
    if (intptr_listadd (&table[i], j, intptr_world)) {
        *added = 0;
        return 1;
    }
    *added = 1;
    return 0;
}

int CCkdtree_node_quadrant_k_nearest (CCkdtree *kt, int ncount, int n, int k,
            CCdatagroup *dat, double *wcoord, int *list, CCrandstate *rstate)
{
    int rval = 0, i, lcount = 0, newtree = 0, iter, nrun;
    int *llist = (int *) NULL;
    CCkdbnds localbnds;
    CCkdtree localkt;
    CCkdtree *thetree;

    if (kt == (CCkdtree *) NULL) {
        rval = CCkdtree_build (&localkt, ncount, dat, wcoord, rstate);
        CCcheck_rval (rval, "CCkdtree_build failed");
        thetree = &localkt;
        newtree = 1;
    } else {
        thetree = kt;
    }

    CC_MALLOC (llist, k, int);

    if (dat->z) nrun = 2;
    else        nrun = 1;

    for (iter = 0; iter < nrun; iter++) {
        if (dat->z) {
            if (iter == 0) {
                localbnds.z[0] = dat->z[n];
                localbnds.z[1] = BIGDOUBLE;
            } else {
                localbnds.z[0] = -BIGDOUBLE;
                localbnds.z[1] = dat->z[n];
            }
        }

        localbnds.x[0] = dat->x[n];
        localbnds.x[1] = BIGDOUBLE;
        localbnds.y[0] = dat->y[n];
        localbnds.y[1] = BIGDOUBLE;
        rval = q_run_it (thetree, dat, wcoord, llist, &lcount, list, n, k,
                         &localbnds);
        CCcheck_rval (rval, "q_run_it failed");

        localbnds.y[0] = -BIGDOUBLE;
        localbnds.y[1] = dat->y[n];
        rval = q_run_it (thetree, dat, wcoord, llist, &lcount, list, n, k,
                         &localbnds);
        CCcheck_rval (rval, "q_run_it failed");

        localbnds.x[0] = -BIGDOUBLE;
        localbnds.x[1] = dat->x[n];
        localbnds.y[0] = dat->y[n];
        localbnds.y[1] = BIGDOUBLE;
        rval = q_run_it (thetree, dat, wcoord, llist, &lcount, list, n, k,
                         &localbnds);
        CCcheck_rval (rval, "q_run_it failed");

        localbnds.y[0] = -BIGDOUBLE;
        localbnds.y[1] = dat->y[n];
        rval = q_run_it (thetree, dat, wcoord, llist, &lcount, list, n, k,
                         &localbnds);
        CCcheck_rval (rval, "q_run_it failed");
    }


    if (dat->z) {
        for (i = lcount; i < (8 * k); i++) list[i] = -1;
    } else {
        for (i = lcount; i < (4 * k); i++) list[i] = -1;
    }

CLEANUP:
    CC_FREE (llist, int);
    if (newtree) CCkdtree_free (&localkt);
    return rval;
}

static int q_run_it (CCkdtree *thetree, CCdatagroup *dat, double *datw,
        int *llist, int *lcount, int *list, int target, int num, CCkdbnds *box)
{
    int rval = 0, i, j;

    rval = run_kdtree_node_k_nearest (thetree, dat, datw, llist, target, num,
                                      box);
    CCcheck_rval (rval, "run_kdtree_node_k_nearest failed");

    for (i = 0; i < num; i++) {
        if (llist[i] != -1) {
            for (j = 0; j < *lcount; j++) {
                if (list[j] == llist[i]) break;
            }
            if (j == *lcount) list[(*lcount)++] = llist[i];
        }
    }

CLEANUP:
    return rval;
}

int CCkdtree_node_k_nearest (CCkdtree *kt, int ncount, int n, int k,
        CCdatagroup *dat, double *wcoord, int *list, CCrandstate *rstate)
{
    int rval = 0, newtree = 0;
    CCkdtree localkt, *thetree;

    if (kt == (CCkdtree *) NULL) {
        rval = CCkdtree_build (&localkt, ncount, dat, wcoord, rstate);
        CCcheck_rval (rval, "CCkdtree_build failed");
        thetree = &localkt;
        newtree = 1;
    } else {
        thetree = kt;
    }

    rval = run_kdtree_node_k_nearest (thetree, dat, wcoord, list, n, k,
                                      (CCkdbnds *) NULL);
    CCcheck_rval (rval, "run_kdtree_node_k_nearest failed");

CLEANUP:
    if (newtree) CCkdtree_free (&localkt);
    return rval;
}

static int run_kdtree_node_k_nearest (CCkdtree *thetree, CCdatagroup *dat,
        double *datw, int *list, int target, int num, CCkdbnds *box)
{
    int rval = 0, i;
    int *heap_names =  (int *) NULL, heap_count = 0, have_heap = 0;
    double diff, worst_on_list = BIGDOUBLE;
    CCkdnode *p, *lastp;
    CCdheap near_heap;
    shortedge *nearlist = (shortedge *) NULL;

    if (num >= NEAR_HEAP_CUTOFF) {
        rval = CCutil_dheap_init (&near_heap, num);
        CCcheck_rval (rval, "CCutil_dheap_init failed");
        have_heap = 1;
        CC_MALLOC (heap_names, num, int);
    } else {
        CC_MALLOC (nearlist, num+1, shortedge);
        for (i = 0; i < num; i++) nearlist[i].length = BIGDOUBLE;
        nearlist[num].length = -BIGDOUBLE;
    }

/*
    To do top down search just use:

        node_k_nearest_work (thetree->root);
*/

    p = thetree->bucketptr[target];
    node_k_nearest_work (thetree, dat, datw, p, &near_heap, heap_names,
                         &heap_count, target, num, nearlist, &worst_on_list,
                         box);
    while (1) {
        lastp = p;
        p = p->father;
        if (p == (CCkdnode *) NULL) break;
        switch (p->cutdim) {
        case 0:
            diff = p->cutval - dat->x[target];
            if (lastp == p->loson) {    /* So target is on low side */
               if (worst_on_list > dtrunc(diff))
                   if (box == (CCkdbnds *) NULL || p->cutval <= box->x[1])
                       node_k_nearest_work (thetree, dat, datw, p->hison,
                              &near_heap, heap_names, &heap_count, target,
                              num, nearlist, &worst_on_list, box);
            } else {
               if (worst_on_list > dtrunc(-diff))
                   if (box == (CCkdbnds *) NULL || p->cutval >= box->x[0])
                       node_k_nearest_work (thetree, dat, datw, p->loson,
                              &near_heap, heap_names, &heap_count, target,
                              num, nearlist, &worst_on_list, box);
            }
            break;
        case 1:
            diff = p->cutval - dat->y[target];
            if (lastp == p->loson) {
               if (worst_on_list > dtrunc(diff))
                   if (box == (CCkdbnds *) NULL || p->cutval <= box->y[1])
                       node_k_nearest_work (thetree, dat, datw, p->hison,
                              &near_heap, heap_names, &heap_count, target,
                              num, nearlist, &worst_on_list, box);
            } else {
               if (worst_on_list > dtrunc(-diff))
                   if (box == (CCkdbnds *) NULL || p->cutval >= box->y[0])
                       node_k_nearest_work (thetree, dat, datw, p->loson,
                              &near_heap, heap_names, &heap_count, target,
                              num, nearlist, &worst_on_list, box);
            }
            break;
        case 2:
            diff = p->cutval - dat->z[target];
            if (lastp == p->loson) {
               if (worst_on_list > dtrunc(diff))
                   if (box == (CCkdbnds *) NULL || p->cutval <= box->z[1])
                       node_k_nearest_work (thetree, dat, datw, p->hison,
                              &near_heap, heap_names, &heap_count, target,
                              num, nearlist, &worst_on_list, box);
            } else {
               if (worst_on_list > dtrunc(-diff))
                   if (box == (CCkdbnds *) NULL || p->cutval >= box->z[0])
                       node_k_nearest_work (thetree, dat, datw, p->loson,
                              &near_heap, heap_names, &heap_count, target,
                              num, nearlist, &worst_on_list, box);
            }

            break;
        case 3:
            if (lastp == p->loson) {
                if (worst_on_list > p->cutval + datw[target])
                    node_k_nearest_work (thetree, dat, datw, p->hison,
                              &near_heap, heap_names, &heap_count, target,
                              num, nearlist, &worst_on_list, box);
            } else {
                node_k_nearest_work (thetree, dat, datw, p->loson, &near_heap,
                              heap_names, &heap_count, target, num, nearlist,
                              &worst_on_list, box);
            }
            break;
        }
        if (datw == (double *) NULL && p->bnds &&
               ball_in_bounds (dat, p->bnds, target, worst_on_list))
              /* Doing extra check for box with quad-nearest appears to slow */
              /* things down.                                                */
            break;
    }

    if (num >= NEAR_HEAP_CUTOFF) {
        if (heap_count < num) {
            if (box == (CCkdbnds *) NULL) {
                fprintf (stderr, "WARNING: There do not exist %d neighbors\n",
                         num);
            }
            for (i = 0; i < heap_count; i++) list[i] = heap_names[i];
            for (; i < num; i++) list[i] = -1;
        } else {
            for (i = 0; i < num; i++) list[i] = heap_names[i];
        }
    } else {
        int ntot = 0;
        for (i = 0; i < num; i++) {
            if (nearlist[i].length < BIGDOUBLE)
                list[ntot++] = nearlist[i].end;
        }
        if (ntot < num) {
            if (box == (CCkdbnds *) NULL) {
                fprintf (stderr, "WARNING: There do not exist %d neighbors\n",
                         num);
            }
            for (i = ntot; i < num; i++) list[i] = -1;
        }
    }

CLEANUP:
    if (have_heap) CCutil_dheap_free (&near_heap);
    CC_IFFREE (heap_names, int);
    CC_IFFREE (nearlist, shortedge);
    return rval;
}

static void node_k_nearest_work (CCkdtree *thetree, CCdatagroup *dat,
        double *datw, CCkdnode *p, CCdheap *near_heap, int *heap_names,
        int *heap_count, int target, int num, shortedge *nearlist,
        double *worst_on_list, CCkdbnds *box)
{
    int i, h, k;
    double val, thisx, thisdist;

    if (p->bucket) {
        if (num >= NEAR_HEAP_CUTOFF) {
            for (i = p->lopt; i <= p->hipt; i++) {
                if (thetree->perm[i] != target) {
                    if (box == (CCkdbnds *) NULL ||
                       (dat->x[thetree->perm[i]] >= box->x[0] &&
                        dat->x[thetree->perm[i]] <= box->x[1] &&
                        dat->y[thetree->perm[i]] >= box->y[0] &&
                        dat->y[thetree->perm[i]] <= box->y[1] &&
                        (!dat->z || (
                          dat->z[thetree->perm[i]] >= box->z[0] &&
                          dat->z[thetree->perm[i]] <= box->z[1]))
                        )) {
                        thisdist = Fedgelen (thetree->perm[i], target);
                        if (*heap_count < num) {
                            near_heap->key[*heap_count] = -thisdist;
                            heap_names[*heap_count] = thetree->perm[i];
                            /* this can't fail since the heap is big enough */
                            CCutil_dheap_insert (near_heap, *heap_count);
                            (*heap_count)++;
                        } else if (*worst_on_list > thisdist) {
                            h = CCutil_dheap_deletemin (near_heap);
                            heap_names[h] = thetree->perm[i];
                            near_heap->key[h] = -thisdist;
                            /* this can't fail since the heap is big enough */
                            CCutil_dheap_insert (near_heap, h);
                            h = CCutil_dheap_findmin (near_heap);
                            *worst_on_list = -near_heap->key[h];
                        }
                    }
                }
            }
        } else {
            for (i = p->lopt; i <= p->hipt; i++) {
                if (thetree->perm[i] != target) {
                    if (box == (CCkdbnds *) NULL ||
                       (dat->x[thetree->perm[i]] >= box->x[0] &&
                        dat->x[thetree->perm[i]] <= box->x[1] &&
                        dat->y[thetree->perm[i]] >= box->y[0] &&
                        dat->y[thetree->perm[i]] <= box->y[1] &&
                        (!dat->z || (
                          dat->z[thetree->perm[i]] >= box->z[0] &&
                          dat->z[thetree->perm[i]] <= box->z[1]))
                        )) {
                        thisdist = Fedgelen (thetree->perm[i], target);
                        if (*worst_on_list > thisdist) {
                            for (k = 0; nearlist[k+1].length > thisdist; k++) {
                                nearlist[k].end = nearlist[k + 1].end;
                                nearlist[k].length = nearlist[k + 1].length;
                            }
                            nearlist[k].length = thisdist;
                            nearlist[k].end = thetree->perm[i];
                            *worst_on_list = nearlist[0].length;
                        }
                    }
                }
            }
        }
    } else {
        val = p->cutval;
        switch (p->cutdim) {
        case 0:
            thisx = dat->x[target];
            if (thisx < val) {
                node_k_nearest_work (thetree, dat, datw, p->loson, near_heap,
                        heap_names, heap_count, target, num, nearlist,
                        worst_on_list, box);
                /* Truncation for floating point coords */
                if (*worst_on_list > dtrunc(val - thisx))
                    if (box == (CCkdbnds *) NULL || val >= box->x[0])
                        node_k_nearest_work (thetree, dat, datw, p->hison,
                               near_heap, heap_names, heap_count, target,
                               num, nearlist, worst_on_list, box);
            } else {
                node_k_nearest_work (thetree, dat, datw, p->hison, near_heap,
                               heap_names, heap_count, target, num, nearlist,
                               worst_on_list, box);
                if (*worst_on_list > dtrunc(thisx - val))
                    if (box == (CCkdbnds *) NULL || val <= box->x[1])
                        node_k_nearest_work (thetree, dat, datw, p->loson,
                               near_heap, heap_names, heap_count, target,
                               num, nearlist, worst_on_list, box);
            }
            break;
        case 1:
            thisx = dat->y[target];
            if (thisx < val) {
                node_k_nearest_work (thetree, dat, datw, p->loson, near_heap,
                               heap_names, heap_count, target, num, nearlist,
                               worst_on_list, box);
                if (*worst_on_list > dtrunc(val - thisx))
                    if (box == (CCkdbnds *) NULL || val >= box->y[0])
                        node_k_nearest_work (thetree, dat, datw, p->hison,
                               near_heap, heap_names, heap_count, target,
                               num, nearlist, worst_on_list, box);
            } else {
                node_k_nearest_work (thetree, dat, datw, p->hison, near_heap,
                               heap_names, heap_count, target, num, nearlist,
                               worst_on_list, box);
                if (*worst_on_list > dtrunc(thisx - val))
                    if (box == (CCkdbnds *) NULL || val <= box->y[1])
                        node_k_nearest_work (thetree, dat, datw, p->loson,
                               near_heap, heap_names, heap_count, target,
                               num, nearlist, worst_on_list, box);
            }
            break;
        case 2:
            thisx = dat->z[target];
            if (thisx < val) {
                node_k_nearest_work (thetree, dat, datw, p->loson, near_heap,
                               heap_names, heap_count, target, num, nearlist,
                               worst_on_list, box);
                if (*worst_on_list > dtrunc(val - thisx))
                    if (box == (CCkdbnds *) NULL || val >= box->z[0])
                        node_k_nearest_work (thetree, dat, datw, p->hison,
                               near_heap, heap_names, heap_count, target,
                               num, nearlist, worst_on_list, box);
            } else {
                node_k_nearest_work (thetree, dat, datw, p->hison, near_heap,
                               heap_names, heap_count, target, num, nearlist,
                               worst_on_list, box);
                if (*worst_on_list > dtrunc(thisx - val))
                    if (box == (CCkdbnds *) NULL || val <= box->z[1])
                        node_k_nearest_work (thetree, dat, datw, p->loson,
                               near_heap, heap_names, heap_count, target,
                               num, nearlist, worst_on_list, box);
            }
            break;
        case 3:
            thisx = datw[target];
            node_k_nearest_work (thetree, dat, datw, p->loson, near_heap,
                               heap_names, heap_count, target, num, nearlist,
                               worst_on_list, box);
            if (*worst_on_list > val + thisx)
                node_k_nearest_work (thetree, dat, datw, p->hison, near_heap,
                               heap_names, heap_count, target, num, nearlist,
                               worst_on_list, box);
            break;
        }
    }
}


int CCkdtree_node_nearest (CCkdtree *kt, int n, CCdatagroup *dat,
        double *wcoord)
{
    int nnode;
    double diff = 0.0, ndist = BIGDOUBLE;
    CCkdnode *p, *lastp;
    CCkdtree *thetree = (CCkdtree *) NULL;

    if (kt == (CCkdtree *) NULL) {
        fprintf (stderr, "ERROR: kt  NULL in CCkdtree_node_nearest)\n");
        return n;
    }

    thetree = kt;

    ndist = BIGDOUBLE;
    nnode = n;

/*
    To do top down search just use:

        node_nearest_work (kt->root);
        return nnode;
*/

    p = kt->bucketptr[n];
    node_nearest_work (thetree, dat, wcoord, p, n, &ndist, &nnode);
    while (1) {
        lastp = p;
        p = p->father;
        if (p == (CCkdnode *) NULL) break;

        switch (p->cutdim) {
        case 0:
            diff = p->cutval - dat->x[n];
            break;
        case 1:
            diff = p->cutval - dat->y[n];
            break;
        case 2:
            diff = p->cutval - dat->z[n];
            break;
        case 3:
            break;
        }

        if (p->cutdim < 3) {
            if (lastp == p->loson) {
                if (ndist > dtrunc(diff))
                   node_nearest_work (thetree, dat, wcoord, p->hison, n,
                                      &ndist, &nnode);
            } else {
                if (ndist > dtrunc(-diff))
                   node_nearest_work (thetree, dat, wcoord, p->loson, n,
                                      &ndist, &nnode);
            }
        } else {
            if (lastp == p->loson) {
                if (ndist > p->cutval + wcoord[n])
                    node_nearest_work (thetree, dat, wcoord, p->hison, n,
                                      &ndist, &nnode);
            } else {
                node_nearest_work (thetree, dat, wcoord, p->loson, n,
                                   &ndist, &nnode);
            }
        }


        if (wcoord == (double *) NULL && p->bnds &&
               ball_in_bounds (dat, p->bnds, n, ndist))
            break;
    }
    return nnode;
}

static int ball_in_bounds (CCdatagroup *dat, CCkdbnds *bnds, int n,
        double dist)
{
    if (dtrunc(dat->x[n] - bnds->x[0]) < dist ||
        dtrunc(bnds->x[1] - dat->x[n]) < dist ||
        dtrunc(dat->y[n] - bnds->y[0]) < dist ||
        dtrunc(bnds->y[1] - dat->y[n]) < dist ||
        (dat->z && (dtrunc(dat->z[n] - bnds->z[0]) < dist ||
                    dtrunc(bnds->z[1] - dat->z[n]) < dist)))
        return 0;
    return 1;
}

static void node_nearest_work (CCkdtree *thetree, CCdatagroup *dat,
        double *datw, CCkdnode *p, int target, double *ndist, int *nnode)
{
    int i;
    double val, thisx = 0.0, thisdist;

    if (!p->empty) {
        if (p->bucket) {
            for (i = p->lopt; i <= p->hipt; i++) {
                if (thetree->perm[i] != target) {
                    thisdist = Fedgelen (thetree->perm[i], target);
                    if (*ndist > thisdist) {
                        *ndist = thisdist;
                        *nnode = thetree->perm[i];
                    }
                }
            }
        } else {
            val = p->cutval;
            switch (p->cutdim) {
            case 0:
                thisx = dat->x[target];
                break;
            case 1:
                thisx = dat->y[target];
                break;
            case 2:
                thisx = dat->z[target];
                break;
            case 3:
                thisx = datw[target];
                break;
            }
            if (p->cutdim < 3) {
                if (thisx < val) {
                    node_nearest_work (thetree, dat, datw, p->loson, target,
                                       ndist, nnode);
                    if (*ndist > dtrunc(val - thisx))
                        node_nearest_work (thetree, dat, datw, p->hison,
                                       target, ndist, nnode);
                } else {
                    node_nearest_work (thetree, dat, datw, p->hison, target,
                                       ndist, nnode);
                    if (*ndist > dtrunc(thisx - val))
                        node_nearest_work (thetree, dat, datw, p->loson,
                                       target, ndist, nnode);
                }
            } else {
                node_nearest_work (thetree, dat, datw, p->loson, target, ndist,
                                       nnode);
                if (*ndist > val + thisx)
                    node_nearest_work (thetree, dat, datw, p->hison, target,
                                       ndist, nnode);
            }
        }
    }
}

int CCkdtree_fixed_radius_nearest (CCkdtree *kt, CCdatagroup *dat,
        double *datw, int n, double rad, int (*doit_fn) (int, int, void *),
        void *pass_param)
{
    int target;
    CCkdnode *p, *lastp;
    double dist, diff, xtarget, ytarget, ztarget = 0.0;

    if (kt == (CCkdtree *) NULL) {
        fprintf (stderr, "ERROR: fixed_radius_nearest needs CCkdtree\n");
        return 0;
    }

    dist = (double) rad;
    target = n;
    xtarget = dat->x[target];
    ytarget = dat->y[target];
    if (dat->z) ztarget = dat->z[target];
        
    p = kt->bucketptr[target];
    if (fixed_radius_nearest_work (kt, p, doit_fn, target, dist, dat, datw,
                                   xtarget, ytarget, ztarget, pass_param))
        return 1;

    if (datw) {
        double wdist = dist - datw[target];
        while (1) {
            lastp = p;
            p = p->father;
            if (p == (CCkdnode *) NULL) return 0;

            if (p->cutdim == 0)      diff = p->cutval - xtarget;
            else if (p->cutdim == 1) diff = p->cutval - ytarget;
            else if (p->cutdim == 2) diff = p->cutval - ztarget;
            else                     diff = p->cutval;

            if (lastp == p->loson) {
                if (wdist > dtrunc(diff)) {
                    if (fixed_radius_nearest_work (kt, p->hison, doit_fn,
                              target, dist, dat, datw, xtarget, ytarget,
                              ztarget, pass_param))
                        return 1;
                }
            } else {
                if (wdist > dtrunc(-diff)) {
                    if (fixed_radius_nearest_work (kt, p->loson, doit_fn,
                              target, dist, dat, datw, xtarget, ytarget,
                              ztarget, pass_param))
                        return 1;
                }
            }

            if (p->bnds &&  /* ball_in_bounds */
              !(dtrunc(xtarget - p->bnds->x[0]) < wdist ||
                dtrunc(p->bnds->x[1] - xtarget) < wdist ||
                dtrunc(ytarget - p->bnds->y[0]) < wdist ||
                dtrunc(p->bnds->y[1] - ytarget) < wdist ||
                (dat->z && (dtrunc(ztarget - p->bnds->z[0]) < wdist ||
                            dtrunc(p->bnds->z[1] - ztarget) < wdist))))
                return 0;
        }
    } else {
        while (1) {
            lastp = p;
            p = p->father;
            if (p == (CCkdnode *) NULL) return 0;

            if (p->cutdim == 0)      diff = p->cutval - xtarget;
            else if (p->cutdim == 1) diff = p->cutval - ytarget;
            else                     diff = p->cutval - ztarget;

            if (lastp == p->loson) {
                if (dist > dtrunc(diff)) {
                    if (fixed_radius_nearest_work (kt, p->hison, doit_fn,
                              target, dist, dat, datw, xtarget, ytarget,
                              ztarget, pass_param))
                        return 1;
                }
            } else {
                if (dist > dtrunc(-diff)) {
                    if (fixed_radius_nearest_work (kt, p->loson, doit_fn,
                              target, dist, dat, datw, xtarget, ytarget,
                              ztarget, pass_param))
                        return 1;
                }
            }
            if (p->bnds &&  /* ball_in_bounds */
                !(dtrunc(xtarget - p->bnds->x[0]) < dist ||
                  dtrunc(p->bnds->x[1] - xtarget) < dist ||
                  dtrunc(ytarget - p->bnds->y[0]) < dist ||
                  dtrunc(p->bnds->y[1] - ytarget) < dist ||
                  (dat->z && (dtrunc(ztarget - p->bnds->z[0]) < dist ||
                              dtrunc(p->bnds->z[1] - ztarget) < dist))))
                return 0;
        }
    }
}

static int fixed_radius_nearest_work (CCkdtree *thetree, CCkdnode *p,
            int (*doit_fn) (int, int, void *), int target, double dist,
            CCdatagroup *dat, double *datw,  double xtarget, double ytarget,
            double ztarget, void *pass_param)
{
    int i, c;
    double val, thisx = 0.0, thisdist;

    if (p->empty) return 0;

    if (p->bucket) {
        for (i = p->lopt; i <= p->hipt; i++) {
            c = thetree->perm[i];
            if (c != target) {
                thisdist = Fedgelen (c, target);
                if (thisdist < dist) {
                    if (doit_fn (target, c, pass_param)) {
                        return 1;
                    }
                }
            }
        }
        return 0;
    } else {
        if (datw) {
            double wdist = dist - datw[target];

            val = p->cutval;
            switch (p->cutdim) {
            case 0:
                thisx = xtarget;
                break;
            case 1:
                thisx = ytarget;
                break;
            case 2:
                thisx = ztarget;
                break;
            case 3:
                if (fixed_radius_nearest_work (thetree, p->loson, doit_fn,
                      target, dist, dat, datw, xtarget, ytarget, ztarget,
                      pass_param)) {
                    return 1;
                }
                if (p->cutval <= wdist) {
                    if (fixed_radius_nearest_work (thetree, p->hison, doit_fn,
                         target, dist, dat, datw, xtarget, ytarget, ztarget,
                         pass_param)) {
                        return 1;
                    }
                }
                return 0;
            default:
                return 0;
            }
            if (thisx < val) {
                if (fixed_radius_nearest_work (thetree, p->loson, doit_fn,
                        target, dist, dat, datw, xtarget, ytarget, ztarget,
                        pass_param)) {
                    return 1;
                }
                if (wdist > dtrunc(val - thisx)) {
                    if (fixed_radius_nearest_work (thetree, p->hison, doit_fn,
                            target, dist, dat, datw, xtarget, ytarget, ztarget,
                            pass_param)) {
                        return 1;
                    }
                }
            } else {
                if (fixed_radius_nearest_work (thetree, p->hison, doit_fn,
                        target, dist, dat, datw, xtarget, ytarget, ztarget,
                        pass_param)) {
                    return 1;
                }
                if (wdist > dtrunc(thisx - val)) {
                    if (fixed_radius_nearest_work (thetree, p->loson, doit_fn,
                            target, dist, dat, datw, xtarget, ytarget, ztarget,
                            pass_param)) {
                        return 1;
                    }
                }
            }
        } else {
            val = p->cutval;
            switch (p->cutdim) {
            case 0:
                thisx = xtarget;
                break;
            case 1:
                thisx = ytarget;
                break;
            case 2:
                thisx = ztarget;
                break;
            case 3:
            default:
                fprintf (stderr, "ERROR: split on w without node weights\n");
                return 1;
            }
            if (thisx < val) {
                if (fixed_radius_nearest_work (thetree, p->loson, doit_fn,
                        target, dist, dat, datw, xtarget, ytarget, ztarget,
                        pass_param)) {
                    return 1;
                }
                if (dist > dtrunc(val - thisx)) {
                    if (fixed_radius_nearest_work (thetree, p->hison, doit_fn,
                            target, dist, dat, datw, xtarget, ytarget, ztarget,
                            pass_param)) {
                        return 1;
                    }
                }
            } else {
                if (fixed_radius_nearest_work (thetree, p->hison, doit_fn,
                        target, dist, dat, datw, xtarget, ytarget, ztarget,
                        pass_param)) {
                    return 1;
                }
                if (dist > dtrunc(thisx - val)) {
                    if (fixed_radius_nearest_work (thetree, p->loson, doit_fn,
                            target, dist, dat, datw, xtarget, ytarget, ztarget,
                            pass_param)) {
                        return 1;
                    }
                }
            }
        }
    }
    return 0;
}

int CCkdtree_nearest_neighbor_tour (CCkdtree *kt, int ncount, int start,
         CCdatagroup *dat, int *outcycle, double *val, CCrandstate *rstate)
{
    int rval = 0, newtree = 0, i, current, next;
    double len;
    CCkdtree localkt, *mykt;

    if (ncount < 3) {
        fprintf (stderr, "Cannot find tour in an %d node graph\n", ncount);
        rval = 1; goto CLEANUP;
    }

    if (kt == (CCkdtree *) NULL) {
        rval = CCkdtree_build (&localkt, ncount, dat, (double *) NULL, rstate);
        CCcheck_rval (rval, "CCkdtree_build failed");
        mykt = &localkt;
        newtree = 1;
    } else {
        mykt = kt;
    }

    len = 0.0;
    current = start;
    if (outcycle != (int *) NULL) outcycle[0] = start;
    for (i = 1; i < ncount; i++) {
        CCkdtree_delete (mykt, current);
        next = CCkdtree_node_nearest (mykt, current, dat, (double *) NULL);
        if (outcycle != (int *) NULL) outcycle [i] = next;
        len += (double) CCutil_dat_edgelen (current, next, dat);
        current = next;
    }
    len += (double) CCutil_dat_edgelen (current, start, dat);
    *val = len;

CLEANUP:
    if (newtree) CCkdtree_free (&localkt);
    if (kt) CCkdtree_undelete_all (kt, ncount);
    return rval;
}

int CCkdtree_nearest_neighbor_2match (CCkdtree *kt, int ncount, int start,
         CCdatagroup *dat, int *outmatch, double *val, CCrandstate *rstate)
{
    int rval = 0, i, j, cur, next;
    int count = 0, cyccount = 0, newtree = 0;
    char *mark = (char *) NULL;
    double len, szeit;
    CCkdtree localkt, *mykt;

    if (ncount < 3) {
        fprintf (stderr, "Cannot find 2-matching in an %d node graph\n",
                 ncount);
        rval = 1; goto CLEANUP;
    }

    if (kt == (CCkdtree *) NULL) {
        rval = CCkdtree_build (&localkt, ncount, dat, (double *) NULL, rstate);
        CCcheck_rval (rval, "CCkdtree_build failed");
        mykt = &localkt;
        newtree = 1;
    } else {
        mykt = kt;
    }

    CC_MALLOC (mark, ncount, char);
    for (i = 0 ; i < ncount; i++) mark[i] = 0;

    printf ("Grow nearest neighbor 2-matching from node %d\n", start);
    fflush (stdout);
    szeit = CCutil_zeit ();
    len = 0.0;

    while (count < ncount) {
        for (j = start; j < ncount && mark[j]; j++);
        if (j == ncount) {
            for (j = 0; j < start && mark[j]; j++);
            if (j == start) {
                fprintf (stderr, "ERROR in near-2match\n");
                rval = 1; goto CLEANUP;
            }
        }
        start = j;
        mark[start] = 1;
        CCkdtree_delete (mykt, start);
        next = CCkdtree_node_nearest (mykt, start, dat, (double *) NULL);
        mark[next] = 1;
        len += (double) CCutil_dat_edgelen (start, next, dat);
        if (outmatch != (int *) NULL) {
            outmatch[2 * count] = start;
            outmatch[(2 * count) + 1] = next;
        }
        count++;
        CCkdtree_delete (mykt, next);
        cur = CCkdtree_node_nearest (mykt, next, dat, (double *) NULL);
        len += (double) CCutil_dat_edgelen (next, cur, dat);
        if (outmatch != (int *) NULL) {
            outmatch[2 * count] = next;
            outmatch[(2 * count) + 1] = cur;
        }
        count++;
        CCkdtree_undelete (mykt, start);
        while (cur != start && count < ncount - 3) {
            mark[cur] = 1;
            CCkdtree_delete (mykt, cur);
            next = CCkdtree_node_nearest (mykt, cur, dat, (double *) NULL);
            len += (double) CCutil_dat_edgelen (cur, next, dat);
            if (outmatch != (int *) NULL) {
                outmatch[2 * count] = cur;
                outmatch[(2 * count) + 1] = next;
            }
            count++;
            cur = next;
        }
        CCkdtree_delete (mykt, start);

        if (cur != start) {   /* Not enough nodes for another circuit */
            while (count < ncount - 1) {
                mark[cur] = 1;
                CCkdtree_delete (mykt, cur);
                next = CCkdtree_node_nearest (mykt, cur, dat, (double *) NULL);
                len += (double) CCutil_dat_edgelen (cur, next, dat);
                if (outmatch != (int *) NULL) {
                    outmatch[2 * count] = cur;
                    outmatch[(2 * count) + 1] = next;
                }
                count++;
                cur = next;
            }
            len += (double) CCutil_dat_edgelen (cur, start, dat);
            if (outmatch != (int *) NULL) {
                outmatch[2 * count] = cur;
                outmatch[(2 * count) + 1] = start;
            }
            count++;
        }
        cyccount++;
    }

    *val = len;
    printf ("%d cycles in 2-matching\n", cyccount);
    printf ("Running time for Nearest Neighbor 2-match: %.2f\n",
                                                  CCutil_zeit () - szeit);
    fflush (stdout);

CLEANUP:
    if (newtree) CCkdtree_free (&localkt);
    if (kt) CCkdtree_undelete_all (kt, ncount);
    CC_IFFREE (mark, char);
    return rval;
}

/****************************************************************************/
/*                                                                          */
/*                    INPUT/OUTPUT ROUTINES                                 */
/*                                                                          */
/*                           TSP CODE                                       */
/*                                                                          */
/* Written by:  Applegate, Bixby, Chvatal, and Cook                         */
/* Date: February 13, 1995                                                  */
/*       September 30, 1997 (dave)                                          */
/*                                                                          */
/*    EXPORTED FUNCTIONS:                                                   */
/*                                                                          */
/*  CC_SFILE *CCutil_sopen (char *f, char *s)                               */
/*      Opens a file for buffered binary I/O.  The buffered binary I/O      */
/*      routines using CC_SFILE's attempt to be machine independent,        */
/*      and only compatible with themselves.  Comparable to the stdio       */
/*      routine fopen().  If the file already exists and is being           */
/*      opened for output, the old file is renamed by prepending an O       */
/*      to is name.                                                         */
/*    f - the filename to open.  "stdin" means descriptor 0, "stdout"       */
/*        descriptor 1, and "stderr" descriptor 2.  "-" means               */
/*        descriptor 0 or 1, depending on wither the file is opened         */
/*        for reading or writing.                                           */
/*    s - the mode to open, either "r" for input, or "w" for output.        */
/*    returns a pointer to the opened file, or NULL if there is an          */
/*        error.                                                            */
/*                                                                          */
/*  CC_SFILE *CCutil_sdopen (int d, char *s)                                */
/*      Opens a descriptor for buffered binary I/O.  The buffered binary    */
/*      I/O routines using CC_SFILE's attempt to be machine independent,    */
/*      and only compatible with themselves.  Comparable to the stdio       */
/*      routine fdopen().                                                   */
/*    d - the descriptor to open.                                           */
/*    s - the mode to open, either "r" for input, "w" for output, or        */
/*        "rw" for both input and output.                                   */
/*    returns a pointer to the opened file, or NULL if there is an          */
/*        error.                                                            */
/*                                                                          */
/*  int CCutil_swrite (CC_SFILE *f, char *buf, int size)                    */
/*      writes to a buffered binary I/O file.                               */
/*    f - the CC_SFILE to write to                                          */
/*    buf - the data to write                                               */
/*    size - the number of bytes to write.                                  */
/*    returns 0 if succesful, nonzero if failure.                           */
/*                                                                          */
/*  int CCutil_swrite_bits (CC_SFILE *f, int x, int xbits)                  */
/*      writes bits to a buffered binary I/O file.                          */
/*    f - the CC_SFILE to write to                                          */
/*    x - an int containing the data to write                               */
/*    xbits - the number of bits to write.  The lowest order xbits          */
/*            bits of x will be written.                                    */
/*    returns 0 if succesful, nonzero if failure.                           */
/*                                                                          */
/*  int CCutil_swrite_ubits (CC_SFILE *f, unsigned int x, int xbits)        */
/*      writes bits to a buffered binary I/O file.                          */
/*    f - the CC_SFILE to write to                                          */
/*    x - an unsigned int int containing the data to write                  */
/*    xbits - the number of bits to write.  The lowest order xbits          */
/*            bits of x will be written.                                    */
/*    returns 0 if succesful, nonzero if failure.                           */
/*                                                                          */
/*  int CCutil_swrite_char (CC_SFILE *f, char x)                            */
/*      writes a char to a buffered binary I/O file.                        */
/*    f - the CC_SFILE to write to                                          */
/*    x - the char to write                                                 */
/*    returns 0 if succesful, nonzero if failure.                           */
/*                                                                          */
/*  int CCutil_swrite_string (CC_SFILE *f, const char *s)                   */
/*      writes a string to a buffered binary I/O file.                      */
/*    f - the CC_SFILE to write to                                          */
/*    s - the string to write.  The array of characters in s up to and      */
/*        including the first NULL are written.                             */
/*    returns 0 if succesful, nonzero if failure.                           */
/*                                                                          */
/*  int CCutil_swrite_short (CC_SFILE *f, short x)                          */
/*      writes a short to a buffered binary I/O file.                       */
/*    f - the CC_SFILE to write to                                          */
/*    x - the short to write                                                */
/*    returns 0 if succesful, nonzero if failure.                           */
/*                                                                          */
/*  int CCutil_swrite_ushort (CC_SFILE *f, unsigned short x)                */
/*      writes an unsigned short to a buffered binary I/O file.             */
/*    f - the CC_SFILE to write to                                          */
/*    x - the unsigned short to write                                       */
/*    returns 0 if succesful, nonzero if failure.                           */
/*                                                                          */
/*  int CCutil_swrite_int (CC_SFILE *f, int x)                              */
/*      writes an int to a buffered binary I/O file.                        */
/*    f - the CC_SFILE to write to                                          */
/*    x - the int to write                                                  */
/*    returns 0 if succesful, nonzero if failure.                           */
/*                                                                          */
/*  int CCutil_swrite_uint (CC_SFILE *f, unsigned int x)                    */
/*      writes an unsigned int to a buffered binary I/O file.               */
/*    f - the CC_SFILE to write to                                          */
/*    x - the unsigned int to write                                         */
/*    returns 0 if succesful, nonzero if failure.                           */
/*                                                                          */
/*  int CCutil_swrite_double (CC_SFILE *f, double x)                        */
/*      writes a double to a buffered binary I/O file.                      */
/*    f - the CC_SFILE to write to                                          */
/*    x - the double to write                                               */
/*    returns 0 if succesful, nonzero if failure.                           */
/*                                                                          */
/*  int CCutil_sread (CC_SFILE *f, char *buf, int size)                     */
/*      reads from a buffered binary I/O file.                              */
/*    f - the CC_SFILE to read from.                                        */
/*    buf - a buffer in which to store the data read.  buf should have      */
/*          space for size characters.                                      */
/*    size - the number of bytes to read.                                   */
/*    returns 0 if succesful, nonzero if failure.                           */
/*                                                                          */
/*  int CCutil_sread_bits (CC_SFILE *f, int *x, int xbits)                  */
/*      reads bits from a buffered binary I/O file.                         */
/*    f - the CC_SFILE to read from.                                        */
/*    x - on return, will contain the bits read (in the low-order           */
/*        xbits bits).                                                      */
/*    xbits - the number of bits read.                                      */
/*    returns 0 if succesful, nonzero if failure.                           */
/*                                                                          */
/*  int CCutil_sread_ubits (CC_SFILE *f, unsigned int *x, int xbits)        */
/*      reads bits from a buffered binary I/O file.                         */
/*    f - the CC_SFILE to read from.                                        */
/*    x - on return, will contain the bits read (in the low-order           */
/*        xbits bits).                                                      */
/*    xbits - the number of bits read.                                      */
/*    returns 0 if succesful, nonzero if failure.                           */
/*                                                                          */
/*  int CCutil_sread_char (CC_SFILE *f, char *x)                            */
/*      reads a char from a buffered binary I/O file.                       */
/*    f - the CC_SFILE to read from.                                        */
/*    x - on return, will contain the char read                             */
/*    returns 0 if succesful, nonzero if failure.                           */
/*                                                                          */
/*  int CCutil_sread_string (CC_SFILE *f, char *x, int maxlen)              */
/*      reads a string from a buffered binary I/O file.                     */
/*    f - the CC_SFILE to read from.                                        */
/*    x - on return, will contain the string read.                          */
/*    maxlen - the maximum number of characters to read.                    */
/*    returns 0 if succesful, nonzero if failure.                           */
/*                                                                          */
/*  int CCutil_sread_short (CC_SFILE *f, short *x)                          */
/*      reads a short from a buffered binary I/O file.                      */
/*    f - the CC_SFILE to read from.                                        */
/*    x - on return, will contain the short read                            */
/*    returns 0 if succesful, nonzero if failure.                           */
/*                                                                          */
/*  int CCutil_sread_ushort (CC_SFILE *f, unsigned short *x)                */
/*      reads an unsigned short from a buffered binary I/O file.            */
/*    f - the CC_SFILE to read from.                                        */
/*    x - on return, will contain the unsigned short read                   */
/*    returns 0 if succesful, nonzero if failure.                           */
/*                                                                          */
/*  int CCutil_sread_short_r (CC_SFILE *f, short *x)                        */
/*      reads a reversed short from a buffered binary I/O file.             */
/*    f - the CC_SFILE to read from.                                        */
/*    x - on return, will contain the short read                            */
/*    returns 0 if succesful, nonzero if failure.                           */
/*    CCutil_sread_short_r is provided for compatability with some          */
/*    binary files written by other tools which use a different byte        */
/*    order.                                                                */
/*                                                                          */
/*  int CCutil_sread_int (CC_SFILE *f, int *x)                              */
/*      reads an int from a buffered binary I/O file.                       */
/*    f - the CC_SFILE to read from.                                        */
/*    x - on return, will contain the int read                              */
/*    returns 0 if succesful, nonzero if failure.                           */
/*                                                                          */
/*  int CCutil_sread_uint (CC_SFILE *f, unsigned int *x)                    */
/*      reads an unsigned int from a buffered binary I/O file.              */
/*    f - the CC_SFILE to read from.                                        */
/*    x - on return, will contain the unsigned int read                     */
/*    returns 0 if succesful, nonzero if failure.                           */
/*                                                                          */
/*  int CCutil_sread_int_r (CC_SFILE *f, int *x)                            */
/*      reads a reversed int from a buffered binary I/O file.               */
/*    f - the CC_SFILE to read from.                                        */
/*    x - on return, will contain the int read                              */
/*    returns 0 if succesful, nonzero if failure.                           */
/*    CCutil_sread_int_r is provided for compatability with some            */
/*    binary files written by other tools which use a different byte        */
/*    order.                                                                */
/*                                                                          */
/*  int CCutil_sread_double (CC_SFILE *f, double *x)                        */
/*      reads a double from a buffered binary I/O file.                     */
/*    f - the CC_SFILE to read from.                                        */
/*    x - on return, will contain the double read                           */
/*    returns 0 if succesful, nonzero if failure.                           */
/*                                                                          */
/*  int CCutil_sread_double_r (CC_SFILE *f, double *x)                      */
/*      reads a reversed double from a buffered binary I/O file.            */
/*    f - the CC_SFILE to read from.                                        */
/*    x - on return, will contain the double read                           */
/*    returns 0 if succesful, nonzero if failure.                           */
/*    CCutil_sread_double_r is provided for compatability with some         */
/*    binary files written by other tools which use a different byte        */
/*    order.                                                                */
/*                                                                          */
/*  int CCutil_sflush (CC_SFILE *f)                                         */
/*      flushes the buffer of a buffered binary I/O file.                   */
/*    f - the CC_SFILE to flush                                             */
/*    returns 0 if succesful, nonzero if failure.                           */
/*                                                                          */
/*  int CCutil_stell (CC_SFILE *f)                                          */
/*      returns the current location in a buffered binary I/O file.         */
/*      Comparable to the stdio function ftell().                           */
/*    f - the CC_SFILE                                                      */
/*    returns the current location, or -1 for failure.                      */
/*                                                                          */
/*  int CCutil_sseek (CC_SFILE *f, int offset)                              */
/*      changes the current location in a buffered binary I/O file.         */
/*      Comparable to the stdio function fseek().                           */
/*    f - the CC_SFILE                                                      */
/*    offset - a value returned by CCutil_stell().                          */
/*    returns 0 for success, nonzero for failure.                           */
/*                                                                          */
/*  int CCutil_srewind (CC_SFILE *f)                                        */
/*      changes the current location in a buffered binary I/O file to       */
/*      the beginning.  Comparable to the stdio function rewind().          */
/*    f - the CC_SFILE                                                      */
/*    returns 0 for success, nonzero for failure.                           */
/*                                                                          */
/*  int CCutil_sclose (CC_SFILE *f)                                         */
/*      closes a CC_SFILE.                                                  */
/*    f - the CC_SFILE to close                                             */
/*    returns 0 for success, nonzero for failure.                           */
/*                                                                          */
/*  int CCutil_sbits (unsigned int x)                                       */
/*      computes the number of bits necessary to represent all              */
/*      nonnegative integers <= x                                           */
/*    x - a number                                                          */
/*    returns the number of bits necessary to represent x.                  */
/*                                                                          */
/*  int CCutil_sdelete_file (const char *fname)                             */
/*      deletes a file.                                                     */
/*    fname - the file to delete                                            */
/*    returns 0 for success, nonzero for failure.                           */
/*                                                                          */
/*  int CCutil_sdelete_file_backup (const char *fname)                      */
/*      deletes the backup file for fname (created if fname was             */
/*      overwritten by CCutil_sopen).                                       */
/*    fname - the file name whose backup is to be deleted.                  */
/*    returns 0 for success, nonzero for failure.                           */
/*                                                                          */
/*  CC_SFILE *CCutil_snet_open (char *h, unsigned short p)                  */
/*      Opens a network connection to a port on a remote host               */
/*    h - the name of the host to connect to                                */
/*    p - the port on the host to connect to                                */
/*    returns a CC_SFILE (opened for input and output) for buffered         */
/*            binary I/O to the specified port on the remote host,          */
/*            or NULL if there is a failure.                                */
/*    Only exists if CC_NETREADY is defined                                 */
/*                                                                          */
/*  CC_SFILE *CCutil_snet_receive (CC_SPORT *s)                             */
/*      Accepts a network connection on a port.                             */
/*    s - the CC_SPORT to accept a connection from.  Must be the            */
/*        returned result of a successfull CCutil_snet_listen call.         */
/*    returns a CC_SFILE (opened for input and output) for buffered         */
/*        binary I/O on the specified port, or NULL if there is a           */
/*        failure.                                                          */
/*    Only exists if CC_NETREADY is defined                                 */
/*                                                                          */
/*  CC_SPORT *CCutil_snet_listen (unsigned short p)                         */
/*      Prepares to accept network connections on a port.                   */
/*    p - the port on which to accept connections.                          */
/*    returns a CC_SPORT for accepting connections on the specified         */
/*        port.  This return value is passed to CCutil_snet_receive to      */
/*        accept a connection.  Returns NULL if there is a failure.         */
/*    Only exists if CC_NETREADY is defined                                 */
/*                                                                          */
/*  void CCutil_snet_unlisten (CC_SPORT *s)                                 */
/*      Ceases accepting network connections from an CC_SPORT.              */
/*    s - the CC_SPORT to close.                                            */
/*    Only exists if CC_NETREADY is defined                                 */
/*                                                                          */
/****************************************************************************/

#include "ccutil.h"

static CC_SFILE
    *sopen_write (const char *f),
    *sopen_read (const char *f),
    *sdopen (int t),
    *sdopen_write (int t),
    *sdopen_read (int t),
    *sdopen_readwrite (int t);

static int
    swrite_buffer (CC_SFILE *f),
    sread_buffer (CC_SFILE *f),
    prepare_write (CC_SFILE *f),
    prepare_read (CC_SFILE *f);

static void
    sinit (CC_SFILE *s);



/* VERSION A3 */

/* CCutil_sopen interprets filenames "stdin" as descriptor 0, "stdout" as
 * descriptor 1, and "stderr" as descriptor 2.  "-" is interpreted as
 * 0 or 1, depending on whether the file is opened for reading or writing.
 *
 * CCutil_sclose doesn't close descriptors 0, 1, and 2.
 */

/* When writing, written data extends from buffer[0] bit 7 through
 * buffer[chars_in_buffer-1] bit bits_in_last_char.  Empty space extends
 * from buffer[chars_in_buffer-1] bit bits_in_last_char-1 through
 * buffer[CC_SBUFFER_SIZE-1] bit 0.
 *
 * When reading, read data extends from buffer[0] bit 7 through
 * buffer[current_buffer_char] bit bits_in_last_char.  unread data
 * extends from buffer[current_buffer_char] bit bits_in_last_char-1
 * through buffer[chars_in_buffer-1] bit 0.  Empty space extends from
 * buffer[chars_in_buffer] bit 7 through buffer[CC_SBUFFER_SIZE-1] bit 0.
 */

/* If the routines detect an error, they return -1.
 */

#define SREAD 1
#define SWRITE 2
#define SRW_EMPTY 3
#define SRW_READ 4
#define SRW_WRITE 5

#define TFILE 1
#define TDESC 2
#define TNET 3

#define NBITMASK(n) ((1<<(n))-1)
#define BITRANGE(x,start,length) (((x) >> (start)) & NBITMASK(length))
#define BITS_PER_CHAR (8)

#ifndef O_BINARY
#define O_BINARY 0
#endif
#ifndef O_EXCL
#define O_EXCL 0
#endif

CC_SFILE *CCutil_sopen (const char *f, const char *s)
{
    if (strcmp (s, "r") == 0) {
        return sopen_read (f);
    } else if (strcmp (s, "w") == 0) {
        return sopen_write (f);
    } else {
        fprintf (stderr, "Need to specify read/write in CCutil_sopen\n");
        return (CC_SFILE *) NULL;
    }
}

CC_SFILE *CCutil_sdopen (int d, const char *s)
{
    if (strcmp (s, "r") == 0) {
        return sdopen_read (d);
    } else if (strcmp (s, "w") == 0) {
        return sdopen_write (d);
    } else if (strcmp (s, "rw") == 0) {
        return sdopen_readwrite (d);
    } else {
        fprintf (stderr, "Need to specify read/write in CCutil_sdopen\n");
        return (CC_SFILE *) NULL;
    }
}

static CC_SFILE *sopen_write (const char *f)
{
    CC_SFILE *s = (CC_SFILE *) NULL;
    int t;
    char fbuf[CC_SFNAME_SIZE];
    char fbuf_N[CC_SFNAME_SIZE + 32];
    char fbuf_Nx[CC_SFNAME_SIZE + 64];

    strncpy (fbuf, f, sizeof (fbuf) - 12);
    fbuf[sizeof (fbuf) - 12] = '\0';
    sprintf (fbuf_N, "N%s", fbuf);
    sprintf (fbuf_Nx, "N%s~", fbuf);


    if (strcmp (f, "stdout") == 0 || strcmp (f, "-") == 0) {
        s = sdopen_write (1);
    } else if (strcmp (f, "stderr") == 0) {
        s = sdopen_write (2);
    } else {
        t = open (fbuf_N, O_WRONLY | O_CREAT | O_BINARY | O_EXCL, 0644);
        if (t == -1 && errno == EEXIST) {
            fprintf (stderr, "%s already exists, renaming to %s\n",
                          fbuf_N, fbuf_Nx);
            if (rename (fbuf_N, fbuf_Nx)) {
                perror (fbuf_Nx);
                fprintf (stderr, "Couldn't rename %s to %s\n", fbuf_N,
                         fbuf_Nx);
                return (CC_SFILE *) NULL;
            }
            t = open (fbuf_N, O_WRONLY | O_CREAT | O_BINARY | O_EXCL, 0644);
        }
        if (t == -1) {
            perror (fbuf_N);
            fprintf (stderr, "Couldn't open %s for output\n", fbuf_N);
            return (CC_SFILE *) NULL;
        }
        s = sdopen_write (t);
        if (!s) {
            close (t);
        } else {
            s->type = TFILE;
        }
    }
    if (s) {
        strncpy (s->fname, fbuf, sizeof (s->fname));
        s->fname[sizeof (s->fname)-1] = '\0';
    }
    return s;
}

static CC_SFILE *sopen_read (const char *f)
{
    CC_SFILE *s = (CC_SFILE *) NULL;
    int t;

    if (strcmp (f, "stdin") == 0 || strcmp (f, "-") == 0) {
        s = sdopen_read (0);
    } else {
        t = open (f, O_RDONLY | O_BINARY, 0644);
        if (t == -1) {
            perror (f);
            fprintf (stderr, "Couldn't open for input\n");
            s = (CC_SFILE *) NULL;
        }
        s = sdopen_read (t);
        if (!s) {
            close (t);
        } else {
            s->type = TFILE;
        }
    }
    if (s) {
        strncpy (s->fname, f, sizeof (s->fname));
        s->fname[sizeof (s->fname)-1] = '\0';
    }
    return s;
}

static CC_SFILE *sdopen (int t)
{
    CC_SFILE *s = (CC_SFILE *) NULL;

    if (t < 0) {
        fprintf (stderr, "Invalid descriptor %d\n", t);
        return (CC_SFILE *) NULL;
    }

    s = CC_SAFE_MALLOC (1, CC_SFILE);
    if (s == (CC_SFILE *) NULL) {
        return (CC_SFILE *) NULL;
    }
    sinit (s);

    s->desc = t;
    s->type = TDESC;
    sprintf (s->fname, "descriptor %d", t);
    return s;
}

static CC_SFILE *sdopen_write (int t)
{
    CC_SFILE *s = (CC_SFILE *) NULL;

    s = sdopen (t);
    if (s) {
        s->status = SWRITE;
    }

    return s;
}

static CC_SFILE *sdopen_read (int t)
{
    CC_SFILE *s = (CC_SFILE *) NULL;

    s = sdopen (t);
    if (s) {
        s->status = SREAD;
    }

    return s;
}

static CC_SFILE *sdopen_readwrite (int t)
{
    CC_SFILE *s = (CC_SFILE *) NULL;

    s = sdopen (t);
    if (s) {
        s->status = SRW_EMPTY;
    }

    return s;
}

int CCutil_swrite (CC_SFILE *f, char *buf, int size)
{
    int i;

    for (i=0; i<size; i++) {
        if (CCutil_swrite_char (f, buf[i])) return -1;
    }
    return 0;
}

int CCutil_swrite_bits (CC_SFILE *f, int x, int xbits)
{
    if (x < 0) {
        fprintf (stderr, "CCutil_swrite_bits cannot write negative numbers\n");
        return -1;
    }
    return CCutil_swrite_ubits (f, (unsigned int) x, xbits);
}

int CCutil_swrite_ubits (CC_SFILE *f, unsigned int x, int xbits)
{
    int getbits;
    unsigned int v;

    if (prepare_write (f)) return -1;

    while (xbits) {
        if (f->bits_in_last_char == 0) {
            if (f->chars_in_buffer == CC_SBUFFER_SIZE) {
                if (swrite_buffer (f)) return -1;
            }
            f->buffer[f->chars_in_buffer++] = 0;
            f->bits_in_last_char = BITS_PER_CHAR;
        }
        getbits = f->bits_in_last_char;
        if (getbits > xbits)
            getbits = xbits;
        xbits -= getbits;
        f->bits_in_last_char -= getbits;
        v = BITRANGE (x, xbits, getbits);
        f->buffer[f->chars_in_buffer - 1] =
            (unsigned int) f->buffer[f->chars_in_buffer - 1] |
            (unsigned int) (v << f->bits_in_last_char);
    }
    return 0;
}

int CCutil_swrite_char (CC_SFILE *f, char x)
{
    unsigned char ux = (unsigned char) x;
    
    if (prepare_write (f)) return -1;

    f->bits_in_last_char = 0;
    if (f->chars_in_buffer + 1 > CC_SBUFFER_SIZE) {
        if (swrite_buffer (f)) return -1;
    }
    f->buffer[f->chars_in_buffer++] = ((unsigned int) ux) & 0xff;
    return 0;
}

int CCutil_swrite_string (CC_SFILE *f, const char *s)
{
    int rval;

    while (*s) {
        rval = CCutil_swrite_char (f, *s);
        if (rval)
            return rval;
        s++;
    }
    CCutil_swrite_char (f, (unsigned char) 0);
    return 0;
}

int CCutil_swrite_short (CC_SFILE *f, short x)
{
    return CCutil_swrite_ushort (f, (unsigned short) x);
}

int CCutil_swrite_ushort (CC_SFILE *f, unsigned short x)
{
    if (prepare_write (f)) return -1;

    f->bits_in_last_char = 0;
    if (f->chars_in_buffer + 2 > CC_SBUFFER_SIZE) {
        if (swrite_buffer (f)) return -1;
    }

    f->buffer[f->chars_in_buffer++] = (((unsigned int) x) >> 8) & 0xff;
    f->buffer[f->chars_in_buffer++] = ((unsigned int) x) & 0xff;
    return 0;
}

int CCutil_swrite_int (CC_SFILE *f, int x)
{
    return CCutil_swrite_uint (f, (unsigned int) x);
}

int CCutil_swrite_uint (CC_SFILE *f, unsigned int x)
{
    if (prepare_write (f)) return -1;

    f->bits_in_last_char = 0;
    if (f->chars_in_buffer + 4 > CC_SBUFFER_SIZE) {
        if (swrite_buffer (f)) return -1;
    }

    f->buffer[f->chars_in_buffer++] = (((unsigned int) x) >> 24) & 0xff;
    f->buffer[f->chars_in_buffer++] = (((unsigned int) x) >> 16) & 0xff;
    f->buffer[f->chars_in_buffer++] = (((unsigned int) x) >> 8) & 0xff;
    f->buffer[f->chars_in_buffer++] = ((unsigned int) x) & 0xff;
    return 0;
}

int CCutil_swrite_double (CC_SFILE *f, double x)
{
    unsigned short e;
    unsigned int m1;
    unsigned int m2;

    e = 128;

    if (x < 0) {
        e = (unsigned int) e + 256;
        x = -x;
    }

    if (x >= 1.0) {
#define MUNCH_HI_EXP(x,e,v,lv) if (x >= v) {e = (unsigned int) e + lv; x *= 1/v;}
        MUNCH_HI_EXP(x,e,18446744073709551616.0,64);
        MUNCH_HI_EXP(x,e,4294967296.0,32);
        MUNCH_HI_EXP(x,e,65536.0, 16);
        MUNCH_HI_EXP(x,e,256.0, 8);
        MUNCH_HI_EXP(x,e,16.0, 4);
        MUNCH_HI_EXP(x,e,4.0, 2);
        MUNCH_HI_EXP(x,e,2.0, 1);
#undef MUNCH_HI_EXP
        x /= 2;
        e = (unsigned int) e + 1;
    } else if (x < 0.5) {
#define MUNCH_LO_EXP(x,e,v,lv) if (x < 1/v) {e = (unsigned int) e - lv; x *= v;}
        MUNCH_LO_EXP(x,e,18446744073709551616.0,64);
        MUNCH_LO_EXP(x,e,4294967296.0,32);
        MUNCH_LO_EXP(x,e,65536.0, 16);
        MUNCH_LO_EXP(x,e,256.0, 8);
        MUNCH_LO_EXP(x,e,16.0, 4);
        MUNCH_LO_EXP(x,e,4.0, 2);
        MUNCH_LO_EXP(x,e,2.0, 1);
#undef MUNCH_LP_EXP
    }
    x *= 4294967296.0;
    m1 = (unsigned int) x;
    m2 = (unsigned int) ((x - m1) * 4294967296.0);
    if (CCutil_swrite_ushort (f, e)) return -1;
    if (CCutil_swrite_uint (f, m1)) return -1;
    if (CCutil_swrite_uint (f, m2)) return -1;
    return 0;
}

int CCutil_sread (CC_SFILE *f, char *buf, int size)
{
    int i;

    for (i=0; i<size; i++) {
        if (CCutil_sread_char (f, &buf[i])) return -1;
    }
    return 0;
}

int CCutil_sread_bits (CC_SFILE *f, int *x, int xbits)
{
    unsigned int ux = 0;
    int rval;
    
    rval = CCutil_sread_ubits (f, &ux, xbits);
    *x = (int) ux;
    return rval;
}

int CCutil_sread_ubits (CC_SFILE *f, unsigned int *x, int xbits)
{
    int getbits;
    unsigned int v;

    if (prepare_read (f)) return -1;

    *x = 0;
    while (xbits) {
        if (f->bits_in_last_char == 0) {
            if (f->current_buffer_char + 1 == f->chars_in_buffer) {
                if (sread_buffer (f)) return -1;
            }
            f->current_buffer_char++;
            f->bits_in_last_char = BITS_PER_CHAR;
        }
        getbits = f->bits_in_last_char;
        if (getbits > xbits)
            getbits = xbits;
        f->bits_in_last_char -= getbits;
        xbits -= getbits;
        v = BITRANGE ((unsigned int) f->buffer[f->current_buffer_char],
                      f->bits_in_last_char, getbits);
        *x |= v << xbits;
    }
    return 0;
}

int CCutil_sread_char (CC_SFILE *f, char *x)
{
    if (prepare_read (f)) return -1;

    f->bits_in_last_char = 0;
    if (f->current_buffer_char + 1 == f->chars_in_buffer) {
        if (sread_buffer (f)) return -1;
    }
    *x = (char) (f->buffer[++f->current_buffer_char]);
    return 0;
}

int CCutil_sread_string (CC_SFILE *f, char *x, int maxlen)
{
    int i, rval;

    maxlen--;
    for (i = 0; i < maxlen;  i++, x++) {
        rval = CCutil_sread_char (f, x);
        if (rval)
            return rval;
        if (*x == 0)
            return 0;
    }
    *x = 0;
    return 0;
}

int CCutil_sread_short (CC_SFILE *f, short *x)
{
    unsigned short ux = 0;
    int rval;

    rval = CCutil_sread_ushort (f, &ux);
    *x = (short) ux;
    return rval;
}

int CCutil_sread_ushort (CC_SFILE *f, unsigned short *x)
{
    if (prepare_read (f)) return -1;

    f->bits_in_last_char = 0;
    if (f->current_buffer_char + 1 == f->chars_in_buffer) {
        if (sread_buffer (f)) return -1;
    }
    *x = ((unsigned int) f->buffer[++f->current_buffer_char]) << 8;
    if (f->current_buffer_char + 1 == f->chars_in_buffer) {
        if (sread_buffer (f)) return -1;
    }
    *x = (unsigned int) *x | ((unsigned int) f->buffer[++f->current_buffer_char]);
    return 0;
}

int CCutil_sread_short_r (CC_SFILE *f, short *x)
{
    unsigned short ux = 0;
    
    if (prepare_read (f)) return -1;

    f->bits_in_last_char = 0;
    if (f->current_buffer_char + 1 == f->chars_in_buffer) {
        if (sread_buffer (f)) return -1;
    }
    ux = ((unsigned short) f->buffer[++f->current_buffer_char]);
    if (f->current_buffer_char + 1 == f->chars_in_buffer) {
        if (sread_buffer (f)) return -1;
    }
    ux = (unsigned int) ux | ((unsigned int) f->buffer[++f->current_buffer_char]) << 8;
    *x = (short) ux;
    return 0;
}

int CCutil_sread_int (CC_SFILE *f, int *x)
{
    unsigned int ux = 0;
    int rval;

    rval = CCutil_sread_uint (f, &ux);
    *x = (int) ux;
    return rval;
}

int CCutil_sread_uint (CC_SFILE *f, unsigned int *x)
{
    if (prepare_read (f)) return -1;

    f->bits_in_last_char = 0;
    if (f->current_buffer_char + 1 == f->chars_in_buffer) {
        if (sread_buffer (f)) return -1;
    }
    *x = ((unsigned int) f->buffer[++f->current_buffer_char]) << 24;
    if (f->current_buffer_char + 1 == f->chars_in_buffer) {
        if (sread_buffer (f)) return -1;
    }
    *x |= ((unsigned int) f->buffer[++f->current_buffer_char]) << 16;
    if (f->current_buffer_char + 1 == f->chars_in_buffer) {
        if (sread_buffer (f)) return -1;
    }
    *x |= ((unsigned int) f->buffer[++f->current_buffer_char]) << 8;
    if (f->current_buffer_char + 1 == f->chars_in_buffer) {
        if (sread_buffer (f)) return -1;
    }
    *x |= ((unsigned int) f->buffer[++f->current_buffer_char]);
    return 0;
}

int CCutil_sread_int_r (CC_SFILE *f, int *x)
{
    unsigned int ux = 0;
    
    if (prepare_read (f)) return -1;

    f->bits_in_last_char = 0;
    if (f->current_buffer_char + 1 == f->chars_in_buffer) {
        if (sread_buffer (f)) return -1;
    }
    ux = ((unsigned int) f->buffer[++f->current_buffer_char]);
    if (f->current_buffer_char + 1 == f->chars_in_buffer) {
        if (sread_buffer (f)) return -1;
    }
    ux |= ((unsigned int) f->buffer[++f->current_buffer_char]) << 8;
    if (f->current_buffer_char + 1 == f->chars_in_buffer) {
        if (sread_buffer (f)) return -1;
    }
    ux |= ((unsigned int) f->buffer[++f->current_buffer_char]) << 16;
    if (f->current_buffer_char + 1 == f->chars_in_buffer) {
        if (sread_buffer (f)) return -1;
    }
    ux |= ((unsigned int) f->buffer[++f->current_buffer_char]) << 24;
    *x = (int) ux;
    return 0;
}

int CCutil_sread_double (CC_SFILE *f, double *x)
{
    unsigned short e;
    unsigned int m1;
    unsigned int m2;

    if (CCutil_sread_ushort (f, &e)) return -1;
    if (CCutil_sread_uint (f, &m1)) return -1;
    if (CCutil_sread_uint (f, &m2)) return -1;

    *x = ((m2 / 4294967296.0) + m1) / 4294967296.0;

    if ((unsigned int) e >= 256) {
        *x = -*x;
        e = (unsigned int) e - 256;
    }

    if ((unsigned int) e > 128) {
#define UNMUNCH_HI_EXP(x,e,v,lv) if ((unsigned int) e >= (unsigned int) (128 + lv)) \
                                     {e = (unsigned int) e - lv; x *= v;}
        UNMUNCH_HI_EXP(*x,e,18446744073709551616.0,64);
        UNMUNCH_HI_EXP(*x,e,4294967296.0,32);
        UNMUNCH_HI_EXP(*x,e,65536.0, 16);
        UNMUNCH_HI_EXP(*x,e,256.0, 8);
        UNMUNCH_HI_EXP(*x,e,16.0, 4);
        UNMUNCH_HI_EXP(*x,e,4.0, 2);
        UNMUNCH_HI_EXP(*x,e,2.0, 1);
#undef UNMUNCH_HI_EXP
    } else if ((unsigned int) e < 128) {
#define UNMUNCH_LO_EXP(x,e,v,lv) if ((unsigned int) e <= (unsigned int) (128 - lv)) \
                                     {e = (unsigned int) e + lv; x *= 1/v;}
        UNMUNCH_LO_EXP(*x,e,18446744073709551616.0,64);
        UNMUNCH_LO_EXP(*x,e,4294967296.0,32);
        UNMUNCH_LO_EXP(*x,e,65536.0, 16);
        UNMUNCH_LO_EXP(*x,e,256.0, 8);
        UNMUNCH_LO_EXP(*x,e,16.0, 4);
        UNMUNCH_LO_EXP(*x,e,4.0, 2);
        UNMUNCH_LO_EXP(*x,e,2.0, 1);
#undef UNMUNCH_LO_EXP
    }

    return 0;
}

int CCutil_sread_double_r (CC_SFILE *f, double *x)
{
    unsigned short e;
    unsigned int m1;
    unsigned int m2;
    short se;
    int sm1;
    int sm2;

    if (CCutil_sread_short_r (f, &se)) return -1;
    if (CCutil_sread_int_r (f, &sm1)) return -1;
    if (CCutil_sread_int_r (f, &sm2)) return -1;
    e = (unsigned short) se;
    m1 = (unsigned int) sm1;
    m2 = (unsigned int) sm2;

    *x = ((m2 / 4294967296.0) + m1) / 4294967296.0;

    if ((unsigned int) e >= 256) {
        *x = -*x;
        e = (unsigned int) e - 256;
    }

    if ((unsigned int) e > 128) {
#define UNMUNCH_HI_EXP(x,e,v,lv) if ((unsigned int) e >= (unsigned int) (128 + lv)) \
                                     {e = (unsigned int) e - lv; x *= v;}
        UNMUNCH_HI_EXP(*x,e,18446744073709551616.0,64);
        UNMUNCH_HI_EXP(*x,e,4294967296.0,32);
        UNMUNCH_HI_EXP(*x,e,65536.0, 16);
        UNMUNCH_HI_EXP(*x,e,256.0, 8);
        UNMUNCH_HI_EXP(*x,e,16.0, 4);
        UNMUNCH_HI_EXP(*x,e,4.0, 2);
        UNMUNCH_HI_EXP(*x,e,2.0, 1);
#undef UNMUNCH_HI_EXP
    } else if ((unsigned int) e < 128) {
#define UNMUNCH_LO_EXP(x,e,v,lv) if ((unsigned int) e <= (unsigned int) (128 - lv)) \
                                     {e = (unsigned int) e + lv; x *= 1/v;}
        UNMUNCH_LO_EXP(*x,e,18446744073709551616.0,64);
        UNMUNCH_LO_EXP(*x,e,4294967296.0,32);
        UNMUNCH_LO_EXP(*x,e,65536.0, 16);
        UNMUNCH_LO_EXP(*x,e,256.0, 8);
        UNMUNCH_LO_EXP(*x,e,16.0, 4);
        UNMUNCH_LO_EXP(*x,e,4.0, 2);
        UNMUNCH_LO_EXP(*x,e,2.0, 1);
#undef UNMUNCH_LO_EXP
    }

    return 0;
}

int CCutil_sflush (CC_SFILE *f)
{
    int rval;
    
    if (f == (CC_SFILE *) NULL) {
        rval = -1;
    } else if (f->status == SREAD || f->status == SRW_READ) {
        f->bits_in_last_char = 0;
        rval = 0;
    } else if (f->status == SWRITE || f->status == SRW_WRITE) {
        rval = swrite_buffer (f);
    } else if (f->status == SRW_EMPTY) {
        rval = 0;
    } else {
        fprintf (stderr, "Buffer %s has invalid status %d\n", f->fname,
                 f->status);
        rval = -1;
    }
        
    return rval;
}

int CCutil_stell (CC_SFILE *f)
{
    if (!f) return -1;
    f->bits_in_last_char = 0;
    if (f->status == SREAD) {
        return f->pos - f->chars_in_buffer + f->current_buffer_char + 1;
    } else if (f->status == SWRITE) {
        return f->pos + f->chars_in_buffer;
    } else if (f->status == SRW_EMPTY || f->status == SRW_READ ||
               f->status == SRW_WRITE) {
        fprintf (stderr, "Cannot CCutil_stell for a r/w CC_SFILE\n");
        return -1;
    } else {
        fprintf (stderr, "Buffer %s has invalid status %d\n", f->fname,
                 f->status);
        return -1;
    }
}

int CCutil_sseek (CC_SFILE *f, int offset)
{
    int curloc;

    if (!f) return -1;
    if (CCutil_sflush (f)) return -1;
    curloc = CCutil_stell (f);
    if (curloc < 0) return curloc;
    if (curloc == offset) return 0;
    if (lseek (f->desc, offset, SEEK_SET) < 0) {
        perror (f->fname);
        fprintf (stderr, "Unable to lseek on %s\n", f->fname);
        return -1;
    }
    f->chars_in_buffer = 0;
    f->current_buffer_char = -1;
    f->pos = offset;

    return 0;
}

int CCutil_srewind (CC_SFILE *f)
{
    return CCutil_sseek (f, 0);
}

int CCutil_sclose (CC_SFILE *f)
{
    int retval = 0;
    char fbuf_O[CC_SFNAME_SIZE + 32];
    char fbuf_N[CC_SFNAME_SIZE + 32];

    if (!f) return -1;

    if ((f->status == SWRITE || f->status == SRW_WRITE) &&
        f->chars_in_buffer) {
        if (swrite_buffer (f)) retval = -1;
    }

    if (f->desc >= 3) {
        if (close (f->desc)) {
            perror ("close");
            fprintf (stderr, "Unable to close swrite file %s\n", f->fname);
            retval = -1;
        }
        if (f->status == SWRITE && f->type == TFILE) {
            sprintf (fbuf_N, "N%s", f->fname);
            sprintf (fbuf_O, "O%s", f->fname);
            rename (f->fname, fbuf_O);
            if (rename (fbuf_N, f->fname)) {
                perror (f->fname);
                fprintf (stderr, "Couldn't rename %s to %s\n",
                                               fbuf_N, f->fname);
                retval = -1;
            }
        }
    }

    CC_FREE (f, CC_SFILE);

    return retval;
}

static int swrite_buffer (CC_SFILE *f)
{
    char *p;
    int nleft;
    int n;

    if (!f) return -1;
    if (f->status != SWRITE && f->status != SRW_WRITE &&
        f->status != SRW_EMPTY) {
        fprintf (stderr, "%s not open for output\n", f->fname);
        return -1;
    }

    p = (char *) f->buffer;
    nleft = f->chars_in_buffer;
    while (nleft) {
        n = (int) write (f->desc, p, nleft);
        if (n == -1) {
            if (errno == EINTR) {
                fprintf (stderr, "swrite_buffer interrupted, retrying\n");
                continue;
            }
            perror ("write");
            fprintf (stderr, "swrite_buffer of %d chars to %s failed\n", nleft,
                     f->fname);
            return -1;
        }
        nleft -= n;
        p += n;
        f->pos += n;
    }
    f->bits_in_last_char = 0;
    f->chars_in_buffer = 0;
    return 0;
}

static int sread_buffer (CC_SFILE *f)
{
    int n;

    if (!f) return -1;
    if (f->status != SREAD && f->status != SRW_READ &&
        f->status != SRW_EMPTY) {
        fprintf (stderr, "%s not open for input\n", f->fname);
        return -1;
    }

    if (f->current_buffer_char + 1 == f->chars_in_buffer) {
        f->chars_in_buffer = 0;
        f->current_buffer_char = -1;
    }
    if (f->chars_in_buffer == CC_SBUFFER_SIZE) {
        fprintf (stderr, "sread_buffer for %s when buffer full\n", f->fname);
        return 0;
    }

  retry:
    n = (int) read (f->desc, (char *) f->buffer + f->chars_in_buffer,
              CC_SBUFFER_SIZE - f->chars_in_buffer);

    if (n == -1) {
        if (errno == EINTR) {
            fprintf (stderr, "sread_buffer interrupted, retrying\n");
            goto retry;
        }
        perror ("read");
        fprintf (stderr, "sread_buffer failed\n");
        return -1;
    }
    if (n == 0) {
        fprintf (stderr, "sread_buffer encountered EOF\n");
        return -1;
    }
    f->pos += n;
    f->chars_in_buffer += n;

    if (f->status == SRW_EMPTY) f->status = SRW_READ;
    
    return 0;
}

static void sinit (CC_SFILE *s)
{
    s->status = 0;
    s->desc = -1;
    s->type = 0;
    s->chars_in_buffer = 0;
    s->current_buffer_char = -1;
    s->bits_in_last_char = 0;
    s->pos = 0;
    s->fname[0] = '\0';
}

int CCutil_sbits (unsigned int x)
{
    int i;
    int ux = x;
    unsigned int b;

    i = 32;
    b = ((unsigned int) 1) << 31;
    while ((ux & b) == 0 && i > 1) {
        b >>= 1;
        i--;
    }
    return i;
}

int CCutil_sdelete_file (const char *fname)
{
    int rval;

    rval = unlink (fname);
    if (rval) {
        perror (fname);
        fprintf (stderr, "unlink: could not delete %s\n", fname);
    }
    return rval;
}

int CCutil_sdelete_file_backup (const char *fname)
{
    int rval;
    char fbuf_O[CC_SFNAME_SIZE + 32];

    sprintf (fbuf_O, "O%s", fname);
    rval = unlink (fbuf_O);

    return rval;
}

static int prepare_write (CC_SFILE *f)
{
    if (!f) return -1;
    if (f->status == SREAD) {
        fprintf (stderr, "%s not open for output\n", f->fname);
        return -1;
    } else if (f->status == SRW_READ) {
        f->chars_in_buffer = 0;
        f->current_buffer_char = -1;
        f->bits_in_last_char = 0;
        f->status = SRW_WRITE;
    } else if (f->status == SRW_EMPTY) {
        f->status = SRW_WRITE;
    } else if (f->status != SWRITE && f->status != SRW_WRITE) {
        fprintf (stderr, "%s has bogus status %d\n", f->fname, f->status);
        return -1;
    }
    
    return 0;
}

static int prepare_read (CC_SFILE *f)
{
    if (!f) return -1;
    if (f->status == SWRITE) {
        fprintf (stderr, "%s not open for input\n", f->fname);
        return -1;
    } else if (f->status == SRW_WRITE) {
        if (CCutil_sflush (f)) return -1;
        f->chars_in_buffer = 0;
        f->current_buffer_char = -1;
        f->bits_in_last_char = 0;
        f->status = SRW_EMPTY;
    } else if (f->status != SREAD && f->status != SRW_READ &&
               f->status != SRW_EMPTY) {
        fprintf (stderr, "%s has bogus status %d\n", f->fname, f->status);
        return -1;
    }
    
    return 0;
}

#ifdef CC_NETREADY

CC_SFILE *CCutil_snet_open (const char *hname, unsigned short p)
{
    struct hostent *h;
    struct sockaddr_in hsock;
    int s;
    CC_SFILE *f = (CC_SFILE *) NULL;

    memset ((void *) &hsock, 0, sizeof (hsock));

    h = gethostbyname (hname);
    if (h == (struct hostent *) NULL) {
        fprintf (stderr, "cannot get host info for %s\n", hname);
        return (CC_SFILE *) NULL;
    }
    memcpy ((void *) &hsock.sin_addr, (void *) h->h_addr, h->h_length);
    hsock.sin_family = AF_INET;
    hsock.sin_port = htons(p);

    s = socket (AF_INET, SOCK_STREAM, 0);
    if (s < 0) {
        perror ("socket");
        fprintf (stderr, "Unable to get socket\n");
        return (CC_SFILE *) NULL;
    }
    if (connect (s, (struct sockaddr *) &hsock, sizeof (hsock)) < 0) {
        perror ("connect");
        fprintf (stderr, "Unable to connect to %s\n", hname);
        return (CC_SFILE *) NULL;
    }

    f = sdopen_readwrite (s);
    if (f == (CC_SFILE *) NULL) {
        fprintf (stderr, "sdopen_readwrite failed\n");
        return (CC_SFILE *) NULL;
    }

    return f;
}

CC_SFILE *CCutil_snet_receive (CC_SPORT *s)
{
    struct sockaddr_in new;
    unsigned int l;  /* Bico 160728 changed from int to unsigned int */
    int t;
    CC_SFILE *f = (CC_SFILE *) NULL;

    memset ((void *) &new, 0, sizeof (new));
    new.sin_family = AF_INET;
    new.sin_addr.s_addr = INADDR_ANY;
    new.sin_port = 0;
    l = sizeof (new);

    t = accept (s->t, (struct sockaddr *) &new, &l);
    if (t < 0) {
        perror ("accept");
        fprintf (stderr, "accept failed\n");
        return (CC_SFILE *) NULL;
    }

    f = sdopen_readwrite (t);
    if (f == (CC_SFILE *) NULL) {
        fprintf (stderr, "sdopen_readwrite failed\n");
        return (CC_SFILE *) NULL;
    }

    return f;
}

CC_SPORT *CCutil_snet_listen (unsigned short p)
{
    int s = -1;
    struct sockaddr_in me;
    CC_SPORT *sp = (CC_SPORT *) NULL;
    
    s = socket (AF_INET, SOCK_STREAM, 0);
    if (s < 0) {
        perror ("socket");
        fprintf (stderr, "Unable to get socket\n");
        goto FAILURE;
    }

    memset ((void *) &me, 0, sizeof (me));

    me.sin_addr.s_addr = INADDR_ANY;
    me.sin_family = AF_INET;
    me.sin_port = htons (p);

    if (bind (s, (struct sockaddr *) &me, sizeof (me)) < 0) {
        perror ("bind");
        fprintf (stderr, "Cannot bind socket\n");
        goto FAILURE;
    }

    if (listen (s, 100) < 0) {
        perror ("listen");
        fprintf (stderr, "Cannot listen to socket\n");
        goto FAILURE;
    }

    sp = CC_SAFE_MALLOC (1, CC_SPORT);
    if (sp == (CC_SPORT *) NULL) {
        fprintf (stderr, "Out of memory in CCutil_snet_listen\n");
        goto FAILURE;
    }

    sp->t = s;
    sp->port = p;

    return sp;

  FAILURE:
    if (s >= 0) close (s);
    CC_IFFREE (sp, CC_SPORT);
    return (CC_SPORT *) NULL;
}

void CCutil_snet_unlisten (CC_SPORT *s)
{
    if (s != (CC_SPORT *) NULL) {
        close (s->t);
        CC_FREE (s, CC_SPORT);
    }
}

#endif /* CC_NETREADY */

/****************************************************************************/
/*                                                                          */
/*                         SORTING ROUTINES                                 */
/*                                                                          */
/*                             TSP CODE                                     */
/*                                                                          */
/*   Written by:  Applegate, Bixby, Chvatal, and Cook                       */
/*   DATE:  February 24, 1994                                               */
/*                                                                          */
/*    EXPORTED FUNCTIONS:                                                   */
/*                                                                          */
/*  void CCutil_int_array_quicksort (int *len, int n)                       */
/*    len - the array to be sorted                                          */
/*    n - the number of elements in len                                     */
/*    Uses quicksort to put len in increasing order.                        */
/*                                                                          */
/*  void CCutil_int_perm_quicksort (int *perm, int *len, int n)             */
/*    perm - must be allocated by the calling routine, it will return       */
/*           a permutation in increasing order of len.                      */
/*    n - the number of elements in perm and len.                           */
/*                                                                          */
/*  void CCutil_rselect (int *arr, int l, int r, int m,                     */
/*      double *coord, CCrandstate *rstate)                                 */
/*    arr - permutation that will be rearranged                             */
/*    l,r - specify the range of arr that we are interested in              */
/*    m - is the index into l,r that is the break point for the perm        */
/*    coord - gives the keys that determine the ordering                    */
/*                                                                          */
/****************************************************************************/

#include "ccutil.h"

static void
    work_int_perm_quicksort (int *perm, int *len, int n),
    select_split (int *arr, int n, double v, int *start, int *end,
           double *coord),
    select_sort (int *arr, int n, double *coord),
    select_sort_dsample (double *samp, int n);

void CCutil_int_array_quicksort (int *len, int n)
{
    int i, j, temp, t;

    if (n <= 1)
        return;

    CC_SWAP (len[0], len[(n - 1)/2], temp);

    i = 0;
    j = n;
    t = len[0];

    while (1) {
        do i++; while (i < n && len[i] < t);
        do j--; while (len[j] > t);
        if (j < i) break;
        CC_SWAP (len[i], len[j], temp);
    }
    CC_SWAP (len[0], len[j], temp);

    CCutil_int_array_quicksort (len, j);
    CCutil_int_array_quicksort (len + i, n - i);
}

void CCutil_int_perm_quicksort (int *perm, int *len, int n)
{
    int i;

    for (i = 0; i < n; i++) perm[i] = i;
    work_int_perm_quicksort (perm, len, n);
}

static void work_int_perm_quicksort (int *perm, int *len, int n)
{
    int i, j, temp, t;

    if (n <= 1)
        return;

    CC_SWAP (perm[0], perm[(n - 1)/2], temp);

    i = 0;
    j = n;
    t = len[perm[0]];

    while (1) {
        do i++; while (i < n && len[perm[i]] < t);
        do j--; while (len[perm[j]] > t);
        if (j < i) break;
        CC_SWAP (perm[i], perm[j], temp);
    }
    CC_SWAP (perm[0], perm[j], temp);

    work_int_perm_quicksort (perm, len, j);
    work_int_perm_quicksort (perm + i, len, n - i);
}

/**********  Median - Select Routines **********/

/* NSAMPLES should be odd */
#define NSAMPLES 3
#define SORTSIZE 20


void CCutil_rselect (int *arr, int l, int r, int m, double *coord,
        CCrandstate *rstate)
{
    double samplevals[NSAMPLES];
    int i;
    int st, en;
    int n;

    arr += l;
    n = r - l + 1;
    m -= l;

    while (n > SORTSIZE) {
        for (i = 0; i < NSAMPLES; i++) {
            samplevals[i] = coord[arr[CCutil_lprand (rstate) % n]];
        }
        select_sort_dsample (samplevals, NSAMPLES);
        select_split (arr, n, samplevals[(NSAMPLES - 1) / 2], &st, &en, coord);
        if (st > m) {
            n = st;
        } else if (en <= m) {
            arr += en;
            n -= en;
            m -= en;
        } else {
            return;
        }
    }

    select_sort (arr, n, coord);
    return;
}

static void select_split (int *arr, int n, double v, int *start, int *end,
                          double *coord)
{
    int i, j, k;
    int t;

    i = 0;
    j = k = n;

    while (i < j) {
        if (coord[arr[i]] < v) {
            i++;
        } else if (coord[arr[i]] == v) {
            j--;
            CC_SWAP (arr[i], arr[j], t);
        } else {
            j--;
            k--;
            t = arr[i];
            arr[i] = arr[j];
            arr[j] = arr[k];
            arr[k] = t;
        }
    }
    *start = j;
    *end = k;
    return;
}

static void select_sort (int *arr, int n, double *coord)
{
    int i, j;
    int t;

    for (i = 1; i < n; i++) {
        t = arr[i];
        for (j = i; j > 0 && coord[arr[j - 1]] > coord[t]; j--) {
            arr[j] = arr[j - 1];
        }
        arr[j] = t;
    }
}

static void select_sort_dsample (double *samp, int n)
{
    int i, j;
    double t;

    for (i = 1; i < n; i++) {
        t = samp[i];
        for (j = i; j > 0 && samp[j - 1] > t; j--) {
            samp[j] = samp[j - 1];
        }
        samp[j] = t;
    }
}

static void iselect_split (int *arr, int n, double v, int *start, int *end,
    int *coord);
static void iselect_sort (int *arr, int n, int *coord);
static void iselect_sort_dsample (int *samp, int n);

void CCutil_iselect (int *arr, int l, int r, int m, int *coord,
        CCrandstate *rstate)
{
    int samplevals[NSAMPLES];
    int i;
    int st, en;
    int n;

    arr += l;
    n = r - l + 1;
    m -= l;

    while (n > SORTSIZE) {
        for (i = 0; i < NSAMPLES; i++) {
            samplevals[i] = coord[arr[CCutil_lprand (rstate) % n]];
        }
        iselect_sort_dsample (samplevals, NSAMPLES);
        iselect_split (arr, n, samplevals[(NSAMPLES - 1) / 2], &st, &en, coord);
        if (st > m) {
            n = st;
        } else if (en <= m) {
            arr += en;
            n -= en;
            m -= en;
        } else {
            return;
        }
    }

    iselect_sort (arr, n, coord);
    return;
}

static void iselect_split (int *arr, int n, double v, int *start, int *end,
                           int *coord)
{
    int i, j, k;
    int t;

    i = 0;
    j = k = n;

    while (i < j) {
        if (coord[arr[i]] < v) {
            i++;
        } else if (coord[arr[i]] == v) {
            j--;
            CC_SWAP (arr[i], arr[j], t);
        } else {
            j--;
            k--;
            t = arr[i];
            arr[i] = arr[j];
            arr[j] = arr[k];
            arr[k] = t;
        }
    }
    *start = j;
    *end = k;
    return;
}

static void iselect_sort (int *arr, int n, int *coord)
{
    int i, j;
    int t;

    for (i = 1; i < n; i++) {
        t = arr[i];
        for (j = i; j > 0 && coord[arr[j - 1]] > coord[t]; j--) {
            arr[j] = arr[j - 1];
        }
        arr[j] = t;
    }
}

static void iselect_sort_dsample (int *samp, int n)
{
    int i, j;
    int t;

    for (i = 1; i < n; i++) {
        t = samp[i];
        for (j = i; j > 0 && samp[j - 1] > t; j--) {
            samp[j] = samp[j - 1];
        }
        samp[j] = t;
    }
}

/****************************************************************************/
/*                                                                          */
/*              MACHINE INDEPENDENT RANDOM NUMBER GENERATOR                 */
/*                                                                          */
/*                             TSP CODE                                     */
/*                                                                          */
/*                                                                          */
/*  Written by:  DIMACS  (modified for TSP)                                 */
/*  Date: February 7, 1995  (cofeb16)                                       */
/*        September 18, 2001  (billenium fix)                               */
/*                                                                          */
/*    EXPORTED FUNCTIONS:                                                   */
/*                                                                          */
/*  void CCutil_sprand (int seed, CCrandstate *r)                           */
/*    - Call once to initialize the generator.                              */
/*                                                                          */
/*  int CCutil_lprand (CCrandstate *r)                                      */
/*    - Returns an integer in the range 0 to CC_PRANDMAX - 1.               */
/*                                                                          */
/*  double CCutil_normrand (CCrandstate *r)                                 */
/*    - Returns a normally-distributed random value with mean 0 and         */
/*      deviation 1.                                                        */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/*    NOTES (from DIMACS):                                                  */
/*        This file contains a set of c-language functions for generating   */
/*    uniform integers.   This is a COMPLETELY PORTABLE generator. It will  */
/*    give IDENTICAL sequences of random numbers for any architecture with  */
/*    at least 30-bit integers, regardless of the integer representation,   */
/*    INT_MAX value, or roundoff/truncation method, etc.                    */
/*        This Truly Remarkable RNG is described more fully in              */
/*    J. Bentley's column, ``The Software Exploratorium ''. It is based on  */
/*    one in Knuth, Vol 2, Section 3.2.2 (Algorithm A).                     */
/*                                                                          */
/*  CCutil_normrand is not from DIMACS or Bentley, but rather just uses     */
/*  the Box-Muller transformation to generate a normally-distributed        */
/*  random variable from two uniform ones.                                  */
/*                                                                          */
/****************************************************************************/


#include "ccutil.h"


void CCutil_sprand (int seed, CCrandstate *r)
{
    int i, ii;
    int last, next;
    int *arr = r->arr;

    seed %= CC_PRANDMAX;
    if (seed < 0) seed += CC_PRANDMAX;

    arr[0] = last = seed;
    next = 1;
    for (i = 1; i < 55; i++) {
        ii = (21 * i) % 55;
        arr[ii] = next;
        next = last - next;
        if (next < 0)
            next += CC_PRANDMAX;
        last = arr[ii];
    }
    r->a = 0;
    r->b = 24;
    for (i = 0; i < 165; i++)
        last = CCutil_lprand (r);
}


int CCutil_lprand (CCrandstate *r)
{
    int t;

    if (r->a-- == 0)
        r->a = 54;
    if (r->b-- == 0)
        r->b = 54;

    t = r->arr[r->a] - r->arr[r->b];

    if (t < 0)
        t += CC_PRANDMAX;

    r->arr[r->a] = t;

    return t;
}


#ifdef      TRY_CODE

/*-----------------------------------------------*/
/* This is a little driver program so you can    */
/* test the code.                                */
/* Typing: a.out 0 3 1                           */
/* should produce                                */
/*     921674862                                 */
/*     250065336                                 */
/*     377506581                                 */
/*  Typing: a.out 1000000 1 2                    */
/*  should produce                               */
/*     57265995                                  */
/*-----------------------------------------------*/

int main (int ac, char **av)
{
    int i;
    int j;
    int n;
    int m;
    int seed;
    CCrandstate rstate;

    if (ac < 4) {
        fprintf (stderr, "Usage: #discard #print #seed\n");
        return 0;
    }
    m = atoi (av[1]);           /* Number to discard initially */
    n = atoi (av[2]);           /* Number to print */
    seed = atoi (av[3]);        /* Seed */

    CCutil_sprand (seed, &rstate);

    for (i = 0; i < m; i++)
        j = CCutil_lprand (&rstate);
    for (i = 0; i < n; i++)
        printf ("%ld\n", CCutil_lprand (&rstate));
    return 0;
}

#endif  /* TRY_CODE */


double CCutil_normrand (CCrandstate *r)
{
    double x1 = ((double) CCutil_lprand(r)) / ((double) CC_PRANDMAX);
    double x2 = ((double) CCutil_lprand(r)) / ((double) CC_PRANDMAX);

    return sqrt (-2*log(x1)) * cos(2*M_PI*x2);

}

/****************************************************************************/
/*                                                                          */
/*               MISCELLANEOUS UTILITY ROUTINES                             */
/*                                                                          */
/*                           TSP CODE                                       */
/*                                                                          */
/*  Written by:  Applegate, Bixby, Chvatal, and Cook                        */
/*  Date: October 12, 1995                                                  */
/*  Date: September 28, 1997                                                */
/*        April 7, 2003 (bico)                                              */
/*                                                                          */
/*    EXPORTED FUNCTIONS:                                                   */
/*                                                                          */
/*  unsigned int CCutil_nextprime (unsigned int x)                          */
/*    FINDS the smallest positive prime >= x                                */
/*                                                                          */
/*  void CCutil_printlabel (void)                                           */
/*    PRINTS information identifying a machine                              */
/*                                                                          */
/*  int CCutil_print_command (int ac, char **av)                            */
/*    PRINTS the command line                                               */
/*                                                                          */
/****************************************************************************/

#include "ccutil.h"


static int
    isprime (unsigned int x);


unsigned int CCutil_nextprime (unsigned int x)
{
    if (x < 3) return 3;
    x |= 1;
    while (!isprime (x)) x += 2;
    return x;
}

static int isprime (unsigned int p)
{
    unsigned int i;

    if ((p&1) == 0) return 0;
    for (i=3; i*i<=p; i+=2) {
        if (p%i == 0) return 0;
    }
    return 1;
}

void CCutil_printlabel (void)
{
#ifdef CC_NETREADY
    char buf[1024];

    gethostname (buf, 1024);
    printf ("Host: %s  Current process id: %d\n", buf, (int) getpid());
    fflush (stdout);
#else
    printf ("No label - need function to non-NETREADY machines\n");
    fflush (stdout);
#endif
}

int CCutil_print_command (int ac, char **av)
{
    int rval = 0;
    int i, cmdlen = 0;
    char *cmdout = (char *) NULL;

    for (i=0; i<ac; i++) {
        cmdlen += strlen(av[i]) + 1;
        if (strlen(av[i]) == 1) cmdlen++;
    }
    cmdout = CC_SAFE_MALLOC (cmdlen, char);
    CCcheck_NULL (cmdout, "out of memory in print_command");

    cmdlen = 0;
    for (i=0; i<ac; i++) {
        if (strlen(av[i]) == 1) cmdout[cmdlen++] = '-';
        strcpy (cmdout + cmdlen, av[i]);
        cmdlen += strlen(av[i]);
        cmdout[cmdlen++] = ' ';
    }
    cmdout[cmdlen-1] = '\0';
    printf ("%s\n", cmdout); fflush (stdout);

CLEANUP:

    CC_IFFREE (cmdout, char);
    return rval;
}

char *CCutil_strchr (char *s, int c)
{
    while (*s) {
        if (*s == c) return s;
        s++;
    }
    return (char *) NULL;
}

const char *CCutil_strrchr_c (const char *s, int c)
{
    const char *l = (char *) NULL;

    while (*s) {
        if (*s == c) l = s;
        s++;
    }
    return l;
}

char *CCutil_strdup (const char *s)
{
    char *p = CC_SAFE_MALLOC (strlen(s)+1, char);

    if (p == (char *) NULL) {
        fprintf (stderr, "Out of memory in CCutil_strdup\n");
        return (char *) NULL;
    }
    strcpy (p, s);
    return p;
}


/****************************************************************************/
/*                                                                          */
/*                        TIMING FUNCTIONS                                  */
/*                                                                          */
/*                            TSP CODE                                      */
/*                                                                          */
/*                                                                          */
/*  Written by:  Applegate, Bixby, Chvatal, and Cook                        */
/*  Date: Summer 1994  (cofeb16)                                            */
/*        December 1997 (dla)                                               */
/*                                                                          */
/*                                                                          */
/*    EXPORTED FUNCTIONS:                                                   */
/*                                                                          */
/*  double CCutil_zeit (void)                                               */
/*        - To measure cpu time.                                            */
/*    To use this, set double t = CCutil_zeit (), run the function you      */
/*    want to time, then compute CCutil_zeit () - t.                        */
/*                                                                          */
/*  double CCutil_real_zeit (void)                                          */
/*    - To measure wall clock time.                                         */
/*                                                                          */
/*    To use this, set double t = CCutil_real_zeit (), run the function     */
/*    you want to time, then compute CCutil_real_zeit () - t.               */
/*                                                                          */
/****************************************************************************/

#include "ccutil.h"

#ifdef HAVE_GETRUSAGE

#ifdef HAVE_SYS_RESOURCE_H
# include <sys/resource.h>
#endif

double CCutil_zeit (void)
{
    struct rusage ru;

    getrusage (RUSAGE_SELF, &ru);

    return ((double) ru.ru_utime.tv_sec) +
            ((double) ru.ru_utime.tv_usec) / 1000000.0;
}
#else /* HAVE_GETRUSAGE */

#ifdef HAVE_TIMES

#ifdef HAVE_SYS_PARAM_H
# include <sys/param.h>
#endif
#ifdef HAVE_SYS_TIMES_H
# include <sys/times.h>
#endif

#ifdef CLK_TCK
#define MACHINE_FREQ CLK_TCK
#else
#define MACHINE_FREQ HZ
#endif

double CCutil_zeit (void)
{
    struct tms now;

    times (&now);
    return ((double) now.tms_utime) / ((double) MACHINE_FREQ);
}
#else /* HAVE_TIMES */

#ifdef HAVE_CLOCK

#ifndef CLOCKS_PER_SEC
#ifdef CLK_TCK
#define CLOCKS_PER_SEC CLK_TCK
#else
#define CLOCKS_PER_SEC 60
#endif
#endif

double CCutil_zeit (void)
{
    return ((double) clock()) / ((double) CLOCKS_PER_SEC);
}

#else /* HAVE_CLOCK */

double CCutil_zeit (void)
{
    return 0.0;
}
#endif /* HAVE_CLOCK */
#endif /* HAVE_TIMES */
#endif /* HAVE_GETRUSAGE */

double CCutil_real_zeit (void)
{
    return (double) time (0);
}

char *CCutil_problabel (const char *probloc)
{
    const char *p;
    const char *problabel = probloc;
    char *probcopy = (char *) NULL;
    char *p2;

    p = CCutil_strrchr_c (problabel, ':');
    if (p != (const char *) NULL) problabel = p+1;
    p = CCutil_strrchr_c (problabel, '/');
    if (p != (const char *) NULL) problabel = p+1;
    probcopy = CCutil_strdup (problabel);
    if (probcopy == (char *) NULL) return (char *) NULL;
    p2 = CCutil_strchr (probcopy, '.');
    if (p2 != (char *) NULL) *p2 = '\0';
    return probcopy;
}
