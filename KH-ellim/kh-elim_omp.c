/*                  
MIT License

Copyright (c) 2022  Keld Helsgaun

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

#include <omp.h>
#include "kh-elim.h"

typedef struct node {
    int id, degree, *cost, *deg;
    double x, y, z;
    edge *first_edge;
    struct node *tour_suc, **adj[2];
    int *level;
    omp_lock_t lock;
} node;

#include "kh-elim_common.c"

int max_threads;
node ***cd_set;
int *cd_set_size;
int cache_mask, **cache_val, **cache_sig;

int main(int argc, char *argv[])
{
    int old_eliminated, i, j, k, m;
    double start_time;

    max_threads = omp_get_max_threads();
    // max_threads = 1; omp_set_num_threads(max_threads);
    if (parseargs(argc, argv))
        exit(1);
    read_edge_set_file(edge_file_name);
    read_tsp_file(tsp_file_name);
    read_tour_file(tour_file_name);
    for (i = 0; i < n; i++) {
        node_set[i].cost = (int *) malloc(max_threads * sizeof(int));
        node_set[i].deg = (int *) malloc(max_threads * sizeof(int));
        node_set[i].level = (int *) malloc(max_threads * sizeof(int));
        node_set[i].adj[0] =
            (node **) malloc(max_threads * sizeof(node *));
        node_set[i].adj[1] =
            (node **) malloc(max_threads * sizeof(node *));
        omp_init_lock(&node_set[i].lock);
    }
    if (C == C_cache) {
        for (i = 0; (1 << i) < (n << 1); i++);
        i = 1 << i;
        cache_mask = i - 1;
        cache_val = (int **) malloc(max_threads * sizeof(int *));
        cache_sig = (int **) malloc(max_threads * sizeof(int *));
        for (j = 0; j < max_threads; j++) {
            cache_val[j] = (int *) malloc(i * sizeof(int));
            cache_sig[j] = (int *) malloc(i * sizeof(int));
            for (k = 0; k < i; k++)
                cache_sig[j][k] = -1;
        }
    }
    printf("Edges:     %d\n", edges);
    edge_set = (edge **) malloc((n - 1) * sizeof(edge *));
    sort_edges();
    cd_set = (node ***) malloc(max_threads * sizeof(node **));
    for (i = 0; i < max_threads; i++)
        cd_set[i] = (node **) malloc((max_cd_set_size + 1) * sizeof(node *));
    cd_set_size = (int *) malloc(max_threads * sizeof(int));
    for (m = 0; m < MAX_N; m++)
        for (k = 0; k <= m; k++)
            if (k == 0 || k == m)
                binomial[m][k] = 1;
            else
                binomial[m][k] = binomial[m - 1][k - 1] +
                                 binomial[m - 1][k];
    start_time = get_time();

    if (Jonker_Volgenant) {
        JV_elimination();
        printf("|");
        fflush(stdout);
    }
    do {
        old_eliminated = eliminated;
        HS_elimination();
        if (Jonker_Volgenant)
            JV_elimination();
        printf("|");
        fflush(stdout);
    } while (repeated && eliminated > old_eliminated);

    printf("\nEliminated: %d\n", eliminated);
    printf("Remaining: %d", edges - eliminated);
    write_output_file(output_file_name);
    printf("\nTime:      %.2f seconds\n", get_time() - start_time);
    return 0;
}

int eliminate_edge(node * a, node * b)
{
    edge *e = 0, *prev;
    int i;

    if (a->degree <= 2 || b->degree <= 2)
        return 0;
    if (a->id < b->id) {
        omp_set_lock(&a->lock);
        omp_set_lock(&b->lock);
    } else {
        omp_set_lock(&b->lock);
        omp_set_lock(&a->lock);
    }
    for (i = 1; i <= 2; i++, swap(&a, &b)) {
        prev = 0;
        for (e = a->first_edge; e && e->to != b; prev = e, e = e->next);
        if (!e)
            break;
        if (prev)
            prev->next = e->next;
        else
            a->first_edge = e->next;
        a->degree--;
    }
    omp_unset_lock(&b->lock);
    omp_unset_lock(&a->lock);
    if (print_E && e != NULL) {
        printf("E");
        fflush(stdout);
    }
    return e != NULL;
}

void JV_elimination()
{
    node *a, *b, *c, *d;
    edge *ea, *ec;
    int Cac, Cbc, i, j, tid;
    
    processed = 0;
#pragma omp parallel for schedule(dynamic) \
        private(a, b, c, d, ea, ec, tid, Cac, Cbc, i, j)
    for (i = 0; i < n; i++) {
        a = &node_set[i];
        tid = omp_get_thread_num();
        for (ea = a->first_edge; ea; ea = ea->next) {
            b = ea->to;
            if (a->id > b->id)
                continue;
            if (a->degree <= 2 || b->degree <= 2 || is_tour_edge(a, b))
                goto Next_b;
            potential_cd_nodes(a, b);
            for (j = 0; j < cd_set_size[tid]; j++) {
                c = cd_set[tid][j];
                Cac = C(a, c);
                Cbc = C(b, c);
                for (ec = c->first_edge; ec; ec = ec->next) {
                    d = ec->to;
                    if (d != a && d != b &&
                        (ea->cost + ec->cost <= Cac + C(d, b) ||
                         ea->cost + ec->cost <= C(a, d) + Cbc))
                        break;
                }
                if (!ec)
                    break;
            }
            if (j < cd_set_size[tid] && eliminate_edge(a, b))
#pragma omp atomic
                eliminated++;
          Next_b:
#pragma omp atomic
            processed++;
            if (processed % 1000 == 999) {
                printf(".");
                fflush(stdout);
            }
        }
    }
}

void HS_elimination()
{
    node *a, *b, *c, *d;
    edge *ea;
    int Cab, Ccd, cd_count, i, j, k, tid;

    processed = 0;
#pragma omp parallel for schedule(dynamic) \
        private(a, b, c, d, ea, Cab, Ccd, cd_count, i, j, k, tid)
    for (i = 0; i < n; i++) {
        a = &node_set[i];
        tid = omp_get_thread_num();
        for (ea = a->first_edge; ea; ea = ea->next) {
            b = ea->to;
            if (a->id > b->id)
                continue;
            triangle_elimination(a, b);
            if (a->degree <= 2 || b->degree <= 2 || is_tour_edge(a, b))
                goto Next_b;
            Cab = ea->cost;
            if (diamond_elimination(a, b, Cab))
                goto Next_b;
            potential_cd_nodes(a, b);
            cd_count = 0;
            for (j = 1; j < cd_set_size[tid]; j++) {
                c = cd_set[tid][j];
                for (k = 0; k < j; k++) {
                    if (++cd_count > max_cd_count)
                        goto Next_b;
                    d = cd_set[tid][k];
                    if (compatible(a, b, c, d, Cab, Ccd = C(c, d)))
                        continue;
                    if (can_HS_eliminate(a, b, c, d, Cab)) {
                        if (eliminate_edge(a, b))
#pragma omp atomic
                            eliminated++;
                        goto Next_b;
                    }
                }
            }
          Next_b:
#pragma omp atomic
            processed++;
            if (processed % 1000 == 999) {
                printf(".");
                fflush(stdout);
            }
        }
    }
}

void triangle_elimination(node * a, node * b)
{
    node *c, *d;
    edge *ea, *eb;

    if (a->degree == 2) {
        for (ea = a->first_edge; ea; ea = ea->next) {
            c = ea->to;
            if (c != b && is_edge(b, c)) {
                eliminate_edge(b, c);
#pragma omp atomic
                eliminated++;
            }
        }
    }
    if (b->degree == 2) {
        for (eb = b->first_edge; eb; eb = eb->next) {
            d = eb->to;
            if (d != a && is_edge(a, d)) {
                eliminate_edge(a, d);
#pragma omp atomic
                eliminated++;
            }
        }
    }
}

int diamond_elimination(node * a, node * b, int Cab)
{
    node *c, *d;
    edge *ea, *eb;
    int diff;
    
    if (a->degree != 3 || b->degree != 3)
        return 0;
    for (ea = a->first_edge; ea; ea = ea->next) {
        c = ea->to;
        if (c != b && c->degree == 3 && is_edge(c, b)) {
            for (eb = b->first_edge; eb; eb = eb->next) {
                d = eb->to;
                if (d != c && d != a && d->degree == 3 && is_edge(d, a)) {
                    diff = ea->cost + eb->cost - C(c, b) - C(d, a);
                    if (diff > 0) {
                        eliminate_edge(c, a);
                        eliminate_edge(d, b);
                    } else if (diff < 0) {
                        eliminate_edge(c, b);
                        eliminate_edge(d, a);
                    }
#pragma omp atomic
                    eliminated += 2;
                    return 1;
                }
            }
        }
    }
    return 0;
}

void potential_cd_nodes(node * a, node * b)
{
    node *c, *d;
    edge *ea, *ec;
    int i, j, cost, tid = omp_get_thread_num();

    for (i = 1; i <= 2; i++, swap(&a, &b)) {
        for (ea = a->first_edge; ea; ea = ea->next) {
            c = ea->to;
            c->cost[tid] = INT_MAX;
            if (quicker_cd_search)
                continue;
            for (ec = c->first_edge; ec; ec = ec->next)
                ec->to->cost[tid] = INT_MAX;
        }
    }
    a->cost[tid] = b->cost[tid] = 0;
    cd_set_size[tid] = 0;
    for (i = 1; i <= 2; i++, swap(&a, &b)) {
        for (ea = a->first_edge; ea; ea = ea->next) {
            c = ea->to;
            if (c->cost[tid] == INT_MAX) {
                cost = c->cost[tid] = ea->cost + C(c, b);
                if (cd_set_size[tid] == max_cd_set_size &&
                    cost >= cd_set[tid][max_cd_set_size - 1]->cost[tid])
                    continue;
                j = cd_set_size[tid];
                for (; j > 0 && cost < cd_set[tid][j - 1]->cost[tid]; j--)
                    cd_set[tid][j] = cd_set[tid][j - 1];
                cd_set[tid][j] = c;
                if (cd_set_size[tid] < max_cd_set_size)
                    cd_set_size[tid]++;
            }
        }
    }
    if (quicker_cd_search)
        return;
    for (i = 1; i <= 2; i++, swap(&a, &b)) {
        for (ea = a->first_edge; ea; ea = ea->next) {
            c = ea->to;
            for (ec = c->first_edge; ec; ec = ec->next) {
                d = ec->to;
                if (d->cost[tid] == INT_MAX) {
                    cost = d->cost[tid] = C(a, d) + C(d, b);
                    if (cd_set_size[tid] == max_cd_set_size &&
                        cost >= cd_set[tid][ max_cd_set_size - 1]->cost[tid])
                        continue;
                    j = cd_set_size[tid];
                    for (; j > 0 && cost < cd_set[tid][j - 1]->cost[tid];
                         j--)
                        cd_set[tid][j] = cd_set[tid][j - 1];
                    cd_set[tid][j] = d;
                    if (cd_set_size[tid] < max_cd_set_size)
                        cd_set_size[tid]++;
                }
            }
        }
    }
}

int extra_node_opt(node * a, node * b,
                   node * c1, node * c, node * c2,
                   node * d1, node * d, node * d2,
                   node * e1, node * e, node * e2,
                   node * f1, node * f, node * f2,
                   int Cab,
                   int Cc1c, int Ccc2,
                   int Cd1d, int Cdd2,
                   int Ce1e, int Cee2,
                   int Cf1f, int Cff2)
{
    node *x, *x1, *x2;
    int Cx1x, Cxx2, i, tid;
    edge *ex1, *ex2;

    if (extra_nodes < 1 ||
        (e && extra_nodes == 1) ||
        (f && extra_nodes == 2))
        return 1;
    tid = omp_get_thread_num();
    for (i = 0; i < cd_set_size[tid]; i++) {
        x = cd_set[tid][i];
        if (x == a || x == b ||
            x == c1 || x == c || x == c2 ||
            x == d1 || x == d || x == d2 ||
            x == e1 || x == e || x == e2 ||
            x == f1 || x == f || x == f2)
            continue;
        for (ex1 = x->first_edge; ex1; ex1 = ex1->next) {
            x1 = ex1->to;
            if (x1 == c || x1 == d || x1 == e || x1 == f ||
                ((x1 == a || x1 == b) &&
                 (x1 == c1 || x1 == c2 || x1 == d1 || x1 == d2 ||
                  x1 == e1 || x1 == e2 || x1 == f1 || x1 == f2)) ||
                ((x1 == c1 || x1 == c2) &&
                 (x1 == a || x1 == b || x1 == d1 || x1 == d2 ||
                  x1 == e1 || x1 == e2 || x1 == f1 || x1 == f2)) ||
                ((x1 == d1 || x1 == d2) &&
                 (x1 == a || x1 == b || x1 == c1 || x1 == c2 ||
                  x1 == e1 || x1 == e2 || x1 == f1 || x1 == f2)) ||
                ((x1 == e1 || x1 == e2) &&
                 (x1 == a || x1 == b || x1 == c1 || x1 == c2 ||
                  x1 == d1 || x1 == d2 || x1 == f1 || x1 == f2)) ||
                ((x1 == f1 || x1 == f2) &&
                 (x1 == a || x1 == b || x1 == c1 || x1 == c2 ||
                  x1 == d1 || x1 == d2 || x1 == e1 || x1 == e2)))
                continue;
            Cx1x = ex1->cost;
            if (!opt22(x1, x, a, b, Cx1x, Cab) ||
                !opt22(x1, x, c1, c, Cx1x, Cc1c) ||
                !opt22(x1, x, c, c2, Cx1x, Ccc2) ||
                !opt22(x1, x, d1, d, Cx1x, Cd1d) ||
                !opt22(x1, x, d, d2, Cx1x, Cdd2) ||
                !opt23(x1, x, c1, c, c2, Cx1x, Cc1c, Ccc2) ||
                !opt23(x1, x, d1, d, d2, Cx1x, Cd1d, Cdd2) ||
                !opt233(x1, x, c1, c, c2, d1, d, d2) ||
                !opt2233(x1, x, a, b, c1, c, c2, d1, d, d2) ||
                (e && (!opt22(x1, x, e1, e, Cx1x, Ce1e) ||
                       !opt22(x1, x, e, e2, Cx1x, Cee2) ||
                       !opt23(x1, x, e1, e, e2, Cx1x, Ce1e, Cee2) ||
                       (!f && !opt22333(x1, x, a, b, c1, c, c2,
                                        d1, d, d2, e1, e, e2)))) ||
                (f && (!opt22(x1, x, f1, f, Cx1x, Cf1f) ||
                       !opt22(x1, x, f, f2, Cx1x, Cff2) ||
                       !opt23(x1, x, f1, f, f2, Cx1x, Cf1f, Cff2) ||
                       !opt223333(x1, x, a, b, c1, c, c2, d1, d, d2,
                                  e1, e, e2, f1, f, f2))))
                continue;
            for (ex2 = ex1->next; ex2; ex2 = ex2->next) {
                x2 = ex2->to;
                if (x2 == c || x2 == d || x2 == e || x2 == f ||
                    ((x2 == a || x2 == b) &&
                     (x1 == a || x1 == b ||
                      x2 == c1 || x2 == c2 || x2 == d1 || x2 == d2 ||
                      x2 == e1 || x2 == e2 || x2 == f1 || x2 == f2)) ||
                    ((x2 == c1 || x2 == c2) &&
                     (x1 == c1 || x1 == c2 ||
                      x2 == a || x2 == b || x2 == d1 || x2 == d2 ||
                      x2 == e1 || x2 == e2 || x2 == f1 || x2 == f2)) ||
                    ((x2 == d1 || x2 == d2) &&
                     (x1 == d1 || x1 == d2 ||
                      x2 == a || x2 == b || x2 == c1 || x2 == c2 ||
                      x2 == e1 || x2 == e2 || x2 == f1 || x2 == f2)) ||
                    ((x2 == e1 || x2 == e2) &&
                     (x1 == e1 || x1 == e2 ||
                      x2 == a || x2 == b || x2 == c1 || x2 == c2 ||
                      x2 == d1 || x2 == d2 || x2 == f1 || x2 == f2)) ||
                    ((x2 == f1 || x2 == f2) &&
                     (x1 == f1 || x1 == f2 ||
                      x2 == a ||  x2 == b || x2 == c1 || x2 == c2 ||
                      x2 == d1 || x2 == d2 || x2 == e1 || x2 == e2)) ||
                    (f ? has_cycle222222(a, b, c1, c2, d1, d2,
                                         e1, e2, f1, f2, x1, x2)
                       : e ? has_cycle22222(a, b, c1, c2, d1, d2,
                                            e1, e2, x1, x2)
                           : has_cycle2222(a, b, c1, c2, d1, d2, x1, x2)))
                    continue;
                Cxx2 = C(x, x2);
                if (opt22(x, x2, a, b, Cxx2, Cab) &&
                    opt22(x, x2, c1, c, Cxx2, Cc1c) &&
                    opt22(x, x2, c, c2, Cxx2, Ccc2) &&
                    opt22(x, x2, d1, d, Cxx2, Cd1d) &&
                    opt22(x, x2, d, d2, Cxx2, Cdd2) &&
                    opt23(a, b, x1, x, x2, Cab, Cx1x, Cxx2) &&
                    opt23(c1, c, x1, x, x2, Cc1c, Cx1x, Cxx2) &&
                    opt23(c, c2, x1, x, x2, Ccc2, Cx1x, Cxx2) &&
                    opt23(d1, d, x1, x, x2, Cd1d, Cx1x, Cxx2) &&
                    opt23(d, d2, x1, x, x2, Cdd2, Cx1x, Cxx2) &&
                    opt23(x, x2, c1, c, c2, Cxx2, Cc1c, Ccc2) &&
                    opt23(x, x2, d1, d, d2, Cxx2, Cd1d, Cdd2) &&
                    opt33(x1, x, x2, c1, c, c2) &&
                    opt33(x1, x, x2, d1, d, d2) &&
                    opt2333(a, b, c1, c, c2, d1, d, d2, x1, x, x2) &&
                    (!e || (opt22(x, x2, e1, e, Cxx2, Ce1e) &&
                            opt22(x, x2, e, e2, Cxx2, Cee2) &&
                            opt23(x, x2, e1, e, e2, Cxx2, Ce1e, Cee2)  &&
                            opt23(e1, e, x1, x, x2, Ce1e, Cx1x, Cxx2) &&
                            opt23(e, e2, x1, x, x2, Cee2, Cx1x, Cxx2) &&
                            opt33(x1, x, x2, e1, e, e2))) &&
                    (!f || (opt22(x, x2, f1, f, Cxx2, Cf1f) &&
                            opt22(x, x2, f, f2, Cxx2, Cff2) &&
                            opt23(x, x2, f1, f, f2, Cxx2, Cf1f, Cff2) &&
                            opt23(f1, f, x1, x, x2, Cf1f, Cx1x, Cxx2) &&
                            opt23(f, f2, x1, x, x2, Cff2, Cx1x, Cxx2) &&
                            opt33(x1, x, x2, f1, f, f2))) &&
                    (!e ||
                     (opt2333(a, b, c1, c, c2, e1, e, e2, x1, x, x2) &&
                      opt2333(a, b, d1, d, d2, e1, e, e2, x1, x, x2) &&
                      opt3333(c1, c, c2, d1, d, d2, e1, e, e2, x1, x, x2) &&
                      opt23333(a, b, c1, c, c2, d1, d, d2, e1, e, e2,
                               x1, x, x2))) &&
                    (!f ||
                     (opt2333(a, b, c1, c, c2, f1, f, f2, x1, x, x2) &&
                      opt2333(a, b, d1, d, d2, f1, f, f2, x1, x, x2) &&
                      opt2333(a, b, e1, e, e2, f1, f, f2, x1, x, x2) &&
                      opt3333(c1, c, c2, d1, d, d2, f1, f, f2, x1, x, x2) &&
                      opt3333(c1, c, c2, e1, e, e2, f1, f, f2, x1, x, x2) &&
                      opt3333(d1, d, d2, e1, e, e2, f1, f, f2, x1, x, x2) &&
                      opt233333(a, b, c1, c, c2, d1, d, d2, e1, e, e2,
                                f1, f, f2, x1, x, x2))) &&
                    extra_edge_opt(a, b, c1, c, c2, d1, d, d2,
                                   x1, x, x2, e1, e, e2) &&
                    (!f || extra_edge_opt(a, b, c1, c, c2, d1, d, d2,
                                          x1, x, x2, f1, f, f2)) &&
                    (f ||
                     (!e ? extra_node_opt(a, b, c1, c, c2, d1, d, d2,
                                          x1, x, x2, 0, 0, 0,
                                          Cab, Cc1c, Ccc2, Cd1d, Cdd2,
                                          Cx1x, Cxx2, 0, 0) :
                           extra_node_opt(a, b, c1, c, c2, d1, d, d2,
                                          e1, e, e2, x1, x, x2,
                                          Cab, Cc1c, Ccc2, Cd1d, Cdd2,
                                          Ce1e, Cee2, Cx1x, Cxx2))))
                    goto Next_x;
            }
        }
        return 0;
    Next_x:;
    }
    return 1;
}

int has_cycle(int count, node ** p)
{
    int i, tid = omp_get_thread_num();
    node *n1, *n2, *start, *prev, *next;

    for (i = 0; i < count; i++)
        p[i]->deg[tid] = 0;
    for (i = 0; i < count; i += 2) {
        n1 = p[i];
        n2 = p[i + 1];
        if (n1->deg[tid] > 1 || n2->deg[tid] > 1)
            return 1;
        n1->adj[n1->deg[tid]++][tid] = n2;
        n2->adj[n2->deg[tid]++][tid] = n1;
    }
    for (i = 0; i < count; i++)
        if (p[i]->deg[tid] != 2)
            break;
    if (i == count)
        return 1;
    for (i = 0; i < count; i++) {
        if (p[i]->deg[tid] != 1)
            continue;
        start = prev = p[i];
        next = start->adj[0][tid];
        start->deg[tid] = 0;
        while (next != start && next->deg[tid] == 2) {
            next->deg[tid] = 0;
            if (next->adj[0][tid] == prev) {
                prev = next;
                next = next->adj[1][tid];
            } else {
                prev = next;
                next = next->adj[0][tid];
            }
        }
        if (next == start)
            return 1;
        next->deg[tid] = 0;
    }
    for (i = 0; i < count; i++)
        if (p[i]->deg[tid] == 2)
            return 1;
    return 0;
}

int C_cache(node * a, node * b)
{
    int i = a->id, j = b->id, k, tid;

    if (i == j)
        return 0;
    if (i > j) {
        k = i;
        i = j;
        j = k;
    }
    k = ((i << 8) + i + j) & cache_mask;
    tid = omp_get_thread_num();
    if (cache_sig[tid][k] == i)
        return cache_val[tid][k];
    cache_sig[tid][k] = i;
    return (cache_val[tid][k] = dist(a, b));
}
