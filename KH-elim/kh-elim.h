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

#ifndef _KH_ELIM_H
#define _KH_ELIM_H

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include <unistd.h>

double ab_stretch = 1.0;
int extra_edges = 0;
int extra_nodes = 0;
int Jonker_Volgenant = 0;
int max_cd_set_size = 10;
int max_cd_count = 10;
int print_E = 0;
int repeated = 0;
int quicker_cd_search = 0;
int strong_3_opt = 0;

char *tsp_file_name; 
char *edge_file_name; 
char *output_file_name; 
char *tour_file_name;

#define MAX_N 20
typedef int aSet[MAX_N];
int binomial[MAX_N][MAX_N];

int concorde_edge_format = 0;

struct node;

typedef struct edge {
    struct node *to;
    int cost;
    struct edge *next;
} edge;

typedef int (*dist_function) (struct node * a, struct node * b);

void read_edge_set_file(char * edge_file_name);
void read_tsp_file(char * edge_file_name);
void read_tour_file(char * edge_file_name);
void write_output_file(char * output_file_name);
void add_edge(struct node * a, struct node * b, int cost);
int eliminate_edge(struct node * a, struct node * b);
int is_edge(struct node * a, struct node * b);
int is_tour_edge(struct node * a, struct node * b);
int compare(const void * ea, const void * eb);
void sort_edges();
void JV_elimination();
void HS_elimination();
void triangle_elimination(struct node * a, struct node * b);
int diamond_elimination(struct node * a, struct node * b, int Cab);
int can_HS_eliminate(struct node * a, struct node * b,
                     struct node * c, struct node * d,
                     int Cab);
int compatible(struct node * a, struct node * b,
               struct node * c, struct node * d,
               int Cab,
               int Ccd);
void swap(struct node ** a, struct node ** b);
void swapC(int * a, int * b);
int fixed(struct node * a, struct node * b);
int min(int a, int b);
void potential_cd_nodes(struct node * a, struct node * b);
int opt22(struct node * a, struct node * b,
          struct node * c, struct node * d,
          int Cab,
          int Ccd);
int opt23_simple(struct node * a, struct node * b,
                 struct node * c1, struct node * c, struct node * c2,
                 int Cab,
                 int Cc1c, int Ccc2);
int opt23(struct node * a, struct node * b,
          struct node * c1, struct node * c, struct node * c2,
          int Cab,
          int Cc1c, int Ccc2);
int opt_excess(struct node * z, 
               struct node * s1, struct node * s2,
               struct node * s3, struct node * s4,
               int Cs1s2, int Cs2s3, int Cs3s4);
void reverse_path(struct node **p, int size);
int is_opt(int paths, struct node *** path, int * size);
int reverse_opt(int paths, struct node ***path, int * size, int n);
int perm_opt(int paths, struct node *** path, int * size, int n);
int opt(int paths, struct node *** path, int * size);
int extra_edge_opt(struct node * a, struct node * b,
                   struct node * c1, struct node * c, struct node * c2,
                   struct node * d1, struct node * d, struct node * d2,
                   struct node * e1, struct node * e, struct node * e2,
                   struct node * f1, struct node * f, struct node * f2);
int extra_extra_edge_24333_opt(struct node * a, struct node * b,
                               struct node * c1, struct node * c2,
                               struct node * c3, struct node * c4,
                               struct node * d1, struct node * d2,
                               struct node * d3,
                               struct node * e1, struct node * e2,
                               struct node * e3,
                               struct node * f1, struct node * f2,
                               struct node * f3,
                               int state); 
int extra_extra_edge_33333_opt(struct node * a1, struct node * a2,
                               struct node * a3,
                               struct node * c1, struct node * c2,
                               struct node * c3,
                               struct node * d1, struct node * d2,
                               struct node * d3,
                               struct node * e1, struct node * e2,
                               struct node * e3,
                               struct node * f1, struct node * f2,
                               struct node * f3,
                               int state);
int extra_node_opt(struct node * a, struct node * b,
                   struct node * c1, struct node * c, struct node * c2,
                   struct node * d1, struct node * d, struct node * d2,
                   struct node * e1, struct node * e, struct node * e2,
                   struct node * f1, struct node * f, struct node * f2,
                   int Cab,
                   int Cc1c, int Ccc2,
                   int Cd1d, int Cdd2,
                   int Ce1e, int Cee2,
                   int Cf1f, int Cff2);
int path_opt(struct node ** perm, int n, int *fix);
int has_cycle(int count, struct node ** p);
int has_cycle222(struct node * a, struct node * b,
                 struct node * c1, struct node * c2,
                 struct node * d1, struct node * d2);
int has_cycle2222(struct node * a, struct node * b,
                  struct node * c1, struct node * c2,
                  struct node * d1, struct node * d2,
                  struct node * e1, struct node * e2);
int has_cycle22222(struct node * a, struct node * b,
                   struct node * c1, struct node * c2,
                   struct node * d1, struct node * d2,
                   struct node * e1, struct node * e2,
                   struct node * f1, struct node * f2);
int has_cycle222222(struct node * a, struct node * b,
                    struct node * c1, struct node * c2,
                    struct node * d1, struct node * d2,
                    struct node * e1, struct node * e2,
                    struct node * f1, struct node * f2,
                    struct node * g1, struct node * g2);
double get_time();

int euc2d(struct node * a, struct node * b);
int euc3d(struct node * a, struct node * b);
int ceil2d(struct node * a, struct node * b);
int geom(struct node * a, struct node * b);
int C_cache(struct node * a, struct node * b);

int parseargs(int argc, char *argv[]);
void usage(char *f);

dist_function dist, C;
int n, edges, eliminated = 0, processed = 0;
struct node *node_set;
edge **edge_set;

#endif
