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

#ifndef _KH_ELIM_COMMON_C
#define _KH_ELIM_COMMON_C

void read_edge_set_file(char * edge_file_name)
{
    FILE *in;
    int from, to, cost = INT_MIN, count, i;
    char line[81];

    in = fopen(edge_file_name, "r");
    if (!in) {
        fprintf(stderr, "Cannot open edge_file: %s\n", edge_file_name);
        exit(1);
    }
    fscanf(in, "%d %d\n", &n, &edges);
    node_set = (node *) calloc(n, sizeof(node));
    for (i = 0; i < n; i++)
        node_set[i].id = i;
    for (i = 0; i < edges; i++) {
        fgets(line, 80, in);
        count = sscanf(line, "%d %d %d\n", &from, &to, &cost);
        if (i == 0)
            concorde_edge_format = count == 3;
         if (concorde_edge_format ? count != 3 : count != 2) {
            printf("Wrong edge_file format: %s\n", line);
            exit(1);
        }
        if (from < 0 || from >= n || to < 0 || to >= n) {
            fprintf(stderr, "Edge (%d, %d) out of range in edge_file: %s\n",
                    from, to, edge_file_name);
            exit(1);
        }
        add_edge(&node_set[from], &node_set[to], cost);
    }
    fclose(in);
}

void read_tsp_file(char * tsp_file_name)
{
    FILE *in;
    int dimension, id, i;
    char line[256], edge_weight_type[80], name[256], *cp;
    double x, y, z;
    edge *e;

    in = fopen(tsp_file_name, "r");
    if (!in) {
        fprintf(stderr, "Cannot open tsp_file: %s\n", tsp_file_name);
        exit(1);
    }
    while (fgets(line, 256, in)) {
        if (!strncmp(line, "EDGE_WEIGHT_TYPE", strlen("EDGE_WEIGHT_TYPE"))) {
            cp = strchr(line, ':');
            sscanf(cp + 1, "%s", edge_weight_type);
            if (!strcmp(edge_weight_type, "EUC_2D"))
                C = dist = euc2d;
            else if (!strcmp(edge_weight_type, "EUC_3D"))
                C = dist = euc3d;
            else if (!strcmp(edge_weight_type, "CEIL_2D"))
                C = dist = ceil2d;
            else if (!strcmp(edge_weight_type, "GEOM")) {
                C = C_cache;
                dist = geom;
            } else {
                fprintf(stderr, "Unknown EDGE_WEIGHT_TYPE: %s\n",
                        edge_weight_type);
                exit(1);
            }
            break;
        }
        if (!strncmp(line, "DIMENSION", strlen("DIMENSION"))) {
            cp = strchr(line, ':');
            sscanf(cp + 1, "%d", &dimension);
        } else if (!strncmp(line, "NAME", strlen("NAME"))) {
            cp = strchr(line, ':');
            sscanf(cp + 1, "%s", name);
        }
    }
    while (fgets(line, 256, in)) {
        if (!strncmp(line, "NODE_COORD_SECTION",
            strlen("NODE_COORD_SECTION")))
            break;
        if (!strncmp(line, "DIMENSION", strlen("DIMENSION"))) {
            cp = strchr(line, ':');
            sscanf(cp + 1, "%d", &dimension);
        } else if (!strncmp(line, "NAME", strlen("NAME"))) {
            cp = strchr(line, ':');
            sscanf(cp + 1, "%s", name);
        }
    }
    if (dimension != n) {
        fprintf(stderr,
                "Node count does not match DIMENSION in TSPLIB file\n");
        exit(1);
    }
    for (i = 0; i < n; i++) {
        if (dist == euc3d)
            fscanf(in, "%d %lf %lf %lf\n", &id, &x, &y, &z);
        else
            fscanf(in, "%d %lf %lf\n", &id, &x, &y);
        node_set[id - 1].x = x;
        node_set[id - 1].y = y;
        node_set[id - 1].z = z;
    }
    for (i = 0; i < n; i++) {
        for (e = node_set[i].first_edge; e; e = e->next) {
            if (e->cost == INT_MIN)
                e->cost = dist(&node_set[i], e->to);
            else if (e->cost != dist(&node_set[i], e->to)) {
                fprintf(stderr, "Wrong cost given for node %d\n", i);
                exit(1);
            }
        }
    }
    fclose(in);
}

void read_tour_file(char * tour_file_name)
{
    FILE *in;
    int from = 0, to = 0, id, i;
    char line[256];

    if (!tour_file_name)
        return;
    in = fopen(tour_file_name, "r");
    if (!in) {
        fprintf(stderr, "Cannot open tour_file: %s\n", tour_file_name);
        exit(1);
    }
    while (fgets(line, 256, in)) {
       if (!strncmp(line, "TOUR_SECTION", strlen("TOUR_SECTION"))) {
           for (i = 0; i < n; i++, to = id) {
               fscanf(in, "%d\n", &id);
               if (id <= 0 || id > n) {
                   fprintf(stderr, "Node %d out of range in tour_file: %s\n",
                           id, tour_file_name);
                   exit(1);
               }
               if (i == 0)
                   from = id;
               else
                   node_set[to - 1].tour_suc = &node_set[id - 1];
           }
           node_set[id - 1].tour_suc = &node_set[from - 1];
       }
    }
    for (id = 0; id < n; id++) {
        if (!is_edge(&node_set[id], node_set[id].tour_suc))
            printf("WARNING: tour edge (%d, %d) is not in edge_set\n",
                   node_set[id].id + 1, node_set[id].tour_suc->id + 1);
    }
    fclose(in);
}

static char *full_name(char * name, int edges)
{
   char *new_name = 0, *buf, *pos;

      if (!(pos = strstr(name, "$"))) {
          new_name = (char *) calloc(strlen(name) + 1, 1);
          strcpy(new_name, name);
          return new_name;
      }
      buf = (char *) malloc(12);
      sprintf(buf, "%d", edges);
      do {
          free(new_name);
          new_name =
             (char *) calloc(strlen(name) + strlen(buf) + 1, 1);
          strncpy(new_name, name, pos - name);
          strcat(new_name, buf);
          strcat(new_name, pos + 1);
          name = new_name;
      }
      while ((pos = strstr(name, "$")));
      free(buf);
      return new_name; 
} 

void write_output_file(char * output_file_name)
{
    FILE *out;
    int i;
    char *full_file_name;
    edge *e;
    node *a;

    if (!output_file_name)
        return;
    full_file_name = full_name(output_file_name, edges - eliminated);
    out = fopen(full_file_name, "w");
    if (!out) {
         fprintf(stderr, "Cannot open output_file: %s\n", output_file_name);
         exit(1);
     }
     fprintf(out, "%d %d\n", n, edges - eliminated);
     for (i = 0, a = node_set; i < n; i++, a++) {
         for (e = a->first_edge; e; e = e->next)
             if (a->id < e->to->id) {
                 if (1 || concorde_edge_format)
                     fprintf(out, "%d %d %d\n", a->id, e->to->id, e->cost);
                 else
                     fprintf(out, "%d %d\n", a->id, e->to->id);
             }
      }
      fclose(out);
      printf(" (written to %s)", full_file_name);
      free(full_file_name);
}

void add_edge(node * a, node * b, int cost)
{
    edge *e;
    int i;

    for (i = 1; i <= 2; i++, swap(&a, &b)) {
        e = (edge *) malloc(sizeof(edge));
        e->to = b;
        e->cost = cost;
        e->next = a->first_edge;
        a->first_edge = e;
        a->degree++;
    }
}

int is_edge(node * a, node * b)
{
    edge *e;

    for (e = a->first_edge; e; e = e->next)
        if (e->to == b)
            return 1;
    return 0;
}

int is_tour_edge(node * a, node * b)
{
    return a->tour_suc == b || b->tour_suc == a;
}

int compare(const void *ea, const void *eb)
{
    return (*(edge **) ea)->cost - (*(edge **) eb)->cost;
}

void sort_edges()
{
    int i, j, count;
    node *a;
    edge *e;

    for (i = 0, a = node_set; i < n; i++, a++) {
        count = 0;
        for (e = a->first_edge; e; e = e->next)
            edge_set[count++] = e;
        qsort(edge_set, count, sizeof(edge *), compare);
        a->first_edge = e = edge_set[0];
        for (j = 1; j < count; j++)
            e = e->next = edge_set[j];
        e->next = 0;
    }
}

int opt24(node * a, node * b,
          node * c1, node * c2, node * c3, node * c4)
{
    node *p1[2] = { a, b };
    node *p2[4] = { c1, c2, c3, c4 };
    node **p[2] = { p1, p2 };
    int size[2] = { 2, 4 };
    return opt(2, p, size);
}

int opt33(node * c1, node * c2, node * c3,
          node * d1, node * d2, node * d3)
{
    node *p1[3] = { c1, c2, c3 };
    node *p2[3] = { d1, d2, d3 };
    node **p[2] = { p1, p2 };
    int size[2] = { 3, 3 };
    return opt(2, p, size);
}

int opt34(node * c1, node * c2, node * c3,
          node * d1, node * d2, node * d3, node * d4)
{
    node *p1[3] = { c1, c2, c3 };
    node *p2[4] = { d1, d2, d3, d4 };
    node **p[2] = { p1, p2 };
    int size[2] = { 3, 4 };
    return opt(2, p, size);
}

int opt222(node * a, node * b,
           node * c1, node * c2,
           node * d1, node * d2)
{
    node *p1[2] = { a, b };
    node *p2[2] = { c1, c2 };
    node *p3[2] = { d1, d2 };
    node **p[3] = { p1, p2, p3 };
    int size[3] = { 2, 2, 2 };
    return opt(3, p, size);
}

int opt232(node * a, node * b,
           node * c1, node * c2, node * c3,
           node * d1, node * d2)
{
    node *p1[2] = { a, b };
    node *p2[3] = { c1, c2, c3 };
    node *p3[3] = { d1, d2 };
    node **p[3] = { p1, p2, p3 };
    int size[3] = { 2, 3, 2 };
    return opt(3, p, size);
}

int opt233(node * a, node * b,
           node * c1, node * c2, node * c3,
           node * d1, node * d2, node * d3)
{
    node *p1[2] = { a, b };
    node *p2[3] = { c1, c2, c3 };
    node *p3[3] = { d1, d2, d3 };
    node **p[3] = { p1, p2, p3 };
    int size[3] = { 2, 3, 3 };
    return opt(3, p, size);
}

int opt243(node * a, node * b,
           node * c1, node * c2, node * c3, node * c4,
           node * d1, node * d2, node * d3)
{
    if (!opt23(c1, c2, c2, c3, c4, C(c1, c2), C(c2, c3), C(c3, c4)) ||
        !opt23(c1, c2, d1, d2, d3, C(c1, c2), C(d1, d2), C(d2, d3)) ||
        !opt23(d1, d2, c1, c2, c3, C(d1, d2), C(c1, c2), C(c2, c3)) ||
        !opt23(d2, d3, c1, c2, c3, C(d2, d3), C(c1, c2), C(c2, c3)) ||
        !opt23(a, b, c1, c2, c3, C(a, b), C(c1, c2), C(c2, c3)))
        return 0;
    node *p1[2] = { a, b };
    node *p2[4] = { c1, c2, c3, c4 };
    node *p3[3] = { d1, d2, d3 }; 
    node **p[3] = { p1, p2, p3 };
    int size[3] = { 2, 4, 3 };
    return opt(3, p, size);
}

int opt244(node * a, node * b,
           node * c1, node * c2, node * c3, node * c4,
           node * d1, node * d2, node * d3, node * d4)
{
    if (!opt23(c1, c2, c2, c3, c4, C(c1, c2), C(c2, c3), C(c3, c4)) ||
        !opt23(c1, c2, d1, d2, d3, C(c1, c2), C(d1, d2), C(d2, d3)) ||
        !opt23(c1, c2, d2, d3, d4, C(c1, c2), C(d2, d3), C(d3, d4)) ||
        !opt23(d1, d2, c1, c2, c3, C(d1, d2), C(c1, c2), C(c2, c3)) ||
        !opt23(d2, d3, c1, c2, c3, C(d2, d3), C(c1, c2), C(c2, c3)) ||
        !opt23(d3, d4, c1, c2, c3, C(d3, d4), C(c1, c2), C(c2, c3)) ||
        !opt23(a, b, c1, c2, c3, C(a, b), C(c1, c2), C(c2, c3)))
        return 0;
    node *p1[2] = { a, b };
    node *p2[4] = { c1, c2, c3, c4 };
    node *p3[4] = { d1, d2, d3, d4 };
    node **p[3] = { p1, p2, p3 };
    int size[3] = { 2, 4, 4 };
    return opt(3, p, size);
}

int opt253(node * a, node * b,
           node * c1, node * c2, node * c3, node * c4, node * c5,
           node * d1, node * d2, node * d3)
{
    if (!opt23(c1, c2, c2, c3, c4, C(c1, c2), C(c2, c3), C(c3, c4)) ||
        !opt23(c1, c2, c3, c4, c5, C(c1, c2), C(c3, c4), C(c4, c5)) ||
        !opt23(c1, c2, d1, d2, d3, C(c1, c2), C(d1, d2), C(d2, d3)) ||
        !opt23(d1, d2, c1, c2, c3, C(d1, d2), C(c1, c2), C(c2, c3)) ||
        !opt23(d2, d3, c1, c2, c3, C(d2, d3), C(c1, c2), C(c2, c3)) ||
        !opt23(a, b, c1, c2, c3, C(a, b), C(c1, c2), C(c2, c3)))
        return 0;
    node *p1[2] = { a, b };
    node *p2[5] = { c1, c2, c3, c4, c5 };
    node *p3[3] = { d1, d2, d3 };
    node **p[3] = { p1, p2, p3 };
    int size[3] = { 2, 5, 3 };
    return opt(3, p, size);
}

int opt333(node * a, node * b, node * b2,
           node * c1, node * c2, node * c3,
           node * d1, node * d2, node * d3)
{
    node *p1[3] = { a, b, b2 };
    node *p2[3] = { c1, c2, c3 };
    node *p3[3] = { d1, d2, d3 };
    node **p[3] = { p1, p2, p3 };
    int size[3] = { 3, 3, 3 };
    return opt(3, p, size);
}

int opt343(node * c1, node * c2, node * c3,
           node * d1, node * d2, node * d3, node * d4,
           node * e1, node * e2, node * e3)
{
    if (!opt23(d1, d2, d2, d3, d4, C(d1, d2), C(d2, d3), C(d3, d4)) ||
        !opt23(d1, d2, c1, c2, c3, C(d1, d2), C(c1, c2), C(c2, c3)) ||
        !opt23(d1, d2, e1, e2, e3, C(d1, d2), C(e1, e2), C(e2, e3)) ||
        !opt23(c1, c2, d1, d2, d3, C(c1, c2), C(d1, d2), C(d2, d3)) ||
        !opt23(c2, c3, d1, d2, d3, C(c2, c3), C(d1, d2), C(d2, d3)) ||
        !opt23(e1, e2, d1, d2, d3, C(e1, e2), C(d1, d2), C(d2, d3)) ||
        !opt23(e2, e3, d1, d2, d3, C(e2, e3), C(d1, d2), C(d2, d3)) ||
        !opt23(c2, c3, d2, d3, d4, C(c2, c3), C(d2, d3), C(d3, d4)) ||
        !opt23(c2, c3, e1, e2, e3, C(c2, c3), C(e1, e2), C(e2, e3)))
        return 0;
    node *p1[3] = { c1, c2, c3 };
    node *p2[4] = { d1, d2, d3, d4 };
    node *p3[3] = { e1, e2, e3 };
    node **p[3] = { p1, p2, p3 };
    int size[3] = { 3, 4, 3 };
    return opt(3, p, size);
}

int opt443(node * c1, node * c2, node * c3, node * c4,
           node * d1, node * d2, node * d3, node * d4,
           node * e1, node * e2, node * e3)
{
    if (!opt23(c1, c2, c2, c3, c4, C(c1, c2), C(c2, c3), C(c3, c4)) ||
        !opt23(c1, c2, d1, d2, d3, C(c1, c2), C(d1, d2), C(d2, d3)) ||
        !opt23(c1, c2, d2, d3, d4, C(c1, c2), C(d2, d3), C(d3, d4)) ||
        !opt23(d1, d2, c1, c2, c3, C(d1, d2), C(c1, c2), C(c2, c3)) ||
        !opt23(d2, d3, c1, c2, c3, C(d2, d3), C(c1, c2), C(c2, c3)) ||
        !opt23(d3, d4, c1, c2, c3, C(d3, d4), C(c1, c2), C(c2, c3)) ||
        !opt23(e1, e2, c1, c2, c3, C(e1, e2), C(c1, c2), C(c2, c3)) ||
        !opt23(e2, e3, c1, c2, c3, C(e2, e3), C(c1, c2), C(c2, c3)))
        return 0;
    node *p1[4] = { c1, c2, c3, c4 };
    node *p2[4] = { d1, d2, d3, d4 };
    node *p3[3] = { e1, e2, e3 };
    node **p[3] = { p1, p2, p3 };
    int size[3] = { 4, 4, 3 };
    return opt(3, p, size);
}

int opt2233(node * a, node * b,
            node * c1, node * c2,
            node * d1, node * d2, node * d3,
            node * e1, node * e2, node * e3)
{
    node *p1[2] = { a, b };
    node *p2[2] = { c1, c2 };
    node *p3[3] = { d1, d2, d3 };
    node *p4[3] = { e1, e2, e3 };
    node **p[4] = { p1, p2, p3, p4 };
    int size[4] = { 2, 2, 3, 3 };
    return opt(4, p, size);
}

int opt2333(node * a, node * b,
            node * c1, node * c2, node * c3,
            node * d1, node * d2, node * d3,
            node * e1, node * e2, node * e3)
{
    node *p1[2] = { a, b };
    node *p2[3] = { c1, c2, c3 };
    node *p3[3] = { d1, d2, d3 };
    node *p4[3] = { e1, e2, e3 };
    node **p[4] = { p1, p2, p3, p4 };
    int size[4] = { 2, 3, 3, 3 };
    return opt(4, p, size);
}

int opt2433(node * a, node * b,
            node * c1, node * c2, node * c3, node * c4,
            node * d1, node * d2, node * d3,
            node * e1, node * e2, node * e3)
{
    if (!opt23(c1, c2, c2, c3, c4, C(c1, c2), C(c2, c3), C(c3, c4)) ||
        !opt23(c1, c2, d1, d2, d3, C(c1, c2), C(d1, d2), C(d2, d3)) ||
        !opt23(c1, c2, e1, e2, e3, C(c1, c2), C(e1, e2), C(e2, e3)) ||
        !opt23(d1, d2, c1, c2, c3, C(d1, d2), C(c1, c2), C(c2, c3)) ||
        !opt23(d2, d3, c1, c2, c3, C(d2, d3), C(c1, c2), C(c2, c3)) ||
        !opt23(e1, e2, c1, c2, c3, C(e1, e2), C(c1, c2), C(c2, c3)) ||
        !opt23(e2, e3, c1, c2, c3, C(e2, e3), C(c1, c2), C(c2, c3)) ||
        !opt23(a, b, c1, c2, c3, C(a, b), C(c1, c2), C(c2, c3)))
        return 0;
    node *p1[2] = { a, b };
    node *p2[4] = { c1, c2, c3, c4 };
    node *p3[3] = { d1, d2, d3 };
    node *p4[3] = { e1, e2, e3 };
    node **p[4] = { p1, p2, p3, p4 };
    int size[4] = { 2, 4, 3, 3 };
    return opt(4, p, size);
}

int opt2443(node * a, node * b,
            node * c1, node * c2, node * c3, node * c4,
            node * d1, node * d2, node * d3, node * d4,
            node * e1, node * e2, node * e3)
{
    if (!opt23(c1, c2, c2, c3, c4, C(c1, c2), C(c2, c3), C(c3, c4)) ||
        !opt23(c1, c2, d1, d2, d3, C(c1, c2), C(d1, d2), C(d2, d3)) ||
        !opt23(c1, c2, d2, d3, d4, C(c1, c2), C(d2, d3), C(d3, d4)) ||
        !opt23(c1, c2, e1, e2, e3, C(c1, c2), C(e1, e2), C(e2, e3)) ||
        !opt23(d1, d2, c1, c2, c3, C(d1, d2), C(c1, c2), C(c2, c3)) ||
        !opt23(d2, d3, c1, c2, c3, C(d2, d3), C(c1, c2), C(c2, c3)) ||
        !opt23(d3, d4, c1, c2, c3, C(d3, d4), C(c1, c2), C(c2, c3)) ||
        !opt23(e1, e2, c1, c2, c3, C(e1, e2), C(c1, c2), C(c2, c3)) ||
        !opt23(e2, e3, c1, c2, c3, C(e2, e3), C(c1, c2), C(c2, c3)) ||
        !opt23(a, b, c1, c2, c3, C(a, b), C(c1, c2), C(c2, c3)))
        return 0;
    node *p1[2] = { a, b };
    node *p2[4] = { c1, c2, c3, c4 };
    node *p3[4] = { d1, d2, d3, d4 };
    node *p4[3] = { e1, e2, e3 };
    node **p[4] = { p1, p2, p3, p4 };
    int size[4] = { 2, 4, 4, 3 };
    return opt(4, p, size);
}

int opt2533(node * a, node * b,
            node * c1, node * c2, node * c3, node * c4, node * c5,
            node * d1, node * d2, node * d3,
            node * e1, node * e2, node * e3)
{
    if (!opt23(c1, c2, c2, c3, c4, C(c1, c2), C(c2, c3), C(c3, c4)) ||
        !opt23(c1, c2, c3, c4, c5, C(c1, c2), C(c3, c4), C(c4, c5)) ||
        !opt23(c1, c2, d1, d2, d3, C(c1, c2), C(d1, d2), C(d2, d3)) ||
        !opt23(c1, c2, e1, e2, e3, C(c1, c2), C(e1, e2), C(e2, e3)) ||
        !opt23(d1, d2, c1, c2, c3, C(d1, d2), C(c1, c2), C(c2, c3)) ||
        !opt23(d2, d3, c1, c2, c3, C(d2, d3), C(c1, c2), C(c2, c3)) ||
        !opt23(e1, e2, c1, c2, c3, C(e1, e2), C(c1, c2), C(c2, c3)) ||
        !opt23(e2, e3, c1, c2, c3, C(e2, e3), C(c1, c2), C(c2, c3)) ||
        !opt23(a, b, c1, c2, c3, C(a, b), C(c1, c2), C(c2, c3)))
        return 0;
    node *p1[2] = { a, b };
    node *p2[5] = { c1, c2, c3, c4, c5 };
    node *p3[3] = { d1, d2, d3 };
    node *p4[3] = { e1, e2, e3 };
    node **p[4] = { p1, p2, p3, p4 };
    int size[4] = { 2, 5, 3, 3 };
    return opt(4, p, size);
}

int opt3333(node * c1, node * c2, node * c3,
            node * d1, node * d2, node * d3,
            node * e1, node * e2, node * e3,
            node * f1, node * f2, node * f3)
{
    node *p1[3] = { c1, c2, c3 };
    node *p2[3] = { d1, d2, d3 };
    node *p3[3] = { e1, e2, e3 };
    node *p4[3] = { f1, f2, f3 };
    node **p[4] = { p1, p2, p3, p4 };
    int size[4] = { 3, 3, 3, 3 };
    return opt(4, p, size);
}

int opt3433(node * c1, node * c2, node * c3,
            node * d1, node * d2, node * d3, node * d4,
            node * e1, node * e2, node * e3,
            node * f1, node * f2, node * f3)
{
    if (!opt23(d1, d2, d2, d3, d4, C(d1, d2), C(d2, d3), C(d3, d4)) ||
        !opt23(d1, d2, c1, c2, c3, C(d1, d2), C(c1, c2), C(c2, c3)) ||
        !opt23(d1, d2, e1, e2, e3, C(d1, d2), C(e1, e2), C(e2, e3)) ||
        !opt23(d1, d2, f1, f2, f3, C(d1, d2), C(f1, f2), C(f2, f3)) ||
        !opt23(c1, c2, d1, d2, d3, C(c1, c2), C(d1, d2), C(d2, d3)) ||
        !opt23(c2, c3, d1, d2, d3, C(c2, c3), C(d1, d2), C(d2, d3)) ||
        !opt23(e1, e2, d1, d2, d3, C(e1, e2), C(d1, d2), C(d2, d3)) ||
        !opt23(e2, e3, d1, d2, d3, C(e2, e3), C(d1, d2), C(d2, d3)) ||
        !opt23(f1, f2, d1, d2, d3, C(f1, f2), C(d1, d2), C(d2, d3)) ||
        !opt23(f2, f3, d1, d2, d3, C(f2, f3), C(d1, d2), C(d2, d3)))
        return 0;
    node *p1[3] = { c1, c2, c3 };
    node *p2[4] = { d1, d2, d3, d4 };
    node *p3[3] = { e1, e2, e3 };
    node *p4[3] = { f1, f2, f3 };
    node **p[4] = { p1, p2, p3, p4 };
    int size[4] = { 3, 4, 3, 3 };
    return opt(4, p, size);
}

int opt22333(node * a, node * b,
             node * c1, node * c2,
             node * d1, node * d2, node * d3,
             node * e1, node * e2, node * e3,
             node * f1, node * f2, node * f3)
{
    node *p1[2] = { a, b };
    node *p2[3] = { c1, c2 };
    node *p3[3] = { d1, d2, d3 };
    node *p4[3] = { e1, e2, e3 };
    node *p5[3] = { f1, f2, f3 };
    node **p[5] = { p1, p2, p3, p4, p5 };
    int size[5] = { 2, 2, 3, 3, 3 };
    return opt(5, p, size);
}

int opt23333(node * a, node * b,
             node * c1, node * c2, node * c3,
             node * d1, node * d2, node * d3,
             node * e1, node * e2, node * e3,
             node * f1, node * f2, node * f3)
{
    node *p1[2] = { a, b };
    node *p2[3] = { c1, c2, c3 };
    node *p3[3] = { d1, d2, d3 };
    node *p4[3] = { e1, e2, e3 };
    node *p5[3] = { f1, f2, f3 };
    node **p[5] = { p1, p2, p3, p4, p5 };
    int size[5] = { 2, 3, 3, 3, 3 };
    return opt(5, p, size);
}

int opt24333(node * a, node * b,
             node * c1, node * c2, node * c3, node * c4,
             node * d1, node * d2, node * d3,
             node * e1, node * e2, node * e3,
             node * f1, node * f2, node * f3)
{
    if (!opt23(c1, c2, c2, c3, c4, C(c1, c2), C(c2, c3), C(c3, c4)) ||
        !opt23(c1, c2, d1, d2, d3, C(c1, c2), C(d1, d2), C(d2, d3)) ||
        !opt23(c1, c2, e1, e2, e3, C(c1, c2), C(e1, e2), C(e2, e3)) ||
        !opt23(c1, c2, f1, f2, f3, C(c1, c2), C(f1, f2), C(f2, f3)) ||
        !opt23(d1, d2, c1, c2, c3, C(d1, d2), C(c1, c2), C(c2, c3)) ||
        !opt23(d2, d3, c1, c2, c3, C(d2, d3), C(c1, c2), C(c2, c3)) ||
        !opt23(e1, e2, c1, c2, c3, C(e1, e2), C(c1, c2), C(c2, c3)) ||
        !opt23(e2, e3, c1, c2, c3, C(e2, e3), C(c1, c2), C(c2, c3)) ||
        !opt23(a, b, c1, c2, c3, C(a, b), C(c1, c2), C(c2, c3)))
        return 0;
    node *p1[2] = { a, b };
    node *p2[4] = { c1, c2, c3, c4 };
    node *p3[3] = { d1, d2, d3 };
    node *p4[3] = { e1, e2, e3 };
    node *p5[3] = { f1, f2, f3 };
    node **p[5] = { p1, p2, p3, p4, p5 };
    int size[5] = { 2, 4, 3, 3, 3 };
    return opt(5, p, size);
}

int opt24433(node * a, node * b,
             node * c1, node * c2, node * c3, node * c4,
             node * d1, node * d2, node * d3, node * d4,
             node * e1, node * e2, node * e3,
             node * f1, node * f2, node * f3)
{
    if (!opt23(c1, c2, c2, c3, c4, C(c1, c2), C(c2, c3), C(c3, c4)) ||
        !opt23(c1, c2, d1, d2, d3, C(c1, c2), C(d1, d2), C(d2, d3)) ||
        !opt23(c1, c2, d2, d3, d4, C(c1, c2), C(d2, d3), C(d3, d4)) ||
        !opt23(c1, c2, e1, e2, e3, C(c1, c2), C(e1, e2), C(e2, e3)) ||
        !opt23(c1, c2, f1, f2, f3, C(c1, c2), C(f1, f2), C(f2, f3)) ||
        !opt23(d1, d2, c1, c2, c3, C(d1, d2), C(c1, c2), C(c2, c3)) ||
        !opt23(d2, d3, c1, c2, c3, C(d2, d3), C(c1, c2), C(c2, c3)) ||
        !opt23(d3, d4, c1, c2, c3, C(d3, d4), C(c1, c2), C(c2, c3)) ||
        !opt23(e1, e2, c1, c2, c3, C(e1, e2), C(c1, c2), C(c2, c3)) ||
        !opt23(e2, e3, c1, c2, c3, C(e2, e3), C(c1, c2), C(c2, c3)) ||
        !opt23(a, b, c1, c2, c3, C(a, b), C(c1, c2), C(c2, c3)))
        return 0;
    node *p1[2] = { a, b };
    node *p2[4] = { c1, c2, c3, c4 };
    node *p3[4] = { d1, d2, d3, d4 };
    node *p4[3] = { e1, e2, e3 };
    node *p5[3] = { f1, f2, f3 };
    node **p[5] = { p1, p2, p3, p4, p5 };
    int size[5] = { 2, 4, 4, 3, 3 };
    return opt(5, p, size);
}

int opt25333(node * a, node * b,
             node * c1, node * c2, node * c3, node * c4, node * c5,
             node * d1, node * d2, node * d3,
             node * e1, node * e2, node * e3,
             node * f1, node * f2, node * f3)
{
    if (!opt23(c1, c2, c2, c3, c4, C(c1, c2), C(c2, c3), C(c3, c4)) ||
        !opt23(c1, c2, c3, c4, c5, C(c1, c2), C(c3, c4), C(c4, c5)) ||
        !opt23(c1, c2, d1, d2, d3, C(c1, c2), C(d1, d2), C(d2, d3)) ||
        !opt23(c1, c2, e1, e2, e3, C(c1, c2), C(e1, e2), C(e2, e3)) ||
        !opt23(c1, c2, f1, f2, f3, C(c1, c2), C(f1, f2), C(f2, f3)) ||
        !opt23(d1, d2, c1, c2, c3, C(d1, d2), C(c1, c2), C(c2, c3)) ||
        !opt23(d2, d3, c1, c2, c3, C(d2, d3), C(c1, c2), C(c2, c3)) ||
        !opt23(e1, e2, c1, c2, c3, C(e1, e2), C(c1, c2), C(c2, c3)) ||
        !opt23(e2, e3, c1, c2, c3, C(e2, e3), C(c1, c2), C(c2, c3)) ||
        !opt23(f1, f2, c1, c2, c3, C(f1, f2), C(c1, c2), C(c2, c3)) ||
        !opt23(f2, f3, c1, c2, c3, C(f2, f3), C(c1, c2), C(c2, c3)) ||
        !opt23(a, b, c1, c2, c3, C(a, b), C(c1, c2), C(c2, c3)))
        return 0;
    node *p1[2] = { a, b };
    node *p2[5] = { c1, c2, c3, c4, c5 };
    node *p3[3] = { d1, d2, d3 };
    node *p4[3] = { e1, e2, e3 };
    node *p5[3] = { f1, f2, f3 };
    node **p[5] = { p1, p2, p3, p4, p5 };
    int size[5] = { 2, 5, 3, 3, 3 };
    return opt(5, p, size);
}

int opt33333(node * a, node * b, node * b2,
             node * c1, node * c2, node * c3,
             node * d1, node * d2, node * d3,
             node * e1, node * e2, node * e3,
             node * f1, node * f2, node * f3)
{
    node *p1[3] = { a, b, b2 };
    node *p2[3] = { c1, c2, c3 };
    node *p3[3] = { d1, d2, d3 };
    node *p4[3] = { e1, e2, e3 };
    node *p5[3] = { f1, f2, f3 };
    node **p[5] = { p1, p2, p3, p4, p5 };
    int size[5] = { 3, 3, 3, 3, 3 };
    return opt(5, p, size);
}

int opt34333(node * c1, node * c2, node * c3,
             node * d1, node * d2, node * d3, node * d4,
             node * e1, node * e2, node * e3,
             node * f1, node * f2, node * f3,
             node * g1, node * g2, node * g3)
{
    if (!opt23(d1, d2, d2, d3, d4, C(d1, d2), C(d2, d3), C(d3, d4)) ||
        !opt23(d1, d2, c1, c2, c3, C(d1, d2), C(c1, c2), C(c2, c3)) ||
        !opt23(d1, d2, e1, e2, e3, C(d1, d2), C(e1, e2), C(e2, e3)) ||
        !opt23(d1, d2, f1, f2, f3, C(d1, d2), C(f1, f2), C(f2, f3)) ||
        !opt23(d1, d2, g1, g2, g3, C(d1, d2), C(g1, g2), C(g2, g3)) ||
        !opt23(c1, c2, d1, d2, d3, C(c1, c2), C(d1, d2), C(d2, d3)) ||
        !opt23(c2, c3, d1, d2, d3, C(c2, c3), C(d1, d2), C(d2, d3)) ||
        !opt23(e1, e2, d1, d2, d3, C(e1, e2), C(d1, d2), C(d2, d3)) ||
        !opt23(e2, e3, d1, d2, d3, C(e2, e3), C(d1, d2), C(d2, d3)) ||
        !opt23(f1, f2, d1, d2, d3, C(f1, f2), C(d1, d2), C(d2, d3)) ||
        !opt23(f2, f3, d1, d2, d3, C(f2, f3), C(d1, d2), C(d2, d3)) ||
        !opt23(g1, g2, d1, d2, d3, C(g1, g2), C(d1, d2), C(d2, d3)) ||
        !opt23(g2, g3, d1, d2, d3, C(g2, g3), C(d1, d2), C(d2, d3)))
        return 0;
    node *p1[3] = { c1, c2, c3 };
    node *p2[4] = { d1, d2, d3, d4 };
    node *p3[3] = { e1, e2, e3 };
    node *p4[3] = { f1, f2, f3 };
    node *p5[3] = { g1, g2, g3 };
    node **p[5] = { p1, p2, p3, p4, p5 };
    int size[5] = { 3, 4, 3, 3, 3 };
    return opt(5, p, size);
}

int opt223333(node * a, node * b,
              node * c1, node * c2,
              node * d1, node * d2, node * d3,
              node * e1, node * e2, node * e3,
              node * f1, node * f2, node * f3,
              node * g1, node * g2, node * g3)
{
    node *p1[2] = { a, b };
    node *p2[2] = { c1, c2 };
    node *p3[3] = { d1, d2, d3 };
    node *p4[3] = { e1, e2, e3 };
    node *p5[3] = { f1, f2, f3 };
    node *p6[3] = { g1, g2, g3 };
    node **p[6] = { p1, p2, p3, p4, p5, p6 };
    int size[6] = { 2, 2, 3, 3, 3, 3 };
    return opt(6, p, size);
}

int opt233333(node * a, node * b,
              node * c1, node * c2, node * c3,
              node * d1, node * d2, node * d3,
              node * e1, node * e2, node * e3,
              node * f1, node * f2, node * f3,
              node * g1, node * g2, node * g3)
{
    node *p1[2] = { a, b };
    node *p2[3] = { c1, c2, c3 };
    node *p3[3] = { d1, d2, d3 };
    node *p4[3] = { e1, e2, e3 };
    node *p5[3] = { f1, f2, f3 };
    node *p6[3] = { g1, g2, g3 };
    node **p[6] = { p1, p2, p3, p4, p5, p6 };
    int size[6] = { 2, 3, 3, 3, 3, 3 };
    return opt(6, p, size);
}

int can_HS_eliminate(node * a, node * b,
                     node * c, node * d,
                     int Cab)
{
    node *c1, *c2, *d1, *d2;
    edge *ec1, *ec2, *ed1, *ed2;
    int Cc1c, Ccc2, Cd1d, Cdd2;

    for (ec1 = c->first_edge; ec1; ec1 = ec1->next) {
        c1 = ec1->to;
        if (c1 == d)
            continue;
        Cc1c = ec1->cost;
        if (!opt22(c1, c, a, b, Cc1c, Cab))
            continue;
        for (ec2 = ec1->next; ec2; ec2 = ec2->next) {
            c2 = ec2->to;
            if (c2 == d ||
                (c2 == a && c1 == b) || (c2 == b && c1 == a))
                continue;
            Ccc2 = ec2->cost;
            if (!opt22(c, c2, a, b, Ccc2, Cab) ||
                !opt23(a, b, c1, c, c2, Cab, Cc1c, Ccc2))
                continue;
            for (ed1 = d->first_edge; ed1; ed1 = ed1->next) {
                d1 = ed1->to;
                if (d1 == c ||
                    ((d1 == a || d1 == b) && (c1 == d1 || c2 == d1)))
                    continue;
                Cd1d = ed1->cost;
                if (!opt22(d, d1, a, b, Cd1d, Cab) ||
                    !opt22(d, d1, c1, c, Cd1d, Cc1c) ||
                    !opt22(d, d1, c, c2, Cd1d, Ccc2) ||
                    !opt23(d, d1, c1, c, c2, Cd1d, Cc1c, Ccc2) ||
                    !opt222(a, b, c1, c, d, d1) ||
                    !opt222(a, b, c, c2, d, d1) ||
                    !opt232(a, b, c1, c,  c2, d1, d))
                    continue;
                for (ed2 = ed1->next; ed2; ed2 = ed2->next) {
                    d2 = ed2->to;
                    if (d2 == c ||
                        (d2 == a && (d1 == b || c1 == d2 || c2 == d2)) ||
                        (d2 == b && (d1 == a || c1 == d2 || c2 == d2)) ||
                        (d2 == c1 && (d1 == c2 || a == d2 || b == d2)) ||
                        (d2 == c2 && (d1 == c1 || a == d2 || b == d2)) ||
                        has_cycle222(a, b, c1, c2, d1, d2))
                        continue;
                    Cdd2 = ed2->cost;
                    if (opt22(d, d2, a, b, Cdd2, Cab) &&
                        opt22(d, d2, c1, c, Cdd2, Cc1c) &&
                        opt22(d, d2, c, c2, Cdd2, Ccc2) &&
                        opt23(a, b, d1, d, d2, Cab, Cd1d, Cdd2) &&
                        opt23(d, d2, c1, c, c2, Cdd2, Cc1c, Ccc2) &&
                        opt23(c1, c, d1, d, d2, Cc1c, Cd1d, Cdd2) &&
                        opt23(c, c2, d1, d, d2, Ccc2, Cd1d, Cdd2) &&
                        opt233(a, b, c1, c, c2, d1, d, d2) &&
                        extra_edge_opt(a, b, c1, c, c2, d1, d, d2,
                                       0, 0, 0, 0, 0, 0) &&
                        extra_node_opt(a, b, c1, c, c2, d1, d, d2,
                                       0, 0, 0, 0, 0, 0,
                                       Cab, Cc1c, Ccc2, Cd1d, Cdd2,
                                       0, 0, 0, 0))
                        return 0;
                }
            }
        }
    }
    return 1;
}

int compatible(node * a, node * b,
               node * c, node * d,
               int Cab,
               int Ccd)
{
    int t = Cab + Ccd;

    return c == a || c == b || d == a || d == b ||
           Ccd > ab_stretch * Cab ||
           (is_edge(c, d) &&
            (t <= C(a, c) + C(d, b) ||
             t <= C(a, d) + C(c, b)));
}

void swap(node ** a, node ** b)
{
    node *temp = *a;
    *a = *b;
    *b = temp;
}

void swapC(int * a, int * b)
{
    int temp = *a;
    *a = *b;
    *b = temp;
}

int fixed(node * a, node * b) 
{
    return (a->degree == 2 || b->degree == 2) && is_edge(a, b);
}

int min(int a, int b)
{
    return a < b ? a : b;
}

int opt22(node * a, node * b,
          node * c, node * d,
          int Cab,
          int Ccd)
{
    return c == a || c == b || d == a || d == b ||
           Cab + Ccd <= C(a, c) + C(d, b) ||
           Cab + Ccd <= C(a, d) + C(c, b);
}

int opt23_simple(node * a, node * b,
                 node * c1, node * c, node * c2,
                 int Cab,
                 int Cc1c, int Ccc2)
{
    return (c1 != a || c2 != b) && (c1 != b || c2 != a) &&
           Cab + Cc1c + Ccc2 <= C(a, c) + C(c, b) + C(c1, c2);
}

int opt23(node * a, node * b,
          node * c1, node * c, node * c2,
          int Cab,
          int Cc1c, int Ccc2)
{
    node *z;
    edge *e;
    int i, j;

    if (!opt23_simple(a, b, c1, c, c2, Cab, Cc1c, Ccc2))
        return 0;
    for (e = c->first_edge; e; e = e->next)
        if (e->to != c1 && e->to != c2 && e->to->degree == 2)
            return 0;
    if (fixed(c1, c2))
        return 0;
    if (!strong_3_opt)
        return 1;
    for (i = 1; i <= 2; i++, swap(&a, &b)) {
        for (j = 1; j <= 2; j++, swap(&c1, &c2), swapC(&Cc1c, &Ccc2)) {
            if (b == c1) {
                if (fixed(a, c2))
                    return 0;
                for (e = b->first_edge; e; e = e->next) {
                    z = e->to;
                    if (z != a && z != c && z != c2 &&
                        !opt_excess(z, a, c1, c, c2, Cab, Cc1c, Ccc2))
                        return 0;
                }
                for (e = c->first_edge; e; e = e->next) {
                    z = e->to;
                    if (z != a && z != b && z != c1 && z != c2 &&
                        !is_edge(b, z) &&
                        !opt_excess(z, a, c1, c, c2, Cab, Cc1c, Ccc2))
                        return 0;
                }
                return 1;
            } else if (c2 == a && fixed(b, c1))
                return 0;
        }
    }
    return 1;
}

int opt_excess(node * z,
               node * s1, node * s2, node * s3, node * s4,
               int Cs1s2, int Cs2s3, int Cs3s4)
{
    node *z1, *z2;
    edge *e1, *e2;
    int Cs2z = C(s2, z), Cs3z = C(s3, z);
    int bound = C(s1, s4) - Cs1s2 - Cs3s4;

    for (e1 = z->first_edge; e1; e1 = e1->next) {
        z1 = e1->to;
        if (z1 == s2 || z1 == s3 ||
            e1->cost - Cs2z - C(s3, z1) > bound ||
            e1->cost - Cs3z - C(s2, z1) > bound ||
            !opt24(z1, z, s1, s2, s3, s4))
            continue;
        for (e2 = e1->next; e2; e2 = e2->next) {
            z2 = e2->to;
            if (z2 != s2 && z2 != s3 &&
                (z1 != s1 || z2 != s4) && (z1 != s4 || z2 != s1) &&
                opt34(z1, z, z2, s1, s2, s3, s4))
                return 1;
        }
    }
    return 0;
}

void reverse_path(node **p, int size)
{
    int i = 0, j = size - 1;
    while (i < j)
        swap(&p[i++], &p[j--]);
}

int is_opt(int paths, node *** path, int * size) {
    int i, j, n, result;
    node **p;
    int *fix;

    n = 0;
    for (i = 0; i < paths; i++)
        n += size[i];
    p = (node **) malloc(n * sizeof(node *));
    fix = (int *) calloc(n, sizeof(int));
    n = 0;
    for (i = 0; i < paths; i++) {
        for (j = 0; j < size[i]; j++)
            p[n++] = path[i][j];
        fix[n - 1] = 1;
    }
    fix[n - 1] = 0;
    result = path_opt(p, n, fix);
    free(fix);
    free(p);
    return result;
}

int reverse_opt(int paths, node ***path, int * size, int n) 
{
    int result;

    if (n <= 1)
        return is_opt(paths, path, size);
    if (reverse_opt(paths, path, size, n - 1))
        return 1;
    reverse_path(path[n - 1], size[n - 1]);
    result = reverse_opt(paths, path, size, n - 1);
    reverse_path(path[n - 1], size[n - 1]);
    return result;
}

int perm_opt(int paths, node *** path, int * size, int n)
{
    int i, s;
    node **t;

    if (n <= 1)
        return reverse_opt(paths, path, size, paths);
    for (i = 1; i < n; i++) {
        t = path[i]; path[i] = path[n - 1]; path[n - 1] = t;
        s = size[i]; size[i] = size[n - 1]; size[n - 1] = s;
        if (perm_opt(paths, path, size, n - 1))
            return 1;
        t = path[i]; path[i] = path[n - 1]; path[n - 1] = t;
        s = size[i]; size[i] = size[n - 1]; size[n - 1] = s;
    }
    return 0;
}

int opt(int paths, node *** path, int * size)
{
    node *p[MAX_N], *a1, *a2, *b1, *b2;
    int i, j, k, n;

    for (i = 1; i < paths; i++) {
        b1 = path[i][0];
        b2 = path[i][size[i] - 1];
        if ((b1 == b2 && size[i] > 1) ||
            (size[i] > 2 && fixed(b1, b2)))
            return 0;
        for (j = 0; j < i; j++) {
            a1 = path[j][0];
            a2 = path[j][size[j] - 1];
            if (j == 0 &&
                ((a1 == a2 && size[j] > 1) || (fixed(a1, a2) && size[j] > 2)))
                return 0;
            n = 0;
            if (a2 == b1 || (a2 != b2 && b1 != a1 && fixed(a2, b1))) {
                if (a1 == b2)
                    return 0;
                for (k = 0; k < size[j]; k++)
                    p[n++] = path[j][k];
                for (k = (a2 == b1 ? 1 : 0); k < size[i]; k++)
                    p[n++] = path[i][k];
            } else if (a2 == b2 || (a2 != b1 && b2 != a1 && fixed(a2, b2))) {
                if (a1 == b1)
                    return 0;
                for (k = 0; k < size[j]; k++)
                    p[n++] = path[j][k];
                for (k = size[i] - (a2 == b2 ? 2 : 1); k >= 0; k--)
                    p[n++] = path[i][k];
            } else if (a1 == b1 || (a1 != b2 && b1 != a2 && fixed(a1, b1))) {
                if (a2 == b2)
                    return 0;
                for (k = size[j] - 1; k >= 0; k--)
                    p[n++] = path[j][k];
                for (k = (a1 == b1 ? 1 : 0); k < size[i]; k++)
                    p[n++] = path[i][k];
            } else if (a1 == b2 || (a1 != b1 && b2 != a2 && fixed(a1, b2))) {
                if (a2 == b1)
                    return 0;
                for (k = size[j] - 1; k >= 0; k--)
                    p[n++] = path[j][k];
                for (k = size[i] - (a1 == b2 ? 2 : 1); k >= 0; k--)
                   p[n++] = path[i][k];
            }
            if (n > 0) {
                path[j] = p;
                size[j] = n;
                path[i] = path[paths - 1];
                size[i] = size[paths - 1];
                return opt(paths - 1, path, size);
            }
        }
    }
    for (i = 0; i < paths; i++) {
        p[2 * i] = path[i][0];
        p[2 * i + 1] = path[i][size[i] - 1];
    }
    if (has_cycle(2 * paths, p))
        return 0;
    return paths > 3 || perm_opt(paths, path, size, paths);
}

int extra_edge_opt(node * a, node * b,
                   node * c1, node * c, node * c2,
                   node * d1, node * d, node * d2,
                   node * e1, node * e, node * e2,
                   node * f1, node * f, node * f2)
{
    int j, i, state = 0;
    edge *ee;

    if (extra_edges < 1)
        return 1;
    for (j = 1; j <= 2; j++, swap(&c1, &d1), swap(&c, &d), swap(&c2, &d2)) {
        for (i = 1; i <= 2; i++, swap(&c1, &c2)) {
            state++;
            if (c1 == a || c1 == b || c1 == d1 || c1 == d2 ||
                c1 == e1 || c1 == e2 || c1 == f1 || c1 == f2)
                continue;
            for (ee = c1->first_edge; ee; ee = ee->next)
                if (ee->to != c && ee->to != c2 &&
                    ee->to != d && ee->to != e && ee->to != f &&
                    opt243(a, b, ee->to, c1, c, c2, d1, d, d2) &&
                    (!e || f || opt2433(a, b, ee->to, c1, c, c2,
                                        d1, d, d2, e1, e, e2)) &&
                    (!f || opt24333(a, b, ee->to, c1, c, c2, d1, d, d2,
                                    e1, e, e2, f1, f, f2)) &&
                    extra_extra_edge_24333_opt(a, b, ee->to, c1, c, c2,
                                               d1, d, d2, e1, e, e2,
                                               f1, f, f2, state))
                    break;
            if (!ee)
                return 0;
        }
    }
    for (i = 1; e && i <= 2; i++, swap(&e1, &e2)) {
        state++;
        if (e1 == a || e1 == b || e1 == c1 || e1 == c2 ||
            e1 == d1 || e1 == d2 || e1 == f1 || e1 == f2)
            continue;
        for (ee = e1->first_edge; ee; ee = ee->next)
            if (ee->to != e && ee->to != e2 &&
                ee->to != c && ee->to != d && ee->to != f && 
                (!f ? opt2433(a, b, ee->to, e1, e, e2, c1, c, c2, d1, d, d2)
                    : opt24333(a, b, ee->to, e1, e, e2, c1, c, c2,
                               d1, d, d2, f1, f, f2)) &&
                extra_extra_edge_24333_opt(a, b, ee->to, e1, e, e2,
                                           c1, c, c2, d1, d, d2, 
                                           f1, f, f2, state))
                break;
        if (!ee)
            return 0;
    }
    for (i = 1; f && i <= 2; i++, swap(&f1, &f2)) {
        state++;
        if (f1 == a || f1 == b || f1 == c1 || f1 == c2 ||
            f1 == d1 || f1 == d2 || f1 == e1 || f1 == e2)
            continue;
        for (ee = f1->first_edge; ee; ee = ee->next)
            if (ee->to != f && ee->to != f2 &&
                ee->to != c && ee->to != d && ee->to != e &&
                opt24333(a, b, ee->to, f1, f, f2, c1, c, c2,
                         d1, d, d2, e1, e, e2) &&
                extra_extra_edge_24333_opt(a, b, ee->to, f1, f, f2,
                                           c1, c, c2, d1, d, d2, 
                                           e1, e, e2, state))
                break;
        if (!ee)
            return 0;
    }
    for (i = 1; i <= 2; i++, swap(&a, &b)) {
        state++;
        if (b == c1 || b == c2 || b == d1 || b == d2 ||
            b == e1 || b == e2 || b == f1 || b == f2)
            continue;
        for (ee = b->first_edge; ee; ee = ee->next)
            if (ee->to != a &&
                ee->to != c && ee->to != d && ee->to != e && ee->to != f &&
                opt333(a, b, ee->to, c1, c, c2, d1, d, d2) &&
                (!e || f || opt3333(a, b, ee->to, c1, c, c2,
                                    d1, d, d2, e1, e, e2)) &&
                (!f || opt33333(a, b, ee->to, c1, c, c2, d1, d, d2,
                                e1, e, e2, f1, f, f2)) &&
                extra_extra_edge_33333_opt(a, b, ee->to, c1, c, c2, d1, d, d2,
                                           e1, e, e2, f1, f, f2, state))
                break;
        if (!ee)
            return 0;
    }
    return 1;
}

int extra_extra_edge_24333_opt(node * a, node * b,
                               node * c1, node * c2, node * c3, node * c4,
                               node * d1, node * d2, node * d3,
                               node * e1, node * e2, node * e3,
                               node * f1, node * f2, node * f3,
                               int state)
{
    int i, my_state = 0;
    edge *ee;

    if (extra_edges < 2)
        return 1;
    for (i = 1; i <= 2; i++, swap(&c1, &c4), swap(&c2, &c3)) {
        if (++my_state < state)
            continue;
        if (c1 == a || c1 == b || c1 == d1 || c1 == d3 ||
            c1 == e1 || c1 == e3 || c1 == f1 || c1 == f3)
            continue;
        for (ee = c1->first_edge; ee; ee = ee->next)
            if (ee->to != c2 && ee->to != c3 && ee->to != c4 &&
                ee->to != d2 && ee->to != e2 && ee->to != f2 &&
                opt253(a, b, ee->to, c1, c2, c3, c4, d1, d2, d3) &&
                (!e2 || f2 || opt2533(a, b, ee->to, c1, c2, c3, c4,
                                      d1, d2, d3, e1, e2, e3)) &&
                (!f2 || opt25333(a, b, ee->to, c1, c2, c3, c4, d1, d2, d3,
                                 e1, e2, e3, f1, f2, f3)))
                break;
        if (!ee)
            return 0;
    }
    for (i = 1; i <= 2; i++, swap(&d1, &d3)) {
        if (++my_state < state)
            continue;
        if (d1 == a || d1 == b || d1 == c1 || d1 == c4 ||
            d1 == e1 || d1 == e3 || d1 == f1 || d1 == f3)
            continue;
        for (ee = d1->first_edge; ee; ee = ee->next)
            if (ee->to != d2 && ee->to != d3 &&
                ee->to != c2 && ee->to != c3 && 
                ee->to != e2 && ee->to != f2 &&
                opt244(a, b, ee->to, d1, d2, d3, c1, c2, c3, c4) &&
                (!e2 || f2 || opt2443(a, b, ee->to, d1, d2, d3,
                                      c1, c2, c3, c4, e1, e2, e3)) &&
                (!f2 || opt24433(a, b, ee->to, d1, d2, d3,
                                 c1, c2, c3, c4, e1, e2, e3, f1, f2, f3)))
                break;
        if (!ee)
            return 0;
    }
    for (i = 1; e2 && i <= 2; i++, swap(&e1, &e3)) {
        if (++my_state < state)
            continue;
        if (e1 == a || e1 == b || e1 == c1 || e1 == c4 ||
            e1 == d1 || e1 == d3 || e1 == f1 || e1 == f3)
            continue;
        for (ee = e1->first_edge; ee; ee = ee->next)
            if (ee->to != e2 && ee->to != e3 &&
                ee->to != c2 && ee->to != c3 &&
                ee->to != d2 && ee->to != f2 &&
                opt443(ee->to, e1, e2, e3, c1, c2, c3, c4, d1, d2, d3) &&
                opt244(a, b, ee->to, e1, e2, e3, c1, c2, c3, c4) &&
                (!f2 || opt24433(a, b, ee->to, e1, e2, e3, c1, c2, c3, c4,
                                 d1, d2, d3, f1, f2, f3)))
                break;
        if (!ee)
            return 0;
    }
    for (i = 1; f2 && i <= 2; i++, swap(&f1, &f3)) {
        if (++my_state < state)
            continue;
        if (f1 == a || f1 == b || f1 == c1 || f1 == c4 ||
            f1 == d1 || f1 == d3 || f1 == e1 || f1 == e3)
            continue;
        for (ee = f1->first_edge; ee; ee = ee->next)
            if (ee->to != f2 && ee->to != f3 &&
                ee->to != c2 && ee->to != c3 &&
                ee->to != d2 && ee->to != e2 &&
                opt24433(a, b, ee->to, f1, f2, f3, c1, c2, c3, c4, 
                         d1, d2, d3, e1, e2, e3))
                break;
        if (!ee)
            return 0;
    }
    for (i = 1; i <= 2; i++, swap(&a, &b)) {
        if (++my_state < state)
            continue;
        if (b == c1 || b == c4 || b == d1 || b == d3 ||
            b == e1 || b == e3 || b == f1 || b == f3)
            continue;
        for (ee = b->first_edge; ee; ee = ee->next)
            if (ee->to != a &&
                ee->to != c2 && ee->to != c3 &&
                ee->to != d2 && ee->to != e2 && ee->to != f2 &&
                opt343(a, b, ee->to, c1, c2, c3, c4, d1, d2, d3) &&
                (!e2 || f2 || opt3433(a, b, ee->to, c1, c2, c3, c4,
                                      d1, d2, d3, e1, e2, e3)) &&
                (!f2 || opt34333(a, b, ee->to, c1, c2, c3, c4, d1, d2, d3,
                                 e1, e2, e3, f1, f2, f3)))
                break;
        if (!ee)
            return 0;
    }
    return 1;
}

int extra_extra_edge_33333_opt(node * a1, node * a2, node * a3,
                               node * c1, node * c2, node * c3,
                               node * d1, node * d2, node * d3,
                               node * e1, node * e2, node * e3,
                               node * f1, node * f2, node * f3,
                               int state)
{
    int i, my_state = 0;
    edge *ee;

    if (extra_edges < 2)
        return 1;
    for (i = 1; i <= 2; i++, swap(&c1, &c3)) {
        if (++my_state < state)
            continue;
        if (c1 == a1 || c1 == a3 || c1 == d1 || c1 == d3 ||
            c1 == e1 || c1 == e3 || c1 == f1 || c1 == f3)
            continue;
        for (ee = c1->first_edge; ee; ee = ee->next) {
            if (ee->to != c2 && ee->to != c3 &&
                ee->to != a2 && ee->to != d2 &&
                ee->to != e2 && ee->to != f2 &&
                opt343(a1, a2, a3, ee->to, c1, c2, c3, d1, d2, d3) &&
                (!e2 || f2 || opt3433(a1, a2, a3, ee->to, c1, c2, c3, 
                                      d1, d2, d3, e1, e2, e3)) &&
                (!f2 || opt34333(a1, a2, a3, ee->to, c1, c2, c3, 
                                 d1, d2, d3, e1, e2, e3, f1, f2, f3)))
                break;
        }
        if (!ee)
            return 0;
    }
    for (i = 1; i <= 2; i++, swap(&d1, &d3)) {
        if (++my_state < state)
            continue;
        if (d1 == a1 || d1 == a3 || d1 == c1 || d1 == c3 ||
            d1 == e1 || d1 == e3 || d1 == f1 || d1 == f3)
            continue;
        for (ee = d1->first_edge; ee; ee = ee->next) {
            if (ee->to != d2 && ee->to != d3 &&
                ee->to != a2 && ee->to != c2 &&
                ee->to != e2 && ee->to != f2 &&
                opt343(a1, a2, a3, ee->to, d1, d2, d3, c1, c2, c3) &&
                (!e2 || f2 || opt3433(a1, a2, a3, ee->to, d1, d2, d3, 
                                      c1, c2, c3, e1, e2, e3)) &&
                (!f2 || opt34333(a1, a2, a3, ee->to, d1, d2, d3, 
                                 c1, c2, c3, e1, e2, e3, f1, f2, f3)))
                break;
        }
        if (!ee)
            return 0;
    }
    for (i = 1; e2 && i <= 2; i++, swap(&e1, &e3)) {
        if (++my_state < state)
            continue;
        if (e1 == a1 || e1 == a3 || e1 == c1 || e1 == c3 ||
            e1 == d1 || e1 == d3 || e1 == f1 || e1 == f3)
            continue;
        for (ee = e1->first_edge; ee; ee = ee->next) {
            if (ee->to != e2 && ee->to != e3 &&
                ee->to != a2 && ee->to != c2 &&
                ee->to != d2 && ee->to != f2 &&
                opt3433(a1, a2, a3, ee->to, e1, e2, e3, 
                        c1, c2, c3, d1, d2, d3) &&
                (!f2 || opt34333(a1, a2, a3, ee->to, e1, e2, e3, 
                                 c1, c2, c3, d1, d2, d3, f1, f2, f3)))
                break;
        }
        if (!ee)
            return 0;
    }
    for (i = 1; f2 && i <= 2; i++, swap(&f1, &f3)) {
        if (++my_state < state)
            continue;
        if (f1 == a1 || f1 == a3 || f1 == c1 || f1 == c3 ||
            f1 == d1 || f1 == d3 || f1 == e1 || f1 == e3)
            continue;
        for (ee = f1->first_edge; ee; ee = ee->next) {
            if (ee->to != e2 && ee->to != f3 &&
                ee->to != a2 && ee->to != c2 &&
                ee->to != d2 && ee->to != e2 &&
                opt34333(a1, a2, a3, ee->to, f1, f2, f3, 
                         c1, c2, c3, d1, d2, d3, e1, e2, e3))
                break;
        }
        if (!ee)
            return 0;
    }
    for (i = 1; i <= 2; i++, swap(&a1, &a3)) {
        if (++my_state < state)
            continue;
        if (a1 == c1 || a1 == c3 || a1 == d1 || a1 == d3 ||
            a1 == e1 || a1 == e3 || a1 == f1 || a1 == f3)
            continue;
        for (ee = a1->first_edge; ee; ee = ee->next) {
            if (ee->to != a2 && ee->to != a3 &&
                ee->to != c2 && ee->to != d2 &&
                ee->to != e2 && ee->to != f2 &&
                opt343(c1, c2, c3, ee->to, a1, a2, a3, d1, d2, d3) &&
                (!e2 || opt3433(c1, c2, c3, ee->to, a1, a2, a3, 
                                d1, d2, d3, e1, e2, e3)) &&
                (!f2 || opt34333(c1, c2, c3, ee->to, a1, a2, a3, 
                                 d1, d2, d3, e1, e2, e3, f1, f2, f3)))
                break;
        }
        if (!ee)
            return 0;
    }
    return 1;
}

static void firstset(int k, aSet S);
static void nextset(int k, aSet S);
static int *Sval (int k, aSet S, int t_indx, int * V, int * b);

/* Bellman-Held-Karp dynamic programming algorithm addapted from
 * Applegate et al.: "A Practical Guide to Discrete Optimization", 2015.
 */

int path_opt(node ** p, int n, int * fix)
{
    int i, j, k, t, minv, v, *valbase, result = 0;
    aSet S, S_minus_t;
    int cost, *dist, *b, *V;

    if (n > MAX_N)
        return 1;
    n--;
    assert((V = (int *) malloc((n - 1) * (1 << (n - 2)) * sizeof(int))));
    assert((b = (int *) malloc(n * sizeof(int))));
    assert((dist = (int *) malloc(n * n * sizeof(int))));
    cost = C(p[0], p[1]);
    for (i = 1; i < n; i++) {
        cost += dist[i * n + i - 1] = dist[(i - 1) * n + i] =
            fix[i] ? INT_MIN / n : C(p[i], p[i + 1]);
        for (j = i - 2; j >= 0; j--)
            dist[i * n + j] = dist[j * n + i] = C(p[i + 1], p[j + 1]);
    }
    b[0] = 0;
    for (i = 1; i < n; i++)
        b[i] = b[i - 1] + (i - 1) * binomial[n - 1][i - 1];
    for (firstset(1, S); S[0] < n - 1; nextset(1, S))
        *Sval(1, S, 0, V, b) = C(p[S[0] + 1], p[0]);
    for (k = 2; k < n; k++) {
        for (firstset(k, S); S[k - 1] < n - 1; nextset(k, S)) {
            for (t = 1; t < k; t++)
                S_minus_t[t - 1] = S[t];
            for (t = 0; t < k; t++) {
                valbase = Sval(k - 1, S_minus_t, 0, V, b);
                minv = INT_MAX;
                for (j = 0; j < k - 1; j++) {
                    v = valbase[j] + dist[S[t] * n + S_minus_t[j]];
                    if (k == n - 1) {
                        if (v + dist[t * n + n - 1] < cost)
                            goto CLEANUP;
                    } else if (v < minv)
                        minv = v;
                }
                *Sval(k, S, t, V, b) = minv;
                S_minus_t[t] = S[t];
            }
        }
    }
    result = 1;
CLEANUP:
    free(V);
    free(b);
    free(dist);
    return result;
}

static void firstset(int k, aSet S)
{
    int i;
    for (i = 0; i < k; i++)
        S[i] = i;
}

static void nextset(int k, aSet S)
{
    int i;
    for (i = 0; i < k - 1 && S[i] + 1 == S[i + 1]; i++)
        S[i] = i;
    S[i]++;
}

static int *Sval(int k, aSet S, int t_indx, int * V, int * b)
{
    unsigned long int loc = 0;
    int i;
    for (i = 0; i < k; i++)
        loc += binomial[S[i]][i + 1];
    return &V[b[k] + k * loc + t_indx];
}

// Simple implementation of Bellman-Held-Karp dynamic programming algorithm 

/*
int path_opt(node ** p, int n, int * fix)
{
    if (n > MAX_N)
        return 1;
    int i, j, S;
    int cost = 0, dist[n][n + 1], f[n][1 << n];

    n--;
    for (i = 1; i <= n; i++) {
        for (j = 0; j < i; j++)
            dist[i][j] = dist[j][i] = 
                i == j + 1 && fix[j] ? INT_MIN / n : C(p[i], p[j]);
        cost += dist[i - 1][i];
    }
    for (i = 0; i < n; i++)
        for (S = 1; S < 1 << n; S++)
            f[i][S] = INT_MAX / 2;
    f[0][1] = 0;
    for (S = 3; S < 1 << n; S += 2)
        for (i = 1; i < n; i++)
            if ((S & 1 << i) != 0)
                for (j = 0; j < n; j++)
                    if (j != i && (S & 1 << j) != 0)
                        f[i][S] = min(f[i][S],
                                      f[j][S ^ (1 << i)] + dist[j][i]);
    for (i = 1; i < n; i++)
        if (f[i][(1 << n) - 1] + dist[i][n] < cost)
            return 0;
    return 1;
}
*/

int has_cycle222(node * a, node * b,
                 node * c1, node * c2,
                 node * d1, node * d2)
{
    node *p[6] = { a, b, c1, c2, d1, d2 };
    return has_cycle(6, p);
}

int has_cycle2222(node * a, node * b,
                  node * c1, node * c2,
                  node * d1, node * d2,
                  node * e1, node * e2)
{
    node *p[8] = { a, b, c1, c2, d1, d2, e1, e2 };
    return has_cycle(8, p);
}

int has_cycle22222(node * a, node * b,
                   node * c1, node * c2,
                   node * d1, node * d2,
                   node * e1, node * e2,
                   node * f1, node * f2)
{
    node *p[10] = { a, b, c1, c2, d1, d2, e1, e2, f1, f2 };
    return has_cycle(10, p);
}

int has_cycle222222(node * a, node * b,
                    node * c1, node * c2,
                    node * d1, node * d2,
                    node * e1, node * e2,
                    node * f1, node * f2,
                    node * g1, node * g2)
{
    node *p[12] = { a, b, c1, c2, d1, d2, e1, e2, f1, f2, g1, g2 };
    return has_cycle(12, p);
}

double get_time()
{
    struct rusage ru;
    getrusage(RUSAGE_SELF, &ru);
    return ru.ru_utime.tv_sec + ru.ru_utime.tv_usec / 1000000.0;
}

int euc2d(node * a, node * b)
{
    double xd = a->x - b->x, yd = a->y - b->y;
    return (int) (sqrt(xd * xd + yd * yd) + 0.5);
}

int euc3d(node * a, node * b)
{
    double xd = a->x - b->x, yd = a->y - b->y, zd = a->z - b->z;
    return (int) (sqrt(xd * xd + yd * yd + zd * zd) + 0.5);
}

int ceil2d(node * a, node * b)
{
    double xd = a->x - b->x, yd = a->y - b->y;
    return (int) ceil(sqrt(xd * xd + yd * yd));
}

#undef M_PI
#define M_PI 3.14159265358979323846264
#define M_RRR 6378388.0

int geom(node * a, node * b)
{
    double lati = M_PI * (a->x / 180.0);
    double latj = M_PI * (b->x / 180.0);
    double longi = M_PI * (a->y / 180.0);
    double longj = M_PI * (b->y / 180.0);
    double q1 = cos(latj) * sin(longi - longj);
    double q3 = sin((longi - longj) / 2.0);
    double q4 = cos((longi - longj) / 2.0);
    double q2 = sin(lati + latj) * q3 * q3 - sin(lati - latj) * q4 * q4;
    double q5 = cos(lati - latj) * q4 * q4 - cos(lati + latj) * q3 * q3;
    return (int) (M_RRR * atan2(sqrt(q1 * q1 + q2 * q2), q5) + 1.0);
}

int parseargs(int argc, char *argv[])
 {
    int c;
    while ((c = getopt(argc, argv, "EJe:n:o:p:qrst:T:x:w:")) != -1) {
        switch(c) {
        case 'E':
            print_E = 1;
            break;
        case 'J':
            Jonker_Volgenant = 1;
            break;
        case 'e':
            extra_edges = atoi(optarg);
            if (extra_edges < 0)
                extra_edges = 0;
            if (extra_edges > 2)
                extra_edges = 2;
            break;
        case 'n':
            extra_nodes = atoi(optarg);
            if (extra_nodes < 0)
                extra_nodes = 0;
            if (extra_nodes > 3)
                extra_nodes = 3;
            break;
        case 'o':
            output_file_name = optarg;
            break;
        case 'p':
            max_cd_set_size = atoi(optarg);
            if (max_cd_set_size <= 0)
                max_cd_set_size = 1;
            break;
        case 'q':
            quicker_cd_search = 1;
            break;
        case 'r':
            repeated = 1;
            break;
        case 's':
            strong_3_opt = 1;
            break;
        case 't':
            tour_file_name = optarg;
            break;
        case 'T':
            tsp_file_name = optarg;
            break;
        case 'w':
            max_cd_count = atoi(optarg);
            break;
        case 'x':
            ab_stretch = atof(optarg);
            break;
        case '?':
        default:
            usage (argv[0]);
            return 1;
        }
    }
    if (optind < argc)
        edge_file_name = argv[optind++];
    if (!tsp_file_name || !edge_file_name) {
        if (!tsp_file_name)
            printf("*** Missing tsp_file\n");
        if (!edge_file_name)
            printf("*** Missing edge_file\n");
        usage(argv[0]);
        return 1;
    }
    return 0;
}

void usage(char *f)
{
    fprintf(stderr, "usage: %s [- see below -] -T tsp_file edge_file\n", f);
    fprintf(stderr, "   -E    print an E for each eliminated edge\n");
    fprintf(stderr, "   -J    use Jonker-Volgenant elimination\n");
    fprintf(stderr, "   -e #  use extra edges, range 0 to 2, default 0\n");
    fprintf(stderr, "   -n #  use extra nodes, range 0 to 3, default 0\n");
    fprintf(stderr, "   -o f  dump the remaining edges to file f\n");
    fprintf(stderr, "   -p #  maximum potential points, default 10\n");
    fprintf(stderr, "   -q #  quicker search for potential points\n");
    fprintf(stderr, "   -r    repeat until no more eliminations can be made\n");
    fprintf(stderr, "   -s    perform strong close point elimination\n");
    fprintf(stderr, "   -t f  tour file (no tour edges will be eliminated)\n");
    fprintf(stderr, "   -T f  TSPLIB file\n");
    fprintf(stderr, "   -w #  maximum witness edges for each edge, default 10\n");
    fprintf(stderr, "   -x d  double, stretch factor, d in [0.5;2.5] is best, default 1.0\n");
}

#endif
