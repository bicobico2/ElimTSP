# ElimTSP
Eliminating edges in the (symmetric) traveling salesman problem.

The codes are based on the paper "Local elimination in the traveling salesman problem" by William Cook, Keld Helsgaun,
Stefan Hougardy, and Rasmus T. Schroeder.

### Executables

**elim:**  Runs the main edge elimination code. Takes as input a TSPLIB file (to specify
   the edge lengths in the complete graph) and a Concorde-style edge file
   for the initial sparse graph.

**verify:** Verifies the output of elim. Takes as input a Hamilton-Tutte tree
   file, a TSPLIB file, and a edge file for the initial graph.

**kh-elim:** Run an experimental elimination code. Takes as input a TSPLIB file
   and an edge file.

Running elim, verify, or kh-elim with no arguments will display a short description
of the run time options.

### Example

elim -z4 -H pr299.ht -E -T pr299.tsp pr299.edg

The -z4 option specifies that Tutte moves up to depth 4 should be considered.
The option -H pr299.ht directs the code to write a Hamilton-Tutte tree to the
file pr299.ht. The -E options directs the code to print an "E" whenever an
edge is eliminated.

verify -e pr299.ver.edg -V pr299.ht -T pr299.tsp pr299.edg

The option -e pr299.ver.edg directs the code to write the verified edge
list to pr299.ver.edg. 
The -V pr299.ht specifies the Hamilton-Tutte tree file pr299.ht (that was
generated by the run of elim).

### Notes 
 1. The files cc_edgelen.c, cc_getdata.c, cc_heldkarp.c, and cc_util.c
    contain functions from the Concorde TSP library.
 2. The results in Table 4 of the paper "Local elimination in the 
    traveling salesman problem" were produced using "elim -A" to loop through
    a sequence of elimination levels.
 3. The results in Table 7 of the paper "Local elimination in the 
    traveling salesman problem" were produced using "kh-elim -Jq" (for a
    single core) and "kh-elim_omp -Jq" (for multiple cores).
