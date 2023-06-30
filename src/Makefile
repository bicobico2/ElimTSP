CC=gcc
CFLAGS = -O3 -g -Wall -Wno-comment -Wno-unused-result -Wstrict-prototypes -Wpointer-arith -Wnested-externs
LDFLAGS = -O3 -g 

ELIMOBJS=elim.o cdtest.o path.o pair.o improve.o graph.o swap.o tree.o
CCUTILOBJS= cc_heldkarp.o cc_edgelen.o cc_getdata.o cc_util.o

elim: elimmain.o $(ELIMOBJS) $(CCUTILOBJS) 
	$(CC) $(CLAGS) -o elim elimmain.o $(ELIMOBJS) $(CCUTILOBJS) -lm

verify: verify.o $(ELIMOBJS) $(CCUTILOBJS)
	$(CC) $(CFLAGS) -o verify verify.o $(ELIMOBJS)  $(CCUTILOBJS) -lm

clean:
	-rm -f *.o elim verify

elimmain.o:  elimmain.c  elim.h  ccutil.h
elim.o:      elim.c      elim.h  ccutil.h
cdtest.o:    cdtest.c    elim.h  ccutil.h
path.o:      path.c      elim.h  ccutil.h
pair.o:      pair.c      elim.h  ccutil.h
improve.o:   improve.c   elim.h  ccutil.h
graph.o:     graph.c     elim.h  ccutil.h
swap.o:      swap.c      elim.h  ccutil.h
tree.o:      tree.c      elim.h  ccutil.h
verify.o:    verify.c    elim.h  ccutil.h

cc_edgelen.o:   cc_edgelen.c   ccutil.h
cc_getdata.o:   cc_getdata.c   ccutil.h
cc_heldkarp.o:  cc_heldkarp.c  ccutil.h
cc_util.o:      cc_util.c      ccutil.h
