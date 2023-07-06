CC = gcc
CFLAGS = -O3 -Wall -Wno-unused-result
 
INSTANCE = pr299
OPTIONS = -o $(INSTANCE).$$.edg -E -J -s -r -x1.5 -p10 -w10 -n1 -e1

all:
	@$(MAKE) kh-elim
	@$(MAKE) kh-elim_omp

kh-elim: kh-elim.c kh-elim_common.c kh-elim.h
	$(CC) $(CFLAGS) -o kh-elim kh-elim.c -lm

kh-elim_omp: kh-elim_omp.c kh-elim_common.c kh-elim.h
	$(CC) -fopenmp $(CFLAGS) -o kh-elim_omp kh-elim_omp.c -lm

clean:
	/bin/rm -f kh-elim kh-elim_omp
	/bin/rm -f  *~

check:
	$(MAKE) all
	@echo
#	kh-elim $(OPTIONS) -T $(INSTANCE).tsp $(INSTANCE).edg
	@echo
	time -p kh-elim_omp $(OPTIONS) -T $(INSTANCE).tsp $(INSTANCE).edg
