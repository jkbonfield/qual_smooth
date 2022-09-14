all: qual_smooth qual_train

# NB: HTSLIB and HTSCODECS dirs are local source trees
HTSLIB=/nfs/users/nfs_j/jkb/work/samtools_master/htslib
HTSCODECS=/nfs/users/nfs_j/jkb/work/samtools_master/htscodecs

CFLAGS=-g -O3 -Wall
CPPFLAGS=-I$(HTSLIB) -I$(HTSCODECS)
LDFLAGS=-L$(HTSLIB) -L$(HTSCODECS)/htscodecs/.libs -Wl,-R$(HTSLIB) -lhts -lhtscodecs -lm -pthread
CC=gcc

.c.o:
	$(CC) -c $< -o $@ $(CFLAGS) $(CPPFLAGS)

OBJS = smooth.o entropy.o
qual_smooth: $(OBJS)
	$(CC) $(CFLAGS) $(CPPFLAGS) $(OBJS) -o $@ $(LDFLAGS)

qual_train: train.o
	$(CC) $(CFLAGS) $(CPPFLAGS) $< -o $@ $(LDFLAGS)

clean:
	rm qual_smooth qual_train *.o
