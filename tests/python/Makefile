LIBS = "/home/ransom/cprogs/Working/presto/lib/"
INC = "/home/ransom/cprogs/Working/presto/include/"
PYINC = "/usr/local/include/python1.5/" 
CC = egcc
F77 = g77
CFLAGS = -I$(PYINC) -I$(INC) -I./ -Wall -ansi -pedantic -pipe -fPIC -O -g

OBJS = cvects.o /home/ransom/cprogs/Working/presto/src/vectors.o
WRAPS = cvects_wrap.o

all: wrap my_range.so

wrap: test_wrap.c
	swig -python test.i

libcvects.so: $(OBJS)
	$(CC) $(CFLAGS) -shared $(OBJS) -o $@

cvects.so: $(OBJS) $(WRAPS)
	$(CC) $(CFLAGS) -shared $(OBJS) $(WRAPS) -o $@

statsum.so: statsum.o statsum_wrap.o
	$(CC) $(CFLAGS) -shared statsum.o statsum_wrap.o -o $@

my_range.so: shorttest.o test_wrap.o
	$(CC) $(CFLAGS) -shared shorttest.o test_wrap.o -o $@

test_cvects: $(OBJS) test_cvects.o
	$(CC) $(CFLAGS) -o $@ $(OBJS) test_cvects.o -L./ -lcvects -lm

clean:
	rm -f *.o *~ *# *.so
	rm -f test_cvects cvects_wrap.o














