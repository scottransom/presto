PYTHONBASE    = /home/ransom/python2.1
DATABASE      = "$(PRESTO)/lib"
INSTALLDIR    = /home/ransom/python
INCLUDE       = -I$(PYTHONBASE)/include/python2.1 -I$(PYTHONBASE)/include/python2.1/Numeric -I../include
#INCLUDE       = -IPYTHONBASE/include/python2.1 -I../include -I/home/ransom/include
CC            = cc
CWARN         = -Xa -vc
CFLAGS        = -g -O -fpic $(INCLUDE) -DDATABASE='$(DATABASE)' -DUSEFFTW
#CFLAGS        = -g -O $(CWARN) $(INCLUDE) -DDATABASE='$(DATABASE)' -DUSEFFTW
LIBPATHS      = -L"$(PRESTO)/lib" -L/home/ransom/lib
LIBS          = -lpresto -lsfftw -lm
SWIG          = swig 
SWIGOPT       = -python -shadow
PRESTOOBJS    = presto_wrap.o

wrap: prestowrap
lib: prestoc.so

prestowrap:
	$(SWIG) $(SWIGOPT) presto.i
	mv presto.py presto.py.bak
	cat presto.py.bak presto_addins.py > presto.py

prestoc.so: $(PRESTOOBJS)
	rm -f prestoc.so
	$(CC) $(CFLAGS) $(LIBPATHS) -shared $(PRESTOOBJS) -o $@ $(LIBS)
#	ld -G $(LIBPATHS) $(PRESTOOBJS) $(LIBS) -o $@ 

update:
	cat presto.py.bak presto_addins.py > presto.py

install:
	cp prestoc.so $(INSTALLDIR)
	cp presto.py $(INSTALLDIR)

clean:
	rm -f *~ *.o *.so *.pyc




