PYTHONBASE    = /usr
DATABASE      = "$(PRESTO)/lib"
INSTALLDIR    = /home/ransom/python
INCLUDE       = -I$(PYTHONBASE)/include/python2.2 -I$(PYTHONBASE)/include/python2.2/Numeric -I../include
#INCLUDE       = -IPYTHONBASE/include/python2.2 -I../include -I/home/ransom/include
CC            = cc
CWARN         = -Xa -vc
CFLAGS        = -g -O -fPIC $(INCLUDE) -DDATABASE='$(DATABASE)' -DUSEFFTW
#CFLAGS        = -g -O $(CWARN) $(INCLUDE) -DDATABASE='$(DATABASE)' -DUSEFFTW
LIBPATHS      = -L"$(PRESTO)/lib" -L/usr/local/lib
LIBS          = -lpresto -lfftw3f -lm
SWIG          = /usr/bin/swig 
SWIGOPT       = -make_default -python -shadow
PRESTOOBJS    = presto_wrap.o

lib: _presto.so
wrap: prestowrap

prestowrap:
	$(SWIG) $(SWIGOPT) presto.i
	mv presto.py presto.py.bak
	cat presto.py.bak presto_addins.py > presto.py

_presto.so: $(PRESTOOBJS)
	rm -f _presto.so
	$(CC) $(CFLAGS) $(LIBPATHS) -shared $(PRESTOOBJS) -o $@ $(LIBS)
#	ld -G $(LIBPATHS) $(PRESTOOBJS) $(LIBS) -o $@ 

update:
	cat presto.py.bak presto_addins.py > presto.py

install:
	cp _presto.so $(INSTALLDIR)
	cp presto.py $(INSTALLDIR)

clean:
	rm -f *~ *.o *.so *.pyc




