SWIG          = swig
SWIGOPT       = -Wall -python -py3

prestowrap:
	$(SWIG) $(SWIGOPT) presto.i
	cp presto.py ../presto_src/prestoswig.py
	cp presto_wrap.c ../presto_src/presto_wrap.c
clean:
	rm -f *~ *.o *.so *.pyc *.py *.c
