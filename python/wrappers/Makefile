SWIG          = swig3.0
SWIGOPT       = -make_default -python -shadow

prestowrap:
	$(SWIG) $(SWIGOPT) presto.i
	cp presto.py ../presto_src/prestoswig.py
	cp presto_wrap.c ../presto_src/presto_wrap.c
clean:
	rm -f *~ *.o *.so *.pyc *.py *.c
