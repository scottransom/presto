SWIG          = swig 
SWIGOPT       = -make_default -python -shadow

prestowrap:
	$(SWIG) $(SWIGOPT) presto.i
	cp presto.py ../presto_src/prestoswig.py
	cp presto_wrap.c ../presto_src/
clean:
	rm -f *~ *.o *.so *.pyc *.py *.c
