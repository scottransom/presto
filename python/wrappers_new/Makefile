SWIG          = swig 
SWIGOPT       = -make_default -python -shadow

prestowrap:
	$(SWIG) $(SWIGOPT) presto.i
	cp newpresto.py ../presto_src/newprestoswig.py
	cp presto_wrap.c ../presto_src/newpresto_wrap.c
clean:
	rm -f *~ *.o *.so *.pyc *.py *.c
