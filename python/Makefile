build:
	python setup.py install --home=${PRESTO}

fftfit:
	cd fftfit_src ; f2py -c fftfit.pyf *.f
	cd fftfit_src ; cp fftfit.so ${PRESTO}/lib/python

clean:
	rm -rf build
	rm -f *~ *.o *.so *.pyc
