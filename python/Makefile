build:
	python setup.py install --home=${PRESTO}
	cd fftfit_src ; f2py -c fftfit.pyf *.f
	cd fftfit_src ; cp fftfit.so ${PRESTO}/lib/python
	python fftfit_src/test_fftfit.py

clean:
	rm -rf build
	rm -f *~ *.o *.so *.pyc
