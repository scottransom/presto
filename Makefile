#
#  Makefile for PRESTO:  Pulsar Search Software
#     for Unix (and hopefully MPI soon)
#           by Scott M. Ransom
#            V 0.94, 16 Nov 01
#
DATE = $(shell date +%d%b%y)

tar:  squeaky tags package

tags:
	cd src ; find . -name "*.[ch]" -print | etags -
	cd src ; find ../include -name "*.[ch]" -print | etags -a -

package:
	cd ..; tar -cf presto$(DATE).tar presto
	gzip -9 ../presto$(DATE).tar
	mv ../presto$(DATE).tar.gz .

squeaky:
	rm -f *~ presto*.tar.gz *#
	cd lib ; rm -f *~ *.o *.a *.so *#
	cd bin ; rm -f *~ *.o *.a *.so *.dat *.fft *.inf *#
	cd src ; rm -f *~ *.o *.a *.so *.dat *.fft *.inf *#
	cd include ; rm -f *~ *.o *.a *.so *.dat *.fft *.inf *#
	cd python ; rm -f *~ *.o *.a *.so *.dat *.fft *.inf *#
	cd docs ; rm -f *~ *.o *.a *.so *.dat *.fft *.inf *#
	cd tests ; rm -f *~ *.o *.a *.so *.dat *.fft *.inf *#
	cd docs ; rm -f *~ *.o *.a *.so *.dat *.fft *.inf *#
	cd clig ; rm -f *~ *.o *.a *.so *.dat *.fft *.inf *#
	cd MPI ; rm -f *~ *.o *.a *.so *.dat *.fft *.inf *#
	cd oldsource ; rm -f *~ *.o *.a *.so *.dat *.fft *.inf *#
	cd python ; make clean
	cd src ; make squeaky














