#
#  Makefile for PRESTO:  Pulsar Search Software
#     for Unix (and hopefully MPI soon)
#           by Scott M. Ransom
#            V 0.92, 23 Jun 99
#
DATE = $(shell date +%d%b%y)

tar:  squeaky package

package:
	@echo "Tarring..."
	@echo ""
	@cd ..; tar -cf presto$(DATE).tar presto
	@echo "Gzipping..."
	@gzip -9 ../presto$(DATE).tar
	@mv ../presto$(DATE).tar.gz .
	@cp presto$(DATE).tar.gz /pulsar/
	@cp /home/ransom/programming/fftwresults/fftw_wisdom.txt ./lib/
	@echo ""
	@echo "Done."
	@echo ""

squeaky:
	@echo "Cleaning house:"
	@echo ""
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
	cd src ; find .. -name "*.[ch]" -print | etags -
	@echo ""
	@echo "All is now squeaky clean and neat."
	@echo ""














