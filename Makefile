#
#  Makefile for PRESTO:  Pulsar Search Software
#     for Unix (and hopefully MPI soon)
#           by Scott M. Ransom
#
DATE = $(shell date +%d%b%y)

tar:  squeaky tags package

tags:
	cd src ; find . -name "*.[ch]" -print | etags -
	cd src ; find ../include -name "*.[ch]" -print | etags -a -

package:
	cd ..; tar --exclude=.svn -cvf presto$(DATE).tar presto
	gzip -9 ../presto$(DATE).tar
	mv ../presto$(DATE).tar.gz .

squeaky:
	rm -f *~ presto*.tar.gz *#
	find . -name "*.[oa]" | xargs rm -f
	find . -name "*.so"   | xargs rm -f
	find . -name "*.dat"  | xargs rm -f
	find . -name "*.fft"  | xargs rm -f
	find . -name "*.inf"  | xargs rm -f
	find . -name "*[~#]"  | xargs rm -f
	cd python ; make clean
	cd src ; make squeaky
