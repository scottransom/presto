#!/usr/bin/env python
import mygetopt, string, sys, os, presto, Numeric

if (len(sys.argv)==1):
    print """
usage:  preppaste [whatever options you send to prepdata] rawfiles
  This program takes Parkes Multibeam rawfiles and pastes them
  together to create a single barycentered and de-dispersed
  datafile.  This is useful if your OS cannot handle files greater
  than 2GB in length or if your DLT has bad blocks that had to be
  skipped during the reading of the tape.
  
  Extra options:
    '-numchan':   The number of channels in the data.
                  Default = 256 for PKMB
    '-padval':    The value to use when filling gaps or padding
                  the final data set.
                  Default = 0.5 * (numchan - 1) for PKMB
    '-blocklen':  The number of bytes per data block.
                  Default = 49792 for PKMB

  Note:  The rawfiles need to be able to be sorted into the
         correct order!  (i.e. alphabetical == chronological)"""
    sys.exit(0)

def pkmb_prepdata(optlist, file, filenum):
    # Runs prepdata and returns a tuple containg number of points
    # written, time per point, and starting time (MJD)
    outfile = optlist['o']+`filenum`
    command = 'prepdata -pkmb -nobary -dm '+optlist['dm']+' -o '+\
              outfile+' '+file
    print ""
    print command
    print ""
    os.system(command)
    inf = presto.read_inffile(outfile)
    return (int(inf.N), inf.dt, inf.mjd_i+inf.mjd_f)

# Get the command line args

options = ['o=', 'pkmb', 'ebpp', 'pad0', 'padavg',
           'numout=', 'nobary', 'DE405', 'dm=', 'padval=']
optlist, files = mygetopt.getopt(sys.argv[1:], options)
files.sort()

if optlist.has_key('pkmb'):
    datatype = 'pkmb'
    numout = int(optlist.get('numout', 0))
    blocklen = int(optlist.get('blocklen', 49792))
    numchan = int(optlist.get('numchan', 256))
    padval = float(optlist.get('padval', 0.5 * (numchan - 1.0)))
    ptsperblock = (blocklen - 640) * 8 / numchan

# Prep the individual data files (makes topocentric data)

print "\n\n**** Making topocentric data files...\n\n"

ns = []
epochs = []
for filenum in xrange(len(files)):
    file = files[filenum]
    if (datatype=='pkmb'):
        (nout, dt, epoch) = pkmb_prepdata(optlist, file, filenum)
    ns.append(nout)
    epochs.append(epoch)
if (numout):
    epochs.append(epochs[0] + numout * dt / 86400.0)
    ns.append(0)
ns = Numeric.asarray(ns)
epochs = Numeric.asarray(epochs)

# Calculate the amount of padding to add after each file

binsneeded = ((epochs[1:] - epochs[:-1]) * 86400.0 /
              dt + 0.5).astype('i')
padbins = binsneeded - ns[:-1]

# Add padding to the data topocentric files we just wrote

print "\n\n**** Adding padding...\n\n"

for filenum in xrange(len(padbins)):
    outfile = optlist['o']+`filenum`+'.dat'
    if (padbins[filenum] > 0):
	command = 'patchdata '+`padbins[filenum]`+' '+`padval`+' >> '+outfile
	print ""
	print command
	print ""
	os.system(command)

# Cat the files together and remove the temp files
# Note:  need to add ability to check length of data written if
#        the user wants to write an amount that is less than
#        the data written in the topocentric steps...
#        This might not be needed since we will run prepdata
#        to barycenter the data.

print "\n\n**** Joining files...\n\n"

outfile = optlist['o']+'0.dat'
for filenum in xrange(1, len(files)):
    infile = optlist['o']+`filenum`+'.dat'
    command = 'cat '+infile+' >> '+outfile
    print ""
    print command
    print ""
    os.system(command)
    os.remove(infile)
    os.remove(optlist['o']+`filenum`+'.inf')

# Adjust infofile for the pasted data set

inf = presto.read_inffile(optlist['o']+'0')
inf.N = os.stat(outfile)[6] / 4
presto.writeinf(inf)

# Run prepdata on the big file to barycenter it

print "\n\n**** Barycentering...\n\n"

if (numout):
    command = ('prepdata -numout '+`numout`+' -o '+
               optlist['o']+' '+outfile)
else:
    command = ('prepdata -o '+optlist['o']+' '+outfile)
print ""
print command
print ""
os.system(command)
os.remove(outfile)
os.remove(optlist['o']+'0.inf')

print "\n\n**** Done.\n\n"
