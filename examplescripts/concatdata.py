from __future__ import print_function
from future import standard_library
standard_library.install_aliases()
from builtins import range
# Binary floating point data file concatenation routine
#
# Written by Scott M. Ransom <ransom@cfa.harvard.edu>
# last revision: 1 Mar 99
#

usage = """

Usage:  concatdata outfile numpts padval infile1 infile2 ...
    This routine helps properly connect binary time-series.
    Mandatory arguments:
      outfile:  The intended name of the freshly concat-ed time-series.
       numpts:  The number of points for the new time-series to contain.
       padval:  The value to pad the time-series with.
      infiles:  A list (in order) of the short data sets to concatenate.
                Each file must have a '.dat' suffix and a corresponding
                '.inf' file describing each must exist.  There
                must be at least 2 files to concat.

"""

from sys import argv, exit
from string import atol, rfind
from presto import *
from infodata import *
from math import modf
from Numeric import *

def addtoMJD(daystoadd, MJDi, MJDf):
    (fdays, idays) = modf(daystoadd)
    MJDf = MJDf + fdays
    if (idays >= 0.0):
        MJDi = MJDi + int(idays + 1.0e-10)
    else:
        MJDi = MJDi + int(idays - 1.0e-10)
    if (MJDf >= 1.0):
        MJDf = MJDf - 1.0
        MJDi = MJDi + 1
    if (MJDf < 0.0):
        MJDf = MJDf + 1.0
        MJDi = MJDi - 1
    return (MJDi, MJDf)

def subtractMJDs(MJDi1, MJDf1, MJDi2, MJDf2):
    # return MJD1 - MJD2
    return MJDi1 - MJDi2 + MJDf1 - MJDf2
    
debug = 1
SBF = 1.0e-4   # Smallest bin fraction to worry about

# Show a usage statement if necessary

if (len(argv)<6):
    print(usage)
    exit(0)

# Get and check the arguments

print('')
print('   Binary Data Concatenation Routine')
print('      Written by Scott M. Ransom')
print('              1 Mar 99\n')

outfilenm = argv[1]
numpts = atol(argv[2])
padval = float(argv[3])
if (numpts < 0):
    print('numpts must be greater than 0.  Exiting.')
    print(usage)
    exit(-1)

print('Creating a %ld point file named \'%s\'.' % (numpts, outfilenm))
print('Using %f for each padding point.\n' % padval)

# Read the important data from the infofiles into lists

infile = []
file_data = []
file_startMJDi = []
file_startMJDf = []
file_endMJDi = []
file_endMJDf = []
file_N = []
padbins = []
print('The input files are:')
for index in range(len(argv)-4):
    infile.append(argv[index+4])
    infile[index] = infile[index][0:rfind(infile[index],'.')]

    # Get the info about the data file

    file_data.append(infodata(infile[index]+".inf"))
    file_data[index].mjd_i = int(file_data[index].epoch)
    file_data[index].mjd_f = file_data[index].epoch - file_data[index].mjd_i
    file_N.append(int(file_data[index].N + 1.0e-10))
    file_startMJDi.append(file_data[index].mjd_i)
    file_startMJDf.append(file_data[index].mjd_f)

    # Calculate the ending MJDs of each segment

    (MJDi, MJDf) = addtoMJD((file_data[index].dt * file_N[index]) \
                            / SECPERDAY, file_startMJDi[index], \
                            file_startMJDf[index])
    file_endMJDi.append(MJDi)
    file_endMJDf.append(MJDf)
    print('  %s.dat:  %9.0f pts at MJD %5d.%015.0f' % \
          (infile[index], file_N[index], \
           file_startMJDi[index], file_startMJDf[index] * 1.0e15))
    if (index > 0):
        if not (dt == file_data[index].dt):
            print('\nCannot concatenate the data.  The input file dt\'s')
            print('   are different.  Exiting.')
            exit(-1)
        else:
            dt = file_data[index].dt            
            
        # Calculate the number of bins of padding to use between
        # each data segment.

        padbins.append(subtractMJDs(file_startMJDi[index], \
                                    file_startMJDf[index], \
                                    file_endMJDi[index-1], \
                                    file_endMJDf[index-1]) \
                       * SECPERDAY / dt)
    else:
        dt = file_data[index].dt
                                            
print('')

# Convert the infodata into Numpy Arrays and determine the number of
# bins to add as padding as well as the shifts needed in the data sets

nf = len(file_data)
padbins = asarray(padbins);

# Calculate the number of whole bins of padding

wholebins = (padbins+SBF).astype('l')

# Calculate the shifts required to keep the latter data segment
# in proper phase with the first segment

shift = padbins - wholebins - 1.0
shift = where(less(shift, -1.0 + SBF), 0.0, shift)
shift = where(greater(shift, 1.0 - SBF), 0.0, shift)
for index in range(len(shift)):
    if (fabs(shift[index]) > SBF):
        file_N[index + 1] = file_N[index + 1] + 1;
shift = where(greater(fabs(shift), SBF), shift, 0.0)
wholebins = wholebins.tolist()

# Calculate the number of bins of padding to tack on the end

endpad = numpts - add.reduce(wholebins) - add.reduce(file_N)
if endpad:
    wholebins.append(endpad)

# Adjust the start MJDs for the shifting of bins in the latter
# data sets.

for index in range(len(shift)):
    if (shift[index] < -SBF):
        (MJDi, MJDf) = addtoMJD((1.0 + shift[index]) * dt / SECPERDAY, \
                               file_startMJDi[index+1], \
                               file_startMJDf[index+1])
        file_startMJDi[index+1] = MJDi
        file_startMJDf[index+1] = MJDf
        
# Show the user what shifts were required

print('The bin shifts requires to align the data files in phase are:')
print('  %s.dat:  %+f bins' % (infile[0], 0.0))
for index in range(len(shift)):
    print('  %s.dat:  %+f bins' % (infile[index+1], shift[index]))
print('')

# Show the user what the output files will consist of

print('The output file will consist of:')

outfile_N = []
commands = []
totalbins = 0
for index in range(len(wholebins)):
    outfile_N.append(file_N[index])
    (MJDi, MJDf) = addtoMJD(totalbins * dt / SECPERDAY, \
                            file_startMJDi[0], \
                            file_startMJDf[0])
    totalbins = totalbins + outfile_N[2 * index]
    print('     data:  %9.0f pts starting at MJD %5d.%015.0f' % \
          (outfile_N[2 * index], MJDi, MJDf * 1.0e15))
    if (index == 0):
        commands.append("  cp %s.dat %s" % (infile[0], outfilenm))
    else:
        commands.append("  shiftdata %f %s.dat >> %s" % \
                       (shift[index-1], infile[index], outfilenm))
    outfile_N.append(wholebins[index])
    (MJDi, MJDf) = addtoMJD(totalbins * dt / SECPERDAY, \
                            file_startMJDi[0], \
                            file_startMJDf[0])
    print('  padding:  %9.0f pts starting at MJD %5d.%015.0f' % \
          (outfile_N[2 * index + 1], MJDi, MJDf * 1.0e15))
    totalbins = totalbins + outfile_N[2 * index + 1]
    commands.append("  patchdata %ld %f >> %s" % \
                   (wholebins[index], padval, outfilenm))

if (len(wholebins) < len(file_N)):
    outfile_N.append(file_N[len(file_N)])
    (MJDi, MJDf) = addtoMJD(totalbins * dt / SECPERDAY, \
                            file_startMJDi[0], \
                            file_startMJDf[0])
    print('     data:  %9.0f pts starting at MJD %5d.%015.0f' % \
          (outfile_N[2 * index + 1], MJDi, MJDf * 1.0e15))
    commands.append("  shiftdata %f %s.dat >> %s" % \
                   (shift[len(file_N)-1], infile[len(file_N)], outfilenm))
print('')

# Show the user the commands we will use to concat everything

print('The commands to perform the concatenation will be:')

for index in range(len(commands)):
    print(commands[index])
print('')
