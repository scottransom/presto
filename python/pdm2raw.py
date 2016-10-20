from __future__ import print_function
from numpyio import fread, fwrite
from sys import argv

print("\nReading info from %s.hdr and" % argv[1])
print("%s.dat\n" % argv[1])
print("Writing %s.raw\n" % argv[1])

HEADERLEN = 640
BLOCKLEN = 49152

# Read the header file

file = open(argv[1]+'.hdr', 'r')
data = fread(file, HEADERLEN+8, 'b')
file.close()
header = data[4:-4]
infile = open(argv[1]+'.dat', 'r')
outfile = open(argv[1]+'.raw', 'w')

# Read and write the raw data

while (1):
    data = fread(infile, BLOCKLEN+8, 'b')
    if (len(data)==BLOCKLEN+8):
        fwrite(outfile, HEADERLEN, header, 'b')
        fwrite(outfile, BLOCKLEN, data[4:-4], 'b')
    else:
        break
print('')
infile.close()
outfile.close()
