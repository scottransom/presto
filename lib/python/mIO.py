## Automatically adapted for numpy Apr 14, 2006 by convertcode.py

from Numeric import *
import numpy as Numeric
N = Numeric
import scipy.io.numpyio as numpyio
import struct

def getsize_type(mtype):
    if mtype in ['B','uchar','byte','unsigned char','integer*1', 'int8']:
        mtype = 'B'
    elif (mtype == 'c') or mtype in ['char','char*1']:
        mtype = 'c'
    elif mtype in ['b','schar', 'signed char']:
        mtype = 'b'
    elif mtype in ['s','short','int16','integer*2']:
        mtype = 's'
    elif mtype in ['l','i','int','long','int32','integer*4']:
        mtype = 'l'
    elif mtype in ['f','float','float32','real*4']:
        mtype = 'f'
    elif mtype in ['d','double','float64','real*8']:
        mtype = 'd'
    elif mtype in ['F','complex float','complex*8','complex64']:
        mtype = 'F'
    elif mtype in ['D','complex*16','complex128','complex','complex double']:
        mtype = 'D'
    else:
        raise TypeError, 'Bad datatype -- ' + mtype

    argout = (array(0,mtype).itemsize,mtype)
    return argout

class binary_file:
    def __init__(self,file_name,permission='r',format='n'):
        if type(file_name) == type(''):
            self.fid = open(file_name,permission)
        elif 'fileno' in file_name.__methods__:  # first argument is an open file
            self.fid = file_name 
        if format in ['native','n']:
            self.bs = 0
            self.format = 'native'
        elif format in ['ieee-le','l']:
            self.bs = not LittleEndian
            self.format = 'ieee-le'
        elif format in ['ieee-be','B']:
            self.bs = LittleEndian
            self.format = 'ieee-be'

        self.seek = self.fid.seek
        self.tell = self.fid.tell
        self.close = self.fid.close
        self.fileno = self.fid.fileno
        self.mode = self.fid.mode
        self.closed = self.fid.closed
        self.name = self.fid.name

    def __del__(self):
        try:
            self.fid.close()
        except:
            pass
    
    def fwrite(self,data,mtype):
        howmany,mtype = getsize_type(mtype)
        data = asarray(data)
        count = product(data.shape)
        val = numpyio.fwrite(self.fid,count,data,mtype,self.bs)
        return val

    def fread(self,count,mtype):
        howmany,mtype = getsize_type(mtype)
        retval = numpyio.fread(self.fid, count, mtype, mtype, self.bs)
        if len(retval) == 1:
            retval = retval[0]
        return retval

    def rewind(self,howmany=None):
        if howmany is None:
            self.seek(0)
        else:
            self.seek(-howmany,1)

    def size(self):
        try:
            sz = self.thesize
        except AttributeError:            
            curpos = self.tell()
            self.fid.seek(0,2)
            sz = self.fid.tell()
            self.fid.seek(curpos)
            self.thesize = sz
        return sz

    def fort_write(self,fmt,*args):
        if self.format == 'ieee-le':
            nfmt = "<i"
        elif self.format == 'ieee-be':
            nfmt = ">i"
        else:
            nfmt = "i"
        if type(fmt) == type(''):
            if self.format == 'ieee-le':
                fmt = "<"+fmt
            elif self.format == 'ieee-be':
                fmt = ">"+fmt
            str = apply(struct.pack,(fmt,)+args)
            strlen = struct.pack(nfmt,len(str))
            self.fid.write(strlen)
            self.fid.write(str)
            self.fid.write(strlen)
        elif type(fmt) == type(array([0])):
            sz,mtype = getsize_type(args[0])
            count = product(fmt.shape)
            strlen = struct.pack(nfmt,count*sz)
            self.fid.write(strlen)
            numpyio.fwrite(self.fid,count,fmt,mtype,self.bs)
            self.fid.write(strlen)
        else:
            raise TypeError, "Unknown type in first argument"

    def fort_read(self,fmt,dtype=None):
        lookup_dict = {'ieee-le':"<",'ieee-be':">",'native':''}
        if dtype is None:
            fmt = lookup_dict[self.format] + fmt
            numbytes = struct.calcsize(fmt)
            nn = struct.calcsize("i");
            self.fid.read(nn)
            data = struct.unpack(fmt,self.fid.read(numbytes))
            self.fid.read(nn)
            return data
        else:  # Ignore format string and read in next record as an array.
            fmt = lookup_dict[self.format] + "i"
            nn = struct.calcsize(fmt)
            nbytes = struct.unpack(fmt,self.fid.read(nn))[0]
            howmany, dtype = getsize_type(dtype)
            ncount = nbytes / howmany
            if ncount*howmany != nbytes:
                self.rewind(4)
                raise ValueError, "A mismatch between the type requested and the data stored."
            retval = numpyio.fread(self.fid, ncount, dtype, dtype, self.bs)
            if len(retval) == 1:
                retval = retval[0]
            self.fid.read(nn)
            return retval
                                      
fopen = binary_file

def read_raw(name,size,type,bs=0):
    if bs:
        fid = fopen(name,'r','B')
    else:
        fid = fopen(name,'r')
    numels = N.product(size)
    data = fid.fread(numels,type)
    data.shape = size
    fid.close()
    return data
    
def write_raw(name,arr,type,bs=0):
    if bs:
        fid = fopen(name,'w','B')
    else:
        fid = fopen(name,'w')
    fid.fwrite(arr,type)
    fid.close()
    return

