#include "multifiles.h"

#ifndef __USE_FILE_OFFSET64
#ifndef __USE_LARGEFILE
#define fseeko fseek
#endif
#endif

void print_multifile(multifile * mfile, int full)
/* Print a multifile structure.  If 'full' is true, */
/* print everything.  Used for debugging.           */
{
    if (full) {
        int ii;

        printf("There are %d files (maxlen = %lld, mode = '%s'):\n",
               mfile->numfiles, mfile->maxfilelen, mfile->mode);
        for (ii = 0; ii < mfile->numfiles; ii++) {
            printf("   %d:  '%s', len = %lld'\n",
                   ii, mfile->filenames[ii], mfile->filelens[ii]);
        }
    }
    printf(" File = %d, Filepos = %lld, Fullpos = %lld, Length = %lld\n",
           mfile->currentfile, mfile->currentpos, mfile->position, mfile->length);
}

multifile *fopen_multifile(int numfiles, char **filenames, char *mode,
                           long long maxlen)
/* Open a multifile for use and return the multifile structure.  */
/*  'numfiles' is the number of files in the multifile.          */
/*  'filenames' is an array of the names of the component files  */
/*  'mode' is the method of opening the file (binary is assumed) */
/* 	"r" : read only, do not create (truncate) the files      */
/*      "r+": read+write, do not create (truncate) the files     */
/*      "w" : read+write, create (truncate) the files            */
/*      "a" : read+write, create files or open at end of files   */
/*  'maxlen' is the maximum length in bytes of each file.  This  */
/*      number is only used if a file is opened fo writing.  The */
/*      default value is DEFAULT_MAXLEN.  If you want to use     */
/*      the default, simply set 'maxlen' to 0.                   */
{
    multifile *mfile;
    int ii, append = 0;

    mfile = (multifile *) malloc(sizeof(multifile));
    mfile->numfiles = numfiles;
    mfile->currentfile = 0;
    mfile->currentpos = 0;
    mfile->length = 0;
    mfile->position = 0;
    if (maxlen > 0)
        mfile->maxfilelen = maxlen;
    else
        mfile->maxfilelen = DEFAULT_MAXLEN;
    mfile->filelens = (long long *) malloc(numfiles * sizeof(long long));
    mfile->filenames = (char **) malloc(numfiles * sizeof(char *));
    mfile->fileptrs = (FILE **) malloc(numfiles * sizeof(FILE *));
    if (strncmp(mode, "r", 1) == 0) {
        if (strncmp(mode, "r+", 2) == 0) {
            strcpy(mfile->mode, "rb+");
        } else {
            strcpy(mfile->mode, "rb");
        }
    } else if (strncmp(mode, "w", 1) == 0) {
        if (strncmp(mode, "w+", 2) == 0) {
            strcpy(mfile->mode, "wb+");
        } else {
            strcpy(mfile->mode, "wb");
        }
    } else if (strncmp(mode, "a", 1) == 0) {
        if (strncmp(mode, "a+", 2) == 0) {
            strcpy(mfile->mode, "ab+");
        } else {
            strcpy(mfile->mode, "ab");
        }
        append = 1;
    } else {
        printf("\n'mode' = '%s' in open_multifile() is not valid.\n\n", mode);
        return NULL;
    }
    for (ii = 0; ii < numfiles; ii++) {
        mfile->filenames[ii] = (char *) calloc(strlen(filenames[ii]) + 1, 1);
        strcpy(mfile->filenames[ii], filenames[ii]);
        mfile->fileptrs[ii] = chkfopen(filenames[ii], mfile->mode);
        mfile->filelens[ii] = chkfilelen(mfile->fileptrs[ii], 1);
        mfile->length += mfile->filelens[ii];
    }
    if (append) {
        mfile->currentfile = numfiles - 1;
        chkfileseek(mfile->fileptrs[mfile->currentfile], 0, 1, SEEK_END);
    }
    return mfile;
}

int fclose_multifile(multifile * mfile)
/* Close a multifile and free its resources. */
{
    int ii, flag = 0;

    for (ii = 0; ii < mfile->numfiles; ii++) {
        flag = fclose(mfile->fileptrs[ii]);
        free(mfile->filenames[ii]);
    }
    free(mfile->filenames);
    free(mfile->fileptrs);
    free(mfile->filelens);
    free(mfile);
    return flag;
}

int fread_multifile(void *data, size_t type, size_t number, multifile * mfile)
/* Read binary data from a multifile.                        */
/*   'data' is an array of the proper type to store the data */
/*   'type' is the size of each nugget of data to read       */
/*   'number' is the number of nuggets to read               */
/*   'mfile' is a pointer to a valid multifile structure     */
{
    int findex;
    size_t readbytes, bytesread, tmpbytesread;

    findex = mfile->currentfile;
    readbytes = number * type;
    bytesread = chkfread((char *) data, 1, readbytes, mfile->fileptrs[findex]);
    mfile->currentpos += bytesread;
    mfile->position += bytesread;
    readbytes -= bytesread;
    while (readbytes) {
        if (feof(mfile->fileptrs[findex])) {
            if (findex == mfile->numfiles - 1) {
                return bytesread / type;
            } else {
                findex++;
                mfile->currentfile++;
                mfile->currentpos = 0;
                tmpbytesread = chkfread((char *) data + bytesread, 1,
                                        readbytes, mfile->fileptrs[findex]);
                bytesread += tmpbytesread;
                mfile->currentpos += tmpbytesread;
                mfile->position += tmpbytesread;
                readbytes -= tmpbytesread;
            }
        } else {
            printf("\nRead error in read_multifile():\n");
            printf("\tTried to read %zd bytes, only read %zd!\n\n",
                   number * type, bytesread);
            return bytesread / type;
        }
    }
    return bytesread / type;
}

int fwrite_multifile(void *data, size_t type, size_t number, multifile * mfile)
/* Write binary data to a multifile.                         */
/*   'data' is an array of data to write                     */
/*   'type' is the size of each nugget of data to write      */
/*   'number' is the number of nuggets to write              */
/*   'mfile' is a pointer to a valid multifile structure     */
{
    int findex;
    unsigned long long bytesleft;
    size_t writebytes, tmpwritebytes, byteswritten, tmpbyteswritten;

    findex = mfile->currentfile;
    writebytes = number * type;
    bytesleft = mfile->maxfilelen - mfile->currentpos;
    tmpwritebytes = (writebytes > bytesleft) ? bytesleft : writebytes;
    byteswritten = chkfwrite((char *) data, 1, tmpwritebytes,
                             mfile->fileptrs[findex]);
    mfile->currentpos += byteswritten;
    if (mfile->currentpos > mfile->filelens[findex]) {
        mfile->length += (mfile->currentpos - mfile->filelens[findex]);
        mfile->filelens[findex] = mfile->currentpos;
    }
    mfile->position += byteswritten;
    writebytes -= byteswritten;
    while (writebytes) {
        if (mfile->currentpos == mfile->maxfilelen) {
            if (findex < mfile->numfiles - 1) {
                findex++;
                mfile->currentfile++;
                mfile->currentpos = 0;
                bytesleft = mfile->maxfilelen;
                tmpwritebytes = (writebytes > bytesleft) ? bytesleft : writebytes;
                tmpbyteswritten = chkfwrite((char *) data + byteswritten,
                                            1, tmpwritebytes,
                                            mfile->fileptrs[findex]);
                byteswritten += tmpbyteswritten;
                mfile->currentpos += tmpbyteswritten;
                if (mfile->currentpos > mfile->filelens[findex]) {
                    mfile->length += (mfile->currentpos - mfile->filelens[findex]);
                    mfile->filelens[findex] = mfile->currentpos;
                }
                mfile->position += tmpbyteswritten;
                writebytes -= tmpbyteswritten;
            } else {
                printf("\nWrite error in write_multifile():\n");
                printf("\tMultifile is maximum length!\n\n");
                return byteswritten / type;
            }
        } else {
            printf("\nWrite error in write_multifile():\n");
            printf("\tTried to write %zd bytes, only wrote %zd!\n\n",
                   number * type, byteswritten);
            return byteswritten / type;
        }
    }
    return byteswritten / type;
}

int fseek_multifile(multifile * mfile, long long offset, int whence)
/* Set the multifile position indicator as in fseek().   */
/*   'offset' is the file offset in bytes.               */
/*   'whence' is either SEEK_SET, SEEK_CUR, or SEEK_END. */
/*   Note:  Return is 0 for success, -1 for failure      */
{
    int findex, rt;
    long long cumlen = 0, maxlen;

    if (whence == SEEK_SET) {
        mfile->position = offset;
    } else if (whence == SEEK_CUR) {
        mfile->position += offset;
    } else if (whence == SEEK_END) {
        mfile->position = mfile->length + offset;
    } else {
        printf("\nSeek error in seek_multifile():\n");
        printf("\tValue for 'whence' was %d!\n\n", whence);
    }
    maxlen = mfile->maxfilelen * mfile->numfiles;
    if (mfile->position > maxlen) {
        printf("\nWarning in fseek_multifile():\n");
        printf(" Seeking beyond max file length.\n");
        mfile->position = maxlen;
        findex = mfile->numfiles - 1;
        mfile->currentfile = findex;
        mfile->currentpos = mfile->maxfilelen;
        if ((rt = fseeko(mfile->fileptrs[findex], 0, SEEK_END)) == -1) {
            perror("\nError in fseek_multifile()");
            printf("\n");
        }
        return rt;
    }
    if (mfile->position <= 0) {
        findex = 0;
        mfile->position = 0;
        mfile->currentfile = 0;
        mfile->currentpos = 0;
    } else {
        /* Allow each file to be extended to mfile->maxfilelen */
        /* if the file has been opened in write mode.          */
        if (mfile->mode[0] == 'w') {
            findex = mfile->position / mfile->maxfilelen;
            mfile->currentfile = findex;
            mfile->currentpos = mfile->position - findex * mfile->maxfilelen;
        } else {
            findex = -1;
            while (cumlen <= mfile->position && findex < mfile->numfiles - 1) {
                findex++;
                cumlen += mfile->filelens[findex];
            }
            /* Extend the last file (for files opened with */
            /* any form of read access).                   */
            if (cumlen < mfile->position) {
                findex = mfile->numfiles - 1;
                mfile->currentfile = findex;
                mfile->currentpos = mfile->position - findex * mfile->maxfilelen;
            } else {
                mfile->currentfile = findex;
                mfile->currentpos =
                    mfile->position - (cumlen - mfile->filelens[findex]);
            }
        }
    }
    if ((rt = fseeko(mfile->fileptrs[findex], mfile->currentpos, SEEK_SET)) == -1) {
        perror("\nError in fseek_multifile()");
        printf("\n");
    }
    return rt;
}

void rewind_multifile(multifile * mfile)
/* Rewind a multifile. */
{
    fseek_multifile(mfile, 0, SEEK_SET);
}

long long ftell_multifile(multifile * mfile)
/* Report the current position of a muiltifile. */
{
    return mfile->position;
}
