#include "presto.h"
#include "mask.h"
#include "sigproc_fb.h"

static unsigned char *cdatabuffer;
static float *fdatabuffer;
static int currentfile = 0, currentblock = 0;
static int numbuffered = 0, numpadded = 0;

extern void add_padding(float *fdata, float *padding, int numchan, int numtopad);

/* Note:  Much of this has been ripped out of SIGPROC      */
/* and then slightly modified.  Thanks Dunc!               */

static void send_string(char *string, FILE * outfile)
{
    int len;
    len = strlen(string);
    chkfwrite(&len, sizeof(int), 1, outfile);
    chkfwrite(string, sizeof(char), len, outfile);
}

static void get_string(FILE * inputfile, int *nbytes, char string[])
{
    int nchar;
    strcpy(string, "ERROR");
    chkfread(&nchar, sizeof(int), 1, inputfile);
    *nbytes = sizeof(int);
    if (feof(inputfile)) exit(0);
    if (nchar > 80 || nchar < 1)
        return;
    chkfread(string, nchar, 1, inputfile);
    string[nchar] = '\0';
    *nbytes += nchar;
}

static int strings_equal(char *string1, char *string2)
{
    if (!strcmp(string1, string2)) {
        return 1;
    } else {
        return 0;
    }
}

static void send_double(char *name, double double_precision, FILE * outfile)
{
    send_string(name, outfile);
    chkfwrite(&double_precision, sizeof(double), 1, outfile);
}

static void send_int(char *name, int integer, FILE * outfile)
{
    send_string(name, outfile);
    chkfwrite(&integer, sizeof(int), 1, outfile);
}

static void send_coords(double raj, double dej, double az, double za, FILE * outfile)
{
    if ((raj != 0.0) || (raj != -1.0))
        send_double("src_raj", raj, outfile);
    if ((dej != 0.0) || (dej != -1.0))
        send_double("src_dej", dej, outfile);
    if ((az != 0.0) || (az != -1.0))
        send_double("az_start", az, outfile);
    if ((za != 0.0) || (za != -1.0))
        send_double("za_start", za, outfile);
}

void get_telescope_name(int telescope_id, struct spectra_info *s)
{
    float default_beam = 1.0;   // deg
    switch (telescope_id) {
    case 0:
        strcpy(s->telescope, "Fake");
        s->beam_FWHM = default_beam;
        break;
    case 1:
        strcpy(s->telescope, "Arecibo");
        // Don't use full 305m size for AO.  Use illuminated size.
        s->beam_FWHM = 2.0 / 3600.0 * beam_halfwidth(s->fctr, 200.0);
        break;
    case 2:
        strcpy(s->telescope, "Ooty");
        s->beam_FWHM = default_beam;
        break;
    case 3:
        strcpy(s->telescope, "Nancay");
        s->beam_FWHM = 2.0 / 3600.0 * beam_halfwidth(s->fctr, 100.0);
        break;
    case 4:
        strcpy(s->telescope, "Parkes");
        s->beam_FWHM = 2.0 / 3600.0 * beam_halfwidth(s->fctr, 64.0);
        break;
    case 5:
        strcpy(s->telescope, "Jodrell");
        s->beam_FWHM = 2.0 / 3600.0 * beam_halfwidth(s->fctr, 76.0);
        break;
    case 6:
        strcpy(s->telescope, "GBT");
        s->beam_FWHM = 2.0 / 3600.0 * beam_halfwidth(s->fctr, 100.0);
        break;
    case 7:
        strcpy(s->telescope, "GMRT");
        s->beam_FWHM = default_beam;
        break;
    case 8:
        strcpy(s->telescope, "Effelsberg");
        s->beam_FWHM = 2.0 / 3600.0 * beam_halfwidth(s->fctr, 100.0);
        break;
    case 9:
        strcpy(s->telescope, "ATA");
        s->beam_FWHM = default_beam;
        break;
    case 10:
        strcpy(s->telescope, "SRT");
        s->beam_FWHM = default_beam;
        break;
    case 11:
        strcpy(s->telescope, "LOFAR");
        s->beam_FWHM = default_beam;
        break;
    case 12:
        strcpy(s->telescope, "VLA");
        s->beam_FWHM = default_beam;
        break;
    case 20:  // May need to change....
        strcpy(s->telescope, "CHIME");
        s->beam_FWHM = 2.0 / 3600.0 * beam_halfwidth(s->fctr, 20.0);
        break;
    case 21:  // May need to change....
        strcpy(s->telescope, "FAST");
        s->beam_FWHM = 2.0 / 3600.0 * beam_halfwidth(s->fctr, 350.0);
        break;
    case 64:
        strcpy(s->telescope, "MeerKAT");
        s->beam_FWHM = default_beam;
        break;
    case 65:
        strcpy(s->telescope, "KAT-7");
        s->beam_FWHM = default_beam;
        break;
    default:
        strcpy(s->telescope, "Unknown");
        s->beam_FWHM = default_beam;
        break;
    }
}

void get_backend_name(int machine_id, struct spectra_info *s)
{
    char *backend, string[80];
    switch (machine_id) {
    case 0:
        strcpy(string, "FAKE");
        break;
    case 1:
        strcpy(string, "PSPM");
        break;
    case 2:
        strcpy(string, "WAPP");
        break;
    case 3:
        strcpy(string, "AOFTM");
        break;
    case 4:
        strcpy(string, "BPP");
        break;
    case 5:
        strcpy(string, "OOTY");
        break;
    case 6:
        strcpy(string, "SCAMP");
        break;
    case 7:
        strcpy(string, "SPIGOT");
        break;
    case 11:
        strcpy(string, "BG/P");
        break;
    case 12:
        strcpy(string, "PDEV");
        break;
    case 20:
        strcpy(string, "CHIME+PSR");
        break;
    case 64:
        strcpy(string, "KAT");
        break;
    case 65:
        strcpy(string, "KAT-DC2");
        break;
    default:
        strcpy(s->backend, "Unknown");
        break;
    }
    backend = (char *) calloc(strlen(string) + 1, 1);
    strcpy(backend, string);
}


void write_filterbank_header(sigprocfb * fb, FILE * outfile)
{
    int ii, jj;

    if (fb->machine_id != 0) {
        send_string("HEADER_START", outfile);
        send_string("rawdatafile", outfile);
        send_string(fb->inpfile, outfile);
        if (!strings_equal(fb->source_name, "")) {
            send_string("source_name", outfile);
            send_string(fb->source_name, outfile);
        }
        send_int("machine_id", fb->machine_id, outfile);
        send_int("telescope_id", fb->telescope_id, outfile);
        send_coords(fb->src_raj, fb->src_dej, fb->az_start, fb->za_start, outfile);
        send_int("data_type", 1, outfile);      /* filterbank data */
        send_double("fch1", fb->fch1, outfile);
        send_double("foff", fb->foff, outfile);
        send_int("nchans", fb->nchans, outfile);
        send_int("nbits", fb->nbits, outfile);
        send_double("tstart", fb->tstart, outfile);
        send_double("tsamp", fb->tsamp, outfile);
        if (fb->sumifs) {
            send_int("nifs", 1, outfile);
        } else {
            jj = 0;
            for (ii = 1; ii <= fb->nifs; ii++)
                if (fb->ifstream[ii - 1] == 'Y')
                    jj++;
            if (jj == 0)
                printf("\nNo valid IF streams selected!\n\n");
            send_int("nifs", jj, outfile);
        }
        send_string("HEADER_END", outfile);
    }
}

/* attempt to read in the general header info from a pulsar data file */
int read_filterbank_header(sigprocfb * fb, FILE * inputfile)
{
    char string[80], message[80];
    int itmp, nbytes = 0, totalbytes;
    int expecting_rawdatafile = 0, expecting_source_name = 0;
    int barycentric, pulsarcentric;
    /* try to read in the first line of the header */
    get_string(inputfile, &nbytes, string);
    if (!strings_equal(string, "HEADER_START")) {
        /* the data file is not in standard format, rewind and return */
        rewind(inputfile);
        return 0;
    }
    /* store total number of bytes read so far */
    totalbytes = nbytes;

    fb->ibeam = 1;              // default value
    fb->signedints = 0;         // default value
    /* loop over and read remaining header lines until HEADER_END reached */
    while (1) {
        get_string(inputfile, &nbytes, string);
        if (strings_equal(string, "HEADER_END"))
            break;
        totalbytes += nbytes;
        if (strings_equal(string, "rawdatafile")) {
            expecting_rawdatafile = 1;
        } else if (strings_equal(string, "source_name")) {
            expecting_source_name = 1;
        } else if (strings_equal(string, "az_start")) {
            chkfread(&(fb->az_start), sizeof(double), 1, inputfile);
            totalbytes += sizeof(double);
        } else if (strings_equal(string, "za_start")) {
            chkfread(&(fb->za_start), sizeof(double), 1, inputfile);
            totalbytes += sizeof(double);
        } else if (strings_equal(string, "src_raj")) {
            chkfread(&(fb->src_raj), sizeof(double), 1, inputfile);
            totalbytes += sizeof(double);
        } else if (strings_equal(string, "src_dej")) {
            chkfread(&(fb->src_dej), sizeof(double), 1, inputfile);
            totalbytes += sizeof(double);
        } else if (strings_equal(string, "tstart")) {
            chkfread(&(fb->tstart), sizeof(double), 1, inputfile);
            totalbytes += sizeof(double);
        } else if (strings_equal(string, "tsamp")) {
            chkfread(&(fb->tsamp), sizeof(double), 1, inputfile);
            totalbytes += sizeof(double);
        } else if (strings_equal(string, "fch1")) {
            chkfread(&(fb->fch1), sizeof(double), 1, inputfile);
            totalbytes += sizeof(double);
        } else if (strings_equal(string, "foff")) {
            chkfread(&(fb->foff), sizeof(double), 1, inputfile);
            totalbytes += sizeof(double);
        } else if (strings_equal(string, "nchans")) {
            chkfread(&(fb->nchans), sizeof(int), 1, inputfile);
            totalbytes += sizeof(int);
        } else if (strings_equal(string, "telescope_id")) {
            chkfread(&(fb->telescope_id), sizeof(int), 1, inputfile);
            totalbytes += sizeof(int);
        } else if (strings_equal(string, "machine_id")) {
            chkfread(&(fb->machine_id), sizeof(int), 1, inputfile);
            totalbytes += sizeof(int);
        } else if (strings_equal(string, "data_type")) {
            chkfread(&itmp, sizeof(int), 1, inputfile);
            totalbytes += sizeof(int);
        } else if (strings_equal(string, "nbits")) {
            chkfread(&(fb->nbits), sizeof(int), 1, inputfile);
            totalbytes += sizeof(int);
        } else if (strings_equal(string, "barycentric")) {
            chkfread(&barycentric, sizeof(int), 1, inputfile);
            totalbytes += sizeof(int);
        } else if (strings_equal(string, "pulsarcentric")) {
            chkfread(&pulsarcentric, sizeof(int), 1, inputfile);
            totalbytes += sizeof(int);
        } else if (strings_equal(string, "nsamples")) {
            /* read this one only for backwards compatibility */
            chkfread(&itmp, sizeof(int), 1, inputfile);
            totalbytes += sizeof(int);
        } else if (strings_equal(string, "nifs")) {
            chkfread(&(fb->nifs), sizeof(int), 1, inputfile);
            totalbytes += sizeof(int);
        } else if (strings_equal(string, "nbeams")) {
            chkfread(&(fb->nbeams), sizeof(int), 1, inputfile);
            totalbytes += sizeof(int);
        } else if (strings_equal(string, "ibeam")) {
            chkfread(&(fb->ibeam), sizeof(int), 1, inputfile);
            totalbytes += sizeof(int);
        } else if (strings_equal(string, "signed")) {
            char tmp;
            chkfread(&tmp, sizeof(char), 1, inputfile);
            fb->signedints = tmp;
            totalbytes += sizeof(int);
        } else if (expecting_rawdatafile) {
            strcpy(fb->inpfile, string);
            expecting_rawdatafile = 0;
        } else if (expecting_source_name) {
            strcpy(fb->source_name, string);
            expecting_source_name = 0;
        } else {
            sprintf(message, "read_filterbank_header - unknown parameter: %s\n",
                    string);
            fprintf(stderr, "ERROR: %s\n", message);
            exit(1);
        }
    }
    /* add on last header string */
    totalbytes += nbytes;
    /* return total number of bytes read */
    fb->headerlen = totalbytes;
    /* Calculate the number of samples in the file */
    fb->N = (chkfilelen(inputfile, 1) - fb->headerlen) / fb->nchans * 8 / fb->nbits;
    return totalbytes;
}



void read_filterbank_files(struct spectra_info *s)
{
    sigprocfb fb;
    int ii;

    // s->num_files and s->filenames are assumed to be set
    s->datatype = SIGPROCFB;
    s->files = (FILE **) malloc(sizeof(FILE *) * s->num_files);
    s->header_offset = gen_ivect(s->num_files);
    // The following two aren't used for filterbank data,
    // but they should be initialized for mpiprepsubband
    s->start_subint = gen_ivect(s->num_files);
    for (ii = 0; ii < s->num_files; ii++)
        s->start_subint[ii] = 0;
    s->num_subint = gen_ivect(s->num_files);
    for (ii = 0; ii < s->num_files; ii++)
        s->start_subint[ii] = 0;
    s->start_spec = (long long *) malloc(sizeof(long long) * s->num_files);
    s->num_spec = (long long *) malloc(sizeof(long long) * s->num_files);
    s->num_pad = (long long *) malloc(sizeof(long long) * s->num_files);
    s->start_MJD = (long double *) malloc(sizeof(long double) * s->num_files);
#if DEBUG_OUT
    printf("Reading '%s'\n", s->filenames[0]);
#endif
    s->files[0] = chkfopen(s->filenames[0], "r");
    // Read the filterbank header into a SIGPROCFB struct
    s->header_offset[0] = read_filterbank_header(&fb, s->files[0]);
    // Make an initial offset into the file
    chkfseek(s->files[0], s->header_offset[0], SEEK_SET);
    strncpy(s->source, fb.source_name, 40);
    if (fb.sumifs) {
        s->summed_polns = 1;
        s->num_polns = 1;
    } else {
        s->num_polns = fb.nifs;
        strncpy(s->poln_order, fb.ifstream, 8);
    }
    // Position info
    {
        int d, h, m;
        double sec;
        h = (int) floor(fb.src_raj / 10000.0);
        m = (int) floor((fb.src_raj - h * 10000) / 100.0);
        sec = fb.src_raj - h * 10000 - m * 100;
        ra_dec_to_string(s->ra_str, h, m, sec);
        s->ra2000 = hms2rad(h, m, sec) * RADTODEG;
        d = (int) floor(fabs(fb.src_dej) / 10000.0);
        m = (int) floor((fabs(fb.src_dej) - d * 10000) / 100.0);
        sec = fabs(fb.src_dej) - d * 10000 - m * 100;
        if (fb.src_dej < 0.0)
            d = -d;
        ra_dec_to_string(s->dec_str, d, m, sec);
        s->dec2000 = dms2rad(d, m, sec) * RADTODEG;
    }
    s->bits_per_sample = fb.nbits;
    s->signedints = fb.signedints;
    s->num_channels = fb.nchans;
    s->samples_per_spectra = s->num_polns * s->num_channels;
    s->bytes_per_spectra = s->bits_per_sample * s->samples_per_spectra / 8;
    s->spectra_per_subint = 2400;        // use this as the blocksize
    s->bytes_per_subint = s->bytes_per_spectra * s->spectra_per_subint;
    s->samples_per_subint = s->spectra_per_subint * s->samples_per_spectra;
    s->min_spect_per_read = 1;  // Can read a single spectra at a time
    // allocate the raw data buffers
    cdatabuffer = gen_bvect(s->bytes_per_subint);
    fdatabuffer = gen_fvect(s->spectra_per_subint * s->num_channels);
    s->padvals = gen_fvect(s->num_channels);
    for (ii = 0; ii < s->num_channels; ii++)
        s->padvals[ii] = 0.0;
    s->dt = fb.tsamp;
    s->time_per_subint = s->spectra_per_subint * s->dt;
    s->T = s->N * s->dt;
    s->df = fabs(fb.foff);
    if (fb.foff < 0.0 && s->apply_flipband == -1)
        s->apply_flipband = 0;  // we do this automatically
    s->BW = s->num_channels * s->df;
    s->lo_freq = fb.fch1 - (s->num_channels - 1) * s->df;
    s->hi_freq = fb.fch1;
    s->fctr = s->lo_freq - 0.5 * s->df + 0.5 * s->BW;
    s->azimuth = fb.az_start;
    s->zenith_ang = fb.za_start;
    s->num_beams = 1;
    s->beamnum = fb.ibeam;
    get_telescope_name(fb.telescope_id, s);
    get_backend_name(fb.machine_id, s);
    s->start_MJD[0] = fb.tstart;
    mjd_to_datestr(s->start_MJD[0], s->date_obs);
    s->start_spec[0] = 0L;
    s->num_spec[0] = fb.N;
    s->N = fb.N;
    s->get_rawblock = &get_filterbank_rawblock;
    s->offset_to_spectra = &offset_to_filterbank_spectra;

    // Step through the other files
    for (ii = 1; ii < s->num_files; ii++) {

#if DEBUG_OUT
        printf("Reading '%s'\n", s->filenames[ii]);
#endif
        s->files[ii] = chkfopen(s->filenames[ii], "r");
        s->header_offset[ii] = read_filterbank_header(&fb, s->files[ii]);
        // Make an initial offset into each file to the spactra
        chkfseek(s->files[ii], s->header_offset[ii], SEEK_SET);
        // Compare key values with s->XXX[0] to see if things are the same
        if (s->num_channels != fb.nchans) {
            fprintf(stderr,
                    "Error:  num chans %d in file #%d does not match original num chans %d!!\n",
                    fb.nchans, ii + 1, s->num_channels);
            exit(1);
        }
        if (s->bits_per_sample != fb.nbits) {
            fprintf(stderr,
                    "Error:  bits per sample %d in file #%d does not match original bits per sample %d!!\n",
                    fb.nbits, ii + 1, s->bits_per_sample);
            exit(1);
        }
        if (s->dt != fb.tsamp) {
            fprintf(stderr,
                    "Error:  sample time %f in file #%d does not match original sample time %f!!\n",
                    fb.tsamp, ii + 1, s->dt);
            exit(1);
        }
        if (s->df != fabs(fb.foff)) {
            fprintf(stderr,
                    "Error:  channel width %f in file #%d does not match original channel width %f!!\n",
                    fabs(fb.foff), ii + 1, s->df);
            exit(1);
        }
        if (s->hi_freq != fb.fch1) {
            fprintf(stderr,
                    "Error:  high chan freq %f in file #%d does not match original high chan freq %f!!\n",
                    fb.fch1, ii + 1, s->hi_freq);
            exit(1);
        }
        s->start_MJD[ii] = fb.tstart;
        s->start_spec[ii] =
            (long long) ((s->start_MJD[ii] - s->start_MJD[0]) * SECPERDAY / s->dt +
                         0.5);
        s->num_pad[ii - 1] = s->start_spec[ii] - s->N;
        s->num_spec[ii] = fb.N;
        s->N += s->num_spec[ii] + s->num_pad[ii - 1];
    }
    s->T = s->N * s->dt;
    s->num_pad[s->num_files - 1] = 0L;
}


long long offset_to_filterbank_spectra(long long specnum, struct spectra_info *s)
// This routine offsets into the filterbank files to the spectra
// 'specnum'.  It returns the current spectra number.
{
    long long offset_spec;
    int filenum = 0;

    if (specnum > s->N) {
        fprintf(stderr, "Error:  offset spectra %lld is > total spectra %lld\n\n",
                specnum, s->N);
        exit(1);
    }
    // Find which file we need
    while (filenum + 1 < s->num_files && specnum > s->start_spec[filenum + 1])
        filenum++;

    // Shift to that file
    currentfile = filenum;

    // Are we in a padding zone?
    if (specnum > (s->start_spec[currentfile] + s->num_spec[currentfile])) {
        // Seek to the end of the file
        chkfileseek(s->files[currentfile], 0, 1, SEEK_END);
        numbuffered = 0;
        numpadded =
            specnum - (s->start_spec[currentfile] + s->num_spec[currentfile]);
        return specnum;
    }
    // Otherwise, seek to the spectra
    // Skip the header first
    chkfseek(s->files[currentfile], s->header_offset[currentfile], SEEK_SET);
    offset_spec = specnum - s->start_spec[currentfile];
    chkfseek(s->files[currentfile], offset_spec * s->bytes_per_spectra, SEEK_CUR);
    numbuffered = 0;
    numpadded = 0;
    return specnum;
}


int get_filterbank_rawblock(float *fdata, struct spectra_info *s, int *padding)
// This routine reads a single block (i.e subint) from the input files
// which contain raw data in SIGPROC filterbank format.  If padding is
// returned as 1, then padding was added and statistics should not be
// calculated.  Return 1 on success.
{
    int numread, numtopad = 0, numtoread;
    float *fdataptr = fdata;

    // If we only made a partial read last time, adjust our pointer
    // offsetting into the output floating-point array
    fdataptr = fdata + numbuffered * s->num_channels;
    numtoread = s->spectra_per_subint - numbuffered;

    /* Make sure our current file number is valid */
    if (currentfile >= s->num_files)
        return 0;

    /* First, attempt to read data from the current file */
    numread = chkfread(cdatabuffer, s->bytes_per_spectra,
                       numtoread, s->files[currentfile]);
    if (s->flip_bytes) {        // byte-swap if necessary
        /* Need to add this later */
    }
    convert_filterbank_block(fdataptr, cdatabuffer, numread, s);

    if (numread == numtoread) { // Got all we needed
        *padding = 0;
        numbuffered = 0;
        currentblock++;
        goto return_block;
    } else {                    // Didn't get all the data we needed
        numbuffered += numread;
        if (feof(s->files[currentfile])) {      // End of file?
            numtopad = s->num_pad[currentfile] - numpadded;
            if (numtopad == 0) {        // No padding needed.  Try reading the next file
                currentfile++;
                return get_filterbank_rawblock(fdata, s, padding);
            } else {            // Pad the data instead
                *padding = 1;
                fdataptr = fdata + numbuffered * s->num_channels;
                if (numtopad >= numtoread - numread) {  // Need lots of padding
                    // Fill the rest of the buffer with padding
                    numtopad = s->spectra_per_subint - numbuffered;
                    add_padding(fdataptr, s->padvals, s->num_channels, numtopad);
                    numpadded += numtopad;
                    numbuffered = 0;
                    currentblock++;
                    // If done with padding reset padding variables
                    if (numpadded == s->num_pad[currentfile]) {
                        numpadded = 0;
                        currentfile++;
                    }
                    goto return_block;
                } else {        // Need < 1 block (or remaining block) of padding
                    add_padding(fdataptr, s->padvals, s->num_channels, numtopad);
                    numbuffered += numtopad;
                    // Done with padding, so reset padding variables
                    numpadded = 0;
                    currentfile++;
                    return get_filterbank_rawblock(fdata, s, padding);
                }
            }
        } else {
            fprintf(stderr,
                    "Error: Problem reading record from filterbank data file:\n");
            fprintf(stderr, "   currentfile = %d, currentblock = %d.  Exiting.\n",
                    currentfile, currentblock);
            exit(1);
        }
    }

  return_block:
    // Apply the corrections that need a full block

    // Invert the band if requested
    if (s->apply_flipband)
        flip_band(fdata, s);

    // Perform Zero-DMing if requested
    if (s->remove_zerodm)
        remove_zerodm(fdata, s);

    return 1;
}


void convert_filterbank_block(float *outdata, unsigned char *indata,
                              int numread, struct spectra_info *s)
/* This routine converts SIGPROC filterbank-format data into PRESTO format */
{
    int ii, jj, spec_ct, offset;

    if (s->bits_per_sample == 32) {
        float *floatdata = (float *) indata;
        for (spec_ct = 0; spec_ct < numread; spec_ct++) {
            offset = spec_ct * s->num_channels;
            for (ii = 0, jj = s->num_channels - 1; ii < s->num_channels; ii++, jj--)
                outdata[ii + offset] = floatdata[jj + offset];
        }
    } else if (s->bits_per_sample == 16) {
        unsigned short *sdata = (unsigned short *) indata;
        for (spec_ct = 0; spec_ct < numread; spec_ct++) {
            offset = spec_ct * s->num_channels;
            for (ii = 0, jj = s->num_channels - 1; ii < s->num_channels; ii++, jj--)
                outdata[ii + offset] = (float) sdata[jj + offset];
        }
    } else if (s->bits_per_sample == 8) {
        if (s->signedints) {
            char *chardata = (char *) indata;
            for (spec_ct = 0; spec_ct < numread; spec_ct++) {
                offset = spec_ct * s->num_channels;
                for (ii = 0, jj = s->num_channels - 1; ii < s->num_channels; ii++, jj--)
                    outdata[ii + offset] = (float) chardata[jj + offset];
            }
        } else {
            unsigned char *chardata = (unsigned char *) indata;
            for (spec_ct = 0; spec_ct < numread; spec_ct++) {
                offset = spec_ct * s->num_channels;
                for (ii = 0, jj = s->num_channels - 1; ii < s->num_channels; ii++, jj--)
                    outdata[ii + offset] = (float) chardata[jj + offset];
            }
        }
    } else if (s->bits_per_sample == 4) {
        unsigned char c, *chardata = (unsigned char *) indata;
        for (spec_ct = 0; spec_ct < numread; spec_ct++) {
            offset = spec_ct * s->num_channels;
            for (ii = 0, jj = s->num_channels / 2 - 1; ii < s->num_channels;
                 ii += 2, jj--) {
                c = chardata[(jj + offset / 2)];
                outdata[(ii + 1) + offset] = (float) (c & 15);
                outdata[ii + offset] = (float) ((c & 240) >> 4);
            }
        }
    } else if (s->bits_per_sample == 2) {
        unsigned char c, *chardata = (unsigned char *) indata;
        float *outdataptr;
        for (spec_ct = 0; spec_ct < numread; spec_ct++) {
            offset = spec_ct * s->num_channels;
            outdataptr = outdata + offset;
            for (ii = 0, jj = s->num_channels / 4 - 1; ii < s->num_channels;
                 ii += 4, jj--) {
                c = chardata[(jj + offset / 4)];
                *outdataptr++ = (float) ((c >> 0x06) & 0x03);
                *outdataptr++ = (float) ((c >> 0x04) & 0x03);
                *outdataptr++ = (float) ((c >> 0x02) & 0x03);
                *outdataptr++ = (float) (c & 0x03);
            }
        }
    } else if (s->bits_per_sample == 1) {
        unsigned char c, *chardata = (unsigned char *) indata;
        float *outdataptr;
        for (spec_ct = 0; spec_ct < numread; spec_ct++) {
            offset = spec_ct * s->num_channels;
            outdataptr = outdata + offset;
            for (ii = 0, jj = s->num_channels / 8 - 1; ii < s->num_channels;
                 ii += 8, jj--) {
                c = chardata[(jj + offset / 8)];
                *outdataptr++ = (float) ((c >> 0x07) & 0x01);
                *outdataptr++ = (float) ((c >> 0x06) & 0x01);
                *outdataptr++ = (float) ((c >> 0x05) & 0x01);
                *outdataptr++ = (float) ((c >> 0x04) & 0x01);
                *outdataptr++ = (float) ((c >> 0x03) & 0x01);
                *outdataptr++ = (float) ((c >> 0x02) & 0x01);
                *outdataptr++ = (float) ((c >> 0x01) & 0x01);
                *outdataptr++ = (float) (c & 0x01);
            }
        }
    } else {
        printf
            ("\nYikes!!! Not supposed to be here in convert_filterbank_block()\n\n");
    }
}
