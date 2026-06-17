/* Stage 3 unit check: prove that the factored read-once/dedisperse-many
 * reader produces bit-identical subbands to the original read_subbands().
 *
 * Pass A drives the new reader: read_clean_rawblock() once per block, then
 * dedisp_subbands() at a single DM.  Pass B rewinds and drives the original
 * read_subbands() at the same DM.  The two subband streams must match byte
 * for byte.  Exit status 0 means identical, 1 means a mismatch was found.
 *
 * If one argument ends in ".mask" it is read as an rfifind mask and applied
 * in both passes, so the shared-once masking path is checked too.
 *
 * Usage:  reader_equiv_check <rawfile> [more raw files ...] [maskfile.mask]
 */
#include "presto.h"
#include "backend_common.h"

int main(int argc, char *argv[])
{
    struct spectra_info s;
    infodata idata;
    mask obsmask;
    int ii, numread, padding = 0, nummasked = 0;
    int numchan, nsub, sps;
    int *idispdts, *maskchans;
    double dm = 242.30;
    double *dispdts;
    float *fdata_new, *fdata_old, *stored = NULL;
    long long nblocks_stored = 0, block_floats;
    rawblock_reader *reader;
    char *maskfile = NULL;
    char *rawfiles[16384];
    int numraw = 0;

    if (argc < 2) {
        fprintf(stderr, "Usage: %s <rawfile> [more raw files ...] [maskfile.mask]\n", argv[0]);
        return 2;
    }

    /* Separate an optional ".mask" argument from the raw data file names. */
    for (ii = 1; ii < argc; ii++) {
        size_t len = strlen(argv[ii]);
        if (len > 5 && strcmp(argv[ii] + len - 5, ".mask") == 0)
            maskfile = argv[ii];
        else
            rawfiles[numraw++] = argv[ii];
    }
    if (numraw == 0) {
        fprintf(stderr, "Error: no raw data files given.\n");
        return 2;
    }

    /* --- Shared observation setup (mirrors prepfold/prepsubband) --- */
    spectra_info_set_defaults(&s);
    s.filenames = rawfiles;
    s.num_files = numraw;
    identify_psrdatatype(&s, 0);
    if (s.datatype != PSRFITS && s.datatype != SIGPROCFB) {
        fprintf(stderr, "Error: this check only supports raw (PSRFITS/filterbank) data.\n");
        return 2;
    }
    read_rawdata_files(&s);
    spectra_info_to_inf(&s, &idata);
    s.clip_sigma = 0.0;                         /* no clipping */
    if (maskfile) {
        read_mask(maskfile, &obsmask);
        printf("reader_equiv_check: applying mask '%s'\n", maskfile);
    } else {
        obsmask.numchan = obsmask.numint = 0;   /* no mask */
    }

    numchan = s.num_channels;
    sps = s.spectra_per_subint;
    nsub = (numchan % 32 == 0) ? 32 : numchan;  /* a divisor of numchan */
    block_floats = (long long) nsub * sps;

    /* --- Per-candidate dispersion delays (same recipe as prepfold) --- */
    dispdts = subband_search_delays(numchan, nsub, dm,
                                    idata.freq, idata.chan_wid, 0.0);
    idispdts = gen_ivect(numchan);
    for (ii = 0; ii < numchan; ii++)
        idispdts[ii] = (int) (dispdts[ii] / s.dt + 0.5);
    vect_free(dispdts);

    maskchans = gen_ivect(numchan);
    fdata_new = gen_fvect(block_floats);
    fdata_old = gen_fvect(block_floats);

    printf("reader_equiv_check: numchan=%d nsub=%d spectra_per_subint=%d dm=%.2f\n",
           numchan, nsub, sps, dm);

    /* --- Pass A: new reader, read once + dedisperse once per block --- */
    reader = rawblock_reader_init(&s, &obsmask);
    /* Prime the reader (fills the "last" slot, returns 0). */
    (void) read_clean_rawblock(reader, &s, maskchans, &nummasked, &obsmask, &padding);
    while (read_clean_rawblock(reader, &s, maskchans, &nummasked, &obsmask, &padding)) {
        dedisp_subbands(reader->current_clean, reader->last_clean,
                        sps, numchan, idispdts, nsub, fdata_new);
        rawblock_reader_advance(reader);
        stored = (float *) realloc(stored,
                                   (size_t) (nblocks_stored + 1) * block_floats * sizeof(float));
        memcpy(stored + nblocks_stored * block_floats, fdata_new,
               block_floats * sizeof(float));
        nblocks_stored++;
    }
    free_rawblock_reader(reader);
    printf("Pass A (new reader): %lld foldable blocks\n", nblocks_stored);

    /* --- Rewind to the start of the observation for Pass B --- */
    offset_to_spectra(0, &s);

    /* --- Pass B: original read_subbands at the same DM --- */
    long long nblocks_old = 0;
    int mismatch = 0;
    while ((numread = read_subbands(fdata_old, idispdts, nsub, &s, 1, &padding,
                                    maskchans, &nummasked, &obsmask))) {
        if (nblocks_old >= nblocks_stored) {
            fprintf(stderr, "Error: read_subbands produced more blocks (%lld) than the new reader (%lld)\n",
                    nblocks_old + 1, nblocks_stored);
            mismatch = 1;
            break;
        }
        if (memcmp(stored + nblocks_old * block_floats, fdata_old,
                   block_floats * sizeof(float)) != 0) {
            fprintf(stderr, "MISMATCH in subband block %lld\n", nblocks_old);
            mismatch = 1;
            break;
        }
        nblocks_old++;
    }
    printf("Pass B (read_subbands): %lld foldable blocks\n", nblocks_old);

    if (!mismatch && nblocks_old != nblocks_stored) {
        fprintf(stderr, "Error: block count mismatch (new=%lld old=%lld)\n",
                nblocks_stored, nblocks_old);
        mismatch = 1;
    }

    free(stored);
    vect_free(idispdts);
    vect_free(maskchans);
    vect_free(fdata_new);
    vect_free(fdata_old);

    if (mismatch) {
        printf("FAIL: subband outputs differ.\n");
        return 1;
    }
    printf("PASS: new reader + dedisp_subbands is bit-identical to read_subbands over %lld blocks.\n",
           nblocks_stored);
    return 0;
}
