#include "presto.h"
#include "math.h"
#include "rednoise_cmd.h"
#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

int main(int argc, char *argv[])
{
    FILE *infile, *outfile;
    char *rootfilenm, *outname;
    Cmdline *cmd;

    /* Call usage() if we have no command line arguments */

    if (argc == 1) {
        Program = argv[0];
        printf("\n");
        usage();
        exit(1);
    }

    /* Parse the command line using the excellent program Clig */

    cmd = parseCmdline(argc, argv);

#ifdef DEBUG
    showOptionValues();
#endif

    printf("\n\n");
    printf("     Rednoise Removal Routine v2.0\n");
    printf("            Feb, 2015\n\n");

    {
        int hassuffix = 0;
        char *suffix;
        hassuffix = split_root_suffix(cmd->argv[0], &rootfilenm, &suffix);
        if (hassuffix) {
            if (strcmp(suffix, "fft") != 0) {
                printf
                    ("\nInput file ('%s') must be a fourier transform ('.fft')!\n\n",
                     cmd->argv[0]);
                free(suffix);
                exit(0);
            }
            free(suffix);
        } else {
            printf("\nInput file ('%s') must be a fourier transform ('.fft')!\n\n",
                   cmd->argv[0]);
            exit(0);
        }
        outname = (char *) calloc(strlen(rootfilenm) + 11, sizeof(char));
        sprintf(outname, "%s_red.fft", rootfilenm);
    }

    {
        long numsamp, binnum = 1, numwrote = 0;
        int bufflen, nblk_old, nblk_new, mid_old, mid_new;
        int ii, ind, old_percent = 0, new_percent = 0;
        float mean_old, mean_new, T, dslope = 1.0, norm;
        float *powbuf, powargr, powargi;
        fcomplex *newbuf, *oldbuf, *inbuf1, *inbuf2, *outbuf, *tempzz;

        /* Read the info file, and write a copy */
        {
            char *newinf;
            infodata idata;
            readinf(&idata, rootfilenm);
            numsamp = idata.N;
            T = numsamp * idata.dt;
            newinf = (char *) calloc(strlen(idata.name) + 5, sizeof(char));
            sprintf(newinf, "%s_red", idata.name);
            sprintf(idata.name, "%s", newinf);
            writeinf(&idata);
            free(newinf);
        }

        /* Open files and create arrays */

        infile = chkfopen(argv[1], "rb");
        outfile = chkfopen(outname, "wb");
        bufflen = cmd->startwidth;
        inbuf1 = gen_cvect(cmd->endwidth);
        inbuf2 = gen_cvect(cmd->endwidth);
        outbuf = gen_cvect(cmd->endwidth);
        powbuf = gen_fvect(cmd->endwidth);
        oldbuf = inbuf1;
        newbuf = inbuf2;

        /* Takes care of the DC offset and Nyquist */
        chkfread(oldbuf, sizeof(fcomplex), 1, infile);
        oldbuf[0].r = 1.0;
        oldbuf[0].i = 0.0;
        chkfwrite(oldbuf, sizeof(fcomplex), 1, outfile);
        numwrote += 1;

        // Calculates the first mean
        nblk_old = chkfread(oldbuf, sizeof(fcomplex), bufflen, infile);
        if (nblk_old != bufflen) {
            perror("\nError in rednoise:  number read != bufflen");
            printf("\n");
            exit(-1);
        }
        // Buffer bin of the ~midpoint of the current block
        mid_old = nblk_old / 2;

        // Compute the powers
        for (ii = 0; ii < nblk_old; ii++)
            powbuf[ii] = POWER(oldbuf[ii].r, oldbuf[ii].i);
        mean_old = median(powbuf, nblk_old) / log(2.0);

        // Write out the first half of the normalized block
        // Note that this does *not* include a slope, but since it
        // is only a few bins, that is probably OK.
        norm = invsqrtf(mean_old);
        for (ii = 0; ii < mid_old; ii++) {
            outbuf[ii].r = oldbuf[ii].r * norm;
            outbuf[ii].i = oldbuf[ii].i * norm;
            //printf("  %10ld %4d %.5g\n", ii+numwrote, ii, 1.0/(norm*norm));
        }
        chkfwrite(outbuf, sizeof(fcomplex), mid_old, outfile);
        numwrote += mid_old;

        // This is the Fourier bin index for the next read
        binnum += nblk_old;
        // This updates the length of the median block logarithmically
        bufflen = cmd->startwidth * log(binnum);

        while ((nblk_new = chkfread(newbuf, sizeof(fcomplex), bufflen, infile))) {
            mid_new = nblk_new / 2;
            for (ii = 0; ii < nblk_new; ii++)
                powbuf[ii] = POWER(newbuf[ii].r, newbuf[ii].i);
            mean_new = median(powbuf, nblk_new) / log(2.0);

            // The slope between the last block median and the current median
            dslope = (mean_new - mean_old) / (0.5 * (nblk_old + nblk_new));
            //printf("\n%d %.5g %.5g %.5g\n", nblk_new, mean_old, mean_new, dslope);

            // Correct the last-half of the old block...
            for (ii = 0, ind = mid_old; ind < nblk_old; ii++, ind++) {
                norm = invsqrtf(mean_old + dslope * ii);
                outbuf[ii].r = oldbuf[ind].r * norm;
                outbuf[ii].i = oldbuf[ind].i * norm;
                //printf("  %10ld %4d %.5g\n", ii+numwrote, ii, 1.0/(norm*norm));
            }
            // ...and the first-half of the new block
            for (ind = 0; ind < mid_new; ii++, ind++) {
                norm = invsqrtf(mean_old + dslope * ii);
                outbuf[ii].r = newbuf[ind].r * norm;
                outbuf[ii].i = newbuf[ind].i * norm;
                //printf("  %10ld %4d %.5g\n", ii+numwrote, ii, 1.0/(norm*norm));
            }
            // Write the normalized amplitudes
            chkfwrite(outbuf, sizeof(fcomplex), ii, outfile);
            numwrote += ii;

            // Update the variables and pointers

            binnum += nblk_new;
            if ((float) binnum / T < cmd->endfreq)
                bufflen = cmd->startwidth * log(binnum);
            else
                bufflen = cmd->endwidth;
            SWAP(oldbuf, newbuf);
            nblk_old = nblk_new;
            mean_old = mean_new;
            mid_old = mid_new;

            /* Print percent complete */

            new_percent = (int) 100 *((binnum * 2.0) / numsamp);
            if (new_percent != old_percent) {
                // printf("\rAmount Complete = %d%%", new_percent);
                old_percent = new_percent;
                fflush(stdout);
            }
        }
        // Deal with the last chunk (assume same slope as before)
        for (ii = 0, ind = mid_old; ind < nblk_old; ii++, ind++) {
            norm = invsqrtf(mean_old + dslope * ii);
            outbuf[ii].r = oldbuf[ind].r * norm;
            outbuf[ii].i = oldbuf[ind].r * norm;
        }
        chkfwrite(outbuf, sizeof(fcomplex), nblk_old - mid_old, outfile);
        numwrote += nblk_old - mid_old;

        printf("\nDone.  Rednoise removed from %ld of %ld points.\n\n",
               numwrote, numsamp / 2);
        vect_free(inbuf1);
        vect_free(inbuf2);
        vect_free(powbuf);
        vect_free(outbuf);
    }

    fclose(infile);
    fclose(outfile);
    free(rootfilenm);
    free(outname);
    exit(0);
}
