#include "presto.h"
#include "vectors.h"
#include "limits.h"
#include "errno.h"

char bands[NUMBANDS][40] = { "Radio", "IR", "Optical", "UV", "X-ray", "Gamma" };

char scopes[NUMSCOPES][40] =
    { "None (Artificial Data Set)", "Arecibo", "Parkes", "VLA",
    "MMT", "Las Campanas 2.5m", "Mt. Hopkins 48in", "Other"
};

void add_to_inf_epoch(infodata * data, double seconds)
{
    data->mjd_f += seconds / SECPERDAY;
    data->mjd_i += 1 ? (data->mjd_f > 1.0) : 0;
    data->mjd_f -= 1.0 ? (data->mjd_f > 1.0) : 0.0;
}

void read_inf_line_valstr(FILE * infofile, char *valstr, char *errdesc)
{
    char line[250], *sptr = NULL;
    int ii, slen;

    sptr = fgets(line, 250, infofile);
    if (sptr != NULL && sptr[0] != '\n' && 0 != (ii = strlen(sptr))) {
        // Check to see if this is a "standard" .inf line
        // which has an '=' in character 40
        if (ii >= 40 && line[40] == '=') {
            sptr = line + 41;
        } else {
            // Else, look for the first '=' back from the end of the line
            while (--ii >= 0) {
                if (sptr[ii] == '=')
                    break;
            }
            if (ii + 1 == 0) {
                sprintf(line,
                        "Error:  no '=' to separate key/val while looking for '%s' in readinf()\n",
                        errdesc);
                perror(line);
                exit(EXIT_FAILURE);
            }
            sptr = line + ii + 1;
        }
        sptr = remove_whitespace(sptr);
        slen = strlen(sptr);
        if (slen) {
            if ((strcmp(errdesc, "data->name") == 0 && slen > 199) ||
                (strcmp(errdesc, "data->telescope") == 0 && slen > 39) ||
                (strcmp(errdesc, "data->band") == 0 && slen > 39) ||
                (strcmp(errdesc, "data->name") != 0 && slen > 99)) {
                sprintf(line,
                        "Error:  value string is too long (%d char) while looking for '%s' in readinf()\n",
                        slen, errdesc);
                perror(line);
                exit(EXIT_FAILURE);
            }
            strcpy(valstr, sptr);
        } else {
            strcpy(valstr, "Unknown");
        }
        return;
    } else {
        if (feof(infofile)) {
            sprintf(line,
                    "Error:  end-of-file while looking for '%s' in readinf()\n",
                    errdesc);
        } else {
            sprintf(line,
                    "Error:  found blank line while looking for '%s' in readinf()\n",
                    errdesc);
        }
        perror(line);
        exit(EXIT_FAILURE);
    }
    // Should never get here....
}

double chk_str2double(char *instr, char *desc)
{
    char tmp[100], *sptr = instr, *endptr;
    double retval;

    retval = strtod(sptr, &endptr);
    if (retval == 0.0 && endptr == instr) {
        sprintf(tmp,
                "Error:  can not convert '%s' to a double (%s) in chk_str2double()\n",
                instr, desc);
        perror(tmp);
        exit(EXIT_FAILURE);
    }
    return retval;
}

long chk_str2long(char *instr, char *desc)
{
    char tmp[100], *sptr = instr, *endptr;
    long retval;

    errno = 0;
    retval = strtol(sptr, &endptr, 10);
    if ((errno == ERANGE && (retval == LONG_MAX || retval == LONG_MIN))
        || (errno != 0 && retval == 0)) {
        sprintf(tmp,
                "Error:  can not convert '%s' to an int/long (%s) in chk_str2long()\n",
                instr, desc);
        perror(tmp);
        exit(EXIT_FAILURE);
    }
    if (endptr == instr) {
        sprintf(tmp,
                "Error:  No digits were found in '%s' for %s in chk_str2long()\n",
                instr, desc);
        perror(tmp);
        exit(EXIT_FAILURE);
    }
    return retval;
}

void readinf(infodata * data, char *filenm)
{
    char tmp1[100], tmp2[100], tmp3[100], *infofilenm, *sptr;
    int ii, retval, noteslen = 0;
    FILE *infofile;

    infofilenm = malloc(strlen(filenm) + 5);
    sprintf(infofilenm, "%s.inf", filenm);
    infofile = chkfopen(infofilenm, "r");
    free(infofilenm);

    read_inf_line_valstr(infofile, data->name, "data->name");
    read_inf_line_valstr(infofile, data->telescope, "data->telescope");
    /* If not using makedata */
    if (strcmp(data->telescope, scopes[0]) != 0) {
        read_inf_line_valstr(infofile, data->instrument, "data->instrument");
        read_inf_line_valstr(infofile, data->object, "data->object");
        read_inf_line_valstr(infofile, tmp1, "RA string");
        ra_dec_from_string(tmp1, &data->ra_h, &data->ra_m, &data->ra_s);
        read_inf_line_valstr(infofile, tmp1, "DEC string");
        ra_dec_from_string(tmp1, &data->dec_d, &data->dec_m, &data->dec_s);
        read_inf_line_valstr(infofile, data->observer, "data->observer");
        read_inf_line_valstr(infofile, tmp1, "MJD string");
        retval = sscanf(tmp1, "%d.%s", &data->mjd_i, tmp2);
        if (retval != 2) {
            sprintf(tmp3, "Error:  can not parse MJD string '%s' in readinf()'\n",
                    tmp1);
            perror(tmp3);
            exit(EXIT_FAILURE);
        }
        sprintf(tmp3, "0.%s", tmp2);
        data->mjd_f = chk_str2double(tmp3, "data->mjd_f");
        read_inf_line_valstr(infofile, tmp1, "data->bary");
        data->bary = chk_str2long(tmp1, "data->bary");
    } else {
        data->mjd_i = -1;
        strcpy(data->object, "fake pulsar");
    }
    read_inf_line_valstr(infofile, tmp1, "data->N");
    data->N = chk_str2double(tmp1, "data->N");
    read_inf_line_valstr(infofile, tmp1, "data->dt");
    data->dt = chk_str2double(tmp1, "data->dt");
    read_inf_line_valstr(infofile, tmp1, "data->numonoff");
    data->numonoff = chk_str2long(tmp1, "data->numonoff");
    if (data->numonoff) {
        ii = 0;
        do {
            read_inf_line_valstr(infofile, tmp1, "on-off pairs");
            retval = sscanf(tmp1, "%lf %*[ ,] %lf",
                            &data->onoff[ii], &data->onoff[ii + 1]);
            if (retval != 2) {
                sprintf(tmp3,
                        "Error:  can not parse on-off pair (%d) in readinf()\n",
                        ii / 2);
                perror(tmp3);
                exit(EXIT_FAILURE);
            }
            ii += 2;
        } while (data->onoff[ii - 1] < data->N - 1 && ii < 2 * MAXNUMONOFF);
        data->numonoff = ii / 2;
        if (data->numonoff >= MAXNUMONOFF) {
            sprintf(tmp3,
                    "Error:  number of onoff pairs (%d) >= MAXNUMONOFF (%d) in readinf().\n",
                    data->numonoff, MAXNUMONOFF);
            perror(tmp3);
            exit(EXIT_FAILURE);
        }
    } else {
        data->numonoff = 1;
        data->onoff[0] = 0;
        data->onoff[1] = data->N - 1;
    }
    /* If not using makedata */
    if (strcmp(data->telescope, scopes[0]) != 0) {
        read_inf_line_valstr(infofile, data->band, "data->band");
        if (strcmp(data->band, bands[0]) == 0) {
            read_inf_line_valstr(infofile, tmp1, "data->fov");
            data->fov = chk_str2double(tmp1, "data->fov");
            read_inf_line_valstr(infofile, tmp1, "data->dm");
            data->dm = chk_str2double(tmp1, "data->dm");
            read_inf_line_valstr(infofile, tmp1, "data->freq");
            data->freq = chk_str2double(tmp1, "data->freq");
            read_inf_line_valstr(infofile, tmp1, "data->freqband");
            data->freqband = chk_str2double(tmp1, "data->freqband");
            read_inf_line_valstr(infofile, tmp1, "data->num_chan");
            data->num_chan = chk_str2long(tmp1, "data->num_chan");
            read_inf_line_valstr(infofile, tmp1, "data->chan_wid");
            data->chan_wid = chk_str2double(tmp1, "data->chan_wid");
        } else if ((strcmp(data->band, bands[4]) == 0) ||
                   (strcmp(data->band, bands[5]) == 0)) {
            read_inf_line_valstr(infofile, tmp1, "data->fov");
            data->fov = chk_str2double(tmp1, "data->fov");
            read_inf_line_valstr(infofile, tmp1, "data->energy");
            data->energy = chk_str2double(tmp1, "data->energy");
            read_inf_line_valstr(infofile, tmp1, "data->energyband");
            data->energyband = chk_str2double(tmp1, "data->energyband");
        } else {
            read_inf_line_valstr(infofile, data->filt, "data->filt");
            read_inf_line_valstr(infofile, tmp1, "data->fov");
            data->fov = chk_str2double(tmp1, "data->fov");
            read_inf_line_valstr(infofile, tmp1, "data->wavelen");
            data->wavelen = chk_str2double(tmp1, "data->wavelen");
            read_inf_line_valstr(infofile, tmp1, "data->waveband");
            data->waveband = chk_str2double(tmp1, "data->waveband");
        }
    }
    read_inf_line_valstr(infofile, data->analyzer, "data->analyzer");
    // The following is description line for the 'Notes' part
    sptr = fgets(tmp1, 100, infofile);
    // Now read all the notes lines
    while (1) {
        sptr = fgets(tmp1, 100, infofile);
        if (noteslen + strlen(tmp1) > 500)
            break;
        if (sptr) {
            if (noteslen == 0)
                strcpy(data->notes + noteslen, rmlead(tmp1));
            else
                strcpy(data->notes + noteslen, tmp1);
            noteslen += strlen(data->notes);
        } else {
            if (feof(infofile))
                break;
        }
    }
    fclose(infofile);
}


void chk_empty(char *instr)
{
    if (strlen(remove_whitespace(instr)) == 0)
        strcpy(instr, "Unknown");
}


void writeinf(infodata * data)
{
    char tmp1[100], tmp2[100], *infofilenm;
    int itmp, ii;
    FILE *infofile;

    infofilenm = malloc(strlen(data->name) + 5);
    sprintf(infofilenm, "%s.inf", data->name);
    infofile = chkfopen(infofilenm, "w");
    free(infofilenm);

    fprintf(infofile, " Data file name without suffix          =  %s\n", data->name);
    chk_empty(data->telescope);
    fprintf(infofile,
            " Telescope used                         =  %s\n", data->telescope);
    if (strcmp(data->telescope, scopes[0]) != 0) {      /* If using makedata */
        chk_empty(data->instrument);
        fprintf(infofile,
                " Instrument used                        =  %s\n", data->instrument);
        chk_empty(data->object);
        fprintf(infofile,
                " Object being observed                  =  %s\n", data->object);
        ra_dec_to_string(tmp1, data->ra_h, data->ra_m, data->ra_s);
        fprintf(infofile, " J2000 Right Ascension (hh:mm:ss.ssss)  =  %s\n", tmp1);
        ra_dec_to_string(tmp1, data->dec_d, data->dec_m, data->dec_s);
        fprintf(infofile, " J2000 Declination     (dd:mm:ss.ssss)  =  %s\n", tmp1);
        chk_empty(data->observer);
        fprintf(infofile,
                " Data observed by                       =  %s\n", data->observer);
        sprintf(tmp1, "%.15f", data->mjd_f);
        sscanf(tmp1, "%d.%s", &itmp, tmp2);
        fprintf(infofile,
                " Epoch of observation (MJD)             =  %d.%s\n",
                data->mjd_i, tmp2);
        fprintf(infofile,
                " Barycentered?           (1 yes, 0 no)  =  %d\n", data->bary);
    }
    fprintf(infofile,
            " Number of bins in the time series      =  %-11.0f\n", data->N);
    fprintf(infofile, " Width of each time series bin (sec)    =  %.15g\n",
            data->dt);
    fprintf(infofile, " Any breaks in the data? (1 yes, 0 no)  =  %d\n",
            data->numonoff > 1 ? 1 : 0);
    if (data->numonoff > 1) {
        for (ii = 0; ii < data->numonoff; ii++) {
            fprintf(infofile,
                    " On/Off bin pair #%3d                   =  %-11.0f, %-11.0f\n",
                    ii + 1, data->onoff[2 * ii], data->onoff[2 * ii + 1]);
        }
    }
    if (strcmp(data->telescope, scopes[0]) != 0) {      /* If using makedata */
        fprintf(infofile,
                " Type of observation (EM band)          =  %s\n", data->band);
        if (strcmp(data->band, bands[0]) == 0) {
            fprintf(infofile,
                    " Beam diameter (arcsec)                 =  %.0f\n", data->fov);
            fprintf(infofile,
                    " Dispersion measure (cm-3 pc)           =  %.12g\n", data->dm);
            fprintf(infofile,
                    " Central freq of low channel (MHz)      =  %.12g\n",
                    data->freq);
            fprintf(infofile, " Total bandwidth (MHz)                  =  %.12g\n",
                    data->freqband);
            fprintf(infofile, " Number of channels                     =  %d\n",
                    data->num_chan);
            fprintf(infofile, " Channel bandwidth (MHz)                =  %.12g\n",
                    data->chan_wid);
        } else if ((strcmp(data->band, bands[4]) == 0)
                   || (strcmp(data->band, bands[5]) == 0)) {
            fprintf(infofile, " Field-of-view diameter (arcsec)        =  %.2f\n",
                    data->fov);
            fprintf(infofile, " Central energy (kev)                   =  %.1f\n",
                    data->energy);
            fprintf(infofile, " Energy bandpass (kev)                  =  %.1f\n",
                    data->energyband);
        } else {
            chk_empty(data->filt);
            fprintf(infofile,
                    " Photometric filter used                =  %s\n", data->filt);
            fprintf(infofile,
                    " Field-of-view diameter (arcsec)        =  %.2f\n", data->fov);
            fprintf(infofile,
                    " Central wavelength (nm)                =  %.1f\n",
                    data->wavelen);
            fprintf(infofile, " Bandpass (nm)                          =  %.1f\n",
                    data->waveband);
        }
    }
    chk_empty(data->analyzer);
    fprintf(infofile,
            " Data analyzed by                       =  %s\n", data->analyzer);
    fprintf(infofile, " Any additional notes:\n    %s\n\n", data->notes);
    fclose(infofile);
}
