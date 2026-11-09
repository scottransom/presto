/*
 * coherentsearch: PGPLOT pages of candidate pulse profiles (CoherentSearch.jl's
 * bin/plotting.jl).  Each candidate is refolded at the full harmonic depth
 * with a wide exact kernel, rotated so its peak sits at phase 0.5, and drawn
 * in a 3 x 5 grid per US-Letter portrait page under a header describing the
 * file and the search.
 */
#include <time.h>
#include "coherentsearch.h"
#include "cpgplot.h"

#define CS_PLOT_NCOLS 3
#define CS_PLOT_NROWS 5

/* Render <stem>.ps into anti-aliased <stem>_01.png, _02.png, ... with
 * Ghostscript, as prepfold does: PGPLOT's own PNG driver draws poorly (no
 * anti-aliasing, a black background).  The PNGs are a convenience, so a
 * missing or failing gs is only noted. */
static void ps_to_png(const char *stem)
{
    if (system("gs --version > /dev/null 2>&1") != 0) {
        fprintf(stderr, "coherentsearch: no Ghostscript ('gs'), so no PNGs (the .ps is fine)\n");
        return;
    }
    char *cmd = malloc(2 * strlen(stem) + 200);
    sprintf(cmd, "gs -q -dSAFER -dBATCH -dNOPAUSE -sDEVICE=png16m -r150 "
            "-dTextAlphaBits=4 -dGraphicsAlphaBits=4 -sOutputFile='%s_%%02d.png' '%s.ps' "
            "> /dev/null 2>&1", stem, stem);
    if (system(cmd))
        fprintf(stderr, "coherentsearch: gs failed to render the PNGs (the .ps is fine).  "
                "On Ubuntu this is often the AppArmor 'gs' profile; see FAQ.md.\n");
    free(cmd);
}

/* Circularly shift prof so that its maximum lands on bin n/2 (phase 0.5). */
static void rotate_to_peak(const double *prof, int n, float *out)
{
    int imax = 0;
    for (int i = 1; i < n; i++)
        if (prof[i] > prof[imax])
            imax = i;
    const int shift = n / 2 - imax;
    for (int i = 0; i < n; i++)
        out[((i + shift) % n + n) % n] = (float) prof[i];
}

void cs_plot_candidates(const cs_fft *fft, const cs_candlist *l, int nharms, double threshold,
                        const char *fftname, const char *stem, int png)
{
    if (l->n == 0)
        return;
    char device[4096];
    snprintf(device, sizeof(device), "%s.ps/VCPS", stem);
    if (cpgopen(device) <= 0) {
        fprintf(stderr, "coherentsearch: could not open PGPLOT device '%s'\n", device);
        return;
    }
    cpgask(0);
    cpgpap(8.0, 11.0 / 8.5);
    cpgscr(16, 0.0, 0.0, 0.5);  /* navy */
    const int perpage = CS_PLOT_NCOLS * CS_PLOT_NROWS;
    const long npages = (l->n + perpage - 1) / perpage;
    const char *base = strrchr(fftname, '/') ? strrchr(fftname, '/') + 1 : fftname;
    char date[64];
    time_t now = time(NULL);
    strftime(date, sizeof(date), "%Y-%m-%dT%H:%M:%S", localtime(&now));
    double *prof = malloc((2 * nharms + 2) * sizeof(double));
    float *y = malloc((2 * nharms + 3) * sizeof(float));
    float *x = malloc((2 * nharms + 3) * sizeof(float));
    /* Grid geometry, in normalised device coordinates. */
    const float left = 0.03f, right = 0.99f, bottom = 0.04f, top = 0.90f;
    const float cw = (right - left) / CS_PLOT_NCOLS, ch = (top - bottom) / CS_PLOT_NROWS;
    for (long page = 0; page < npages; page++) {
        cpgpage();
        /* Header */
        cpgsci(1);
        cpgsvp(0.0, 1.0, 0.0, 1.0);
        cpgswin(0.0, 1.0, 0.0, 1.0);
        cpgsch(0.9);
        cpgslw(3);
        char txt[512];
        snprintf(txt, sizeof(txt), "%s      Candidate pulse profiles", base);
        cpgtext(0.03, 0.975, txt);
        cpgslw(1);
        cpgsch(0.75);
        snprintf(txt, sizeof(txt),
                 "object=%s   DM=%.2f   T=%.1f s   N=%lld   nharms=%d   threshold=%.1f", fft->object[0] ? fft->object : "?", fft->dm, fft->T,
                 fft->N, nharms, threshold);
        cpgtext(0.03, 0.955, txt);
        snprintf(txt, sizeof(txt), "%s        Page %ld of %ld", date, page + 1, npages);
        cpgtext(0.03, 0.935, txt);
        const long lo = page * perpage;
        const long hi = lo + perpage < l->n ? lo + perpage : l->n;
        const int lastrow = (int) ((hi - lo + CS_PLOT_NCOLS - 1) / CS_PLOT_NCOLS) - 1;
        for (long idx = lo; idx < hi; idx++) {
            const cs_cand *c = l->c + idx;
            const int slot = (int) (idx - lo), row = slot / CS_PLOT_NCOLS,
                col = slot % CS_PLOT_NCOLS;
            /* Fold at the full requested depth, not the decimated count that
             * found it, so the profile approximates a time-domain fold. */
            const int nb = cs_candidate_profile(fft, c->r, nharms, prof);
            rotate_to_peak(prof, nb, y);
            y[nb] = y[0];                   /* close the cycle at phase 1 */
            float ymin = y[0], ymax = y[0];
            for (int i = 0; i <= nb; i++) {
                x[i] = (float) i / nb;
                ymin = y[i] < ymin ? y[i] : ymin;
                ymax = y[i] > ymax ? y[i] : ymax;
            }
            float pad = 0.05f * (ymax - ymin > 0 ? ymax - ymin : 1.0f);
            /* Room above the (centred) peak for a known pulsar's name. */
            const float toppad = c->note[0] ? 5.0f * pad : pad;
            const float x0 = left + col * cw + 0.025f, x1 = left + (col + 1) * cw - 0.015f;
            const float y1 = top - row * ch - 0.035f, y0 = top - (row + 1) * ch + 0.01f;
            cpgsvp(x0, x1, y0, y1);
            cpgswin(0.0, 1.0, ymin - pad, ymax + toppad);
            cpgsci(1);
            cpgsch(0.7);
            cpgbox(row == lastrow ? "BCNST" : "BCST", 0.5, 5, "BC", 0.0, 0);
            if (row == lastrow)
                cpgmtxt("B", 2.3, 0.5, 0.5, "Pulse phase");
            int k = c->nharm > 0 ? nharms / c->nharm : 0;
            cpgsch(0.6);
            snprintf(txt, sizeof(txt), "#%ld   S/N %.2f   #Harm %d (k=%d)", idx + 1,
                     c->metric, c->nharm, k);
            cpgmtxt("T", 1.8, 0.5, 0.5, txt);
            snprintf(txt, sizeof(txt), "%.10f Hz   %.6f ms", c->freq, 1000.0 / c->freq);
            cpgmtxt("T", 0.6, 0.5, 0.5, txt);
            cpgsci(16);
            cpgslw(2);
            cpgline(nb + 1, x, y);
            cpgslw(1);
            if (c->note[0]) {       /* a known pulsar: its name, bold and red */
                cpgsci(2);
                cpgslw(5);
                cpgsch(0.65);
                cpgmtxt("T", -1.3, 0.03, 0.0, c->note);
                cpgslw(1);
            }
        }
    }
    cpgclos();
    if (png)
        ps_to_png(stem);
    free(prof);
    free(x);
    free(y);
}
