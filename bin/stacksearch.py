#!/usr/bin/env python
import os
from sys import argv, exit, stdout
import numpy as np
import presto.presto as pp
from presto.psr_utils import *

import glob
import time

import matplotlib.pyplot as plt


class cand:
    "Class to hold pulsar candidate info"

    def makeharm(self, harm_num, harm_power, harm_sigma):
        if harm_num > len(self.harmonics):
            self.harmonics.append(harmonic(self.f, harm_power, harm_sigma, harm_num))
        else:
            self.harmonics[harm_num - 1] = harmonic(
                self.f, harm_power, harm_sigma, harm_num
            )

    def __init__(self, freq, cand_power, cand_sigma, file_list=[], T=[]):
        self.f = freq
        self.p = 1000.0 / freq
        self.power = cand_power
        self.sigma = cand_sigma
        self.harmonics = []
        self.makeharm(1, self.power, self.sigma)  # Self is first harmonic
        self.T = T
        self.F_ERROR = (
            20.0 / T
        )  # 20 fourier bins of error -- paulr asks if this should be settable???
        if len(file_list) > 0:
            self.files_dict = self.zoom_powers(
                self.f, file_list
            )  # Dictionary of power:fft_file -- should be helpful for rfi checking
        else:
            self.files_dict = {}

    def zoom_powers(self, freq=-1, fft_files=[]):
        "Returns a dictionary of the maximum power:fft file around a central frequency freq.  fft_files is a sequence of open file objects."
        freq_dict = {}
        if freq < 0:
            freq = self.f
        for fft_file in fft_files:
            # Jump to the right frequency minus one bin (8 bytes per complex sample)
            fft_file.seek(8 * (int(self.T * (freq)) - 1))
            zoomed_powers = pp.spectralpower(
                np.fromfile(fft_file, np.complex64, count=3)
            )
            freq_dict[max(zoomed_powers)] = fft_file.name
        return freq_dict


class harmonic:
    "Class to store and print harmonics"

    def __init__(self, base_freq, harm_power, harm_sigma, harm_num):
        self.base = base_freq
        self.harmnum = harm_num
        self.power = harm_power
        self.sigma = harm_sigma

    def __str__(self):
        return (
            "   "
            + ("%i" % self.harmnum).ljust(8)
            + ("%1.5f" % self.power).rjust(10)
            + ("%1.2f" % self.sigma).rjust(10)
        )


def get_options():
    # Some default and initialization values

    usage_string = """\nStacks a series of fft files (of the same length) and outputs pulsar candidates.
Can also read in a single file of prestacked powers.

Usage: stacksearch.py [options...]  files...
Options:
    -h                         Print this help and exit
    -n nfreqs                  Number of frequencies to read from fft files. Default: determine from first input file .inf
    -o powfile                 Name of file to store stacked powers. Default: none
    -ncands num_candidates     Maximum number of candidates to gather. Default: -1 (no limit)
    -candfile candfile_name    Name for text file of candidates. Default: name of powfile .cands or 'candidates'
    -sigma sigma_cutoff        Lowest sigma that will be considered a viable candidate. Default: 4.5
    -noadd                     Use -noadd with a powers file as last argument to search prestacked data.
                                   Automatically set to true if there is only one input file and its extension is .pow"""

    nfreqs = 0  # number of points to read from fft file. Equal to half the number of samples in the corresponding time series
    powfile = ""
    ncands = -1
    noadd = 0
    SIGMA_CUTOFF = 3.0
    candfile_name = ""
    fft_list = []
    debug = False

    # Parse command line
    if len(argv) == 1:
        print(usage_string)
        exit(0)
    else:
        x = 1
        while x < len(argv):
            if argv[x].startswith("-"):
                if argv[x] == "-n":
                    nfreqs = int(argv[x + 1])
                    x += 2
                elif argv[x] == "-o":
                    powfile = argv[x + 1]
                    x += 2
                elif argv[x] == "-ncands":
                    ncands = int(argv[x + 1])
                    x += 2
                elif argv[x] == "-noadd":
                    noadd = 1
                    x += 1
                elif argv[x] == "-candfile":
                    candfile_name = argv[x + 1]
                    x += 2
                elif argv[x] == "-sigma":
                    SIGMA_CUTOFF = float(argv[x + 1])
                    x += 2
                elif argv[x] == "-h":
                    print(usage_string)
                    exit(0)
                elif argv[x] == "-d":
                    debug = True
                    x += 1
                else:
                    print("Unrecognized option: %s" % argv[x])
                    print(usage_string)
                    exit(0)

            # Store remaining (required) arguments
            else:
                input_files = argv[x:]
                break

        # Be smart about running on a .pow file
        if len(input_files) == 1 and input_files[0].endswith(".pow"):
            noadd = 1

        nfiles = len(input_files)
        print("Reading %d files" % nfiles)
        # Read necessary data from .inf file
        inf_file = open(".".join(argv[x].split(".")[:-1]) + ".inf", "r")
        inf_file_lines = inf_file.readlines()
        DM = 0.0  # Default to 0 DM, if not found in the .inf file
        for line in inf_file_lines:
            if line.startswith(" Width of each time series bin"):
                SAMPLE_TIME = float(line.split()[-1])
                print("Sample time is", SAMPLE_TIME, "s")
            if line.startswith(" Dispersion measure (cm-3 pc)"):
                DM = float(line.split()[-1])
                print("Files dedisperded at DM=%.2f" % DM)
            if nfreqs == 0:
                if line.startswith(" Number of bins in the time series"):
                    nfreqs = int(line.split()[-1]) // 2
                    print("Number of frequencies in fft files is %d" % nfreqs)
            if noadd:
                if line.startswith("    Number of files stacked:"):
                    nfiles = int(line.split()[-1])
                if line.startswith("    File:"):
                    fft_list += [line.split()[-1]]

        # Establish the remaining default values
        if not noadd:
            fft_list = input_files
        inf_file.close()
        T = nfreqs * SAMPLE_TIME * 2
        F_ERROR = 20 / T  # 20 fourier bins of error
        if candfile_name == "":  # Set name for candidate list output file
            if noadd:
                candfile_name = ".".join(input_files[0].split(".")[:-1])
            elif powfile != "":
                candfile_name = powfile
            else:
                candfile_name = "candidates"

        print(f"candfile_name = {candfile_name}")
        # make dictionary of all the relevant variables
        options = {}
        options["nfreqs"] = nfreqs  # number of frequencies in the fft's
        options["powfile"] = powfile  # name of file to write stored powers to
        options["ncands"] = ncands  # max number of candidates
        options["noadd"] = noadd  # flag for using prestacked powers
        options["candfile_name"] = (
            candfile_name  # name of file to write candidate list to
        )
        options["SIGMA_CUTOFF"] = SIGMA_CUTOFF  # min sigma
        options["T"] = T  # length of time series (s)
        options["F_ERROR"] = F_ERROR  # frequency error margin
        options["fft_list"] = fft_list  # names of stacked fft files
        options["input_files"] = input_files  # names of input files from command line
        options["inf_file_lines"] = inf_file_lines  # from readlines() on .inf file
        options["nfiles"] = nfiles  # number of stacked fft files
        options["DM"] = DM  # DM of the fft files
        options["DEBUG"] = debug  # Turn on debug plotting
        print("T=", options["T"])
        return options


class Candidates:  # This class holds the power spectrum and candidate list
    def __init__(self, options):
        # initialize the values of options that get used a lot
        powfile = options["powfile"]
        nfreqs = options["nfreqs"]
        self.nfreqs = nfreqs
        noadd = options["noadd"]
        SIGMA_CUTOFF = options["SIGMA_CUTOFF"]
        F_ERROR = options["F_ERROR"]
        self.F_ERROR = F_ERROR
        ncands = options["ncands"]
        self.DEBUG = options["DEBUG"]
        self.candfile_name = options["candfile_name"]

        # First, read in or create stacked power spectrum
        if noadd:
            print("Reading pre-added powers from", options["input_files"][0], "... ")
            self.powers = np.fromfile(
                options["input_files"][0], np.float32, count=nfreqs
            )
            print("Done.")
            print("Number of files stacked: " + str(options["nfiles"]))
            print("The stacked files were:")
            for file_name in options["fft_list"]:
                print(file_name)
        else:
            self.powers = self.stack_files(options["input_files"])  # Get power spectrum
            self.sumharms()  # Harmonic summing
            self.write_stacked_powers(
                self.powers, powfile, options["inf_file_lines"], options["input_files"]
            )

        # Now we extract just the powers above the sigma threshold
        print("Culling the spectrum powers above the sigma threshold:")
        # power_cutoff = powersum_at_sigma(SIGMA_CUTOFF, options['nfiles']) #THIS SEEMS TOO LARGE...
        power_cutoff = powersum_at_sigma(
            SIGMA_CUTOFF, 1
        )  # paulr - seems to be ignoring the number of summed power spectrum bins
        print("Power cut-off is %.2f" % power_cutoff)
        self.hi_pow_args = np.compress(
            self.powers > power_cutoff, np.arange(len(self.powers))
        )
        self.hi_pows = np.take(self.powers, self.hi_pow_args)
        print("\n", len(self.hi_pows), "peaks above sigma threshold of", SIGMA_CUTOFF)

        # open the fft files so we can quickly check powers with zoom_powers
        fft_file_list = [open(name, "rb") for name in options["fft_list"]]

        # Now we build the candidate list from the culled powers
        self.candlist = []
        try:
            for culled_index in np.argsort(self.hi_pows)[::-1]:
                dupe = 0
                maxbin = self.hi_pow_args[culled_index]
                maxfreq = (
                    maxbin / options["T"]
                ) / 8  # Divide by n here to offset the effect of harmonic summing by n...
                maxpow = float(
                    self.hi_pows[culled_index]
                )  # WHY IS THIS NEEDED???????????????
                # sigma = presto.candidate_sigma(maxpow, options['nfiles'], 1)
                sigma = pp.candidate_sigma(
                    maxpow, 1, 1
                )  # paulr - seems to be ignoring the number of summed power spectrum bins
                if (
                    maxfreq > 25.0 and maxfreq < 2000
                ):  # paulr - These should not be hard coded!
                    # This loops through the previous candidates in the list and checks for duplicates within an error margin
                    for candidate in self.candlist:
                        if abs(candidate.f - maxfreq) < F_ERROR:
                            dupe = 1
                            break
                    if not dupe:
                        self.candlist += [
                            cand(maxfreq, maxpow, sigma, fft_file_list, options["T"])
                        ]

                print(
                    "\rCandidates found:",
                    len(self.candlist),
                    "Most recent sigma: ",
                    sigma,
                    end=" ",
                )
                stdout.flush()
                if len(self.candlist) > ncands and ncands > 0:
                    break
        except KeyboardInterrupt:
            print("\n\nKeyboardInterrupt! Finishing...")
        for file in fft_file_list:
            file.close()

    def find_harms(self):
        # Sorting candlist by frequency and finds all harmonics, adds them to their base candidates harmlist, and deletes them from candlist
        ii = 0  # Index of current base candidate
        n = 2  # Harmonic being searched
        self.candlist.sort(
            key=lambda x: x.f
        )  # candidate list is now sorted by frequency.  makes searching for harms much simpler
        print("Searching for harmonics... ")
        while "TRUE":
            if ii >= len(self.candlist):
                break
            freq = self.candlist[ii].f  # Base frequency
            temp_harmlist = (
                []
            )  # List to hold all the harmonics found each time we loop, so we can pick the one with the greatest power
            jj = ii + 1  # Index of harmonic candidate
            while "TRUE":  # Search through higher frequencies in the list
                if jj >= len(self.candlist):
                    break
                if (
                    abs(self.candlist[jj].f - n * freq) <= n * self.F_ERROR
                ):  # Do I need to multiply error by n??
                    temp_harmlist += [self.candlist[jj]]
                    del self.candlist[jj]
                else:
                    jj += 1
            if len(temp_harmlist) > 0:
                harmonic = temp_harmlist[
                    np.argmax([candidate.power for candidate in temp_harmlist])
                ]  # Pick the harmonic with the highest power
                self.candlist[ii].makeharm(n, harmonic.power, harmonic.sigma)
                n += 1
            else:  # Move on to next base freq if a harmonic isn't found
                ii += 1
                # n=2
        print("Done.")

    def print_cands(self):
        print(f"Writing candidates to {self.candfile_name}.cands")
        candfile = open(self.candfile_name + ".cands", "w")
        candfile.write(
            "#Cand num:      power:     sigma:     period(ms):     frequency(Hz):        DM:\n"
        )
        z = 0
        for y in np.argsort([candidate.power for candidate in self.candlist])[::-1]:
            candfile.write(
                "Cand   "
                + str(z).ljust(5)
                + "  %5.2f     %3.1f     %.5f           %.10f          %.2f\n"
                % (
                    self.candlist[y].power,
                    self.candlist[y].sigma,
                    self.candlist[y].p,
                    self.candlist[y].f,
                    options["DM"],
                )
            )
            z += 1
        candfile.close()

    def normalize_power(self, pow, nfreqs, segwidth):
        # This divides the summed power spectrum by the local median power (within "swidth" bins)
        steps = np.arange(0, nfreqs, segwidth)
        for step in steps:
            a = step
            b = step + segwidth
            # Last chunk may be shorter than segwidth
            piece = pow[a:b]
            median = np.median(piece)
            pow[a:b] /= median

    def stack_files(self, fft_files):
        print(
            "Reads files from the list of files and stacks n_frequencies points into an array"
        )
        # Empty array to hold stacked powers
        stackedpow_norm = np.zeros(self.nfreqs, dtype=np.float32)
        nstack = len(fft_files)
        print("Stacking %i files of %i frequencies:" % (nstack, self.nfreqs))

        for fft_file in fft_files:
            print("Stacking", fft_file, "...", end=" ")
            stdout.flush()
            ft = np.fromfile(fft_file, np.complex64, count=self.nfreqs)
            currentpows = pp.spectralpower(ft)
            currentpows[0] = 0.0  # Zap DC bin since it is often huge
            # Normalize powers to 1
            if self.DEBUG:
                print(f"Prenorm mean {currentpows.mean()}")
            self.normalize_power(currentpows, self.nfreqs, self.nfreqs // 128)
            if self.DEBUG:
                print(f"Postnorm mean {currentpows.mean()}")
            stackedpow_norm = np.add(stackedpow_norm, currentpows, stackedpow_norm)
            del currentpows
            print("Done.")

        # Divide by number of power spectra summed
        stackedpow_norm /= nstack
        if self.DEBUG:
            fig, ax = plt.subplots()
            ax.hist(stackedpow_norm, range=(0.0, 20.0), bins=200)
            plt.show()

        return stackedpow_norm

    def write_stacked_powers(self, stackedpow, powfile, inf_file_lines, input_files):
        print("Writing stacked powers...")
        outfile = open(powfile + ".pow", "wb")
        stackedpow.tofile(outfile)
        outfile.close()
        inf_file = open(powfile + ".inf", "w")
        inf_file_lines[0] = (
            " Data file name without suffix          =  "
            + ".".join(input_files[0].split(".")[:-1])
            + "\n"
        )
        inf_file_lines += "    Number of files stacked: " + str(len(input_files))
        for file_name in input_files:
            inf_file_lines += "\n    File:  " + file_name
        inf_file.writelines(inf_file_lines)
        inf_file.close()
        print("Done.")

    def sumharms(self):
        print("Performing harmonic summing...")
        self.powers = (
            self.powers
            + zoomharm(2, self.powers)
            + zoomharm(4, self.powers)
            + zoomharm(8, self.powers)
        )
        print("Done.")


def zoomharm(harmonic, arr):
    return np.resize(
        np.ravel(
            np.transpose(
                np.resize(arr[: len(arr) // harmonic], (harmonic, len(arr) // harmonic))
            )
        )[(harmonic - 1) // 2 :],
        (len(arr),),
    )


# Launch the stacksearch:
if __name__ == "__main__":

    options = get_options()
    cands = Candidates(options)
    cands.find_harms()
    cands.print_cands()
