#!/usr/bin/env python
from numpy import *
from Pgplot import *

class observation:
    def __init__(self, dt, f_ctr, BW, numchan):
        # dt in sec, f_ctr and in MHz
        self.dt = dt
        self.f_ctr = f_ctr
        self.BW = BW
        self.numchan = numchan
        self.chanwidth = BW/numchan
    def guess_dDM(self, DM):
        """
        guess_dDM(self, DM):
            Choose a reasonable dDM by setting the maximum smearing across the
                'BW' to equal the sampling time 'dt'.
        """
        return self.dt*0.0001205*self.f_ctr**3.0/(0.5*self.BW)

class dedisp_method:
    def __init__(self, obs, downsamp, loDM, hiDM, dDM, numDMs=0,
                 numsub=0, smearfact=2.0):
        self.obs = obs
        self.downsamp = downsamp
        self.loDM = loDM
        self.dDM = dDM
        self.numsub = numsub
        self.BW_smearing = BW_smear(dDM, self.obs.BW, self.obs.f_ctr)
        self.numprepsub = 0
        if (numsub):  # Calculate the maximum subband smearing we can handle
            DMs_per_prepsub = 2
            while(1):
                next_dsubDM = (DMs_per_prepsub+2) * dDM
                next_ss = subband_smear(next_dsubDM, numsub, self.obs.BW, self.obs.f_ctr)
                # The 0.8 is a small fudge factor to make sure that the subband
                # smearing is always the smallest contribution
                if (next_ss > 0.8*min(self.BW_smearing, 1000.0*obs.dt*downsamp)):
                    self.dsubDM = DMs_per_prepsub*dDM
                    self.DMs_per_prepsub = DMs_per_prepsub
                    self.sub_smearing = subband_smear(self.dsubDM, numsub,
                                                      self.obs.BW, self.obs.f_ctr)
                    break
                DMs_per_prepsub += 2
        else:
            self.dsubDM = dDM
            self.sub_smearing = subband_smear(self.dsubDM, numsub, self.obs.BW, self.obs.f_ctr)
        # Calculate the nominal DM to move to the next method
        cross_DM = self.DM_for_smearfact(smearfact)
        if (cross_DM > hiDM):
            cross_DM = hiDM
        if (numDMs==0):
            self.numDMs = int(ceil((cross_DM-loDM)/dDM))
            if (numsub):
                self.numprepsub = int(ceil(self.numDMs*dDM / self.dsubDM))
                self.numDMs = self.numprepsub * DMs_per_prepsub
        else:
            self.numDMs = numDMs
        self.hiDM = loDM + self.numDMs*dDM
        self.DMs = arange(self.numDMs, dtype='d')*dDM + loDM
    def chan_smear(self, DM):
        """
        Return the smearing (in ms) in each channel at the specified DM
        """
        try:
            DM = where(DM==0.0, self.dDM/2.0, DM)
        except TypeError:
            if (DM==0.0): DM = self.dDM/2.0
        return dm_smear(DM, self.obs.chanwidth, self.obs.f_ctr)
    def total_smear(self, DM):
        """
        Return the total smearing in ms due to the sampling rate,
        the smearing over each channel, the smearing over each subband
        (if numsub > 0) and the smearing over the full BW assuming the
        worst-case DM error.
        """
        return sqrt((1000.0*self.obs.dt)**2.0 +
                    (1000.0*self.obs.dt*self.downsamp)**2.0 +
                    self.BW_smearing**2.0 +
                    self.sub_smearing**2.0 +
                    self.chan_smear(DM)**2.0)
    def DM_for_smearfact(self, smearfact):
        """
        Return the DM where the smearing in a single channel is a factor smearfact
        larger than all the other smaring causes combined.
        """
        other_smear = sqrt((1000.0*self.obs.dt)**2.0 +
                           (1000.0*self.obs.dt*self.downsamp)**2.0 +
                           self.BW_smearing**2.0 +
                           self.sub_smearing**2.0)
        return smearfact*0.001*other_smear/self.obs.chanwidth*0.0001205*self.obs.f_ctr**3.0
    def DM_for_newparams(self, dDM, downsamp):
        """
        Return the DM where the smearing in a single channel is causes the same smearing
        as the effects of the new dosnsampling rate and dDM.
        """
        other_smear = sqrt((1000.0*self.obs.dt)**2.0 +
                           (1000.0*self.obs.dt*downsamp)**2.0 +
                           BW_smear(dDM, self.obs.BW, self.obs.f_ctr)**2.0 +
                           self.sub_smearing**2.0)
        return 0.001*other_smear/self.obs.chanwidth*0.0001205*self.obs.f_ctr**3.0
    def plot(self, work_fract):
        DMspan = self.DMs[-1]-self.DMs[0]
        loDM  = self.DMs[0]  + DMspan*0.02
        hiDM  = self.DMs[-1] - DMspan*0.02
        midDM = self.DMs[0]  + DMspan*0.5
        dt_ms = 1000.0*self.obs.dt*self.downsamp
        plotxy(log10(self.total_smear(self.DMs)), self.DMs, width=4)
        ppgplot.pgslw(1)
        ppgplot.pgsch(0.8)
        ppgplot.pgptxt(midDM, log10(1.1*self.total_smear(midDM)), 90.0, 0.0,
                       "%d (%.1f%%)" % (self.numDMs, 100.0*work_fract))
        # Sample time
        plotxy(log10(zeros(self.numDMs)+dt_ms), self.DMs, width=1, color='green')
        ppgplot.pgptxt(loDM, log10(0.85*dt_ms), 0.0, 0.0,
                       "%g" % dt_ms)
        # DM stepsize smearing
        plotxy(log10(zeros(self.numDMs)+self.BW_smearing), self.DMs, width=1, color='red')
        ppgplot.pgptxt(hiDM, log10(0.85*self.BW_smearing), 0.0, 1.0,
                       "%g" % self.dDM)
        # channel smearing
        plotxy(log10(self.chan_smear(self.DMs)), self.DMs, width=1, color='blue')
        # subband smearing
        if (self.numsub):
            plotxy(log10(zeros(self.numDMs)+self.sub_smearing),
                   self.DMs, width=1, color='purple')
            ppgplot.pgptxt(midDM, log10(0.85*self.sub_smearing), 0.0, 0.5,
                           "%g (%d)" % (self.dsubDM, self.numprepsub))
        ppgplot.pgsci(1)
        ppgplot.pgsch(1.0)
    def __str__(self):
        if (self.numsub):
            return "%9.3f  %9.3f  %6.2f    %4d  %6.2f  %6d  %6d  %6d " % \
                   (self.loDM, self.hiDM, self.dDM, self.downsamp, self.dsubDM,
                    self.numDMs, self.DMs_per_prepsub, self.numprepsub)
        else:
            return "%9.3f  %9.3f  %6.2f    %4d  %6d" % \
                   (self.loDM, self.hiDM, self.dDM, self.downsamp, self.numDMs)
        
def dm_smear(DM, BW, f_ctr):
    """
    dm_smear(DM, BW, f_ctr):
        Return the smearing in ms caused by a 'DM' over a bandwidth
        of 'BW' MHz centered at 'f_ctr' MHz.
    """
    return 1000.0*DM*BW/(0.0001205*f_ctr**3.0)

def BW_smear(DMstep, BW, f_ctr):
    """
    BW_smear(DMstep, BW, f_ctr):
        Return the smearing in ms caused by a search using a DM stepsize of
        'DMstep' over a bandwidth of 'BW' MHz centered at 'f_ctr' MHz.
    """
    maxDMerror = 0.5*DMstep
    return dm_smear(maxDMerror, BW, f_ctr)

def guess_DMstep(DM, dt, BW, f_ctr):
    """
    guess_DMstep(DM, dt, BW, f_ctr):
        Choose a reasonable DMstep by setting the maximum smearing across the
        'BW' to equal the sampling time 'dt'.
    """
    return dt*0.0001205*f_ctr**3.0/(0.5*BW)

def subband_smear(subDMstep, numsub, BW, f_ctr):
    """
    subband_smear(subDMstep, numsub, BW, f_ctr):
        Return the smearing in ms caused by a search using a subband
        DM stepsize of 'subDMstep' over a total bandwidth of 'BW' MHz
        centered at 'f_ctr' MHz, and having numsub subbands.
    """
    if (numsub==0): return 0.0
    subBW = BW/numsub
    maxsubDMerror = 0.5*subDMstep
    return dm_smear(maxsubDMerror, subBW, f_ctr)

def total_smear(DM, DMstep, dt, f_ctr, BW, numchan, subDMstep, numsub=0):
    """
    total_smear(DM, DMstep, dt, f_ctr, BW, numchan, subDMstep, numsub=0):
        Return the total smearing in ms due to the sampling rate,
        the smearing over each channel, the smearing over each subband
        (if numsub > 0) and the smearing over the full BW assuming the
        worst-case DM error.
    """
    return sqrt(2 * (1000.0*dt)**2.0 +
                dm_smear(DM, BW/numchan, f_ctr)**2.0 +
                subband_smear(subDMstep, numsub, BW, f_ctr)**2.0 + 
                BW_smear(DMstep, BW, f_ctr)**2.0)

def dm_steps(loDM, hiDM, obs, numsub=0, ok_smearing=0.0, device="/XWIN"):
    """
    dm_steps(loDM, hiDM, obs, numsub=0, ok_smearing=0.0):
        Return the optimal DM stepsizes (and subband DM stepsizes if
        numsub>0) to keep the total smearing below 'ok_smearing' (in ms),
        for the DMs between loDM and hiDM.  If 'ok_smearing'=0.0, then
        use the best values based only on the data.
    """
    # Allowable DM stepsizes
    allow_dDMs = [0.01, 0.02, 0.03, 0.05, 0.1, 0.2, 0.3, 0.5, 1.0,
                  2.0, 3.0, 5.0, 10.0, 20.0, 30.0, 50.0, 100.0, 200.0, 300.0]
    # Allowable number of downsampling factors
    allow_downsamps = [1, 2, 4, 8, 16, 32, 64]

    # Initial values
    index_downsamps = index_dDMs = 0
    downsamp = allow_downsamps[index_downsamps]
    dDM = allow_dDMs[index_dDMs]
    dtms = 1000.0*obs.dt
    
    # Fudge factor that "softens" the boundary defining
    # if 2 time scales are equal or not
    ff = 1.2
    
    # This is the array that will hold the de-dispersion plans
    methods = []

    # Minimum possible smearing
    min_tot_smearing = total_smear(loDM+0.5*dDM, dDM, obs.dt, obs.f_ctr,
                                   obs.BW, obs.numchan, allow_dDMs[0], 0)
    # Minimum channel smearing
    min_chan_smearing = dm_smear(loDM+0.5*dDM, obs.chanwidth, obs.f_ctr)
    # Minimum smearing across the obs.BW
    min_BW_smearing = BW_smear(dDM, obs.BW, obs.f_ctr)

    print
    print "Minimum total smearing     : %.3g ms" % min_tot_smearing
    print "--------------------------------------------"
    print "Minimum channel smearing   : %.3g ms" % min_chan_smearing
    print "Minimum smearing across BW : %.3g ms" % min_BW_smearing
    print "Minimum sample time        : %.3g ms" % dtms
    print

    ok_smearing = max([ok_smearing, min_chan_smearing, min_BW_smearing, dtms])
    print "Setting the new 'best' resolution to : %.3g ms" % ok_smearing

    # See if the data is too high time resolution for our needs
    if (ff*min_chan_smearing > dtms or
        ok_smearing > dtms):
        if (ok_smearing > ff*min_chan_smearing):
            print "   Note: ok_smearing > dt (i.e. data is higher resolution than needed)"
            okval = ok_smearing
        else:
            print "   Note: min_chan_smearing > dt (i.e. data is higher resolution than needed)"
            okval = ff*min_chan_smearing

        while (dtms*allow_downsamps[index_downsamps+1] < okval):
            index_downsamps += 1    
        downsamp = allow_downsamps[index_downsamps]
        print "         New dt is %d x %.12g ms = %.12g ms" % \
              (downsamp, dtms, dtms*downsamp)

    # Calculate the appropriate initial dDM 
    dDM = guess_DMstep(loDM, obs.dt*downsamp, obs.BW, obs.f_ctr)
    print "Best guess for optimal initial dDM is %.3f" % dDM
    while (allow_dDMs[index_dDMs+1] < ff*dDM):
        index_dDMs += 1

    # Create the first method
    methods = [dedisp_method(obs, downsamp, loDM, hiDM,
                             allow_dDMs[index_dDMs], numsub=numsub)]
    numDMs = [methods[-1].numDMs]
   
    # Calculate the next methods
    while(methods[-1].hiDM < hiDM):

        # Determine the new downsample factor
        index_downsamps += 1
        downsamp = allow_downsamps[index_downsamps]
        eff_dt = dtms*downsamp

        # Determine the new DM step
        while (BW_smear(allow_dDMs[index_dDMs+1], obs.BW, obs.f_ctr) < ff*eff_dt):
            index_dDMs += 1
        dDM = allow_dDMs[index_dDMs]

        # Get the next method
        methods.append(dedisp_method(obs, downsamp, methods[-1].hiDM,
                                     hiDM, dDM, numsub=numsub))
        numDMs.append(methods[-1].numDMs)
            
    # Calculate the DMs to search and the smearing at each
    total_numDMs = sum(numDMs)
    DMs = zeros(total_numDMs, dtype='d')
    total_smears = zeros(total_numDMs, dtype='d')

    # Calculate the DMs and optimal smearing for all the DMs
    for ii, offset in enumerate(add.accumulate([0]+numDMs[:-1])):
        DMs[offset:offset+numDMs[ii]] = methods[ii].DMs
        total_smears[offset:offset+numDMs[ii]] = methods[ii].total_smear(methods[ii].DMs)

    # Calculate the predicted amount of time that will be spent in searching
    # this batch of DMs as a fraction of the total
    work_fracts = [meth.numDMs/float(meth.downsamp) for meth in methods]
    work_fracts = asarray(work_fracts)/sum(work_fracts)

    # The optimal smearing
    tot_smear = total_smear(DMs, allow_dDMs[0], obs.dt, obs.f_ctr,
                            obs.BW, obs.numchan, allow_dDMs[0], 0)
    # Plot them
    plotxy(log10(tot_smear), DMs, color='orange', logy=1, rangex=[loDM, hiDM],
           rangey=[log10(0.3*min(tot_smear)), log10(2.5*max(tot_smear))],
           labx="Dispersion Measure (pc/cm\u3\d)", laby="Smearing (ms)",
           device=device)
    ppgplot.pgsch(1.1)
    ppgplot.pgsci(1)
    if (numsub):
        ppgplot.pgmtxt("t", 1.5, 0.6/10.0, 0.5, "\(2156)\dctr\u = %g MHz" % obs.f_ctr)
        if (dtms < 0.1):
            ppgplot.pgmtxt("t", 1.5, 2.8/10.0, 0.5, "dt = %g \gms" % (dtms*1000))
        else:
            ppgplot.pgmtxt("t", 1.5, 2.8/10.0, 0.5, "dt = %g ms" % dtms)
        ppgplot.pgmtxt("t", 1.5, 5.0/10.0, 0.5, "BW = %g MHz" % obs.BW)
        ppgplot.pgmtxt("t", 1.5, 7.2/10.0, 0.5, "N\dchan\u = %d" % obs.numchan)
        ppgplot.pgmtxt("t", 1.5, 9.4/10.0, 0.5, "N\dsub\u = %d" % numsub)
    else:
        ppgplot.pgmtxt("t", 1.5, 1.0/8.0, 0.5, "\(2156)\dctr\u = %g MHz" % obs.f_ctr)
        if (dtms < 0.1):
            ppgplot.pgmtxt("t", 1.5, 3.0/8.0, 0.5, "dt = %g \gms" % (dtms*1000))
        else:
            ppgplot.pgmtxt("t", 1.5, 3.0/8.0, 0.5, "dt = %g ms" % dtms)
        ppgplot.pgmtxt("t", 1.5, 5.0/8.0, 0.5, "BW = %g MHz" % obs.BW)
        ppgplot.pgmtxt("t", 1.5, 7.0/8.0, 0.5, "N\dchan\u = %d" % obs.numchan)

    ppgplot.pgsch(1.0)
    dy = -1.4
    ppgplot.pgsci(1)
    ppgplot.pgmtxt("b", 6*dy, 0.97, 1.0, "Total Smearing")
    ppgplot.pgsci(8)
    ppgplot.pgmtxt("b", 5*dy, 0.97, 1.0, "Optimal Smearing")
    ppgplot.pgsci(4)
    ppgplot.pgmtxt("b", 4*dy, 0.97, 1.0, "Channel Smearing")
    ppgplot.pgsci(3)
    ppgplot.pgmtxt("b", 3*dy, 0.97, 1.0, "Sample Time (ms)")
    ppgplot.pgsci(2)
    ppgplot.pgmtxt("b", 2*dy, 0.97, 1.0, "DM Stepsize Smearing")
    if (numsub):
        ppgplot.pgsci(12)
        ppgplot.pgmtxt("b", 1*dy, 0.97, 1.0, "Subband Stepsize Smearing (# passes)")
    ppgplot.pgsci(11)
    
    if (numsub):
        print "\n  Low DM    High DM     dDM  DownSamp  dsubDM   #DMs  DMs/call  calls  WorkFract"
    else:
        print "\n  Low DM    High DM     dDM  DownSamp   #DMs  WorkFract"
    for method, fract in zip(methods, work_fracts):
        print method, "  %.4g" % fract
        method.plot(fract)
    print "\n\n"
    closeplot()
    
def usage():
    print """
usage:  DDplan.py [options]
  [-h, --help]                    : Display this help
  [-o outfile, --outfile=outfile] : Output .eps plot file (default is xwin)
  [-l loDM, --loDM=loDM]          : Low DM to search   (default = 0 pc cm-3)
  [-d hiDM, --hiDM=HIDM]          : High DM to search  (default = 1000 pc cm-3)
  [-f fctr, --fctr=fctr]          : Center frequency   (default = 1400MHz)
  [-b BW, --bw=bandwidth]         : Bandwidth in MHz   (default = 300MHz)
  [-n #chan, --numchan=#chan]     : Number of channels (default = 1024)
  [-t dt, --dt=dt]                : Sample time (s)    (default = 0.000064 s)
  [-s subbands, --subbands=nsub]  : Number of subbands (default = #chan) 
  [-r resolution, --res=res]      : Acceptable time resolution (ms)
  The program generates a good plan for de-dispersing raw data.  It
  trades a small amount of sensitivity in order to save computation costs.

"""    

if __name__=='__main__':
    import getopt, sys

    try:
        opts, args = getopt.getopt(sys.argv[1:], "ho:l:d:f:b:n:t:s:r:",
                                   ["help", "output=", "loDM=", "hiDM=",
                                    "fctr=", "bw=", "numchan=", "dt=",
                                    "subbands=", "res="])

    except getopt.GetoptError:
        # print help information and exit:
        usage()
        sys.exit(2)
    if len(sys.argv)==1:
        usage()
        sys.exit(2)
    # The defaults are close to the future ALFA survey
    loDM, hiDM = 0.0, 1000.0
    fctr = 1400.0
    BW = 300.0
    numchan = 1024
    numsubbands = 0
    dt = 0.000064
    ok_smearing = 0.0
    device = "/xwin"

    for o, a in opts:
        if o in ("-h", "--help"):
            usage()
            sys.exit()
        if o in ("-o", "--output"):
            device = a
            if not (device.endswith(".eps") or device.endswith(".ps")):
                device += ".eps"
            if not (device.endswith("/CPS") or device.endswith("/VCPS") or
                    device.endswith("/cps") or device.endswith("/vcps")):
                device += "/CPS"
        if o in ("-l", "--loDM"):
            loDM = float(a)
        if o in ("-d", "--hiDM"):
            hiDM = float(a)
        if o in ("-f", "--fctr"):
            fctr = float(a)
        if o in ("-b", "--bw"):
            BW = float(a)
        if o in ("-n", "--numchan"):
            numchan = int(a)
        if o in ("-t", "--dt"):
            dt = float(a)
        if o in ("-s", "--subbands"):
            numsubbands = int(a)
        if o in ("-r", "--res"):
            ok_smearing = float(a)

    # The following is an instance of an "observation" class
    obs = observation(dt, fctr, BW, numchan)
    # The following function creates the de-dispersion plan
    # The ok_smearing values is optional and allows you to raise the floor
    # and provide a level of smearing that you are willing to accept (in ms)
    dm_steps(loDM, hiDM, obs, numsub=numsubbands,
             ok_smearing=ok_smearing, device=device)
    
    # The following is an instance of an "observation" class
    # Here's one for a "best" resolution GBT search using the SPIGOT
    # Check out how many DMs you need!  Cool.  ;-)
    #                    dt     f_ctr   BW  numchan
    #obs = observation(0.00008192, 350.0, 50.0, 2048)
    #dm_steps(0.0, 500.0, obs, ok_smearing=0.3) # Create an X-window is the default
    #dm_steps(0.0, 500.0, obs, ok_smearing=0.3, device="GBT_350_DD.ps/CPS")

    # Example for Parkes survey
    #obs = observation(0.000250, 1374.0, 288.0, 96)
    #dm_steps(0.0, 500.0, obs) # Create an X-window is the default
    #dm_steps(0.0, 1500.0, obs, device="PKMB_DD.ps/CPS")
