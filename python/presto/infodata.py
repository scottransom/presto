from builtins import object
## Automatically adapted for numpy Apr 14, 2006 by convertcode.py

class infodata(object):
    def __init__(self, filenm):
        self.breaks = 0
        for line in open(filenm):
            if line.startswith(" Data file name"):
                self.basenm = line.split("=")[-1].strip()
                continue
            if line.startswith(" Telescope"):
                self.telescope = line.split("=")[-1].strip()
                continue
            if line.startswith(" Instrument"):
                self.instrument = line.split("=")[-1].strip()
                continue
            if line.startswith(" Object being observed"):
                self.object = line.split("=")[-1].strip()
                continue
            if line.startswith(" J2000 Right Ascension"):
                self.RA = line.split("=")[-1].strip()
                continue
            if line.startswith(" J2000 Declination"):
                self.DEC = line.split("=")[-1].strip()
                continue
            if line.startswith(" Data observed by"):
                self.observer = line.split("=")[-1].strip()
                continue
            if line.startswith(" Epoch"):
                self.epoch = float(line.split("=")[-1].strip())
                continue
            if line.startswith(" Barycentered?"):
                self.bary = int(line.split("=")[-1].strip())
                continue
            if line.startswith(" Number of bins"):
                self.N = int(line.split("=")[-1].strip())
                continue
            if line.startswith(" Width of each time series bin"):
                self.dt = float(line.split("=")[-1].strip())
                continue
            if line.startswith(" Any breaks in the data?"):
                self.breaks = int(line.split("=")[-1].strip())
                if self.breaks:
                    self.onoff = []
                continue
            if line.startswith(" On/Off bin pair"):
                vals = line.split("=")[-1].strip().split(",")
                self.onoff.append((int(vals[0]), int(vals[1])))
                continue
            if line.startswith(" Type of observation"):
                self.waveband = line.split("=")[-1].strip()
                continue
            if line.startswith(" Beam diameter"):
                self.beam_diam = float(line.split("=")[-1].strip())
                continue
            if line.startswith(" Dispersion measure"):
                self.DM = float(line.split("=")[-1].strip())
                continue
            if line.startswith(" Central freq of low channel"):
                self.lofreq = float(line.split("=")[-1].strip())
                continue
            if line.startswith(" Total bandwidth"):
                self.BW = float(line.split("=")[-1].strip())
                continue
            if line.startswith(" Number of channels"):
                self.numchan = int(line.split("=")[-1].strip())
                continue
            if line.startswith(" Channel bandwidth"):
                self.chan_width = float(line.split("=")[-1].strip())
                continue
            if line.startswith(" Data analyzed by"):
                self.analyzer = line.split("=")[-1].strip()
                continue
            
    def to_file(self, inffn, notes=None):
        if not inffn.endswith(".inf"):
            raise ValueError("PRESTO info files must end with '.inf'. "
                             "Got: %s" % inffn)
        with open(inffn, 'w') as ff:
            if hasattr(self, 'basenm'):
                ff.write(" Data file name without suffix          =  %s\n" %
                         self.basenm)
            if hasattr(self, 'telescope'):
                ff.write(" Telescope used                         =  %s\n" %
                         self.telescope)
            if hasattr(self, 'instrument'):
                ff.write(" Instrument used                        =  %s\n" %
                         self.instrument)
            if hasattr(self, 'object'):
                ff.write(" Object being observed                  =  %s\n" %
                         self.object)
            if hasattr(self, 'RA'):
                ff.write(" J2000 Right Ascension (hh:mm:ss.ssss)  =  %s\n" %
                         self.RA)
            if hasattr(self, 'DEC'):
                ff.write(" J2000 Declination     (dd:mm:ss.ssss)  =  %s\n" %
                         self.DEC)
            if hasattr(self, 'observer'):
                ff.write(" Data observed by                       =  %s\n" %
                         self.observer)
            if hasattr(self, 'epoch'):
                ff.write(" Epoch of observation (MJD)             =  %05.15f\n" %
                         self.epoch)
            if hasattr(self, 'bary'):
                ff.write(" Barycentered?           (1=yes, 0=no)  =  %d\n" %
                         self.bary)
            if hasattr(self, 'N'):
                ff.write(" Number of bins in the time series      =  %-11.0f\n" %
                         self.N)
            if hasattr(self, 'dt'):
                ff.write(" Width of each time series bin (sec)    =  %.15g\n" %
                         self.dt)
            if hasattr(self, 'breaks') and self.breaks:
                ff.write(" Any breaks in the data? (1 yes, 0 no)  =  1\n")
                if hasattr(self, 'onoff'):
                    for ii, (on, off) in enumerate(self.onoff, 1):
                        ff.write(" On/Off bin pair #%3d                   =  %-11.0f, %-11.0f\n" %
                                 (ii, on, off))
            else:
                ff.write(" Any breaks in the data? (1 yes, 0 no)  =  0\n")
            if hasattr(self, 'DM'):
                ff.write(" Dispersion measure (cm-3 pc)           =  %.12g\n" %
                         self.DM)
            if hasattr(self, 'lofreq'):
                ff.write(" Central freq of low channel (Mhz)      =  %.12g\n" %
                         self.lofreq)
            if hasattr(self, 'BW'):
                ff.write(" Total bandwidth (Mhz)                  =  %.12g\n" %
                         self.BW)
            if hasattr(self, 'numchan'):
                ff.write(" Number of channels                     =  %d\n" %
                         self.numchan)
            if hasattr(self, 'chan_width'):
                ff.write(" Channel bandwidth (Mhz)                =  %.12g\n" %
                         self.chan_width)
            if hasattr(self, 'analyzer'):
                ff.write(" Data analyzed by                       =  %s\n" %
                         self.analyzer)
            if hasattr(self, 'deorbited'):
                ff.write(" Orbit removed?          (1=yes, 0=no)  =  %d\n" %
                         self.deorbited)
            ff.write(" Any additional notes:\n")
            if notes is not None:
                ff.write("    %s\n" % notes.strip())
