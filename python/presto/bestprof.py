from builtins import object
## Automatically adapted for numpy Apr 14, 2006 by convertcode.py

import numpy as num

def get_epochs(line):
    i, f = line.split("=")[-1].split(".")
    f = "0."+f
    epochi = float(i)
    epochf = float(f)
    # Check to see if it is very close to 1 sec
    # If it is, assume the epoch was _exactly_ at the second
    fsec = epochf*86400.0 + 1e-10
    if (num.fabs(fsec - int(fsec)) < 1e-6):
        # print "Looks like an exact second"
        epochf = float(int(fsec))/86400.0
    return epochi, epochf

class bestprof(object):
    def __init__(self, filenm):
        infile = open(filenm)
        self.topo = 0
        self.profile = []
        for line in infile.readlines():
            if line[0]=="#":
                if line.startswith("# Input file"):
                    self.datnm = line.split("=")[-1][:-1]
                    continue
                if line.startswith("# Candidate"):
                    if line.startswith("# Candidate        =  PSR_"):
                        self.psr = line.split("=")[-1].split("_")[1][:-1]
                        continue
                    else:
                        self.psr = None
                if line.startswith("# T_sample"):
                    self.dt = float(line.split("=")[-1])
                    continue
                if line.startswith("# Data Folded"):
                    self.N = float(line.split("=")[-1])
                    continue 
                if line.startswith("# Data Avg"):
                    self.data_avg = float(line.split("=")[-1])
                    continue 
                if line.startswith("# Data StdDev"):
                    self.data_std = float(line.split("=")[-1])
                    continue 
                if line.startswith("# Profile Avg"):
                    self.prof_avg = float(line.split("=")[-1])
                    continue 
                if line.startswith("# Profile StdDev"):
                    self.prof_std = float(line.split("=")[-1])
                    continue 
                if line.startswith("# Reduced chi-sqr"):
                    self.chi_sqr = float(line.split("=")[-1])
                    continue 
                if line.startswith("# Epoch_topo"):
                    try:
                        self.epochi, self.epochf = get_epochs(line)
                        self.epochi_topo, self.epochf_topo = self.epochi, self.epochf
                        self.topo = 1
                    except ValueError:
                        pass
                    continue
                if line.startswith("# Epoch_bary"):
                    try:
                        self.epochi_bary, self.epochf_bary = get_epochs(line)
                        if not self.topo:
                            self.epochi, self.epochf = self.epochi_bary, self.epochf_bary
                    except ValueError:
                        pass
                    continue
                if line.startswith("# P_topo"):
                    try:
                        self.p0_topo = float(line.split("=")[-1].split("+")[0])/1000.0
                        self.p0err_topo = float(line.split("=")[-1].split("+")[1][2:])/1000.0
                        if self.topo:
                            self.p0, self.p0err = self.p0_topo, self.p0err_topo
                    except:
                        pass
                    continue
                if line.startswith("# P_bary"):
                    try:
                        self.p0_bary = float(line.split("=")[-1].split("+")[0])/1000.0
                        self.p0err_bary = float(line.split("=")[-1].split("+")[1][2:])/1000.0
                        if not self.topo:
                            self.p0, self.p0err = self.p0_bary, self.p0err_bary
                    except:
                        pass
                    continue
                if line.startswith("# P'_topo"):
                    try:
                        self.p1_topo = float(line.split("=")[-1].split("+")[0])
                        self.p1err_topo = float(line.split("=")[-1].split("+")[1][2:])
                        if self.topo:
                            self.p1, self.p1err = self.p1_topo, self.p1err_topo
                    except:
                        pass
                    continue
                if line.startswith("# P'_bary"):
                    try:
                        self.p1_bary = float(line.split("=")[-1].split("+")[0])
                        self.p1err_bary = float(line.split("=")[-1].split("+")[1][2:])
                        if not self.topo:
                            self.p1, self.p1err = self.p1_bary, self.p1err_bary
                    except:
                        pass
                    continue
                if line.startswith("# P''_topo"):
                    try:
                        self.p2_topo = float(line.split("=")[-1].split("+")[0])
                        self.p2err_topo = float(line.split("=")[-1].split("+")[1][2:])
                        if self.topo:
                            self.p2, self.p2err = self.p2_topo, self.p2err_topo
                    except:
                        pass
                    continue
                if line.startswith("# P''_bary"):
                    try:
                        self.p2_bary = float(line.split("=")[-1].split("+")[0])
                        self.p2err_bary = float(line.split("=")[-1].split("+")[1][2:])
                        if not self.topo:
                            self.p2, self.p2err = self.p2_bary, self.p2err_bary
                    except:
                        pass
                    continue
            else:
                self.profile.append(float(line.split()[-1]))
        infile.close()
        self.T = self.dt*self.N
        self.proflen = len(self.profile)
    def normalize(self):
        normprof = num.asarray(self.profile)
        normprof -= min(normprof)
        normprof /= max(normprof)
        return normprof

