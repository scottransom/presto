# This file was created automatically by SWIG.
import prestoc
class orbitparams:
    __setmethods__ = {}
    for _s in []: __setmethods__.update(_s.__setmethods__)
    def __setattr__(self,name,value):
        if (name == "this"):
            if isinstance(value,orbitparams):
                self.__dict__[name] = value.this
                if hasattr(value,"thisown"): self.__dict__["thisown"] = value.thisown
                del value.thisown
                return
        method = orbitparams.__setmethods__.get(name,None)
        if method: return method(self,value)
        self.__dict__[name] = value

    __getmethods__ = {}
    for _s in []: __getmethods__.update(_s.__getmethods__)
    def __getattr__(self,name):
        method = orbitparams.__getmethods__.get(name,None)
        if method: return method(self)
        raise AttributeError,name

    __setmethods__["p"] = prestoc.orbitparams_p_set
    __getmethods__["p"] = prestoc.orbitparams_p_get
    __setmethods__["e"] = prestoc.orbitparams_e_set
    __getmethods__["e"] = prestoc.orbitparams_e_get
    __setmethods__["x"] = prestoc.orbitparams_x_set
    __getmethods__["x"] = prestoc.orbitparams_x_get
    __setmethods__["w"] = prestoc.orbitparams_w_set
    __getmethods__["w"] = prestoc.orbitparams_w_get
    __setmethods__["t"] = prestoc.orbitparams_t_set
    __getmethods__["t"] = prestoc.orbitparams_t_get
    __setmethods__["pd"] = prestoc.orbitparams_pd_set
    __getmethods__["pd"] = prestoc.orbitparams_pd_get
    __setmethods__["wd"] = prestoc.orbitparams_wd_set
    __getmethods__["wd"] = prestoc.orbitparams_wd_get
    def __init__(self,*args):
        self.this = apply(prestoc.new_orbitparams,args)
        self.thisown = 1
    def __del__(self,prestoc=prestoc):
        if getattr(self,'thisown',0):
            prestoc.delete_orbitparams(self)
    def __repr__(self):
        return "<C orbitparams instance at %s>" % (self.this,)

class orbitparamsPtr(orbitparams):
    def __init__(self,this):
        self.this = this
        if not hasattr(self,"thisown"): self.thisown = 0
        self.__class__ = orbitparams
prestoc.orbitparams_swigregister(orbitparamsPtr)
class psrdata:
    __setmethods__ = {}
    for _s in []: __setmethods__.update(_s.__setmethods__)
    def __setattr__(self,name,value):
        if (name == "this"):
            if isinstance(value,psrdata):
                self.__dict__[name] = value.this
                if hasattr(value,"thisown"): self.__dict__["thisown"] = value.thisown
                del value.thisown
                return
        method = psrdata.__setmethods__.get(name,None)
        if method: return method(self,value)
        self.__dict__[name] = value

    __getmethods__ = {}
    for _s in []: __getmethods__.update(_s.__getmethods__)
    def __getattr__(self,name):
        method = psrdata.__getmethods__.get(name,None)
        if method: return method(self)
        raise AttributeError,name

    __setmethods__["ra2000"] = prestoc.psrdata_ra2000_set
    __getmethods__["ra2000"] = prestoc.psrdata_ra2000_get
    __setmethods__["ra1950"] = prestoc.psrdata_ra1950_set
    __getmethods__["ra1950"] = prestoc.psrdata_ra1950_get
    __setmethods__["rae"] = prestoc.psrdata_rae_set
    __getmethods__["rae"] = prestoc.psrdata_rae_get
    __setmethods__["dec2000"] = prestoc.psrdata_dec2000_set
    __getmethods__["dec2000"] = prestoc.psrdata_dec2000_get
    __setmethods__["dec1950"] = prestoc.psrdata_dec1950_set
    __getmethods__["dec1950"] = prestoc.psrdata_dec1950_get
    __setmethods__["dece"] = prestoc.psrdata_dece_set
    __getmethods__["dece"] = prestoc.psrdata_dece_get
    __setmethods__["dmin__"] = prestoc.psrdata_dmin___set
    __getmethods__["dmin__"] = prestoc.psrdata_dmin___get
    __setmethods__["dmax__"] = prestoc.psrdata_dmax___set
    __getmethods__["dmax__"] = prestoc.psrdata_dmax___get
    __setmethods__["dist"] = prestoc.psrdata_dist_set
    __getmethods__["dist"] = prestoc.psrdata_dist_get
    __setmethods__["ldeg"] = prestoc.psrdata_ldeg_set
    __getmethods__["ldeg"] = prestoc.psrdata_ldeg_get
    __setmethods__["bdeg"] = prestoc.psrdata_bdeg_set
    __getmethods__["bdeg"] = prestoc.psrdata_bdeg_get
    __setmethods__["pmra"] = prestoc.psrdata_pmra_set
    __getmethods__["pmra"] = prestoc.psrdata_pmra_get
    __setmethods__["pmrae"] = prestoc.psrdata_pmrae_set
    __getmethods__["pmrae"] = prestoc.psrdata_pmrae_get
    __setmethods__["pmdec"] = prestoc.psrdata_pmdec_set
    __getmethods__["pmdec"] = prestoc.psrdata_pmdec_get
    __setmethods__["pmdece"] = prestoc.psrdata_pmdece_set
    __getmethods__["pmdece"] = prestoc.psrdata_pmdece_get
    __setmethods__["posepoch"] = prestoc.psrdata_posepoch_set
    __getmethods__["posepoch"] = prestoc.psrdata_posepoch_get
    __setmethods__["p"] = prestoc.psrdata_p_set
    __getmethods__["p"] = prestoc.psrdata_p_get
    __setmethods__["pe"] = prestoc.psrdata_pe_set
    __getmethods__["pe"] = prestoc.psrdata_pe_get
    __setmethods__["pdot"] = prestoc.psrdata_pdot_set
    __getmethods__["pdot"] = prestoc.psrdata_pdot_get
    __setmethods__["pdote"] = prestoc.psrdata_pdote_set
    __getmethods__["pdote"] = prestoc.psrdata_pdote_get
    __setmethods__["f2"] = prestoc.psrdata_f2_set
    __getmethods__["f2"] = prestoc.psrdata_f2_get
    __setmethods__["f2e"] = prestoc.psrdata_f2e_set
    __getmethods__["f2e"] = prestoc.psrdata_f2e_get
    __setmethods__["f3"] = prestoc.psrdata_f3_set
    __getmethods__["f3"] = prestoc.psrdata_f3_get
    __setmethods__["f3e"] = prestoc.psrdata_f3e_set
    __getmethods__["f3e"] = prestoc.psrdata_f3e_get
    __setmethods__["epoch"] = prestoc.psrdata_epoch_set
    __getmethods__["epoch"] = prestoc.psrdata_epoch_get
    __setmethods__["dm"] = prestoc.psrdata_dm_set
    __getmethods__["dm"] = prestoc.psrdata_dm_get
    __setmethods__["dme"] = prestoc.psrdata_dme_set
    __getmethods__["dme"] = prestoc.psrdata_dme_get
    __setmethods__["rm"] = prestoc.psrdata_rm_set
    __getmethods__["rm"] = prestoc.psrdata_rm_get
    __setmethods__["rme"] = prestoc.psrdata_rme_set
    __getmethods__["rme"] = prestoc.psrdata_rme_get
    __setmethods__["we"] = prestoc.psrdata_we_set
    __getmethods__["we"] = prestoc.psrdata_we_get
    __setmethods__["w50"] = prestoc.psrdata_w50_set
    __getmethods__["w50"] = prestoc.psrdata_w50_get
    __setmethods__["w10"] = prestoc.psrdata_w10_set
    __getmethods__["w10"] = prestoc.psrdata_w10_get
    __setmethods__["s400"] = prestoc.psrdata_s400_set
    __getmethods__["s400"] = prestoc.psrdata_s400_get
    __setmethods__["s600"] = prestoc.psrdata_s600_set
    __getmethods__["s600"] = prestoc.psrdata_s600_get
    __setmethods__["s1400"] = prestoc.psrdata_s1400_set
    __getmethods__["s1400"] = prestoc.psrdata_s1400_get
    __setmethods__["tau"] = prestoc.psrdata_tau_set
    __getmethods__["tau"] = prestoc.psrdata_tau_get
    __setmethods__["t408"] = prestoc.psrdata_t408_set
    __getmethods__["t408"] = prestoc.psrdata_t408_get
    __setmethods__["distmod"] = prestoc.psrdata_distmod_set
    __getmethods__["distmod"] = prestoc.psrdata_distmod_get
    __setmethods__["lum"] = prestoc.psrdata_lum_set
    __getmethods__["lum"] = prestoc.psrdata_lum_get
    __setmethods__["bsurf"] = prestoc.psrdata_bsurf_set
    __getmethods__["bsurf"] = prestoc.psrdata_bsurf_get
    __setmethods__["age"] = prestoc.psrdata_age_set
    __getmethods__["age"] = prestoc.psrdata_age_get
    __setmethods__["edot"] = prestoc.psrdata_edot_set
    __getmethods__["edot"] = prestoc.psrdata_edot_get
    __setmethods__["pb"] = prestoc.psrdata_pb_set
    __getmethods__["pb"] = prestoc.psrdata_pb_get
    __setmethods__["pbe"] = prestoc.psrdata_pbe_set
    __getmethods__["pbe"] = prestoc.psrdata_pbe_get
    __setmethods__["a1"] = prestoc.psrdata_a1_set
    __getmethods__["a1"] = prestoc.psrdata_a1_get
    __setmethods__["a1e"] = prestoc.psrdata_a1e_set
    __getmethods__["a1e"] = prestoc.psrdata_a1e_get
    __setmethods__["om"] = prestoc.psrdata_om_set
    __getmethods__["om"] = prestoc.psrdata_om_get
    __setmethods__["ome"] = prestoc.psrdata_ome_set
    __getmethods__["ome"] = prestoc.psrdata_ome_get
    __setmethods__["omdot"] = prestoc.psrdata_omdot_set
    __getmethods__["omdot"] = prestoc.psrdata_omdot_get
    __setmethods__["omdote"] = prestoc.psrdata_omdote_set
    __getmethods__["omdote"] = prestoc.psrdata_omdote_get
    __setmethods__["e"] = prestoc.psrdata_e_set
    __getmethods__["e"] = prestoc.psrdata_e_get
    __setmethods__["ee"] = prestoc.psrdata_ee_set
    __getmethods__["ee"] = prestoc.psrdata_ee_get
    __setmethods__["t0"] = prestoc.psrdata_t0_set
    __getmethods__["t0"] = prestoc.psrdata_t0_get
    __setmethods__["t0e"] = prestoc.psrdata_t0e_set
    __getmethods__["t0e"] = prestoc.psrdata_t0e_get
    __setmethods__["gamma"] = prestoc.psrdata_gamma_set
    __getmethods__["gamma"] = prestoc.psrdata_gamma_get
    __setmethods__["gammae"] = prestoc.psrdata_gammae_set
    __getmethods__["gammae"] = prestoc.psrdata_gammae_get
    __setmethods__["pbdot"] = prestoc.psrdata_pbdot_set
    __getmethods__["pbdot"] = prestoc.psrdata_pbdot_get
    __setmethods__["pbdote"] = prestoc.psrdata_pbdote_set
    __getmethods__["pbdote"] = prestoc.psrdata_pbdote_get
    __setmethods__["si"] = prestoc.psrdata_si_set
    __getmethods__["si"] = prestoc.psrdata_si_get
    __setmethods__["sie"] = prestoc.psrdata_sie_set
    __getmethods__["sie"] = prestoc.psrdata_sie_get
    __setmethods__["r__"] = prestoc.psrdata_r___set
    __getmethods__["r__"] = prestoc.psrdata_r___get
    __setmethods__["re"] = prestoc.psrdata_re_set
    __getmethods__["re"] = prestoc.psrdata_re_get
    __setmethods__["pb2"] = prestoc.psrdata_pb2_set
    __getmethods__["pb2"] = prestoc.psrdata_pb2_get
    __setmethods__["pb2e"] = prestoc.psrdata_pb2e_set
    __getmethods__["pb2e"] = prestoc.psrdata_pb2e_get
    __setmethods__["a12"] = prestoc.psrdata_a12_set
    __getmethods__["a12"] = prestoc.psrdata_a12_get
    __setmethods__["a12e"] = prestoc.psrdata_a12e_set
    __getmethods__["a12e"] = prestoc.psrdata_a12e_get
    __setmethods__["om2"] = prestoc.psrdata_om2_set
    __getmethods__["om2"] = prestoc.psrdata_om2_get
    __setmethods__["om2e"] = prestoc.psrdata_om2e_set
    __getmethods__["om2e"] = prestoc.psrdata_om2e_get
    __setmethods__["omdot2"] = prestoc.psrdata_omdot2_set
    __getmethods__["omdot2"] = prestoc.psrdata_omdot2_get
    __setmethods__["omdot2e"] = prestoc.psrdata_omdot2e_set
    __getmethods__["omdot2e"] = prestoc.psrdata_omdot2e_get
    __setmethods__["e2"] = prestoc.psrdata_e2_set
    __getmethods__["e2"] = prestoc.psrdata_e2_get
    __setmethods__["e2e"] = prestoc.psrdata_e2e_set
    __getmethods__["e2e"] = prestoc.psrdata_e2e_get
    __setmethods__["t02"] = prestoc.psrdata_t02_set
    __getmethods__["t02"] = prestoc.psrdata_t02_get
    __setmethods__["t02e"] = prestoc.psrdata_t02e_set
    __getmethods__["t02e"] = prestoc.psrdata_t02e_get
    __setmethods__["gamma2"] = prestoc.psrdata_gamma2_set
    __getmethods__["gamma2"] = prestoc.psrdata_gamma2_get
    __setmethods__["gamma2e"] = prestoc.psrdata_gamma2e_set
    __getmethods__["gamma2e"] = prestoc.psrdata_gamma2e_get
    __setmethods__["pbdot2"] = prestoc.psrdata_pbdot2_set
    __getmethods__["pbdot2"] = prestoc.psrdata_pbdot2_get
    __setmethods__["pbdot2e"] = prestoc.psrdata_pbdot2e_set
    __getmethods__["pbdot2e"] = prestoc.psrdata_pbdot2e_get
    __setmethods__["si2"] = prestoc.psrdata_si2_set
    __getmethods__["si2"] = prestoc.psrdata_si2_get
    __setmethods__["si2e"] = prestoc.psrdata_si2e_set
    __getmethods__["si2e"] = prestoc.psrdata_si2e_get
    __setmethods__["r2"] = prestoc.psrdata_r2_set
    __getmethods__["r2"] = prestoc.psrdata_r2_get
    __setmethods__["r2e"] = prestoc.psrdata_r2e_set
    __getmethods__["r2e"] = prestoc.psrdata_r2e_get
    __setmethods__["nscode"] = prestoc.psrdata_nscode_set
    __getmethods__["nscode"] = prestoc.psrdata_nscode_get
    __setmethods__["ndflag"] = prestoc.psrdata_ndflag_set
    __getmethods__["ndflag"] = prestoc.psrdata_ndflag_get
    __setmethods__["ntauflag"] = prestoc.psrdata_ntauflag_set
    __getmethods__["ntauflag"] = prestoc.psrdata_ntauflag_get
    __setmethods__["ntype"] = prestoc.psrdata_ntype_set
    __getmethods__["ntype"] = prestoc.psrdata_ntype_get
    __setmethods__["modcode"] = prestoc.psrdata_modcode_set
    __getmethods__["modcode"] = prestoc.psrdata_modcode_get
    __setmethods__["limcode"] = prestoc.psrdata_limcode_set
    __getmethods__["limcode"] = prestoc.psrdata_limcode_get
    __setmethods__["ibin"] = prestoc.psrdata_ibin_set
    __getmethods__["ibin"] = prestoc.psrdata_ibin_get
    __setmethods__["jname"] = prestoc.psrdata_jname_set
    __getmethods__["jname"] = prestoc.psrdata_jname_get
    __setmethods__["bname"] = prestoc.psrdata_bname_set
    __getmethods__["bname"] = prestoc.psrdata_bname_get
    __setmethods__["lcode"] = prestoc.psrdata_lcode_set
    __getmethods__["lcode"] = prestoc.psrdata_lcode_get
    __setmethods__["ucode"] = prestoc.psrdata_ucode_set
    __getmethods__["ucode"] = prestoc.psrdata_ucode_get
    def __init__(self,*args):
        self.this = apply(prestoc.new_psrdata,args)
        self.thisown = 1
    def __del__(self,prestoc=prestoc):
        if getattr(self,'thisown',0):
            prestoc.delete_psrdata(self)
    def __repr__(self):
        return "<C psrdata instance at %s>" % (self.this,)

class psrdataPtr(psrdata):
    def __init__(self,this):
        self.this = this
        if not hasattr(self,"thisown"): self.thisown = 0
        self.__class__ = psrdata
prestoc.psrdata_swigregister(psrdataPtr)
class psrparams:
    __setmethods__ = {}
    for _s in []: __setmethods__.update(_s.__setmethods__)
    def __setattr__(self,name,value):
        if (name == "this"):
            if isinstance(value,psrparams):
                self.__dict__[name] = value.this
                if hasattr(value,"thisown"): self.__dict__["thisown"] = value.thisown
                del value.thisown
                return
        method = psrparams.__setmethods__.get(name,None)
        if method: return method(self,value)
        self.__dict__[name] = value

    __getmethods__ = {}
    for _s in []: __getmethods__.update(_s.__getmethods__)
    def __getattr__(self,name):
        method = psrparams.__getmethods__.get(name,None)
        if method: return method(self)
        raise AttributeError,name

    __setmethods__["jname"] = prestoc.psrparams_jname_set
    __getmethods__["jname"] = prestoc.psrparams_jname_get
    __setmethods__["bname"] = prestoc.psrparams_bname_set
    __getmethods__["bname"] = prestoc.psrparams_bname_get
    __setmethods__["ntype"] = prestoc.psrparams_ntype_set
    __getmethods__["ntype"] = prestoc.psrparams_ntype_get
    __setmethods__["ra2000"] = prestoc.psrparams_ra2000_set
    __getmethods__["ra2000"] = prestoc.psrparams_ra2000_get
    __setmethods__["dec2000"] = prestoc.psrparams_dec2000_set
    __getmethods__["dec2000"] = prestoc.psrparams_dec2000_get
    __setmethods__["dm"] = prestoc.psrparams_dm_set
    __getmethods__["dm"] = prestoc.psrparams_dm_get
    __setmethods__["dist"] = prestoc.psrparams_dist_set
    __getmethods__["dist"] = prestoc.psrparams_dist_get
    __setmethods__["fwhm"] = prestoc.psrparams_fwhm_set
    __getmethods__["fwhm"] = prestoc.psrparams_fwhm_get
    __setmethods__["timepoch"] = prestoc.psrparams_timepoch_set
    __getmethods__["timepoch"] = prestoc.psrparams_timepoch_get
    __setmethods__["p"] = prestoc.psrparams_p_set
    __getmethods__["p"] = prestoc.psrparams_p_get
    __setmethods__["pd"] = prestoc.psrparams_pd_set
    __getmethods__["pd"] = prestoc.psrparams_pd_get
    __setmethods__["pdd"] = prestoc.psrparams_pdd_set
    __getmethods__["pdd"] = prestoc.psrparams_pdd_get
    __setmethods__["f"] = prestoc.psrparams_f_set
    __getmethods__["f"] = prestoc.psrparams_f_get
    __setmethods__["fd"] = prestoc.psrparams_fd_set
    __getmethods__["fd"] = prestoc.psrparams_fd_get
    __setmethods__["fdd"] = prestoc.psrparams_fdd_set
    __getmethods__["fdd"] = prestoc.psrparams_fdd_get
    __setmethods__["orb"] = prestoc.psrparams_orb_set
    __getmethods__["orb"] = prestoc.psrparams_orb_get
    def __init__(self,*args):
        self.this = apply(prestoc.new_psrparams,args)
        self.thisown = 1
    def __del__(self,prestoc=prestoc):
        if getattr(self,'thisown',0):
            prestoc.delete_psrparams(self)
    def __repr__(self):
        return "<C psrparams instance at %s>" % (self.this,)

class psrparamsPtr(psrparams):
    def __init__(self,this):
        self.this = this
        if not hasattr(self,"thisown"): self.thisown = 0
        self.__class__ = psrparams
prestoc.psrparams_swigregister(psrparamsPtr)
class DoubleArray:
    __setmethods__ = {}
    for _s in []: __setmethods__.update(_s.__setmethods__)
    def __setattr__(self,name,value):
        if (name == "this"):
            if isinstance(value,DoubleArray):
                self.__dict__[name] = value.this
                if hasattr(value,"thisown"): self.__dict__["thisown"] = value.thisown
                del value.thisown
                return
        method = DoubleArray.__setmethods__.get(name,None)
        if method: return method(self,value)
        self.__dict__[name] = value

    __getmethods__ = {}
    for _s in []: __getmethods__.update(_s.__getmethods__)
    def __getattr__(self,name):
        method = DoubleArray.__getmethods__.get(name,None)
        if method: return method(self)
        raise AttributeError,name

    __setmethods__["dptr"] = prestoc.DoubleArray_dptr_set
    __getmethods__["dptr"] = prestoc.DoubleArray_dptr_get
    def __init__(self,*args):
        self.this = apply(prestoc.new_DoubleArray,args)
        self.thisown = 1
    def __del__(self,prestoc=prestoc):
        if getattr(self,'thisown',0):
            prestoc.delete_DoubleArray(self)
    def __getitem__(*args): return apply(prestoc.DoubleArray___getitem__,args)
    def __setitem__(*args): return apply(prestoc.DoubleArray___setitem__,args)
    def __repr__(self):
        return "<C DoubleArray instance at %s>" % (self.this,)

class DoubleArrayPtr(DoubleArray):
    def __init__(self,this):
        self.this = this
        if not hasattr(self,"thisown"): self.thisown = 0
        self.__class__ = DoubleArray
prestoc.DoubleArray_swigregister(DoubleArrayPtr)
class infodata:
    __setmethods__ = {}
    for _s in []: __setmethods__.update(_s.__setmethods__)
    def __setattr__(self,name,value):
        if (name == "this"):
            if isinstance(value,infodata):
                self.__dict__[name] = value.this
                if hasattr(value,"thisown"): self.__dict__["thisown"] = value.thisown
                del value.thisown
                return
        method = infodata.__setmethods__.get(name,None)
        if method: return method(self,value)
        self.__dict__[name] = value

    __getmethods__ = {}
    for _s in []: __getmethods__.update(_s.__getmethods__)
    def __getattr__(self,name):
        method = infodata.__getmethods__.get(name,None)
        if method: return method(self)
        raise AttributeError,name

    __setmethods__["ra_s"] = prestoc.infodata_ra_s_set
    __getmethods__["ra_s"] = prestoc.infodata_ra_s_get
    __setmethods__["dec_s"] = prestoc.infodata_dec_s_set
    __getmethods__["dec_s"] = prestoc.infodata_dec_s_get
    __setmethods__["N"] = prestoc.infodata_N_set
    __getmethods__["N"] = prestoc.infodata_N_get
    __setmethods__["dt"] = prestoc.infodata_dt_set
    __getmethods__["dt"] = prestoc.infodata_dt_get
    __setmethods__["fov"] = prestoc.infodata_fov_set
    __getmethods__["fov"] = prestoc.infodata_fov_get
    __setmethods__["mjd_f"] = prestoc.infodata_mjd_f_set
    __getmethods__["mjd_f"] = prestoc.infodata_mjd_f_get
    __setmethods__["dm"] = prestoc.infodata_dm_set
    __getmethods__["dm"] = prestoc.infodata_dm_get
    __setmethods__["freq"] = prestoc.infodata_freq_set
    __getmethods__["freq"] = prestoc.infodata_freq_get
    __setmethods__["freqband"] = prestoc.infodata_freqband_set
    __getmethods__["freqband"] = prestoc.infodata_freqband_get
    __setmethods__["chan_wid"] = prestoc.infodata_chan_wid_set
    __getmethods__["chan_wid"] = prestoc.infodata_chan_wid_get
    __setmethods__["wavelen"] = prestoc.infodata_wavelen_set
    __getmethods__["wavelen"] = prestoc.infodata_wavelen_get
    __setmethods__["waveband"] = prestoc.infodata_waveband_set
    __getmethods__["waveband"] = prestoc.infodata_waveband_get
    __setmethods__["energy"] = prestoc.infodata_energy_set
    __getmethods__["energy"] = prestoc.infodata_energy_get
    __setmethods__["energyband"] = prestoc.infodata_energyband_set
    __getmethods__["energyband"] = prestoc.infodata_energyband_get
    __getmethods__["onoff"] = prestoc.infodata_onoff_get
    __setmethods__["num_chan"] = prestoc.infodata_num_chan_set
    __getmethods__["num_chan"] = prestoc.infodata_num_chan_get
    __setmethods__["mjd_i"] = prestoc.infodata_mjd_i_set
    __getmethods__["mjd_i"] = prestoc.infodata_mjd_i_get
    __setmethods__["ra_h"] = prestoc.infodata_ra_h_set
    __getmethods__["ra_h"] = prestoc.infodata_ra_h_get
    __setmethods__["ra_m"] = prestoc.infodata_ra_m_set
    __getmethods__["ra_m"] = prestoc.infodata_ra_m_get
    __setmethods__["dec_d"] = prestoc.infodata_dec_d_set
    __getmethods__["dec_d"] = prestoc.infodata_dec_d_get
    __setmethods__["dec_m"] = prestoc.infodata_dec_m_set
    __getmethods__["dec_m"] = prestoc.infodata_dec_m_get
    __setmethods__["bary"] = prestoc.infodata_bary_set
    __getmethods__["bary"] = prestoc.infodata_bary_get
    __setmethods__["numonoff"] = prestoc.infodata_numonoff_set
    __getmethods__["numonoff"] = prestoc.infodata_numonoff_get
    __setmethods__["notes"] = prestoc.infodata_notes_set
    __getmethods__["notes"] = prestoc.infodata_notes_get
    __setmethods__["name"] = prestoc.infodata_name_set
    __getmethods__["name"] = prestoc.infodata_name_get
    __setmethods__["object"] = prestoc.infodata_object_set
    __getmethods__["object"] = prestoc.infodata_object_get
    __setmethods__["instrument"] = prestoc.infodata_instrument_set
    __getmethods__["instrument"] = prestoc.infodata_instrument_get
    __setmethods__["observer"] = prestoc.infodata_observer_set
    __getmethods__["observer"] = prestoc.infodata_observer_get
    __setmethods__["analyzer"] = prestoc.infodata_analyzer_set
    __getmethods__["analyzer"] = prestoc.infodata_analyzer_get
    __setmethods__["telescope"] = prestoc.infodata_telescope_set
    __getmethods__["telescope"] = prestoc.infodata_telescope_get
    __setmethods__["band"] = prestoc.infodata_band_set
    __getmethods__["band"] = prestoc.infodata_band_get
    __setmethods__["filt"] = prestoc.infodata_filt_set
    __getmethods__["filt"] = prestoc.infodata_filt_get
    def __init__(self,*args):
        self.this = apply(prestoc.new_infodata,args)
        self.thisown = 1
    def __del__(self,prestoc=prestoc):
        if getattr(self,'thisown',0):
            prestoc.delete_infodata(self)
    def __repr__(self):
        return "<C infodata instance at %s>" % (self.this,)

class infodataPtr(infodata):
    def __init__(self,this):
        self.this = this
        if not hasattr(self,"thisown"): self.thisown = 0
        self.__class__ = infodata
prestoc.infodata_swigregister(infodataPtr)
class makedata:
    __setmethods__ = {}
    for _s in []: __setmethods__.update(_s.__setmethods__)
    def __setattr__(self,name,value):
        if (name == "this"):
            if isinstance(value,makedata):
                self.__dict__[name] = value.this
                if hasattr(value,"thisown"): self.__dict__["thisown"] = value.thisown
                del value.thisown
                return
        method = makedata.__setmethods__.get(name,None)
        if method: return method(self,value)
        self.__dict__[name] = value

    __getmethods__ = {}
    for _s in []: __getmethods__.update(_s.__getmethods__)
    def __getattr__(self,name):
        method = makedata.__getmethods__.get(name,None)
        if method: return method(self)
        raise AttributeError,name

    __setmethods__["basefilenm"] = prestoc.makedata_basefilenm_set
    __getmethods__["basefilenm"] = prestoc.makedata_basefilenm_get
    __setmethods__["description"] = prestoc.makedata_description_set
    __getmethods__["description"] = prestoc.makedata_description_get
    __setmethods__["N"] = prestoc.makedata_N_set
    __getmethods__["N"] = prestoc.makedata_N_get
    __setmethods__["next2_to_n"] = prestoc.makedata_next2_to_n_set
    __getmethods__["next2_to_n"] = prestoc.makedata_next2_to_n_get
    __setmethods__["dt"] = prestoc.makedata_dt_set
    __getmethods__["dt"] = prestoc.makedata_dt_get
    __setmethods__["T"] = prestoc.makedata_T_set
    __getmethods__["T"] = prestoc.makedata_T_get
    __setmethods__["ptype"] = prestoc.makedata_ptype_set
    __getmethods__["ptype"] = prestoc.makedata_ptype_get
    __setmethods__["pnum"] = prestoc.makedata_pnum_set
    __getmethods__["pnum"] = prestoc.makedata_pnum_get
    __setmethods__["fwhm"] = prestoc.makedata_fwhm_set
    __getmethods__["fwhm"] = prestoc.makedata_fwhm_get
    __setmethods__["round"] = prestoc.makedata_round_set
    __getmethods__["round"] = prestoc.makedata_round_get
    __setmethods__["roundnum"] = prestoc.makedata_roundnum_set
    __getmethods__["roundnum"] = prestoc.makedata_roundnum_get
    __setmethods__["f"] = prestoc.makedata_f_set
    __getmethods__["f"] = prestoc.makedata_f_get
    __setmethods__["fd"] = prestoc.makedata_fd_set
    __getmethods__["fd"] = prestoc.makedata_fd_get
    __setmethods__["fdd"] = prestoc.makedata_fdd_set
    __getmethods__["fdd"] = prestoc.makedata_fdd_get
    __setmethods__["p"] = prestoc.makedata_p_set
    __getmethods__["p"] = prestoc.makedata_p_get
    __setmethods__["pd"] = prestoc.makedata_pd_set
    __getmethods__["pd"] = prestoc.makedata_pd_get
    __setmethods__["pdd"] = prestoc.makedata_pdd_set
    __getmethods__["pdd"] = prestoc.makedata_pdd_get
    __setmethods__["r"] = prestoc.makedata_r_set
    __getmethods__["r"] = prestoc.makedata_r_get
    __setmethods__["z"] = prestoc.makedata_z_set
    __getmethods__["z"] = prestoc.makedata_z_get
    __setmethods__["w"] = prestoc.makedata_w_set
    __getmethods__["w"] = prestoc.makedata_w_get
    __setmethods__["amp"] = prestoc.makedata_amp_set
    __getmethods__["amp"] = prestoc.makedata_amp_get
    __setmethods__["phs"] = prestoc.makedata_phs_set
    __getmethods__["phs"] = prestoc.makedata_phs_get
    __setmethods__["dc"] = prestoc.makedata_dc_set
    __getmethods__["dc"] = prestoc.makedata_dc_get
    __setmethods__["binary"] = prestoc.makedata_binary_set
    __getmethods__["binary"] = prestoc.makedata_binary_get
    __setmethods__["orb"] = prestoc.makedata_orb_set
    __getmethods__["orb"] = prestoc.makedata_orb_get
    __setmethods__["ampmod"] = prestoc.makedata_ampmod_set
    __getmethods__["ampmod"] = prestoc.makedata_ampmod_get
    __setmethods__["ampmoda"] = prestoc.makedata_ampmoda_set
    __getmethods__["ampmoda"] = prestoc.makedata_ampmoda_get
    __setmethods__["ampmodf"] = prestoc.makedata_ampmodf_set
    __getmethods__["ampmodf"] = prestoc.makedata_ampmodf_get
    __setmethods__["ampmodp"] = prestoc.makedata_ampmodp_set
    __getmethods__["ampmodp"] = prestoc.makedata_ampmodp_get
    __setmethods__["noisetype"] = prestoc.makedata_noisetype_set
    __getmethods__["noisetype"] = prestoc.makedata_noisetype_get
    __setmethods__["noise"] = prestoc.makedata_noise_set
    __getmethods__["noise"] = prestoc.makedata_noise_get
    __setmethods__["noisesig"] = prestoc.makedata_noisesig_set
    __getmethods__["noisesig"] = prestoc.makedata_noisesig_get
    __setmethods__["numonoff"] = prestoc.makedata_numonoff_set
    __getmethods__["numonoff"] = prestoc.makedata_numonoff_get
    __setmethods__["onoff"] = prestoc.makedata_onoff_set
    __getmethods__["onoff"] = prestoc.makedata_onoff_get
    def __init__(self,*args):
        self.this = apply(prestoc.new_makedata,args)
        self.thisown = 1
    def __del__(self,prestoc=prestoc):
        if getattr(self,'thisown',0):
            prestoc.delete_makedata(self)
    def __repr__(self):
        return "<C makedata instance at %s>" % (self.this,)

class makedataPtr(makedata):
    def __init__(self,this):
        self.this = this
        if not hasattr(self,"thisown"): self.thisown = 0
        self.__class__ = makedata
prestoc.makedata_swigregister(makedataPtr)
class rderivs:
    __setmethods__ = {}
    for _s in []: __setmethods__.update(_s.__setmethods__)
    def __setattr__(self,name,value):
        if (name == "this"):
            if isinstance(value,rderivs):
                self.__dict__[name] = value.this
                if hasattr(value,"thisown"): self.__dict__["thisown"] = value.thisown
                del value.thisown
                return
        method = rderivs.__setmethods__.get(name,None)
        if method: return method(self,value)
        self.__dict__[name] = value

    __getmethods__ = {}
    for _s in []: __getmethods__.update(_s.__getmethods__)
    def __getattr__(self,name):
        method = rderivs.__getmethods__.get(name,None)
        if method: return method(self)
        raise AttributeError,name

    __setmethods__["pow"] = prestoc.rderivs_pow_set
    __getmethods__["pow"] = prestoc.rderivs_pow_get
    __setmethods__["phs"] = prestoc.rderivs_phs_set
    __getmethods__["phs"] = prestoc.rderivs_phs_get
    __setmethods__["dpow"] = prestoc.rderivs_dpow_set
    __getmethods__["dpow"] = prestoc.rderivs_dpow_get
    __setmethods__["dphs"] = prestoc.rderivs_dphs_set
    __getmethods__["dphs"] = prestoc.rderivs_dphs_get
    __setmethods__["d2pow"] = prestoc.rderivs_d2pow_set
    __getmethods__["d2pow"] = prestoc.rderivs_d2pow_get
    __setmethods__["d2phs"] = prestoc.rderivs_d2phs_set
    __getmethods__["d2phs"] = prestoc.rderivs_d2phs_get
    __setmethods__["locpow"] = prestoc.rderivs_locpow_set
    __getmethods__["locpow"] = prestoc.rderivs_locpow_get
    def __init__(self,*args):
        self.this = apply(prestoc.new_rderivs,args)
        self.thisown = 1
    def __del__(self,prestoc=prestoc):
        if getattr(self,'thisown',0):
            prestoc.delete_rderivs(self)
    def __repr__(self):
        return "<C rderivs instance at %s>" % (self.this,)

class rderivsPtr(rderivs):
    def __init__(self,this):
        self.this = this
        if not hasattr(self,"thisown"): self.thisown = 0
        self.__class__ = rderivs
prestoc.rderivs_swigregister(rderivsPtr)
class fourierprops:
    __setmethods__ = {}
    for _s in []: __setmethods__.update(_s.__setmethods__)
    def __setattr__(self,name,value):
        if (name == "this"):
            if isinstance(value,fourierprops):
                self.__dict__[name] = value.this
                if hasattr(value,"thisown"): self.__dict__["thisown"] = value.thisown
                del value.thisown
                return
        method = fourierprops.__setmethods__.get(name,None)
        if method: return method(self,value)
        self.__dict__[name] = value

    __getmethods__ = {}
    for _s in []: __getmethods__.update(_s.__getmethods__)
    def __getattr__(self,name):
        method = fourierprops.__getmethods__.get(name,None)
        if method: return method(self)
        raise AttributeError,name

    __setmethods__["r"] = prestoc.fourierprops_r_set
    __getmethods__["r"] = prestoc.fourierprops_r_get
    __setmethods__["rerr"] = prestoc.fourierprops_rerr_set
    __getmethods__["rerr"] = prestoc.fourierprops_rerr_get
    __setmethods__["z"] = prestoc.fourierprops_z_set
    __getmethods__["z"] = prestoc.fourierprops_z_get
    __setmethods__["zerr"] = prestoc.fourierprops_zerr_set
    __getmethods__["zerr"] = prestoc.fourierprops_zerr_get
    __setmethods__["w"] = prestoc.fourierprops_w_set
    __getmethods__["w"] = prestoc.fourierprops_w_get
    __setmethods__["werr"] = prestoc.fourierprops_werr_set
    __getmethods__["werr"] = prestoc.fourierprops_werr_get
    __setmethods__["pow"] = prestoc.fourierprops_pow_set
    __getmethods__["pow"] = prestoc.fourierprops_pow_get
    __setmethods__["powerr"] = prestoc.fourierprops_powerr_set
    __getmethods__["powerr"] = prestoc.fourierprops_powerr_get
    __setmethods__["sig"] = prestoc.fourierprops_sig_set
    __getmethods__["sig"] = prestoc.fourierprops_sig_get
    __setmethods__["rawpow"] = prestoc.fourierprops_rawpow_set
    __getmethods__["rawpow"] = prestoc.fourierprops_rawpow_get
    __setmethods__["phs"] = prestoc.fourierprops_phs_set
    __getmethods__["phs"] = prestoc.fourierprops_phs_get
    __setmethods__["phserr"] = prestoc.fourierprops_phserr_set
    __getmethods__["phserr"] = prestoc.fourierprops_phserr_get
    __setmethods__["cen"] = prestoc.fourierprops_cen_set
    __getmethods__["cen"] = prestoc.fourierprops_cen_get
    __setmethods__["cenerr"] = prestoc.fourierprops_cenerr_set
    __getmethods__["cenerr"] = prestoc.fourierprops_cenerr_get
    __setmethods__["pur"] = prestoc.fourierprops_pur_set
    __getmethods__["pur"] = prestoc.fourierprops_pur_get
    __setmethods__["purerr"] = prestoc.fourierprops_purerr_set
    __getmethods__["purerr"] = prestoc.fourierprops_purerr_get
    __setmethods__["locpow"] = prestoc.fourierprops_locpow_set
    __getmethods__["locpow"] = prestoc.fourierprops_locpow_get
    def __init__(self,*args):
        self.this = apply(prestoc.new_fourierprops,args)
        self.thisown = 1
    def __del__(self,prestoc=prestoc):
        if getattr(self,'thisown',0):
            prestoc.delete_fourierprops(self)
    def __repr__(self):
        return "<C fourierprops instance at %s>" % (self.this,)

class fourierpropsPtr(fourierprops):
    def __init__(self,this):
        self.this = this
        if not hasattr(self,"thisown"): self.thisown = 0
        self.__class__ = fourierprops
prestoc.fourierprops_swigregister(fourierpropsPtr)
class binaryprops:
    __setmethods__ = {}
    for _s in []: __setmethods__.update(_s.__setmethods__)
    def __setattr__(self,name,value):
        if (name == "this"):
            if isinstance(value,binaryprops):
                self.__dict__[name] = value.this
                if hasattr(value,"thisown"): self.__dict__["thisown"] = value.thisown
                del value.thisown
                return
        method = binaryprops.__setmethods__.get(name,None)
        if method: return method(self,value)
        self.__dict__[name] = value

    __getmethods__ = {}
    for _s in []: __getmethods__.update(_s.__getmethods__)
    def __getattr__(self,name):
        method = binaryprops.__getmethods__.get(name,None)
        if method: return method(self)
        raise AttributeError,name

    __setmethods__["ppsr"] = prestoc.binaryprops_ppsr_set
    __getmethods__["ppsr"] = prestoc.binaryprops_ppsr_get
    __setmethods__["fpsr"] = prestoc.binaryprops_fpsr_set
    __getmethods__["fpsr"] = prestoc.binaryprops_fpsr_get
    __setmethods__["rpsr"] = prestoc.binaryprops_rpsr_set
    __getmethods__["rpsr"] = prestoc.binaryprops_rpsr_get
    __setmethods__["pbin"] = prestoc.binaryprops_pbin_set
    __getmethods__["pbin"] = prestoc.binaryprops_pbin_get
    __setmethods__["rbin"] = prestoc.binaryprops_rbin_set
    __getmethods__["rbin"] = prestoc.binaryprops_rbin_get
    __setmethods__["z"] = prestoc.binaryprops_z_set
    __getmethods__["z"] = prestoc.binaryprops_z_get
    __setmethods__["asinic"] = prestoc.binaryprops_asinic_set
    __getmethods__["asinic"] = prestoc.binaryprops_asinic_get
    __setmethods__["rdetect"] = prestoc.binaryprops_rdetect_set
    __getmethods__["rdetect"] = prestoc.binaryprops_rdetect_get
    __setmethods__["nfftbins"] = prestoc.binaryprops_nfftbins_set
    __getmethods__["nfftbins"] = prestoc.binaryprops_nfftbins_get
    __setmethods__["lowbin"] = prestoc.binaryprops_lowbin_set
    __getmethods__["lowbin"] = prestoc.binaryprops_lowbin_get
    __setmethods__["ppsrerr"] = prestoc.binaryprops_ppsrerr_set
    __getmethods__["ppsrerr"] = prestoc.binaryprops_ppsrerr_get
    __setmethods__["fpsrerr"] = prestoc.binaryprops_fpsrerr_set
    __getmethods__["fpsrerr"] = prestoc.binaryprops_fpsrerr_get
    __setmethods__["rpsrerr"] = prestoc.binaryprops_rpsrerr_set
    __getmethods__["rpsrerr"] = prestoc.binaryprops_rpsrerr_get
    __setmethods__["pbinerr"] = prestoc.binaryprops_pbinerr_set
    __getmethods__["pbinerr"] = prestoc.binaryprops_pbinerr_get
    __setmethods__["rbinerr"] = prestoc.binaryprops_rbinerr_set
    __getmethods__["rbinerr"] = prestoc.binaryprops_rbinerr_get
    __setmethods__["zerr"] = prestoc.binaryprops_zerr_set
    __getmethods__["zerr"] = prestoc.binaryprops_zerr_get
    __setmethods__["asinicerr"] = prestoc.binaryprops_asinicerr_set
    __getmethods__["asinicerr"] = prestoc.binaryprops_asinicerr_get
    __setmethods__["rdetecterr"] = prestoc.binaryprops_rdetecterr_set
    __getmethods__["rdetecterr"] = prestoc.binaryprops_rdetecterr_get
    __setmethods__["sig"] = prestoc.binaryprops_sig_set
    __getmethods__["sig"] = prestoc.binaryprops_sig_get
    __setmethods__["phs"] = prestoc.binaryprops_phs_set
    __getmethods__["phs"] = prestoc.binaryprops_phs_get
    __setmethods__["phserr"] = prestoc.binaryprops_phserr_set
    __getmethods__["phserr"] = prestoc.binaryprops_phserr_get
    __setmethods__["cen"] = prestoc.binaryprops_cen_set
    __getmethods__["cen"] = prestoc.binaryprops_cen_get
    __setmethods__["cenerr"] = prestoc.binaryprops_cenerr_set
    __getmethods__["cenerr"] = prestoc.binaryprops_cenerr_get
    __setmethods__["pur"] = prestoc.binaryprops_pur_set
    __getmethods__["pur"] = prestoc.binaryprops_pur_get
    __setmethods__["purerr"] = prestoc.binaryprops_purerr_set
    __getmethods__["purerr"] = prestoc.binaryprops_purerr_get
    __setmethods__["pow"] = prestoc.binaryprops_pow_set
    __getmethods__["pow"] = prestoc.binaryprops_pow_get
    __setmethods__["powerr"] = prestoc.binaryprops_powerr_set
    __getmethods__["powerr"] = prestoc.binaryprops_powerr_get
    def __init__(self,*args):
        self.this = apply(prestoc.new_binaryprops,args)
        self.thisown = 1
    def __del__(self,prestoc=prestoc):
        if getattr(self,'thisown',0):
            prestoc.delete_binaryprops(self)
    def __repr__(self):
        return "<C binaryprops instance at %s>" % (self.this,)

class binarypropsPtr(binaryprops):
    def __init__(self,this):
        self.this = this
        if not hasattr(self,"thisown"): self.thisown = 0
        self.__class__ = binaryprops
prestoc.binaryprops_swigregister(binarypropsPtr)
class rawbincand:
    __setmethods__ = {}
    for _s in []: __setmethods__.update(_s.__setmethods__)
    def __setattr__(self,name,value):
        if (name == "this"):
            if isinstance(value,rawbincand):
                self.__dict__[name] = value.this
                if hasattr(value,"thisown"): self.__dict__["thisown"] = value.thisown
                del value.thisown
                return
        method = rawbincand.__setmethods__.get(name,None)
        if method: return method(self,value)
        self.__dict__[name] = value

    __getmethods__ = {}
    for _s in []: __getmethods__.update(_s.__getmethods__)
    def __getattr__(self,name):
        method = rawbincand.__getmethods__.get(name,None)
        if method: return method(self)
        raise AttributeError,name

    __setmethods__["full_N"] = prestoc.rawbincand_full_N_set
    __getmethods__["full_N"] = prestoc.rawbincand_full_N_get
    __setmethods__["full_T"] = prestoc.rawbincand_full_T_set
    __getmethods__["full_T"] = prestoc.rawbincand_full_T_get
    __setmethods__["full_lo_r"] = prestoc.rawbincand_full_lo_r_set
    __getmethods__["full_lo_r"] = prestoc.rawbincand_full_lo_r_get
    __setmethods__["mini_N"] = prestoc.rawbincand_mini_N_set
    __getmethods__["mini_N"] = prestoc.rawbincand_mini_N_get
    __setmethods__["mini_r"] = prestoc.rawbincand_mini_r_set
    __getmethods__["mini_r"] = prestoc.rawbincand_mini_r_get
    __setmethods__["mini_power"] = prestoc.rawbincand_mini_power_set
    __getmethods__["mini_power"] = prestoc.rawbincand_mini_power_get
    __setmethods__["mini_numsum"] = prestoc.rawbincand_mini_numsum_set
    __getmethods__["mini_numsum"] = prestoc.rawbincand_mini_numsum_get
    __setmethods__["mini_sigma"] = prestoc.rawbincand_mini_sigma_set
    __getmethods__["mini_sigma"] = prestoc.rawbincand_mini_sigma_get
    __setmethods__["psr_p"] = prestoc.rawbincand_psr_p_set
    __getmethods__["psr_p"] = prestoc.rawbincand_psr_p_get
    __setmethods__["orb_p"] = prestoc.rawbincand_orb_p_set
    __getmethods__["orb_p"] = prestoc.rawbincand_orb_p_get
    def __init__(self,*args):
        self.this = apply(prestoc.new_rawbincand,args)
        self.thisown = 1
    def __del__(self,prestoc=prestoc):
        if getattr(self,'thisown',0):
            prestoc.delete_rawbincand(self)
    def __repr__(self):
        return "<C rawbincand instance at %s>" % (self.this,)

class rawbincandPtr(rawbincand):
    def __init__(self,this):
        self.this = this
        if not hasattr(self,"thisown"): self.thisown = 0
        self.__class__ = rawbincand
prestoc.rawbincand_swigregister(rawbincandPtr)
class foldstats:
    __setmethods__ = {}
    for _s in []: __setmethods__.update(_s.__setmethods__)
    def __setattr__(self,name,value):
        if (name == "this"):
            if isinstance(value,foldstats):
                self.__dict__[name] = value.this
                if hasattr(value,"thisown"): self.__dict__["thisown"] = value.thisown
                del value.thisown
                return
        method = foldstats.__setmethods__.get(name,None)
        if method: return method(self,value)
        self.__dict__[name] = value

    __getmethods__ = {}
    for _s in []: __getmethods__.update(_s.__getmethods__)
    def __getattr__(self,name):
        method = foldstats.__getmethods__.get(name,None)
        if method: return method(self)
        raise AttributeError,name

    __setmethods__["numdata"] = prestoc.foldstats_numdata_set
    __getmethods__["numdata"] = prestoc.foldstats_numdata_get
    __setmethods__["data_avg"] = prestoc.foldstats_data_avg_set
    __getmethods__["data_avg"] = prestoc.foldstats_data_avg_get
    __setmethods__["data_var"] = prestoc.foldstats_data_var_set
    __getmethods__["data_var"] = prestoc.foldstats_data_var_get
    __setmethods__["numprof"] = prestoc.foldstats_numprof_set
    __getmethods__["numprof"] = prestoc.foldstats_numprof_get
    __setmethods__["prof_avg"] = prestoc.foldstats_prof_avg_set
    __getmethods__["prof_avg"] = prestoc.foldstats_prof_avg_get
    __setmethods__["prof_var"] = prestoc.foldstats_prof_var_set
    __getmethods__["prof_var"] = prestoc.foldstats_prof_var_get
    __setmethods__["redchi"] = prestoc.foldstats_redchi_set
    __getmethods__["redchi"] = prestoc.foldstats_redchi_get
    def __init__(self,*args):
        self.this = apply(prestoc.new_foldstats,args)
        self.thisown = 1
    def __del__(self,prestoc=prestoc):
        if getattr(self,'thisown',0):
            prestoc.delete_foldstats(self)
    def __repr__(self):
        return "<C foldstats instance at %s>" % (self.this,)

class foldstatsPtr(foldstats):
    def __init__(self,this):
        self.this = this
        if not hasattr(self,"thisown"): self.thisown = 0
        self.__class__ = foldstats
prestoc.foldstats_swigregister(foldstatsPtr)
tofloatvector = prestoc.tofloatvector

float_to_complex = prestoc.float_to_complex

complex_to_float = prestoc.complex_to_float

SQRT2 = prestoc.SQRT2
PI = prestoc.PI
TWOPI = prestoc.TWOPI
DEGTORAD = prestoc.DEGTORAD
RADTODEG = prestoc.RADTODEG
PIBYTWO = prestoc.PIBYTWO
SOL = prestoc.SOL
SECPERJULYR = prestoc.SECPERJULYR
SECPERDAY = prestoc.SECPERDAY
LOWACC = prestoc.LOWACC
HIGHACC = prestoc.HIGHACC
INTERBIN = prestoc.INTERBIN
INTERPOLATE = prestoc.INTERPOLATE
NO_CHECK_ALIASED = prestoc.NO_CHECK_ALIASED
CHECK_ALIASED = prestoc.CHECK_ALIASED
CONV = prestoc.CONV
CORR = prestoc.CORR
INPLACE_CONV = prestoc.INPLACE_CONV
INPLACE_CORR = prestoc.INPLACE_CORR
FFTDK = prestoc.FFTDK
FFTD = prestoc.FFTD
FFTK = prestoc.FFTK
NOFFTS = prestoc.NOFFTS
RAW = prestoc.RAW
PREPPED = prestoc.PREPPED
FFT = prestoc.FFT
SAME = prestoc.SAME
power_arr = prestoc.power_arr

phase_arr = prestoc.phase_arr

dpower_arr = prestoc.dpower_arr

dphase_arr = prestoc.dphase_arr

num_psrs_in_database = prestoc.num_psrs_in_database

get_psrdata = prestoc.get_psrdata

get_psrdata_by_num = prestoc.get_psrdata_by_num

return_psrparams_at_epoch = prestoc.return_psrparams_at_epoch

get_psr_at_epoch = prestoc.get_psr_at_epoch

readinf = prestoc.readinf

writeinf = prestoc.writeinf

read_mak_input = prestoc.read_mak_input

read_mak_file = prestoc.read_mak_file

write_mak_file = prestoc.write_mak_file

read_int = prestoc.read_int

read_double = prestoc.read_double

frotate = prestoc.frotate

drotate = prestoc.drotate

dorbint = prestoc.dorbint

keplars_eqn = prestoc.keplars_eqn

lin_interp_E = prestoc.lin_interp_E

E_to_phib = prestoc.E_to_phib

E_to_v = prestoc.E_to_v

E_to_p = prestoc.E_to_p

E_to_z = prestoc.E_to_z

E_to_phib_BT = prestoc.E_to_phib_BT

r_resp_halfwidth = prestoc.r_resp_halfwidth

z_resp_halfwidth = prestoc.z_resp_halfwidth

w_resp_halfwidth = prestoc.w_resp_halfwidth

binary_velocity = prestoc.binary_velocity

bin_resp_halfwidth = prestoc.bin_resp_halfwidth

gen_r_response = prestoc.gen_r_response

gen_z_response = prestoc.gen_z_response

gen_w_response = prestoc.gen_w_response

gen_bin_response = prestoc.gen_bin_response

get_numphotons = prestoc.get_numphotons

get_localpower = prestoc.get_localpower

get_localpower3d = prestoc.get_localpower3d

get_derivs3d = prestoc.get_derivs3d

calc_props = prestoc.calc_props

calc_binprops = prestoc.calc_binprops

calc_rzwerrs = prestoc.calc_rzwerrs

candidate_sigma = prestoc.candidate_sigma

power_for_sigma = prestoc.power_for_sigma

chisqr = prestoc.chisqr

print_candidate = prestoc.print_candidate

print_bin_candidate = prestoc.print_bin_candidate

read_rzw_cand = prestoc.read_rzw_cand

read_bin_cand = prestoc.read_bin_cand

read_rawbin_cand = prestoc.read_rawbin_cand

get_rzw_cand = prestoc.get_rzw_cand

get_bin_cand = prestoc.get_bin_cand

get_rawbin_cand = prestoc.get_rawbin_cand

chkfilelen = prestoc.chkfilelen

read_fcomplex_file = prestoc.read_fcomplex_file

read_float_file = prestoc.read_float_file

prune_powers = prestoc.prune_powers

median = prestoc.median

dms2rad = prestoc.dms2rad

hms2rad = prestoc.hms2rad

sphere_ang_diff = prestoc.sphere_ang_diff

spread_with_pad = prestoc.spread_with_pad

spread_no_pad = prestoc.spread_no_pad

paddata = prestoc.paddata

place_complex_kernel = prestoc.place_complex_kernel

place_real_kernel = prestoc.place_real_kernel

chop_complex_ends = prestoc.chop_complex_ends

chop_real_ends = prestoc.chop_real_ends

complex_corr_conv = prestoc.complex_corr_conv

real_corr_conv = prestoc.real_corr_conv

corr_complex = prestoc.corr_complex

stretch_fft = prestoc.stretch_fft

corr_loc_pow = prestoc.corr_loc_pow

rz_interp = prestoc.rz_interp

max_r_arr = prestoc.max_r_arr

max_rz_arr = prestoc.max_rz_arr

fold_errors = prestoc.fold_errors

foldfile = prestoc.foldfile

fold = prestoc.fold

simplefold = prestoc.simplefold

doppler = prestoc.doppler

search_minifft = prestoc.search_minifft

print_rawbincand = prestoc.print_rawbincand

barycenter = prestoc.barycenter

fftwcall = prestoc.fftwcall

tablesixstepfft = prestoc.tablesixstepfft

realfft = prestoc.realfft

tree_max_dm = prestoc.tree_max_dm

smearing_from_bw = prestoc.smearing_from_bw

delay_from_dm = prestoc.delay_from_dm

dm_from_delay = prestoc.dm_from_delay

dedisp_delays = prestoc.dedisp_delays

subband_search_delays = prestoc.subband_search_delays

nice_output_1 = prestoc.nice_output_1

nice_output_2 = prestoc.nice_output_2

corr_rz_plane = prestoc.corr_rz_plane

corr_rz_interp = prestoc.corr_rz_interp




#-------------- Extra Stuff to Make Things Easier -----------------

import math, umath, Numeric, Pgplot, string, numpyio, miscutils

def read_foldstats(file):
   stats = foldstats()
   stats.numdata = read_double(file)
   stats.data_avg = read_double(file)
   stats.data_var = read_double(file)
   stats.numprof = read_double(file)
   stats.prof_avg = read_double(file)
   stats.prof_var = read_double(file)
   stats.redchi = read_double(file)
   return stats
  
class pfd:
   def __init__(self, filename):
      infile = open(filename, "rb")
      self.npart = int(read_double(infile))
      self.nsub = int(read_double(infile))
      self.proflen = int(read_double(infile))
      self.stats = []
      self.profs = Numeric.zeros(self.npart * self.nsub *
                                 self.proflen, 'd')
      if (self.nsub > 1):
         self.profs.shape = (self.npart, self.nsub, self.proflen)
      else:
         self.profs.shape = (self.npart, self.proflen)
      for ii in xrange(self.npart):
         if (self.nsub > 1):
            self.stats.append([])
            for jj in xrange(self.nsub):
               self.stats[ii].append(read_foldstats(infile))
               self.profs[ii][jj] = self.profs[ii][jj] + \
                                    numpyio.fread(infile,
                                                  self.proflen, 'd')
         else:
            self.stats.append(read_foldstats(infile))
            self.profs[ii] = self.profs[ii]+ \
                             numpyio.fread(infile, self.proflen, 'd')
      infile.close()
   
def val_with_err(value, error, len=0, digits=2):
   """
   val_with_err(value, error, len=0, digits=2):
       Returns a string of length len (auto if 0) with 'value'
          rounded to the appropriate decimal place and the
          'error' in parenthesis as in scientific journals.
          The error has 'digits' decimal places.    
       Notes:
          'len' should be ~20 to show full double precision          
             if the base 10 exponent of the error needs to be shown.       
          If len == 0, left-justified minimum length string is returned.
          If len > 0, the string returned is right justified.       
          If len < 0, the string returned is left justified.       
   """
   slen = 40
   if abs(len) > slen: slen = abs(len)
   if digits==2:
      return nice_output_2(' '*slen, value, error, len)
   else:
      return nice_output_1(' '*slen, value, error, len)

def read_inffile(filename):
   """
   read_inffile(filename):
       Return an infodata 'C' structure containing the data from the
       'inf' file in 'filename'.  'filename' should not include the
       '.inf' suffix.
   """
   id = infodata()
   print "Reading information from", "\""+filename+".inf\""
   readinf(id, filename)
   return id

def read_makfile(filename):
   """
   read_makfile(filename):
       Return an makedata 'C' structure containing the data from the
       'mak' in 'filename'.  'filename' should not include the
       '.mak' suffix.
   """
   md = makedata()
   read_mak_file(filename, md)
   return md

def psrepoch(psrname, epoch):
   """
   psrepoch(psrname, epoch):
       Return a psrparams 'C' structure which includes data for
           PSR 'psrname' (a string of the B1950 or J2000 name of the
           pulsar -- without PSR, J, or B included) at epoch 'epoch'
           (in MJD format).
   """
   pp = psrparams()
   num = return_psrparams_at_epoch(pp, psrname, epoch)
   print 'Retrieved data at MJD %f for %s' % (epoch, pp.jname)
   print 'The pulsar was #%d in the database.' % num
   return pp

def collect_psrdata():
    """
    collect_psrdata():
        Return a list of all of the pulsars in the Taylor et al.
            pulsar database including their characteristics.
    """
    pdata = []
    np = num_psrs_in_database()
    print 'There are %d pulsars in the database.' % np
    for i in range(0, np):
        pdata.append(psrdata())
        get_psrdata_by_num(pdata[i], i)
    return pdata

def read_rzwcands(filename):
    """
    read_rzwcands(filename):
        Return a list of all of the rzw search candidates from
            the file 'filename'.
    """
    infile = open(filename, "r")
    cands = []
    nextcand = fourierprops()
    while (read_rzw_cand(infile, nextcand)):
       cands.append(nextcand)
       nextcand = fourierprops()
    infile.close()
    return cands

def read_rawbincands(filename):
    """
    read_rawbincands(filename):
        Return a list of all of the raw binary search candidates
            from the file 'filename'.
    """
    infile = open(filename, "r")
    cands = []
    nextcand = rawbincand()
    while (read_rawbin_cand(infile, nextcand)):
       cands.append(nextcand)
       nextcand = rawbincand()
    infile.close()
    return cands

def next2_to_n(x):
    """
    next2_to_n(x):
        Return the first value of 2^n >= x.
    """
    i = 1L
    while (i < x): i = i << 1
    return i

def rfft(data, sign=-1):
   """
   rfft(data, sign=-1):
       Return the FFT of the real-valued 'data'.
       Note:  This only returns the positive frequency half of the FFT,
              since the other half is symmetric.  The Nyquist frequency
              is stored in the complex part of frequency 0 as per
              Numerical Recipes.
       The optional value 'sign' should be positive or negative 1.
   """
   # Default to sign = -1 if the user gives a bad value
   tmp = Numeric.array(data, copy=1)
   if (sign == -1 or sign != 1):
      tmp = tofloatvector(tmp)
      realfft(tmp, len(tmp), -1)
      float_to_complex(tmp)
   else:
      complex_to_float(tmp)
      realfft(tmp, len(tmp), 1)
   return tmp

def spectralpower(fftarray):
    """
    spectralpower(fftarray):
        Return the power spectrum of a complex FFT 'fftarray'.
    """
    fftarray = Numeric.asarray(fftarray)
    if fftarray.typecode()=='F':
       return power_arr(fftarray, len(fftarray))
    elif fftarray.typecode()=='D':
       return dpower_arr(fftarray, len(fftarray))
    else:
       print 'fftarray must be complex in spectralpower()'
       return None
    
def spectralphase(fftarray):
    """
    spectralphase(fftarray):
        Return the spectral phase (deg) of a complex FFT 'fftarray'.
    """
    fftarray = Numeric.asarray(fftarray)
    if fftarray.typecode()=='F':
       return phase_arr(fftarray, len(fftarray))
    elif fftarray.typecode()=='D':
       return dphase_arr(fftarray, len(fftarray))
    else:
       print 'fftarray must be complex in spectralpower()'
       return None

def maximize_rz(data, r, z, norm = None):
   """
   maximize_rz(data, r, z, norm = None):
       Optimize the detection of a signal at location 'r', 'z' in
           the F-Fdot plane.  The routine returns a list containing
           the optimized values of the maximum normalized power, rmax,
           zmax, and an rderivs structure for the peak.
   """
   rd = rderivs()
   (maxpow, rmax, zmax) = max_rz_arr(data, len(data), r, z, rd)
   if not norm:
      maxpow = maxpow / rd.locpow
   else:
      maxpow = maxpow / norm
   return [maxpow, rmax, zmax, rd]

def maximize_r(data, r, norm = None):
   """
   maximize_r(data, r, norm = None):
       Optimize the detection of a signal at Fourier frequency 'r' in
           a FFT 'data'.  The routine returns a list containing
           the optimized values of the maximum normalized power, rmax,
           and an rderivs structure for the peak.
   """
   rd = rderivs()
   (maxpow, rmax) = max_r_arr(data, len(data), r, rd)
   if not norm:
      maxpow = maxpow / rd.locpow
   else:
      maxpow = maxpow / norm
   return [maxpow, rmax, rd]

def search_fft(data, numcands, norm='default'):
   """
   search_fft(data, numcands):
      Search a short FFT and return a list containing the powers and
      Fourier frequencies of the 'numcands' highest candidates in 'data'.
      'norm' is the value to multiply each pow power by to get
         a normalized power spectrum (defaults to  1.0/(Freq 0) value)
   """
   if (norm=='default'): norm = 1.0/data[0].real
   hp = Numeric.zeros(numcands, 'f')
   hf = Numeric.zeros(numcands, 'f')
   search_minifft(data, len(data), norm, numcands, hp, hf) 
   cands = []
   for i in range(numcands):
      cands.append([hp[i],hf[i]])
   return cands

def ffdot_plane(data, r, dr, numr, z, dz, numz):
   """
   ffdot_plane(data, r, dr, numr, z, dz, numz):
       Generate an F-Fdot plane centered on the point 'r', 'z'.
       There will be a total of 'numr' x 'numz' points in the array.
       The F-Fdot plane will be interpolated such the points are
       separated by 'dr' in the 'r' (i.e. f) direction and 'dz'
       in the 'z' (i.e. fdot) direction.  'data' is the input FFT.
       Note:  'dr' much be the reciprocal of an integer
              (i.e. 1 / numbetween).  Also, 'r' is considered to be
              the average frequency (r = ro + z / 2).
   """
   numbetween = int(1.0 / dr)
   startbin = int(r - (numr * dr) / 2)
   loz = z - (numz * dz) / 2
   hiz = loz + (numz - 1) * dz
   maxabsz = max(abs(loz), abs(hiz))
   kern_half_width = z_resp_halfwidth(maxabsz, LOWACC)
   fftlen = next2_to_n(numr + 2 * numbetween * kern_half_width)
   (ffdraw, nextbin) = corr_rz_plane(data, len(data), numbetween,
                                     startbin, loz, hiz, numz,
                                     fftlen, LOWACC)
   return Numeric.array(ffdraw[:,0:numr], copy=1)

def estimate_rz(psr, T, show=0, device='/XWIN'):
    """
    estimate_rz(psr, T, eo=0.0, show=0, device='/XWIN'):
        Return estimates of a pulsar's average Fourier freq ('r')
        relative to its nominal Fourier freq as well as its
        Fourier f-dot ('z') in bins, of a pulsar.
           'psr' is a psrparams structure describing the pulsar.
           'T' is the length of the observation in sec.
           'show' if true, displays plots of 'r' and 'z'.
           'device' if the device to plot to if 'show' is true.
    """
    from Scientific.Statistics import average
    startE = keplars_eqn(psr.orb.t, psr.orb.p, psr.orb.e, 1.0E-15)
    numorbpts = int(T / psr.orb.p + 1.0) * 1024 + 1
    dt = T / (numorbpts - 1)
    E = dorbint(startE, numorbpts, dt, psr.orb)
    z = z_from_e(E, psr, T)
    r = T/p_from_e(E, psr) - T/psr.p
    if show:
        times = Numeric.arange(numorbpts) * dt
        Pgplot.plotxy(r, times, labx = 'Time', \
                      laby = 'Fourier Frequency (r)', device=device)
        if device=='/XWIN':
           print 'Press enter to continue:'
           i = raw_input()
        Pgplot.nextplotpage()
        Pgplot.plotxy(z, times, labx = 'Time',
                      laby = 'Fourier Frequency Derivative (z)', device=device)
        Pgplot.closeplot()
    return (average(r), average(z))
    
def alias(r, rny):
    """
    alias_to_r(r, rny):
        Convert an aliased Fourier frequency into the 'true' Fourier
        frequency of a signal.  Or vise-versa -- the transformation is
        symmetric about the Nyquist Freq.
           'r' is the signal's Fourier frequency to convert.
           'rny' is the Nyquist frequency (in bins).  For an FFT
              of real data, 'rny' = number of data points FFT'd / 2.
    """
    return 2.0 * rny - r

def show_ffdot_plane(data, r, z, dr = 0.125, dz = 0.5,
                     numr = 300, numz = 300, T = None, 
                     contours = None, title = None, 
                     image = "astro", device = "/XWIN", norm = 1.0):
   """
   show_ffdot_plane(data, r, z):
       Show a color plot of the F-Fdot plane centered on the point 'r', 'z'.
   """
   ffdp = ffdot_plane(data, r, dr, numr, z, dz, numz)
   ffdpow = spectralpower(ffdp.flat)
   ffdpow.shape = (numz, numr)
   startbin = int(r - (numr * dr) / 2)
   startz = int(z - (numz * dz) / 2)
   x = Numeric.arange(numr, typecode="d") * dr + startbin
   y = Numeric.arange(numz, typecode="d") * dz + startz
   highpt = Numeric.argmax(ffdpow.flat)
   hir = highpt % numr
   hiz = highpt / numr
   print ""
   print "Fourier Freqs from ", min(x), "to", max(x), "."
   print "Fourier Fdots from ", min(y), "to", max(y), "."
   print "Maximum normalized power is ", ffdpow[hiz][hir]
   print "The max value is located at:  r =", startbin + hir * dr, \
         "  z =", startz + hiz * dz
   print ""
   if not T:
      Pgplot.plot2d(ffdpow, x, y, labx = "Fourier Frequency (bins)", \
                    laby = "Fourier Frequency Derivative", \
                    title = title, image = image, \
                    contours = contours, device = device)
   else:
      Pgplot.plot2d(ffdpow, x/T, y/(T**2.0), labx = "Frequency (hz)", \
                    laby = "Frequency Derivative (Hz/sec)", \
                    rangex2 = [x[0], x[-1]], rangey2 = [y[0], y[-1]], \
                    labx2 = "Fourier Frequency", \
                    laby2 = "Fourier Frequency Derivative", \
                    title = title, image = image, \
                    contours = contours, device = device)


def v_from_e(e, psr):
   """
   v_from_e(e, psr):
       Return a vector of velocities (km/s) from a vector of Eccentric
       anomalys.
           'e' is the vector of Eccentric anomalys.
           'psr' is a psrparams instance containing info about the pulsar.
   """
   oldw = psr.orb.w
   psr.orb.w = psr.orb.w * DEGTORAD
   v = Numeric.array(e, copy=1)
   E_to_v(v, len(v), psr.orb)
   psr.orb.w = oldw
   return v

def d_from_e(e, psr):
   """
   d_from_e(e, psr):
       Return a vector of time delays (s) from a vector of Eccentric
       anomalys.
           'e' is the vector of Eccentric anomalys.
           'psr' is a psrparams instance containing info about the pulsar.
   """
   oldw = psr.orb.w
   psr.orb.w = psr.orb.w * DEGTORAD
   d = Numeric.array(e, copy=1)
   E_to_phib(d, len(d), psr.orb)
   psr.orb.w = oldw
   return d

def p_from_e(e, psr):
   """
   p_from_e(e, psr):
       Return a vector of pulsar periods (s) from a vector of Eccentric
       anomalys.
           'e' is the vector of Eccentric anomalys.
           'psr' is a psrparams instance containing info about the pulsar.
   """
   oldw = psr.orb.w
   psr.orb.w = psr.orb.w * DEGTORAD
   p = Numeric.array(e, copy=1)
   E_to_p(p, len(p), psr.p, psr.orb)
   psr.orb.w = oldw
   return p

def z_from_e(e, psr, T):
   """
   z_from_e(e, psr):
       Return a vector of Fourier F-dots (bins) from a vector of Eccentric
       anomalys.
           'e' is the vector of Eccentric anomalys.
           'psr' is a psrparams instance containing info about the pulsar.
           'T' is the total length of the observation (s).
   """
   oldw = psr.orb.w
   psr.orb.w = psr.orb.w * DEGTORAD
   z = Numeric.array(e, copy=1)
   E_to_z(z, len(z), psr.p, T, psr.orb)
   psr.orb.w = oldw
   return z

def pcorr(data, kernel, numbetween, lo, hi):
   """
   pcorr(data, kernel, numbetween, lo, hi):
       Perform a correlation with the raw complex vectors 'data' and
       'kernel'.  The returned vector should start at frequency
       'lo' (must be an integer), and go up to but not include 'hi'
       (also an integer).
   """
   kern_half_width = len(kernel)/(2 * numbetween)
   result = Numeric.zeros((hi-lo)*numbetween, 'F')
   corr_complex(data, len(data), RAW,
                kernel, len(kernel), RAW,
                result, len(result), lo,
                numbetween, kern_half_width, CORR)
   return result

def p_to_f(p, pd, pdd):
   """
   p_to_f(p, pd, pdd):
      Convert period, period derivative and period second
      derivative to the equivalent frequency counterparts.
      Will also convert from f to p.
   """
   f = 1.0 / p
   fd = -pd / (p * p)
   if (pdd==0.0):
      fdd = 0.0
   else:
      fdd = 2.0 * pd * pd / (p**3.0) - pdd / (p * p)
   return [f, fd, fdd]

def bary_to_topo(pb, pbd, pbdd, infofilenm, ephem="DE200"):
   """
   bary_to_topo(pb, pbd, pbdd, infofilenm, ephem="DE200"):
      Use least squares to calculate topocentric period
      period derivative, and period second derivative
      for the corresponding barycentric values.  The data
      for the observation must be found in the info file.
   """
   from LinearAlgebra import linear_least_squares
   if infofilenm[-4:]==".inf":  infofilenm = infofilenm[:-4]
   obs = read_inffile(infofilenm)
   T = obs.N * obs.dt
   dt = 10.0
   tto = obs.mjd_i + obs.mjd_f
   tts = Numeric.arange(tto, tto + (T + dt) / SECPERDAY, dt / SECPERDAY)
   nn = len(tts)
   bts = Numeric.zeros(nn, 'd')
   vel = Numeric.zeros(nn, 'd')
   ra = miscutils.coord_to_string(obs.ra_h, obs.ra_m, obs.ra_s)
   dec = miscutils.coord_to_string(obs.dec_d, obs.dec_m, obs.dec_s)
   if (obs.telescope == 'Parkes'):  tel = 'PK'
   elif (obs.telescope == 'Effelsberg'):  tel = 'EB'
   elif (obs.telescope == 'Arecibo'):  tel = 'AO'
   elif (obs.telescope == 'MMT'):  tel = 'MT'
   else:
      print "Telescope not recognized."
      return 0
   barycenter(tts, bts, vel, nn, ra, dec, tel, ephem)
   print "Topocentric start time = %17.11f" % tts[0]
   print "Barycentric start time = %17.11f" % bts[0]
   avgvel = Numeric.add.reduce(vel) / nn
   print "Average Earth velocity = %10.5e c" % (avgvel)
   tts = Numeric.arange(nn, typecode='d') * dt
   bts = (bts - bts[0]) * SECPERDAY
   [fb, fbd, fbdd] = p_to_f(pb, pbd, pbdd)
   b = fb * bts + fbd * bts**2.0 / 2.0 + fbdd * bts**3.0 / 6.0
   a = Numeric.transpose(Numeric.asarray([tts, tts**2.0, tts**3.0]))
   [ft, ftd, ftdd], residuals, rank, sv = linear_least_squares(a,b)
   [pt, ptd, ptdd] = p_to_f(ft, ftd, ftdd)
   print "    Topocentric period = %15.12f" % pt
   print "     Topocentric p-dot = %15.9e" % ptd
   print "  Topocentric p-dotdot = %15.9e" % ptdd
   print "     Quick Topo period = %15.12f" % (pb * (1.0 + avgvel))
   print "      Quick Topo p-dot = %15.9e" % (pbd * (1.0 + avgvel))
   print "   Quick Topo p-dotdot = %15.9e" % (pbdd * (1.0 + avgvel))
   return [pt, ptd, ptdd]


def measure_phase(profile, template, sigma, fwhm):
    """
    measure_phase(profile, template, sigma, fwhm):
       TOA measurement technique from J. H. Taylor's talk
       _Pulsar_Timing_and_Relativistic_Gravity_.  Routine
       takes two profiles, the first measured and the
       second a high S/N template and determines the phase
       offset of 'profile' from 'template'.  Both profiles
       must have the same number of points.  'sigma' denotes
       the RMS noise level of the 'profile'.  'fwhm' is the
       approximate width of the template pulse (0-1).  The phase
       returned is cyclic (i.e. from 0-1).  The routine
       returns a tuple comtaining (tau, tau_err, b, b_err, a).
       Where 'tau' is the phase, 'b' is the scaling factor,
       and 'a' is the DC offset.  The error values are
       estimates of the 1 sigma errors.
    """
    from simple_roots import newton_raphson
    N = len(profile)
    if not (N == len(template)):
       print "Lengths of 'profile' and 'template' must"
       print "  be equal in measure_phase()."
       return 0.0
    ft = rfft(profile)
    p0 = ft[0].real
    # Nyquist freq
    ft[0] = complex(ft[0].imag, 0.0)
    P_k = abs(ft)
    frotate(P_k, len(ft), 1)
    Theta_k = umath.arctan2(-ft.imag, ft.real)
    frotate(Theta_k, len(ft), 1)
    ft = rfft(template)
    s0 = ft[0].real
    # Nyquist freq
    ft[0] = complex(ft[0].imag, 0.0)
    S_k = abs(ft)
    frotate(S_k, len(ft), 1)
    Phi_k = umath.arctan2(-ft.imag, ft.real)
    frotate(Phi_k, len(ft), 1)
    # Estimate of the noise sigma (This needs to be checked)
    # Note:  Checked 10 Jul 2000.  Looks OK.
    sig = sigma * math.sqrt(N)
    k = Numeric.arange(len(ft), typecode='d') + 1.0
    def fn(tau, k=k, p=P_k, s=S_k, theta=Theta_k, phi=Phi_k):
       # Since Nyquist freq always has phase = 0.0
       k[-1] = 0.0
       return Numeric.add.reduce(k * p * s *
                                 umath.sin(phi - theta + k * tau))
    def dfn(tau, k=k, p=P_k, s=S_k, theta=Theta_k, phi=Phi_k):
       # Since Nyquist freq always has phase = 0.0
       k[-1] = 0.0
       return Numeric.add.reduce(k * k * p * s *
                                 umath.cos(phi - theta + k * tau))
    numphases = 200
    ddchidt = Numeric.zeros(numphases, 'd')
    phases = Numeric.arange(numphases, typecode='d') / \
             float(numphases-1) * TWOPI - PI
    for i in Numeric.arange(numphases):
       ddchidt[i] = dfn(phases[i])
    maxdphase = phases[Numeric.argmax(ddchidt)] + \
                0.5 * TWOPI / (numphases - 1.0)
    # Solve for tau
    tau = newton_raphson(fn, dfn, maxdphase - 0.5 * fwhm * TWOPI,
                         maxdphase + 0.5 * fwhm * TWOPI)
    # Solve for b
    c = P_k * S_k * umath.cos(Phi_k - Theta_k + k * tau)
    d = Numeric.add.reduce(S_k**2.0)
    b = Numeric.add.reduce(c) / d
    # tau sigma
    tau_err = sig * umath.sqrt(1.0 / (2.0 * b *
                                      Numeric.add.reduce(k**2.0 * c)))
    # b sigma  (Note:  This seems to be an underestimate...)
    b_err = sig * umath.sqrt(1.0 / (2.0 * d))
    # Solve for a
    a = (p0 - b * s0) / float(N)
    return (tau / TWOPI, tau_err / TWOPI, b, b_err, a)
