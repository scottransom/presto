# This file was created automatically by SWIG.
import prestoc
class orbitparams:
    def __init__(self,*args):
        self.this = apply(prestoc.new_orbitparams,args)
        self.thisown = 1

    def __del__(self,prestoc=prestoc):
        if getattr(self,'thisown',0):
            prestoc.delete_orbitparams(self)
    __setmethods__ = {}
    for _s in []: __setmethods__.update(_s.__setmethods__)
    __setmethods__.update({
        "p" : prestoc.orbitparams_p_set,
        "e" : prestoc.orbitparams_e_set,
        "x" : prestoc.orbitparams_x_set,
        "w" : prestoc.orbitparams_w_set,
        "t" : prestoc.orbitparams_t_set,
        "pd" : prestoc.orbitparams_pd_set,
        "wd" : prestoc.orbitparams_wd_set,
    })
    def __setattr__(self,name,value):
        if (name == "this") or (name == "thisown"): self.__dict__[name] = value; return
        method = orbitparams.__setmethods__.get(name,None)
        if method: return method(self,value)
        self.__dict__[name] = value
    __getmethods__ = {}
    for _s in []: __getmethods__.update(_s.__getmethods__)
    __getmethods__.update({
        "p" : prestoc.orbitparams_p_get,
        "e" : prestoc.orbitparams_e_get,
        "x" : prestoc.orbitparams_x_get,
        "w" : prestoc.orbitparams_w_get,
        "t" : prestoc.orbitparams_t_get,
        "pd" : prestoc.orbitparams_pd_get,
        "wd" : prestoc.orbitparams_wd_get,
    })
    def __getattr__(self,name):
        method = orbitparams.__getmethods__.get(name,None)
        if method: return method(self)
        raise AttributeError,name
    def __repr__(self):
        return "<C orbitparams instance at %s>" % (self.this,)
class orbitparamsPtr(orbitparams):
    def __init__(self,this):
        self.this = this
        self.thisown = 0
        self.__class__ = orbitparams



class PSRDATA:
    def __init__(self,*args):
        self.this = apply(prestoc.new_PSRDATA,args)
        self.thisown = 1

    def __del__(self,prestoc=prestoc):
        if getattr(self,'thisown',0):
            prestoc.delete_psrdata(self)
    __setmethods__ = {}
    for _s in []: __setmethods__.update(_s.__setmethods__)
    __setmethods__.update({
        "ra2000" : prestoc.PSRDATA_ra2000_set,
        "ra1950" : prestoc.PSRDATA_ra1950_set,
        "rae" : prestoc.PSRDATA_rae_set,
        "dec2000" : prestoc.PSRDATA_dec2000_set,
        "dec1950" : prestoc.PSRDATA_dec1950_set,
        "dece" : prestoc.PSRDATA_dece_set,
        "dmin__" : prestoc.PSRDATA_dmin___set,
        "dmax__" : prestoc.PSRDATA_dmax___set,
        "dist" : prestoc.PSRDATA_dist_set,
        "ldeg" : prestoc.PSRDATA_ldeg_set,
        "bdeg" : prestoc.PSRDATA_bdeg_set,
        "pmra" : prestoc.PSRDATA_pmra_set,
        "pmrae" : prestoc.PSRDATA_pmrae_set,
        "pmdec" : prestoc.PSRDATA_pmdec_set,
        "pmdece" : prestoc.PSRDATA_pmdece_set,
        "posepoch" : prestoc.PSRDATA_posepoch_set,
        "p" : prestoc.PSRDATA_p_set,
        "pe" : prestoc.PSRDATA_pe_set,
        "pdot" : prestoc.PSRDATA_pdot_set,
        "pdote" : prestoc.PSRDATA_pdote_set,
        "f2" : prestoc.PSRDATA_f2_set,
        "f2e" : prestoc.PSRDATA_f2e_set,
        "f3" : prestoc.PSRDATA_f3_set,
        "f3e" : prestoc.PSRDATA_f3e_set,
        "epoch" : prestoc.PSRDATA_epoch_set,
        "dm" : prestoc.PSRDATA_dm_set,
        "dme" : prestoc.PSRDATA_dme_set,
        "rm" : prestoc.PSRDATA_rm_set,
        "rme" : prestoc.PSRDATA_rme_set,
        "we" : prestoc.PSRDATA_we_set,
        "w50" : prestoc.PSRDATA_w50_set,
        "w10" : prestoc.PSRDATA_w10_set,
        "s400" : prestoc.PSRDATA_s400_set,
        "s600" : prestoc.PSRDATA_s600_set,
        "s1400" : prestoc.PSRDATA_s1400_set,
        "tau" : prestoc.PSRDATA_tau_set,
        "t408" : prestoc.PSRDATA_t408_set,
        "distmod" : prestoc.PSRDATA_distmod_set,
        "lum" : prestoc.PSRDATA_lum_set,
        "bsurf" : prestoc.PSRDATA_bsurf_set,
        "age" : prestoc.PSRDATA_age_set,
        "edot" : prestoc.PSRDATA_edot_set,
        "pb" : prestoc.PSRDATA_pb_set,
        "pbe" : prestoc.PSRDATA_pbe_set,
        "a1" : prestoc.PSRDATA_a1_set,
        "a1e" : prestoc.PSRDATA_a1e_set,
        "om" : prestoc.PSRDATA_om_set,
        "ome" : prestoc.PSRDATA_ome_set,
        "omdot" : prestoc.PSRDATA_omdot_set,
        "omdote" : prestoc.PSRDATA_omdote_set,
        "e" : prestoc.PSRDATA_e_set,
        "ee" : prestoc.PSRDATA_ee_set,
        "t0" : prestoc.PSRDATA_t0_set,
        "t0e" : prestoc.PSRDATA_t0e_set,
        "gamma" : prestoc.PSRDATA_gamma_set,
        "gammae" : prestoc.PSRDATA_gammae_set,
        "pbdot" : prestoc.PSRDATA_pbdot_set,
        "pbdote" : prestoc.PSRDATA_pbdote_set,
        "si" : prestoc.PSRDATA_si_set,
        "sie" : prestoc.PSRDATA_sie_set,
        "r__" : prestoc.PSRDATA_r___set,
        "re" : prestoc.PSRDATA_re_set,
        "pb2" : prestoc.PSRDATA_pb2_set,
        "pb2e" : prestoc.PSRDATA_pb2e_set,
        "a12" : prestoc.PSRDATA_a12_set,
        "a12e" : prestoc.PSRDATA_a12e_set,
        "om2" : prestoc.PSRDATA_om2_set,
        "om2e" : prestoc.PSRDATA_om2e_set,
        "omdot2" : prestoc.PSRDATA_omdot2_set,
        "omdot2e" : prestoc.PSRDATA_omdot2e_set,
        "e2" : prestoc.PSRDATA_e2_set,
        "e2e" : prestoc.PSRDATA_e2e_set,
        "t02" : prestoc.PSRDATA_t02_set,
        "t02e" : prestoc.PSRDATA_t02e_set,
        "gamma2" : prestoc.PSRDATA_gamma2_set,
        "gamma2e" : prestoc.PSRDATA_gamma2e_set,
        "pbdot2" : prestoc.PSRDATA_pbdot2_set,
        "pbdot2e" : prestoc.PSRDATA_pbdot2e_set,
        "si2" : prestoc.PSRDATA_si2_set,
        "si2e" : prestoc.PSRDATA_si2e_set,
        "r2" : prestoc.PSRDATA_r2_set,
        "r2e" : prestoc.PSRDATA_r2e_set,
        "nscode" : prestoc.PSRDATA_nscode_set,
        "ndflag" : prestoc.PSRDATA_ndflag_set,
        "ntauflag" : prestoc.PSRDATA_ntauflag_set,
        "ntype" : prestoc.PSRDATA_ntype_set,
        "modcode" : prestoc.PSRDATA_modcode_set,
        "limcode" : prestoc.PSRDATA_limcode_set,
        "ibin" : prestoc.PSRDATA_ibin_set,
        "jname" : prestoc.PSRDATA_jname_set,
        "bname" : prestoc.PSRDATA_bname_set,
        "lcode" : prestoc.PSRDATA_lcode_set,
        "ucode" : prestoc.PSRDATA_ucode_set,
    })
    def __setattr__(self,name,value):
        if (name == "this") or (name == "thisown"): self.__dict__[name] = value; return
        method = PSRDATA.__setmethods__.get(name,None)
        if method: return method(self,value)
        self.__dict__[name] = value
    __getmethods__ = {}
    for _s in []: __getmethods__.update(_s.__getmethods__)
    __getmethods__.update({
        "ra2000" : prestoc.PSRDATA_ra2000_get,
        "ra1950" : prestoc.PSRDATA_ra1950_get,
        "rae" : prestoc.PSRDATA_rae_get,
        "dec2000" : prestoc.PSRDATA_dec2000_get,
        "dec1950" : prestoc.PSRDATA_dec1950_get,
        "dece" : prestoc.PSRDATA_dece_get,
        "dmin__" : prestoc.PSRDATA_dmin___get,
        "dmax__" : prestoc.PSRDATA_dmax___get,
        "dist" : prestoc.PSRDATA_dist_get,
        "ldeg" : prestoc.PSRDATA_ldeg_get,
        "bdeg" : prestoc.PSRDATA_bdeg_get,
        "pmra" : prestoc.PSRDATA_pmra_get,
        "pmrae" : prestoc.PSRDATA_pmrae_get,
        "pmdec" : prestoc.PSRDATA_pmdec_get,
        "pmdece" : prestoc.PSRDATA_pmdece_get,
        "posepoch" : prestoc.PSRDATA_posepoch_get,
        "p" : prestoc.PSRDATA_p_get,
        "pe" : prestoc.PSRDATA_pe_get,
        "pdot" : prestoc.PSRDATA_pdot_get,
        "pdote" : prestoc.PSRDATA_pdote_get,
        "f2" : prestoc.PSRDATA_f2_get,
        "f2e" : prestoc.PSRDATA_f2e_get,
        "f3" : prestoc.PSRDATA_f3_get,
        "f3e" : prestoc.PSRDATA_f3e_get,
        "epoch" : prestoc.PSRDATA_epoch_get,
        "dm" : prestoc.PSRDATA_dm_get,
        "dme" : prestoc.PSRDATA_dme_get,
        "rm" : prestoc.PSRDATA_rm_get,
        "rme" : prestoc.PSRDATA_rme_get,
        "we" : prestoc.PSRDATA_we_get,
        "w50" : prestoc.PSRDATA_w50_get,
        "w10" : prestoc.PSRDATA_w10_get,
        "s400" : prestoc.PSRDATA_s400_get,
        "s600" : prestoc.PSRDATA_s600_get,
        "s1400" : prestoc.PSRDATA_s1400_get,
        "tau" : prestoc.PSRDATA_tau_get,
        "t408" : prestoc.PSRDATA_t408_get,
        "distmod" : prestoc.PSRDATA_distmod_get,
        "lum" : prestoc.PSRDATA_lum_get,
        "bsurf" : prestoc.PSRDATA_bsurf_get,
        "age" : prestoc.PSRDATA_age_get,
        "edot" : prestoc.PSRDATA_edot_get,
        "pb" : prestoc.PSRDATA_pb_get,
        "pbe" : prestoc.PSRDATA_pbe_get,
        "a1" : prestoc.PSRDATA_a1_get,
        "a1e" : prestoc.PSRDATA_a1e_get,
        "om" : prestoc.PSRDATA_om_get,
        "ome" : prestoc.PSRDATA_ome_get,
        "omdot" : prestoc.PSRDATA_omdot_get,
        "omdote" : prestoc.PSRDATA_omdote_get,
        "e" : prestoc.PSRDATA_e_get,
        "ee" : prestoc.PSRDATA_ee_get,
        "t0" : prestoc.PSRDATA_t0_get,
        "t0e" : prestoc.PSRDATA_t0e_get,
        "gamma" : prestoc.PSRDATA_gamma_get,
        "gammae" : prestoc.PSRDATA_gammae_get,
        "pbdot" : prestoc.PSRDATA_pbdot_get,
        "pbdote" : prestoc.PSRDATA_pbdote_get,
        "si" : prestoc.PSRDATA_si_get,
        "sie" : prestoc.PSRDATA_sie_get,
        "r__" : prestoc.PSRDATA_r___get,
        "re" : prestoc.PSRDATA_re_get,
        "pb2" : prestoc.PSRDATA_pb2_get,
        "pb2e" : prestoc.PSRDATA_pb2e_get,
        "a12" : prestoc.PSRDATA_a12_get,
        "a12e" : prestoc.PSRDATA_a12e_get,
        "om2" : prestoc.PSRDATA_om2_get,
        "om2e" : prestoc.PSRDATA_om2e_get,
        "omdot2" : prestoc.PSRDATA_omdot2_get,
        "omdot2e" : prestoc.PSRDATA_omdot2e_get,
        "e2" : prestoc.PSRDATA_e2_get,
        "e2e" : prestoc.PSRDATA_e2e_get,
        "t02" : prestoc.PSRDATA_t02_get,
        "t02e" : prestoc.PSRDATA_t02e_get,
        "gamma2" : prestoc.PSRDATA_gamma2_get,
        "gamma2e" : prestoc.PSRDATA_gamma2e_get,
        "pbdot2" : prestoc.PSRDATA_pbdot2_get,
        "pbdot2e" : prestoc.PSRDATA_pbdot2e_get,
        "si2" : prestoc.PSRDATA_si2_get,
        "si2e" : prestoc.PSRDATA_si2e_get,
        "r2" : prestoc.PSRDATA_r2_get,
        "r2e" : prestoc.PSRDATA_r2e_get,
        "nscode" : prestoc.PSRDATA_nscode_get,
        "ndflag" : prestoc.PSRDATA_ndflag_get,
        "ntauflag" : prestoc.PSRDATA_ntauflag_get,
        "ntype" : prestoc.PSRDATA_ntype_get,
        "modcode" : prestoc.PSRDATA_modcode_get,
        "limcode" : prestoc.PSRDATA_limcode_get,
        "ibin" : prestoc.PSRDATA_ibin_get,
        "jname" : prestoc.PSRDATA_jname_get,
        "bname" : prestoc.PSRDATA_bname_get,
        "lcode" : prestoc.PSRDATA_lcode_get,
        "ucode" : prestoc.PSRDATA_ucode_get,
    })
    def __getattr__(self,name):
        method = PSRDATA.__getmethods__.get(name,None)
        if method: return method(self)
        raise AttributeError,name
    def __repr__(self):
        return "<C PSRDATA instance at %s>" % (self.this,)
class PSRDATAPtr(PSRDATA):
    def __init__(self,this):
        self.this = this
        self.thisown = 0
        self.__class__ = PSRDATA



class PSRPARAMS:
    def __init__(self,*args):
        self.this = apply(prestoc.new_PSRPARAMS,args)
        self.thisown = 1

    def __del__(self,prestoc=prestoc):
        if getattr(self,'thisown',0):
            prestoc.delete_psrparams(self)
    __setmethods__ = {}
    for _s in []: __setmethods__.update(_s.__setmethods__)
    __setmethods__.update({
        "jname" : prestoc.PSRPARAMS_jname_set,
        "bname" : prestoc.PSRPARAMS_bname_set,
        "ntype" : prestoc.PSRPARAMS_ntype_set,
        "ra2000" : prestoc.PSRPARAMS_ra2000_set,
        "dec2000" : prestoc.PSRPARAMS_dec2000_set,
        "dm" : prestoc.PSRPARAMS_dm_set,
        "dist" : prestoc.PSRPARAMS_dist_set,
        "fwhm" : prestoc.PSRPARAMS_fwhm_set,
        "timepoch" : prestoc.PSRPARAMS_timepoch_set,
        "p" : prestoc.PSRPARAMS_p_set,
        "pd" : prestoc.PSRPARAMS_pd_set,
        "pdd" : prestoc.PSRPARAMS_pdd_set,
        "f" : prestoc.PSRPARAMS_f_set,
        "fd" : prestoc.PSRPARAMS_fd_set,
        "fdd" : prestoc.PSRPARAMS_fdd_set,
        "orb" : prestoc.PSRPARAMS_orb_set,
    })
    def __setattr__(self,name,value):
        if (name == "this") or (name == "thisown"): self.__dict__[name] = value; return
        method = PSRPARAMS.__setmethods__.get(name,None)
        if method: return method(self,value)
        self.__dict__[name] = value
    __getmethods__ = {}
    for _s in []: __getmethods__.update(_s.__getmethods__)
    __getmethods__.update({
        "jname" : prestoc.PSRPARAMS_jname_get,
        "bname" : prestoc.PSRPARAMS_bname_get,
        "ntype" : prestoc.PSRPARAMS_ntype_get,
        "ra2000" : prestoc.PSRPARAMS_ra2000_get,
        "dec2000" : prestoc.PSRPARAMS_dec2000_get,
        "dm" : prestoc.PSRPARAMS_dm_get,
        "dist" : prestoc.PSRPARAMS_dist_get,
        "fwhm" : prestoc.PSRPARAMS_fwhm_get,
        "timepoch" : prestoc.PSRPARAMS_timepoch_get,
        "p" : prestoc.PSRPARAMS_p_get,
        "pd" : prestoc.PSRPARAMS_pd_get,
        "pdd" : prestoc.PSRPARAMS_pdd_get,
        "f" : prestoc.PSRPARAMS_f_get,
        "fd" : prestoc.PSRPARAMS_fd_get,
        "fdd" : prestoc.PSRPARAMS_fdd_get,
        "orb" : lambda x : orbitparamsPtr(prestoc.PSRPARAMS_orb_get(x)),
    })
    def __getattr__(self,name):
        method = PSRPARAMS.__getmethods__.get(name,None)
        if method: return method(self)
        raise AttributeError,name
    def __repr__(self):
        return "<C PSRPARAMS instance at %s>" % (self.this,)
class PSRPARAMSPtr(PSRPARAMS):
    def __init__(self,this):
        self.this = this
        self.thisown = 0
        self.__class__ = PSRPARAMS



class DoubleArray:
    def __init__(self,*args):
        self.this = apply(prestoc.new_DoubleArray,args)
        self.thisown = 1

    def __del__(self,prestoc=prestoc):
        if getattr(self,'thisown',0):
            prestoc.delete_DoubleArray(self)
    def __getitem__(*args):
        val = apply(prestoc.DoubleArray___getitem__,args)
        return val
    def __setitem__(*args):
        val = apply(prestoc.DoubleArray___setitem__,args)
        return val
    __setmethods__ = {}
    for _s in []: __setmethods__.update(_s.__setmethods__)
    __setmethods__.update({
        "dptr" : prestoc.DoubleArray_dptr_set,
    })
    def __setattr__(self,name,value):
        if (name == "this") or (name == "thisown"): self.__dict__[name] = value; return
        method = DoubleArray.__setmethods__.get(name,None)
        if method: return method(self,value)
        self.__dict__[name] = value
    __getmethods__ = {}
    for _s in []: __getmethods__.update(_s.__getmethods__)
    __getmethods__.update({
        "dptr" : prestoc.DoubleArray_dptr_get,
    })
    def __getattr__(self,name):
        method = DoubleArray.__getmethods__.get(name,None)
        if method: return method(self)
        raise AttributeError,name
    def __repr__(self):
        return "<C DoubleArray instance at %s>" % (self.this,)
class DoubleArrayPtr(DoubleArray):
    def __init__(self,this):
        self.this = this
        self.thisown = 0
        self.__class__ = DoubleArray



class INFODATA:
    def __init__(self,*args):
        self.this = apply(prestoc.new_INFODATA,args)
        self.thisown = 1

    def __del__(self,prestoc=prestoc):
        if getattr(self,'thisown',0):
            prestoc.delete_infodata(self)
    __setmethods__ = {}
    for _s in []: __setmethods__.update(_s.__setmethods__)
    __setmethods__.update({
        "name" : prestoc.INFODATA_name_set,
        "object" : prestoc.INFODATA_object_set,
        "ra_h" : prestoc.INFODATA_ra_h_set,
        "ra_m" : prestoc.INFODATA_ra_m_set,
        "ra_s" : prestoc.INFODATA_ra_s_set,
        "dec_d" : prestoc.INFODATA_dec_d_set,
        "dec_m" : prestoc.INFODATA_dec_m_set,
        "dec_s" : prestoc.INFODATA_dec_s_set,
        "telescope" : prestoc.INFODATA_telescope_set,
        "instrument" : prestoc.INFODATA_instrument_set,
        "observer" : prestoc.INFODATA_observer_set,
        "N" : prestoc.INFODATA_N_set,
        "dt" : prestoc.INFODATA_dt_set,
        "numonoff" : prestoc.INFODATA_numonoff_set,
        "fov" : prestoc.INFODATA_fov_set,
        "mjd_i" : prestoc.INFODATA_mjd_i_set,
        "mjd_f" : prestoc.INFODATA_mjd_f_set,
        "bary" : prestoc.INFODATA_bary_set,
        "band" : prestoc.INFODATA_band_set,
        "dm" : prestoc.INFODATA_dm_set,
        "freq" : prestoc.INFODATA_freq_set,
        "freqband" : prestoc.INFODATA_freqband_set,
        "num_chan" : prestoc.INFODATA_num_chan_set,
        "chan_wid" : prestoc.INFODATA_chan_wid_set,
        "filt" : prestoc.INFODATA_filt_set,
        "wavelen" : prestoc.INFODATA_wavelen_set,
        "waveband" : prestoc.INFODATA_waveband_set,
        "energy" : prestoc.INFODATA_energy_set,
        "energyband" : prestoc.INFODATA_energyband_set,
        "analyzer" : prestoc.INFODATA_analyzer_set,
        "notes" : prestoc.INFODATA_notes_set,
    })
    def __setattr__(self,name,value):
        if (name == "this") or (name == "thisown"): self.__dict__[name] = value; return
        method = INFODATA.__setmethods__.get(name,None)
        if method: return method(self,value)
        self.__dict__[name] = value
    __getmethods__ = {}
    for _s in []: __getmethods__.update(_s.__getmethods__)
    __getmethods__.update({
        "name" : prestoc.INFODATA_name_get,
        "object" : prestoc.INFODATA_object_get,
        "ra_h" : prestoc.INFODATA_ra_h_get,
        "ra_m" : prestoc.INFODATA_ra_m_get,
        "ra_s" : prestoc.INFODATA_ra_s_get,
        "dec_d" : prestoc.INFODATA_dec_d_get,
        "dec_m" : prestoc.INFODATA_dec_m_get,
        "dec_s" : prestoc.INFODATA_dec_s_get,
        "telescope" : prestoc.INFODATA_telescope_get,
        "instrument" : prestoc.INFODATA_instrument_get,
        "observer" : prestoc.INFODATA_observer_get,
        "N" : prestoc.INFODATA_N_get,
        "dt" : prestoc.INFODATA_dt_get,
        "numonoff" : prestoc.INFODATA_numonoff_get,
        "onoff" : prestoc.INFODATA_onoff_get,
        "fov" : prestoc.INFODATA_fov_get,
        "mjd_i" : prestoc.INFODATA_mjd_i_get,
        "mjd_f" : prestoc.INFODATA_mjd_f_get,
        "bary" : prestoc.INFODATA_bary_get,
        "band" : prestoc.INFODATA_band_get,
        "dm" : prestoc.INFODATA_dm_get,
        "freq" : prestoc.INFODATA_freq_get,
        "freqband" : prestoc.INFODATA_freqband_get,
        "num_chan" : prestoc.INFODATA_num_chan_get,
        "chan_wid" : prestoc.INFODATA_chan_wid_get,
        "filt" : prestoc.INFODATA_filt_get,
        "wavelen" : prestoc.INFODATA_wavelen_get,
        "waveband" : prestoc.INFODATA_waveband_get,
        "energy" : prestoc.INFODATA_energy_get,
        "energyband" : prestoc.INFODATA_energyband_get,
        "analyzer" : prestoc.INFODATA_analyzer_get,
        "notes" : prestoc.INFODATA_notes_get,
    })
    def __getattr__(self,name):
        method = INFODATA.__getmethods__.get(name,None)
        if method: return method(self)
        raise AttributeError,name
    def __repr__(self):
        return "<C INFODATA instance at %s>" % (self.this,)
class INFODATAPtr(INFODATA):
    def __init__(self,this):
        self.this = this
        self.thisown = 0
        self.__class__ = INFODATA



class MAKEDATA:
    def __init__(self,*args):
        self.this = apply(prestoc.new_MAKEDATA,args)
        self.thisown = 1

    def __del__(self,prestoc=prestoc):
        if getattr(self,'thisown',0):
            prestoc.delete_makedata(self)
    __setmethods__ = {}
    for _s in []: __setmethods__.update(_s.__setmethods__)
    __setmethods__.update({
        "basefilenm" : prestoc.MAKEDATA_basefilenm_set,
        "description" : prestoc.MAKEDATA_description_set,
        "N" : prestoc.MAKEDATA_N_set,
        "next2_to_n" : prestoc.MAKEDATA_next2_to_n_set,
        "dt" : prestoc.MAKEDATA_dt_set,
        "T" : prestoc.MAKEDATA_T_set,
        "ptype" : prestoc.MAKEDATA_ptype_set,
        "pnum" : prestoc.MAKEDATA_pnum_set,
        "fwhm" : prestoc.MAKEDATA_fwhm_set,
        "round" : prestoc.MAKEDATA_round_set,
        "roundnum" : prestoc.MAKEDATA_roundnum_set,
        "f" : prestoc.MAKEDATA_f_set,
        "fd" : prestoc.MAKEDATA_fd_set,
        "fdd" : prestoc.MAKEDATA_fdd_set,
        "p" : prestoc.MAKEDATA_p_set,
        "pd" : prestoc.MAKEDATA_pd_set,
        "pdd" : prestoc.MAKEDATA_pdd_set,
        "r" : prestoc.MAKEDATA_r_set,
        "z" : prestoc.MAKEDATA_z_set,
        "w" : prestoc.MAKEDATA_w_set,
        "amp" : prestoc.MAKEDATA_amp_set,
        "phs" : prestoc.MAKEDATA_phs_set,
        "dc" : prestoc.MAKEDATA_dc_set,
        "binary" : prestoc.MAKEDATA_binary_set,
        "orb" : prestoc.MAKEDATA_orb_set,
        "ampmod" : prestoc.MAKEDATA_ampmod_set,
        "ampmoda" : prestoc.MAKEDATA_ampmoda_set,
        "ampmodf" : prestoc.MAKEDATA_ampmodf_set,
        "ampmodp" : prestoc.MAKEDATA_ampmodp_set,
        "noisetype" : prestoc.MAKEDATA_noisetype_set,
        "noise" : prestoc.MAKEDATA_noise_set,
        "noisesig" : prestoc.MAKEDATA_noisesig_set,
        "numonoff" : prestoc.MAKEDATA_numonoff_set,
        "onoff" : prestoc.MAKEDATA_onoff_set,
    })
    def __setattr__(self,name,value):
        if (name == "this") or (name == "thisown"): self.__dict__[name] = value; return
        method = MAKEDATA.__setmethods__.get(name,None)
        if method: return method(self,value)
        self.__dict__[name] = value
    __getmethods__ = {}
    for _s in []: __getmethods__.update(_s.__getmethods__)
    __getmethods__.update({
        "basefilenm" : prestoc.MAKEDATA_basefilenm_get,
        "description" : prestoc.MAKEDATA_description_get,
        "N" : prestoc.MAKEDATA_N_get,
        "next2_to_n" : prestoc.MAKEDATA_next2_to_n_get,
        "dt" : prestoc.MAKEDATA_dt_get,
        "T" : prestoc.MAKEDATA_T_get,
        "ptype" : prestoc.MAKEDATA_ptype_get,
        "pnum" : prestoc.MAKEDATA_pnum_get,
        "fwhm" : prestoc.MAKEDATA_fwhm_get,
        "round" : prestoc.MAKEDATA_round_get,
        "roundnum" : prestoc.MAKEDATA_roundnum_get,
        "f" : prestoc.MAKEDATA_f_get,
        "fd" : prestoc.MAKEDATA_fd_get,
        "fdd" : prestoc.MAKEDATA_fdd_get,
        "p" : prestoc.MAKEDATA_p_get,
        "pd" : prestoc.MAKEDATA_pd_get,
        "pdd" : prestoc.MAKEDATA_pdd_get,
        "r" : prestoc.MAKEDATA_r_get,
        "z" : prestoc.MAKEDATA_z_get,
        "w" : prestoc.MAKEDATA_w_get,
        "amp" : prestoc.MAKEDATA_amp_get,
        "phs" : prestoc.MAKEDATA_phs_get,
        "dc" : prestoc.MAKEDATA_dc_get,
        "binary" : prestoc.MAKEDATA_binary_get,
        "orb" : lambda x : orbitparamsPtr(prestoc.MAKEDATA_orb_get(x)),
        "ampmod" : prestoc.MAKEDATA_ampmod_get,
        "ampmoda" : prestoc.MAKEDATA_ampmoda_get,
        "ampmodf" : prestoc.MAKEDATA_ampmodf_get,
        "ampmodp" : prestoc.MAKEDATA_ampmodp_get,
        "noisetype" : prestoc.MAKEDATA_noisetype_get,
        "noise" : prestoc.MAKEDATA_noise_get,
        "noisesig" : prestoc.MAKEDATA_noisesig_get,
        "numonoff" : prestoc.MAKEDATA_numonoff_get,
        "onoff" : lambda x : DoubleArrayPtr(prestoc.MAKEDATA_onoff_get(x)),
    })
    def __getattr__(self,name):
        method = MAKEDATA.__getmethods__.get(name,None)
        if method: return method(self)
        raise AttributeError,name
    def __repr__(self):
        return "<C MAKEDATA instance at %s>" % (self.this,)
class MAKEDATAPtr(MAKEDATA):
    def __init__(self,this):
        self.this = this
        self.thisown = 0
        self.__class__ = MAKEDATA



class RDERIVS:
    def __init__(self,*args):
        self.this = apply(prestoc.new_RDERIVS,args)
        self.thisown = 1

    def __del__(self,prestoc=prestoc):
        if getattr(self,'thisown',0):
            prestoc.delete_rderivs(self)
    __setmethods__ = {}
    for _s in []: __setmethods__.update(_s.__setmethods__)
    __setmethods__.update({
        "pow" : prestoc.RDERIVS_pow_set,
        "phs" : prestoc.RDERIVS_phs_set,
        "dpow" : prestoc.RDERIVS_dpow_set,
        "dphs" : prestoc.RDERIVS_dphs_set,
        "d2pow" : prestoc.RDERIVS_d2pow_set,
        "d2phs" : prestoc.RDERIVS_d2phs_set,
        "locpow" : prestoc.RDERIVS_locpow_set,
    })
    def __setattr__(self,name,value):
        if (name == "this") or (name == "thisown"): self.__dict__[name] = value; return
        method = RDERIVS.__setmethods__.get(name,None)
        if method: return method(self,value)
        self.__dict__[name] = value
    __getmethods__ = {}
    for _s in []: __getmethods__.update(_s.__getmethods__)
    __getmethods__.update({
        "pow" : prestoc.RDERIVS_pow_get,
        "phs" : prestoc.RDERIVS_phs_get,
        "dpow" : prestoc.RDERIVS_dpow_get,
        "dphs" : prestoc.RDERIVS_dphs_get,
        "d2pow" : prestoc.RDERIVS_d2pow_get,
        "d2phs" : prestoc.RDERIVS_d2phs_get,
        "locpow" : prestoc.RDERIVS_locpow_get,
    })
    def __getattr__(self,name):
        method = RDERIVS.__getmethods__.get(name,None)
        if method: return method(self)
        raise AttributeError,name
    def __repr__(self):
        return "<C RDERIVS instance at %s>" % (self.this,)
class RDERIVSPtr(RDERIVS):
    def __init__(self,this):
        self.this = this
        self.thisown = 0
        self.__class__ = RDERIVS



class FOURIERPROPS:
    def __init__(self,*args):
        self.this = apply(prestoc.new_FOURIERPROPS,args)
        self.thisown = 1

    def __del__(self,prestoc=prestoc):
        if getattr(self,'thisown',0):
            prestoc.delete_fourierprops(self)
    __setmethods__ = {}
    for _s in []: __setmethods__.update(_s.__setmethods__)
    __setmethods__.update({
        "r" : prestoc.FOURIERPROPS_r_set,
        "rerr" : prestoc.FOURIERPROPS_rerr_set,
        "z" : prestoc.FOURIERPROPS_z_set,
        "zerr" : prestoc.FOURIERPROPS_zerr_set,
        "w" : prestoc.FOURIERPROPS_w_set,
        "werr" : prestoc.FOURIERPROPS_werr_set,
        "pow" : prestoc.FOURIERPROPS_pow_set,
        "powerr" : prestoc.FOURIERPROPS_powerr_set,
        "sig" : prestoc.FOURIERPROPS_sig_set,
        "rawpow" : prestoc.FOURIERPROPS_rawpow_set,
        "phs" : prestoc.FOURIERPROPS_phs_set,
        "phserr" : prestoc.FOURIERPROPS_phserr_set,
        "cen" : prestoc.FOURIERPROPS_cen_set,
        "cenerr" : prestoc.FOURIERPROPS_cenerr_set,
        "pur" : prestoc.FOURIERPROPS_pur_set,
        "purerr" : prestoc.FOURIERPROPS_purerr_set,
        "locpow" : prestoc.FOURIERPROPS_locpow_set,
    })
    def __setattr__(self,name,value):
        if (name == "this") or (name == "thisown"): self.__dict__[name] = value; return
        method = FOURIERPROPS.__setmethods__.get(name,None)
        if method: return method(self,value)
        self.__dict__[name] = value
    __getmethods__ = {}
    for _s in []: __getmethods__.update(_s.__getmethods__)
    __getmethods__.update({
        "r" : prestoc.FOURIERPROPS_r_get,
        "rerr" : prestoc.FOURIERPROPS_rerr_get,
        "z" : prestoc.FOURIERPROPS_z_get,
        "zerr" : prestoc.FOURIERPROPS_zerr_get,
        "w" : prestoc.FOURIERPROPS_w_get,
        "werr" : prestoc.FOURIERPROPS_werr_get,
        "pow" : prestoc.FOURIERPROPS_pow_get,
        "powerr" : prestoc.FOURIERPROPS_powerr_get,
        "sig" : prestoc.FOURIERPROPS_sig_get,
        "rawpow" : prestoc.FOURIERPROPS_rawpow_get,
        "phs" : prestoc.FOURIERPROPS_phs_get,
        "phserr" : prestoc.FOURIERPROPS_phserr_get,
        "cen" : prestoc.FOURIERPROPS_cen_get,
        "cenerr" : prestoc.FOURIERPROPS_cenerr_get,
        "pur" : prestoc.FOURIERPROPS_pur_get,
        "purerr" : prestoc.FOURIERPROPS_purerr_get,
        "locpow" : prestoc.FOURIERPROPS_locpow_get,
    })
    def __getattr__(self,name):
        method = FOURIERPROPS.__getmethods__.get(name,None)
        if method: return method(self)
        raise AttributeError,name
    def __repr__(self):
        return "<C FOURIERPROPS instance at %s>" % (self.this,)
class FOURIERPROPSPtr(FOURIERPROPS):
    def __init__(self,this):
        self.this = this
        self.thisown = 0
        self.__class__ = FOURIERPROPS



class BINARYPROPS:
    def __init__(self,*args):
        self.this = apply(prestoc.new_BINARYPROPS,args)
        self.thisown = 1

    def __del__(self,prestoc=prestoc):
        if getattr(self,'thisown',0):
            prestoc.delete_binaryprops(self)
    __setmethods__ = {}
    for _s in []: __setmethods__.update(_s.__setmethods__)
    __setmethods__.update({
        "ppsr" : prestoc.BINARYPROPS_ppsr_set,
        "fpsr" : prestoc.BINARYPROPS_fpsr_set,
        "rpsr" : prestoc.BINARYPROPS_rpsr_set,
        "pbin" : prestoc.BINARYPROPS_pbin_set,
        "rbin" : prestoc.BINARYPROPS_rbin_set,
        "z" : prestoc.BINARYPROPS_z_set,
        "asinic" : prestoc.BINARYPROPS_asinic_set,
        "rdetect" : prestoc.BINARYPROPS_rdetect_set,
        "nfftbins" : prestoc.BINARYPROPS_nfftbins_set,
        "lowbin" : prestoc.BINARYPROPS_lowbin_set,
        "ppsrerr" : prestoc.BINARYPROPS_ppsrerr_set,
        "fpsrerr" : prestoc.BINARYPROPS_fpsrerr_set,
        "rpsrerr" : prestoc.BINARYPROPS_rpsrerr_set,
        "pbinerr" : prestoc.BINARYPROPS_pbinerr_set,
        "rbinerr" : prestoc.BINARYPROPS_rbinerr_set,
        "zerr" : prestoc.BINARYPROPS_zerr_set,
        "asinicerr" : prestoc.BINARYPROPS_asinicerr_set,
        "rdetecterr" : prestoc.BINARYPROPS_rdetecterr_set,
        "sig" : prestoc.BINARYPROPS_sig_set,
        "phs" : prestoc.BINARYPROPS_phs_set,
        "phserr" : prestoc.BINARYPROPS_phserr_set,
        "cen" : prestoc.BINARYPROPS_cen_set,
        "cenerr" : prestoc.BINARYPROPS_cenerr_set,
        "pur" : prestoc.BINARYPROPS_pur_set,
        "purerr" : prestoc.BINARYPROPS_purerr_set,
        "pow" : prestoc.BINARYPROPS_pow_set,
        "powerr" : prestoc.BINARYPROPS_powerr_set,
    })
    def __setattr__(self,name,value):
        if (name == "this") or (name == "thisown"): self.__dict__[name] = value; return
        method = BINARYPROPS.__setmethods__.get(name,None)
        if method: return method(self,value)
        self.__dict__[name] = value
    __getmethods__ = {}
    for _s in []: __getmethods__.update(_s.__getmethods__)
    __getmethods__.update({
        "ppsr" : prestoc.BINARYPROPS_ppsr_get,
        "fpsr" : prestoc.BINARYPROPS_fpsr_get,
        "rpsr" : prestoc.BINARYPROPS_rpsr_get,
        "pbin" : prestoc.BINARYPROPS_pbin_get,
        "rbin" : prestoc.BINARYPROPS_rbin_get,
        "z" : prestoc.BINARYPROPS_z_get,
        "asinic" : prestoc.BINARYPROPS_asinic_get,
        "rdetect" : prestoc.BINARYPROPS_rdetect_get,
        "nfftbins" : prestoc.BINARYPROPS_nfftbins_get,
        "lowbin" : prestoc.BINARYPROPS_lowbin_get,
        "ppsrerr" : prestoc.BINARYPROPS_ppsrerr_get,
        "fpsrerr" : prestoc.BINARYPROPS_fpsrerr_get,
        "rpsrerr" : prestoc.BINARYPROPS_rpsrerr_get,
        "pbinerr" : prestoc.BINARYPROPS_pbinerr_get,
        "rbinerr" : prestoc.BINARYPROPS_rbinerr_get,
        "zerr" : prestoc.BINARYPROPS_zerr_get,
        "asinicerr" : prestoc.BINARYPROPS_asinicerr_get,
        "rdetecterr" : prestoc.BINARYPROPS_rdetecterr_get,
        "sig" : prestoc.BINARYPROPS_sig_get,
        "phs" : prestoc.BINARYPROPS_phs_get,
        "phserr" : prestoc.BINARYPROPS_phserr_get,
        "cen" : prestoc.BINARYPROPS_cen_get,
        "cenerr" : prestoc.BINARYPROPS_cenerr_get,
        "pur" : prestoc.BINARYPROPS_pur_get,
        "purerr" : prestoc.BINARYPROPS_purerr_get,
        "pow" : prestoc.BINARYPROPS_pow_get,
        "powerr" : prestoc.BINARYPROPS_powerr_get,
    })
    def __getattr__(self,name):
        method = BINARYPROPS.__getmethods__.get(name,None)
        if method: return method(self)
        raise AttributeError,name
    def __repr__(self):
        return "<C BINARYPROPS instance at %s>" % (self.this,)
class BINARYPROPSPtr(BINARYPROPS):
    def __init__(self,this):
        self.this = this
        self.thisown = 0
        self.__class__ = BINARYPROPS



class RAWBINCAND:
    def __init__(self,*args):
        self.this = apply(prestoc.new_RAWBINCAND,args)
        self.thisown = 1

    def __del__(self,prestoc=prestoc):
        if getattr(self,'thisown',0):
            prestoc.delete_rawbincand(self)
    __setmethods__ = {}
    for _s in []: __setmethods__.update(_s.__setmethods__)
    __setmethods__.update({
        "full_N" : prestoc.RAWBINCAND_full_N_set,
        "full_T" : prestoc.RAWBINCAND_full_T_set,
        "full_lo_r" : prestoc.RAWBINCAND_full_lo_r_set,
        "mini_N" : prestoc.RAWBINCAND_mini_N_set,
        "mini_r" : prestoc.RAWBINCAND_mini_r_set,
        "mini_power" : prestoc.RAWBINCAND_mini_power_set,
        "mini_numsum" : prestoc.RAWBINCAND_mini_numsum_set,
        "mini_sigma" : prestoc.RAWBINCAND_mini_sigma_set,
        "psr_p" : prestoc.RAWBINCAND_psr_p_set,
        "orb_p" : prestoc.RAWBINCAND_orb_p_set,
    })
    def __setattr__(self,name,value):
        if (name == "this") or (name == "thisown"): self.__dict__[name] = value; return
        method = RAWBINCAND.__setmethods__.get(name,None)
        if method: return method(self,value)
        self.__dict__[name] = value
    __getmethods__ = {}
    for _s in []: __getmethods__.update(_s.__getmethods__)
    __getmethods__.update({
        "full_N" : prestoc.RAWBINCAND_full_N_get,
        "full_T" : prestoc.RAWBINCAND_full_T_get,
        "full_lo_r" : prestoc.RAWBINCAND_full_lo_r_get,
        "mini_N" : prestoc.RAWBINCAND_mini_N_get,
        "mini_r" : prestoc.RAWBINCAND_mini_r_get,
        "mini_power" : prestoc.RAWBINCAND_mini_power_get,
        "mini_numsum" : prestoc.RAWBINCAND_mini_numsum_get,
        "mini_sigma" : prestoc.RAWBINCAND_mini_sigma_get,
        "psr_p" : prestoc.RAWBINCAND_psr_p_get,
        "orb_p" : prestoc.RAWBINCAND_orb_p_get,
    })
    def __getattr__(self,name):
        method = RAWBINCAND.__getmethods__.get(name,None)
        if method: return method(self)
        raise AttributeError,name
    def __repr__(self):
        return "<C RAWBINCAND instance at %s>" % (self.this,)
class RAWBINCANDPtr(RAWBINCAND):
    def __init__(self,this):
        self.this = this
        self.thisown = 0
        self.__class__ = RAWBINCAND



class foldstats:
    def __init__(self,*args):
        self.this = apply(prestoc.new_foldstats,args)
        self.thisown = 1

    def __del__(self,prestoc=prestoc):
        if getattr(self,'thisown',0):
            prestoc.delete_foldstats(self)
    __setmethods__ = {}
    for _s in []: __setmethods__.update(_s.__setmethods__)
    __setmethods__.update({
        "numdata" : prestoc.foldstats_numdata_set,
        "data_avg" : prestoc.foldstats_data_avg_set,
        "data_var" : prestoc.foldstats_data_var_set,
        "numprof" : prestoc.foldstats_numprof_set,
        "prof_avg" : prestoc.foldstats_prof_avg_set,
        "prof_var" : prestoc.foldstats_prof_var_set,
        "redchi" : prestoc.foldstats_redchi_set,
    })
    def __setattr__(self,name,value):
        if (name == "this") or (name == "thisown"): self.__dict__[name] = value; return
        method = foldstats.__setmethods__.get(name,None)
        if method: return method(self,value)
        self.__dict__[name] = value
    __getmethods__ = {}
    for _s in []: __getmethods__.update(_s.__getmethods__)
    __getmethods__.update({
        "numdata" : prestoc.foldstats_numdata_get,
        "data_avg" : prestoc.foldstats_data_avg_get,
        "data_var" : prestoc.foldstats_data_var_get,
        "numprof" : prestoc.foldstats_numprof_get,
        "prof_avg" : prestoc.foldstats_prof_avg_get,
        "prof_var" : prestoc.foldstats_prof_var_get,
        "redchi" : prestoc.foldstats_redchi_get,
    })
    def __getattr__(self,name):
        method = foldstats.__getmethods__.get(name,None)
        if method: return method(self)
        raise AttributeError,name
    def __repr__(self):
        return "<C foldstats instance at %s>" % (self.this,)
class foldstatsPtr(foldstats):
    def __init__(self,this):
        self.this = this
        self.thisown = 0
        self.__class__ = foldstats





#-------------- FUNCTION WRAPPERS ------------------

tofloatvector = prestoc.tofloatvector

float_to_complex = prestoc.float_to_complex

complex_to_float = prestoc.complex_to_float

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



#-------------- VARIABLE WRAPPERS ------------------

SQRT2 = prestoc.SQRT2
PI = prestoc.PI
TWOPI = prestoc.TWOPI
DEGTORAD = prestoc.DEGTORAD
RADTODEG = prestoc.RADTODEG
PIBYTWO = prestoc.PIBYTWO
SOL = prestoc.SOL
SECPERJULYR = prestoc.SECPERJULYR
SECPERDAY = prestoc.SECPERDAY
ARCSEC2RAD = prestoc.ARCSEC2RAD
SEC2RAD = prestoc.SEC2RAD
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
