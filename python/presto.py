# This file was created automatically by SWIG.
import prestoc
import new
class orbitparams:
    def __init__(self,*args):
        self.this = apply(prestoc.new_orbitparams,args)
        self.thisown = 1

    def __del__(self,prestoc=prestoc):
        if self.thisown == 1 :
            prestoc.delete_orbitparams(self)
    __setmethods__ = {
        "p" : prestoc.orbitparams_p_set,
        "e" : prestoc.orbitparams_e_set,
        "x" : prestoc.orbitparams_x_set,
        "w" : prestoc.orbitparams_w_set,
        "t" : prestoc.orbitparams_t_set,
        "pd" : prestoc.orbitparams_pd_set,
        "wd" : prestoc.orbitparams_wd_set,
    }
    def __setattr__(self,name,value):
        if (name == "this") or (name == "thisown"): self.__dict__[name] = value; return
        method = orbitparams.__setmethods__.get(name,None)
        if method: return method(self,value)
        self.__dict__[name] = value
    __getmethods__ = {
        "p" : prestoc.orbitparams_p_get,
        "e" : prestoc.orbitparams_e_get,
        "x" : prestoc.orbitparams_x_get,
        "w" : prestoc.orbitparams_w_get,
        "t" : prestoc.orbitparams_t_get,
        "pd" : prestoc.orbitparams_pd_get,
        "wd" : prestoc.orbitparams_wd_get,
    }
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



class psrdata:
    def __init__(self,*args):
        self.this = apply(prestoc.new_psrdata,args)
        self.thisown = 1

    def __del__(self,prestoc=prestoc):
        if self.thisown == 1 :
            prestoc.delete_psrdata(self)
    __setmethods__ = {
        "ra2000" : prestoc.psrdata_ra2000_set,
        "ra1950" : prestoc.psrdata_ra1950_set,
        "rae" : prestoc.psrdata_rae_set,
        "dec2000" : prestoc.psrdata_dec2000_set,
        "dec1950" : prestoc.psrdata_dec1950_set,
        "dece" : prestoc.psrdata_dece_set,
        "dmin__" : prestoc.psrdata_dmin___set,
        "dmax__" : prestoc.psrdata_dmax___set,
        "dist" : prestoc.psrdata_dist_set,
        "ldeg" : prestoc.psrdata_ldeg_set,
        "bdeg" : prestoc.psrdata_bdeg_set,
        "pmra" : prestoc.psrdata_pmra_set,
        "pmrae" : prestoc.psrdata_pmrae_set,
        "pmdec" : prestoc.psrdata_pmdec_set,
        "pmdece" : prestoc.psrdata_pmdece_set,
        "posepoch" : prestoc.psrdata_posepoch_set,
        "p" : prestoc.psrdata_p_set,
        "pe" : prestoc.psrdata_pe_set,
        "pdot" : prestoc.psrdata_pdot_set,
        "pdote" : prestoc.psrdata_pdote_set,
        "f2" : prestoc.psrdata_f2_set,
        "f2e" : prestoc.psrdata_f2e_set,
        "f3" : prestoc.psrdata_f3_set,
        "f3e" : prestoc.psrdata_f3e_set,
        "epoch" : prestoc.psrdata_epoch_set,
        "dm" : prestoc.psrdata_dm_set,
        "dme" : prestoc.psrdata_dme_set,
        "rm" : prestoc.psrdata_rm_set,
        "rme" : prestoc.psrdata_rme_set,
        "we" : prestoc.psrdata_we_set,
        "w50" : prestoc.psrdata_w50_set,
        "w10" : prestoc.psrdata_w10_set,
        "s400" : prestoc.psrdata_s400_set,
        "s600" : prestoc.psrdata_s600_set,
        "s1400" : prestoc.psrdata_s1400_set,
        "tau" : prestoc.psrdata_tau_set,
        "t408" : prestoc.psrdata_t408_set,
        "distmod" : prestoc.psrdata_distmod_set,
        "lum" : prestoc.psrdata_lum_set,
        "bsurf" : prestoc.psrdata_bsurf_set,
        "age" : prestoc.psrdata_age_set,
        "edot" : prestoc.psrdata_edot_set,
        "pb" : prestoc.psrdata_pb_set,
        "pbe" : prestoc.psrdata_pbe_set,
        "a1" : prestoc.psrdata_a1_set,
        "a1e" : prestoc.psrdata_a1e_set,
        "om" : prestoc.psrdata_om_set,
        "ome" : prestoc.psrdata_ome_set,
        "omdot" : prestoc.psrdata_omdot_set,
        "omdote" : prestoc.psrdata_omdote_set,
        "e" : prestoc.psrdata_e_set,
        "ee" : prestoc.psrdata_ee_set,
        "t0" : prestoc.psrdata_t0_set,
        "t0e" : prestoc.psrdata_t0e_set,
        "gamma" : prestoc.psrdata_gamma_set,
        "gammae" : prestoc.psrdata_gammae_set,
        "pbdot" : prestoc.psrdata_pbdot_set,
        "pbdote" : prestoc.psrdata_pbdote_set,
        "si" : prestoc.psrdata_si_set,
        "sie" : prestoc.psrdata_sie_set,
        "r__" : prestoc.psrdata_r___set,
        "re" : prestoc.psrdata_re_set,
        "pb2" : prestoc.psrdata_pb2_set,
        "pb2e" : prestoc.psrdata_pb2e_set,
        "a12" : prestoc.psrdata_a12_set,
        "a12e" : prestoc.psrdata_a12e_set,
        "om2" : prestoc.psrdata_om2_set,
        "om2e" : prestoc.psrdata_om2e_set,
        "omdot2" : prestoc.psrdata_omdot2_set,
        "omdot2e" : prestoc.psrdata_omdot2e_set,
        "e2" : prestoc.psrdata_e2_set,
        "e2e" : prestoc.psrdata_e2e_set,
        "t02" : prestoc.psrdata_t02_set,
        "t02e" : prestoc.psrdata_t02e_set,
        "gamma2" : prestoc.psrdata_gamma2_set,
        "gamma2e" : prestoc.psrdata_gamma2e_set,
        "pbdot2" : prestoc.psrdata_pbdot2_set,
        "pbdot2e" : prestoc.psrdata_pbdot2e_set,
        "si2" : prestoc.psrdata_si2_set,
        "si2e" : prestoc.psrdata_si2e_set,
        "r2" : prestoc.psrdata_r2_set,
        "r2e" : prestoc.psrdata_r2e_set,
        "nscode" : prestoc.psrdata_nscode_set,
        "ndflag" : prestoc.psrdata_ndflag_set,
        "ntauflag" : prestoc.psrdata_ntauflag_set,
        "ntype" : prestoc.psrdata_ntype_set,
        "modcode" : prestoc.psrdata_modcode_set,
        "limcode" : prestoc.psrdata_limcode_set,
        "ibin" : prestoc.psrdata_ibin_set,
        "jname" : prestoc.psrdata_jname_set,
        "bname" : prestoc.psrdata_bname_set,
        "lcode" : prestoc.psrdata_lcode_set,
        "ucode" : prestoc.psrdata_ucode_set,
    }
    def __setattr__(self,name,value):
        if (name == "this") or (name == "thisown"): self.__dict__[name] = value; return
        method = psrdata.__setmethods__.get(name,None)
        if method: return method(self,value)
        self.__dict__[name] = value
    __getmethods__ = {
        "ra2000" : prestoc.psrdata_ra2000_get,
        "ra1950" : prestoc.psrdata_ra1950_get,
        "rae" : prestoc.psrdata_rae_get,
        "dec2000" : prestoc.psrdata_dec2000_get,
        "dec1950" : prestoc.psrdata_dec1950_get,
        "dece" : prestoc.psrdata_dece_get,
        "dmin__" : prestoc.psrdata_dmin___get,
        "dmax__" : prestoc.psrdata_dmax___get,
        "dist" : prestoc.psrdata_dist_get,
        "ldeg" : prestoc.psrdata_ldeg_get,
        "bdeg" : prestoc.psrdata_bdeg_get,
        "pmra" : prestoc.psrdata_pmra_get,
        "pmrae" : prestoc.psrdata_pmrae_get,
        "pmdec" : prestoc.psrdata_pmdec_get,
        "pmdece" : prestoc.psrdata_pmdece_get,
        "posepoch" : prestoc.psrdata_posepoch_get,
        "p" : prestoc.psrdata_p_get,
        "pe" : prestoc.psrdata_pe_get,
        "pdot" : prestoc.psrdata_pdot_get,
        "pdote" : prestoc.psrdata_pdote_get,
        "f2" : prestoc.psrdata_f2_get,
        "f2e" : prestoc.psrdata_f2e_get,
        "f3" : prestoc.psrdata_f3_get,
        "f3e" : prestoc.psrdata_f3e_get,
        "epoch" : prestoc.psrdata_epoch_get,
        "dm" : prestoc.psrdata_dm_get,
        "dme" : prestoc.psrdata_dme_get,
        "rm" : prestoc.psrdata_rm_get,
        "rme" : prestoc.psrdata_rme_get,
        "we" : prestoc.psrdata_we_get,
        "w50" : prestoc.psrdata_w50_get,
        "w10" : prestoc.psrdata_w10_get,
        "s400" : prestoc.psrdata_s400_get,
        "s600" : prestoc.psrdata_s600_get,
        "s1400" : prestoc.psrdata_s1400_get,
        "tau" : prestoc.psrdata_tau_get,
        "t408" : prestoc.psrdata_t408_get,
        "distmod" : prestoc.psrdata_distmod_get,
        "lum" : prestoc.psrdata_lum_get,
        "bsurf" : prestoc.psrdata_bsurf_get,
        "age" : prestoc.psrdata_age_get,
        "edot" : prestoc.psrdata_edot_get,
        "pb" : prestoc.psrdata_pb_get,
        "pbe" : prestoc.psrdata_pbe_get,
        "a1" : prestoc.psrdata_a1_get,
        "a1e" : prestoc.psrdata_a1e_get,
        "om" : prestoc.psrdata_om_get,
        "ome" : prestoc.psrdata_ome_get,
        "omdot" : prestoc.psrdata_omdot_get,
        "omdote" : prestoc.psrdata_omdote_get,
        "e" : prestoc.psrdata_e_get,
        "ee" : prestoc.psrdata_ee_get,
        "t0" : prestoc.psrdata_t0_get,
        "t0e" : prestoc.psrdata_t0e_get,
        "gamma" : prestoc.psrdata_gamma_get,
        "gammae" : prestoc.psrdata_gammae_get,
        "pbdot" : prestoc.psrdata_pbdot_get,
        "pbdote" : prestoc.psrdata_pbdote_get,
        "si" : prestoc.psrdata_si_get,
        "sie" : prestoc.psrdata_sie_get,
        "r__" : prestoc.psrdata_r___get,
        "re" : prestoc.psrdata_re_get,
        "pb2" : prestoc.psrdata_pb2_get,
        "pb2e" : prestoc.psrdata_pb2e_get,
        "a12" : prestoc.psrdata_a12_get,
        "a12e" : prestoc.psrdata_a12e_get,
        "om2" : prestoc.psrdata_om2_get,
        "om2e" : prestoc.psrdata_om2e_get,
        "omdot2" : prestoc.psrdata_omdot2_get,
        "omdot2e" : prestoc.psrdata_omdot2e_get,
        "e2" : prestoc.psrdata_e2_get,
        "e2e" : prestoc.psrdata_e2e_get,
        "t02" : prestoc.psrdata_t02_get,
        "t02e" : prestoc.psrdata_t02e_get,
        "gamma2" : prestoc.psrdata_gamma2_get,
        "gamma2e" : prestoc.psrdata_gamma2e_get,
        "pbdot2" : prestoc.psrdata_pbdot2_get,
        "pbdot2e" : prestoc.psrdata_pbdot2e_get,
        "si2" : prestoc.psrdata_si2_get,
        "si2e" : prestoc.psrdata_si2e_get,
        "r2" : prestoc.psrdata_r2_get,
        "r2e" : prestoc.psrdata_r2e_get,
        "nscode" : prestoc.psrdata_nscode_get,
        "ndflag" : prestoc.psrdata_ndflag_get,
        "ntauflag" : prestoc.psrdata_ntauflag_get,
        "ntype" : prestoc.psrdata_ntype_get,
        "modcode" : prestoc.psrdata_modcode_get,
        "limcode" : prestoc.psrdata_limcode_get,
        "ibin" : prestoc.psrdata_ibin_get,
        "jname" : prestoc.psrdata_jname_get,
        "bname" : prestoc.psrdata_bname_get,
        "lcode" : prestoc.psrdata_lcode_get,
        "ucode" : prestoc.psrdata_ucode_get,
    }
    def __getattr__(self,name):
        method = psrdata.__getmethods__.get(name,None)
        if method: return method(self)
        raise AttributeError,name
    def __repr__(self):
        return "<C psrdata instance at %s>" % (self.this,)
class psrdataPtr(psrdata):
    def __init__(self,this):
        self.this = this
        self.thisown = 0
        self.__class__ = psrdata



class psrparams:
    def __init__(self,*args):
        self.this = apply(prestoc.new_psrparams,args)
        self.thisown = 1

    def __del__(self,prestoc=prestoc):
        if self.thisown == 1 :
            prestoc.delete_psrparams(self)
    __setmethods__ = {
        "jname" : prestoc.psrparams_jname_set,
        "bname" : prestoc.psrparams_bname_set,
        "ntype" : prestoc.psrparams_ntype_set,
        "ra2000" : prestoc.psrparams_ra2000_set,
        "dec2000" : prestoc.psrparams_dec2000_set,
        "dm" : prestoc.psrparams_dm_set,
        "dist" : prestoc.psrparams_dist_set,
        "fwhm" : prestoc.psrparams_fwhm_set,
        "timepoch" : prestoc.psrparams_timepoch_set,
        "p" : prestoc.psrparams_p_set,
        "pd" : prestoc.psrparams_pd_set,
        "pdd" : prestoc.psrparams_pdd_set,
        "f" : prestoc.psrparams_f_set,
        "fd" : prestoc.psrparams_fd_set,
        "fdd" : prestoc.psrparams_fdd_set,
        "orb" : prestoc.psrparams_orb_set,
    }
    def __setattr__(self,name,value):
        if (name == "this") or (name == "thisown"): self.__dict__[name] = value; return
        method = psrparams.__setmethods__.get(name,None)
        if method: return method(self,value)
        self.__dict__[name] = value
    __getmethods__ = {
        "jname" : prestoc.psrparams_jname_get,
        "bname" : prestoc.psrparams_bname_get,
        "ntype" : prestoc.psrparams_ntype_get,
        "ra2000" : prestoc.psrparams_ra2000_get,
        "dec2000" : prestoc.psrparams_dec2000_get,
        "dm" : prestoc.psrparams_dm_get,
        "dist" : prestoc.psrparams_dist_get,
        "fwhm" : prestoc.psrparams_fwhm_get,
        "timepoch" : prestoc.psrparams_timepoch_get,
        "p" : prestoc.psrparams_p_get,
        "pd" : prestoc.psrparams_pd_get,
        "pdd" : prestoc.psrparams_pdd_get,
        "f" : prestoc.psrparams_f_get,
        "fd" : prestoc.psrparams_fd_get,
        "fdd" : prestoc.psrparams_fdd_get,
        "orb" : lambda x : orbitparamsPtr(prestoc.psrparams_orb_get(x)),
    }
    def __getattr__(self,name):
        method = psrparams.__getmethods__.get(name,None)
        if method: return method(self)
        raise AttributeError,name
    def __repr__(self):
        return "<C psrparams instance at %s>" % (self.this,)
class psrparamsPtr(psrparams):
    def __init__(self,this):
        self.this = this
        self.thisown = 0
        self.__class__ = psrparams



class DoubleArray:
    def __init__(self,*args):
        self.this = apply(prestoc.new_DoubleArray,args)
        self.thisown = 1

    def __del__(self,prestoc=prestoc):
        if self.thisown == 1 :
            prestoc.delete_DoubleArray(self)
    __setmethods__ = {
        "dptr" : prestoc.DoubleArray_dptr_set,
    }
    def __setattr__(self,name,value):
        if (name == "this") or (name == "thisown"): self.__dict__[name] = value; return
        method = DoubleArray.__setmethods__.get(name,None)
        if method: return method(self,value)
        self.__dict__[name] = value
    __getmethods__ = {
        "dptr" : prestoc.DoubleArray_dptr_get,
    }
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


DoubleArray.__getitem__ = new.instancemethod(prestoc.DoubleArray___getitem__, None, DoubleArray)
DoubleArray.__setitem__ = new.instancemethod(prestoc.DoubleArray___setitem__, None, DoubleArray)

class infodata:
    def __init__(self,*args):
        self.this = apply(prestoc.new_infodata,args)
        self.thisown = 1

    def __del__(self,prestoc=prestoc):
        if self.thisown == 1 :
            prestoc.delete_infodata(self)
    __setmethods__ = {
        "name" : prestoc.infodata_name_set,
        "object" : prestoc.infodata_object_set,
        "ra_h" : prestoc.infodata_ra_h_set,
        "ra_m" : prestoc.infodata_ra_m_set,
        "ra_s" : prestoc.infodata_ra_s_set,
        "dec_d" : prestoc.infodata_dec_d_set,
        "dec_m" : prestoc.infodata_dec_m_set,
        "dec_s" : prestoc.infodata_dec_s_set,
        "telescope" : prestoc.infodata_telescope_set,
        "instrument" : prestoc.infodata_instrument_set,
        "observer" : prestoc.infodata_observer_set,
        "N" : prestoc.infodata_N_set,
        "dt" : prestoc.infodata_dt_set,
        "numonoff" : prestoc.infodata_numonoff_set,
        "onoff" : prestoc.infodata_onoff_set,
        "fov" : prestoc.infodata_fov_set,
        "mjd_i" : prestoc.infodata_mjd_i_set,
        "mjd_f" : prestoc.infodata_mjd_f_set,
        "bary" : prestoc.infodata_bary_set,
        "band" : prestoc.infodata_band_set,
        "dm" : prestoc.infodata_dm_set,
        "freq" : prestoc.infodata_freq_set,
        "freqband" : prestoc.infodata_freqband_set,
        "num_chan" : prestoc.infodata_num_chan_set,
        "chan_wid" : prestoc.infodata_chan_wid_set,
        "filt" : prestoc.infodata_filt_set,
        "wavelen" : prestoc.infodata_wavelen_set,
        "waveband" : prestoc.infodata_waveband_set,
        "energy" : prestoc.infodata_energy_set,
        "energyband" : prestoc.infodata_energyband_set,
        "analyzer" : prestoc.infodata_analyzer_set,
        "notes" : prestoc.infodata_notes_set,
    }
    def __setattr__(self,name,value):
        if (name == "this") or (name == "thisown"): self.__dict__[name] = value; return
        method = infodata.__setmethods__.get(name,None)
        if method: return method(self,value)
        self.__dict__[name] = value
    __getmethods__ = {
        "name" : prestoc.infodata_name_get,
        "object" : prestoc.infodata_object_get,
        "ra_h" : prestoc.infodata_ra_h_get,
        "ra_m" : prestoc.infodata_ra_m_get,
        "ra_s" : prestoc.infodata_ra_s_get,
        "dec_d" : prestoc.infodata_dec_d_get,
        "dec_m" : prestoc.infodata_dec_m_get,
        "dec_s" : prestoc.infodata_dec_s_get,
        "telescope" : prestoc.infodata_telescope_get,
        "instrument" : prestoc.infodata_instrument_get,
        "observer" : prestoc.infodata_observer_get,
        "N" : prestoc.infodata_N_get,
        "dt" : prestoc.infodata_dt_get,
        "numonoff" : prestoc.infodata_numonoff_get,
        "onoff" : lambda x : DoubleArrayPtr(prestoc.infodata_onoff_get(x)),
        "fov" : prestoc.infodata_fov_get,
        "mjd_i" : prestoc.infodata_mjd_i_get,
        "mjd_f" : prestoc.infodata_mjd_f_get,
        "bary" : prestoc.infodata_bary_get,
        "band" : prestoc.infodata_band_get,
        "dm" : prestoc.infodata_dm_get,
        "freq" : prestoc.infodata_freq_get,
        "freqband" : prestoc.infodata_freqband_get,
        "num_chan" : prestoc.infodata_num_chan_get,
        "chan_wid" : prestoc.infodata_chan_wid_get,
        "filt" : prestoc.infodata_filt_get,
        "wavelen" : prestoc.infodata_wavelen_get,
        "waveband" : prestoc.infodata_waveband_get,
        "energy" : prestoc.infodata_energy_get,
        "energyband" : prestoc.infodata_energyband_get,
        "analyzer" : prestoc.infodata_analyzer_get,
        "notes" : prestoc.infodata_notes_get,
    }
    def __getattr__(self,name):
        method = infodata.__getmethods__.get(name,None)
        if method: return method(self)
        raise AttributeError,name
    def __repr__(self):
        return "<C infodata instance at %s>" % (self.this,)
class infodataPtr(infodata):
    def __init__(self,this):
        self.this = this
        self.thisown = 0
        self.__class__ = infodata



class makedata:
    def __init__(self,*args):
        self.this = apply(prestoc.new_makedata,args)
        self.thisown = 1

    def __del__(self,prestoc=prestoc):
        if self.thisown == 1 :
            prestoc.delete_makedata(self)
    __setmethods__ = {
        "basefilenm" : prestoc.makedata_basefilenm_set,
        "description" : prestoc.makedata_description_set,
        "N" : prestoc.makedata_N_set,
        "next2_to_n" : prestoc.makedata_next2_to_n_set,
        "dt" : prestoc.makedata_dt_set,
        "T" : prestoc.makedata_T_set,
        "ptype" : prestoc.makedata_ptype_set,
        "pnum" : prestoc.makedata_pnum_set,
        "fwhm" : prestoc.makedata_fwhm_set,
        "round" : prestoc.makedata_round_set,
        "roundnum" : prestoc.makedata_roundnum_set,
        "f" : prestoc.makedata_f_set,
        "fd" : prestoc.makedata_fd_set,
        "fdd" : prestoc.makedata_fdd_set,
        "p" : prestoc.makedata_p_set,
        "pd" : prestoc.makedata_pd_set,
        "pdd" : prestoc.makedata_pdd_set,
        "r" : prestoc.makedata_r_set,
        "z" : prestoc.makedata_z_set,
        "w" : prestoc.makedata_w_set,
        "amp" : prestoc.makedata_amp_set,
        "phs" : prestoc.makedata_phs_set,
        "dc" : prestoc.makedata_dc_set,
        "binary" : prestoc.makedata_binary_set,
        "orb" : prestoc.makedata_orb_set,
        "ampmod" : prestoc.makedata_ampmod_set,
        "ampmoda" : prestoc.makedata_ampmoda_set,
        "ampmodf" : prestoc.makedata_ampmodf_set,
        "ampmodp" : prestoc.makedata_ampmodp_set,
        "noisetype" : prestoc.makedata_noisetype_set,
        "noise" : prestoc.makedata_noise_set,
        "noisesig" : prestoc.makedata_noisesig_set,
        "numonoff" : prestoc.makedata_numonoff_set,
        "onoff" : prestoc.makedata_onoff_set,
    }
    def __setattr__(self,name,value):
        if (name == "this") or (name == "thisown"): self.__dict__[name] = value; return
        method = makedata.__setmethods__.get(name,None)
        if method: return method(self,value)
        self.__dict__[name] = value
    __getmethods__ = {
        "basefilenm" : prestoc.makedata_basefilenm_get,
        "description" : prestoc.makedata_description_get,
        "N" : prestoc.makedata_N_get,
        "next2_to_n" : prestoc.makedata_next2_to_n_get,
        "dt" : prestoc.makedata_dt_get,
        "T" : prestoc.makedata_T_get,
        "ptype" : prestoc.makedata_ptype_get,
        "pnum" : prestoc.makedata_pnum_get,
        "fwhm" : prestoc.makedata_fwhm_get,
        "round" : prestoc.makedata_round_get,
        "roundnum" : prestoc.makedata_roundnum_get,
        "f" : prestoc.makedata_f_get,
        "fd" : prestoc.makedata_fd_get,
        "fdd" : prestoc.makedata_fdd_get,
        "p" : prestoc.makedata_p_get,
        "pd" : prestoc.makedata_pd_get,
        "pdd" : prestoc.makedata_pdd_get,
        "r" : prestoc.makedata_r_get,
        "z" : prestoc.makedata_z_get,
        "w" : prestoc.makedata_w_get,
        "amp" : prestoc.makedata_amp_get,
        "phs" : prestoc.makedata_phs_get,
        "dc" : prestoc.makedata_dc_get,
        "binary" : prestoc.makedata_binary_get,
        "orb" : lambda x : orbitparamsPtr(prestoc.makedata_orb_get(x)),
        "ampmod" : prestoc.makedata_ampmod_get,
        "ampmoda" : prestoc.makedata_ampmoda_get,
        "ampmodf" : prestoc.makedata_ampmodf_get,
        "ampmodp" : prestoc.makedata_ampmodp_get,
        "noisetype" : prestoc.makedata_noisetype_get,
        "noise" : prestoc.makedata_noise_get,
        "noisesig" : prestoc.makedata_noisesig_get,
        "numonoff" : prestoc.makedata_numonoff_get,
        "onoff" : lambda x : DoubleArrayPtr(prestoc.makedata_onoff_get(x)),
    }
    def __getattr__(self,name):
        method = makedata.__getmethods__.get(name,None)
        if method: return method(self)
        raise AttributeError,name
    def __repr__(self):
        return "<C makedata instance at %s>" % (self.this,)
class makedataPtr(makedata):
    def __init__(self,this):
        self.this = this
        self.thisown = 0
        self.__class__ = makedata



class rderivs:
    def __init__(self,*args):
        self.this = apply(prestoc.new_rderivs,args)
        self.thisown = 1

    def __del__(self,prestoc=prestoc):
        if self.thisown == 1 :
            prestoc.delete_rderivs(self)
    __setmethods__ = {
        "pow" : prestoc.rderivs_pow_set,
        "phs" : prestoc.rderivs_phs_set,
        "dpow" : prestoc.rderivs_dpow_set,
        "dphs" : prestoc.rderivs_dphs_set,
        "d2pow" : prestoc.rderivs_d2pow_set,
        "d2phs" : prestoc.rderivs_d2phs_set,
        "locpow" : prestoc.rderivs_locpow_set,
    }
    def __setattr__(self,name,value):
        if (name == "this") or (name == "thisown"): self.__dict__[name] = value; return
        method = rderivs.__setmethods__.get(name,None)
        if method: return method(self,value)
        self.__dict__[name] = value
    __getmethods__ = {
        "pow" : prestoc.rderivs_pow_get,
        "phs" : prestoc.rderivs_phs_get,
        "dpow" : prestoc.rderivs_dpow_get,
        "dphs" : prestoc.rderivs_dphs_get,
        "d2pow" : prestoc.rderivs_d2pow_get,
        "d2phs" : prestoc.rderivs_d2phs_get,
        "locpow" : prestoc.rderivs_locpow_get,
    }
    def __getattr__(self,name):
        method = rderivs.__getmethods__.get(name,None)
        if method: return method(self)
        raise AttributeError,name
    def __repr__(self):
        return "<C rderivs instance at %s>" % (self.this,)
class rderivsPtr(rderivs):
    def __init__(self,this):
        self.this = this
        self.thisown = 0
        self.__class__ = rderivs



class fourierprops:
    def __init__(self,*args):
        self.this = apply(prestoc.new_fourierprops,args)
        self.thisown = 1

    def __del__(self,prestoc=prestoc):
        if self.thisown == 1 :
            prestoc.delete_fourierprops(self)
    __setmethods__ = {
        "r" : prestoc.fourierprops_r_set,
        "rerr" : prestoc.fourierprops_rerr_set,
        "z" : prestoc.fourierprops_z_set,
        "zerr" : prestoc.fourierprops_zerr_set,
        "w" : prestoc.fourierprops_w_set,
        "werr" : prestoc.fourierprops_werr_set,
        "pow" : prestoc.fourierprops_pow_set,
        "powerr" : prestoc.fourierprops_powerr_set,
        "sig" : prestoc.fourierprops_sig_set,
        "rawpow" : prestoc.fourierprops_rawpow_set,
        "phs" : prestoc.fourierprops_phs_set,
        "phserr" : prestoc.fourierprops_phserr_set,
        "cen" : prestoc.fourierprops_cen_set,
        "cenerr" : prestoc.fourierprops_cenerr_set,
        "pur" : prestoc.fourierprops_pur_set,
        "purerr" : prestoc.fourierprops_purerr_set,
        "locpow" : prestoc.fourierprops_locpow_set,
    }
    def __setattr__(self,name,value):
        if (name == "this") or (name == "thisown"): self.__dict__[name] = value; return
        method = fourierprops.__setmethods__.get(name,None)
        if method: return method(self,value)
        self.__dict__[name] = value
    __getmethods__ = {
        "r" : prestoc.fourierprops_r_get,
        "rerr" : prestoc.fourierprops_rerr_get,
        "z" : prestoc.fourierprops_z_get,
        "zerr" : prestoc.fourierprops_zerr_get,
        "w" : prestoc.fourierprops_w_get,
        "werr" : prestoc.fourierprops_werr_get,
        "pow" : prestoc.fourierprops_pow_get,
        "powerr" : prestoc.fourierprops_powerr_get,
        "sig" : prestoc.fourierprops_sig_get,
        "rawpow" : prestoc.fourierprops_rawpow_get,
        "phs" : prestoc.fourierprops_phs_get,
        "phserr" : prestoc.fourierprops_phserr_get,
        "cen" : prestoc.fourierprops_cen_get,
        "cenerr" : prestoc.fourierprops_cenerr_get,
        "pur" : prestoc.fourierprops_pur_get,
        "purerr" : prestoc.fourierprops_purerr_get,
        "locpow" : prestoc.fourierprops_locpow_get,
    }
    def __getattr__(self,name):
        method = fourierprops.__getmethods__.get(name,None)
        if method: return method(self)
        raise AttributeError,name
    def __repr__(self):
        return "<C fourierprops instance at %s>" % (self.this,)
class fourierpropsPtr(fourierprops):
    def __init__(self,this):
        self.this = this
        self.thisown = 0
        self.__class__ = fourierprops



class binaryprops:
    def __init__(self,*args):
        self.this = apply(prestoc.new_binaryprops,args)
        self.thisown = 1

    def __del__(self,prestoc=prestoc):
        if self.thisown == 1 :
            prestoc.delete_binaryprops(self)
    __setmethods__ = {
        "ppsr" : prestoc.binaryprops_ppsr_set,
        "fpsr" : prestoc.binaryprops_fpsr_set,
        "rpsr" : prestoc.binaryprops_rpsr_set,
        "pbin" : prestoc.binaryprops_pbin_set,
        "rbin" : prestoc.binaryprops_rbin_set,
        "z" : prestoc.binaryprops_z_set,
        "asinic" : prestoc.binaryprops_asinic_set,
        "rdetect" : prestoc.binaryprops_rdetect_set,
        "nfftbins" : prestoc.binaryprops_nfftbins_set,
        "lowbin" : prestoc.binaryprops_lowbin_set,
        "ppsrerr" : prestoc.binaryprops_ppsrerr_set,
        "fpsrerr" : prestoc.binaryprops_fpsrerr_set,
        "rpsrerr" : prestoc.binaryprops_rpsrerr_set,
        "pbinerr" : prestoc.binaryprops_pbinerr_set,
        "rbinerr" : prestoc.binaryprops_rbinerr_set,
        "zerr" : prestoc.binaryprops_zerr_set,
        "asinicerr" : prestoc.binaryprops_asinicerr_set,
        "rdetecterr" : prestoc.binaryprops_rdetecterr_set,
        "sig" : prestoc.binaryprops_sig_set,
        "phs" : prestoc.binaryprops_phs_set,
        "phserr" : prestoc.binaryprops_phserr_set,
        "cen" : prestoc.binaryprops_cen_set,
        "cenerr" : prestoc.binaryprops_cenerr_set,
        "pur" : prestoc.binaryprops_pur_set,
        "purerr" : prestoc.binaryprops_purerr_set,
        "pow" : prestoc.binaryprops_pow_set,
        "powerr" : prestoc.binaryprops_powerr_set,
    }
    def __setattr__(self,name,value):
        if (name == "this") or (name == "thisown"): self.__dict__[name] = value; return
        method = binaryprops.__setmethods__.get(name,None)
        if method: return method(self,value)
        self.__dict__[name] = value
    __getmethods__ = {
        "ppsr" : prestoc.binaryprops_ppsr_get,
        "fpsr" : prestoc.binaryprops_fpsr_get,
        "rpsr" : prestoc.binaryprops_rpsr_get,
        "pbin" : prestoc.binaryprops_pbin_get,
        "rbin" : prestoc.binaryprops_rbin_get,
        "z" : prestoc.binaryprops_z_get,
        "asinic" : prestoc.binaryprops_asinic_get,
        "rdetect" : prestoc.binaryprops_rdetect_get,
        "nfftbins" : prestoc.binaryprops_nfftbins_get,
        "lowbin" : prestoc.binaryprops_lowbin_get,
        "ppsrerr" : prestoc.binaryprops_ppsrerr_get,
        "fpsrerr" : prestoc.binaryprops_fpsrerr_get,
        "rpsrerr" : prestoc.binaryprops_rpsrerr_get,
        "pbinerr" : prestoc.binaryprops_pbinerr_get,
        "rbinerr" : prestoc.binaryprops_rbinerr_get,
        "zerr" : prestoc.binaryprops_zerr_get,
        "asinicerr" : prestoc.binaryprops_asinicerr_get,
        "rdetecterr" : prestoc.binaryprops_rdetecterr_get,
        "sig" : prestoc.binaryprops_sig_get,
        "phs" : prestoc.binaryprops_phs_get,
        "phserr" : prestoc.binaryprops_phserr_get,
        "cen" : prestoc.binaryprops_cen_get,
        "cenerr" : prestoc.binaryprops_cenerr_get,
        "pur" : prestoc.binaryprops_pur_get,
        "purerr" : prestoc.binaryprops_purerr_get,
        "pow" : prestoc.binaryprops_pow_get,
        "powerr" : prestoc.binaryprops_powerr_get,
    }
    def __getattr__(self,name):
        method = binaryprops.__getmethods__.get(name,None)
        if method: return method(self)
        raise AttributeError,name
    def __repr__(self):
        return "<C binaryprops instance at %s>" % (self.this,)
class binarypropsPtr(binaryprops):
    def __init__(self,this):
        self.this = this
        self.thisown = 0
        self.__class__ = binaryprops



class rawbincand:
    def __init__(self,*args):
        self.this = apply(prestoc.new_rawbincand,args)
        self.thisown = 1

    def __del__(self,prestoc=prestoc):
        if self.thisown == 1 :
            prestoc.delete_rawbincand(self)
    __setmethods__ = {
        "full_N" : prestoc.rawbincand_full_N_set,
        "full_T" : prestoc.rawbincand_full_T_set,
        "full_lo_r" : prestoc.rawbincand_full_lo_r_set,
        "mini_N" : prestoc.rawbincand_mini_N_set,
        "mini_r" : prestoc.rawbincand_mini_r_set,
        "mini_power" : prestoc.rawbincand_mini_power_set,
        "mini_numsum" : prestoc.rawbincand_mini_numsum_set,
        "mini_sigma" : prestoc.rawbincand_mini_sigma_set,
        "psr_p" : prestoc.rawbincand_psr_p_set,
        "orb_p" : prestoc.rawbincand_orb_p_set,
    }
    def __setattr__(self,name,value):
        if (name == "this") or (name == "thisown"): self.__dict__[name] = value; return
        method = rawbincand.__setmethods__.get(name,None)
        if method: return method(self,value)
        self.__dict__[name] = value
    __getmethods__ = {
        "full_N" : prestoc.rawbincand_full_N_get,
        "full_T" : prestoc.rawbincand_full_T_get,
        "full_lo_r" : prestoc.rawbincand_full_lo_r_get,
        "mini_N" : prestoc.rawbincand_mini_N_get,
        "mini_r" : prestoc.rawbincand_mini_r_get,
        "mini_power" : prestoc.rawbincand_mini_power_get,
        "mini_numsum" : prestoc.rawbincand_mini_numsum_get,
        "mini_sigma" : prestoc.rawbincand_mini_sigma_get,
        "psr_p" : prestoc.rawbincand_psr_p_get,
        "orb_p" : prestoc.rawbincand_orb_p_get,
    }
    def __getattr__(self,name):
        method = rawbincand.__getmethods__.get(name,None)
        if method: return method(self)
        raise AttributeError,name
    def __repr__(self):
        return "<C rawbincand instance at %s>" % (self.this,)
class rawbincandPtr(rawbincand):
    def __init__(self,this):
        self.this = this
        self.thisown = 0
        self.__class__ = rawbincand



class foldstats:
    def __init__(self,*args):
        self.this = apply(prestoc.new_foldstats,args)
        self.thisown = 1

    def __del__(self,prestoc=prestoc):
        if self.thisown == 1 :
            prestoc.delete_foldstats(self)
    __setmethods__ = {
        "numdata" : prestoc.foldstats_numdata_set,
        "data_avg" : prestoc.foldstats_data_avg_set,
        "data_var" : prestoc.foldstats_data_var_set,
        "numprof" : prestoc.foldstats_numprof_set,
        "prof_avg" : prestoc.foldstats_prof_avg_set,
        "prof_var" : prestoc.foldstats_prof_var_set,
        "redchi" : prestoc.foldstats_redchi_set,
    }
    def __setattr__(self,name,value):
        if (name == "this") or (name == "thisown"): self.__dict__[name] = value; return
        method = foldstats.__setmethods__.get(name,None)
        if method: return method(self,value)
        self.__dict__[name] = value
    __getmethods__ = {
        "numdata" : prestoc.foldstats_numdata_get,
        "data_avg" : prestoc.foldstats_data_avg_get,
        "data_var" : prestoc.foldstats_data_var_get,
        "numprof" : prestoc.foldstats_numprof_get,
        "prof_avg" : prestoc.foldstats_prof_avg_get,
        "prof_var" : prestoc.foldstats_prof_var_get,
        "redchi" : prestoc.foldstats_redchi_get,
    }
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

sigma_from_sumpows = prestoc.sigma_from_sumpows

sumpows_from_sigma = prestoc.sumpows_from_sigma

nice_output_1 = prestoc.nice_output_1

nice_output_2 = prestoc.nice_output_2

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

selectkth = prestoc.selectkth

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

corr_rz_plane = prestoc.corr_rz_plane

corr_rz_interp = prestoc.corr_rz_interp

tree_max_dm = prestoc.tree_max_dm

smearing_from_bw = prestoc.smearing_from_bw

delay_from_dm = prestoc.delay_from_dm

dm_from_delay = prestoc.dm_from_delay

dedisp_delays = prestoc.dedisp_delays

subband_search_delays = prestoc.subband_search_delays



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

import math, Numeric, Pgplot, string

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
    from Statistics import average
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

def ra_dec_to_string(h_or_d, m, s):
   """
   ra_dec_to_string(h_or_d, m, s):
      Return a formatted string of RA or DEC values as
      'hh:mm:ss.ssss' if RA, or 'dd:mm:ss.ssss' if DEC.
   """
   if (s >= 10.0):
      return "%.2d:%.2d:%.4f" % (h_or_d, m, s)
   else:
      return "%.2d:%.2d:0%.4f" % (h_or_d, m, s)

def ra_to_hours(ra_string):
   """
   ra_to_hours(ar_string):
      Given a string containing RA information as
      'hh:mm:ss.ssss', return the equivalent decimal
      hours.
   """
   h, m, s = string.split(ra_string, ":")
   h = int(h)
   m = int(m)
   s = float(s)
   return 12.0/PI * hms2rad(h, m, s)

def dec_to_deg(dec_string):
   """
   dec_to_deg(dec_string):
      Given a string containing DEC information as
      'dd:mm:ss.ssss', return the equivalent decimal
      degrees.
   """
   d, m, s = string.split(dec_string, ":")
   d = int(d)
   m = int(m)
   s = float(s)
   return RADTODEG * dms2rad(d, m, s)

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
   ra = ra_dec_to_string(obs.ra_h, obs.ra_m, obs.ra_s)
   dec = ra_dec_to_string(obs.dec_d, obs.dec_m, obs.dec_s)
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
