# This file was created automatically by SWIG.
# Don't modify this file, modify the SWIG interface instead.
# This file is compatible with both classic and new-style classes.
import _presto
def _swig_setattr(self,class_type,name,value):
    if (name == "this"):
        if isinstance(value, class_type):
            self.__dict__[name] = value.this
            if hasattr(value,"thisown"): self.__dict__["thisown"] = value.thisown
            del value.thisown
            return
    method = class_type.__swig_setmethods__.get(name,None)
    if method: return method(self,value)
    self.__dict__[name] = value

def _swig_getattr(self,class_type,name):
    method = class_type.__swig_getmethods__.get(name,None)
    if method: return method(self)
    raise AttributeError,name

import types
try:
    _object = types.ObjectType
    _newclass = 1
except AttributeError:
    class _object : pass
    _newclass = 0


SQRT2 = _presto.SQRT2
PI = _presto.PI
TWOPI = _presto.TWOPI
DEGTORAD = _presto.DEGTORAD
RADTODEG = _presto.RADTODEG
PIBYTWO = _presto.PIBYTWO
SOL = _presto.SOL
SECPERJULYR = _presto.SECPERJULYR
SECPERDAY = _presto.SECPERDAY
ARCSEC2RAD = _presto.ARCSEC2RAD
SEC2RAD = _presto.SEC2RAD
LOWACC = _presto.LOWACC
HIGHACC = _presto.HIGHACC
INTERBIN = _presto.INTERBIN
INTERPOLATE = _presto.INTERPOLATE
NO_CHECK_ALIASED = _presto.NO_CHECK_ALIASED
CHECK_ALIASED = _presto.CHECK_ALIASED
CONV = _presto.CONV
CORR = _presto.CORR
INPLACE_CONV = _presto.INPLACE_CONV
INPLACE_CORR = _presto.INPLACE_CORR
FFTDK = _presto.FFTDK
FFTD = _presto.FFTD
FFTK = _presto.FFTK
NOFFTS = _presto.NOFFTS
RAW = _presto.RAW
PREPPED = _presto.PREPPED
FFT = _presto.FFT
SAME = _presto.SAME
power_arr = _presto.power_arr

phase_arr = _presto.phase_arr

dpower_arr = _presto.dpower_arr

dphase_arr = _presto.dphase_arr

class orbitparams(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, orbitparams, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, orbitparams, name)
    __swig_setmethods__["p"] = _presto.orbitparams_p_set
    __swig_getmethods__["p"] = _presto.orbitparams_p_get
    if _newclass:p = property(_presto.orbitparams_p_get,_presto.orbitparams_p_set)
    __swig_setmethods__["e"] = _presto.orbitparams_e_set
    __swig_getmethods__["e"] = _presto.orbitparams_e_get
    if _newclass:e = property(_presto.orbitparams_e_get,_presto.orbitparams_e_set)
    __swig_setmethods__["x"] = _presto.orbitparams_x_set
    __swig_getmethods__["x"] = _presto.orbitparams_x_get
    if _newclass:x = property(_presto.orbitparams_x_get,_presto.orbitparams_x_set)
    __swig_setmethods__["w"] = _presto.orbitparams_w_set
    __swig_getmethods__["w"] = _presto.orbitparams_w_get
    if _newclass:w = property(_presto.orbitparams_w_get,_presto.orbitparams_w_set)
    __swig_setmethods__["t"] = _presto.orbitparams_t_set
    __swig_getmethods__["t"] = _presto.orbitparams_t_get
    if _newclass:t = property(_presto.orbitparams_t_get,_presto.orbitparams_t_set)
    __swig_setmethods__["pd"] = _presto.orbitparams_pd_set
    __swig_getmethods__["pd"] = _presto.orbitparams_pd_get
    if _newclass:pd = property(_presto.orbitparams_pd_get,_presto.orbitparams_pd_set)
    __swig_setmethods__["wd"] = _presto.orbitparams_wd_set
    __swig_getmethods__["wd"] = _presto.orbitparams_wd_get
    if _newclass:wd = property(_presto.orbitparams_wd_get,_presto.orbitparams_wd_set)
    def __init__(self,*args):
        self.this = apply(_presto.new_orbitparams,args)
        self.thisown = 1
    def __del__(self, destroy= _presto.delete_orbitparams):
        try:
            if self.thisown: destroy(self)
        except: pass
    def __repr__(self):
        return "<C orbitparams instance at %s>" % (self.this,)

class orbitparamsPtr(orbitparams):
    def __init__(self,this):
        self.this = this
        if not hasattr(self,"thisown"): self.thisown = 0
        self.__class__ = orbitparams
_presto.orbitparams_swigregister(orbitparamsPtr)
tofloatvector = _presto.tofloatvector

float_to_complex = _presto.float_to_complex

complex_to_float = _presto.complex_to_float


class DoubleArray(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, DoubleArray, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, DoubleArray, name)
    __swig_setmethods__["dptr"] = _presto.DoubleArray_dptr_set
    __swig_getmethods__["dptr"] = _presto.DoubleArray_dptr_get
    if _newclass:dptr = property(_presto.DoubleArray_dptr_get,_presto.DoubleArray_dptr_set)
    def __init__(self,*args):
        self.this = apply(_presto.new_DoubleArray,args)
        self.thisown = 1
    def __del__(self, destroy= _presto.delete_DoubleArray):
        try:
            if self.thisown: destroy(self)
        except: pass
    def __getitem__(*args): return apply(_presto.DoubleArray___getitem__,args)
    def __setitem__(*args): return apply(_presto.DoubleArray___setitem__,args)
    def __repr__(self):
        return "<C DoubleArray instance at %s>" % (self.this,)

class DoubleArrayPtr(DoubleArray):
    def __init__(self,this):
        self.this = this
        if not hasattr(self,"thisown"): self.thisown = 0
        self.__class__ = DoubleArray
_presto.DoubleArray_swigregister(DoubleArrayPtr)

class infodata(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, infodata, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, infodata, name)
    __swig_setmethods__["ra_s"] = _presto.infodata_ra_s_set
    __swig_getmethods__["ra_s"] = _presto.infodata_ra_s_get
    if _newclass:ra_s = property(_presto.infodata_ra_s_get,_presto.infodata_ra_s_set)
    __swig_setmethods__["dec_s"] = _presto.infodata_dec_s_set
    __swig_getmethods__["dec_s"] = _presto.infodata_dec_s_get
    if _newclass:dec_s = property(_presto.infodata_dec_s_get,_presto.infodata_dec_s_set)
    __swig_setmethods__["N"] = _presto.infodata_N_set
    __swig_getmethods__["N"] = _presto.infodata_N_get
    if _newclass:N = property(_presto.infodata_N_get,_presto.infodata_N_set)
    __swig_setmethods__["dt"] = _presto.infodata_dt_set
    __swig_getmethods__["dt"] = _presto.infodata_dt_get
    if _newclass:dt = property(_presto.infodata_dt_get,_presto.infodata_dt_set)
    __swig_setmethods__["fov"] = _presto.infodata_fov_set
    __swig_getmethods__["fov"] = _presto.infodata_fov_get
    if _newclass:fov = property(_presto.infodata_fov_get,_presto.infodata_fov_set)
    __swig_setmethods__["mjd_f"] = _presto.infodata_mjd_f_set
    __swig_getmethods__["mjd_f"] = _presto.infodata_mjd_f_get
    if _newclass:mjd_f = property(_presto.infodata_mjd_f_get,_presto.infodata_mjd_f_set)
    __swig_setmethods__["dm"] = _presto.infodata_dm_set
    __swig_getmethods__["dm"] = _presto.infodata_dm_get
    if _newclass:dm = property(_presto.infodata_dm_get,_presto.infodata_dm_set)
    __swig_setmethods__["freq"] = _presto.infodata_freq_set
    __swig_getmethods__["freq"] = _presto.infodata_freq_get
    if _newclass:freq = property(_presto.infodata_freq_get,_presto.infodata_freq_set)
    __swig_setmethods__["freqband"] = _presto.infodata_freqband_set
    __swig_getmethods__["freqband"] = _presto.infodata_freqband_get
    if _newclass:freqband = property(_presto.infodata_freqband_get,_presto.infodata_freqband_set)
    __swig_setmethods__["chan_wid"] = _presto.infodata_chan_wid_set
    __swig_getmethods__["chan_wid"] = _presto.infodata_chan_wid_get
    if _newclass:chan_wid = property(_presto.infodata_chan_wid_get,_presto.infodata_chan_wid_set)
    __swig_setmethods__["wavelen"] = _presto.infodata_wavelen_set
    __swig_getmethods__["wavelen"] = _presto.infodata_wavelen_get
    if _newclass:wavelen = property(_presto.infodata_wavelen_get,_presto.infodata_wavelen_set)
    __swig_setmethods__["waveband"] = _presto.infodata_waveband_set
    __swig_getmethods__["waveband"] = _presto.infodata_waveband_get
    if _newclass:waveband = property(_presto.infodata_waveband_get,_presto.infodata_waveband_set)
    __swig_setmethods__["energy"] = _presto.infodata_energy_set
    __swig_getmethods__["energy"] = _presto.infodata_energy_get
    if _newclass:energy = property(_presto.infodata_energy_get,_presto.infodata_energy_set)
    __swig_setmethods__["energyband"] = _presto.infodata_energyband_set
    __swig_getmethods__["energyband"] = _presto.infodata_energyband_get
    if _newclass:energyband = property(_presto.infodata_energyband_get,_presto.infodata_energyband_set)
    __swig_setmethods__["onoff"] = _presto.infodata_onoff_set
    __swig_getmethods__["onoff"] = _presto.infodata_onoff_get
    if _newclass:onoff = property(_presto.infodata_onoff_get,_presto.infodata_onoff_set)
    __swig_setmethods__["num_chan"] = _presto.infodata_num_chan_set
    __swig_getmethods__["num_chan"] = _presto.infodata_num_chan_get
    if _newclass:num_chan = property(_presto.infodata_num_chan_get,_presto.infodata_num_chan_set)
    __swig_setmethods__["mjd_i"] = _presto.infodata_mjd_i_set
    __swig_getmethods__["mjd_i"] = _presto.infodata_mjd_i_get
    if _newclass:mjd_i = property(_presto.infodata_mjd_i_get,_presto.infodata_mjd_i_set)
    __swig_setmethods__["ra_h"] = _presto.infodata_ra_h_set
    __swig_getmethods__["ra_h"] = _presto.infodata_ra_h_get
    if _newclass:ra_h = property(_presto.infodata_ra_h_get,_presto.infodata_ra_h_set)
    __swig_setmethods__["ra_m"] = _presto.infodata_ra_m_set
    __swig_getmethods__["ra_m"] = _presto.infodata_ra_m_get
    if _newclass:ra_m = property(_presto.infodata_ra_m_get,_presto.infodata_ra_m_set)
    __swig_setmethods__["dec_d"] = _presto.infodata_dec_d_set
    __swig_getmethods__["dec_d"] = _presto.infodata_dec_d_get
    if _newclass:dec_d = property(_presto.infodata_dec_d_get,_presto.infodata_dec_d_set)
    __swig_setmethods__["dec_m"] = _presto.infodata_dec_m_set
    __swig_getmethods__["dec_m"] = _presto.infodata_dec_m_get
    if _newclass:dec_m = property(_presto.infodata_dec_m_get,_presto.infodata_dec_m_set)
    __swig_setmethods__["bary"] = _presto.infodata_bary_set
    __swig_getmethods__["bary"] = _presto.infodata_bary_get
    if _newclass:bary = property(_presto.infodata_bary_get,_presto.infodata_bary_set)
    __swig_setmethods__["numonoff"] = _presto.infodata_numonoff_set
    __swig_getmethods__["numonoff"] = _presto.infodata_numonoff_get
    if _newclass:numonoff = property(_presto.infodata_numonoff_get,_presto.infodata_numonoff_set)
    __swig_setmethods__["notes"] = _presto.infodata_notes_set
    __swig_getmethods__["notes"] = _presto.infodata_notes_get
    if _newclass:notes = property(_presto.infodata_notes_get,_presto.infodata_notes_set)
    __swig_setmethods__["name"] = _presto.infodata_name_set
    __swig_getmethods__["name"] = _presto.infodata_name_get
    if _newclass:name = property(_presto.infodata_name_get,_presto.infodata_name_set)
    __swig_setmethods__["object"] = _presto.infodata_object_set
    __swig_getmethods__["object"] = _presto.infodata_object_get
    if _newclass:object = property(_presto.infodata_object_get,_presto.infodata_object_set)
    __swig_setmethods__["instrument"] = _presto.infodata_instrument_set
    __swig_getmethods__["instrument"] = _presto.infodata_instrument_get
    if _newclass:instrument = property(_presto.infodata_instrument_get,_presto.infodata_instrument_set)
    __swig_setmethods__["observer"] = _presto.infodata_observer_set
    __swig_getmethods__["observer"] = _presto.infodata_observer_get
    if _newclass:observer = property(_presto.infodata_observer_get,_presto.infodata_observer_set)
    __swig_setmethods__["analyzer"] = _presto.infodata_analyzer_set
    __swig_getmethods__["analyzer"] = _presto.infodata_analyzer_get
    if _newclass:analyzer = property(_presto.infodata_analyzer_get,_presto.infodata_analyzer_set)
    __swig_setmethods__["telescope"] = _presto.infodata_telescope_set
    __swig_getmethods__["telescope"] = _presto.infodata_telescope_get
    if _newclass:telescope = property(_presto.infodata_telescope_get,_presto.infodata_telescope_set)
    __swig_setmethods__["band"] = _presto.infodata_band_set
    __swig_getmethods__["band"] = _presto.infodata_band_get
    if _newclass:band = property(_presto.infodata_band_get,_presto.infodata_band_set)
    __swig_setmethods__["filt"] = _presto.infodata_filt_set
    __swig_getmethods__["filt"] = _presto.infodata_filt_get
    if _newclass:filt = property(_presto.infodata_filt_get,_presto.infodata_filt_set)
    def __init__(self,*args):
        self.this = apply(_presto.new_infodata,args)
        self.thisown = 1
    def __del__(self, destroy= _presto.delete_infodata):
        try:
            if self.thisown: destroy(self)
        except: pass
    def __repr__(self):
        return "<C infodata instance at %s>" % (self.this,)

class infodataPtr(infodata):
    def __init__(self,this):
        self.this = this
        if not hasattr(self,"thisown"): self.thisown = 0
        self.__class__ = infodata
_presto.infodata_swigregister(infodataPtr)

readinf = _presto.readinf

writeinf = _presto.writeinf

class makedata(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, makedata, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, makedata, name)
    __swig_setmethods__["basefilenm"] = _presto.makedata_basefilenm_set
    __swig_getmethods__["basefilenm"] = _presto.makedata_basefilenm_get
    if _newclass:basefilenm = property(_presto.makedata_basefilenm_get,_presto.makedata_basefilenm_set)
    __swig_setmethods__["description"] = _presto.makedata_description_set
    __swig_getmethods__["description"] = _presto.makedata_description_get
    if _newclass:description = property(_presto.makedata_description_get,_presto.makedata_description_set)
    __swig_setmethods__["N"] = _presto.makedata_N_set
    __swig_getmethods__["N"] = _presto.makedata_N_get
    if _newclass:N = property(_presto.makedata_N_get,_presto.makedata_N_set)
    __swig_setmethods__["next2_to_n"] = _presto.makedata_next2_to_n_set
    __swig_getmethods__["next2_to_n"] = _presto.makedata_next2_to_n_get
    if _newclass:next2_to_n = property(_presto.makedata_next2_to_n_get,_presto.makedata_next2_to_n_set)
    __swig_setmethods__["dt"] = _presto.makedata_dt_set
    __swig_getmethods__["dt"] = _presto.makedata_dt_get
    if _newclass:dt = property(_presto.makedata_dt_get,_presto.makedata_dt_set)
    __swig_setmethods__["T"] = _presto.makedata_T_set
    __swig_getmethods__["T"] = _presto.makedata_T_get
    if _newclass:T = property(_presto.makedata_T_get,_presto.makedata_T_set)
    __swig_setmethods__["ptype"] = _presto.makedata_ptype_set
    __swig_getmethods__["ptype"] = _presto.makedata_ptype_get
    if _newclass:ptype = property(_presto.makedata_ptype_get,_presto.makedata_ptype_set)
    __swig_setmethods__["pnum"] = _presto.makedata_pnum_set
    __swig_getmethods__["pnum"] = _presto.makedata_pnum_get
    if _newclass:pnum = property(_presto.makedata_pnum_get,_presto.makedata_pnum_set)
    __swig_setmethods__["fwhm"] = _presto.makedata_fwhm_set
    __swig_getmethods__["fwhm"] = _presto.makedata_fwhm_get
    if _newclass:fwhm = property(_presto.makedata_fwhm_get,_presto.makedata_fwhm_set)
    __swig_setmethods__["round"] = _presto.makedata_round_set
    __swig_getmethods__["round"] = _presto.makedata_round_get
    if _newclass:round = property(_presto.makedata_round_get,_presto.makedata_round_set)
    __swig_setmethods__["roundnum"] = _presto.makedata_roundnum_set
    __swig_getmethods__["roundnum"] = _presto.makedata_roundnum_get
    if _newclass:roundnum = property(_presto.makedata_roundnum_get,_presto.makedata_roundnum_set)
    __swig_setmethods__["f"] = _presto.makedata_f_set
    __swig_getmethods__["f"] = _presto.makedata_f_get
    if _newclass:f = property(_presto.makedata_f_get,_presto.makedata_f_set)
    __swig_setmethods__["fd"] = _presto.makedata_fd_set
    __swig_getmethods__["fd"] = _presto.makedata_fd_get
    if _newclass:fd = property(_presto.makedata_fd_get,_presto.makedata_fd_set)
    __swig_setmethods__["fdd"] = _presto.makedata_fdd_set
    __swig_getmethods__["fdd"] = _presto.makedata_fdd_get
    if _newclass:fdd = property(_presto.makedata_fdd_get,_presto.makedata_fdd_set)
    __swig_setmethods__["p"] = _presto.makedata_p_set
    __swig_getmethods__["p"] = _presto.makedata_p_get
    if _newclass:p = property(_presto.makedata_p_get,_presto.makedata_p_set)
    __swig_setmethods__["pd"] = _presto.makedata_pd_set
    __swig_getmethods__["pd"] = _presto.makedata_pd_get
    if _newclass:pd = property(_presto.makedata_pd_get,_presto.makedata_pd_set)
    __swig_setmethods__["pdd"] = _presto.makedata_pdd_set
    __swig_getmethods__["pdd"] = _presto.makedata_pdd_get
    if _newclass:pdd = property(_presto.makedata_pdd_get,_presto.makedata_pdd_set)
    __swig_setmethods__["r"] = _presto.makedata_r_set
    __swig_getmethods__["r"] = _presto.makedata_r_get
    if _newclass:r = property(_presto.makedata_r_get,_presto.makedata_r_set)
    __swig_setmethods__["z"] = _presto.makedata_z_set
    __swig_getmethods__["z"] = _presto.makedata_z_get
    if _newclass:z = property(_presto.makedata_z_get,_presto.makedata_z_set)
    __swig_setmethods__["w"] = _presto.makedata_w_set
    __swig_getmethods__["w"] = _presto.makedata_w_get
    if _newclass:w = property(_presto.makedata_w_get,_presto.makedata_w_set)
    __swig_setmethods__["amp"] = _presto.makedata_amp_set
    __swig_getmethods__["amp"] = _presto.makedata_amp_get
    if _newclass:amp = property(_presto.makedata_amp_get,_presto.makedata_amp_set)
    __swig_setmethods__["phs"] = _presto.makedata_phs_set
    __swig_getmethods__["phs"] = _presto.makedata_phs_get
    if _newclass:phs = property(_presto.makedata_phs_get,_presto.makedata_phs_set)
    __swig_setmethods__["dc"] = _presto.makedata_dc_set
    __swig_getmethods__["dc"] = _presto.makedata_dc_get
    if _newclass:dc = property(_presto.makedata_dc_get,_presto.makedata_dc_set)
    __swig_setmethods__["binary"] = _presto.makedata_binary_set
    __swig_getmethods__["binary"] = _presto.makedata_binary_get
    if _newclass:binary = property(_presto.makedata_binary_get,_presto.makedata_binary_set)
    __swig_setmethods__["orb"] = _presto.makedata_orb_set
    __swig_getmethods__["orb"] = _presto.makedata_orb_get
    if _newclass:orb = property(_presto.makedata_orb_get,_presto.makedata_orb_set)
    __swig_setmethods__["ampmod"] = _presto.makedata_ampmod_set
    __swig_getmethods__["ampmod"] = _presto.makedata_ampmod_get
    if _newclass:ampmod = property(_presto.makedata_ampmod_get,_presto.makedata_ampmod_set)
    __swig_setmethods__["ampmoda"] = _presto.makedata_ampmoda_set
    __swig_getmethods__["ampmoda"] = _presto.makedata_ampmoda_get
    if _newclass:ampmoda = property(_presto.makedata_ampmoda_get,_presto.makedata_ampmoda_set)
    __swig_setmethods__["ampmodf"] = _presto.makedata_ampmodf_set
    __swig_getmethods__["ampmodf"] = _presto.makedata_ampmodf_get
    if _newclass:ampmodf = property(_presto.makedata_ampmodf_get,_presto.makedata_ampmodf_set)
    __swig_setmethods__["ampmodp"] = _presto.makedata_ampmodp_set
    __swig_getmethods__["ampmodp"] = _presto.makedata_ampmodp_get
    if _newclass:ampmodp = property(_presto.makedata_ampmodp_get,_presto.makedata_ampmodp_set)
    __swig_setmethods__["noisetype"] = _presto.makedata_noisetype_set
    __swig_getmethods__["noisetype"] = _presto.makedata_noisetype_get
    if _newclass:noisetype = property(_presto.makedata_noisetype_get,_presto.makedata_noisetype_set)
    __swig_setmethods__["noise"] = _presto.makedata_noise_set
    __swig_getmethods__["noise"] = _presto.makedata_noise_get
    if _newclass:noise = property(_presto.makedata_noise_get,_presto.makedata_noise_set)
    __swig_setmethods__["noisesig"] = _presto.makedata_noisesig_set
    __swig_getmethods__["noisesig"] = _presto.makedata_noisesig_get
    if _newclass:noisesig = property(_presto.makedata_noisesig_get,_presto.makedata_noisesig_set)
    __swig_setmethods__["numonoff"] = _presto.makedata_numonoff_set
    __swig_getmethods__["numonoff"] = _presto.makedata_numonoff_get
    if _newclass:numonoff = property(_presto.makedata_numonoff_get,_presto.makedata_numonoff_set)
    __swig_setmethods__["onoff"] = _presto.makedata_onoff_set
    __swig_getmethods__["onoff"] = _presto.makedata_onoff_get
    if _newclass:onoff = property(_presto.makedata_onoff_get,_presto.makedata_onoff_set)
    def __init__(self,*args):
        self.this = apply(_presto.new_makedata,args)
        self.thisown = 1
    def __del__(self, destroy= _presto.delete_makedata):
        try:
            if self.thisown: destroy(self)
        except: pass
    def __repr__(self):
        return "<C makedata instance at %s>" % (self.this,)

class makedataPtr(makedata):
    def __init__(self,this):
        self.this = this
        if not hasattr(self,"thisown"): self.thisown = 0
        self.__class__ = makedata
_presto.makedata_swigregister(makedataPtr)

read_mak_input = _presto.read_mak_input

read_mak_file = _presto.read_mak_file

write_mak_file = _presto.write_mak_file

class rderivs(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, rderivs, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, rderivs, name)
    __swig_setmethods__["pow"] = _presto.rderivs_pow_set
    __swig_getmethods__["pow"] = _presto.rderivs_pow_get
    if _newclass:pow = property(_presto.rderivs_pow_get,_presto.rderivs_pow_set)
    __swig_setmethods__["phs"] = _presto.rderivs_phs_set
    __swig_getmethods__["phs"] = _presto.rderivs_phs_get
    if _newclass:phs = property(_presto.rderivs_phs_get,_presto.rderivs_phs_set)
    __swig_setmethods__["dpow"] = _presto.rderivs_dpow_set
    __swig_getmethods__["dpow"] = _presto.rderivs_dpow_get
    if _newclass:dpow = property(_presto.rderivs_dpow_get,_presto.rderivs_dpow_set)
    __swig_setmethods__["dphs"] = _presto.rderivs_dphs_set
    __swig_getmethods__["dphs"] = _presto.rderivs_dphs_get
    if _newclass:dphs = property(_presto.rderivs_dphs_get,_presto.rderivs_dphs_set)
    __swig_setmethods__["d2pow"] = _presto.rderivs_d2pow_set
    __swig_getmethods__["d2pow"] = _presto.rderivs_d2pow_get
    if _newclass:d2pow = property(_presto.rderivs_d2pow_get,_presto.rderivs_d2pow_set)
    __swig_setmethods__["d2phs"] = _presto.rderivs_d2phs_set
    __swig_getmethods__["d2phs"] = _presto.rderivs_d2phs_get
    if _newclass:d2phs = property(_presto.rderivs_d2phs_get,_presto.rderivs_d2phs_set)
    __swig_setmethods__["locpow"] = _presto.rderivs_locpow_set
    __swig_getmethods__["locpow"] = _presto.rderivs_locpow_get
    if _newclass:locpow = property(_presto.rderivs_locpow_get,_presto.rderivs_locpow_set)
    def __init__(self,*args):
        self.this = apply(_presto.new_rderivs,args)
        self.thisown = 1
    def __del__(self, destroy= _presto.delete_rderivs):
        try:
            if self.thisown: destroy(self)
        except: pass
    def __repr__(self):
        return "<C rderivs instance at %s>" % (self.this,)

class rderivsPtr(rderivs):
    def __init__(self,this):
        self.this = this
        if not hasattr(self,"thisown"): self.thisown = 0
        self.__class__ = rderivs
_presto.rderivs_swigregister(rderivsPtr)

class fourierprops(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, fourierprops, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, fourierprops, name)
    __swig_setmethods__["r"] = _presto.fourierprops_r_set
    __swig_getmethods__["r"] = _presto.fourierprops_r_get
    if _newclass:r = property(_presto.fourierprops_r_get,_presto.fourierprops_r_set)
    __swig_setmethods__["rerr"] = _presto.fourierprops_rerr_set
    __swig_getmethods__["rerr"] = _presto.fourierprops_rerr_get
    if _newclass:rerr = property(_presto.fourierprops_rerr_get,_presto.fourierprops_rerr_set)
    __swig_setmethods__["z"] = _presto.fourierprops_z_set
    __swig_getmethods__["z"] = _presto.fourierprops_z_get
    if _newclass:z = property(_presto.fourierprops_z_get,_presto.fourierprops_z_set)
    __swig_setmethods__["zerr"] = _presto.fourierprops_zerr_set
    __swig_getmethods__["zerr"] = _presto.fourierprops_zerr_get
    if _newclass:zerr = property(_presto.fourierprops_zerr_get,_presto.fourierprops_zerr_set)
    __swig_setmethods__["w"] = _presto.fourierprops_w_set
    __swig_getmethods__["w"] = _presto.fourierprops_w_get
    if _newclass:w = property(_presto.fourierprops_w_get,_presto.fourierprops_w_set)
    __swig_setmethods__["werr"] = _presto.fourierprops_werr_set
    __swig_getmethods__["werr"] = _presto.fourierprops_werr_get
    if _newclass:werr = property(_presto.fourierprops_werr_get,_presto.fourierprops_werr_set)
    __swig_setmethods__["pow"] = _presto.fourierprops_pow_set
    __swig_getmethods__["pow"] = _presto.fourierprops_pow_get
    if _newclass:pow = property(_presto.fourierprops_pow_get,_presto.fourierprops_pow_set)
    __swig_setmethods__["powerr"] = _presto.fourierprops_powerr_set
    __swig_getmethods__["powerr"] = _presto.fourierprops_powerr_get
    if _newclass:powerr = property(_presto.fourierprops_powerr_get,_presto.fourierprops_powerr_set)
    __swig_setmethods__["sig"] = _presto.fourierprops_sig_set
    __swig_getmethods__["sig"] = _presto.fourierprops_sig_get
    if _newclass:sig = property(_presto.fourierprops_sig_get,_presto.fourierprops_sig_set)
    __swig_setmethods__["rawpow"] = _presto.fourierprops_rawpow_set
    __swig_getmethods__["rawpow"] = _presto.fourierprops_rawpow_get
    if _newclass:rawpow = property(_presto.fourierprops_rawpow_get,_presto.fourierprops_rawpow_set)
    __swig_setmethods__["phs"] = _presto.fourierprops_phs_set
    __swig_getmethods__["phs"] = _presto.fourierprops_phs_get
    if _newclass:phs = property(_presto.fourierprops_phs_get,_presto.fourierprops_phs_set)
    __swig_setmethods__["phserr"] = _presto.fourierprops_phserr_set
    __swig_getmethods__["phserr"] = _presto.fourierprops_phserr_get
    if _newclass:phserr = property(_presto.fourierprops_phserr_get,_presto.fourierprops_phserr_set)
    __swig_setmethods__["cen"] = _presto.fourierprops_cen_set
    __swig_getmethods__["cen"] = _presto.fourierprops_cen_get
    if _newclass:cen = property(_presto.fourierprops_cen_get,_presto.fourierprops_cen_set)
    __swig_setmethods__["cenerr"] = _presto.fourierprops_cenerr_set
    __swig_getmethods__["cenerr"] = _presto.fourierprops_cenerr_get
    if _newclass:cenerr = property(_presto.fourierprops_cenerr_get,_presto.fourierprops_cenerr_set)
    __swig_setmethods__["pur"] = _presto.fourierprops_pur_set
    __swig_getmethods__["pur"] = _presto.fourierprops_pur_get
    if _newclass:pur = property(_presto.fourierprops_pur_get,_presto.fourierprops_pur_set)
    __swig_setmethods__["purerr"] = _presto.fourierprops_purerr_set
    __swig_getmethods__["purerr"] = _presto.fourierprops_purerr_get
    if _newclass:purerr = property(_presto.fourierprops_purerr_get,_presto.fourierprops_purerr_set)
    __swig_setmethods__["locpow"] = _presto.fourierprops_locpow_set
    __swig_getmethods__["locpow"] = _presto.fourierprops_locpow_get
    if _newclass:locpow = property(_presto.fourierprops_locpow_get,_presto.fourierprops_locpow_set)
    def __init__(self,*args):
        self.this = apply(_presto.new_fourierprops,args)
        self.thisown = 1
    def __del__(self, destroy= _presto.delete_fourierprops):
        try:
            if self.thisown: destroy(self)
        except: pass
    def __repr__(self):
        return "<C fourierprops instance at %s>" % (self.this,)

class fourierpropsPtr(fourierprops):
    def __init__(self,this):
        self.this = this
        if not hasattr(self,"thisown"): self.thisown = 0
        self.__class__ = fourierprops
_presto.fourierprops_swigregister(fourierpropsPtr)

class binaryprops(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, binaryprops, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, binaryprops, name)
    __swig_setmethods__["ppsr"] = _presto.binaryprops_ppsr_set
    __swig_getmethods__["ppsr"] = _presto.binaryprops_ppsr_get
    if _newclass:ppsr = property(_presto.binaryprops_ppsr_get,_presto.binaryprops_ppsr_set)
    __swig_setmethods__["fpsr"] = _presto.binaryprops_fpsr_set
    __swig_getmethods__["fpsr"] = _presto.binaryprops_fpsr_get
    if _newclass:fpsr = property(_presto.binaryprops_fpsr_get,_presto.binaryprops_fpsr_set)
    __swig_setmethods__["rpsr"] = _presto.binaryprops_rpsr_set
    __swig_getmethods__["rpsr"] = _presto.binaryprops_rpsr_get
    if _newclass:rpsr = property(_presto.binaryprops_rpsr_get,_presto.binaryprops_rpsr_set)
    __swig_setmethods__["pbin"] = _presto.binaryprops_pbin_set
    __swig_getmethods__["pbin"] = _presto.binaryprops_pbin_get
    if _newclass:pbin = property(_presto.binaryprops_pbin_get,_presto.binaryprops_pbin_set)
    __swig_setmethods__["rbin"] = _presto.binaryprops_rbin_set
    __swig_getmethods__["rbin"] = _presto.binaryprops_rbin_get
    if _newclass:rbin = property(_presto.binaryprops_rbin_get,_presto.binaryprops_rbin_set)
    __swig_setmethods__["z"] = _presto.binaryprops_z_set
    __swig_getmethods__["z"] = _presto.binaryprops_z_get
    if _newclass:z = property(_presto.binaryprops_z_get,_presto.binaryprops_z_set)
    __swig_setmethods__["asinic"] = _presto.binaryprops_asinic_set
    __swig_getmethods__["asinic"] = _presto.binaryprops_asinic_get
    if _newclass:asinic = property(_presto.binaryprops_asinic_get,_presto.binaryprops_asinic_set)
    __swig_setmethods__["rdetect"] = _presto.binaryprops_rdetect_set
    __swig_getmethods__["rdetect"] = _presto.binaryprops_rdetect_get
    if _newclass:rdetect = property(_presto.binaryprops_rdetect_get,_presto.binaryprops_rdetect_set)
    __swig_setmethods__["nfftbins"] = _presto.binaryprops_nfftbins_set
    __swig_getmethods__["nfftbins"] = _presto.binaryprops_nfftbins_get
    if _newclass:nfftbins = property(_presto.binaryprops_nfftbins_get,_presto.binaryprops_nfftbins_set)
    __swig_setmethods__["lowbin"] = _presto.binaryprops_lowbin_set
    __swig_getmethods__["lowbin"] = _presto.binaryprops_lowbin_get
    if _newclass:lowbin = property(_presto.binaryprops_lowbin_get,_presto.binaryprops_lowbin_set)
    __swig_setmethods__["ppsrerr"] = _presto.binaryprops_ppsrerr_set
    __swig_getmethods__["ppsrerr"] = _presto.binaryprops_ppsrerr_get
    if _newclass:ppsrerr = property(_presto.binaryprops_ppsrerr_get,_presto.binaryprops_ppsrerr_set)
    __swig_setmethods__["fpsrerr"] = _presto.binaryprops_fpsrerr_set
    __swig_getmethods__["fpsrerr"] = _presto.binaryprops_fpsrerr_get
    if _newclass:fpsrerr = property(_presto.binaryprops_fpsrerr_get,_presto.binaryprops_fpsrerr_set)
    __swig_setmethods__["rpsrerr"] = _presto.binaryprops_rpsrerr_set
    __swig_getmethods__["rpsrerr"] = _presto.binaryprops_rpsrerr_get
    if _newclass:rpsrerr = property(_presto.binaryprops_rpsrerr_get,_presto.binaryprops_rpsrerr_set)
    __swig_setmethods__["pbinerr"] = _presto.binaryprops_pbinerr_set
    __swig_getmethods__["pbinerr"] = _presto.binaryprops_pbinerr_get
    if _newclass:pbinerr = property(_presto.binaryprops_pbinerr_get,_presto.binaryprops_pbinerr_set)
    __swig_setmethods__["rbinerr"] = _presto.binaryprops_rbinerr_set
    __swig_getmethods__["rbinerr"] = _presto.binaryprops_rbinerr_get
    if _newclass:rbinerr = property(_presto.binaryprops_rbinerr_get,_presto.binaryprops_rbinerr_set)
    __swig_setmethods__["zerr"] = _presto.binaryprops_zerr_set
    __swig_getmethods__["zerr"] = _presto.binaryprops_zerr_get
    if _newclass:zerr = property(_presto.binaryprops_zerr_get,_presto.binaryprops_zerr_set)
    __swig_setmethods__["asinicerr"] = _presto.binaryprops_asinicerr_set
    __swig_getmethods__["asinicerr"] = _presto.binaryprops_asinicerr_get
    if _newclass:asinicerr = property(_presto.binaryprops_asinicerr_get,_presto.binaryprops_asinicerr_set)
    __swig_setmethods__["rdetecterr"] = _presto.binaryprops_rdetecterr_set
    __swig_getmethods__["rdetecterr"] = _presto.binaryprops_rdetecterr_get
    if _newclass:rdetecterr = property(_presto.binaryprops_rdetecterr_get,_presto.binaryprops_rdetecterr_set)
    __swig_setmethods__["sig"] = _presto.binaryprops_sig_set
    __swig_getmethods__["sig"] = _presto.binaryprops_sig_get
    if _newclass:sig = property(_presto.binaryprops_sig_get,_presto.binaryprops_sig_set)
    __swig_setmethods__["phs"] = _presto.binaryprops_phs_set
    __swig_getmethods__["phs"] = _presto.binaryprops_phs_get
    if _newclass:phs = property(_presto.binaryprops_phs_get,_presto.binaryprops_phs_set)
    __swig_setmethods__["phserr"] = _presto.binaryprops_phserr_set
    __swig_getmethods__["phserr"] = _presto.binaryprops_phserr_get
    if _newclass:phserr = property(_presto.binaryprops_phserr_get,_presto.binaryprops_phserr_set)
    __swig_setmethods__["cen"] = _presto.binaryprops_cen_set
    __swig_getmethods__["cen"] = _presto.binaryprops_cen_get
    if _newclass:cen = property(_presto.binaryprops_cen_get,_presto.binaryprops_cen_set)
    __swig_setmethods__["cenerr"] = _presto.binaryprops_cenerr_set
    __swig_getmethods__["cenerr"] = _presto.binaryprops_cenerr_get
    if _newclass:cenerr = property(_presto.binaryprops_cenerr_get,_presto.binaryprops_cenerr_set)
    __swig_setmethods__["pur"] = _presto.binaryprops_pur_set
    __swig_getmethods__["pur"] = _presto.binaryprops_pur_get
    if _newclass:pur = property(_presto.binaryprops_pur_get,_presto.binaryprops_pur_set)
    __swig_setmethods__["purerr"] = _presto.binaryprops_purerr_set
    __swig_getmethods__["purerr"] = _presto.binaryprops_purerr_get
    if _newclass:purerr = property(_presto.binaryprops_purerr_get,_presto.binaryprops_purerr_set)
    __swig_setmethods__["pow"] = _presto.binaryprops_pow_set
    __swig_getmethods__["pow"] = _presto.binaryprops_pow_get
    if _newclass:pow = property(_presto.binaryprops_pow_get,_presto.binaryprops_pow_set)
    __swig_setmethods__["powerr"] = _presto.binaryprops_powerr_set
    __swig_getmethods__["powerr"] = _presto.binaryprops_powerr_get
    if _newclass:powerr = property(_presto.binaryprops_powerr_get,_presto.binaryprops_powerr_set)
    def __init__(self,*args):
        self.this = apply(_presto.new_binaryprops,args)
        self.thisown = 1
    def __del__(self, destroy= _presto.delete_binaryprops):
        try:
            if self.thisown: destroy(self)
        except: pass
    def __repr__(self):
        return "<C binaryprops instance at %s>" % (self.this,)

class binarypropsPtr(binaryprops):
    def __init__(self,this):
        self.this = this
        if not hasattr(self,"thisown"): self.thisown = 0
        self.__class__ = binaryprops
_presto.binaryprops_swigregister(binarypropsPtr)

class rawbincand(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, rawbincand, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, rawbincand, name)
    __swig_setmethods__["full_N"] = _presto.rawbincand_full_N_set
    __swig_getmethods__["full_N"] = _presto.rawbincand_full_N_get
    if _newclass:full_N = property(_presto.rawbincand_full_N_get,_presto.rawbincand_full_N_set)
    __swig_setmethods__["full_T"] = _presto.rawbincand_full_T_set
    __swig_getmethods__["full_T"] = _presto.rawbincand_full_T_get
    if _newclass:full_T = property(_presto.rawbincand_full_T_get,_presto.rawbincand_full_T_set)
    __swig_setmethods__["full_lo_r"] = _presto.rawbincand_full_lo_r_set
    __swig_getmethods__["full_lo_r"] = _presto.rawbincand_full_lo_r_get
    if _newclass:full_lo_r = property(_presto.rawbincand_full_lo_r_get,_presto.rawbincand_full_lo_r_set)
    __swig_setmethods__["mini_N"] = _presto.rawbincand_mini_N_set
    __swig_getmethods__["mini_N"] = _presto.rawbincand_mini_N_get
    if _newclass:mini_N = property(_presto.rawbincand_mini_N_get,_presto.rawbincand_mini_N_set)
    __swig_setmethods__["mini_r"] = _presto.rawbincand_mini_r_set
    __swig_getmethods__["mini_r"] = _presto.rawbincand_mini_r_get
    if _newclass:mini_r = property(_presto.rawbincand_mini_r_get,_presto.rawbincand_mini_r_set)
    __swig_setmethods__["mini_power"] = _presto.rawbincand_mini_power_set
    __swig_getmethods__["mini_power"] = _presto.rawbincand_mini_power_get
    if _newclass:mini_power = property(_presto.rawbincand_mini_power_get,_presto.rawbincand_mini_power_set)
    __swig_setmethods__["mini_numsum"] = _presto.rawbincand_mini_numsum_set
    __swig_getmethods__["mini_numsum"] = _presto.rawbincand_mini_numsum_get
    if _newclass:mini_numsum = property(_presto.rawbincand_mini_numsum_get,_presto.rawbincand_mini_numsum_set)
    __swig_setmethods__["mini_sigma"] = _presto.rawbincand_mini_sigma_set
    __swig_getmethods__["mini_sigma"] = _presto.rawbincand_mini_sigma_get
    if _newclass:mini_sigma = property(_presto.rawbincand_mini_sigma_get,_presto.rawbincand_mini_sigma_set)
    __swig_setmethods__["psr_p"] = _presto.rawbincand_psr_p_set
    __swig_getmethods__["psr_p"] = _presto.rawbincand_psr_p_get
    if _newclass:psr_p = property(_presto.rawbincand_psr_p_get,_presto.rawbincand_psr_p_set)
    __swig_setmethods__["orb_p"] = _presto.rawbincand_orb_p_set
    __swig_getmethods__["orb_p"] = _presto.rawbincand_orb_p_get
    if _newclass:orb_p = property(_presto.rawbincand_orb_p_get,_presto.rawbincand_orb_p_set)
    def __init__(self,*args):
        self.this = apply(_presto.new_rawbincand,args)
        self.thisown = 1
    def __del__(self, destroy= _presto.delete_rawbincand):
        try:
            if self.thisown: destroy(self)
        except: pass
    def __repr__(self):
        return "<C rawbincand instance at %s>" % (self.this,)

class rawbincandPtr(rawbincand):
    def __init__(self,this):
        self.this = this
        if not hasattr(self,"thisown"): self.thisown = 0
        self.__class__ = rawbincand
_presto.rawbincand_swigregister(rawbincandPtr)

class foldstats(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, foldstats, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, foldstats, name)
    __swig_setmethods__["numdata"] = _presto.foldstats_numdata_set
    __swig_getmethods__["numdata"] = _presto.foldstats_numdata_get
    if _newclass:numdata = property(_presto.foldstats_numdata_get,_presto.foldstats_numdata_set)
    __swig_setmethods__["data_avg"] = _presto.foldstats_data_avg_set
    __swig_getmethods__["data_avg"] = _presto.foldstats_data_avg_get
    if _newclass:data_avg = property(_presto.foldstats_data_avg_get,_presto.foldstats_data_avg_set)
    __swig_setmethods__["data_var"] = _presto.foldstats_data_var_set
    __swig_getmethods__["data_var"] = _presto.foldstats_data_var_get
    if _newclass:data_var = property(_presto.foldstats_data_var_get,_presto.foldstats_data_var_set)
    __swig_setmethods__["numprof"] = _presto.foldstats_numprof_set
    __swig_getmethods__["numprof"] = _presto.foldstats_numprof_get
    if _newclass:numprof = property(_presto.foldstats_numprof_get,_presto.foldstats_numprof_set)
    __swig_setmethods__["prof_avg"] = _presto.foldstats_prof_avg_set
    __swig_getmethods__["prof_avg"] = _presto.foldstats_prof_avg_get
    if _newclass:prof_avg = property(_presto.foldstats_prof_avg_get,_presto.foldstats_prof_avg_set)
    __swig_setmethods__["prof_var"] = _presto.foldstats_prof_var_set
    __swig_getmethods__["prof_var"] = _presto.foldstats_prof_var_get
    if _newclass:prof_var = property(_presto.foldstats_prof_var_get,_presto.foldstats_prof_var_set)
    __swig_setmethods__["redchi"] = _presto.foldstats_redchi_set
    __swig_getmethods__["redchi"] = _presto.foldstats_redchi_get
    if _newclass:redchi = property(_presto.foldstats_redchi_get,_presto.foldstats_redchi_set)
    def __init__(self,*args):
        self.this = apply(_presto.new_foldstats,args)
        self.thisown = 1
    def __del__(self, destroy= _presto.delete_foldstats):
        try:
            if self.thisown: destroy(self)
        except: pass
    def __repr__(self):
        return "<C foldstats instance at %s>" % (self.this,)

class foldstatsPtr(foldstats):
    def __init__(self,this):
        self.this = this
        if not hasattr(self,"thisown"): self.thisown = 0
        self.__class__ = foldstats
_presto.foldstats_swigregister(foldstatsPtr)

read_int = _presto.read_int

read_double = _presto.read_double

frotate = _presto.frotate

drotate = _presto.drotate

dorbint = _presto.dorbint

keplars_eqn = _presto.keplars_eqn

lin_interp_E = _presto.lin_interp_E

E_to_phib = _presto.E_to_phib

E_to_v = _presto.E_to_v

E_to_p = _presto.E_to_p

E_to_z = _presto.E_to_z

E_to_phib_BT = _presto.E_to_phib_BT

r_resp_halfwidth = _presto.r_resp_halfwidth

z_resp_halfwidth = _presto.z_resp_halfwidth

w_resp_halfwidth = _presto.w_resp_halfwidth

binary_velocity = _presto.binary_velocity

bin_resp_halfwidth = _presto.bin_resp_halfwidth

gen_r_response = _presto.gen_r_response

gen_z_response = _presto.gen_z_response

gen_w_response = _presto.gen_w_response

gen_bin_response = _presto.gen_bin_response

get_numphotons = _presto.get_numphotons

get_localpower = _presto.get_localpower

get_localpower3d = _presto.get_localpower3d

get_derivs3d = _presto.get_derivs3d

calc_props = _presto.calc_props

calc_binprops = _presto.calc_binprops

calc_rzwerrs = _presto.calc_rzwerrs

candidate_sigma = _presto.candidate_sigma

power_for_sigma = _presto.power_for_sigma

chisqr = _presto.chisqr

print_candidate = _presto.print_candidate

print_bin_candidate = _presto.print_bin_candidate

read_rzw_cand = _presto.read_rzw_cand

read_bin_cand = _presto.read_bin_cand

read_rawbin_cand = _presto.read_rawbin_cand

get_rzw_cand = _presto.get_rzw_cand

get_bin_cand = _presto.get_bin_cand

get_rawbin_cand = _presto.get_rawbin_cand

chkfilelen = _presto.chkfilelen

read_fcomplex_file = _presto.read_fcomplex_file

read_float_file = _presto.read_float_file

prune_powers = _presto.prune_powers

median = _presto.median

dms2rad = _presto.dms2rad

hms2rad = _presto.hms2rad

sphere_ang_diff = _presto.sphere_ang_diff

spread_with_pad = _presto.spread_with_pad

spread_no_pad = _presto.spread_no_pad

paddata = _presto.paddata

place_complex_kernel = _presto.place_complex_kernel

place_real_kernel = _presto.place_real_kernel

chop_complex_ends = _presto.chop_complex_ends

chop_real_ends = _presto.chop_real_ends

complex_corr_conv = _presto.complex_corr_conv

real_corr_conv = _presto.real_corr_conv

corr_complex = _presto.corr_complex

stretch_fft = _presto.stretch_fft

corr_loc_pow = _presto.corr_loc_pow

rz_interp = _presto.rz_interp

max_r_arr = _presto.max_r_arr

max_rz_arr = _presto.max_rz_arr

fold_errors = _presto.fold_errors

foldfile = _presto.foldfile

fold = _presto.fold

simplefold = _presto.simplefold

doppler = _presto.doppler

search_minifft = _presto.search_minifft

print_rawbincand = _presto.print_rawbincand

barycenter = _presto.barycenter

fftwcall = _presto.fftwcall

tablesixstepfft = _presto.tablesixstepfft

realfft = _presto.realfft

tree_max_dm = _presto.tree_max_dm

smearing_from_bw = _presto.smearing_from_bw

delay_from_dm = _presto.delay_from_dm

dm_from_delay = _presto.dm_from_delay

dedisp_delays = _presto.dedisp_delays

subband_search_delays = _presto.subband_search_delays

subband_delays = _presto.subband_delays

nice_output_1 = _presto.nice_output_1

nice_output_2 = _presto.nice_output_2

corr_rz_plane = _presto.corr_rz_plane

corr_rz_interp = _presto.corr_rz_interp




#-------------- Extra Stuff to Make Things Easier -----------------

import math, umath, Numeric, Pgplot, string, numpyio, miscutils

def read_foldstats(file, byteswap=0):
   stats = foldstats()
   stats.numdata = read_double(file, byteswap)
   stats.data_avg = read_double(file, byteswap)
   stats.data_var = read_double(file, byteswap)
   stats.numprof = read_double(file, byteswap)
   stats.prof_avg = read_double(file, byteswap)
   stats.prof_var = read_double(file, byteswap)
   stats.redchi = read_double(file, byteswap)
   return stats
  
class pfd:
   def __init__(self, filename):
      infile = open(filename, "rb")
      byteswap = 0
      testswap = read_double(infile, byteswap)
      if (testswap < 1.0 or testswap > 10000.0):
         byteswap = 1
      infile.seek(0)
      self.npart = int(read_double(infile, byteswap))
      self.nsub = int(read_double(infile, byteswap))
      self.proflen = int(read_double(infile, byteswap))
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
               self.stats[ii].append(read_foldstats(infile, byteswap))
               self.profs[ii][jj] = self.profs[ii][jj] + \
                                    numpyio.fread(infile,
                                                  self.proflen, 'd',
                                                  'd', byteswap)
         else:
            self.stats.append(read_foldstats(infile, byteswap))
            self.profs[ii] = self.profs[ii]+ \
                             numpyio.fread(infile, self.proflen, 'd',
                                           'd', byteswap)
      infile.close()
   
def val_with_err(value, error, len=0, digits=2, latex=0):
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
          If latex=1, the string is converted into LaTeX markup.
   """
   slen = 40
   outstr = ' '*slen
   if abs(len) > slen:
      slen = abs(len)
   if digits==2:
      slen = nice_output_2(outstr, value, error, len)
   else:
      slen = nice_output_1(outstr, value, error, len)
   if len <= 0 and outstr[0]==' ':  # Not quite sure why this is necessary...
      outstr = outstr[1:slen]
      if len < 0:
         outstr += ' '
   outstr = outstr[:slen]
   if latex:
      if outstr.find("x10") > 0:
         outstr = outstr.replace("x10^", "$\times$10$^{")+"}$"
   return outstr

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
