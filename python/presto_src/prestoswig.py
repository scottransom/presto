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
del types


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
    def __repr__(self):
        return "<C orbitparams instance at %s>" % (self.this,)
    __swig_setmethods__["p"] = _presto.orbitparams_p_set
    __swig_getmethods__["p"] = _presto.orbitparams_p_get
    if _newclass:p = property(_presto.orbitparams_p_get, _presto.orbitparams_p_set)
    __swig_setmethods__["e"] = _presto.orbitparams_e_set
    __swig_getmethods__["e"] = _presto.orbitparams_e_get
    if _newclass:e = property(_presto.orbitparams_e_get, _presto.orbitparams_e_set)
    __swig_setmethods__["x"] = _presto.orbitparams_x_set
    __swig_getmethods__["x"] = _presto.orbitparams_x_get
    if _newclass:x = property(_presto.orbitparams_x_get, _presto.orbitparams_x_set)
    __swig_setmethods__["w"] = _presto.orbitparams_w_set
    __swig_getmethods__["w"] = _presto.orbitparams_w_get
    if _newclass:w = property(_presto.orbitparams_w_get, _presto.orbitparams_w_set)
    __swig_setmethods__["t"] = _presto.orbitparams_t_set
    __swig_getmethods__["t"] = _presto.orbitparams_t_get
    if _newclass:t = property(_presto.orbitparams_t_get, _presto.orbitparams_t_set)
    __swig_setmethods__["pd"] = _presto.orbitparams_pd_set
    __swig_getmethods__["pd"] = _presto.orbitparams_pd_get
    if _newclass:pd = property(_presto.orbitparams_pd_get, _presto.orbitparams_pd_set)
    __swig_setmethods__["wd"] = _presto.orbitparams_wd_set
    __swig_getmethods__["wd"] = _presto.orbitparams_wd_get
    if _newclass:wd = property(_presto.orbitparams_wd_get, _presto.orbitparams_wd_set)
    def __init__(self, *args):
        _swig_setattr(self, orbitparams, 'this', _presto.new_orbitparams(*args))
        _swig_setattr(self, orbitparams, 'thisown', 1)
    def __del__(self, destroy=_presto.delete_orbitparams):
        try:
            if self.thisown: destroy(self)
        except: pass

class orbitparamsPtr(orbitparams):
    def __init__(self, this):
        _swig_setattr(self, orbitparams, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, orbitparams, 'thisown', 0)
        _swig_setattr(self, orbitparams,self.__class__,orbitparams)
_presto.orbitparams_swigregister(orbitparamsPtr)
tofloatvector = _presto.tofloatvector

float_to_complex = _presto.float_to_complex

complex_to_float = _presto.complex_to_float


class DoubleArray(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, DoubleArray, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, DoubleArray, name)
    def __repr__(self):
        return "<C DoubleArray instance at %s>" % (self.this,)
    __swig_setmethods__["dptr"] = _presto.DoubleArray_dptr_set
    __swig_getmethods__["dptr"] = _presto.DoubleArray_dptr_get
    if _newclass:dptr = property(_presto.DoubleArray_dptr_get, _presto.DoubleArray_dptr_set)
    def __init__(self, *args):
        _swig_setattr(self, DoubleArray, 'this', _presto.new_DoubleArray(*args))
        _swig_setattr(self, DoubleArray, 'thisown', 1)
    def __del__(self, destroy=_presto.delete_DoubleArray):
        try:
            if self.thisown: destroy(self)
        except: pass
    def __getitem__(*args): return _presto.DoubleArray___getitem__(*args)
    def __setitem__(*args): return _presto.DoubleArray___setitem__(*args)

class DoubleArrayPtr(DoubleArray):
    def __init__(self, this):
        _swig_setattr(self, DoubleArray, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, DoubleArray, 'thisown', 0)
        _swig_setattr(self, DoubleArray,self.__class__,DoubleArray)
_presto.DoubleArray_swigregister(DoubleArrayPtr)

class infodata(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, infodata, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, infodata, name)
    def __repr__(self):
        return "<C infodata instance at %s>" % (self.this,)
    __swig_setmethods__["ra_s"] = _presto.infodata_ra_s_set
    __swig_getmethods__["ra_s"] = _presto.infodata_ra_s_get
    if _newclass:ra_s = property(_presto.infodata_ra_s_get, _presto.infodata_ra_s_set)
    __swig_setmethods__["dec_s"] = _presto.infodata_dec_s_set
    __swig_getmethods__["dec_s"] = _presto.infodata_dec_s_get
    if _newclass:dec_s = property(_presto.infodata_dec_s_get, _presto.infodata_dec_s_set)
    __swig_setmethods__["N"] = _presto.infodata_N_set
    __swig_getmethods__["N"] = _presto.infodata_N_get
    if _newclass:N = property(_presto.infodata_N_get, _presto.infodata_N_set)
    __swig_setmethods__["dt"] = _presto.infodata_dt_set
    __swig_getmethods__["dt"] = _presto.infodata_dt_get
    if _newclass:dt = property(_presto.infodata_dt_get, _presto.infodata_dt_set)
    __swig_setmethods__["fov"] = _presto.infodata_fov_set
    __swig_getmethods__["fov"] = _presto.infodata_fov_get
    if _newclass:fov = property(_presto.infodata_fov_get, _presto.infodata_fov_set)
    __swig_setmethods__["mjd_f"] = _presto.infodata_mjd_f_set
    __swig_getmethods__["mjd_f"] = _presto.infodata_mjd_f_get
    if _newclass:mjd_f = property(_presto.infodata_mjd_f_get, _presto.infodata_mjd_f_set)
    __swig_setmethods__["dm"] = _presto.infodata_dm_set
    __swig_getmethods__["dm"] = _presto.infodata_dm_get
    if _newclass:dm = property(_presto.infodata_dm_get, _presto.infodata_dm_set)
    __swig_setmethods__["freq"] = _presto.infodata_freq_set
    __swig_getmethods__["freq"] = _presto.infodata_freq_get
    if _newclass:freq = property(_presto.infodata_freq_get, _presto.infodata_freq_set)
    __swig_setmethods__["freqband"] = _presto.infodata_freqband_set
    __swig_getmethods__["freqband"] = _presto.infodata_freqband_get
    if _newclass:freqband = property(_presto.infodata_freqband_get, _presto.infodata_freqband_set)
    __swig_setmethods__["chan_wid"] = _presto.infodata_chan_wid_set
    __swig_getmethods__["chan_wid"] = _presto.infodata_chan_wid_get
    if _newclass:chan_wid = property(_presto.infodata_chan_wid_get, _presto.infodata_chan_wid_set)
    __swig_setmethods__["wavelen"] = _presto.infodata_wavelen_set
    __swig_getmethods__["wavelen"] = _presto.infodata_wavelen_get
    if _newclass:wavelen = property(_presto.infodata_wavelen_get, _presto.infodata_wavelen_set)
    __swig_setmethods__["waveband"] = _presto.infodata_waveband_set
    __swig_getmethods__["waveband"] = _presto.infodata_waveband_get
    if _newclass:waveband = property(_presto.infodata_waveband_get, _presto.infodata_waveband_set)
    __swig_setmethods__["energy"] = _presto.infodata_energy_set
    __swig_getmethods__["energy"] = _presto.infodata_energy_get
    if _newclass:energy = property(_presto.infodata_energy_get, _presto.infodata_energy_set)
    __swig_setmethods__["energyband"] = _presto.infodata_energyband_set
    __swig_getmethods__["energyband"] = _presto.infodata_energyband_get
    if _newclass:energyband = property(_presto.infodata_energyband_get, _presto.infodata_energyband_set)
    __swig_setmethods__["onoff"] = _presto.infodata_onoff_set
    __swig_getmethods__["onoff"] = _presto.infodata_onoff_get
    if _newclass:onoff = property(_presto.infodata_onoff_get, _presto.infodata_onoff_set)
    __swig_setmethods__["num_chan"] = _presto.infodata_num_chan_set
    __swig_getmethods__["num_chan"] = _presto.infodata_num_chan_get
    if _newclass:num_chan = property(_presto.infodata_num_chan_get, _presto.infodata_num_chan_set)
    __swig_setmethods__["mjd_i"] = _presto.infodata_mjd_i_set
    __swig_getmethods__["mjd_i"] = _presto.infodata_mjd_i_get
    if _newclass:mjd_i = property(_presto.infodata_mjd_i_get, _presto.infodata_mjd_i_set)
    __swig_setmethods__["ra_h"] = _presto.infodata_ra_h_set
    __swig_getmethods__["ra_h"] = _presto.infodata_ra_h_get
    if _newclass:ra_h = property(_presto.infodata_ra_h_get, _presto.infodata_ra_h_set)
    __swig_setmethods__["ra_m"] = _presto.infodata_ra_m_set
    __swig_getmethods__["ra_m"] = _presto.infodata_ra_m_get
    if _newclass:ra_m = property(_presto.infodata_ra_m_get, _presto.infodata_ra_m_set)
    __swig_setmethods__["dec_d"] = _presto.infodata_dec_d_set
    __swig_getmethods__["dec_d"] = _presto.infodata_dec_d_get
    if _newclass:dec_d = property(_presto.infodata_dec_d_get, _presto.infodata_dec_d_set)
    __swig_setmethods__["dec_m"] = _presto.infodata_dec_m_set
    __swig_getmethods__["dec_m"] = _presto.infodata_dec_m_get
    if _newclass:dec_m = property(_presto.infodata_dec_m_get, _presto.infodata_dec_m_set)
    __swig_setmethods__["bary"] = _presto.infodata_bary_set
    __swig_getmethods__["bary"] = _presto.infodata_bary_get
    if _newclass:bary = property(_presto.infodata_bary_get, _presto.infodata_bary_set)
    __swig_setmethods__["numonoff"] = _presto.infodata_numonoff_set
    __swig_getmethods__["numonoff"] = _presto.infodata_numonoff_get
    if _newclass:numonoff = property(_presto.infodata_numonoff_get, _presto.infodata_numonoff_set)
    __swig_setmethods__["notes"] = _presto.infodata_notes_set
    __swig_getmethods__["notes"] = _presto.infodata_notes_get
    if _newclass:notes = property(_presto.infodata_notes_get, _presto.infodata_notes_set)
    __swig_setmethods__["name"] = _presto.infodata_name_set
    __swig_getmethods__["name"] = _presto.infodata_name_get
    if _newclass:name = property(_presto.infodata_name_get, _presto.infodata_name_set)
    __swig_setmethods__["object"] = _presto.infodata_object_set
    __swig_getmethods__["object"] = _presto.infodata_object_get
    if _newclass:object = property(_presto.infodata_object_get, _presto.infodata_object_set)
    __swig_setmethods__["instrument"] = _presto.infodata_instrument_set
    __swig_getmethods__["instrument"] = _presto.infodata_instrument_get
    if _newclass:instrument = property(_presto.infodata_instrument_get, _presto.infodata_instrument_set)
    __swig_setmethods__["observer"] = _presto.infodata_observer_set
    __swig_getmethods__["observer"] = _presto.infodata_observer_get
    if _newclass:observer = property(_presto.infodata_observer_get, _presto.infodata_observer_set)
    __swig_setmethods__["analyzer"] = _presto.infodata_analyzer_set
    __swig_getmethods__["analyzer"] = _presto.infodata_analyzer_get
    if _newclass:analyzer = property(_presto.infodata_analyzer_get, _presto.infodata_analyzer_set)
    __swig_setmethods__["telescope"] = _presto.infodata_telescope_set
    __swig_getmethods__["telescope"] = _presto.infodata_telescope_get
    if _newclass:telescope = property(_presto.infodata_telescope_get, _presto.infodata_telescope_set)
    __swig_setmethods__["band"] = _presto.infodata_band_set
    __swig_getmethods__["band"] = _presto.infodata_band_get
    if _newclass:band = property(_presto.infodata_band_get, _presto.infodata_band_set)
    __swig_setmethods__["filt"] = _presto.infodata_filt_set
    __swig_getmethods__["filt"] = _presto.infodata_filt_get
    if _newclass:filt = property(_presto.infodata_filt_get, _presto.infodata_filt_set)
    def __init__(self, *args):
        _swig_setattr(self, infodata, 'this', _presto.new_infodata(*args))
        _swig_setattr(self, infodata, 'thisown', 1)
    def __del__(self, destroy=_presto.delete_infodata):
        try:
            if self.thisown: destroy(self)
        except: pass

class infodataPtr(infodata):
    def __init__(self, this):
        _swig_setattr(self, infodata, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, infodata, 'thisown', 0)
        _swig_setattr(self, infodata,self.__class__,infodata)
_presto.infodata_swigregister(infodataPtr)


readinf = _presto.readinf

writeinf = _presto.writeinf
class makedata(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, makedata, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, makedata, name)
    def __repr__(self):
        return "<C makedata instance at %s>" % (self.this,)
    __swig_setmethods__["basefilenm"] = _presto.makedata_basefilenm_set
    __swig_getmethods__["basefilenm"] = _presto.makedata_basefilenm_get
    if _newclass:basefilenm = property(_presto.makedata_basefilenm_get, _presto.makedata_basefilenm_set)
    __swig_setmethods__["description"] = _presto.makedata_description_set
    __swig_getmethods__["description"] = _presto.makedata_description_get
    if _newclass:description = property(_presto.makedata_description_get, _presto.makedata_description_set)
    __swig_setmethods__["N"] = _presto.makedata_N_set
    __swig_getmethods__["N"] = _presto.makedata_N_get
    if _newclass:N = property(_presto.makedata_N_get, _presto.makedata_N_set)
    __swig_setmethods__["next2_to_n"] = _presto.makedata_next2_to_n_set
    __swig_getmethods__["next2_to_n"] = _presto.makedata_next2_to_n_get
    if _newclass:next2_to_n = property(_presto.makedata_next2_to_n_get, _presto.makedata_next2_to_n_set)
    __swig_setmethods__["dt"] = _presto.makedata_dt_set
    __swig_getmethods__["dt"] = _presto.makedata_dt_get
    if _newclass:dt = property(_presto.makedata_dt_get, _presto.makedata_dt_set)
    __swig_setmethods__["T"] = _presto.makedata_T_set
    __swig_getmethods__["T"] = _presto.makedata_T_get
    if _newclass:T = property(_presto.makedata_T_get, _presto.makedata_T_set)
    __swig_setmethods__["ptype"] = _presto.makedata_ptype_set
    __swig_getmethods__["ptype"] = _presto.makedata_ptype_get
    if _newclass:ptype = property(_presto.makedata_ptype_get, _presto.makedata_ptype_set)
    __swig_setmethods__["pnum"] = _presto.makedata_pnum_set
    __swig_getmethods__["pnum"] = _presto.makedata_pnum_get
    if _newclass:pnum = property(_presto.makedata_pnum_get, _presto.makedata_pnum_set)
    __swig_setmethods__["fwhm"] = _presto.makedata_fwhm_set
    __swig_getmethods__["fwhm"] = _presto.makedata_fwhm_get
    if _newclass:fwhm = property(_presto.makedata_fwhm_get, _presto.makedata_fwhm_set)
    __swig_setmethods__["round"] = _presto.makedata_round_set
    __swig_getmethods__["round"] = _presto.makedata_round_get
    if _newclass:round = property(_presto.makedata_round_get, _presto.makedata_round_set)
    __swig_setmethods__["roundnum"] = _presto.makedata_roundnum_set
    __swig_getmethods__["roundnum"] = _presto.makedata_roundnum_get
    if _newclass:roundnum = property(_presto.makedata_roundnum_get, _presto.makedata_roundnum_set)
    __swig_setmethods__["f"] = _presto.makedata_f_set
    __swig_getmethods__["f"] = _presto.makedata_f_get
    if _newclass:f = property(_presto.makedata_f_get, _presto.makedata_f_set)
    __swig_setmethods__["fd"] = _presto.makedata_fd_set
    __swig_getmethods__["fd"] = _presto.makedata_fd_get
    if _newclass:fd = property(_presto.makedata_fd_get, _presto.makedata_fd_set)
    __swig_setmethods__["fdd"] = _presto.makedata_fdd_set
    __swig_getmethods__["fdd"] = _presto.makedata_fdd_get
    if _newclass:fdd = property(_presto.makedata_fdd_get, _presto.makedata_fdd_set)
    __swig_setmethods__["p"] = _presto.makedata_p_set
    __swig_getmethods__["p"] = _presto.makedata_p_get
    if _newclass:p = property(_presto.makedata_p_get, _presto.makedata_p_set)
    __swig_setmethods__["pd"] = _presto.makedata_pd_set
    __swig_getmethods__["pd"] = _presto.makedata_pd_get
    if _newclass:pd = property(_presto.makedata_pd_get, _presto.makedata_pd_set)
    __swig_setmethods__["pdd"] = _presto.makedata_pdd_set
    __swig_getmethods__["pdd"] = _presto.makedata_pdd_get
    if _newclass:pdd = property(_presto.makedata_pdd_get, _presto.makedata_pdd_set)
    __swig_setmethods__["r"] = _presto.makedata_r_set
    __swig_getmethods__["r"] = _presto.makedata_r_get
    if _newclass:r = property(_presto.makedata_r_get, _presto.makedata_r_set)
    __swig_setmethods__["z"] = _presto.makedata_z_set
    __swig_getmethods__["z"] = _presto.makedata_z_get
    if _newclass:z = property(_presto.makedata_z_get, _presto.makedata_z_set)
    __swig_setmethods__["w"] = _presto.makedata_w_set
    __swig_getmethods__["w"] = _presto.makedata_w_get
    if _newclass:w = property(_presto.makedata_w_get, _presto.makedata_w_set)
    __swig_setmethods__["amp"] = _presto.makedata_amp_set
    __swig_getmethods__["amp"] = _presto.makedata_amp_get
    if _newclass:amp = property(_presto.makedata_amp_get, _presto.makedata_amp_set)
    __swig_setmethods__["phs"] = _presto.makedata_phs_set
    __swig_getmethods__["phs"] = _presto.makedata_phs_get
    if _newclass:phs = property(_presto.makedata_phs_get, _presto.makedata_phs_set)
    __swig_setmethods__["dc"] = _presto.makedata_dc_set
    __swig_getmethods__["dc"] = _presto.makedata_dc_get
    if _newclass:dc = property(_presto.makedata_dc_get, _presto.makedata_dc_set)
    __swig_setmethods__["binary"] = _presto.makedata_binary_set
    __swig_getmethods__["binary"] = _presto.makedata_binary_get
    if _newclass:binary = property(_presto.makedata_binary_get, _presto.makedata_binary_set)
    __swig_setmethods__["orb"] = _presto.makedata_orb_set
    __swig_getmethods__["orb"] = _presto.makedata_orb_get
    if _newclass:orb = property(_presto.makedata_orb_get, _presto.makedata_orb_set)
    __swig_setmethods__["ampmod"] = _presto.makedata_ampmod_set
    __swig_getmethods__["ampmod"] = _presto.makedata_ampmod_get
    if _newclass:ampmod = property(_presto.makedata_ampmod_get, _presto.makedata_ampmod_set)
    __swig_setmethods__["ampmoda"] = _presto.makedata_ampmoda_set
    __swig_getmethods__["ampmoda"] = _presto.makedata_ampmoda_get
    if _newclass:ampmoda = property(_presto.makedata_ampmoda_get, _presto.makedata_ampmoda_set)
    __swig_setmethods__["ampmodf"] = _presto.makedata_ampmodf_set
    __swig_getmethods__["ampmodf"] = _presto.makedata_ampmodf_get
    if _newclass:ampmodf = property(_presto.makedata_ampmodf_get, _presto.makedata_ampmodf_set)
    __swig_setmethods__["ampmodp"] = _presto.makedata_ampmodp_set
    __swig_getmethods__["ampmodp"] = _presto.makedata_ampmodp_get
    if _newclass:ampmodp = property(_presto.makedata_ampmodp_get, _presto.makedata_ampmodp_set)
    __swig_setmethods__["noisetype"] = _presto.makedata_noisetype_set
    __swig_getmethods__["noisetype"] = _presto.makedata_noisetype_get
    if _newclass:noisetype = property(_presto.makedata_noisetype_get, _presto.makedata_noisetype_set)
    __swig_setmethods__["noise"] = _presto.makedata_noise_set
    __swig_getmethods__["noise"] = _presto.makedata_noise_get
    if _newclass:noise = property(_presto.makedata_noise_get, _presto.makedata_noise_set)
    __swig_setmethods__["noisesig"] = _presto.makedata_noisesig_set
    __swig_getmethods__["noisesig"] = _presto.makedata_noisesig_get
    if _newclass:noisesig = property(_presto.makedata_noisesig_get, _presto.makedata_noisesig_set)
    __swig_setmethods__["numonoff"] = _presto.makedata_numonoff_set
    __swig_getmethods__["numonoff"] = _presto.makedata_numonoff_get
    if _newclass:numonoff = property(_presto.makedata_numonoff_get, _presto.makedata_numonoff_set)
    __swig_setmethods__["onoff"] = _presto.makedata_onoff_set
    __swig_getmethods__["onoff"] = _presto.makedata_onoff_get
    if _newclass:onoff = property(_presto.makedata_onoff_get, _presto.makedata_onoff_set)
    def __init__(self, *args):
        _swig_setattr(self, makedata, 'this', _presto.new_makedata(*args))
        _swig_setattr(self, makedata, 'thisown', 1)
    def __del__(self, destroy=_presto.delete_makedata):
        try:
            if self.thisown: destroy(self)
        except: pass

class makedataPtr(makedata):
    def __init__(self, this):
        _swig_setattr(self, makedata, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, makedata, 'thisown', 0)
        _swig_setattr(self, makedata,self.__class__,makedata)
_presto.makedata_swigregister(makedataPtr)

class psrparams(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, psrparams, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, psrparams, name)
    def __repr__(self):
        return "<C psrparams instance at %s>" % (self.this,)
    __swig_setmethods__["jname"] = _presto.psrparams_jname_set
    __swig_getmethods__["jname"] = _presto.psrparams_jname_get
    if _newclass:jname = property(_presto.psrparams_jname_get, _presto.psrparams_jname_set)
    __swig_setmethods__["bname"] = _presto.psrparams_bname_set
    __swig_getmethods__["bname"] = _presto.psrparams_bname_get
    if _newclass:bname = property(_presto.psrparams_bname_get, _presto.psrparams_bname_set)
    __swig_setmethods__["alias"] = _presto.psrparams_alias_set
    __swig_getmethods__["alias"] = _presto.psrparams_alias_get
    if _newclass:alias = property(_presto.psrparams_alias_get, _presto.psrparams_alias_set)
    __swig_setmethods__["ra2000"] = _presto.psrparams_ra2000_set
    __swig_getmethods__["ra2000"] = _presto.psrparams_ra2000_get
    if _newclass:ra2000 = property(_presto.psrparams_ra2000_get, _presto.psrparams_ra2000_set)
    __swig_setmethods__["dec2000"] = _presto.psrparams_dec2000_set
    __swig_getmethods__["dec2000"] = _presto.psrparams_dec2000_get
    if _newclass:dec2000 = property(_presto.psrparams_dec2000_get, _presto.psrparams_dec2000_set)
    __swig_setmethods__["dm"] = _presto.psrparams_dm_set
    __swig_getmethods__["dm"] = _presto.psrparams_dm_get
    if _newclass:dm = property(_presto.psrparams_dm_get, _presto.psrparams_dm_set)
    __swig_setmethods__["timepoch"] = _presto.psrparams_timepoch_set
    __swig_getmethods__["timepoch"] = _presto.psrparams_timepoch_get
    if _newclass:timepoch = property(_presto.psrparams_timepoch_get, _presto.psrparams_timepoch_set)
    __swig_setmethods__["p"] = _presto.psrparams_p_set
    __swig_getmethods__["p"] = _presto.psrparams_p_get
    if _newclass:p = property(_presto.psrparams_p_get, _presto.psrparams_p_set)
    __swig_setmethods__["pd"] = _presto.psrparams_pd_set
    __swig_getmethods__["pd"] = _presto.psrparams_pd_get
    if _newclass:pd = property(_presto.psrparams_pd_get, _presto.psrparams_pd_set)
    __swig_setmethods__["pdd"] = _presto.psrparams_pdd_set
    __swig_getmethods__["pdd"] = _presto.psrparams_pdd_get
    if _newclass:pdd = property(_presto.psrparams_pdd_get, _presto.psrparams_pdd_set)
    __swig_setmethods__["f"] = _presto.psrparams_f_set
    __swig_getmethods__["f"] = _presto.psrparams_f_get
    if _newclass:f = property(_presto.psrparams_f_get, _presto.psrparams_f_set)
    __swig_setmethods__["fd"] = _presto.psrparams_fd_set
    __swig_getmethods__["fd"] = _presto.psrparams_fd_get
    if _newclass:fd = property(_presto.psrparams_fd_get, _presto.psrparams_fd_set)
    __swig_setmethods__["fdd"] = _presto.psrparams_fdd_set
    __swig_getmethods__["fdd"] = _presto.psrparams_fdd_get
    if _newclass:fdd = property(_presto.psrparams_fdd_get, _presto.psrparams_fdd_set)
    __swig_setmethods__["orb"] = _presto.psrparams_orb_set
    __swig_getmethods__["orb"] = _presto.psrparams_orb_get
    if _newclass:orb = property(_presto.psrparams_orb_get, _presto.psrparams_orb_set)
    def __init__(self, *args):
        _swig_setattr(self, psrparams, 'this', _presto.new_psrparams(*args))
        _swig_setattr(self, psrparams, 'thisown', 1)
    def __del__(self, destroy=_presto.delete_psrparams):
        try:
            if self.thisown: destroy(self)
        except: pass

class psrparamsPtr(psrparams):
    def __init__(self, this):
        _swig_setattr(self, psrparams, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, psrparams, 'thisown', 0)
        _swig_setattr(self, psrparams,self.__class__,psrparams)
_presto.psrparams_swigregister(psrparamsPtr)


read_mak_input = _presto.read_mak_input

read_mak_file = _presto.read_mak_file

write_mak_file = _presto.write_mak_file
class rderivs(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, rderivs, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, rderivs, name)
    def __repr__(self):
        return "<C rderivs instance at %s>" % (self.this,)
    __swig_setmethods__["pow"] = _presto.rderivs_pow_set
    __swig_getmethods__["pow"] = _presto.rderivs_pow_get
    if _newclass:pow = property(_presto.rderivs_pow_get, _presto.rderivs_pow_set)
    __swig_setmethods__["phs"] = _presto.rderivs_phs_set
    __swig_getmethods__["phs"] = _presto.rderivs_phs_get
    if _newclass:phs = property(_presto.rderivs_phs_get, _presto.rderivs_phs_set)
    __swig_setmethods__["dpow"] = _presto.rderivs_dpow_set
    __swig_getmethods__["dpow"] = _presto.rderivs_dpow_get
    if _newclass:dpow = property(_presto.rderivs_dpow_get, _presto.rderivs_dpow_set)
    __swig_setmethods__["dphs"] = _presto.rderivs_dphs_set
    __swig_getmethods__["dphs"] = _presto.rderivs_dphs_get
    if _newclass:dphs = property(_presto.rderivs_dphs_get, _presto.rderivs_dphs_set)
    __swig_setmethods__["d2pow"] = _presto.rderivs_d2pow_set
    __swig_getmethods__["d2pow"] = _presto.rderivs_d2pow_get
    if _newclass:d2pow = property(_presto.rderivs_d2pow_get, _presto.rderivs_d2pow_set)
    __swig_setmethods__["d2phs"] = _presto.rderivs_d2phs_set
    __swig_getmethods__["d2phs"] = _presto.rderivs_d2phs_get
    if _newclass:d2phs = property(_presto.rderivs_d2phs_get, _presto.rderivs_d2phs_set)
    __swig_setmethods__["locpow"] = _presto.rderivs_locpow_set
    __swig_getmethods__["locpow"] = _presto.rderivs_locpow_get
    if _newclass:locpow = property(_presto.rderivs_locpow_get, _presto.rderivs_locpow_set)
    def __init__(self, *args):
        _swig_setattr(self, rderivs, 'this', _presto.new_rderivs(*args))
        _swig_setattr(self, rderivs, 'thisown', 1)
    def __del__(self, destroy=_presto.delete_rderivs):
        try:
            if self.thisown: destroy(self)
        except: pass

class rderivsPtr(rderivs):
    def __init__(self, this):
        _swig_setattr(self, rderivs, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, rderivs, 'thisown', 0)
        _swig_setattr(self, rderivs,self.__class__,rderivs)
_presto.rderivs_swigregister(rderivsPtr)

class fourierprops(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, fourierprops, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, fourierprops, name)
    def __repr__(self):
        return "<C fourierprops instance at %s>" % (self.this,)
    __swig_setmethods__["r"] = _presto.fourierprops_r_set
    __swig_getmethods__["r"] = _presto.fourierprops_r_get
    if _newclass:r = property(_presto.fourierprops_r_get, _presto.fourierprops_r_set)
    __swig_setmethods__["rerr"] = _presto.fourierprops_rerr_set
    __swig_getmethods__["rerr"] = _presto.fourierprops_rerr_get
    if _newclass:rerr = property(_presto.fourierprops_rerr_get, _presto.fourierprops_rerr_set)
    __swig_setmethods__["z"] = _presto.fourierprops_z_set
    __swig_getmethods__["z"] = _presto.fourierprops_z_get
    if _newclass:z = property(_presto.fourierprops_z_get, _presto.fourierprops_z_set)
    __swig_setmethods__["zerr"] = _presto.fourierprops_zerr_set
    __swig_getmethods__["zerr"] = _presto.fourierprops_zerr_get
    if _newclass:zerr = property(_presto.fourierprops_zerr_get, _presto.fourierprops_zerr_set)
    __swig_setmethods__["w"] = _presto.fourierprops_w_set
    __swig_getmethods__["w"] = _presto.fourierprops_w_get
    if _newclass:w = property(_presto.fourierprops_w_get, _presto.fourierprops_w_set)
    __swig_setmethods__["werr"] = _presto.fourierprops_werr_set
    __swig_getmethods__["werr"] = _presto.fourierprops_werr_get
    if _newclass:werr = property(_presto.fourierprops_werr_get, _presto.fourierprops_werr_set)
    __swig_setmethods__["pow"] = _presto.fourierprops_pow_set
    __swig_getmethods__["pow"] = _presto.fourierprops_pow_get
    if _newclass:pow = property(_presto.fourierprops_pow_get, _presto.fourierprops_pow_set)
    __swig_setmethods__["powerr"] = _presto.fourierprops_powerr_set
    __swig_getmethods__["powerr"] = _presto.fourierprops_powerr_get
    if _newclass:powerr = property(_presto.fourierprops_powerr_get, _presto.fourierprops_powerr_set)
    __swig_setmethods__["sig"] = _presto.fourierprops_sig_set
    __swig_getmethods__["sig"] = _presto.fourierprops_sig_get
    if _newclass:sig = property(_presto.fourierprops_sig_get, _presto.fourierprops_sig_set)
    __swig_setmethods__["rawpow"] = _presto.fourierprops_rawpow_set
    __swig_getmethods__["rawpow"] = _presto.fourierprops_rawpow_get
    if _newclass:rawpow = property(_presto.fourierprops_rawpow_get, _presto.fourierprops_rawpow_set)
    __swig_setmethods__["phs"] = _presto.fourierprops_phs_set
    __swig_getmethods__["phs"] = _presto.fourierprops_phs_get
    if _newclass:phs = property(_presto.fourierprops_phs_get, _presto.fourierprops_phs_set)
    __swig_setmethods__["phserr"] = _presto.fourierprops_phserr_set
    __swig_getmethods__["phserr"] = _presto.fourierprops_phserr_get
    if _newclass:phserr = property(_presto.fourierprops_phserr_get, _presto.fourierprops_phserr_set)
    __swig_setmethods__["cen"] = _presto.fourierprops_cen_set
    __swig_getmethods__["cen"] = _presto.fourierprops_cen_get
    if _newclass:cen = property(_presto.fourierprops_cen_get, _presto.fourierprops_cen_set)
    __swig_setmethods__["cenerr"] = _presto.fourierprops_cenerr_set
    __swig_getmethods__["cenerr"] = _presto.fourierprops_cenerr_get
    if _newclass:cenerr = property(_presto.fourierprops_cenerr_get, _presto.fourierprops_cenerr_set)
    __swig_setmethods__["pur"] = _presto.fourierprops_pur_set
    __swig_getmethods__["pur"] = _presto.fourierprops_pur_get
    if _newclass:pur = property(_presto.fourierprops_pur_get, _presto.fourierprops_pur_set)
    __swig_setmethods__["purerr"] = _presto.fourierprops_purerr_set
    __swig_getmethods__["purerr"] = _presto.fourierprops_purerr_get
    if _newclass:purerr = property(_presto.fourierprops_purerr_get, _presto.fourierprops_purerr_set)
    __swig_setmethods__["locpow"] = _presto.fourierprops_locpow_set
    __swig_getmethods__["locpow"] = _presto.fourierprops_locpow_get
    if _newclass:locpow = property(_presto.fourierprops_locpow_get, _presto.fourierprops_locpow_set)
    def __init__(self, *args):
        _swig_setattr(self, fourierprops, 'this', _presto.new_fourierprops(*args))
        _swig_setattr(self, fourierprops, 'thisown', 1)
    def __del__(self, destroy=_presto.delete_fourierprops):
        try:
            if self.thisown: destroy(self)
        except: pass

class fourierpropsPtr(fourierprops):
    def __init__(self, this):
        _swig_setattr(self, fourierprops, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, fourierprops, 'thisown', 0)
        _swig_setattr(self, fourierprops,self.__class__,fourierprops)
_presto.fourierprops_swigregister(fourierpropsPtr)

class binaryprops(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, binaryprops, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, binaryprops, name)
    def __repr__(self):
        return "<C binaryprops instance at %s>" % (self.this,)
    __swig_setmethods__["ppsr"] = _presto.binaryprops_ppsr_set
    __swig_getmethods__["ppsr"] = _presto.binaryprops_ppsr_get
    if _newclass:ppsr = property(_presto.binaryprops_ppsr_get, _presto.binaryprops_ppsr_set)
    __swig_setmethods__["fpsr"] = _presto.binaryprops_fpsr_set
    __swig_getmethods__["fpsr"] = _presto.binaryprops_fpsr_get
    if _newclass:fpsr = property(_presto.binaryprops_fpsr_get, _presto.binaryprops_fpsr_set)
    __swig_setmethods__["rpsr"] = _presto.binaryprops_rpsr_set
    __swig_getmethods__["rpsr"] = _presto.binaryprops_rpsr_get
    if _newclass:rpsr = property(_presto.binaryprops_rpsr_get, _presto.binaryprops_rpsr_set)
    __swig_setmethods__["pbin"] = _presto.binaryprops_pbin_set
    __swig_getmethods__["pbin"] = _presto.binaryprops_pbin_get
    if _newclass:pbin = property(_presto.binaryprops_pbin_get, _presto.binaryprops_pbin_set)
    __swig_setmethods__["rbin"] = _presto.binaryprops_rbin_set
    __swig_getmethods__["rbin"] = _presto.binaryprops_rbin_get
    if _newclass:rbin = property(_presto.binaryprops_rbin_get, _presto.binaryprops_rbin_set)
    __swig_setmethods__["z"] = _presto.binaryprops_z_set
    __swig_getmethods__["z"] = _presto.binaryprops_z_get
    if _newclass:z = property(_presto.binaryprops_z_get, _presto.binaryprops_z_set)
    __swig_setmethods__["asinic"] = _presto.binaryprops_asinic_set
    __swig_getmethods__["asinic"] = _presto.binaryprops_asinic_get
    if _newclass:asinic = property(_presto.binaryprops_asinic_get, _presto.binaryprops_asinic_set)
    __swig_setmethods__["rdetect"] = _presto.binaryprops_rdetect_set
    __swig_getmethods__["rdetect"] = _presto.binaryprops_rdetect_get
    if _newclass:rdetect = property(_presto.binaryprops_rdetect_get, _presto.binaryprops_rdetect_set)
    __swig_setmethods__["nfftbins"] = _presto.binaryprops_nfftbins_set
    __swig_getmethods__["nfftbins"] = _presto.binaryprops_nfftbins_get
    if _newclass:nfftbins = property(_presto.binaryprops_nfftbins_get, _presto.binaryprops_nfftbins_set)
    __swig_setmethods__["lowbin"] = _presto.binaryprops_lowbin_set
    __swig_getmethods__["lowbin"] = _presto.binaryprops_lowbin_get
    if _newclass:lowbin = property(_presto.binaryprops_lowbin_get, _presto.binaryprops_lowbin_set)
    __swig_setmethods__["ppsrerr"] = _presto.binaryprops_ppsrerr_set
    __swig_getmethods__["ppsrerr"] = _presto.binaryprops_ppsrerr_get
    if _newclass:ppsrerr = property(_presto.binaryprops_ppsrerr_get, _presto.binaryprops_ppsrerr_set)
    __swig_setmethods__["fpsrerr"] = _presto.binaryprops_fpsrerr_set
    __swig_getmethods__["fpsrerr"] = _presto.binaryprops_fpsrerr_get
    if _newclass:fpsrerr = property(_presto.binaryprops_fpsrerr_get, _presto.binaryprops_fpsrerr_set)
    __swig_setmethods__["rpsrerr"] = _presto.binaryprops_rpsrerr_set
    __swig_getmethods__["rpsrerr"] = _presto.binaryprops_rpsrerr_get
    if _newclass:rpsrerr = property(_presto.binaryprops_rpsrerr_get, _presto.binaryprops_rpsrerr_set)
    __swig_setmethods__["pbinerr"] = _presto.binaryprops_pbinerr_set
    __swig_getmethods__["pbinerr"] = _presto.binaryprops_pbinerr_get
    if _newclass:pbinerr = property(_presto.binaryprops_pbinerr_get, _presto.binaryprops_pbinerr_set)
    __swig_setmethods__["rbinerr"] = _presto.binaryprops_rbinerr_set
    __swig_getmethods__["rbinerr"] = _presto.binaryprops_rbinerr_get
    if _newclass:rbinerr = property(_presto.binaryprops_rbinerr_get, _presto.binaryprops_rbinerr_set)
    __swig_setmethods__["zerr"] = _presto.binaryprops_zerr_set
    __swig_getmethods__["zerr"] = _presto.binaryprops_zerr_get
    if _newclass:zerr = property(_presto.binaryprops_zerr_get, _presto.binaryprops_zerr_set)
    __swig_setmethods__["asinicerr"] = _presto.binaryprops_asinicerr_set
    __swig_getmethods__["asinicerr"] = _presto.binaryprops_asinicerr_get
    if _newclass:asinicerr = property(_presto.binaryprops_asinicerr_get, _presto.binaryprops_asinicerr_set)
    __swig_setmethods__["rdetecterr"] = _presto.binaryprops_rdetecterr_set
    __swig_getmethods__["rdetecterr"] = _presto.binaryprops_rdetecterr_get
    if _newclass:rdetecterr = property(_presto.binaryprops_rdetecterr_get, _presto.binaryprops_rdetecterr_set)
    __swig_setmethods__["sig"] = _presto.binaryprops_sig_set
    __swig_getmethods__["sig"] = _presto.binaryprops_sig_get
    if _newclass:sig = property(_presto.binaryprops_sig_get, _presto.binaryprops_sig_set)
    __swig_setmethods__["phs"] = _presto.binaryprops_phs_set
    __swig_getmethods__["phs"] = _presto.binaryprops_phs_get
    if _newclass:phs = property(_presto.binaryprops_phs_get, _presto.binaryprops_phs_set)
    __swig_setmethods__["phserr"] = _presto.binaryprops_phserr_set
    __swig_getmethods__["phserr"] = _presto.binaryprops_phserr_get
    if _newclass:phserr = property(_presto.binaryprops_phserr_get, _presto.binaryprops_phserr_set)
    __swig_setmethods__["cen"] = _presto.binaryprops_cen_set
    __swig_getmethods__["cen"] = _presto.binaryprops_cen_get
    if _newclass:cen = property(_presto.binaryprops_cen_get, _presto.binaryprops_cen_set)
    __swig_setmethods__["cenerr"] = _presto.binaryprops_cenerr_set
    __swig_getmethods__["cenerr"] = _presto.binaryprops_cenerr_get
    if _newclass:cenerr = property(_presto.binaryprops_cenerr_get, _presto.binaryprops_cenerr_set)
    __swig_setmethods__["pur"] = _presto.binaryprops_pur_set
    __swig_getmethods__["pur"] = _presto.binaryprops_pur_get
    if _newclass:pur = property(_presto.binaryprops_pur_get, _presto.binaryprops_pur_set)
    __swig_setmethods__["purerr"] = _presto.binaryprops_purerr_set
    __swig_getmethods__["purerr"] = _presto.binaryprops_purerr_get
    if _newclass:purerr = property(_presto.binaryprops_purerr_get, _presto.binaryprops_purerr_set)
    __swig_setmethods__["pow"] = _presto.binaryprops_pow_set
    __swig_getmethods__["pow"] = _presto.binaryprops_pow_get
    if _newclass:pow = property(_presto.binaryprops_pow_get, _presto.binaryprops_pow_set)
    __swig_setmethods__["powerr"] = _presto.binaryprops_powerr_set
    __swig_getmethods__["powerr"] = _presto.binaryprops_powerr_get
    if _newclass:powerr = property(_presto.binaryprops_powerr_get, _presto.binaryprops_powerr_set)
    def __init__(self, *args):
        _swig_setattr(self, binaryprops, 'this', _presto.new_binaryprops(*args))
        _swig_setattr(self, binaryprops, 'thisown', 1)
    def __del__(self, destroy=_presto.delete_binaryprops):
        try:
            if self.thisown: destroy(self)
        except: pass

class binarypropsPtr(binaryprops):
    def __init__(self, this):
        _swig_setattr(self, binaryprops, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, binaryprops, 'thisown', 0)
        _swig_setattr(self, binaryprops,self.__class__,binaryprops)
_presto.binaryprops_swigregister(binarypropsPtr)

class rawbincand(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, rawbincand, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, rawbincand, name)
    def __repr__(self):
        return "<C rawbincand instance at %s>" % (self.this,)
    __swig_setmethods__["full_N"] = _presto.rawbincand_full_N_set
    __swig_getmethods__["full_N"] = _presto.rawbincand_full_N_get
    if _newclass:full_N = property(_presto.rawbincand_full_N_get, _presto.rawbincand_full_N_set)
    __swig_setmethods__["full_T"] = _presto.rawbincand_full_T_set
    __swig_getmethods__["full_T"] = _presto.rawbincand_full_T_get
    if _newclass:full_T = property(_presto.rawbincand_full_T_get, _presto.rawbincand_full_T_set)
    __swig_setmethods__["full_lo_r"] = _presto.rawbincand_full_lo_r_set
    __swig_getmethods__["full_lo_r"] = _presto.rawbincand_full_lo_r_get
    if _newclass:full_lo_r = property(_presto.rawbincand_full_lo_r_get, _presto.rawbincand_full_lo_r_set)
    __swig_setmethods__["mini_N"] = _presto.rawbincand_mini_N_set
    __swig_getmethods__["mini_N"] = _presto.rawbincand_mini_N_get
    if _newclass:mini_N = property(_presto.rawbincand_mini_N_get, _presto.rawbincand_mini_N_set)
    __swig_setmethods__["mini_r"] = _presto.rawbincand_mini_r_set
    __swig_getmethods__["mini_r"] = _presto.rawbincand_mini_r_get
    if _newclass:mini_r = property(_presto.rawbincand_mini_r_get, _presto.rawbincand_mini_r_set)
    __swig_setmethods__["mini_power"] = _presto.rawbincand_mini_power_set
    __swig_getmethods__["mini_power"] = _presto.rawbincand_mini_power_get
    if _newclass:mini_power = property(_presto.rawbincand_mini_power_get, _presto.rawbincand_mini_power_set)
    __swig_setmethods__["mini_numsum"] = _presto.rawbincand_mini_numsum_set
    __swig_getmethods__["mini_numsum"] = _presto.rawbincand_mini_numsum_get
    if _newclass:mini_numsum = property(_presto.rawbincand_mini_numsum_get, _presto.rawbincand_mini_numsum_set)
    __swig_setmethods__["mini_sigma"] = _presto.rawbincand_mini_sigma_set
    __swig_getmethods__["mini_sigma"] = _presto.rawbincand_mini_sigma_get
    if _newclass:mini_sigma = property(_presto.rawbincand_mini_sigma_get, _presto.rawbincand_mini_sigma_set)
    __swig_setmethods__["psr_p"] = _presto.rawbincand_psr_p_set
    __swig_getmethods__["psr_p"] = _presto.rawbincand_psr_p_get
    if _newclass:psr_p = property(_presto.rawbincand_psr_p_get, _presto.rawbincand_psr_p_set)
    __swig_setmethods__["orb_p"] = _presto.rawbincand_orb_p_set
    __swig_getmethods__["orb_p"] = _presto.rawbincand_orb_p_get
    if _newclass:orb_p = property(_presto.rawbincand_orb_p_get, _presto.rawbincand_orb_p_set)
    def __init__(self, *args):
        _swig_setattr(self, rawbincand, 'this', _presto.new_rawbincand(*args))
        _swig_setattr(self, rawbincand, 'thisown', 1)
    def __del__(self, destroy=_presto.delete_rawbincand):
        try:
            if self.thisown: destroy(self)
        except: pass

class rawbincandPtr(rawbincand):
    def __init__(self, this):
        _swig_setattr(self, rawbincand, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, rawbincand, 'thisown', 0)
        _swig_setattr(self, rawbincand,self.__class__,rawbincand)
_presto.rawbincand_swigregister(rawbincandPtr)

class foldstats(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, foldstats, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, foldstats, name)
    def __repr__(self):
        return "<C foldstats instance at %s>" % (self.this,)
    __swig_setmethods__["numdata"] = _presto.foldstats_numdata_set
    __swig_getmethods__["numdata"] = _presto.foldstats_numdata_get
    if _newclass:numdata = property(_presto.foldstats_numdata_get, _presto.foldstats_numdata_set)
    __swig_setmethods__["data_avg"] = _presto.foldstats_data_avg_set
    __swig_getmethods__["data_avg"] = _presto.foldstats_data_avg_get
    if _newclass:data_avg = property(_presto.foldstats_data_avg_get, _presto.foldstats_data_avg_set)
    __swig_setmethods__["data_var"] = _presto.foldstats_data_var_set
    __swig_getmethods__["data_var"] = _presto.foldstats_data_var_get
    if _newclass:data_var = property(_presto.foldstats_data_var_get, _presto.foldstats_data_var_set)
    __swig_setmethods__["numprof"] = _presto.foldstats_numprof_set
    __swig_getmethods__["numprof"] = _presto.foldstats_numprof_get
    if _newclass:numprof = property(_presto.foldstats_numprof_get, _presto.foldstats_numprof_set)
    __swig_setmethods__["prof_avg"] = _presto.foldstats_prof_avg_set
    __swig_getmethods__["prof_avg"] = _presto.foldstats_prof_avg_get
    if _newclass:prof_avg = property(_presto.foldstats_prof_avg_get, _presto.foldstats_prof_avg_set)
    __swig_setmethods__["prof_var"] = _presto.foldstats_prof_var_set
    __swig_getmethods__["prof_var"] = _presto.foldstats_prof_var_get
    if _newclass:prof_var = property(_presto.foldstats_prof_var_get, _presto.foldstats_prof_var_set)
    __swig_setmethods__["redchi"] = _presto.foldstats_redchi_set
    __swig_getmethods__["redchi"] = _presto.foldstats_redchi_get
    if _newclass:redchi = property(_presto.foldstats_redchi_get, _presto.foldstats_redchi_set)
    def __init__(self, *args):
        _swig_setattr(self, foldstats, 'this', _presto.new_foldstats(*args))
        _swig_setattr(self, foldstats, 'thisown', 1)
    def __del__(self, destroy=_presto.delete_foldstats):
        try:
            if self.thisown: destroy(self)
        except: pass

class foldstatsPtr(foldstats):
    def __init__(self, this):
        _swig_setattr(self, foldstats, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, foldstats, 'thisown', 0)
        _swig_setattr(self, foldstats,self.__class__,foldstats)
_presto.foldstats_swigregister(foldstatsPtr)


get_psr_at_epoch = _presto.get_psr_at_epoch

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

