# This file was created automatically by SWIG.
import prestoc
class orbitparamsPtr :
    def __init__(self,this):
        self.this = this
        self.thisown = 0
    def __del__(self):
        if self.thisown == 1 :
            prestoc.delete_orbitparams(self.this)
    def __setattr__(self,name,value):
        if name == "p" :
            prestoc.orbitparams_p_set(self.this,value)
            return
        if name == "e" :
            prestoc.orbitparams_e_set(self.this,value)
            return
        if name == "x" :
            prestoc.orbitparams_x_set(self.this,value)
            return
        if name == "w" :
            prestoc.orbitparams_w_set(self.this,value)
            return
        if name == "t" :
            prestoc.orbitparams_t_set(self.this,value)
            return
        if name == "pd" :
            prestoc.orbitparams_pd_set(self.this,value)
            return
        if name == "wd" :
            prestoc.orbitparams_wd_set(self.this,value)
            return
        self.__dict__[name] = value
    def __getattr__(self,name):
        if name == "p" : 
            return prestoc.orbitparams_p_get(self.this)
        if name == "e" : 
            return prestoc.orbitparams_e_get(self.this)
        if name == "x" : 
            return prestoc.orbitparams_x_get(self.this)
        if name == "w" : 
            return prestoc.orbitparams_w_get(self.this)
        if name == "t" : 
            return prestoc.orbitparams_t_get(self.this)
        if name == "pd" : 
            return prestoc.orbitparams_pd_get(self.this)
        if name == "wd" : 
            return prestoc.orbitparams_wd_get(self.this)
        raise AttributeError,name
    def __repr__(self):
        return "<C orbitparams instance>"
class orbitparams(orbitparamsPtr):
    def __init__(self) :
        self.this = prestoc.new_orbitparams()
        self.thisown = 1




class psrdataPtr :
    def __init__(self,this):
        self.this = this
        self.thisown = 0
    def __del__(self):
        if self.thisown == 1 :
            prestoc.delete_psrdata(self.this)
    def __setattr__(self,name,value):
        if name == "ra2000" :
            prestoc.psrdata_ra2000_set(self.this,value)
            return
        if name == "ra1950" :
            prestoc.psrdata_ra1950_set(self.this,value)
            return
        if name == "rae" :
            prestoc.psrdata_rae_set(self.this,value)
            return
        if name == "dec2000" :
            prestoc.psrdata_dec2000_set(self.this,value)
            return
        if name == "dec1950" :
            prestoc.psrdata_dec1950_set(self.this,value)
            return
        if name == "dece" :
            prestoc.psrdata_dece_set(self.this,value)
            return
        if name == "dmin__" :
            prestoc.psrdata_dmin___set(self.this,value)
            return
        if name == "dmax__" :
            prestoc.psrdata_dmax___set(self.this,value)
            return
        if name == "dist" :
            prestoc.psrdata_dist_set(self.this,value)
            return
        if name == "ldeg" :
            prestoc.psrdata_ldeg_set(self.this,value)
            return
        if name == "bdeg" :
            prestoc.psrdata_bdeg_set(self.this,value)
            return
        if name == "pmra" :
            prestoc.psrdata_pmra_set(self.this,value)
            return
        if name == "pmrae" :
            prestoc.psrdata_pmrae_set(self.this,value)
            return
        if name == "pmdec" :
            prestoc.psrdata_pmdec_set(self.this,value)
            return
        if name == "pmdece" :
            prestoc.psrdata_pmdece_set(self.this,value)
            return
        if name == "posepoch" :
            prestoc.psrdata_posepoch_set(self.this,value)
            return
        if name == "p" :
            prestoc.psrdata_p_set(self.this,value)
            return
        if name == "pe" :
            prestoc.psrdata_pe_set(self.this,value)
            return
        if name == "pdot" :
            prestoc.psrdata_pdot_set(self.this,value)
            return
        if name == "pdote" :
            prestoc.psrdata_pdote_set(self.this,value)
            return
        if name == "f2" :
            prestoc.psrdata_f2_set(self.this,value)
            return
        if name == "f2e" :
            prestoc.psrdata_f2e_set(self.this,value)
            return
        if name == "f3" :
            prestoc.psrdata_f3_set(self.this,value)
            return
        if name == "f3e" :
            prestoc.psrdata_f3e_set(self.this,value)
            return
        if name == "epoch" :
            prestoc.psrdata_epoch_set(self.this,value)
            return
        if name == "dm" :
            prestoc.psrdata_dm_set(self.this,value)
            return
        if name == "dme" :
            prestoc.psrdata_dme_set(self.this,value)
            return
        if name == "rm" :
            prestoc.psrdata_rm_set(self.this,value)
            return
        if name == "rme" :
            prestoc.psrdata_rme_set(self.this,value)
            return
        if name == "we" :
            prestoc.psrdata_we_set(self.this,value)
            return
        if name == "w50" :
            prestoc.psrdata_w50_set(self.this,value)
            return
        if name == "w10" :
            prestoc.psrdata_w10_set(self.this,value)
            return
        if name == "s400" :
            prestoc.psrdata_s400_set(self.this,value)
            return
        if name == "s600" :
            prestoc.psrdata_s600_set(self.this,value)
            return
        if name == "s1400" :
            prestoc.psrdata_s1400_set(self.this,value)
            return
        if name == "tau" :
            prestoc.psrdata_tau_set(self.this,value)
            return
        if name == "t408" :
            prestoc.psrdata_t408_set(self.this,value)
            return
        if name == "distmod" :
            prestoc.psrdata_distmod_set(self.this,value)
            return
        if name == "lum" :
            prestoc.psrdata_lum_set(self.this,value)
            return
        if name == "bsurf" :
            prestoc.psrdata_bsurf_set(self.this,value)
            return
        if name == "age" :
            prestoc.psrdata_age_set(self.this,value)
            return
        if name == "edot" :
            prestoc.psrdata_edot_set(self.this,value)
            return
        if name == "pb" :
            prestoc.psrdata_pb_set(self.this,value)
            return
        if name == "pbe" :
            prestoc.psrdata_pbe_set(self.this,value)
            return
        if name == "a1" :
            prestoc.psrdata_a1_set(self.this,value)
            return
        if name == "a1e" :
            prestoc.psrdata_a1e_set(self.this,value)
            return
        if name == "om" :
            prestoc.psrdata_om_set(self.this,value)
            return
        if name == "ome" :
            prestoc.psrdata_ome_set(self.this,value)
            return
        if name == "omdot" :
            prestoc.psrdata_omdot_set(self.this,value)
            return
        if name == "omdote" :
            prestoc.psrdata_omdote_set(self.this,value)
            return
        if name == "e" :
            prestoc.psrdata_e_set(self.this,value)
            return
        if name == "ee" :
            prestoc.psrdata_ee_set(self.this,value)
            return
        if name == "t0" :
            prestoc.psrdata_t0_set(self.this,value)
            return
        if name == "t0e" :
            prestoc.psrdata_t0e_set(self.this,value)
            return
        if name == "gamma" :
            prestoc.psrdata_gamma_set(self.this,value)
            return
        if name == "gammae" :
            prestoc.psrdata_gammae_set(self.this,value)
            return
        if name == "pbdot" :
            prestoc.psrdata_pbdot_set(self.this,value)
            return
        if name == "pbdote" :
            prestoc.psrdata_pbdote_set(self.this,value)
            return
        if name == "si" :
            prestoc.psrdata_si_set(self.this,value)
            return
        if name == "sie" :
            prestoc.psrdata_sie_set(self.this,value)
            return
        if name == "r__" :
            prestoc.psrdata_r___set(self.this,value)
            return
        if name == "re" :
            prestoc.psrdata_re_set(self.this,value)
            return
        if name == "pb2" :
            prestoc.psrdata_pb2_set(self.this,value)
            return
        if name == "pb2e" :
            prestoc.psrdata_pb2e_set(self.this,value)
            return
        if name == "a12" :
            prestoc.psrdata_a12_set(self.this,value)
            return
        if name == "a12e" :
            prestoc.psrdata_a12e_set(self.this,value)
            return
        if name == "om2" :
            prestoc.psrdata_om2_set(self.this,value)
            return
        if name == "om2e" :
            prestoc.psrdata_om2e_set(self.this,value)
            return
        if name == "omdot2" :
            prestoc.psrdata_omdot2_set(self.this,value)
            return
        if name == "omdot2e" :
            prestoc.psrdata_omdot2e_set(self.this,value)
            return
        if name == "e2" :
            prestoc.psrdata_e2_set(self.this,value)
            return
        if name == "e2e" :
            prestoc.psrdata_e2e_set(self.this,value)
            return
        if name == "t02" :
            prestoc.psrdata_t02_set(self.this,value)
            return
        if name == "t02e" :
            prestoc.psrdata_t02e_set(self.this,value)
            return
        if name == "gamma2" :
            prestoc.psrdata_gamma2_set(self.this,value)
            return
        if name == "gamma2e" :
            prestoc.psrdata_gamma2e_set(self.this,value)
            return
        if name == "pbdot2" :
            prestoc.psrdata_pbdot2_set(self.this,value)
            return
        if name == "pbdot2e" :
            prestoc.psrdata_pbdot2e_set(self.this,value)
            return
        if name == "si2" :
            prestoc.psrdata_si2_set(self.this,value)
            return
        if name == "si2e" :
            prestoc.psrdata_si2e_set(self.this,value)
            return
        if name == "r2" :
            prestoc.psrdata_r2_set(self.this,value)
            return
        if name == "r2e" :
            prestoc.psrdata_r2e_set(self.this,value)
            return
        if name == "nscode" :
            prestoc.psrdata_nscode_set(self.this,value)
            return
        if name == "ndflag" :
            prestoc.psrdata_ndflag_set(self.this,value)
            return
        if name == "ntauflag" :
            prestoc.psrdata_ntauflag_set(self.this,value)
            return
        if name == "ntype" :
            prestoc.psrdata_ntype_set(self.this,value)
            return
        if name == "modcode" :
            prestoc.psrdata_modcode_set(self.this,value)
            return
        if name == "limcode" :
            prestoc.psrdata_limcode_set(self.this,value)
            return
        if name == "ibin" :
            prestoc.psrdata_ibin_set(self.this,value)
            return
        if name == "jname" :
            prestoc.psrdata_jname_set(self.this,value)
            return
        if name == "bname" :
            prestoc.psrdata_bname_set(self.this,value)
            return
        if name == "lcode" :
            prestoc.psrdata_lcode_set(self.this,value)
            return
        if name == "ucode" :
            prestoc.psrdata_ucode_set(self.this,value)
            return
        self.__dict__[name] = value
    def __getattr__(self,name):
        if name == "ra2000" : 
            return prestoc.psrdata_ra2000_get(self.this)
        if name == "ra1950" : 
            return prestoc.psrdata_ra1950_get(self.this)
        if name == "rae" : 
            return prestoc.psrdata_rae_get(self.this)
        if name == "dec2000" : 
            return prestoc.psrdata_dec2000_get(self.this)
        if name == "dec1950" : 
            return prestoc.psrdata_dec1950_get(self.this)
        if name == "dece" : 
            return prestoc.psrdata_dece_get(self.this)
        if name == "dmin__" : 
            return prestoc.psrdata_dmin___get(self.this)
        if name == "dmax__" : 
            return prestoc.psrdata_dmax___get(self.this)
        if name == "dist" : 
            return prestoc.psrdata_dist_get(self.this)
        if name == "ldeg" : 
            return prestoc.psrdata_ldeg_get(self.this)
        if name == "bdeg" : 
            return prestoc.psrdata_bdeg_get(self.this)
        if name == "pmra" : 
            return prestoc.psrdata_pmra_get(self.this)
        if name == "pmrae" : 
            return prestoc.psrdata_pmrae_get(self.this)
        if name == "pmdec" : 
            return prestoc.psrdata_pmdec_get(self.this)
        if name == "pmdece" : 
            return prestoc.psrdata_pmdece_get(self.this)
        if name == "posepoch" : 
            return prestoc.psrdata_posepoch_get(self.this)
        if name == "p" : 
            return prestoc.psrdata_p_get(self.this)
        if name == "pe" : 
            return prestoc.psrdata_pe_get(self.this)
        if name == "pdot" : 
            return prestoc.psrdata_pdot_get(self.this)
        if name == "pdote" : 
            return prestoc.psrdata_pdote_get(self.this)
        if name == "f2" : 
            return prestoc.psrdata_f2_get(self.this)
        if name == "f2e" : 
            return prestoc.psrdata_f2e_get(self.this)
        if name == "f3" : 
            return prestoc.psrdata_f3_get(self.this)
        if name == "f3e" : 
            return prestoc.psrdata_f3e_get(self.this)
        if name == "epoch" : 
            return prestoc.psrdata_epoch_get(self.this)
        if name == "dm" : 
            return prestoc.psrdata_dm_get(self.this)
        if name == "dme" : 
            return prestoc.psrdata_dme_get(self.this)
        if name == "rm" : 
            return prestoc.psrdata_rm_get(self.this)
        if name == "rme" : 
            return prestoc.psrdata_rme_get(self.this)
        if name == "we" : 
            return prestoc.psrdata_we_get(self.this)
        if name == "w50" : 
            return prestoc.psrdata_w50_get(self.this)
        if name == "w10" : 
            return prestoc.psrdata_w10_get(self.this)
        if name == "s400" : 
            return prestoc.psrdata_s400_get(self.this)
        if name == "s600" : 
            return prestoc.psrdata_s600_get(self.this)
        if name == "s1400" : 
            return prestoc.psrdata_s1400_get(self.this)
        if name == "tau" : 
            return prestoc.psrdata_tau_get(self.this)
        if name == "t408" : 
            return prestoc.psrdata_t408_get(self.this)
        if name == "distmod" : 
            return prestoc.psrdata_distmod_get(self.this)
        if name == "lum" : 
            return prestoc.psrdata_lum_get(self.this)
        if name == "bsurf" : 
            return prestoc.psrdata_bsurf_get(self.this)
        if name == "age" : 
            return prestoc.psrdata_age_get(self.this)
        if name == "edot" : 
            return prestoc.psrdata_edot_get(self.this)
        if name == "pb" : 
            return prestoc.psrdata_pb_get(self.this)
        if name == "pbe" : 
            return prestoc.psrdata_pbe_get(self.this)
        if name == "a1" : 
            return prestoc.psrdata_a1_get(self.this)
        if name == "a1e" : 
            return prestoc.psrdata_a1e_get(self.this)
        if name == "om" : 
            return prestoc.psrdata_om_get(self.this)
        if name == "ome" : 
            return prestoc.psrdata_ome_get(self.this)
        if name == "omdot" : 
            return prestoc.psrdata_omdot_get(self.this)
        if name == "omdote" : 
            return prestoc.psrdata_omdote_get(self.this)
        if name == "e" : 
            return prestoc.psrdata_e_get(self.this)
        if name == "ee" : 
            return prestoc.psrdata_ee_get(self.this)
        if name == "t0" : 
            return prestoc.psrdata_t0_get(self.this)
        if name == "t0e" : 
            return prestoc.psrdata_t0e_get(self.this)
        if name == "gamma" : 
            return prestoc.psrdata_gamma_get(self.this)
        if name == "gammae" : 
            return prestoc.psrdata_gammae_get(self.this)
        if name == "pbdot" : 
            return prestoc.psrdata_pbdot_get(self.this)
        if name == "pbdote" : 
            return prestoc.psrdata_pbdote_get(self.this)
        if name == "si" : 
            return prestoc.psrdata_si_get(self.this)
        if name == "sie" : 
            return prestoc.psrdata_sie_get(self.this)
        if name == "r__" : 
            return prestoc.psrdata_r___get(self.this)
        if name == "re" : 
            return prestoc.psrdata_re_get(self.this)
        if name == "pb2" : 
            return prestoc.psrdata_pb2_get(self.this)
        if name == "pb2e" : 
            return prestoc.psrdata_pb2e_get(self.this)
        if name == "a12" : 
            return prestoc.psrdata_a12_get(self.this)
        if name == "a12e" : 
            return prestoc.psrdata_a12e_get(self.this)
        if name == "om2" : 
            return prestoc.psrdata_om2_get(self.this)
        if name == "om2e" : 
            return prestoc.psrdata_om2e_get(self.this)
        if name == "omdot2" : 
            return prestoc.psrdata_omdot2_get(self.this)
        if name == "omdot2e" : 
            return prestoc.psrdata_omdot2e_get(self.this)
        if name == "e2" : 
            return prestoc.psrdata_e2_get(self.this)
        if name == "e2e" : 
            return prestoc.psrdata_e2e_get(self.this)
        if name == "t02" : 
            return prestoc.psrdata_t02_get(self.this)
        if name == "t02e" : 
            return prestoc.psrdata_t02e_get(self.this)
        if name == "gamma2" : 
            return prestoc.psrdata_gamma2_get(self.this)
        if name == "gamma2e" : 
            return prestoc.psrdata_gamma2e_get(self.this)
        if name == "pbdot2" : 
            return prestoc.psrdata_pbdot2_get(self.this)
        if name == "pbdot2e" : 
            return prestoc.psrdata_pbdot2e_get(self.this)
        if name == "si2" : 
            return prestoc.psrdata_si2_get(self.this)
        if name == "si2e" : 
            return prestoc.psrdata_si2e_get(self.this)
        if name == "r2" : 
            return prestoc.psrdata_r2_get(self.this)
        if name == "r2e" : 
            return prestoc.psrdata_r2e_get(self.this)
        if name == "nscode" : 
            return prestoc.psrdata_nscode_get(self.this)
        if name == "ndflag" : 
            return prestoc.psrdata_ndflag_get(self.this)
        if name == "ntauflag" : 
            return prestoc.psrdata_ntauflag_get(self.this)
        if name == "ntype" : 
            return prestoc.psrdata_ntype_get(self.this)
        if name == "modcode" : 
            return prestoc.psrdata_modcode_get(self.this)
        if name == "limcode" : 
            return prestoc.psrdata_limcode_get(self.this)
        if name == "ibin" : 
            return prestoc.psrdata_ibin_get(self.this)
        if name == "jname" : 
            return prestoc.psrdata_jname_get(self.this)
        if name == "bname" : 
            return prestoc.psrdata_bname_get(self.this)
        if name == "lcode" : 
            return prestoc.psrdata_lcode_get(self.this)
        if name == "ucode" : 
            return prestoc.psrdata_ucode_get(self.this)
        raise AttributeError,name
    def __repr__(self):
        return "<C psrdata instance>"
class psrdata(psrdataPtr):
    def __init__(self) :
        self.this = prestoc.new_psrdata()
        self.thisown = 1




class psrparamsPtr :
    def __init__(self,this):
        self.this = this
        self.thisown = 0
    def __del__(self):
        if self.thisown == 1 :
            prestoc.delete_psrparams(self.this)
    def __setattr__(self,name,value):
        if name == "jname" :
            prestoc.psrparams_jname_set(self.this,value)
            return
        if name == "bname" :
            prestoc.psrparams_bname_set(self.this,value)
            return
        if name == "ntype" :
            prestoc.psrparams_ntype_set(self.this,value)
            return
        if name == "ra2000" :
            prestoc.psrparams_ra2000_set(self.this,value)
            return
        if name == "dec2000" :
            prestoc.psrparams_dec2000_set(self.this,value)
            return
        if name == "dm" :
            prestoc.psrparams_dm_set(self.this,value)
            return
        if name == "dist" :
            prestoc.psrparams_dist_set(self.this,value)
            return
        if name == "fwhm" :
            prestoc.psrparams_fwhm_set(self.this,value)
            return
        if name == "timepoch" :
            prestoc.psrparams_timepoch_set(self.this,value)
            return
        if name == "p" :
            prestoc.psrparams_p_set(self.this,value)
            return
        if name == "pd" :
            prestoc.psrparams_pd_set(self.this,value)
            return
        if name == "pdd" :
            prestoc.psrparams_pdd_set(self.this,value)
            return
        if name == "f" :
            prestoc.psrparams_f_set(self.this,value)
            return
        if name == "fd" :
            prestoc.psrparams_fd_set(self.this,value)
            return
        if name == "fdd" :
            prestoc.psrparams_fdd_set(self.this,value)
            return
        if name == "orb" :
            prestoc.psrparams_orb_set(self.this,value.this)
            return
        self.__dict__[name] = value
    def __getattr__(self,name):
        if name == "jname" : 
            return prestoc.psrparams_jname_get(self.this)
        if name == "bname" : 
            return prestoc.psrparams_bname_get(self.this)
        if name == "ntype" : 
            return prestoc.psrparams_ntype_get(self.this)
        if name == "ra2000" : 
            return prestoc.psrparams_ra2000_get(self.this)
        if name == "dec2000" : 
            return prestoc.psrparams_dec2000_get(self.this)
        if name == "dm" : 
            return prestoc.psrparams_dm_get(self.this)
        if name == "dist" : 
            return prestoc.psrparams_dist_get(self.this)
        if name == "fwhm" : 
            return prestoc.psrparams_fwhm_get(self.this)
        if name == "timepoch" : 
            return prestoc.psrparams_timepoch_get(self.this)
        if name == "p" : 
            return prestoc.psrparams_p_get(self.this)
        if name == "pd" : 
            return prestoc.psrparams_pd_get(self.this)
        if name == "pdd" : 
            return prestoc.psrparams_pdd_get(self.this)
        if name == "f" : 
            return prestoc.psrparams_f_get(self.this)
        if name == "fd" : 
            return prestoc.psrparams_fd_get(self.this)
        if name == "fdd" : 
            return prestoc.psrparams_fdd_get(self.this)
        if name == "orb" : 
            return orbitparamsPtr(prestoc.psrparams_orb_get(self.this))
        raise AttributeError,name
    def __repr__(self):
        return "<C psrparams instance>"
class psrparams(psrparamsPtr):
    def __init__(self) :
        self.this = prestoc.new_psrparams()
        self.thisown = 1




class DoubleArrayPtr :
    def __init__(self,this):
        self.this = this
        self.thisown = 0
    def __del__(self):
        if self.thisown == 1 :
            prestoc.delete_DoubleArray(self.this)
    def __getitem__(self,arg0):
        val = prestoc.DoubleArray___getitem__(self.this,arg0)
        return val
    def __setitem__(self,arg0,arg1):
        val = prestoc.DoubleArray___setitem__(self.this,arg0,arg1)
        return val
    def __setattr__(self,name,value):
        if name == "dptr" :
            prestoc.DoubleArray_dptr_set(self.this,value)
            return
        self.__dict__[name] = value
    def __getattr__(self,name):
        if name == "dptr" : 
            return prestoc.DoubleArray_dptr_get(self.this)
        raise AttributeError,name
    def __repr__(self):
        return "<C DoubleArray instance>"
class DoubleArray(DoubleArrayPtr):
    def __init__(self,arg0) :
        self.this = prestoc.new_DoubleArray(arg0)
        self.thisown = 1




class infodataPtr :
    def __init__(self,this):
        self.this = this
        self.thisown = 0
    def __del__(self):
        if self.thisown == 1 :
            prestoc.delete_infodata(self.this)
    def __setattr__(self,name,value):
        if name == "name" :
            prestoc.infodata_name_set(self.this,value)
            return
        if name == "object" :
            prestoc.infodata_object_set(self.this,value)
            return
        if name == "ra_h" :
            prestoc.infodata_ra_h_set(self.this,value)
            return
        if name == "ra_m" :
            prestoc.infodata_ra_m_set(self.this,value)
            return
        if name == "ra_s" :
            prestoc.infodata_ra_s_set(self.this,value)
            return
        if name == "dec_d" :
            prestoc.infodata_dec_d_set(self.this,value)
            return
        if name == "dec_m" :
            prestoc.infodata_dec_m_set(self.this,value)
            return
        if name == "dec_s" :
            prestoc.infodata_dec_s_set(self.this,value)
            return
        if name == "telescope" :
            prestoc.infodata_telescope_set(self.this,value)
            return
        if name == "instrument" :
            prestoc.infodata_instrument_set(self.this,value)
            return
        if name == "observer" :
            prestoc.infodata_observer_set(self.this,value)
            return
        if name == "N" :
            prestoc.infodata_N_set(self.this,value)
            return
        if name == "dt" :
            prestoc.infodata_dt_set(self.this,value)
            return
        if name == "numonoff" :
            prestoc.infodata_numonoff_set(self.this,value)
            return
        if name == "onoff" :
            prestoc.infodata_onoff_set(self.this,value.this)
            return
        if name == "fov" :
            prestoc.infodata_fov_set(self.this,value)
            return
        if name == "mjd_i" :
            prestoc.infodata_mjd_i_set(self.this,value)
            return
        if name == "mjd_f" :
            prestoc.infodata_mjd_f_set(self.this,value)
            return
        if name == "bary" :
            prestoc.infodata_bary_set(self.this,value)
            return
        if name == "band" :
            prestoc.infodata_band_set(self.this,value)
            return
        if name == "dm" :
            prestoc.infodata_dm_set(self.this,value)
            return
        if name == "freq" :
            prestoc.infodata_freq_set(self.this,value)
            return
        if name == "freqband" :
            prestoc.infodata_freqband_set(self.this,value)
            return
        if name == "num_chan" :
            prestoc.infodata_num_chan_set(self.this,value)
            return
        if name == "chan_wid" :
            prestoc.infodata_chan_wid_set(self.this,value)
            return
        if name == "filt" :
            prestoc.infodata_filt_set(self.this,value)
            return
        if name == "wavelen" :
            prestoc.infodata_wavelen_set(self.this,value)
            return
        if name == "waveband" :
            prestoc.infodata_waveband_set(self.this,value)
            return
        if name == "energy" :
            prestoc.infodata_energy_set(self.this,value)
            return
        if name == "energyband" :
            prestoc.infodata_energyband_set(self.this,value)
            return
        if name == "analyzer" :
            prestoc.infodata_analyzer_set(self.this,value)
            return
        if name == "notes" :
            prestoc.infodata_notes_set(self.this,value)
            return
        self.__dict__[name] = value
    def __getattr__(self,name):
        if name == "name" : 
            return prestoc.infodata_name_get(self.this)
        if name == "object" : 
            return prestoc.infodata_object_get(self.this)
        if name == "ra_h" : 
            return prestoc.infodata_ra_h_get(self.this)
        if name == "ra_m" : 
            return prestoc.infodata_ra_m_get(self.this)
        if name == "ra_s" : 
            return prestoc.infodata_ra_s_get(self.this)
        if name == "dec_d" : 
            return prestoc.infodata_dec_d_get(self.this)
        if name == "dec_m" : 
            return prestoc.infodata_dec_m_get(self.this)
        if name == "dec_s" : 
            return prestoc.infodata_dec_s_get(self.this)
        if name == "telescope" : 
            return prestoc.infodata_telescope_get(self.this)
        if name == "instrument" : 
            return prestoc.infodata_instrument_get(self.this)
        if name == "observer" : 
            return prestoc.infodata_observer_get(self.this)
        if name == "N" : 
            return prestoc.infodata_N_get(self.this)
        if name == "dt" : 
            return prestoc.infodata_dt_get(self.this)
        if name == "numonoff" : 
            return prestoc.infodata_numonoff_get(self.this)
        if name == "onoff" : 
            return DoubleArrayPtr(prestoc.infodata_onoff_get(self.this))
        if name == "fov" : 
            return prestoc.infodata_fov_get(self.this)
        if name == "mjd_i" : 
            return prestoc.infodata_mjd_i_get(self.this)
        if name == "mjd_f" : 
            return prestoc.infodata_mjd_f_get(self.this)
        if name == "bary" : 
            return prestoc.infodata_bary_get(self.this)
        if name == "band" : 
            return prestoc.infodata_band_get(self.this)
        if name == "dm" : 
            return prestoc.infodata_dm_get(self.this)
        if name == "freq" : 
            return prestoc.infodata_freq_get(self.this)
        if name == "freqband" : 
            return prestoc.infodata_freqband_get(self.this)
        if name == "num_chan" : 
            return prestoc.infodata_num_chan_get(self.this)
        if name == "chan_wid" : 
            return prestoc.infodata_chan_wid_get(self.this)
        if name == "filt" : 
            return prestoc.infodata_filt_get(self.this)
        if name == "wavelen" : 
            return prestoc.infodata_wavelen_get(self.this)
        if name == "waveband" : 
            return prestoc.infodata_waveband_get(self.this)
        if name == "energy" : 
            return prestoc.infodata_energy_get(self.this)
        if name == "energyband" : 
            return prestoc.infodata_energyband_get(self.this)
        if name == "analyzer" : 
            return prestoc.infodata_analyzer_get(self.this)
        if name == "notes" : 
            return prestoc.infodata_notes_get(self.this)
        raise AttributeError,name
    def __repr__(self):
        return "<C infodata instance>"
class infodata(infodataPtr):
    def __init__(self) :
        self.this = prestoc.new_infodata()
        self.thisown = 1




class makedataPtr :
    def __init__(self,this):
        self.this = this
        self.thisown = 0
    def __del__(self):
        if self.thisown == 1 :
            prestoc.delete_makedata(self.this)
    def __setattr__(self,name,value):
        if name == "basefilenm" :
            prestoc.makedata_basefilenm_set(self.this,value)
            return
        if name == "description" :
            prestoc.makedata_description_set(self.this,value)
            return
        if name == "N" :
            prestoc.makedata_N_set(self.this,value)
            return
        if name == "next2_to_n" :
            prestoc.makedata_next2_to_n_set(self.this,value)
            return
        if name == "dt" :
            prestoc.makedata_dt_set(self.this,value)
            return
        if name == "T" :
            prestoc.makedata_T_set(self.this,value)
            return
        if name == "ptype" :
            prestoc.makedata_ptype_set(self.this,value)
            return
        if name == "pnum" :
            prestoc.makedata_pnum_set(self.this,value)
            return
        if name == "fwhm" :
            prestoc.makedata_fwhm_set(self.this,value)
            return
        if name == "round" :
            prestoc.makedata_round_set(self.this,value)
            return
        if name == "roundnum" :
            prestoc.makedata_roundnum_set(self.this,value)
            return
        if name == "f" :
            prestoc.makedata_f_set(self.this,value)
            return
        if name == "fd" :
            prestoc.makedata_fd_set(self.this,value)
            return
        if name == "fdd" :
            prestoc.makedata_fdd_set(self.this,value)
            return
        if name == "p" :
            prestoc.makedata_p_set(self.this,value)
            return
        if name == "pd" :
            prestoc.makedata_pd_set(self.this,value)
            return
        if name == "pdd" :
            prestoc.makedata_pdd_set(self.this,value)
            return
        if name == "r" :
            prestoc.makedata_r_set(self.this,value)
            return
        if name == "z" :
            prestoc.makedata_z_set(self.this,value)
            return
        if name == "w" :
            prestoc.makedata_w_set(self.this,value)
            return
        if name == "amp" :
            prestoc.makedata_amp_set(self.this,value)
            return
        if name == "phs" :
            prestoc.makedata_phs_set(self.this,value)
            return
        if name == "dc" :
            prestoc.makedata_dc_set(self.this,value)
            return
        if name == "binary" :
            prestoc.makedata_binary_set(self.this,value)
            return
        if name == "orb" :
            prestoc.makedata_orb_set(self.this,value.this)
            return
        if name == "ampmod" :
            prestoc.makedata_ampmod_set(self.this,value)
            return
        if name == "ampmoda" :
            prestoc.makedata_ampmoda_set(self.this,value)
            return
        if name == "ampmodf" :
            prestoc.makedata_ampmodf_set(self.this,value)
            return
        if name == "ampmodp" :
            prestoc.makedata_ampmodp_set(self.this,value)
            return
        if name == "noisetype" :
            prestoc.makedata_noisetype_set(self.this,value)
            return
        if name == "noise" :
            prestoc.makedata_noise_set(self.this,value)
            return
        if name == "noisesig" :
            prestoc.makedata_noisesig_set(self.this,value)
            return
        if name == "numonoff" :
            prestoc.makedata_numonoff_set(self.this,value)
            return
        if name == "onoff" :
            prestoc.makedata_onoff_set(self.this,value.this)
            return
        self.__dict__[name] = value
    def __getattr__(self,name):
        if name == "basefilenm" : 
            return prestoc.makedata_basefilenm_get(self.this)
        if name == "description" : 
            return prestoc.makedata_description_get(self.this)
        if name == "N" : 
            return prestoc.makedata_N_get(self.this)
        if name == "next2_to_n" : 
            return prestoc.makedata_next2_to_n_get(self.this)
        if name == "dt" : 
            return prestoc.makedata_dt_get(self.this)
        if name == "T" : 
            return prestoc.makedata_T_get(self.this)
        if name == "ptype" : 
            return prestoc.makedata_ptype_get(self.this)
        if name == "pnum" : 
            return prestoc.makedata_pnum_get(self.this)
        if name == "fwhm" : 
            return prestoc.makedata_fwhm_get(self.this)
        if name == "round" : 
            return prestoc.makedata_round_get(self.this)
        if name == "roundnum" : 
            return prestoc.makedata_roundnum_get(self.this)
        if name == "f" : 
            return prestoc.makedata_f_get(self.this)
        if name == "fd" : 
            return prestoc.makedata_fd_get(self.this)
        if name == "fdd" : 
            return prestoc.makedata_fdd_get(self.this)
        if name == "p" : 
            return prestoc.makedata_p_get(self.this)
        if name == "pd" : 
            return prestoc.makedata_pd_get(self.this)
        if name == "pdd" : 
            return prestoc.makedata_pdd_get(self.this)
        if name == "r" : 
            return prestoc.makedata_r_get(self.this)
        if name == "z" : 
            return prestoc.makedata_z_get(self.this)
        if name == "w" : 
            return prestoc.makedata_w_get(self.this)
        if name == "amp" : 
            return prestoc.makedata_amp_get(self.this)
        if name == "phs" : 
            return prestoc.makedata_phs_get(self.this)
        if name == "dc" : 
            return prestoc.makedata_dc_get(self.this)
        if name == "binary" : 
            return prestoc.makedata_binary_get(self.this)
        if name == "orb" : 
            return orbitparamsPtr(prestoc.makedata_orb_get(self.this))
        if name == "ampmod" : 
            return prestoc.makedata_ampmod_get(self.this)
        if name == "ampmoda" : 
            return prestoc.makedata_ampmoda_get(self.this)
        if name == "ampmodf" : 
            return prestoc.makedata_ampmodf_get(self.this)
        if name == "ampmodp" : 
            return prestoc.makedata_ampmodp_get(self.this)
        if name == "noisetype" : 
            return prestoc.makedata_noisetype_get(self.this)
        if name == "noise" : 
            return prestoc.makedata_noise_get(self.this)
        if name == "noisesig" : 
            return prestoc.makedata_noisesig_get(self.this)
        if name == "numonoff" : 
            return prestoc.makedata_numonoff_get(self.this)
        if name == "onoff" : 
            return DoubleArrayPtr(prestoc.makedata_onoff_get(self.this))
        raise AttributeError,name
    def __repr__(self):
        return "<C makedata instance>"
class makedata(makedataPtr):
    def __init__(self) :
        self.this = prestoc.new_makedata()
        self.thisown = 1




class rderivsPtr :
    def __init__(self,this):
        self.this = this
        self.thisown = 0
    def __del__(self):
        if self.thisown == 1 :
            prestoc.delete_rderivs(self.this)
    def __setattr__(self,name,value):
        if name == "pow" :
            prestoc.rderivs_pow_set(self.this,value)
            return
        if name == "phs" :
            prestoc.rderivs_phs_set(self.this,value)
            return
        if name == "dpow" :
            prestoc.rderivs_dpow_set(self.this,value)
            return
        if name == "dphs" :
            prestoc.rderivs_dphs_set(self.this,value)
            return
        if name == "d2pow" :
            prestoc.rderivs_d2pow_set(self.this,value)
            return
        if name == "d2phs" :
            prestoc.rderivs_d2phs_set(self.this,value)
            return
        if name == "locpow" :
            prestoc.rderivs_locpow_set(self.this,value)
            return
        self.__dict__[name] = value
    def __getattr__(self,name):
        if name == "pow" : 
            return prestoc.rderivs_pow_get(self.this)
        if name == "phs" : 
            return prestoc.rderivs_phs_get(self.this)
        if name == "dpow" : 
            return prestoc.rderivs_dpow_get(self.this)
        if name == "dphs" : 
            return prestoc.rderivs_dphs_get(self.this)
        if name == "d2pow" : 
            return prestoc.rderivs_d2pow_get(self.this)
        if name == "d2phs" : 
            return prestoc.rderivs_d2phs_get(self.this)
        if name == "locpow" : 
            return prestoc.rderivs_locpow_get(self.this)
        raise AttributeError,name
    def __repr__(self):
        return "<C rderivs instance>"
class rderivs(rderivsPtr):
    def __init__(self) :
        self.this = prestoc.new_rderivs()
        self.thisown = 1




class fourierpropsPtr :
    def __init__(self,this):
        self.this = this
        self.thisown = 0
    def __del__(self):
        if self.thisown == 1 :
            prestoc.delete_fourierprops(self.this)
    def __setattr__(self,name,value):
        if name == "r" :
            prestoc.fourierprops_r_set(self.this,value)
            return
        if name == "rerr" :
            prestoc.fourierprops_rerr_set(self.this,value)
            return
        if name == "z" :
            prestoc.fourierprops_z_set(self.this,value)
            return
        if name == "zerr" :
            prestoc.fourierprops_zerr_set(self.this,value)
            return
        if name == "w" :
            prestoc.fourierprops_w_set(self.this,value)
            return
        if name == "werr" :
            prestoc.fourierprops_werr_set(self.this,value)
            return
        if name == "pow" :
            prestoc.fourierprops_pow_set(self.this,value)
            return
        if name == "powerr" :
            prestoc.fourierprops_powerr_set(self.this,value)
            return
        if name == "sig" :
            prestoc.fourierprops_sig_set(self.this,value)
            return
        if name == "rawpow" :
            prestoc.fourierprops_rawpow_set(self.this,value)
            return
        if name == "phs" :
            prestoc.fourierprops_phs_set(self.this,value)
            return
        if name == "phserr" :
            prestoc.fourierprops_phserr_set(self.this,value)
            return
        if name == "cen" :
            prestoc.fourierprops_cen_set(self.this,value)
            return
        if name == "cenerr" :
            prestoc.fourierprops_cenerr_set(self.this,value)
            return
        if name == "pur" :
            prestoc.fourierprops_pur_set(self.this,value)
            return
        if name == "purerr" :
            prestoc.fourierprops_purerr_set(self.this,value)
            return
        if name == "locpow" :
            prestoc.fourierprops_locpow_set(self.this,value)
            return
        self.__dict__[name] = value
    def __getattr__(self,name):
        if name == "r" : 
            return prestoc.fourierprops_r_get(self.this)
        if name == "rerr" : 
            return prestoc.fourierprops_rerr_get(self.this)
        if name == "z" : 
            return prestoc.fourierprops_z_get(self.this)
        if name == "zerr" : 
            return prestoc.fourierprops_zerr_get(self.this)
        if name == "w" : 
            return prestoc.fourierprops_w_get(self.this)
        if name == "werr" : 
            return prestoc.fourierprops_werr_get(self.this)
        if name == "pow" : 
            return prestoc.fourierprops_pow_get(self.this)
        if name == "powerr" : 
            return prestoc.fourierprops_powerr_get(self.this)
        if name == "sig" : 
            return prestoc.fourierprops_sig_get(self.this)
        if name == "rawpow" : 
            return prestoc.fourierprops_rawpow_get(self.this)
        if name == "phs" : 
            return prestoc.fourierprops_phs_get(self.this)
        if name == "phserr" : 
            return prestoc.fourierprops_phserr_get(self.this)
        if name == "cen" : 
            return prestoc.fourierprops_cen_get(self.this)
        if name == "cenerr" : 
            return prestoc.fourierprops_cenerr_get(self.this)
        if name == "pur" : 
            return prestoc.fourierprops_pur_get(self.this)
        if name == "purerr" : 
            return prestoc.fourierprops_purerr_get(self.this)
        if name == "locpow" : 
            return prestoc.fourierprops_locpow_get(self.this)
        raise AttributeError,name
    def __repr__(self):
        return "<C fourierprops instance>"
class fourierprops(fourierpropsPtr):
    def __init__(self) :
        self.this = prestoc.new_fourierprops()
        self.thisown = 1




class binarypropsPtr :
    def __init__(self,this):
        self.this = this
        self.thisown = 0
    def __del__(self):
        if self.thisown == 1 :
            prestoc.delete_binaryprops(self.this)
    def __setattr__(self,name,value):
        if name == "ppsr" :
            prestoc.binaryprops_ppsr_set(self.this,value)
            return
        if name == "fpsr" :
            prestoc.binaryprops_fpsr_set(self.this,value)
            return
        if name == "rpsr" :
            prestoc.binaryprops_rpsr_set(self.this,value)
            return
        if name == "pbin" :
            prestoc.binaryprops_pbin_set(self.this,value)
            return
        if name == "rbin" :
            prestoc.binaryprops_rbin_set(self.this,value)
            return
        if name == "z" :
            prestoc.binaryprops_z_set(self.this,value)
            return
        if name == "asinic" :
            prestoc.binaryprops_asinic_set(self.this,value)
            return
        if name == "rdetect" :
            prestoc.binaryprops_rdetect_set(self.this,value)
            return
        if name == "nfftbins" :
            prestoc.binaryprops_nfftbins_set(self.this,value)
            return
        if name == "lowbin" :
            prestoc.binaryprops_lowbin_set(self.this,value)
            return
        if name == "ppsrerr" :
            prestoc.binaryprops_ppsrerr_set(self.this,value)
            return
        if name == "fpsrerr" :
            prestoc.binaryprops_fpsrerr_set(self.this,value)
            return
        if name == "rpsrerr" :
            prestoc.binaryprops_rpsrerr_set(self.this,value)
            return
        if name == "pbinerr" :
            prestoc.binaryprops_pbinerr_set(self.this,value)
            return
        if name == "rbinerr" :
            prestoc.binaryprops_rbinerr_set(self.this,value)
            return
        if name == "zerr" :
            prestoc.binaryprops_zerr_set(self.this,value)
            return
        if name == "asinicerr" :
            prestoc.binaryprops_asinicerr_set(self.this,value)
            return
        if name == "rdetecterr" :
            prestoc.binaryprops_rdetecterr_set(self.this,value)
            return
        if name == "sig" :
            prestoc.binaryprops_sig_set(self.this,value)
            return
        if name == "phs" :
            prestoc.binaryprops_phs_set(self.this,value)
            return
        if name == "phserr" :
            prestoc.binaryprops_phserr_set(self.this,value)
            return
        if name == "cen" :
            prestoc.binaryprops_cen_set(self.this,value)
            return
        if name == "cenerr" :
            prestoc.binaryprops_cenerr_set(self.this,value)
            return
        if name == "pur" :
            prestoc.binaryprops_pur_set(self.this,value)
            return
        if name == "purerr" :
            prestoc.binaryprops_purerr_set(self.this,value)
            return
        if name == "pow" :
            prestoc.binaryprops_pow_set(self.this,value)
            return
        if name == "powerr" :
            prestoc.binaryprops_powerr_set(self.this,value)
            return
        self.__dict__[name] = value
    def __getattr__(self,name):
        if name == "ppsr" : 
            return prestoc.binaryprops_ppsr_get(self.this)
        if name == "fpsr" : 
            return prestoc.binaryprops_fpsr_get(self.this)
        if name == "rpsr" : 
            return prestoc.binaryprops_rpsr_get(self.this)
        if name == "pbin" : 
            return prestoc.binaryprops_pbin_get(self.this)
        if name == "rbin" : 
            return prestoc.binaryprops_rbin_get(self.this)
        if name == "z" : 
            return prestoc.binaryprops_z_get(self.this)
        if name == "asinic" : 
            return prestoc.binaryprops_asinic_get(self.this)
        if name == "rdetect" : 
            return prestoc.binaryprops_rdetect_get(self.this)
        if name == "nfftbins" : 
            return prestoc.binaryprops_nfftbins_get(self.this)
        if name == "lowbin" : 
            return prestoc.binaryprops_lowbin_get(self.this)
        if name == "ppsrerr" : 
            return prestoc.binaryprops_ppsrerr_get(self.this)
        if name == "fpsrerr" : 
            return prestoc.binaryprops_fpsrerr_get(self.this)
        if name == "rpsrerr" : 
            return prestoc.binaryprops_rpsrerr_get(self.this)
        if name == "pbinerr" : 
            return prestoc.binaryprops_pbinerr_get(self.this)
        if name == "rbinerr" : 
            return prestoc.binaryprops_rbinerr_get(self.this)
        if name == "zerr" : 
            return prestoc.binaryprops_zerr_get(self.this)
        if name == "asinicerr" : 
            return prestoc.binaryprops_asinicerr_get(self.this)
        if name == "rdetecterr" : 
            return prestoc.binaryprops_rdetecterr_get(self.this)
        if name == "sig" : 
            return prestoc.binaryprops_sig_get(self.this)
        if name == "phs" : 
            return prestoc.binaryprops_phs_get(self.this)
        if name == "phserr" : 
            return prestoc.binaryprops_phserr_get(self.this)
        if name == "cen" : 
            return prestoc.binaryprops_cen_get(self.this)
        if name == "cenerr" : 
            return prestoc.binaryprops_cenerr_get(self.this)
        if name == "pur" : 
            return prestoc.binaryprops_pur_get(self.this)
        if name == "purerr" : 
            return prestoc.binaryprops_purerr_get(self.this)
        if name == "pow" : 
            return prestoc.binaryprops_pow_get(self.this)
        if name == "powerr" : 
            return prestoc.binaryprops_powerr_get(self.this)
        raise AttributeError,name
    def __repr__(self):
        return "<C binaryprops instance>"
class binaryprops(binarypropsPtr):
    def __init__(self) :
        self.this = prestoc.new_binaryprops()
        self.thisown = 1




class rawbincandPtr :
    def __init__(self,this):
        self.this = this
        self.thisown = 0
    def __del__(self):
        if self.thisown == 1 :
            prestoc.delete_rawbincand(self.this)
    def __setattr__(self,name,value):
        if name == "full_N" :
            prestoc.rawbincand_full_N_set(self.this,value)
            return
        if name == "full_T" :
            prestoc.rawbincand_full_T_set(self.this,value)
            return
        if name == "full_lo_r" :
            prestoc.rawbincand_full_lo_r_set(self.this,value)
            return
        if name == "mini_N" :
            prestoc.rawbincand_mini_N_set(self.this,value)
            return
        if name == "mini_r" :
            prestoc.rawbincand_mini_r_set(self.this,value)
            return
        if name == "mini_power" :
            prestoc.rawbincand_mini_power_set(self.this,value)
            return
        if name == "mini_numsum" :
            prestoc.rawbincand_mini_numsum_set(self.this,value)
            return
        if name == "mini_sigma" :
            prestoc.rawbincand_mini_sigma_set(self.this,value)
            return
        if name == "psr_p" :
            prestoc.rawbincand_psr_p_set(self.this,value)
            return
        if name == "orb_p" :
            prestoc.rawbincand_orb_p_set(self.this,value)
            return
        self.__dict__[name] = value
    def __getattr__(self,name):
        if name == "full_N" : 
            return prestoc.rawbincand_full_N_get(self.this)
        if name == "full_T" : 
            return prestoc.rawbincand_full_T_get(self.this)
        if name == "full_lo_r" : 
            return prestoc.rawbincand_full_lo_r_get(self.this)
        if name == "mini_N" : 
            return prestoc.rawbincand_mini_N_get(self.this)
        if name == "mini_r" : 
            return prestoc.rawbincand_mini_r_get(self.this)
        if name == "mini_power" : 
            return prestoc.rawbincand_mini_power_get(self.this)
        if name == "mini_numsum" : 
            return prestoc.rawbincand_mini_numsum_get(self.this)
        if name == "mini_sigma" : 
            return prestoc.rawbincand_mini_sigma_get(self.this)
        if name == "psr_p" : 
            return prestoc.rawbincand_psr_p_get(self.this)
        if name == "orb_p" : 
            return prestoc.rawbincand_orb_p_get(self.this)
        raise AttributeError,name
    def __repr__(self):
        return "<C rawbincand instance>"
class rawbincand(rawbincandPtr):
    def __init__(self) :
        self.this = prestoc.new_rawbincand()
        self.thisown = 1






#-------------- FUNCTION WRAPPERS ------------------

tofloatvector = prestoc.tofloatvector

float_to_complex = prestoc.float_to_complex

complex_to_float = prestoc.complex_to_float

power_arr = prestoc.power_arr

phase_arr = prestoc.phase_arr

dpower_arr = prestoc.dpower_arr

dphase_arr = prestoc.dphase_arr

num_psrs_in_database = prestoc.num_psrs_in_database

def get_psrdata(arg0,arg1):
    val = prestoc.get_psrdata(arg0.this,arg1)
    return val

def get_psrdata_by_num(arg0,arg1):
    val = prestoc.get_psrdata_by_num(arg0.this,arg1)
    return val

def return_psrparams_at_epoch(arg0,arg1,arg2):
    val = prestoc.return_psrparams_at_epoch(arg0.this,arg1,arg2)
    return val

def get_psr_at_epoch(arg0,arg1,arg2,arg3):
    val = prestoc.get_psr_at_epoch(arg0,arg1,arg2,arg3.this)
    return val

def readinf(arg0,arg1):
    val = prestoc.readinf(arg0.this,arg1)
    return val

def writeinf(arg0):
    val = prestoc.writeinf(arg0.this)
    return val

def read_mak_input(arg0):
    val = prestoc.read_mak_input(arg0.this)
    return val

def read_mak_file(arg0,arg1):
    val = prestoc.read_mak_file(arg0,arg1.this)
    return val

def write_mak_file(arg0):
    val = prestoc.write_mak_file(arg0.this)
    return val

frotate = prestoc.frotate

drotate = prestoc.drotate

def dorbint(arg0,arg1,arg2,arg3):
    val = prestoc.dorbint(arg0,arg1,arg2,arg3.this)
    return val

keplars_eqn = prestoc.keplars_eqn

lin_interp_E = prestoc.lin_interp_E

def E_to_phib(arg0,arg1,arg2):
    val = prestoc.E_to_phib(arg0,arg1,arg2.this)
    return val

def E_to_v(arg0,arg1,arg2):
    val = prestoc.E_to_v(arg0,arg1,arg2.this)
    return val

def E_to_p(arg0,arg1,arg2,arg3):
    val = prestoc.E_to_p(arg0,arg1,arg2,arg3.this)
    return val

def E_to_z(arg0,arg1,arg2,arg3,arg4):
    val = prestoc.E_to_z(arg0,arg1,arg2,arg3,arg4.this)
    return val

def E_to_phib_BT(arg0,arg1,arg2):
    val = prestoc.E_to_phib_BT(arg0,arg1,arg2.this)
    return val

r_resp_halfwidth = prestoc.r_resp_halfwidth

z_resp_halfwidth = prestoc.z_resp_halfwidth

w_resp_halfwidth = prestoc.w_resp_halfwidth

def bin_resp_halfwidth(arg0,arg1,arg2):
    val = prestoc.bin_resp_halfwidth(arg0,arg1,arg2.this)
    return val

gen_r_response = prestoc.gen_r_response

gen_z_response = prestoc.gen_z_response

gen_w_response = prestoc.gen_w_response

def gen_bin_response(arg0,arg1,arg2,arg3,arg4,arg5):
    val = prestoc.gen_bin_response(arg0,arg1,arg2,arg3,arg4.this,arg5)
    return val

get_numphotons = prestoc.get_numphotons

get_localpower = prestoc.get_localpower

get_localpower3d = prestoc.get_localpower3d

def get_derivs3d(arg0,arg1,arg2,arg3,arg4,arg5,arg6):
    val = prestoc.get_derivs3d(arg0,arg1,arg2,arg3,arg4,arg5,arg6.this)
    return val

def calc_props(arg0,arg1,arg2,arg3,arg4):
    val = prestoc.calc_props(arg0.this,arg1,arg2,arg3,arg4.this)
    return val

def calc_binprops(arg0,arg1,arg2,arg3,arg4):
    val = prestoc.calc_binprops(arg0.this,arg1,arg2,arg3,arg4.this)
    return val

def calc_rzwerrs(arg0,arg1,arg2):
    val = prestoc.calc_rzwerrs(arg0.this,arg1,arg2)
    return val

sigma_from_sumpows = prestoc.sigma_from_sumpows

sumpows_from_sigma = prestoc.sumpows_from_sigma

nice_output_1 = prestoc.nice_output_1

nice_output_2 = prestoc.nice_output_2

def print_candidate(arg0,arg1,arg2,arg3,arg4):
    val = prestoc.print_candidate(arg0.this,arg1,arg2,arg3,arg4)
    return val

def print_bin_candidate(arg0,arg1):
    val = prestoc.print_bin_candidate(arg0.this,arg1)
    return val

def read_rzw_cand(arg0,arg1):
    val = prestoc.read_rzw_cand(arg0,arg1.this)
    return val

def get_rzw_cand(arg0,arg1,arg2):
    val = prestoc.get_rzw_cand(arg0,arg1,arg2.this)
    return val

def get_bin_cand(arg0,arg1,arg2):
    val = prestoc.get_bin_cand(arg0,arg1,arg2.this)
    return val

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

def max_r_arr(arg0,arg1,arg2,arg3):
    val = prestoc.max_r_arr(arg0,arg1,arg2,arg3.this)
    return val

def max_rz_arr(arg0,arg1,arg2,arg3,arg4):
    val = prestoc.max_rz_arr(arg0,arg1,arg2,arg3,arg4.this)
    return val

foldfile = prestoc.foldfile

fold = prestoc.fold

simplefold = prestoc.simplefold

doppler = prestoc.doppler

def search_minifft(arg0,arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10):
    val = prestoc.search_minifft(arg0,arg1,arg2.this,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10)
    return val

def print_rawbincand(arg0):
    val = prestoc.print_rawbincand(arg0.this)
    return val

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

import math, Numeric, Pgplot

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


   
