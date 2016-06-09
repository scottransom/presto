#! /usr/bin/env python

"""
   Original code found in presto. Written by Scott M. Ransom.
   Modified to return topocentric and corresponding barycentric
   times. 
"""
from presto.prestoswig import *
import numpy as Num
import Pgplot
import psr_utils
import glob

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

def bary_to_topo(infofilenm, ephem="DE200"):
   """
   bary_to_topo(infofilenm, ephem="DE200"):
      Returns the barycentric and topocentric times evert 10 seconds.
      The data for the observation must be found in the info file.
   """
   if infofilenm[-4:]==".inf":  infofilenm = infofilenm[:-4]
   obs = read_inffile(infofilenm)
   T = obs.N * obs.dt
   dt = 10.0
   tto = obs.mjd_i + obs.mjd_f
   tts = Num.arange(tto, tto + (T + dt) / SECPERDAY, dt / SECPERDAY)
   nn = len(tts)
   bts = Num.zeros(nn, 'd')
   vel = Num.zeros(nn, 'd')
   ra = psr_utils.coord_to_string(obs.ra_h, obs.ra_m, obs.ra_s)
   dec = psr_utils.coord_to_string(obs.dec_d, obs.dec_m, obs.dec_s)
   if (obs.telescope == 'Parkes'):  tel = 'PK'
   elif (obs.telescope == 'Effelsberg'):  tel = 'EB'
   elif (obs.telescope == 'Arecibo'):  tel = 'AO'
   elif (obs.telescope == 'MMT'):  tel = 'MT'
   elif (obs.telescope == 'GBT'):  tel = 'GB'
   else:
      print "Telescope not recognized."
      return 0
   barycenter(tts, bts, vel, nn, ra, dec, tel, ephem)
   avgvel = Num.add.reduce(vel) / nn
   tts = Num.arange(nn, dtype='d') * dt
   bts = (bts - bts[0]) * SECPERDAY
   return tts, bts
