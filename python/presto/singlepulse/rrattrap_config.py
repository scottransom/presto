#! /usr/bin/env python
from builtins import range

CLOSE_DM = 2 # pc cm-3
# MIN_GROUP, DM_THRESH, TIME_THRESH will change later on depending on the DDplan.
MIN_GROUP = 50 #minimum group size that is not considered noise
TIME_THRESH = 0.1
DM_THRESH = 0.5
MIN_SIGMA = 8.0
PLOT = True
PLOTTYPE = 'pgplot' # 'pgplot' or 'matplotlib'
RANKS_TO_WRITE = [2,0,3,4,5,6]
RANKS_TO_PLOT = [2,3,4,5,6]

########This is specific to PALFA. You can make a similar plan depending 
########on the survey's DDplan. The values are the lodm, factor that multiplies the 
######## dm_thresh and time_thresh (I recommend choosing the downsampling factor) and
######## the minimum group size.
######### lodm, dmt(factor), min_group  
DMPLAN = [(0.0, 1, 45),             
          (212.8, 2, 40),
          (443.2, 3, 35),
          (534.4, 5, 30),
          (876.4, 6, 30),
          (990.4, 10, 30)]


def use_dmplan(dm):
    """ Sets a factor which multiplies the DMthreshold and time_threshold. The factor is 
        the downsampling rate. 
        This makes the DM_THRESH and TIME_THRESH depend on the DM instead of having fixed 
        values throughout. This helps at higher DMs where the DM step size is > 0.5 pc cm-3. 
    """
    for i in range(len(DMPLAN)-1):
        if dm<=DMPLAN[i+1][0] and dm>DMPLAN[i][0]:
            dmt, min_group = DMPLAN[i][1], DMPLAN[i][2]
        if dm > DMPLAN[i+1][0]:
            dmt, min_group = DMPLAN[i+1][1], DMPLAN[i+1][2]
    return dmt, min_group
