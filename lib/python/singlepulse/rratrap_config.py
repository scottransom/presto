#! /usr/bin/env python

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

def use_ddplan(dm):
    """ Sets a factor which multiplies the DMthreshold and time_threshold. The factor is 
        the downsampling rate. 
        This makes the DM_THRESH and TIME_THRESH depend on the DM instead of having fixed 
        values throughout. This helps at higher DMs where the DM step size is > 0.5 pc cm-3. 
        This is specific to PALFA. You can make a similar plan depending on the survey's DDplan.
    """
    if (dm <=212.8):
        dmt = 1
        min_group = 45
    elif (dm >212.8) and (dm <=443.2):
        dmt = 2
        min_group = 40
    elif (dm >443.2) and (dm <=543.4):
        dmt = 3
        min_group = 35
    elif (dm >543.4) and (dm <=876.4):
        dmt = 5
        min_group = 30
    elif (dm >876.4) and (dm <=990.4):
        dmt = 6
        min_group = 30
    else:
        dmt = 10
        min_group = 30

    return dmt, min_group
