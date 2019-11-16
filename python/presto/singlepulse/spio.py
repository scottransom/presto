from __future__ import print_function
from builtins import range
import numpy as _np
import fileinput as _fileinput
from scipy.special import erf


def get_textfile(txtfile):
    """ Read in the groups.txt file.
    Contains information about the DM, time, box car width, signal to noise, sample number and rank    of groups. 
    Input: name of the text file.
    Output: Array of strings in the file line by line.
    """
    return  _np.loadtxt(txtfile,dtype = 'str',delimiter='\n')


def group_info(rank, txtfile):
    """
    Extracts out relevant information from the groups.txt file as strings. 
    Input: integer:rank of the groups to read.
           string :name of the text file to read from.
    Output: 
           list of single pulse parameters as strings.
    """
    files = get_textfile(txtfile)
    lis=_np.where(files == '\tRank:             %i.000000'%rank)[0]#Checks for this contidion and gives its indices where true.
    # Extract the Max_ sigma value for the required parameters
    parameters=[]
    for i in range(len(lis)):
        temp_list = files[lis[i]-1].split()
        max_sigma = temp_list[2]
        max_sigma = float(max_sigma)
        max_sigma = '%.2f'%max_sigma
        # Extract the number of pulses for this group
        temp_list = files[lis[i]-6].split()
        number_of_pulses = int(temp_list[2])
        # Slice off a mini array to get the parameters from
        temp_lines = files[(lis[i]+1):(lis[i]+number_of_pulses+1)]
        # Get the parameters as strings containing the max_sigma
        parameters.append(temp_lines[_np.array([max_sigma in line for line in temp_lines])])
    return parameters

def split_parameters(rank, txtfile):
    """
    Splits the string into individual parameters and converts them into floats/int. 
    Input: integer:rank of the groups to read.
           string :name of the text file to read from.
    Output: 
           list of single pulse parameters as floats and integers.
    """
    parameters = group_info(rank, txtfile)
    final_parameters=[]
    for i in range(len(parameters)):
    # If there is a degeneracy in max_sigma values, Picks the first one.(Can be updated to get the best pick) 
        correct_values = parameters[i][0].split()
        correct_values[0] = float(correct_values[0])
        correct_values[1] = float(correct_values[1])
        correct_values[1] = float('%.2f'%correct_values[1])
        correct_values[2] = float(correct_values[2])
        correct_values[3] = int(correct_values[3])
        correct_values[4] = int(correct_values[4])
        final_parameters.append(correct_values)
    return final_parameters

def read_RRATrap_info(groupfile, group_to_read, rank):
    """
    Generates DM, time, signal-to-noise and pulse width arrays for all the single pulse events
    that belong to this group.
    Inputs:
        groups.txt file.
        list of line numbers in the file that read something like : Rank:             6.000000
        rank of the group according to the RRATrap classification.
    Outputs:
        dm_list, dm_arr : DMs of the single pulse events in this group
        time_list : Times of the single pulse events in this group 
        sigma_arr : signal-to-noise of the single pulse events in this group 
        width_arr : box car widths of the single pulse events in this group 
    """
    files = get_textfile(groupfile)
    if files[rank-1] != "Number of rank %i groups: 0 "%rank: # checks whether there are %i \
                                                             # ranked groups in the file.
        print(files[rank-1])
        print("Making arrays for DM vs Signal to Noise...")
        temp_list = files[group_to_read-6].split()
        npulses = int(temp_list[2])
        temp_lines = files[(group_to_read+3):(group_to_read+npulses+1)]
        arr = _np.split(temp_lines, len(temp_lines))
        dm_list = []
        time_list = []
        for i in range(len(arr)):
            dm_val= float(arr[i][0].split()[0])
            time_val = float(arr[i][0].split()[2])
            dm_list.append(dm_val)
            time_list.append(time_val)
        arr_2 = _np.array([arr[i][0].split() for i in range(len(arr))], dtype = _np.float32)
        dm_arr = _np.array([arr_2[i][0] for i in range(len(arr))], dtype = _np.float32)
        sigma_arr = _np.array([arr_2[i][1] for i in range(len(arr))], dtype = _np.float32)
        width_arr = _np.array([arr_2[i][4] for i in range(len(arr))], dtype = _np.int8)
    else:
        print("No Rank %i groups for this candidate."%rank)
    
    return dm_list, time_list, dm_arr, sigma_arr, width_arr


def old_read_sp_files(sp_files):
    """*** OLD VERSION ***
        Read all *.singlepulse files in the current directory.
	Return 5 arrays (properties of all single pulses):
		DM, sigma, time, sample, downfact.
    """
    tmp_sp_params = _np.array(_np.empty((1,0)), 
                             dtype=_np.dtype([('dm', 'float64'),
                                             ('sigma','float32'),
                                             ('time','float64'),
                                             ('sample','uint32'),
                                             ('downfact','uint8')]))

    for file in sp_files:
       if os.path.getsize(file):
           curr = _np.atleast_2d(_np.loadtxt(file, dtype=_np.dtype([('dm', 'float64'),('sigma','float32'),('time','float64'),('sample','uint32'),('downfact','uint8')])))
           tmp_sp_params = _np.concatenate([tmp_sp_params, curr], axis=1)

    return tmp_sp_params

def read_sp_files(files):
    """Read all *.singlepulse files in the current directory in a DM range.
        Return 5 arrays (properties of all single pulses):
                DM, sigma, time, sample, downfact."""
    finput = _fileinput.input(files)
    data = _np.loadtxt(finput,
                       dtype=_np.dtype([('dm', 'float32'),
                                        ('sigma','float32'),
                                        ('time','float32'),
                                        ('sample','uint32'),
                                        ('downfact','uint8')]))
    return _np.atleast_2d(data)

def read_tarfile(filenames, names, tar):
    """Read in the .singlepulse.tgz file instead of individual .singlepulse files.
        Return an array of (properties of all single pulses):
              DM, sigma, time, sample, downfact. 
        Input: filenames: names of all the singlepulse files.
               names: subset of filenames. Names of the singlepulse files to be 
               plotted in DM vs time.
               tar: tar file (.singlepulse.tgz)."""  
    members = []
    for name in names:
        if name in filenames:
            member = tar.getmember(name)
            members.append(member)
        else:
            pass
    fileinfo = []
    filearr = []
    for mem in members:
        file = tar.extractfile(mem)
        for line in file.readlines():
            fileinfo.append(line)
        filearr+=(fileinfo[1:])  #Removes the text labels ("DM", "sigma" etc) of the singlepulse properties. Only keeps the values. 
        fileinfo = []
    temp_list = []
    for i in range(len(filearr)):
        temp_line = filearr[i].split()
        temp_list.append(temp_line)
    main_array = _np.asarray(temp_list)
    main_array = _np.split(main_array, 5, axis=1)
    main_array[0] = main_array[0].astype(_np.float16)
    main_array[1] = main_array[1].astype(_np.float16)
    main_array[2] = main_array[2].astype(_np.float16)
    main_array[3] = main_array[3].astype(_np.int)
    main_array[4] = main_array[4].astype(_np.int)
    return main_array

def pick_DM_for_singlepulse_files(filenm):
    """
       Extract DM values from filenames from .singlepulse files.
       Input: string : .singlepulse file mane 
       Output: float: DM value of that file.
    """
    return float(filenm[filenm.find('DM')+2:filenm.find('.singlepulse')])

def ddm_response(ddm, width_ms, lofreq, hifreq):
    if _np.isscalar(ddm):
        ddm = _np.array([ddm])
        scal = True
    else:
        ddm = _np.array([ddm])
        scal = False
    band_MHz = _np.array((lofreq, hifreq))
    zeta = 6.91e-3 * ddm * _np.diff(band_MHz)[0] / (width_ms * (_np.mean(band_MHz)/1000.)**3)
    result = _np.zeros_like(ddm)
    where_nonzero = _np.where(zeta != 0)
    result[where_nonzero] = 0.5*_np.sqrt(_np.pi)*erf(zeta[where_nonzero])/zeta[where_nonzero]
    result[zeta == 0] = 1.
    if scal: return result[0]
    else: return result

def theoritical_dmspan(maxsigma, minsigma, width_ms, lofreq, hifreq):
    # since the sigma threshold = 5
    sigma_limit = minsigma/maxsigma
    # spans over a dm range of 1000 (500*2)  
    ddm = _np.linspace(0, 5000, 50001)
    # makes a normalized gaussian of sigma values
    sigma_range = ddm_response(ddm, width_ms, lofreq, hifreq)
    # Returns te index where sigma_limit is closest to one of the values in sigma_range
    ind = (_np.abs(sigma_range-sigma_limit)).argmin()
    return 2*ddm[ind]

def gen_arrays(dm, sp_files, tar, threshold):    
    """
    Extract dms, times and signal to noise from each singlepulse file as 1D arrays.
    Input: 
           dm: The dm array of the main pulse. Used to decide the DM range in the DM vs time plot and pick out singlepulse files with those DMs.
           threshold: Min signal to noise of the single pulse event that is plotted.
           sp_files: all the .singlepulse file names.
           tar: Instead of the providing individual singlepulse files, you can provide the .singlepulse.tgz tarball.
    Output:
           Arrays: dms, times, sigmas of the singlepulse events and an array of dm_vs_times file names.
           
    Options: Either a tarball of singlepulse files or individual singlepulse files can be supplied.
             Faster when individual singlepulse files are supplied.   
    """
    max_dm = _np.ceil(_np.max(dm)).astype('int')
    min_dm = _np.min(dm).astype('int')
    diff_dm = max_dm-min_dm
    ddm = min_dm-diff_dm
    hidm = max_dm+diff_dm
    if (ddm <= 0):
        ddm = 0
    name_DMs = _np.asarray([pick_DM_for_singlepulse_files(sp_files[x]) for x in range(len(sp_files))])
    inds = name_DMs.argsort()
    name_DMs = name_DMs[inds]
    sp_files = _np.asarray(sp_files)[inds]
    loidx = _np.argmin(_np.abs(name_DMs-ddm))
    hiidx = _np.argmin(_np.abs(name_DMs-hidm))
    print(loidx, hiidx)
    singlepulsefiles = list(sp_files[loidx:hiidx])

    if tar is not None:
        data = read_tarfile(sp_files, singlepulsefiles, tar)
        dms = _np.reshape(data[0],(len(data[0]),))
        times = _np.reshape(data[2],(len(data[1]),))
        sigmas = _np.reshape(data[1],(len(data[2]),))
        widths = _np.reshape(data[4],(len(data[4]),))
    else:
        data = read_sp_files(singlepulsefiles)[0]
        dms = data['dm']
        times = data['time']
        sigmas = data['sigma']
        widths = data['downfact']
    
    dms = _np.delete(dms, (0), axis = 0)
    times = _np.delete(times, (0), axis = 0)
    sigmas = _np.delete(sigmas, (0), axis = 0)
    widths = _np.delete(widths, (0), axis = 0)
    return dms, times, sigmas, widths, singlepulsefiles

def read_spd(spd_file, tar = None):
    """ 
       Reads in all the .spd and the .singlepulse.tgz info that can reproduce the sp plots.
       Inputs: spd_file: .spd file
               .singlepulse.tgz: if not supplied, it will only output .spd info. 
                                 Default: not supplied. 
       Output: An object that has all the relevant information to remake the plot. 
    """
    sp = spd(spd_file)
    if tar is not None:
        dmVt_dms, dmVt_times, dmVt_sigmas, dmVt_widths, dmVt_files = gen_arrays(sp.dmVt_this_dms, sp.spfiles, tar, threshold=5)
        sp.dmVt_dms = dmVt_dms
        sp.dmVt_times = dmVt_times
        sp.dmVt_sigmas = dmVt_sigmas
        sp.dmVt_widths = dmVt_widths
        return sp
    else:
        return sp 


