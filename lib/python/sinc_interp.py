import Numeric, umath, scipy.special, FFT

def sinc(xs):
    """
    sinc(xs):
        Return the sinc function [i.e. sin(pi * xs)/(pi * xs)]
            for the values xs.
    """
    pxs = umath.pi*xs
    return Numeric.where(umath.fabs(pxs)<1e-3, 1.0-pxs*pxs/6.0, umath.sin(pxs)/pxs)

def kaiser_window(xs, halfwidth, alpha):
    """
    kaiser_window(xs, halfwidth, alpha):
        Return the kaiser window function for the values 'xs' when the
            the half-width of the window should be 'haldwidth' with
            the folloff parameter 'alpha'.  The following values are
            particularly interesting:

            alpha
            -----
            0           Rectangular Window
            5           Similar to Hamming window
            6           Similar to Hanning window
            8.6         Almost identical to the Blackman window 
    """
    win = i0(alpha*umath.sqrt(1.0-(xs/halfwidth)**2.0))/i0(alpha)
    return Numeric.where(umath.fabs(xs)<=halfwidth, win, 0.0)

def hanning_window(xs, halfwidth):
    """
    hanning_window(xs, halfwidth):
        Return the Hanning window of halfwidth 'halfwidth' evaluated at
            the values 'xs'.
    """
    win =  0.5 + 0.5*umath.cos(umath.pi*xs/halfwidth)
    return Numeric.where(umath.fabs(xs)<=halfwidth, win, 0.0)

def hamming_window(xs, halfwidth):
    """
    hamming_window(xs, halfwidth):
        Return the Hamming window of halfwidth 'halfwidth' evaluated at
            the values 'xs'.
    """
    win =  0.54 + 0.46*umath.cos(umath.pi*xs/halfwidth)
    return Numeric.where(umath.fabs(xs)<=halfwidth, win, 0.0)

def blackman_window(xs, halfwidth):
    """
    blackman_window(xs, halfwidth):
        Return the Blackman window of halfwidth 'halfwidth' evaluated at
            the values 'xs'.
    """
    rat = umath.pi*xs/halfwidth
    win =  0.42 + 0.5*umath.cos(rat) + 0.08*umath.cos(2.0*rat) 
    return Numeric.where(umath.fabs(xs)<=halfwidth, win, 0.0)

def rectangular_window(xs, halfwidth):
    """
    rectangular_window(xs, halfwidth):
        Return a rectangular window of halfwidth 'halfwidth' evaluated at
            the values 'xs'.
    """
    return Numeric.where(umath.fabs(xs)<=halfwidth, 1.0, 0.0)
    
_window_function = {"rectangular": rectangular_window,
                    "none": rectangular_window,
                    "hanning": hanning_window,
                    "hamming": hamming_window,
                    "blackman": blackman_window,
                    "kaiser": kaiser_window}

def windowed_sinc_interp(data, newx, halfwidth=None,
                         window='hanning', alpha=6.0):
    """
    windowed_sinc_interp(data, newx, halfwidth=None,
                         window='hanning', alpha=6.0):
        Return a single windowed-sinc-interpolated point from the data.
    """
    if umath.fabs(round(newx)-newx) < 1e-5:
        return data[int(round(newx))]
    num_pts = (int(umath.floor(newx)), len(data)-int(umath.ceil(newx))-1)
    if halfwidth is None:
        halfwidth = min(num_pts)
    lo_pt = int(umath.floor(newx)) - halfwidth
    if lo_pt < 0:
        lo_pt < 0
        print "Warning:  trying to access below the lowest index!"
    hi_pt = lo_pt + 2*halfwidth
    if hi_pt >= len(data):
        hi_pt = len(data)-1
        print "Warning:  trying to access above the highest index!"
    halfwidth = (hi_pt-lo_pt)/2
    pts = Numeric.arange(2*halfwidth)+lo_pt
    xs = newx - pts
    if window.lower() is "kaiser":
        win = _window_function[window](xs, len(data)/2, alpha)
    else:
        win = _window_function[window](xs, len(data)/2)
    return umath.add.reduce(Numeric.take(data, pts) * win * sinc(xs))

def periodic_interp(data, zoomfact, window='hanning', alpha=6.0):
    """
    periodic_interp(data, zoomfact, window='hanning', alpha=6.0):
        Return a periodic, windowed, sinc-interpolation of the data which
            is oversampled by a factor of 'zoomfact'.
    """
    zoomfact = int(zoomfact)
    if (zoomfact < 1):
        print "zoomfact must be >= 1."
        return 0.0
    elif zoomfact==1:
        return data
    newN = len(data)*zoomfact
    # Space out the data
    comb = Numeric.zeros((zoomfact, len(data)), typecode='d')
    comb[0] += data
    comb = Numeric.reshape(Numeric.transpose(comb), (newN,))
    # Compute the offsets
    xs = Numeric.zeros(newN, typecode='d')
    xs[:newN/2+1] = Numeric.arange(newN/2+1, typecode='d')/zoomfact
    xs[-newN/2:]  = xs[::-1][newN/2-1:-1]
    # Calculate the sinc times window for the kernel
    if window.lower()=="kaiser":
        win = _window_function[window](xs, len(data)/2, alpha)
    else:
        win = _window_function[window](xs, len(data)/2)
    kernel = win * sinc(xs)
    if (0):
        plotxy(sinc(xs), color='yellow')
        plotxy(win)
        plotxy(kernel, color='red')
        closeplot()
    return FFT.inverse_real_fft(FFT.real_fft(kernel) * FFT.real_fft(comb))


if __name__=='__main__':
    from psr_utils import *
    from Pgplot import *
    from RandomArray import normal
    #from scipy import interpolate
    from spline import *
    
    fwhm = 0.01
    ctr_phase = 0.505
    noise_sigma = 0.2

    # The theoretical profile with noise
    Ntheo = 1000
    theo = gaussian_profile(Ntheo, ctr_phase, fwhm) + normal(0.0, noise_sigma, Ntheo)
    theo_phases = Numeric.arange(Ntheo, typecode='d')/Ntheo

    # The "sampled" data
    Ndata = 100
    data = theo[::Ntheo/Ndata]
    data_phases = theo_phases[::Ntheo/Ndata]

    # The values to interpolate
    Ncalc = 30
    lo_calc = ctr_phase-0.05
    hi_calc = ctr_phase+0.05
    calc_phases = span(lo_calc, hi_calc, Ncalc)
    plotxy(theo, theo_phases, rangex=[lo_calc-0.2, hi_calc+0.2])
    plotxy(data, data_phases, line=None, symbol=3, color='green')

    # Do the interpolation one point at a time
    halfwidth = Ndata/2-5
    calc_vals = []
    for phs in calc_phases:
        calc_vals.append(windowed_sinc_interp(data, phs*len(data), halfwidth))
    plotxy(calc_vals, calc_phases, line=None, symbol=3, color='red')

    # Interpolate the full profile using convolution
    zoomfact = 10
    newvals = periodic_interp(data, 10)
    new_phases = Numeric.arange(Ndata*zoomfact, typecode='d')/(Ndata*zoomfact)
    plotxy(newvals, new_phases, line=1, symbol=None, color='yellow')

    # Interpolate using cubic splines
    if (0):
        sdata = interpolate.splrep(data, data_phases, s=0)
        svals = interpolate.splrep(new_phases, sdata, der=0)
        plotxy(svals, new_phases, line=1, symbol=None, color='cyan')
    else:
        sdata = Spline(data_phases, data)
        svals = sdata(new_phases)
        plotxy(svals, new_phases, line=1, symbol=None, color='cyan')
    
    closeplot()
           
