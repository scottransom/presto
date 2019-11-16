#!/usr/bin/env python
from __future__ import print_function
from builtins import range
from builtins import object
import os
import sys
from presto.psr_utils import gaussian_profile, read_profile
from matplotlib.patches import Rectangle
from presto.bestprof import bestprof
import matplotlib.pyplot as plt
import numpy as Num
from presto import mpfit


class GaussianSelector(object):
    def __init__(self, ax, profile, errs, profnm, minspanx=None,
                 minspany=None, useblit=True):
        self.ax = ax.axes
        self.profile = profile
        self.proflen = len(profile)
        self.profnm = profnm
        self.phases = Num.arange(self.proflen, dtype='d')/self.proflen
        self.errs = errs
        self.visible = True
        self.DCguess = sorted(profile)[len(profile) // 10 + 1]
        self.init_params = [self.DCguess]
        self.numgaussians = 0
        self.canvas = ax.figure.canvas
        self.canvas.mpl_connect('motion_notify_event', self.onmove)
        self.canvas.mpl_connect('button_press_event', self.press)
        self.canvas.mpl_connect('button_release_event', self.release)
        self.canvas.mpl_connect('draw_event', self.update_background)
        self.background = None
        self.rectprops = dict(facecolor='white', edgecolor = 'black',
                              alpha=0.5, fill=False)
        self.to_draw = Rectangle((0,0), 0, 1, visible=False, **self.rectprops)
        self.ax.add_patch(self.to_draw)
        self.useblit = useblit
        self.minspanx = minspanx
        self.minspany = minspany
        # will save the data (position at mouseclick)
        self.eventpress = None            
        # will save the data (pos. at mouserelease)
        self.eventrelease = None          
        self.plot_gaussians(self.init_params)

    def update_background(self, event):
        'force an update of the background'
        if self.useblit:
            self.background = self.canvas.copy_from_bbox(self.ax.bbox)
        
    def ignore(self, event):
        'return True if event should be ignored'
        # If no button was pressed yet ignore the event if it was out
        # of the axes
        if self.eventpress == None:
            return event.inaxes!= self.ax       
        # If a button was pressed, check if the release-button is the
        # same.
        return (event.inaxes!=self.ax or
                event.button != self.eventpress.button)
      
    def press(self, event):
        'on button press event'
        # Is the correct button pressed within the correct axes?
        if self.ignore(event): return         
        # make the drawed box/line visible get the click-coordinates,
        # button, ...
        self.eventpress = event               
        if event.button==1:
            self.to_draw.set_visible(self.visible)
            self.eventpress.ydata = self.DCguess

    def release(self, event):
        'on button release event'
        if self.eventpress is None or self.ignore(event): return
        # release coordinates, button, ...
        self.eventrelease = event             
        if event.button==1:
            # make the box/line invisible again
            self.to_draw.set_visible(False)       
            self.canvas.draw()
            xmin, ymin = self.eventpress.xdata, self.eventpress.ydata
            xmax, ymax = self.eventrelease.xdata, self.eventrelease.ydata
            # calculate dimensions of box 
            if xmin>xmax: xmin, xmax = xmax, xmin        
            if ymin>ymax: ymin, ymax = ymax, ymin
            spanx = xmax - xmin                  
            spany = ymax - ymin
            xproblems = self.minspanx is not None and spanx<self.minspanx
            yproblems = self.minspany is not None and spany<self.minspany
        # call desired function
        self.onselect()
        self.eventpress = None                # reset the variables to their
        self.eventrelease = None              #   inital values

    def update(self):
        'draw using newfangled blit or oldfangled draw depending on useblit'
        if self.useblit:
            if self.background is not None:
                self.canvas.restore_region(self.background)
            self.ax.draw_artist(self.to_draw)
            self.canvas.blit(self.ax.bbox)
        else:
            self.canvas.draw_idle()

    def onmove(self, event):
        if self.eventpress is None or self.ignore(event): return
        x,y = event.xdata, event.ydata # actual position with button still pressed
        minx, maxx = self.eventpress.xdata, x # click-x and actual mouse-x
        miny, maxy = self.eventpress.ydata, y # click-y and actual mouse-y
        if minx>maxx: minx, maxx = maxx, minx # get them in the right order
        if miny>maxy: miny, maxy = maxy, miny
        self.to_draw.set_x(minx)             # set lower left of box
        self.to_draw.set_y(miny)
        self.to_draw.set_width(maxx-minx)     # set width and height of box
        self.to_draw.set_height(maxy-miny)
        self.update()
    
    def plot_gaussians(self, params):
        plt.subplot(211)
        plt.cla()
        # Re-plot the original profile
        plt.plot(self.phases, self.profile, c='black', lw=3, alpha=0.3)
        plt.xlabel('Pulse Phase')
        plt.ylabel('Pulse Amplitude')
        DC = params[0]
        # Plot the individual gaussians
        for ii in range(self.numgaussians):
            phase, FWHM, amp = params[1+ii*3:4+ii*3]
            plt.plot(self.phases, DC + amp*gaussian_profile(self.proflen, phase, FWHM))

    def onselect(self):
        event1 = self.eventpress
        event2 = self.eventrelease
        # Left mouse button = add a gaussian
        if event1.button == event2.button == 1:
            x1, y1 = event1.xdata, event1.ydata
            x2, y2 = event2.xdata, event2.ydata
            phase = 0.5*(x1+x2)
            FWHM = Num.fabs(x2-x1)
            amp = Num.fabs(1.05*(y2-self.init_params[0])*(x2-x1))
            self.init_params += [phase, FWHM, amp]
            self.numgaussians += 1
            self.plot_gaussians(self.init_params)
            plt.draw()
        # Middle mouse button = fit the gaussians
        elif event1.button == event2.button == 2:
            fit_params, fit_errs, chi_sq, dof = \
                        fit_gaussians(self.profile, self.init_params,
                                      Num.zeros(self.proflen)+self.errs,
                                      self.profnm)
            # Save the fit parameters so the caller can retrieve them if needed
            self.fit_params = fit_params
            self.fit_errs = fit_errs
            # scaled uncertainties
            #scaled_fit_errs = fit_errs * Num.sqrt(chi_sq / dof)

            # Plot the best-fit profile
            self.plot_gaussians(fit_params)
            fitprof = gen_gaussians(fit_params, self.proflen)
            self.fitprof = fitprof
            plt.plot(self.phases, fitprof, c='black', lw=1)
            plt.draw()
            
            # Plot the residuals
            plt.subplot(212)
            plt.cla()
            residuals = self.profile - fitprof
            plt.errorbar(self.phases, residuals, self.errs,fmt='.')
            plt.grid(True)
            plt.xlabel('Pulse Phase')
            plt.ylabel('Data-Fit Residuals')
            plt.draw()
        # Right mouse button = remove last gaussian
        elif event1.button == event2.button == 3:
            if self.numgaussians:
                self.init_params = self.init_params[:-3]
                self.numgaussians -= 1
                self.plot_gaussians(self.init_params)
                plt.draw()
                plt.subplot(212)
                plt.cla()
                plt.xlabel('Pulse Phase')
                plt.ylabel('Data-Fit Residuals')
                plt.draw()

def gen_gaussians(params, N):
    """
    gen_gaussians(params, N):
        Return a model of a DC-component + M gaussians
            params is a sequence of 1+M*3 values
                the first value is the DC component.  Each remaining
                group of three represents the gaussians phase (0-1),
                FWHM (0-1), and amplitude (>0.0).
            N is the number of points in the model.
    """
    numgaussians = (len(params)-1) // 3
    model = Num.zeros(N, dtype='d') + params[0]
    for ii in range(numgaussians):
        phase, FWHM, amp = params[1+ii*3:4+ii*3]
        model += amp * gaussian_profile(N, phase, FWHM)
    return model

def fit_function(params, fjac=None, data=None, errs=None):
    return [0, (data - gen_gaussians(params, len(data))) / errs]

def fit_gaussians(data, initial_params, errs, profnm):
    numparams = len(initial_params)
    numgaussians = (len(initial_params)-1) // 3
    # Generate the parameter structure
    parinfo = []
    params0 = []
    for ii in range(numparams):
        params0.append(initial_params[ii])
        parinfo.append({'value':initial_params[ii], 'fixed':0,
                        'limited':[0,0], 'limits':[0.,0.]})
    other_args = {'data':data, 'errs':errs}
    # Now fit it
    mpfit_out = mpfit.mpfit(fit_function, params0, functkw=other_args,
                            parinfo=parinfo, quiet=1)
    fit_params = mpfit_out.params
    fit_errs = mpfit_out.perror
    # degrees of freedom
    dof = len(data) - len(fit_params)
    # chi-squared for the model fit
    chi_sq = mpfit_out.fnorm
    print("------------------------------------------------------------------")
    print("Multi-Gaussian Fit by pygaussfit.py of '%s'"%profnm)
    print("------------------------------------------------------------------")
    print("mpfit status:", mpfit_out.status)
    print("gaussians:", numgaussians)
    print("DOF:", dof)
    print("chi-sq: %.2f" % chi_sq)
    print("reduced chi-sq: %.2f" % (chi_sq/dof))
    residuals = data - gen_gaussians(fit_params, len(data))
    print("residuals  mean: %.3g" % Num.mean(residuals))
    print("residuals stdev: %.3g" % Num.std(residuals))
    print("--------------------------------------")
    print(" const = %.5f +/- %.5f" % (fit_params[0], fit_errs[0]))
    for ii in range(numgaussians):
        print(" phas%d = %.5f +/- %.5f" % (ii+1, fit_params[1+ii*3], fit_errs[1+ii*3]))
        print(" fwhm%d = %.5f +/- %.5f" % (ii+1, fit_params[2+ii*3], fit_errs[2+ii*3]))
        print(" ampl%d = %.5f +/- %.5f" % (ii+1, fit_params[3+ii*3], fit_errs[3+ii*3]))
    print("--------------------------------------")
    return fit_params, fit_errs, chi_sq, dof

if __name__ == '__main__':

    if len(sys.argv)==1:
        from numpy.random import normal

        print("""usage:  python pygaussfit.py input_file [prof_stdev]

Left mouse draws a region roughly boxing where you'll place a gaussian.
    Draw several to fit multiple gaussians.
Middle mouse performs the fit.
Right mouse removes the last gaussian from the fit.

The input_file should simply be an ASCII file with columns for pulse phase
and amplitude.  Comments with "#" are allowed.  .bestprof files work.

Paste the full resulting STDOUT to a '.gaussians' file for use
in get_TOAs.py or sum_profiles.py with the '-g' parameter as a template.""")


        N = 128
        DC = 600.0
        noise_stdev = 8.0
        
        params = [DC]
        params += [0.25, 0.1, 30.0]
        params += [0.3, 0.2, 15.0]
        params += [0.8, 0.05, 20.0]
        
        prof = normal(0.0, noise_stdev, N) + gen_gaussians(params, N)
        filenm = "test"
    else:
        if sys.argv[1].endswith(".pfd"):
            print("Input is PFD")
            # Input is pfd file
            pfdfn = sys.argv[1]
            # Check for bestprof
            if not os.path.exists(pfdfn+".bestprof"):
                print("Creating bestprof file")
                # Create bestprof file with show_pfd
                devnull = open(os.devnull, 'w')
                subprocess.call(['show_pfd', '-noxwin', pfdfn], 
                                stdout=devnull)
                devnull.close()
            filenm = pfdfn+".bestprof"
        else:
            filenm = sys.argv[1]
        prof = read_profile(filenm, normalize=0)
        if len(sys.argv)>=3:
            noise_stdev = float(sys.argv[2])
        else:
            try:
                bprof = bestprof(sys.argv[1])
                noise_stdev = bprof.prof_std
            except:
                noise_stdev = 1.0
    fig = plt.figure()
    dataplot = fig.add_subplot(211)
    interactor = GaussianSelector(dataplot, prof, noise_stdev, filenm)
    plt.show()
    
