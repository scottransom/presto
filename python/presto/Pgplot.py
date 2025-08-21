"""
Routine for easy to use 1-D and 2-D plotting using 'PGPLOT'
      and the Python 'PPGPLOT' package

Written by Scott M. Ransom (ransom@cfa.harvard.edu)
         last revision: Aug 2025

'PGPLOT' was writtten by Tim Pearson <tjp@astro.caltech.edu>,
and can be found at http://astro.caltech.edu/~tjp/pgplot/

'PPGPLOT' was written by Nick Patavalis <npat@ariadne.di.uoa.gr>,
and can be found at http://ariadne.di.uoa.gr/ppgplot/
_or_ an updated version is available in the same directory
where this file was found:  ftp://cfa-ftp.harvard.edu/pub/ransom
"""

import sys
import numpy as np
from presto import ppgplot as ppg


# Check if string in Py2 and Py3 compatible way
def isstr(var):
    return isinstance(var, str if sys.version_info[0] >= 3 else basestring)


# True if we have an /XWIN or /XSERVE device open yet
ppgplot_dev_open_ = False

# True if we have already scaled and/or prepped the current page
ppgplot_dev_prep_ = False

# Default plotting device
ppgplot_device_ = "/XWIN"

# Default font to use
ppgplot_font_ = 1

# Default font size to use
ppgplot_font_size_ = 1.0

# Default line style to use
ppgplot_linestyle_ = 1

# Default line width to use
ppgplot_linewidth_ = 2

# Default symbol to plot
ppgplot_symbol_ = None

# Default label color for contour lines
ppgplot_labels_ = None

# Default label interval for labeling contour lines
ppgplot_labelint_ = 20

# Default minimum label interval for labeling contour lines
ppgplot_labelmin_ = 20

# Default line colors to use
ppgplot_color_ = "white"

# Default color palette for IMAG routines
ppgplot_palette_ = "rainbow"

# The set of colors for PGPLOT
ppgplot_colors_ = {
    "black": 0,
    "Black": 0,
    "BLACK": 0,
    "white": 1,
    "White": 1,
    "WHITE": 1,
    "red": 2,
    "Red": 2,
    "RED": 2,
    "green": 3,
    "Green": 3,
    "GREEN": 3,
    "blue": 4,
    "Blue": 4,
    "BLUE": 4,
    "cyan": 5,
    "Cyan": 5,
    "CYAN": 5,
    "magenta": 6,
    "Magenta": 6,
    "MAGENTA": 6,
    "yellow": 7,
    "Yellow": 7,
    "YELLOW": 7,
    "orange": 8,
    "Orange": 8,
    "ORANGE": 8,
    "green2": 9,
    "Green2": 9,
    "GREEN2": 9,
    "green3": 10,
    "Green3": 10,
    "GREEN3": 10,
    "blue2": 11,
    "Blue2": 14,
    "BLUE2": 11,
    "purple": 12,
    "Purple": 12,
    "PURPLE": 12,
    "pink": 13,
    "Pink": 13,
    "PINK": 13,
    "darkgray": 14,
    "DarkGray": 14,
    "DARKGRAY": 14,
    "dark gray": 14,
    "Dark Gray": 14,
    "DARK GRAY": 14,
    "lightgray": 15,
    "LightGray": 15,
    "LIGHTGRAY": 15,
    "light gray": 15,
    "Light Gray": 15,
    "LIGHT GRAY": 15,
}


# Data should be a Numpy Array
def scalerange(data):
    """
    scalerange(data):
        Adjust the range to be plotted so that it fits nicely on the page.
        Return a list with adjusted minimum and maximum values from 'data'.
    """
    min = np.minimum.reduce(data)
    max = np.maximum.reduce(data)
    extra = 0.1 * (max - min)
    return [min - extra, max + extra]


# Reset global variables to defaults
def resetdefaults():
    """
    resetdefaults():
        Reset global plotting variables to default values.
    """
    global \
        ppgplot_font_, \
        ppgplot_linestyle_, \
        ppgplot_linewidth_, \
        ppgplot_color_, \
        ppgplot_font_size_
    ppg.pgscf(ppgplot_font_)
    ppg.pgsls(ppgplot_linestyle_)
    ppg.pgslw(ppgplot_linewidth_)
    ppg.pgsci(ppgplot_colors_[ppgplot_color_])
    ppg.pgsch(ppgplot_font_size_)


# Go to a subsequent plotting page
def nextplotpage(reset=0):
    """
    nextplotpage():
        Advance the plotting device to a new page.
        The optional entry is:
            reset: reset defaults or not (default = 0 (no))
    """
    global ppgplot_dev_open_, ppgplot_dev_prep_
    if ppgplot_dev_open_:
        ppg.pgpage()
        ppgplot_dev_prep_ = False
    else:
        print("Can't go to the next plot page unless a plotting device is open.")
    if reset:
        resetdefaults()


# Reset the color indices to the default values
def reset_colors():
    lo_col_ind, hi_col_ind = ppg.pgqcol()
    ppg.pgscir(lo_col_ind, hi_col_ind)
    ppg.pgscr(0, 0.00, 0.00, 0.00)  # Black (background)
    ppg.pgscr(1, 1.00, 1.00, 1.00)  # White (default)
    ppg.pgscr(2, 1.00, 0.00, 0.00)  # Red
    ppg.pgscr(3, 0.00, 1.00, 0.00)  # Green
    ppg.pgscr(4, 0.00, 0.00, 1.00)  # Blue
    ppg.pgscr(5, 0.00, 1.00, 1.00)  # Cyan (Green + Blue)
    ppg.pgscr(6, 1.00, 0.00, 1.00)  # Magenta (Red + Blue)
    ppg.pgscr(7, 1.00, 1.00, 0.00)  # Yellow  (Red + Green)
    ppg.pgscr(8, 1.00, 0.50, 0.00)  # Red + Yellow (Orange)
    ppg.pgscr(9, 0.50, 1.00, 0.00)  # Green + Yellow
    ppg.pgscr(10, 0.00, 1.00, 0.50)  # Green + Cyan
    ppg.pgscr(11, 0.00, 0.50, 1.00)  # Blue + Cyan
    ppg.pgscr(12, 0.50, 0.00, 1.00)  # Blue + Magenta
    ppg.pgscr(13, 1.00, 0.00, 0.50)  # Red + Magenta
    ppg.pgscr(14, 0.33, 0.33, 0.33)  # Dark Gray
    ppg.pgscr(15, 0.66, 0.66, 0.66)  # Light Gray
    for ci in range(16, hi_col_ind + 1):
        ppg.pgscr(ci, 0.00, 0.00, 0.00)  # Black (background)


# Open a plotting device
def prepplot(
    rangex,
    rangey,
    title=None,
    labx=None,
    laby=None,
    rangex2=None,
    rangey2=None,
    labx2=None,
    laby2=None,
    logx=False,
    logy=False,
    logx2=False,
    logy2=False,
    font=ppgplot_font_,
    fontsize=ppgplot_font_size_,
    id=False,
    aspect=1,
    vstd=True,
    ticks="in",
    panels=[1, 1],
    device=ppgplot_device_,
):
    """
    prepplot(rangex, rangey, ...)
        Open a PGPLOT device for plotting.
            'rangex' and 'rangey' are sequence objects giving min and
                max values for each axis.
        The optional entries are:
            title:    graph title                 (default = None)
            labx:     label for the x-axis        (default = None)
            laby:     label for the y-axis        (default = None)
            rangex2:  ranges for 2nd x-axis       (default = None)
            rangey2:  ranges for 2nd y-axis       (default = None)
            labx2:    label for the 2nd x-axis    (default = None)
            laby2:    label for the 2nd y-axis    (default = None)
            logx:     make the 1st x-axis log     (default = False)
            logy:     make the 1st y-axis log     (default = False)
            logx2:    make the 2nd x-axis log     (default = False)
            logy2:    make the 2nd y-axis log     (default = False)
            font:     PGPLOT font to use          (default = 1 (normal))
            fontsize: PGPLOT font size to use     (default = 1.0 (normal))
            id:       Show ID line on plot        (default = False)
            aspect:   Aspect ratio                (default = 1 (square))
            vstd:     Standard viewport?          (default = True)
            ticks:    Ticks point in or out       (default = 'in')
            panels:   Number of subpanels [r,c]   (default = [1,1])
            device:   PGPLOT device to use        (default = '/XWIN')
        Note:  Many default values are defined in global variables
            with names like ppgplot_font_ or ppgplot_device_.
    """
    global ppgplot_dev_open_, ppgplot_dev_prep_
    # Check if we will use second X or Y axes
    # Note:  if using a 2nd X axis, the range should correspond
    #   to the minimum and maximum values of the 1st X axis.  If
    #   using a 2nd Y axis, the range should correspond to the
    #   scalerange() values of the 1st Y axis.
    otherxaxis = True
    if rangex2 is None:
        rangex2 = rangex
        otherxaxis = False
    otheryaxis = True
    if rangey2 is None:
        rangey2 = rangey
        otheryaxis = False
    # Open the plot device
    if not ppgplot_dev_open_:
        ppg.pgopen(device)
        # Let the routines know that we already have a device open
        ppgplot_dev_open_ = True
        # Set the aspect ratio
        ppg.pgpap(0.0, aspect)
        if panels != [1, 1]:
            # Set the number of panels
            ppg.pgsubp(panels[0], panels[1])
            ppg.pgpage()
    # Choose the font
    ppg.pgscf(font)
    # Choose the font size
    ppg.pgsch(fontsize)
    # Choose the font size
    ppg.pgslw(ppgplot_linewidth_)
    # Plot the 2nd axis if needed first
    if otherxaxis or otheryaxis:
        # ppg.pgvstd()
        ppg.pgswin(rangex2[0], rangex2[1], rangey2[0], rangey2[1])
        # Decide how the axes will be drawn
        env = "CMST" if ticks == "in" else "CMSTI"
        lxenv = "L" if logx2 else ""
        lyenv = "L" if logy2 else ""
        if otherxaxis and otheryaxis:
            ppg.pgbox(env + lxenv, 0.0, 0, env + lyenv, 0.0, 0)
        elif otheryaxis:
            ppg.pgbox("", 0.0, 0, env + lyenv, 0.0, 0)
        else:
            ppg.pgbox(env + lxenv, 0.0, 0, "", 0.0, 0)
    # Now setup the primary axis
    # ppg.pgvstd()
    ppg.pgswin(rangex[0], rangex[1], rangey[0], rangey[1])
    # Decide how the axes will be drawn
    env = "ST" if ticks == "in" else "STI"
    lxenv = "L" if logx else ""
    lyenv = "L" if logy else ""
    if otherxaxis and otheryaxis:
        ppg.pgbox("BN" + env + lxenv, 0.0, 0, "BN" + env + lyenv, 0.0, 0)
    elif otheryaxis:
        ppg.pgbox("BCN" + env + lxenv, 0.0, 0, "BN" + env + lyenv, 0.0, 0)
    elif otherxaxis:
        ppg.pgbox("BN" + env + lxenv, 0.0, 0, "BCN" + env + lyenv, 0.0, 0)
    else:
        ppg.pgbox("BCN" + env + lxenv, 0.0, 0, "BCN" + env + lyenv, 0.0, 0)
    # Add labels
    if title is not None:
        ppg.pgmtxt("T", 3.2, 0.5, 0.5, title)
    ppg.pgmtxt("B", 3.0, 0.5, 0.5, labx)
    ppg.pgmtxt("L", 2.6, 0.5, 0.5, laby)
    if otherxaxis:
        ppg.pgmtxt("T", 2.0, 0.5, 0.5, labx2)
    if otheryaxis:
        ppg.pgmtxt("R", 3.0, 0.5, 0.5, laby2)
    # Add ID line if required
    if id:
        ppg.pgiden()
    # Let the routines know that we have already prepped the device
    ppgplot_dev_prep_ = True


# Close plotting device
def closeplot():
    """
    closeplot():
        Close the currently open plotting device
    """
    global ppgplot_dev_open_, ppgplot_dev_prep_
    ppg.pgend()
    ppgplot_dev_open_ = False
    ppgplot_dev_prep_ = False


# Plot simple XY line plots with optional arguments and keywords
def plotxy(
    y,
    x=None,
    title=None,
    rangex=None,
    rangey=None,
    labx="",
    laby="",
    rangex2=None,
    rangey2=None,
    labx2="",
    laby2="",
    symbol=ppgplot_symbol_,
    line=ppgplot_linestyle_,
    width=ppgplot_linewidth_,
    color=ppgplot_color_,
    font=ppgplot_font_,
    logx=False,
    logy=False,
    logx2=False,
    logy2=False,
    errx=None,
    erry=None,
    id=False,
    noscale=False,
    aspect=0.7727,
    fontsize=ppgplot_font_size_,
    ticks="in",
    panels=[1, 1],
    device=ppgplot_device_,
    setup=True,
):
    """
    plotxy(y, ...)
        An interface to make various XY style plots using PGPLOT.
            'y' is the 1D sequence object to plot.
        The optional entries are:
            x:        x values                  (default = 0, 1, ...)
            title:    graph title               (default = None)
            rangex:   ranges for the x-axis     (default = automatic)
            rangey:   ranges for the y-axis     (default = automatic)
            labx:     label for the x-axis      (default = None)
            laby:     label for the y-axis      (default = None)
            rangex2:  ranges for 2nd x-axis     (default = None)
            rangey2:  ranges for 2nd y-axis     (default = None)
            labx2:    label for the 2nd x-axis  (default = None)
            laby2:    label for the 2nd y-axis  (default = None)
            logx:     make the 1st x-axis log   (default = False)
            logy:     make the 1st y-axis log   (default = False)
            logx2:    make the 2nd x-axis log   (default = False)
            logy2:    make the 2nd y-axis log   (default = False)
            errx:     symmetric x errors        (default = None)
            erry:     symmetric y errors        (default = None)
            symbol:   symbol for points         (default = None)
            line:     line style                (default = 1 (solid))
            width:    line width                (default = 1 (thin))
            color:    line and/or symbol color  (default = 'white')
            font:     PGPLOT font to use        (default = 1 (normal))
            fontsize: PGPLOT font size to use   (default = 1.0 (normal))
            id:       show ID line on plot      (default = False)
            noscale:  turn off auto scaling     (default = False)
            aspect:   aspect ratio              (default = 0.7727 (rect))
            ticks:    Ticks point in or out     (default = 'in')
            panels:   Number of subpanels [r,c] (default = [1,1])
            device:   PGPLOT device to use      (default = '/XWIN')
            setup:    Auto-setup the plot       (default = True)
        Note:  Many default values are defined in global variables
            with names like ppgplot_font_ or ppgplot_device_.
    """
    # Make sure the input data is an array
    y = np.asarray(y)
    # Announce the global variables we will be using
    global ppgplot_dev_open_, ppgplot_dev_prep_, ppgplot_colors_
    # Define the X axis limits if needed
    x = np.arange(len(y), dtype="f") if x is None else np.asarray(x)
    # Determine the scaling to use for the first axis
    if rangex is None:
        rangex = [x.min(), x.max()]
    if rangey is None:
        if noscale:
            rangey = [y.min(), y.max()]
        else:
            rangey = scalerange(y)
    # Prep the plotting device...
    if not ppgplot_dev_prep_ and setup:
        prepplot(
            rangex,
            rangey,
            title,
            labx,
            laby,
            rangex2,
            rangey2,
            labx2,
            laby2,
            logx,
            logy,
            logx2,
            logy2,
            font,
            fontsize,
            id,
            aspect,
            ticks,
            panels,
            device=device,
        )
    # Choose the line color
    ppg.pgsci(ppgplot_colors_[color] if isstr(color) else color)
    # Plot symbols (and errors) if requested
    if symbol is not None:
        ppg.pgpt(x, y, symbol)
    # Error bars
    if errx is not None:
        if not logx:
            errx = np.asarray(errx)
            ppg.pgerrx(x + errx, x - errx, y, 1.0)
        else:
            errx = 10.0 ** np.asarray(errx)
            ppg.pgerrx(np.log10(10.0**x + errx), np.log10(10.0**x - errx), y, 1.0)
    if erry is not None:
        if not logy:
            erry = np.asarray(erry)
            ppg.pgerry(x, y + erry, y - erry, 1.0)
        else:
            erry = 10.0 ** np.asarray(erry)
            ppg.pgerry(x, np.log10(10.0**y + erry), np.log10(10.0**y - erry), 1.0)
    # Plot connecting lines if requested
    if line is not None:
        # Choose the line style
        ppg.pgsls(line)
        # Choose the line width
        ppg.pgslw(width)
        ppg.pgline(x, y)


# Make an X-Y plot of binned data (i.e. useful for histograms)
def plotbinned(
    y,
    x=None,
    title=None,
    labx="Bins",
    laby="Counts",
    rangex=None,
    rangey=None,
    labx2="",
    laby2="",
    rangex2=None,
    rangey2=None,
    line=ppgplot_linestyle_,
    width=ppgplot_linewidth_,
    color=ppgplot_color_,
    font=ppgplot_font_,
    logx=False,
    logy=False,
    logx2=False,
    logy2=False,
    erry=None,
    id=False,
    noscale=False,
    aspect=0.7727,
    fontsize=ppgplot_font_size_,
    ticks="out",
    panels=[1, 1],
    device=ppgplot_device_,
    setup=True,
):
    """
    plotbinned(y, ...):
        Plot x-y data that is binned.  This routine differs from
            plotxy() in that instead of each point being connected
            by diagonal lines, each point is actually a flat-line
            with the width of a bin.
            'y' is the numerical sequence of binned data to plot.
        The optional entries are:
            x:        x-centers of each bin.    (default = auto)
            title:    graph title               (default = None)
            labx:     label for the x-axis      (default = 'Bins')
            laby:     label for the y-axis      (default = 'Counts')
            rangex:   ranges for the x-axis     (default = automatic)
            rangey:   ranges for the y-axis     (default = automatic)
            labx2:    label for the 2nd x-axis  (default = None)
            laby2:    label for the 2nd y-axis  (default = None)
            rangex2:  ranges for 2nd x-axis     (default = None)
            rangey2:  ranges for 2nd y-axis     (default = None)
            logx:     make the 1st x-axis log   (default = False)
            logy:     make the 1st y-axis log   (default = False)
            logx2:    make the 2nd x-axis log   (default = False)
            logy2:    make the 2nd y-axis log   (default = False)
            erry:     symmetric y errors        (default = None)
            line:     line style                (default = 1 (solid))
            width:    line width                (default = 1 (thin))
            color:    line and/or symbol color  (default = 'white')
            font:     PGPLOT font to use        (default = 1 (normal))
            fontsize: PGPLOT font size to use   (default = 1.0 (normal))
            id:       show ID line on plot      (default = False)
            aspect:   aspect ratio              (default = 0.7727 (rect))
            ticks:    Ticks point in or out     (default = 'in')
            panels:   Number of subpanels [r,c] (default = [1,1])
            device:   PGPLOT device to use      (default = '/XWIN')
            setup:    Auto-setup the plot       (default = True)
        Note:  Many default values are defined in global variables
            with names like ppgplot_font_ or ppgplot_device_.
    """
    # Make sure our entry sequences are Num arrays
    np.asarray(y)
    np.asarray(x)
    if x is None:
        x = np.arange(len(y)) + 0.5
    dx = x[1] - x[0]
    # Correct for the fact that 'x' are the bin centers
    x = x - 0.5 * dx
    # Make the repeat array
    r = np.zeros(len(x), dtype=np.int32) + 2
    ny = np.repeat(y, r)
    r[0] = 1
    nx = np.repeat(x, r)
    # Add the right side of the right-most bin
    nx = np.concatenate((nx, np.zeros(1) + nx[-1] + dx))
    plotxy(
        ny,
        nx,
        title,
        labx=labx,
        laby=laby,
        line=line,
        labx2=labx2,
        laby2=laby2,
        rangex2=rangex2,
        rangey2=rangey2,
        logx=logx,
        logy=logy,
        logx2=logx2,
        logy2=logy2,
        noscale=noscale,
        width=width,
        color=color,
        font=font,
        fontsize=fontsize,
        id=id,
        aspect=aspect,
        rangex=rangex,
        rangey=rangey,
        ticks=ticks,
        panels=panels,
        device=device,
        setup=setup,
    )
    if erry is not None:
        ppg.pgerry(np.arange(len(y)) + 0.5, y + erry, y - erry, 1.0)


# Show a 2D color intensity plot with optional arguments and keywords
def plot2d(
    z,
    x=None,
    y=None,
    title=None,
    rangex=None,
    rangey=None,
    rangez=None,
    labx="",
    laby="",
    rangex2=None,
    rangey2=None,
    labx2="",
    laby2="",
    image=ppgplot_palette_,
    contours=None,
    logx=False,
    logy=False,
    logx2=False,
    logy2=False,
    line=ppgplot_linestyle_,
    width=ppgplot_linewidth_,
    color=ppgplot_color_,
    labels=ppgplot_labels_,
    labelint=ppgplot_labelint_,
    labelmin=ppgplot_labelmin_,
    font=ppgplot_font_,
    id=False,
    noscale=False,
    aspect=1,
    fontsize=ppgplot_font_size_,
    ticks="out",
    panels=[1, 1],
    device=ppgplot_device_,
    setup=True,
):
    """
    plot2d(z, ...)
        An interface to make various 2D plots using PGPLOT.
            'z' is the 2D Numpy array to be plotted.
        The optional entries are:
            x:         x values                    (default = 0, 1, ...)
            y:         y values                    (default = 0, 1, ...)
            title:     graph title                 (default = None)
            rangex:    range for the x-axis        (default = automatic)
            rangey:    range for the y-axis        (default = automatic)
            rangez:    range for the z-axis        (default = automatic)
            labx:      label for the x-axis        (default = None)
            laby:      label for the y-axis        (default = None)
            rangex2:   range for 2nd x-axis        (default = None)
            rangey2:   range for 2nd y-axis        (default = None)
            labx2:     label for the 2nd x-axis    (default = None)
            laby2:     label for the 2nd y-axis    (default = None)
            logx:      make the 1st x-axis log     (default = False)
            logy:      make the 1st y-axis log     (default = False)
            logx2:     make the 2nd x-axis log     (default = False)
            logy2:     make the 2nd y-axis log     (default = False)
            image:     color palette for image     (default = 'rainbow')
            contours:  list of contour values      (default = None)
            line:      contour line style          (default = 1 (solid))
            width:     contour line width          (default = 1 (thin))
            color:     contour line color          (default = 'white')
            labels:    color of contour labels     (default = None)
            labelint:  contour label spacing       (default = 20)
            labelmin:  min contour label spacing   (default = 20)
            font:      PGPLOT font to use          (default = 1 (normal))
            fontsize:  PGPLOT font size to use     (default = 1.0 (normal))
            id:        show ID line on plot        (default = False)
            noscale:   turn off auto scaling       (default = False)
            aspect:    Aspect ratio                (default = 1 (square))
            ticks:     Ticks point in or out       (default = 'out')
            panels:    Number of subpanels [r,c]   (default = [1,1])
            device:    PGPLOT device to use        (default = '/XWIN')
            setup:    Auto-setup the plot          (default = True)
        Note:  Many default values are defined in global variables
            with names like ppgplot_font_ or ppgplot_device_.
    """
    # Make sure the input data is a 2D array
    z = np.asarray(z)
    if not len(z.shape) == 2:
        print("Input data array must be 2 dimensional.")
        return
    # Announce the global variables we will be using
    global ppgplot_dev_open_, ppgplot_dev_prep_, pgpalette
    # Define the X and Y axis limits if needed
    x = np.arange(z.shape[1], dtype="f") if x is None else np.asarray(x)
    y = np.arange(z.shape[0], dtype="f") if y is None else np.asarray(y)
    # Determine the scaling to use for the axes
    if rangex is None:
        dx = x[-1] - x[-2]
        rangex = [x[0], x[-1] + dx]
    if rangey is None:
        dy = y[-1] - y[-2]
        rangey = [y[0], y[-1] + dy]
    if rangez is None:
        rangez = [np.minimum.reduce(np.ravel(z)), np.maximum.reduce(np.ravel(z))]
    # Prep the plotting device...
    if not ppgplot_dev_prep_ and setup:
        prepplot(
            rangex,
            rangey,
            title,
            labx,
            laby,
            rangex2,
            rangey2,
            labx2,
            laby2,
            logx,
            logy,
            logx2,
            logy2,
            font,
            fontsize,
            id,
            aspect,
            ticks,
            panels,
            device=device,
        )
    if image is not None:
        # Set the color indices and the color table
        lo_col_ind, hi_col_ind = ppg.pgqcol()
        lo_col_ind = lo_col_ind + 2
        ppg.pgscir(lo_col_ind, hi_col_ind)
        pgpalette.setpalette(image)
        ppg.pgctab(pgpalette.l, pgpalette.r, pgpalette.g, pgpalette.b)
        # Construct the image
        ppg.pgimag_s(z, 0.0, 0.0, rangex[0], rangey[0], rangex[1], rangey[1])
        reset_colors()
    if contours is not None:
        contours = np.asarray(contours)
        # Choose the line style
        ppg.pgsls(line)
        # Choose the line width
        ppg.pgslw(width)
        # Choose the line color for the contourlines
        ppg.pgsci(ppgplot_colors_[color] if isstr(color) else color)
        # Construct the contours
        ppg.pgcont_s(
            z, len(contours), contours, rangex[0], rangey[0], rangex[1], rangey[1]
        )
        # Label the contours if requested
        if labels is not None:
            # Choose the line color for the contourlines
            ppg.pgsci(ppgplot_colors_[labels] if isstr(labels) else labels)
            for i in range(len(contours)):
                ppg.pgconl_s(z, contours[i], str(contours[i]), labelint, labelmin)


#########################################################################


class Palette(object):
    # Set the color palette
    def setpalette(self, palette):
        """
        setpalette(self, palette):
            Set the color palette for imag-style routines
        """
        if palette == "rainbow":
            self.l = np.array(
                [0.0, 0.015, 0.225, 0.4, 0.59, 0.6, 0.775, 0.955, 0.965, 1.0]
            )
            self.r = np.array([1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.947, 1.0, 1.0])
            self.g = np.array([0.0, 0.0, 1.0, 1.0, 1.0, 0.946, 0.0, 0.8, 0.844, 1.0])
            self.b = np.array([0.0, 0.0, 0.0, 0.0, 0.95, 1.0, 1.0, 1.0, 1.0, 1.0])
        elif palette == "antirainbow":
            self.l = np.array(
                [0.0, 0.035, 0.045, 0.225, 0.4, 0.41, 0.6, 0.775, 0.985, 1.0]
            )
            self.r = np.array([1.0, 1.0, 0.947, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0])
            self.g = np.array([1.0, 0.844, 0.8, 0.0, 0.946, 1.0, 1.0, 1.0, 0.0, 0.0])
            self.b = np.array([1.0, 1.0, 1.0, 1.0, 1.0, 0.95, 0.0, 0.0, 0.0, 0.0])
        elif palette == "astro":
            self.l = np.array([0.0, 0.167, 0.333, 0.5, 0.667, 0.833, 1.0])
            self.r = np.array([0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0])
            self.g = np.array([0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 1.0])
            self.b = np.array([0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0])
        elif palette == "hue":
            self.l = np.array([0.0, 0.167, 0.333, 0.5, 0.667, 0.833, 1.0])
            self.r = np.array([1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0])
            self.g = np.array([0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0])
            self.b = np.array([0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0])
        elif palette == "heat":
            self.l = np.array([0.0, 0.48, 0.7, 0.75, 1.0])
            self.r = np.array([0.0, 1.0, 1.0, 1.0, 1.0])
            self.g = np.array([0.0, 0.0, 0.423, 0.519, 1.0])
            self.b = np.array([0.0, 0.0, 0.0, 0.0, 1.0])
        elif palette == "gamma":
            self.l = np.array([0.0, 0.33, 0.66, 1.0])
            self.r = np.array([0.3, 1.0, 0.0, 0.0])
            self.g = np.array([0.0, 0.3, 1.0, 0.0])
            self.b = np.array([0.0, 0.0, 0.3, 1.0])
        elif palette == "antigray" or palette == "antigrey":
            self.l = np.array([0.0, 1.0])
            self.r = np.array([1.0, 0.0])
            self.g = np.array([1.0, 0.0])
            self.b = np.array([1.0, 0.0])
        elif palette == "apjgray" or palette == "apjgrey":
            self.l = np.array([0.0, 1.0])
            self.r = np.array([1.0, 0.25])
            self.g = np.array([1.0, 0.25])
            self.b = np.array([1.0, 0.25])
        else:
            self.l = np.array([0.0, 1.0])
            self.r = np.array([0.0, 1.0])
            self.g = np.array([0.0, 1.0])
            self.b = np.array([0.0, 1.0])


pgpalette = Palette()

#
# Demo code
#
if __name__ == "__main__":
    from math import *
    from numpy import *

    def distance(width):
        """
        distance(width):
            Return a 'width' x 'width' Numpy array with each
                point set to the geometric distance from the array's center.
        """
        x = np.arange(-width / 2.0 + 0.5, width / 2.0 + 0.5, 1.0) ** 2
        x = np.resize(x, (width, width))
        return np.sqrt(x + np.transpose(x))

    # Do a couple 1-D plots

    x = arange(0.0, 10.0, 0.05)
    xcm = x * 2.54
    rx2 = [min(xcm), max(xcm)]
    ry2 = [-0.25, 0.25]
    y = cos(x)
    f = exp(-0.1 * x)

    # Show the simplest calling sequence
    plotxy(y)
    closeplot()

    # Show something a little more useful
    plotxy(
        y,
        x,
        rangex2=rx2,
        rangey2=ry2,
        labx="inches",
        laby="foobar activity",
        labx2="cm",
        laby2="aged foobar activity",
        id=1,
    )

    # Over-plot the following
    plotxy(y * f, x, color="red", line=2, width=6)
    closeplot()

    # Show a couple 2-D examples

    a = exp(-0.02 * distance(200))
    ca = a * cos(0.04 * distance(200))

    # Show the simplest calling sequence
    plot2d(a)
    closeplot()

    # Show 3 related plots which are a little more useful
    plot2d(
        ca,
        x,
        x,
        title="Contours",
        labx="x",
        laby="y",
        image=None,
        contours=[0.0, 0.4, 0.8],
        labels="yellow",
        color="red",
        labelint=40,
        labelmin=20,
    )
    closeplot()

    # Show the same thing but with an image
    plot2d(ca, x, x, title="Image", labx="x", laby="y", image="heat")
    closeplot()

    # Show the same thing but with an image and contours
    plot2d(
        ca,
        x,
        x,
        title="Image+Contours",
        labx="x",
        laby="y",
        image="heat",
        contours=[0.0, 0.4, 0.8],
    )
    closeplot()
