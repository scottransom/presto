from __future__ import print_function
from future import standard_library
standard_library.install_aliases()
from builtins import range
def catvar(col):
    ret = []
    global a, b, c, d, e
    for i in range(shape(a)[0]):
        ret.append(a[i][col])
    for i in range(shape(b)[0]):
        ret.append(b[i][col])
    for i in range(shape(c)[0]):
        ret.append(c[i][col])
    for i in range(shape(d)[0]):
        ret.append(d[i][col])
    for i in range(shape(e)[0]):
        ret.append(e[i][col])
    return ret

def readsaves(file='montebinopt_saves.txt'):
    f = open(file, 'r')
    result = []
    while 1:
        try:
            result.append(load(f))
        except ValueError:
            continue
        except EOFError:
            break
    f.close()
    return result

def saveresults(file="testresults.txt"):
    from pickle import *
    global psrp, orbp, orbx, orbe, orbw, orbt
    global widthp, widthx, widtht, widthe, widthw
    global mf, z
    vars = ('psrp', 'orbp', 'orbx', 'orbe', 'orbw', 'orbt',
            'widthp', 'widthx', 'widtht', 'widthe', 'widthw',
            'mf', 'z')
    f = open(file, 'w')
    for var in vars:
        print('Saving ', var, '...')
        exec('dump(%s, f)' % (var))
    f.close()
    print('Saved em.')

def quadratic(parameters, x):
    a = parameters[0]
    b = parameters[1]
    c = parameters[2]
    return (a * x + b) * x + c

def linear(parameters, x):
    m = parameters[0]
    b = parameters[1]
    return m * x + b

def genfits():
    from LeastSquares import leastSquaresFit
    global psrp, orbp, orbx, orbe, orbw, orbt
    global widthp, widthx, widtht, widthe, widthw
    global mf, z
    yvals = {'Orb p':widthp, 'Orb x':widthx, 'Orb t':widtht}
    xvals = {'Orb p':orbp/z, 'Orb x':1.0/z, 'Orb t':1.0/z}
    xtits = {'Orb p':'Orb p/z', 'Orb x':'1.0/z', 'Orb t':'1.0/z'}
    for fitvar in ['Orb p', 'Orb x', 'Orb t']:
        vals = []
        for i in range(len(xvals[fitvar])):
            vals.append((xvals[fitvar][i], yvals[fitvar][i]))
        fit = leastSquaresFit(linear, (1.0, 0.0), vals)
        print('%s width = %10.7f * %s + %10.7f  (Acc: %f)' % (fitvar,
                                                              fit[0][0],
                                                              xtits[fitvar],
                                                              fit[0][1],
                                                              fit[1]))
        plotxy(yvals[fitvar], xvals[fitvar],
               laby=fitvar+' Width (Fractional)',
               labx=xtits[fitvar], line=None, font=2,
               symbol=2, color='red', device='fits.ps/CPS')
        plotxy([fit[0][0]*min(xvals[fitvar])+fit[0][1],
                fit[0][0]*max(xvals[fitvar])+fit[0][1]],
               [min(xvals[fitvar]), max(xvals[fitvar])],
               line=1, symbol=None, color='blue')
        nextplotpage(1)
    closeplot()

def genlogfits():
    from LeastSquares import leastSquaresFit
    global psrp, orbp, orbx, orbe, orbw, orbt
    global widthp, widthx, widtht, widthe, widthw
    global mf, z
    yvals = {'Orb p':log(widthp), 'Orb x':log(widthx), 'Orb t':log(widtht)}
    xvals = {'Orb p':log(orbp/z), 'Orb x':log(1.0/z), 'Orb t':log(1.0/z)}
    xtits = {'Orb p':'log(Orb p/z)', 'Orb x':'log(1.0/z)',
             'Orb t':'log(1.0/z)'}
    for fitvar in ['Orb p', 'Orb x', 'Orb t']:
        vals = []
        for i in range(len(xvals[fitvar])):
            vals.append((xvals[fitvar][i], yvals[fitvar][i]))
        fit = leastSquaresFit(linear, (1.0, 0.0), vals)
        print('log(%s) width = %10.7f * %s + %10.7f  (Acc: %f)' % (fitvar,
                                                                   fit[0][0],
                                                                   xtits[fitvar],
                                                                   fit[0][1],
                                                                   fit[1]))
        plotxy(yvals[fitvar], xvals[fitvar],
               laby='log('+fitvar+') Width (Fractional)',
               labx=xtits[fitvar], line=None, font=2,
               symbol=2, color='red', device='logfits.ps/CPS')
        plotxy([fit[0][0]*min(xvals[fitvar])+fit[0][1],
                fit[0][0]*max(xvals[fitvar])+fit[0][1]],
               [min(xvals[fitvar]), max(xvals[fitvar])],
               line=1, symbol=None, color='blue')
        nextplotpage(1)
    closeplot()

if __name__ == '__main__':
    from math import *
    from Numeric import *
    from stats import *
    from Statistics import *
    from Pgplot import *
    from miscutils import *

    def help(funct):
        """
        help(funct):
        Print the documentation string of a function or method.
        """
        print(eval(funct + '.__doc__'))
    
    from pickle import *
    vars = ('psrp', 'orbp', 'orbx', 'orbe', 'orbw', 'orbt',
            'widthp', 'widthx', 'widtht', 'widthe', 'widthw',
            'mf', 'z')
    f = open("testresults.txt")
    for var in vars:
        print('Loading ', var, '...')
        exec(var + ' = asarray(load(f))')
    f.close()
    print('Got em.')
    
