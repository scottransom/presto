# This start up file assumes the shell is 'tcsh'.
# Note: to 'uncomment' a line, remove the initial '#' character with an editor.

#______________________________________________________________________________
# This sets the default file creation mask.  '022' means 'not group writable
# and not world writable, but all other permissions exist'.
umask 022

# This prevents 'core dumps' from crashed programs from being created.
limit	coredumpsize 0

if ($?prompt) then
	set	prompt = "%m %h> "
	set	notify = 
	set	history = 100
	set	savehist = 20

# Here are some suggestions of useful aliases.  You may comment
# them out, add others of your own, etc., as you see fit.
	alias	h	history
	alias	m	more
	alias	lo	exit
	alias	vt100	"set term = vt100"
	alias df 'df -k'
	
	set noclobber
	set	lineedit


setenv MAIL /usr/spool/mail/$USER
setenv PRINTER lp342

#Edit the first 'setenv' line below to indicate the workstation you normally use
# (i.e. the one running the X server) and then uncomment this section
# if you want to be prompted for your $DISPLAY variable upon remote login.

#	setenv	HOMEHOST	myhost
#        if (($TERM == xterm || $TERM == xterms) && ($?DISPLAY == 0)) then
#            setenv REMOTEX 1
#            echo "Enter value for DISPLAY (host:#)(<CR> for "$HOMEHOST":0.0):"
#            set a=($<)
#            if ($a == "")  then
#                    setenv  DISPLAY $HOMEHOST":0.0"
#            else
#                    setenv  DISPLAY $a
#            endif
#       endif
#


endif

#______ LAYERED SOFTWARE _______________________________________________________

# Uncomment the lines corresponding to the package(s) you wish to run, making
# sure not to uncomment the description line, of course.

# Initializations

setenv  MANPATH /usr/man:/usr/local/man

set     path=( . \
                $HOME/bin \
                /usr/bin \
                /usr/local/bin \
        )
 


# TeX
setenv TEXINPUTS .:/usr/local/lib/tex/inputs
setenv TEXINPUTS ${TEXINPUTS}:/usr/local/lib/tex/inputs/latex
setenv TEXINPUTS ${TEXINPUTS}:/usr/local/lib/tex/inputs/latex209
setenv TEXINPUTS ${TEXINPUTS}:/usr/local/lib/tex/inputs/aasmacros
setenv TEXFORMATS /usr/local/lib/tex/formats
setenv TEXFONTS .:/usr/local/lib/tex/fonts

# xdvi fonts
setenv XDVIFONTS $TEXFONTS/pk300\:$TEXFONTS/pk118
#______________________________________________________________________________
# Mongo setup
# not available yet on the ALphas
# use the following setenv to send hardcopy files to printer $PRINTER
#setenv MONGOFILES /usr/local/lib/mongo/mongofiles.lw
# use the following setenv to save files in postscript format in the
# current directory.  the file names will be something like mgoXXXXX.ps
# where XXXXX stands for the process id.
#setenv MONGOFILES /usr/local/lib/mongo/mongofiles.save_ps
#______________________________________________________________________________
#______________________________________________________________________________
# PGPLOT (Caltech VLBI Graphics Subroutine Library)
setenv PGPLOT_FONT /usr/local/pgplot/grfont.dat
#______________________________________________________________________________

# GAG Software
set path=($path /home0/gag/exe/bin)
#
# IDL 
#
source /usr/local/rsi/idl_4/bin/idl_setup
setenv ASTRO_DATA $IDL_DIR/user_contrib/astron/data
setenv IDL_DIR /usr/local/rsi/idl_4                            # latest IDL version 
setenv IDL_PATH \+$IDL_DIR/lib:$IDL_DIR/user_contrib/astron/pro# astro lib
setenv IDL_DEVICE X  
alias idl $IDL_DIR/bin/idl                                     # generic IDL version
#
