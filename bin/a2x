#!/bin/sh
# This is a2x
#
# Copyright (C) 1994-1996 Christoph Beck, <beck@jack.rhein-main.de>
# All rights reserved.
#
# This software may be distributed under the terms of the GNU General
# Public License.
#
# See the manual page for details of installation and usage.
#
# Any bugfixes, improvements, as well as comments are welcome.
# Email to beck@jack.rhein-main.de
#
# Enjoy!
#
#------------------- BEGIN USER DEFINITIONS ---------------------------
#
# where did you copy a2x.ps and a2xshell.ps ?
A2X_PATH=$PRESTO/lib
#A2X_PATH=.
#
# Your printer (use 'gs -h' for a list of devices) or ps for PostScript
#Device=deskjet
Device=ps
#
# one of a4, letter, legal
Paper=letter
#
# Press RET to get next page for these devices:
InteractiveDevices="x11 linux"
#
#-------------------- END USER DEFINITIONS ----------------------------

files=""
Out="-"
A2Xoptions=""
ShowpageHook="-dNOPAUSE"
Jstr=""
#
while test -n "$1"
do
  case $1 in
    -[h?])
cat << EOHELP

----------------------------------------------+---------------------------
a2x: ascii --> something                      |  a2x version 1.0  12/23/96
(c) Christoph Beck, <beck@jack.rhein-main.de> |        All rights reserved
----------------------------------------------+---------------------------

  usage: a2x [-dDEV]
 	     [-a4|-letter|-legal]
             [-title[=name]] [-date] [-num] [-pagenum]
             [-p] [-l]
             [-cn]
             [-o|-e|-O|-E]
	     [-nm|-sn]
	     [-bs] [-tn]
	     [-man]
	     [-ffont]
	     [-ger|-us]
	     [-pagecount]
             [files]

If no files are given, a2x reads from stdin. A2x always writes to stdout,
except for gs-devices like x11 (may be set in a2x).

	... hit <return> to get more ...
EOHELP
  read
cat << EOHELP

Options:

  -dDEV     causes a2x (actually gs) to produce output for device DEV. 
            Default output device is $Device (change this in a2x, if you like;
            use gs -h for possible devices). Use ps for PostScript device.
  -a4,-letter,-legal choose paper sizes.
            The default paper size is $Paper (you may change this in a2x, too).
  -title[=name] prints name (which defaults to the file's name)
            on top (left) of each page.
  -date     causes the current date being printed on top (right) of each page.
  -p        print in portrait mode (default).
  -1up,-2up 1,2 column portrait,landscape modes (make life easy).
  -l        print in landscape mode.
  -man      print (preformatted) man page (try a2x -man -2up ... pretty!).
  -cn       n columns per page (default 1-col-portrait, 2-col-landscape).
  -tn       set tabwidth to n (default 8).
  -bs       turn on backspace handling.
  -sn       select fontsize (default 10pt-portrait, 7pt-landscape).
  -ffont    select font (default is Courier).
  -nm       set lines per column to m.

	... hit <return> to get more ...
EOHELP
  read
cat << EOHELP

  -e,-o     don't print odd,even pages.
  -E,-O     print only even,odd columns.
  -R        print columns per page in reverse order (from right to left).
  -D        Duplex (double sided) output, if your PS device supports this.
  -D -T     Tumble (double sided, pages tumbled) output, if your PS ...
  -num      print column numbers on bottom of each column (default off).
	    ** in versions up to 0.9 -num behaved like -pagenum, see below **
  -pagenum  print page numbers on bottom of each page (default off).
  -iso      ISO-Latin1 character encoding (default, requires PS level 2).
  -ger,-us  german (default), us font encoding.
  -pagecount just prints the number of pages to stdout (and nothing else).










	... hit <return> to get more ...
EOHELP
  read
cat << EOHELP

Hints:

-2up makes the distance between the left and right column equal to the
sum of the left and right margins of the whole page, -l -c2 does not.

Think of columns as logical pages when using -2up and -num; think of
them as just columns composing a page when using -cn, -p, -l and -pagenum.

You can do double sided prints with -o, -r, eg
$ nroff -man a2x.1 | a2x -man -o | lpr	# odd pages
$ nroff -man a2x.1 | a2x -man -r | lpr	# even pages
or, if your printer supports Duplex, just
$ nroff -man a2x.1 | a2x -man -D | lpr	# duplex
to get the same.

You can do a half sized (eg a5) double sided printing with -O, -E, -R:
$ a2x -a4 -2up -num -O    file | lpr	# prints logical pages [1|3] [5|7] ...
Refeed paper to your printer (same order, upside-down) and do
$ a2x -a4 -2up -num -E -R file | lpr	# prints logical pages [4|2] [8|6] ...
Cut between columns, merge and you're done.

Have a good time.

EOHELP
	exit;;
    -date) A2Xoptions="$A2Xoptions --date=`date +%D`";;
    -d[a-z0-9]*) Device=`expr $1 : '-d\(.*\)'`
	for i in $InteractiveDevices
	do
	  if test $Device = $i; then Interactive=TRUE; fi
	done;;
    -fax) Device=faxg3; A2Xoptions="-p -s12 $A2Xoptions"; Out=page%d.$$;;
    -a4|-letter|-legal) Paper=`expr $1 : '-\(.*\)'`;;
    -title=*) Jstr="$Jstr`expr $1 : '-title=\(.*\)'`" # to be used in header
        A2Xoptions="$A2Xoptions -title";;
    -*) A2Xoptions="$A2Xoptions $1";;
     *) if test -r $1; then files="$files $1"; fi;;
  esac
  shift
done


Papersize="$Paper"size
A2Xoptions="--paper=$Papersize $A2Xoptions"

#
# (1) PostScript
#
if test $Device = "ps"
then
  A2Xargs=""; for j in $A2Xoptions; do A2Xargs="$A2Xargs ($j)"; done
  if test -z "$files"  # read from stdin, write to stdout
  then
    if test -z "$Jstr"; then Job=""; else Job="(--job=$Jstr)"; fi
    cat $A2X_PATH/a2x.ps
    echo "a2xdict begin [ (--eof=1) $A2Xargs $Job () ] A2x"
    cat - | sed "s///g"
    echo 
    echo "end"
  else                 # read from file(s)
    cat $A2X_PATH/a2x.ps
    echo "a2xdict begin"
    for i in $files
    do
      if test -z "$Jstr"; then Job="(--job=$i)"; else Job="(--job=$Jstr)"; fi
      echo  "[ (--eof=1) $A2Xargs $Job () ] A2x"
      cat $i | sed "s///g"
      echo 
    done
    echo "end"
  fi   
  exit
fi

#
# (2) GhostScript
#
GSoptions="-sDEVICE=$Device -sPAPERSIZE=$Paper -sOutputFile=$Out -q"
if test -z $Interactive; then GSoptions="$GSoptions -sNOPAUSE"; fi

if test -z "$files" # read from stdin, write to stdout
then
  if test -z "$Jstr"; then Job=""; else Job="(--job=$Jstr)"; fi
  gs $GSoptions $A2X_PATH/a2x.ps \
     -- $A2X_PATH/a2xshell.ps $A2Xoptions $Job %stdin
else		    # read from file(s)
  for i in $files
  do
    if test -z "$Jstr"; then Job="--job=$i"; else Job="--job=$Jstr"; fi
    gs $GSoptions $A2X_PATH/a2x.ps -- $A2X_PATH/a2xshell.ps $A2Xoptions $Job $i
  done
fi

#------------------------------------------------------------------------------
