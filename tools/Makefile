# Make file for qas
# By Yuan-Chung Cheng <yccheng@berkeley.edu>
#

CC = gcc
CFLAGS = -O2 -finline-functions -Wall -Winline -funroll-loops -fexpensive-optimizations
#CFLAGS = -g -Wno-long-double -DDEBUG

LIBDIR = -L/usr/local/lib -L/sw/lib -L../lib
EXTRALIBS = -lqdas -lpthread -lgsl -lgslcblas -lm 
#GSL_LAPACK = -lg2c -lblas -llapack
INCLUDEDIR = -I/usr/local/include -I/usr/X11R6/include -I/sw/include -I../lib

OBJECTS = 
BINARIES = cspline_max normalize 1cpeps 2c3ppes-hj 2c3ppes gau_max 2c3ppes-pi 2c3ppes-htgau 2dspec 2danalysis

INSTALL = /usr/bin/install -c
BINDIR = /usr/local/chem/qdas
INSTALL_BINARIES = cspline_max normalize gau_max 2dspec 2danalysis

all: $(OBJECTS) $(BINARIES)

%.o : %.c %.h
	$(CC) -c $(CFLAGS) $(INCLUDEDIR) $< -o $@

cspline_max: cspline_max.c
	$(CC) $(CFLAGS) $(INCLUDEDIR) $(LIBDIR) -o cspline_max cspline_max.c $(EXTRALIBS)

gau_max: gau_max.c
	$(CC) $(CFLAGS) $(INCLUDEDIR) $(LIBDIR) -o gau_max gau_max.c $(EXTRALIBS)

normalize: normalize.c
	$(CC) $(CFLAGS) $(INCLUDEDIR) $(LIBDIR) -o normalize normalize.c $(EXTRALIBS)

1cpeps: 1cpeps.c $(OBJECTS)
	$(CC) $(CFLAGS) $(INCLUDEDIR) $(LIBDIR) -o 1cpeps 1cpeps.c $(EXTRALIBS)

2c3ppes: 2c3ppes.c $(OBJECTS)
	$(CC) $(CFLAGS) $(INCLUDEDIR) $(LIBDIR) -o 2c3ppes 2c3ppes.c $(EXTRALIBS)

2c3ppes-hj: 2c3ppes-hj.c $(OBJECTS)
	$(CC) $(CFLAGS) $(INCLUDEDIR) $(LIBDIR) -o 2c3ppes-hj 2c3ppes-hj.c $(EXTRALIBS)

2c3ppes-pi: 2c3ppes-pi.c $(OBJECTS)
	$(CC) $(CFLAGS) $(INCLUDEDIR) $(LIBDIR) -o 2c3ppes-pi 2c3ppes-pi.c $(EXTRALIBS)

2c3ppes-htgau: 2c3ppes-htgau.c $(OBJECTS)
	$(CC) $(CFLAGS) $(INCLUDEDIR) $(LIBDIR) -o 2c3ppes-htgau 2c3ppes-htgau.c $(EXTRALIBS)

2dspec: 2dspec.c $(OBJECTS)
	$(CC) $(CFLAGS) $(INCLUDEDIR) $(LIBDIR) -o 2dspec 2dspec.c $(EXTRALIBS) -lfftw3 -lm

2danalysis: 2danalysis.c $(OBJECTS)
	$(CC) $(CFLAGS) $(INCLUDEDIR) $(LIBDIR) -o 2danalysis 2danalysis.c $(EXTRALIBS) -lfftw3 -lm

install: all
	mkdir -p $(BINDIR)
	$(INSTALL) -m 755 $(INSTALL_BINARIES) $(BINDIR)

clean:
	rm -f $(OBJECTS) $(BINARIES)
	rm -f *.a *.o *.so *~ *.mo core a.out gmon.out 

#
# $Log$
# Revision 1.11  2007/05/10 01:47:21  platin
#   - flip w2 output in 2dspec; this makes it easier for analysis.
#   - add 2danalysis program that at this point only useful for
#     computing volume under a 2d curve.
#
# Revision 1.10  2007/02/02 20:01:56  platin
#
#   - bug fix.
#
# Revision 1.9  2007/02/02 08:07:49  platin
#
#   - basic implement for tools to generate 2D spectra.
#
# Revision 1.8  2007/01/12 07:19:59  platin
#
#   - major changes in 2c3ppes* files.
#   - add a 2c3ppes program that uses HTGAU bath.
#
# Revision 1.7  2007/01/03 07:30:00  platin
#
#   - add a new package in which pulse duration integral is carried out
#     to include pulse overlap effect using impulsive response
#     function formalism.
#
# Revision 1.6  2006/12/05 07:33:18  platin
# *** empty log message ***
#
# Revision 1.5  2006/11/21 22:56:47  platin
#
#   - 2c3ppes: generate two-color three-pulse photon-echo signal using
#     response function formalism in the implusive limit.
#
# Revision 1.4  2006/11/17 04:50:09  platin
#
#   - 2c3pes generates full 2-C three-pulse photon-echo signals using
#     implusive response function technique.
#
# Revision 1.3  2006/10/31 23:26:22  platin
#
#   - new aux. program used to compute 2C photon-echo signals in Hohjai's
#     scrambled time-ordering case. i.e. 750, 800, 750.
#
# Revision 1.2  2006/07/15 07:45:59  platin
#
#   - add a simple tool that can be used to compute 1C-PEPS data using
#     Cho's expression.
#
# Revision 1.1  2006/07/12 17:32:42  platin
#
#   - cspline_max: use cspline interpolation to estimate maximum points.
#   - normalize: normalize data set so that maximum of y is 1.0.
#
#
