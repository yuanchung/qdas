# Make file for qas
# By Yuan-Chung Cheng <yccheng@berkeley.edu>
#

CC = gcc
CFLAGS = -O3 -finline-functions -Wall -Winline -funroll-loops -fexpensive-optimizations
#CFLAGS = -g -DBATH_DEBUG

#LIBDIR = -L/home/platin/x-prog/ycclib -L/usr/local/lib -L/sw/lib
LIBDIR = -L/usr/local/lib64 -L.
EXTRALIBS = -lqdas -lpthread -lgsl -lgslcblas -lm
#GSL_LAPACK = -lg2c -lblas -llapack
INCLUDEDIR = -I/usr/local/include -I/usr/X11R6/include
#INCLUDEDIR = -I/usr/local/include -I/usr/X11R6/include -I/sw/include -I/home/ntucc/platin/opt/include

OBJECTS = aux.o spectrum.o bath_js03.o bath_htgau.o bath_mt99ohm.o bath_js03art.o bath_ohmart.o bath.o bathtype_mt99art.o params.o dmop.o parser.o qng2d.o gauss-hermite.o bath_mboart.o bath_mt99art.o spectrum2d.o linpolar.o samplelist.o

all: $(OBJECTS) libqdas

%.o : %.c %.h
	$(CC) -c $(CFLAGS) $(INCLUDEDIR) $< -o $@

libqdas: $(OBJECTS)
	rm -f libqdas.a
	ar -crusv libqdas.a $(OBJECTS)

test_js03art: bath_js03art.c bath_js03art.h bathtype_mt99art.o 
	$(CC) -static $(CFLAGS) $(INCLUDEDIR) -DMAIN -o bath_js03art bath_js03art.c $(LIBDIR) $(EXTRALIBS)

test_ohmart: bath_ohmart.c bath_ohmart.h bathtype_mt99art.o 
	$(CC) $(CFLAGS) $(INCLUDEDIR) -DMAIN -o bath_ohmart bath_ohmart.c $(LIBDIR) $(EXTRALIBS)

test_mboart: bath_mboart.c bath_mboart.h bathtype_mt99art.o 
	$(CC) $(CFLAGS) $(INCLUDEDIR) -DMAIN -o bath_mboart bath_mboart.c $(LIBDIR) $(EXTRALIBS)

test_mt99art: bath_mt99art.c bath_mt99art.h bathtype_mt99art.o 
	$(CC) $(CFLAGS) $(INCLUDEDIR) -DMAIN -o bath_mt99art bath_mt99art.c $(LIBDIR) $(EXTRALIBS)

test_qng2d: qng2d.c qng2d.h
	$(CC) $(CFLAGS) $(INCLUDEDIR) -DMAIN -DDEBUG_QNG2D -o qng2d qng2d.c $(LIBDIR) $(EXTRALIBS)

test_gauss-hermite: gauss-hermite.c gauss-hermite.h aux.o
	$(CC) $(CFLAGS) $(INCLUDEDIR) -DMAIN -DDEBUG_GAUSS_HERMITE -o gauss-hermite gauss-hermite.c $(LIBDIR) $(EXTRALIBS) aux.o

test_linpolar: linpolar.c linpolar.h aux.o
	$(CC) $(CFLAGS) $(INCLUDEDIR) -DMAIN -o linpolar linpolar.c $(LIBDIR) -lgsl -lgslcblas -lm aux.o

clean:
	rm -f *.a *.o *.so *~ *.mo core a.out gmon.out bath_js03art bath_ohmart bath_mboart bath_mt99art qng2d linpolar

#
# $Log$
# Revision 1.19  2007/07/27 22:50:57  platin
#   - minor correction; test two model spectral density in bath_mt99art.c.
#
# Revision 1.18  2007/06/20 18:14:03  platin
#
#   - added time stamp function.
#
# Revision 1.17  2007/06/13 17:22:03  platin
#
#   - basic support for 2D spectrum data type.
#
# Revision 1.16  2007/03/15 08:20:50  platin
#
#   - add direct support for MT99 artificial bath.
#
# Revision 1.15  2007/03/15 02:36:26  platin
#
#   - add multi-mode brownian oscillator model.
#
# Revision 1.14  2007/03/09 06:01:52  platin
#
#   - import new Gauss-Hermite Quadrature module that can be used
#     for static disorder average.
#
# Revision 1.13  2007/01/03 07:30:00  platin
#
#   - add a new package in which pulse duration integral is carried out
#     to include pulse overlap effect using impulsive response
#     function formalism.
#
# Revision 1.12  2006/11/08 06:23:50  platin
#
#   - minor changes and adding a artificial ohmic bath module that
#     uses more efficient artificial fit to Ohmic correlation
#     functions.
#
# Revision 1.11  2006/11/02 21:56:34  platin
#
#   - minor change
#
# Revision 1.10  2006/11/02 21:55:44  platin
#
#   - remove unnessary .h dependency.
#
# Revision 1.9  2006/11/02 19:29:31  platin
#
#   - implement the bathtype interface.
#
# Revision 1.8  2006/10/30 06:31:31  platin
#
#   - add parser library to avoid external dependence.
#
# Revision 1.7  2006/08/15 23:05:06  platin
#
#   - remove ad-hoc fix in tnl-* dynamics; now both programs correctly
#     handle renormalization term (Hren) and eliminate erroneous dynamics,
#     such as nonpositive dynamics and non-zero long time coherence terms...
#
# Revision 1.6  2006/08/03 02:15:27  platin
#
#   - fix the pulse period setup; this should handle the FWHM correctly.
#   - in the js03art module, save some time when gamma is effectively zero.
#
# Revision 1.5  2006/08/01 19:10:51  platin
#
#   - include artificial JS03 bath model that uses four terms for
#     M&T's parameterization of bath correlation functions.
#
# Revision 1.4  2006/07/10 20:25:41  platin
#
#   - a possible sign problem in dm3pes.c,
#   - fist import of the tnl-dm2pes module, which calculates 3P PE signal
#     using non-Markovian TNL method from Meier and Tannor.
#
# Revision 1.3  2006/05/26 23:13:54  platin
#
#   - bug fix and minor changes in makefiles.
#
# Revision 1.2  2006/05/26 19:19:30  platin
#
#   - add dynamics module.
#   - revise the params.c, use a more reasonable scheme to handle
#     parameters for different modules.
#
# Revision 1.1.1.1  2006/05/24 00:42:19  platin
#
#   - initial import of the qdas package.
#   - qdas stands for quantum dynamics and spectroscopy.
#   - basic code inherited from the "lineshape" package.
#
#
#
