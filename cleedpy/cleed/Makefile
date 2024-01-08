#============================================================================
#  GH/26.03.2020
#  makefile for the nonsymmetry LEED program
#============================================================================
.SUFFIXES: .o .c .h
#=============================================================================
MYPRG = .
INCLUDEDIR = $(MYPRG)/include
LIBDIR = .
#

# GH 26/03/2020 added for openblas
INCLUDEDIRBLAS = /opt/OpenBLAS/include/
LIBDIRBLAS = /opt/OpenBLAS/lib/
FLAGSBLAS = -L$(LIBDIRBLAS) -Wl,-rpath,$(LIBDIRBLAS) -lopenblas

#==flags=====================================================================
CCOMP = gcc
OPT = -O3 -march=native
# Use this for homemade BLAS:
# OPT = -O2
CFLAGSSUB = -g -c $(OPT) -I$(INCLUDEDIR) -I$(INCLUDEDIRBLAS) -L$(LIBDIR)
CFLAGS = -g $(OPT) -I$(INCLUDEDIR) -I$(INCLUDEDIRBLAS)
LDFLAGS = -lm


# GH 24/03/2020 added for "homemade" CBLAS
# INCLUDEDIRBLAS = $(MYPRG)/CBLAS/include
# "homemade" CBLAS library
# CBLASLIB = $(MYPRG)/CBLAS/lib/cblas_LINUX.a

# system CBLAS library (not recommended)
# CBLASLIB = /lib64/libblas.so /lib64/libgslcblas.so

#============================================================================
#  header files
#
FILES.h = $(INCLUDEDIR)/stddef.h \
          $(INCLUDEDIR)/leed.h
#===========================================================================
#  subroutines  for matrix operations
#
# complex library
CPLOBJ =  cribasfun.o

# basic qm functions
QMOBJ  =  qmbessm.o \
          qmcgc.o   \
          qmhank.o  \
          qmylm.o

# matrix library
MATOBJ =  matabs.o     \
          matalloc.o   \
          matarralloc.o\
          matarrfree.o \
          matcblasaux_lp.o \
          matcgau.o    \
          matcheck.o   \
          matclu.o     \
          matcol.o     \
          matconj.o    \
          matcop.o     \
          matdiff.o    \
          matext.o     \
          matfree.o    \
          matins.o     \
          matinv_lp.o     \
          matmul_lp.o     \
          matnattovht.o\
          matread.o    \
          matrlu.o     \
          matrow.o     \
          matrm.o      \
          matscal.o    \
          matshow.o    \
          matshowabs.o \
          matshowpar.o \
          matsqmod.o   \
          mattrace.o   \
          mattrans.o   \
          matwrite.o

# input for LEED programs
INPOBJ =  linpbullayer.o \
          linpdebtemp.o  \
          linpovllayer.o \
          linpphasend.o  \
          linprdbulnd.o  \
          linprdpar.o    \
          linprdovlnd.o  \
          linpshowbop.o

          # lreadpar.o
          # lwritepar.o
          # linpmatlm.o

# output for LEED programs
OUTOBJ =  loutbmlist.o \
          louthead.o   \
          loutint.o

# beams:
BEAMOBJ = lbmgen.o    \
          lbmselect.o \
          lbmset.o

          # lbmrotmat.o \

# parameter control:
PCOBJ =   lpcmktlnd.o \
          lpctemtl.o \
          lpcupdate.o

          # lpcmktl.o

# layer doubling:
LDOBJ  = lld2n.o \
         lld2lay.o \
         lld2layrpm.o \
         lldpotstep.o \
         lldpotstep0.o

# multiple scattering:
MSOBJ  = lmsbravlnd.o \
         lmscomplnd.o \
         lmslsumii.o  \
         lmslsumij.o  \
         lmspartinv.o \
         lmstmatii.o  \
         lmstmatndii.o \
         lmstmatij.o  \
         lmsymat.o    \
         lmsymatset.o \
         lmsymmat.o   \
         lmsypy.o

         # lmsbravl.o
         # lmscompl.o
         # lmsltok.o

PRGOBJ = prgcheck.o

TMAOBJ = lpccumtl.o \
         lpcmkms.o


LEEDOBJ = $(CPLOBJ) \
        $(QMOBJ) \
        $(MATOBJ) \
        $(INPOBJ) \
        $(OUTOBJ) \
        $(BEAMOBJ) \
        $(PCOBJ) \
        $(LDOBJ) \
        $(MSOBJ) \
        $(PRGOBJ)

TMATOBJ = $(CPLOBJ) \
        $(QMOBJ) \
        $(MATOBJ) \
        $(INPOBJ) \
        $(PCOBJ) \
        $(PRGOBJ)\
        $(OUTOBJ) \
        $(BEAMOBJ) \
        $(LDOBJ) \
        $(MSOBJ) \
        $(TMAOBJ)

TMATOBJ_531 = $(CPLOBJ) \
        $(QMOBJ) \
        $(MATOBJ) \
        $(INPOBJ) \
        $(PCOBJ) \
        $(PRGOBJ)\
        $(OUTOBJ_531) \
        $(BEAMOBJ) \
        $(LDOBJ) \
        $(MSOBJ) \
        $(TMAOBJ)


#============================================================================
.c.o:   ;  $(CCOMP) $(CFLAGSSUB) $*.c
#============================================================================
test_nsym_lp:	test_nsym_lp.o $(TMATOBJ)
	$(CCOMP) -o $(MYPRG)/test_nsym_lp $(CFLAGS) $(TMATOBJ) test_nsym_lp.c $(LDFLAGS) $(FLAGSBLAS)
	rm -f test_nsym_lp.o
#============================================================================
