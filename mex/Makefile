# MEX
#
MATLABROOT?=/usr/local/MATLAB_R2018a
ARCH?=glnxa64
EXT?=mexa64
MEXCFLAGS?=-I"$(MATLABROOT)/extern/include"
MEXLFLAGS?= -L"$(MATLABROOT)/bin/$(ARCH)" -lmex -lmx
FFTW_LIBS ?= -lfftw3

MEXSOURCES = $(wildcard *.c)
MEXTARGETS = $(patsubst %.c,%.$(EXT),$(MEXSOURCES))
OCTAVEMEXTARGETS = $(patsubst %.c,%.mex,$(MEXSOURCES))

export TARGETUP = ../libltfat/build/libphaseretd.a ../libltfat/build/libltfatd.a

CFLAGS =-std=c99 -Wall -Wextra -fPIC -shared

ifeq ($(COMPTARGET), debug)
	CFLAGS+=-O0 -g
else
	CFLAGS+=-O2 -DNDEBUG
endif

# MATLAB MEX files, only works if there are no spaces in the path to Matlab
matlab: CFLAGS+=-DLTFAT_BUILD_STATIC -DMATLAB_MEXCMD_RELEASE=R2017b
matlab: CFLAGS+=-DMATLAB_MEX_FILE -DLTFAT_LARGEARRAYS
matlab: CFLAGS+=-I../libltfat/modules/libltfat/include -I../libltfat/modules/libphaseret/include
# matlab: MBIN = $(shell matlab -nojvm -nodisplay -nosplash -r "disp([matlabroot,filesep,'bin',filesep,computer('arch'),filesep]);exit;" | grep /bin/)
# matlab: MWFFTWLIB = $(shell ls $(MBIN)*fftw3.* | head -n 1)
matlab: $(MEXTARGETS)

# This appends the correct suffix anyway
# and also only supresses warnings from mex (about old compiler) but not from
# the compiler itself.
%.$(EXT) : %.c
	$(CC) $(CFLAGS) $(MEXCFLAGS) $< -o $@ $(TARGETUP) $(MEXLFLAGS) $(FFTW_LIBS) -lm

# OCTAVE MEX files
# := $(shell ...) runs shell immediatelly even if we do not use the octave
# target at all
octave: MKOCTFILE ?= mkoctfile
octave: FFTW_LIBS = $(shell $(MKOCTFILE) -p FFTW_LIBS)
octave: export CFLAGS = $(shell $(MKOCTFILE) -p CFLAGS) $(CFLAGS)
octave: export LFLAGS = $(shell $(MKOCTFILE) -p LFLAGS) $(TARGETUP) $(FFTW_LIBS) -lm
octave: $(OCTAVEMEXTARGETS)

%.mex: %.c
	$(MKOCTFILE) --mex -o $@ $< 

.PHONY: clean matlab octave

clean:
	@rm -rf *.mex*
	@rm -rf *.o
