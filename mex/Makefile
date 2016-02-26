# MEX
MEXSOURCES = $(wildcard *.c)
MEXTARGETS = $(patsubst %.c,%.mexa64,$(MEXSOURCES))
OCTAVEMEXTARGETS = $(patsubst %.c,%.mex,$(MEXSOURCES))

TARGETUP = ../$(TARGET)

MYCFLAGS = -std=c11 -pedantic -Wall -Wextra -I../include -I../include/phaseret -fPIC -Wl,--no-undefined
MYLFLAGS = -l:$(TARGETUP)

ifeq ($(CONFIG), debug)
	CFLAGS += -O0 -g
else
	CFLAGS += -O2 -DNDEBUG
endif

# MATLAB MEX files, only works if there are no spaces in the path to Matlab
matlab: MEX ?= mex
# matlab: MBIN = $(shell matlab -nojvm -nodisplay -nosplash -r "disp([matlabroot,filesep,'bin',filesep,computer('arch'),filesep]);exit;" | grep /bin/)
# matlab: MWFFTWLIB = $(shell ls $(MBIN)*fftw3.* | head -n 1)
matlab: $(TARGETUP) $(MEXTARGETS)

# This appends the correct suffix anyway
# and also only supresses warnings from mex (about old compiler) but not from
# the compiler itself.
%.mexa64 : %.c
	$(MEX) CFLAGS='$(MYCFLAGS)' $<  $(MYLFLAGS) -lfftw3 > /dev/null

# OCTAVE MEX files
# := $(shell ...) runs shell immediatelly even if we do not use the octave
# target at all
octave: MKOCTFILE ?= mkoctfile
octave: FFTW_LIBS = $(shell $(MKOCTFILE) -p FFTW_LIBS)
octave: export CFLAGS = $(shell $(MKOCTFILE) -p CFLAGS) $(MYCFLAGS)
octave: export LFLAGS = $(shell $(MKOCTFILE) -p LFLAGS) $(FFTW_LIBS) $(MYLFLAGS)
octave: $(TARGETUP) $(OCTAVEMEXTARGETS)

%.mex: %.c
	$(MKOCTFILE) --mex -o $@ $<

.PHONY: clean matlab octave

clean:
	@rm -rf *.mex*
