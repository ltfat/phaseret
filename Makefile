CC ?= gcc
AR ?= ar

VERSION := $(shell cat phaseret_version)
PACKAGE = phaseret-$(VERSION)

PREFIX ?= /usr/local
LIBDIR = $(PREFIX)/lib
INCDIR = $(PREFIX)/include

# LIB
LIBSOURCES = $(wildcard src/*.c)
LIBOBJECTS = $(patsubst %.c,%.o,$(LIBSOURCES))
$(info $(LIBSOBJECTS))

TARGET=build/libphaseret.a
SO_TARGET=$(patsubst %.a,%.so,$(TARGET))

# MEX
MEXSOURCES = $(wildcard mex/*.c)
MEXTARGETS = $(patsubst %.c,%.mexa64,$(MEXSOURCES))
OCTAVEMEXTARGETS = $(patsubst %.c,%.mex,$(MEXSOURCES))

CFLAGS=-std=c11 -pedantic -Wall -Wextra -DNDEBUG -I./include $(OPTFLAGS)
LIBS=-lm -lfftw3

all: lib matlab

lib: $(TARGET) $(SO_TARGET)

dev: CFLAGS=-std=c11 -g -O0 -Wall -Wall -Wextra -I./src $(OPTFLAGS)
dev: all octave

$(TARGET): CFLAGS += -fPIC
$(TARGET): build $(LIBOBJECTS)
	ar rvu $@ $(LIBOBJECTS)
	ranlib $@

$(SO_TARGET): $(TARGET) $(LIBOBJECTS)
	$(CC) -shared -Wl,--no-undefined -o $@ $(LIBOBJECTS) $(LIBS)

%.o: %.c
	$(CC) -c $(CFLAGS) $< -o $@ 

# MATLAB MEX files
matlab: MEX ?= mex
matlab: CFLAGS = -std=c11 -pedantic -Wall -Wextra -DNDEBUG -I./src -fPIC
matlab: LDFLAGS = -l:$(TARGET) $(LIBS)
matlab: $(TARGET) $(MEXTARGETS)

# This appends the correct suffix anyway
# and also only supresses warnings from mex (about old compiler) but not from 
# the compiler itself.
%.mexa64 : %.c
	$(MEX) CFLAGS='$(CFLAGS)' $< $(LDFLAGS) -outdir mex -largeArrayDims > /dev/null

# OCTAVE MEX files
# := $(shell ...) runs shell immediatelly even if we do not use the octave 
# target at all
octave: MKOCTFILE = mkoctfile
octave: FFTW_LIBS = $(shell $(MKOCTFILE) -p FFTW_LIBS)
octave: export CFLAGS = $(shell $(MKOCTFILE) -p CFLAGS) -std=c11 -DNDEBUG -I./src -Wall -Wextra -fPIC
octave: export LFLAGS = $(shell $(MKOCTFILE) -p LFLAGS) $(FFTW_LIBS)
octave: $(TARGET) $(OCTAVEMEXTARGETS) 

%.mex: %.c
	$(MKOCTFILE) --mex -o $@ $<

build:
	@mkdir -p build
	@mkdir -p bin

clean:
	@rm -rf build $(LIBOBJECTS)
	@rm -rf mex/*.mex*

.PHONY: doc doxy mat2doc mat2docmat clean octave matlab cleandoc cleandoxy cleanmat2doc
doc: doxy mat2doc

cleanmat2doc:
	@rm -rf mat2doc_publish

cleandoxy:
	@rm -rf html

cleandoc: cleanmat2doc cleandoxy

doxy:
	doxygen doxygen/doxyconfig

mat2doc: mat2docmat
	mat2doc . html

mat2docmat:
	mat2doc . mat

install: all
	install -d $(DESTDIR)/$(LIBDIR)
	install $(TARGET) $(DESTDIR)/$(LIBDIR)
	mkdir -p $(DESTDIR)/$(INCDIR)
	cp -r include/phaseret.h $(DESTDIR)/$(INCDIR)

uninstall:
	rm -f $(DESTDIR)/$(LIBDIR)
	rm -f $(DESTDIR)/$(INCDIR)/phaseret.h

