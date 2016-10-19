# This is the main makefile for libphaseret
#
# Builds three static and three shared libraries (default prefix is build/):
#
# 	libphaseret.a(.so)     Contains double, single and common code
# 	libphaseretd.a(.so)    Contains double and common code
# 	libphaseretf.a(.so)    Contains single and common code
#
# make CROSS=x86_64-w64-mingw32.static-
# or
# make CROSS=x86_64-w64-mingw32.static- NOBLASLAPACK=1
#
# Single, flattened header file build/phaseret.h can be obtained by
# make build/phaseret.h
# 
#

include ostools.mk

ifdef CROSS
	CC = $(CROSS)gcc
	AR = $(CROSS)ar
	RANLIB =ranlib
	buildprefix ?= build/$(CROSS)
	objprefix ?= obj/$(CROSS)
	MINGW=1
else
	CC ?= gcc
	AR ?= ar
	RANLIB?=ranlib
	buildprefix ?= build
	objprefix ?= obj
endif

# VERSION := $(shell cat phaseret_version)
PACKAGE = phaseret-$(VERSION)

PREFIX ?= /usr/local
LIBDIR = $(PREFIX)/lib
INCDIR = $(PREFIX)/include

# Base CFLAGS
CFLAGS+=-Wall -Wextra -pedantic -std=c99 -Iinclude -I../libltfat/include $(OPTCFLAGS)
CXXFLAGS+=-Wall -Wextra -pedantic -fno-exceptions -fno-rtti -std=c++11 -Iinclude -I../libltfat/include $(OPTCFLAGS)

# The following adds parameters to CFLAGS
COMPTARGET?=release
include comptarget.mk

# LIB
SOURCES = $(addprefix src/, dgtrealwrapper.c gla.c legla.c pghi.c rtisila.c \
		                    rtpghi.c spsi.c utils.c gsrtisila.c gsrtisilapghi.c )
SOURCES_TYPECONSTANT = $(wildcard src/*_typeconstant.c )
DOBJECTS = $(addprefix $(objprefix)/double/,$(notdir $(SOURCES:.c=.o)))
SOBJECTS = $(addprefix $(objprefix)/single/,$(notdir $(SOURCES:.c=.o)))
COMMONDOBJECTS = $(addprefix $(objprefix)/common/d,$(notdir $(SOURCES_TYPECONSTANT:.c=.o)))
COMMONSOBJECTS = $(addprefix $(objprefix)/common/s,$(notdir $(SOURCES_TYPECONSTANT:.c=.o)))

DSTATIC = libphaseretd.a
SSTATIC = libphaseretf.a
DSSTATIC = libphaseret.a

ifdef MINGW
	DSHARED = $(patsubst %.a,%.dll,$(DSTATIC))
	SSHARED = $(patsubst %.a,%.dll,$(SSTATIC))
	DSSHARED = $(patsubst %.a,%.dll,$(DSSTATIC))
	EXTRALFLAGS = -Wl,--out-implib,$@.a -static-libgcc
	CFLAGS += -DPHASERET_BUILD_SHARED
	CXXFLAGS += -DPHASERET_BUILD_SHARED
else
	CFLAGS += -fPIC
	CXXFLAGS += -fPIC
	DSHARED = $(patsubst %.a,%.so,$(DSTATIC))
	SSHARED = $(patsubst %.a,%.so,$(SSTATIC))
	DSSHARED = $(patsubst %.a,%.so,$(DSSTATIC))	
endif

ifdef USECPP
ifeq ($(USECPP),1)
	CC = $(CXX)
	CFLAGS = $(CXXFLAGS)
endif
endif

ifdef NOBLASLAPACK
	CFLAGS += -DNOBLASLAPACK
endif

# Define targets
DTARGET=$(buildprefix)/$(DSTATIC)
STARGET=$(buildprefix)/$(SSTATIC)
DSTARGET=$(buildprefix)/$(DSSTATIC)
SO_DTARGET=$(buildprefix)/$(DSHARED)
SO_STARGET=$(buildprefix)/$(SSHARED)
SO_DSTARGET=$(buildprefix)/$(DSSHARED)

#FFTWLIB ?= -lfftw3 -lfftw3f
LTFATFLIB ?= -lltfatf
LTFATDLIB ?= -lltfatd
LTFATLIB ?= -lltfat

#LIBS=$(FFTWLIB) -lm
LIBS=-lm

DDEP = $(buildprefix) $(objprefix)/double $(objprefix)/common
SDEP = $(buildprefix) $(objprefix)/single $(objprefix)/common

lib: static shared

$(DSTARGET): $(DDEP) $(SDEP) $(COMMONDOBJECTS) $(DOBJECTS) $(SOBJECTS)
	$(AR) rv $@ $(COMMONDOBJECTS) $(DOBJECTS) $(SOBJECTS)
	$(RANLIB) $@ 

$(DTARGET): $(DDEP) $(COMMONDOBJECTS) $(DOBJECTS)
	$(AR) rv $@ $(COMMONDOBJECTS) $(DOBJECTS)
	$(RANLIB) $@

$(STARGET): $(SDEP) $(COMMONSOBJECTS) $(SOBJECTS)
	$(AR) rv $@ $(COMMONSOBJECTS) $(SOBJECTS)
	$(RANLIB) $@

$(SO_DSTARGET): $(DDEP) $(SDEP) $(COMMONDOBJECTS) $(DOBJECTS) $(SOBJECTS)
	$(CC) -shared -Wl,--no-undefined -o $@ $(COMMONDOBJECTS) $(DOBJECTS) $(SOBJECTS) $(EXTRALFLAGS) $(LTFATLIB) $(LIBS)

$(SO_DTARGET): $(DDEP) $(COMMONDOBJECTS) $(DOBJECTS) 
	$(CC) -shared -Wl,--no-undefined -o $@ $(COMMONDOBJECTS) $(DOBJECTS) $(EXTRALFLAGS) $(LTFATDLIB) $(LIBS)

$(SO_STARGET): $(SDEP) $(COMMONSOBJECTS) $(SOBJECTS) 
	$(CC) -shared -Wl,--no-undefined -o $@ $(COMMONSOBJECTS) $(SOBJECTS) $(EXTRALFLAGS) $(LTFATFLIB) $(LIBS)

$(objprefix)/common/d%.o: src/%.c
	$(CC) $(CFLAGS) $(OPTCFLAGS) -DLTFAT_DOUBLE -c $< -o $@ 

$(objprefix)/double/%.o: src/%.c
	$(CC) $(CFLAGS) $(OPTCFLAGS) -DLTFAT_DOUBLE  -c $< -o $@

$(objprefix)/common/s%.o: src/%.c
	$(CC) $(CFLAGS) $(OPTCFLAGS) -DLTFAT_SINGLE -c $< -o $@

$(objprefix)/single/%.o: src/%.c
	$(CC) $(CFLAGS) $(OPTCFLAGS) -DLTFAT_SINGLE  -c $< -o $@

$(buildprefix):
	@$(MKDIR) $(buildprefix)

$(objprefix)/common:
	@$(MKDIR) $(objprefix)$(PS)common

$(objprefix)/double:
	@$(MKDIR) $(objprefix)$(PS)double

$(objprefix)/single:
	@$(MKDIR) $(objprefix)$(PS)single

matlab: $(TARGET)
	$(MAKE) -C mex matlab

octave: $(TARGET)
	$(MAKE) -C mex octave

cleanlib:
	$(RMDIR) $(buildprefix)
	$(RMDIR) $(objprefix)

clean: cleanlib

static: $(DSTARGET) $(DTARGET) $(STARGET) 

shared: $(SO_DSTARGET) $(SO_DTARGET) $(SO_STARGET) 

all: lib matlab octave

.PHONY: all doc doxy mat2doc mat2docmat clean cleanlib octave matlab cleandoc cleandoxy cleanmat2doc static shared munit
doc: doxy mat2doc

cleanmat2doc:
	@rm -rf mat2doc_publish

cleandoxy:
	@rm -rf html latex

cleandoc: cleanmat2doc cleandoxy

doxy:
	doxygen doxygen/doxyconfig

mat2doc: mat2docmat
	mat2doc . html

mat2docmat:
	mat2doc . mat

munit:
	$(MAKE) clean
	$(MAKE) $(SO_DSTARGET)
	$(MAKE) $(buildprefix)/phaseret.h USECPP=0

$(buildprefix)/phaseret.h: $(buildprefix)
	$(CC) -E -P -DNOSYSTEMHEADERS -nostdinc include/phaseret.h -Iinclude -I../libltfat/include -o $(buildprefix)/phaseret.h
	sed -i '1 i #ifndef _PHASERET_H' $(buildprefix)/phaseret.h
	sed -i '1 a #define _PHASERET_H' $(buildprefix)/phaseret.h
	sed -i '2 a #include <ltfat.h>' $(buildprefix)/phaseret.h
	sed -i '$$ a #endif' $(buildprefix)/phaseret.h

install: lib
	install -d $(DESTDIR)$(LIBDIR)
	install $(DTARGET) $(STARGET) $(DSTARGET) $(SO_DTARGET) $(SO_STARGET) $(SO_DSTARGET) $(DESTDIR)$(LIBDIR)
	mkdir -p $(DESTDIR)$(INCDIR)
	cp -r include/* $(DESTDIR)$(INCDIR)

uninstall:
	rm -f $(LIBDIR)/$(DSTATIC) $(LIBDIR)/$(SSTATIC) $(LIBDIR)/$(DSSTATIC)
	rm -f $(LIBDIR)/$(DSHARED) $(LIBDIR)/$(SSHARED) $(LIBDIR)/$(DSSHARED)
	rm -f $(DESTDIR)$(INCDIR)/phaseret.h
	rm -rf $(DESTDIR)$(INCDIR)/phaseret

