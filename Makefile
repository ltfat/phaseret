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

include ostools.mk

# VERSION := $(shell cat phaseret_version)
PACKAGE = phaseret-$(VERSION)

PREFIX ?= /usr/local
LIBDIR = $(PREFIX)/lib
INCDIR = $(PREFIX)/include

# LIB
LIBSOURCES = $(wildcard src/*.c)
LIBOBJECTS = $(addprefix $(objprefix)/,$(notdir $(LIBSOURCES:.c=.o)))

# Base CFLAGS
CFLAGS+=-Wall -std=c99 -Iinclude -Ithirdparty $(OPTCFLAGS)
CXXFLAGS+=-Wall -std=c++11 -Iinclude -Ithirdparty $(OPTCFLAGS)

# The following adds parameters to CFLAGS
include comptarget.mk

STATIC = libphaseret.a
ifdef MINGW
	SHARED = libphaseret.dll
	EXTRALFLAGS = -Wl,--out-implib,$@.a -static-libgcc
else
	CFLAGS += -fPIC
	SHARED = libphaseret.so
endif

TARGET=$(buildprefix)/$(STATIC)
SO_TARGET=$(buildprefix)/$(SHARED)

FFTWLIB ?= -lfftw3

LIBS=-lltfat  $(FFTWLIB) -lm

lib: $(TARGET) $(SO_TARGET)

all: lib matlab octave

$(TARGET): $(objprefix) $(buildprefix) $(LIBOBJECTS)
	$(AR) rvu $@ $(LIBOBJECTS)
	$(RANLIB) $@

$(SO_TARGET): $(objprefix) $(buildprefix) $(LIBOBJECTS)
	$(CC) -shared -Wl,--no-undefined -o $@ $(LIBOBJECTS) $(EXTRALFLAGS) $(LIBS)

$(objprefix)/%.o: src/%.c $(objprefix)
	$(CC) -c $(CFLAGS) $< -o $@

$(buildprefix):
	$(MKDIR) $(buildprefix)

$(objprefix):
	$(MKDIR) $(objprefix)

matlab: $(TARGET)
	$(MAKE) -C mex matlab

octave: $(TARGET)
	$(MAKE) -C mex octave

cleanlib:
	$(RMDIR) $(buildprefix)
	$(RMDIR) $(objprefix)

clean: cleanlib
	$(MAKE) -C mex clean

static: $(TARGET)

shared: $(SO_TARGET)

.PHONY: doc doxy mat2doc mat2docmat clean cleanlib octave matlab cleandoc cleandoxy cleanmat2doc static shared
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

$(buildprefix)/phaseret.h: $(buildprefix)
	$(CC) -E -P -DNOSYSTEMHEADERS -nostdinc include/phaseret.h -o $(buildprefix)/phaseret.h
	sed -i '1 i #ifndef _PHASERET_H' $(buildprefix)/phaseret.h
	sed -i '1 a #define _PHASERET_H' $(buildprefix)/phaseret.h
	sed -i '2 a #include <ltfat.h>' $(buildprefix)/phaseret.h
	sed -i '$$ a #endif' $(buildprefix)/phaseret.h

install: lib
	install -d $(DESTDIR)$(LIBDIR)
	install $(TARGET) $(DESTDIR)$(LIBDIR)
	install $(SO_TARGET) $(DESTDIR)$(LIBDIR)
	mkdir -p $(DESTDIR)$(INCDIR)
	cp -r include/* $(DESTDIR)$(INCDIR)

uninstall:
	rm -f $(DESTDIR)$(LIBDIR)/$(STATIC)
	rm -f $(DESTDIR)$(LIBDIR)/$(SHARED)
	rm -f $(DESTDIR)$(INCDIR)/phaseret.h
	rm -rf $(DESTDIR)$(INCDIR)/phaseret

