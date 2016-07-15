ifdef CROSS
CC = $(CROSS)gcc
AR = $(CROSS)ar
RANLIB =ranlib
MINGW=1
else
CC ?= gcc
AR ?= ar
RANLIB?=ranlib
endif

include ostools.mk

# VERSION := $(shell cat phaseret_version)
PACKAGE = phaseret-$(VERSION)

PREFIX ?= /usr/local
LIBDIR = $(PREFIX)/lib
INCDIR = $(PREFIX)/include

# LIB
LIBSOURCES = $(wildcard src/*.c)
LIBOBJECTS = $(addprefix obj/,$(notdir $(LIBSOURCES:.c=.o)))

CFLAGS=-std=c99 -pedantic -O2 -Wall -Wextra -DNDEBUG -I./include -Ithirdparty $(OPTFLAGS)

STATIC = libphaseret.a
ifdef MINGW
	SHARED = libphaseret.dll
	EXTRALFLAGS = -Wl,--out-implib,$@.a -static-libgcc
else
	CFLAGS += -fPIC
	SHARED = libphaseret.so
endif
TARGET=$(addprefix build/,$(STATIC))
SO_TARGET= $(addprefix build/,$(SHARED))

FFTWLIB ?= -lfftw3

LIBS=-lltfat  $(FFTWLIB) -lm

lib: $(TARGET) $(SO_TARGET)

all: lib matlab octave

dev: CFLAGS=-std=c99 -g -O0 -Wall -Wall -Wextra -I./include $(OPTFLAGS)
dev: all 

$(TARGET): obj build $(LIBOBJECTS)
	$(AR) rvu $@ $(LIBOBJECTS)
	$(RANLIB) $@

$(SO_TARGET): obj build $(LIBOBJECTS)
	$(CC) -shared -Wl,--no-undefined -o $@ $(LIBOBJECTS) $(EXTRALFLAGS) $(LIBS)

obj/%.o: src/%.c obj
	$(CC) -c $(CFLAGS) $< -o $@

build:
	$(MKDIR) build

obj:
	$(MKDIR) obj

matlab: $(TARGET)
	$(MAKE) -C mex matlab

octave: $(TARGET)
	$(MAKE) -C mex octave

cleanlib:
	$(RMDIR) build
	$(RMDIR) obj

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

build/phaseret.h: build 
	$(CC) -E -P -DNOSYSTEMHEADERS -nostdinc include/phaseret.h -o build/phaseret.h
	sed -i '1 i\#include <ltfat.h>' build/phaseret.h

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

