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

CFLAGS=-std=c11 -pedantic -Wall -Wextra -DNDEBUG -I./include/phaseret -Ithirdparty $(OPTFLAGS)

STATIC = libphaseret.a
ifdef MINGW
	SHARED = libphaseret.dll
	EXTRALFLAGS = -Wl,--out-implib,$@.a -static-libgcc
else
	CFLAGS += -fPIC
	SHARED = libphaseret.so
endif
export TARGET=$(addprefix build/,$(STATIC))
SO_TARGET= $(addprefix build/,$(SHARED))

FFTWLIB ?= -lfftw3

export LIBS=-lm $(FFTWLIB)

lib: $(TARGET) $(SO_TARGET)

all: lib matlab octave

dev: CFLAGS=-std=c11 -g -O0 -Wall -Wall -Wextra -I./include/phaseret $(OPTFLAGS)
dev: all octave

$(TARGET): obj build $(LIBOBJECTS)
	$(AR) rvu $@ $(LIBOBJECTS)
	$(RANLIB) $@

$(SO_TARGET): $(TARGET) $(LIBOBJECTS)
	$(CC) -shared -Wl,--no-undefined -o $@ $(LIBOBJECTS) $(EXTRALFLAGS) $(LIBS)

obj/%.o: src/%.c obj
	$(CC) -c $(CFLAGS) $< -o $@

build:
	$(MKDIR) build
	$(MKDIR) bin

obj:
	$(MKDIR) obj

matlab: $(TARGET)
	$(MAKE) -C mex matlab

octave: $(TARGET)
	$(MAKE) -C mex octave

clean:
	$(RMDIR) build
	$(RMDIR) bin
	$(RMDIR) obj
	$(MAKE) -C mex clean

.PHONY: doc doxy mat2doc mat2docmat clean octave matlab cleandoc cleandoxy cleanmat2doc
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

