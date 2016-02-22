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

export TARGET=build/libphaseret.a
SO_TARGET=$(patsubst %.a,%.so,$(TARGET))

CFLAGS=-std=c11 -pedantic -Wall -Wextra -DNDEBUG -I./include $(OPTFLAGS)
export LIBS=-lm -lfftw3

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

build:
	@mkdir -p build
	@mkdir -p bin

matlab:
	make -C mex matlab

octave:
	make -C mex octave

clean:
	@rm -rf build $(LIBOBJECTS)
	make -C mex clean

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

