CC ?= gcc
AR ?= ar

VERSION := $(shell cat phaseret_version)
PACKAGE = phaseret-$(VERSION)

PREFIX ?= /usr/local
LIBDIR = $(PREFIX)/lib
INCDIR = $(PREFIX)/include

# LIB
LIBSOURCES = $(wildcard src/*.c)
LIBOBJECTS = $(addprefix obj/,$(notdir $(LIBSOURCES:.c=.o)))

STATIC = libphaseret.a
SHARED = libphaseret.so
export TARGET=$(addprefix build/,$(STATIC))
SO_TARGET= $(addprefix build/,$(SHARED))

CFLAGS=-std=c11 -pedantic -Wall -Wextra -DNDEBUG -I./include/phaseret $(OPTFLAGS)
export LIBS=-lm -lfftw3

lib: $(TARGET) $(SO_TARGET)

all: lib matlab octave

dev: CFLAGS=-std=c11 -g -O0 -Wall -Wall -Wextra -I./include/phaseret $(OPTFLAGS)
dev: all octave

$(TARGET): CFLAGS += -fPIC
$(TARGET): obj build $(LIBOBJECTS)
	ar rvu $@ $(LIBOBJECTS)
	ranlib $@

$(SO_TARGET): $(TARGET) $(LIBOBJECTS)
	$(CC) -shared -Wl,--no-undefined -o $@ $(LIBOBJECTS) $(LIBS)

obj/%.o: src/%.c obj
	$(CC) -c $(CFLAGS) $< -o $@

build:
	@mkdir -p build
	@mkdir -p bin

obj:
	@mkdir -p obj

matlab: $(TARGET)
	make -C mex matlab

octave: $(TARGET)
	make -C mex octave

clean:
	@rm -rf build bin obj
	make -C mex clean

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

