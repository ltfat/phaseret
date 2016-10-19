ifeq ($(COMPTARGET),debug)
	# Do debug stuff here
	CFLAGS +=-O0 -g
	CXXFLAGS +=-O0 -g
endif

ifeq ($(COMPTARGET),release)
	CFLAGS +=-O2 -DNDEBUG
	CXXFLAGS +=-O2 -DNDEBUG
endif

ifeq ($(COMPTARGET),fulloptim)
	CFLAGS +=-O3 -DNDEBUG
	CXXFLAGS +=-O3 -DNDEBUG
endif

