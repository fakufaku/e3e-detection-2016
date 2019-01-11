#
# Makefile to compile the pyramic demo and tests
#
# Author: Robin Scheibler

CC := g++
DEBUG := -g -Wall -O0

# Speed options
SPEEDFLAGS = -O3 -ffast-math -ftree-vectorize -funroll-loops
SPEEDFLAGS_ARM = -mcpu=cortex-a9 -ftree-loop-ivcanon -mfpu=neon -mfloat-abi=hard
ifneq ($(NO_PYRAMIC),1)
	SPEEDFLAGS:=$(SPEEDFLAGS) $(SPEEDFLAGS_ARM)
endif

# Compiler flags
CPPFLAGS := -std=c++14 -lfftw3f
ifeq ($(DEBUG),1)
	CPPFLAGS:=$(CPPFLAGS) $(DEBUG)
else
	CPPFLAGS:=$(CPPFLAGS) $(SPEEDFLAGS)
endif

# Linker flags
LDFLAGS := -L "./lib"
ifeq ($(NO_PYRAMIC),1)
	LIB := -lfftw3f -lpthread
else
	LIB := -lfftw3f -lpyramicio -lpthread
endif
INC := -I./include

SRCDIR := src
PYM_SRCDIR := src_pyramic
BUILDDIR := build
SRCEXT := cpp

SOURCES := $(shell find $(SRCDIR) -type f | grep \.$(SRCEXT)$$) #stft.cpp srpphat.cpp windows.cpp
OBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))

PYM_SOURCES := $(shell find $(PYM_SRCDIR) -type f | grep \.$(SRCEXT)$$) #stft.cpp srpphat.cpp windows.cpp
PYM_OBJECTS := $(patsubst $(PYM_SRCDIR)/%,$(BUILDDIR)/%,$(PYM_SOURCES:.$(SRCEXT)=.o))

TESTS := $(shell find tests -type f | grep \.$(SRCEXT)$$ | cut -f 1 -d '.')
PYM_TESTS := $(shell find tests_pyramic -type f | grep \.$(SRCEXT)$$ | cut -f 1 -d '.')
DEMOS := $(shell find main -type f | grep \.$(SRCEXT)$$ | cut -f 1 -d '.')
PYM_DEMOS := $(shell find demos -type f | grep \.$(SRCEXT)$$ | cut -f 1 -d '.')

ifneq ($(NO_PYRAMIC),1)
	SOURCES := $(SOURCES) $(PYM_SOURCES)
	OBJECTS := $(OBJECTS) $(PYM_OBJECTS)
	TESTS := $(TESTS) $(PYM_TESTS)
	DEMOS := $(DEMOS) $(PYM_DEMOS)
endif

hello:
	@echo NO_PYRAMIC=$(NO_PYRAMIC)
	@echo $(OBJECTS)
	@echo $(TESTS)
	@echo $(DEMOS)

$(BUILDDIR)/%.o: $(SRCDIR)/%.$(SRCEXT)
	mkdir -p build
	$(CC) -c -o $@ $< $(LDFLAGS) $(INC) $(LIB) $(CPPFLAGS)

$(BUILDDIR)/%.o: $(PYM_SRCDIR)/%.$(SRCEXT)
	mkdir -p build
	$(CC) -c -o $@ $< $(LDFLAGS) $(INC) $(LIB) $(CPPFLAGS)

$(TESTS): $(OBJECTS)
	mkdir -p bin/$(shell dirname $@)
	$(CC) $@.cpp -o bin/$@ $^ $(LDFLAGS) $(INC) $(LIB) $(CPPFLAGS)

$(DEMOS): $(OBJECTS)
	mkdir -p bin/$(shell dirname $@)
	$(CC) $@.cpp -o bin/$@ $^ $(LDFLAGS) $(INC) $(LIB) $(CPPFLAGS)

objects: $(OBJECTS)
tests: $(TESTS)
demos: $(DEMOS)

all: demos tests

clean:
	rm -f bin/* build/*
