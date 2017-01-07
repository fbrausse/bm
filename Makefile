PKGS    := cairo sdl2 glew gl
CXXFLAGS = -std=c++11 $(shell pkg-config --cflags $(PKGS))
LDFLAGS  = $(shell pkg-config --ldflags $(PKGS))
LDLIBS   = $(shell pkg-config --libs $(PKGS))

bm3: bm3.cc
