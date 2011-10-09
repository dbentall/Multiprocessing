CXX = g++
CC = gcc
CFLAGS = -c
INCPATHS = -I/v/frooms/fft/include -I/usr/X11R6/include -I/usr/local/include\
-Iimlib-1.9.15/Imlib
LIBPATHS = -L/usr/X11R6/lib -L/usr/local/lib -L
LIBS = -lfftw3 -lX11 -lXext -ljpeg -lpng -lz -lm -lImlib

all: OMPTest

OMPTest: openMPTest.o wiener.o
	$(CXX) openMPTest.o wiener.o -o OMPTest $(LIBPATHS) $(LIBS)

openMPTest.o: openMPTest.cpp
	$(CXX) -fopenmp $(CFLAGS) openMPTest.cpp

wiener.o: wiener.c
	$(CC) $(CFLAGS) wiener.c $(INCPATHS)

clean:
	rm -rf *o OMPTest