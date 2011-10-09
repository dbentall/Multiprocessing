CXX = g++
CC = gcc
CFLAGS = -fopenmp
INCPATHS = -I/usr/X11R6/include -I/usr/local/include\
-Iimlib-1.9.15/Imlib
LIBPATHS = -L/usr/X11R6/lib -L/usr/local/lib -Liml
LIBS = -lfftw3 -lrfftw -lX11 -lXext -ljpeg -lpng -lz -lm -lImlib

all: OMPTest

OMPTest: openMPTest.o wiener.o
	$(CXX) $(CFLAGS) openMPTest.o wiener.o -o OMPTest $(LIBPATHS) $(LIBS)

openMPTest.o: openMPTest.cpp
	$(CXX) $(CFLAGS) -c openMPTest.cpp

wiener.o: wiener.c
	$(CC) $(CFLAGS) -c wiener.c $(INCPATHS)

clean:
	rm -rf *o OMPTest