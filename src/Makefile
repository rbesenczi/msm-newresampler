include ${FSLCONFDIR}/default.mk

PROJNAME = newresampler
SOFILES  = libfsl-newresampler.so
LIBS     = -lfsl-miscmaths -lfsl-giftiio -lfsl-newmeshclass
OBJS     = meshException.h point.o mpoint.o triangle.o\
		   mesh.o node.o octree.o resampler.o
USRCXXFLAGS += -Wall

all: libfsl-newresampler.so

libfsl-newresampler.so: ${OBJS}
	${CXX} ${CXXFLAGS} -shared -o $@ $^ ${LDFLAGS}
