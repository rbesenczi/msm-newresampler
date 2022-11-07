include ${FSLCONFDIR}/default.mk

PROJNAME = newresampler-demo
XFILES   = metric-resample surface-resample 
LIBS     = -lfsl-newresampler -lboost_program_options
USRCXXFLAGS += -Wall

all: ${XFILES}

%: %.o
	${CXX} ${CXXFLAGS} -o $@ $^ ${LDFLAGS}
