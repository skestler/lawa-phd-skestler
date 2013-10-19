include $(LAWA_HOME)/Makefile.common

SOURCEFILES  = $(wildcard *.cpp)
TARGETS      = $(patsubst %.cpp,%,$(SOURCEFILES))

LDFLAGS += -L$(LAWA_HOME) -llawa -lextensionsflens -lextensionssparsegrid
EX_LDFLAGS += -bind_at_load

ifdef GMP_HOME
  CXXFLAGS += -DMPFR -I$(MPFR_HOME)/include
  LDFLAGS += -L$(MPFR_HOME)/lib -lgmpxx -lgmp
endif

all: $(TARGETS) 

clean:
	$(RM) *.o
	$(RM) $(TARGETS)
	$(RMDIR) *.dSYM
