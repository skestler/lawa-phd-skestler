include $(LAWA_HOME)/Makefile.common

DYLIB    = libextensionsgmpmath.$(DYLIB_EXT)

SOURCES  := $(wildcard *.cc)
OBJECTS  := $(SOURCES:.cc=.o)
DEPS     := $(SOURCES:.cc=.d)

all: $(DYLIB)

$(DYLIB): $(OBJECTS)
	$(CXX) $(CXXDYLIB) $(CXXFLAGS) $(LDFLAGS) -o $(DYLIB) *.o $(LDLIBS)
	$(CP) $(DYLIB) $(LAWA_HOME)

clean:
	$(RM) *.o *.d
	$(RMDIR) *.dSYM

# automatic dependencies -------------------------------------------------------
%.d: %.cc
	$(CXX) -MM $(CPPFLAGS) $(CXXFLAGS) $< \
		| sed 's/$*.o/& $(subst /,\/,$@)/g' > $@

-include $(DEPS)
