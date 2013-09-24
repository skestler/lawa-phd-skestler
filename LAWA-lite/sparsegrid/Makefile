include $(LAWA_HOME)/Makefile.common

DYLIB    = libextensionssparsegrid.$(DYLIB_EXT)

SOURCES  := sandia_rules.cpp sparse_grid_mixed.cpp sparse_grid_open.cpp
OBJECTS  := $(SOURCES:.cpp=.o)
DEPS     := $(SOURCES:.cpp=.d)

all: $(DYLIB)

$(DYLIB): $(OBJECTS)
	$(CXX) $(CXXDYLIB) $(CXXFLAGS) $(LDFLAGS) -o $(DYLIB) $(OBJECTS) $(LDLIBS)
	$(CP) $(DYLIB) $(LAWA_HOME)

clean:
	$(RM) *.o *.d $(DYLIB)
	$(RMDIR) *.dSYM


# automatic dependencies -------------------------------------------------------
%.d: %.cpp
	$(CXX) -MM $(CPPFLAGS) $(CXXFLAGS) $< \
		| sed 's/$*.o/& $(subst /,\/,$@)/g' > $@

-include $(DEPS)
