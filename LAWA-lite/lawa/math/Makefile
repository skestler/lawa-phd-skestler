include $(LAWA_HOME)/Makefile.common

SOURCES  := $(wildcard *.cpp)
OBJECTS  := $(SOURCES:.cpp=.o)
DEPS     := $(SOURCES:.cpp=.d)

all: $(OBJECTS)

clean:
	$(RM) *.o *.d $(DYLIB)
	$(RMDIR) *.dSYM

# automatic dependencies -------------------------------------------------------
%.d: %.cpp
	$(CXX) -MM $(CPPFLAGS) $(CXXFLAGS) $< \
		| sed 's/$*.o/& $(subst /,\/,$@)/g' > $@

-include $(DEPS)
