include $(LAWA_HOME)/Makefile.common
SUBDIRS := lawa extensions tutorials applications validate

all:
	$(SILENT) for i in $(SUBDIRS); \
		do echo ""; echo "processing dir $$i"; \
		$(MAKE) -C $$i all; \
	done;
	
clean:
	$(SILENT) for i in $(SUBDIRS) ; \
		do echo ""; echo "processing dir $$i"; \
		$(MAKE) -C $$i clean; \
	done;
	$(RM) libextensionsflens.$(DYLIB_EXT) \
              libextensionssparsegrid.$(DYLIB_EXT) \
              liblawa.$(DYLIB_EXT)
	$(RM) *.tmp
# DO NOT DELETE
