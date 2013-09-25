SUBDIRS := unbounded_domains new_eval_scheme canutoproject
ifdef BOOST_HOME 
ifdef MPREC_INC_HOME 
ifdef MPREC_LIB_HOME
    SUBDIRS += finance
endif
endif
endif

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
# DO NOT DELETE
