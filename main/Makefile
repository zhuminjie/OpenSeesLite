include ../Makefile.def

all: main sens

main: OpenSeesMain.o

	$(RM) $(RMFLAGS) OpenSeesMain;
	@( \
	for f in $(DIRS); \
	do \
		$(CD) $$f; \
		$(MAKE); \
		$(CD) ..; \
	done );
	$(CD) $(FE)/../main;
	$(LINKER) $(LINKFLAGS) OpenSeesMain.o \
	$(FE_LIBRARY) $(MACHINE_LINKLIBS) \
	$(MACHINE_NUMERICAL_LIBS) $(MACHINE_SPECIFIC_LIBS) \
	 -o OpenSeesMain


sens: OpenSeesSens.o

	$(RM) $(RMFLAGS) OpenSeesSens;
	@( \
	for f in $(DIRS); \
	do \
		$(CD) $$f; \
		$(MAKE); \
		$(CD) ..; \
	done );
	$(CD) $(FE)/../main;
	$(LINKER) $(LINKFLAGS) OpenSeesSens.o \
	$(FE_LIBRARY) $(MACHINE_LINKLIBS) \
	$(MACHINE_NUMERICAL_LIBS) $(MACHINE_SPECIFIC_LIBS) \
	 -o OpenSeesSens


wipe: 
	@( \
	for f in $(DIRS); \
	do \
		$(CD) $$f; \
		$(ECHO) Making lib in $$f; \
		$(MAKE) wipe; \
		$(MAKE) clean; \
		$(CD) ..; \
	done );
	@$(RM) $(RMFLAGS) $(WIPE_LIBS) *.o *~ core 

