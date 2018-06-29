include ../Makefile.def

all: main

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
