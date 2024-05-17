.PHONY: Makefile
LIBS = $(filter-out Makefile, $(wildcard *))

.PHONY: $(MAKECMDGOALS)
$(MAKECMDGOALS): $(LIBS)

.PHONY: $(LIBS)
$(LIBS):
	$(MAKE) -C $@ $(MAKECMDGOALS)

