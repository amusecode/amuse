.PHONY: Makefile
AMUSE_LIBS = $(filter-out Makefile sapporo_light sapporo_2, $(wildcard *))


# We build sapporo_light separately and only on demand, because it needs CUDA and we
# may not have that.
.PHONY: %-sapporo_light
%-sapporo_light:
	$(MAKE) -C sapporo_light $*


# We build sapporo2 separately and only on demand, because it needs CUDA and we
# may not have that.
.PHONY: %-sapporo2
%-sapporo2:
	$(MAKE) -C sapporo_2 $*


.PHONY: install uninstall clean distclean
install uninstall clean distclean: $(AMUSE_LIBS)

.PHONY: $(AMUSE_LIBS)
$(AMUSE_LIBS):
	$(MAKE) -C $@ $(MAKECMDGOALS)

