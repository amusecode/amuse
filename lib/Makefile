.PHONY: Makefile
LIBS = $(filter-out Makefile sapporo_light, $(wildcard *))


# We build sapporo_light separately and only on demand, because it needs CUDA and we
# may not have that.
.PHONY: %-sapporo_light
%-sapporo_light:
	$(MAKE) -C sapporo_light $*


.PHONY: install uninstall clean distclean
install uninstall clean distclean: $(LIBS)

.PHONY: $(LIBS)
$(LIBS):
	$(MAKE) -C $@ $(MAKECMDGOALS)

