.PHONY: Makefile
LIBS = $(filter-out Makefile sapporo_light, $(wildcard *))


# We build sapporo_light separately and only on demand, because it needs CUDA and we
# may not have that.
.PHONY: install-sapporo_light
install-sapporo_light:
	$(MAKE) -C sapporo_light install


.PHONY: $(MAKECMDGOALS)
$(MAKECMDGOALS): $(LIBS)

.PHONY: $(LIBS)
$(LIBS):
	$(MAKE) -C $@ $(MAKECMDGOALS)

