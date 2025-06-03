CODES = $(filter-out Makefile interface __init__.py __pycache__, $(wildcard *))


.PHONY: clean distclean
clean distclean: $(CODES)

.PHONY: $(CODES)
$(CODES):
	@# make clean fails for codes that haven't been updated yet, so ignore any errors
	@# for now and don't print them either, there are many.
	@$(MAKE) -C $@ $(MAKECMDGOALS) >/dev/null 2>&1 || true

