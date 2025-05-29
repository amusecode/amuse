ifneq (,$(filter-out clean distclean, $(MAKECMDGOALS)))
include support/config.mk
endif

include support/shared/venv.mk

LD_LIBRARY_PATH += $(ENV_LIBRARY_PATH)
PYTHONPATH += $(CURDIR)

export

# Detecting dependencies
support/config.mk: support/config.mk.in support/configure
	cd support && ./configure


# Only the compile tests and the speed test need compiling, forward there
# Target 'all' will build all the workers for all the tests
ifneq (,$(filter-out distclean, $(MAKECMDGOALS)))
.PHONY: $(MAKECMDGOALS)
$(MAKECMDGOALS):
	$(MAKE) -C compile_tests $(MAKECMDGOALS)
	$(MAKE) -C reports $(MAKECMDGOALS)
endif


.PHONY: distclean
distclean:
	$(MAKE) -C compile_tests distclean
	rm -f support/config.mk support/config.status support/config.log support/configure~
	rm -rf support/autom4te.cache

