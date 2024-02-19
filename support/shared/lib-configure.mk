# Run configure if needed
ifeq (,$(filter clean distclean, $(MAKECMDGOALS)))

include support/config.mk

support/config.mk: support/config.mk.in
	cd support && ./configure $(CONFIGOPTS)

endif

