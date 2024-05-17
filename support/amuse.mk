ifeq (,$(filter clean distclean, $(MAKECMDGOALS)))


# Detect what features we have
support/features.mk: support/features.mk.in support/configure
	cd support && ./configure $(CONFIGOPTS)

include support/features.mk


# See if we can install the framework
FW_REQ_FEATURES := gcc g++ gfortran python install mpi
FW_MISSING_FEATURES := $(filter-out $(FEATURES), $(FW_REQ_FEATURES))

# Lists of enabled and disabled packages, filled by the .amuse_dep.mk files
ifeq (,$(FW_MISSING_FEATURES))
ENABLED_PACKAGES := \namuse-framework
DISABLED_PACKAGES :=
else
ENABLED_PACKAGES :=
DISABLED_PACKAGES := \namuse-framework (missing features: $(COLOR_RED)$(FW_MISSING_FEATURES)$(COLOR_END))
endif


# Detect community codes and their required features
COMMUNITY_CODES := $(filter-out Makefile, $(subst src/amuse/community/,,$(wildcard src/amuse/community/*)))

# List of worker metadata files for all community codes
COMM_DEPS := $(foreach d,$(wildcard src/amuse/community/*), $(wildcard $d/packages/*.amuse_deps))

# List of makefile snippets produced from the metadata files
COMM_DEPS_MK := $(patsubst src/amuse/community/%, support/comm_deps_mk/%.mk, $(COMM_DEPS))

# Include the makefile snippets to complete the packages targets
include $(COMM_DEPS_MK)

# Creates a makefile snippet that adds targets for the package, and adds those targets
# as a prerequisites of the 'develop-packages' and 'install-packages' targets, but only
# if FEATURES contains all the dependencies of the worker. As a result, when making
# 'install-packages' we will only try to build the packages that we have the
# dependencies for and expect to build cleanly.
support/comm_deps_mk/%.mk: src/amuse/community/%
	@mkdir -p $(dir $@)
	@rm -f $@
	@echo include support/format.mk >>$@
	@echo >>$@
	@echo PACKAGE_NAME = $(notdir $(patsubst %.amuse_deps,%,$*)) >>$@
	@echo >>$@
	@echo 'REQUIRED_FEATURES := $(file < $<)' >>$@
	@echo 'MISSING_FEATURES := $$(filter-out $$(FEATURES), $$(REQUIRED_FEATURES))' >>$@
	@echo >>$@
	@echo 'ifeq (,$$(MISSING_FEATURES))' >>$@
	@echo >>$@
	@echo 'ENABLED_PACKAGES += \\n$(notdir $(patsubst %.amuse_deps,%,$*))' >>$@
	@echo >>$@
	@echo '%-$(notdir $(patsubst %.amuse_deps,%,$*)):' >>$@
	@echo '\tmake -C $(dir $<)' $$\@ >>$@
	@echo >>$@
	@echo 'develop-packages: develop-$(notdir $(patsubst %.amuse_deps,%,$*))' >>$@
	@echo >>$@
	@echo 'install-packages: install-$(notdir $(patsubst %.amuse_deps,%,$*))' >>$@
	@echo >>$@
	@echo else >>$@
	@echo >>$@
	@echo 'DISABLED_PACKAGES += \\n$(notdir $(patsubst %.amuse_deps,%,$*)) (missing features: $(COLOR_RED)$$(MISSING_FEATURES)$(COLOR_END))' >>$@
	@echo endif >>$@


endif               # target is not clean or distclean

