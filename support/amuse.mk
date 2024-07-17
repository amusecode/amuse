ifeq (,$(filter clean distclean, $(MAKECMDGOALS)))


# Detect what features we have
support/features.mk: support/features.mk.in support/configure
	cd support && ./configure $(CONFIGOPTS) || cat config.log

include support/features.mk
include support/format.mk


# See if we can install the framework
FW_REQ_FEATURES := c c++ fortran python install mpi
FW_MISSING_FEATURES := $(filter-out $(FEATURES), $(FW_REQ_FEATURES))

# See if we can install Sapporo light
SAPPORO_LIGHT_REQ_FEATURES := c c++ install cuda
SAPPORO_LIGHT_MISSING_FEATURES := $(filter-out $(FEATURES), $(SAPPORO_LIGHT_REQ_FEATURES))

# Lists of enabled and disabled packages, filled by the .amuse_dep.mk files
ifeq (,$(FW_MISSING_FEATURES))
ENABLED_PACKAGES := \namuse-framework
DISABLED_PACKAGES :=
else
ENABLED_PACKAGES :=
DISABLED_PACKAGES := \namuse-framework (missing features: $(COLOR_RED)$(FW_MISSING_FEATURES)$(COLOR_END))
endif

ifeq (,$(SAPPORO_LIGHT_MISSING_FEATURES))
ENABLED_PACKAGES += \nsapporo_light
else
DISABLED_PACKAGES += \nsapporo_light (missing features: $(COLOR_RED)$(SAPPORO_LIGHT_MISSING_FEATURES)$(COLOR_END))
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
#
# Note: I tried using the file function instead of shell cat, but on GNU make 4.3 at
# least, it doesn't strip the final newline as the docs say. The shell function strips
# all newlines, which is better anyway as it's more robust.
support/comm_deps_mk/%.mk: src/amuse/community/%
	@mkdir -p $(dir $@)
	@rm -f $@
	@printf '%s\n' 'include support/format.mk' >>$@
	@printf '\n' >>$@
	@printf '%s\n' 'PACKAGE_NAME = $(notdir $(patsubst %.amuse_deps,%,$*))' >>$@
	@printf '\n' >>$@
	@printf '%s\n' 'REQUIRED_FEATURES := $(shell cat $<)' >>$@
	@printf '%s\n' 'MISSING_FEATURES := $$(filter-out $$(FEATURES), $$(REQUIRED_FEATURES))' >>$@
	@printf '\n' >>$@
	@printf '%s\n' 'ifeq (,$$(MISSING_FEATURES))' >>$@
	@printf '\n' >>$@
	@printf '%s\n' 'ENABLED_PACKAGES += \n$(notdir $(patsubst %.amuse_deps,%,$*))' >>$@
	@printf '\n' >>$@
	@printf '%s\n' '%-$(notdir $(patsubst %.amuse_deps,%,$*)):' >>$@
	@printf '\t%s\n' 'make -C $(dir $<)/.. $$@' >>$@
	@printf '\n' >>$@
	@printf '%s\n' 'develop-packages: develop-$(notdir $(patsubst %.amuse_deps,%,$*))' >>$@
	@printf '\n' >>$@
	@printf '%s\n' 'install-packages: install-$(notdir $(patsubst %.amuse_deps,%,$*))' >>$@
	@printf '\n' >>$@
	@printf '%s\n' 'else' >>$@
	@printf '\n' >>$@
	@printf '%s\n' 'DISABLED_PACKAGES += \n$(notdir $(patsubst %.amuse_deps,%,$*)) (missing features: $(COLOR_RED)$$(MISSING_FEATURES)$(COLOR_END))' >>$@
	@printf '%s\n' 'endif' >>$@


endif               # target is not clean or distclean

