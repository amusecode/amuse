.DEFAULT_GOAL = configure

ifeq (,$(MAKE_RESTARTS))
$(info Reconfiguring...)
$(shell rm -f support/features.mk)
endif


ifeq (,$(filter clean distclean, $(MAKECMDGOALS)))
 
include support/amuse.mk
include support/dependencies.mk
include support/format.mk
include support/shared/venv.mk
# help depends on the others, so it must be last
include support/help.mk

endif


.PHONY: configure
configure:
	@echo
	@echo
	@printf '%b\n' '$(COLOR_CYAN)*** Configuration complete ***$(COLOR_END)'
	@echo
	@printf 'Detected features: $(FEATURES)\n'
	@echo
	@printf '%b\n' '$(COLOR_GREEN)** Enabled packages **$(COLOR_END)'
	@printf '$(ENABLED_PACKAGES)\n'
	@echo
	@printf '%b\n' '$(COLOR_RED)** Disabled packages **$(COLOR_END)'
	@printf '$(DISABLED_PACKAGES)\n'
	@echo
	@printf '%b\n' '$(COLOR_CYAN)*** Next steps ***$(COLOR_END)'
ifeq (,$(ENV_TYPE))
	@printf '%b\n' '$(ENVIRONMENT_HELP)'
else
ifneq (,$(NEED_PIP_WHEEL))
	@printf '%b\n' '$(NO_PIP_WHEEL_MESSAGE)'
else
ifneq (,$(FW_MISSING_FEATURES)$(DISABLED_PACKAGES))
	@printf '%b\n' '$(DISABLED_PACKAGES_MESSAGE)'
endif
ifeq (,$(NEED_PIP_WHEEL))
ifneq (,$(ENABLED_PACKAGES))
	@printf '%b\n' '$(INSTALL_HELP)'
endif
endif
endif
endif

ifeq (,$(FW_MISSING_FEATURES))

# Build sapporo_light separately so that we can build without CUDA
ifneq (,$(filter cuda, $(FEATURES)))
HAVE_CUDA := yes
export HAVE_CUDA

.PHONY: install-sapporo_light
install-sapporo_light:
	$(MAKE) -C lib install-sapporo_light
endif


.PHONY: install-libs
install-libs:
	$(MAKE) -C lib install

.PHONY: install-amuse-framework
install-amuse-framework: install-libs
	support/shared/uninstall.sh amuse-framework
	cd src && pip --no-cache-dir --debug install .

.PHONY: develop-amuse-framework
develop-amuse-framework: install-libs
	support/shared/uninstall.sh amuse-framework
	cd src && pip install -e .

.PHONY: package-amuse-framework
package-amuse-framework: install-libs
	cd src && python3 -m pip install -vv --no-cache-dir --no-deps --no-build-isolation --prefix ${PREFIX} .

endif


.PHONY: test-amuse-framework
test-amuse-framework:
	$(MAKE) -C src/tests all
	cd src/tests && pytest --pyargs core_tests compile_tests


.PHONY: install-packages

.PHONY: develop-packages

.PHONY: install
install: install-amuse-framework install-packages


.PHONY: clean
clean:
	$(MAKE) -C support clean
	$(MAKE) -C lib clean
	$(MAKE) -C src/amuse/community clean
	$(MAKE) -C src/tests clean

.PHONY: distclean
distclean:
	$(MAKE) -C support distclean
	$(MAKE) -C lib distclean
	$(MAKE) -C src/amuse/community distclean
	$(MAKE) -C src/tests distclean

