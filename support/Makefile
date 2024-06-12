.PHONY: clean
clean: clean_features

.PHONY: distclean
distclean: clean clean_comm_deps

.PHONY: clean_features
clean_features:
	rm -f features.mk config.status config.log
	rm -rf autom4te.cache

.PHONY: clean_comm_deps
clean_comm_deps:
	rm -rf comm_deps_mk/*

