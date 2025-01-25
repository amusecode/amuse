.PHONY: Makefile
TESTS = stopping_condition async_implementation c_implementation
TESTS += fortran_implementation grid_implementation

.PHONY: $(MAKECMDGOALS)
$(MAKECMDGOALS): $(TESTS)

.PHONY: $(TESTS)
$(TESTS):
	$(MAKE) -C $@ $(MAKECMDGOALS)

