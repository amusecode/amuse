.PHONY: Makefile
TESTS = stopping_condition $(wildcard *_implementation) 

.PHONY: $(MAKECMDGOALS)
$(MAKECMDGOALS): $(TESTS)

.PHONY: $(TESTS)
$(TESTS):
	$(MAKE) -C $@ $(MAKECMDGOALS)

