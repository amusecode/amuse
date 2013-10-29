# standard amuse configuration include
# config.mk will be made after ./configure has run
AMUSE_DIR?=../../..
-include $(AMUSE_DIR)/config.mk

MPICXX   ?= mpicxx

CFLAGS   += -Wall -g
CXXFLAGS += $(CFLAGS) 
LDFLAGS  += -lm $(MUSE_LD_FLAGS)

OBJS = interface.o

CODE_GENERATOR = $(AMUSE_DIR)/build.py

RM ?= rm
JAVAC ?= javac
JAR ?= jar
MV ?= mv
CP ?= cp

UPDATE_URL = https://github.com/NLeSC/AMUSE-Distributed/archive/develop.zip
UPDATE_ARCHIVE_DIR = AMUSE-Distributed-develop

#this target intentionally left blank
ALL = 

all: $(ALL)
	#dummy target, use distributed_worker_java to build

java_check:
#Makefile requires these checks to be at the beginning of the line
#Also: make sure the indentation is TABS, not spaces
ifeq ($(JAVA_ENABLED),no)
	$(error "Java support not configured in AMUSE")
endif
ifeq (,$(findstring $(JAVA_VERSION), 1.7))
	$(error "Java 1.7 required for Asterisk, but AMUSE configured with Java $(JAVA_VERSION)")
endif

update:
	$(RM) -rf code.zip $(UPDATE_ARCHIVE_DIR)
	$(MV) -T src src.old
	wget $(UPDATE_URL) -O code.zip
	unzip code.zip -x "*/.*"
	$(MV) $(UPDATE_ARCHIVE_DIR) src
	$(CP) Makefile.code src/Makefile
	$(RM) code.zip

clean:
	$(RM) -f *.so *.o *.pyc 
	$(RM) *~ distributed distributed_worker_java 
	$(RM) -f *.jar CodeInterface.java Worker.java binev.data
	$(RM) -rf tmp
	make -C src clean

CodeInterface.java: interface.py
	$(CODE_GENERATOR) --type=java --mode=interface interface.py DistributedAmuseInterface -o $@

Worker.java: interface.py
	$(CODE_GENERATOR) --type=java --mode=class interface.py DistributedAmuseInterface -o $@

src/dist/distributed-amuse.jar: $(shell find src -name "*.java")
	@echo Compiling Distributed Code
	make -C src

worker.jar: Worker.java CodeInterface.java Code.java src/dist/distributed-amuse.jar
	$(RM) -r tmp
	mkdir tmp
	$(JAVAC) -g -cp "src/dist/*" -d tmp Worker.java CodeInterface.java Code.java
	@echo Building jar file
	$(JAR) -cf worker.jar -C tmp .
	$(RM) -rf tmp

distributed_worker_java: interface.py java_check worker.jar src/dist/distributed-amuse.jar
	$(CODE_GENERATOR) -x --type=java --mode=script interface.py DistributedAmuseInterface -o $@
