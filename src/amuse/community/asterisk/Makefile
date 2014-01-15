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

UPDATE_URL = https://github.com/NLeSC/asterisk/archive/develop.zip
UPDATE_ARCHIVE_DIR = asterisk-develop

#this target intentionally left blank
ALL = 

all: $(ALL)
	#done

code: asterisk_worker_java

java_check:
#Makefile requires these checks to be at the beginning of the line
#Also: make sure the indentation is TABS, not spaces
ifeq ($(JAVA_ENABLED),no)
	$(error "Java support not configured in AMUSE")
endif
ifeq (,$(findstring $(JAVA_VERSION),1.6 1.7))
	$(error "Java 1.6 or 1.7 required for Asterisk, but AMUSE configured with Java $(JAVA_VERSION)")
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
	$(RM) -f *.so *.o *.pyc worker_code.cc worker_code.h 
	$(RM) *~ asterisk_worker asterisk_worker_java worker_code.cc
	$(RM) -f *.jar CodeInterface.java Worker.java binev.data
	$(RM) -rf tmp
	make -C src clean

worker_code.cc: interface.py
	$(CODE_GENERATOR) --type=c interface.py AsteriskInterface -o $@

worker_code.h: interface.py
	$(CODE_GENERATOR) --type=H interface.py AteriskInterface -o $@

CodeInterface.java: interface.py
	$(CODE_GENERATOR) --type=java --mode=interface interface.py AsteriskInterface -o $@

Worker.java: interface.py
	$(CODE_GENERATOR) --type=java --mode=class interface.py AsteriskInterface -o $@

src/lib/asterisk.jar: $(shell find src)
	@echo Compiling Asterisk Lib
	make -C src

worker.jar: Worker.java CodeInterface.java Code.java CodeParticle.java src/lib/asterisk.jar
	$(RM) -r tmp
	mkdir tmp
	$(JAVAC) -g -cp "src/lib/*:src/lib/jogl/*:src/dist/*" -d tmp *.java
	@echo Building jar file
	$(JAR) -cf worker.jar -C tmp .
	$(RM) -rf tmp

asterisk_worker_java: java_check src/lib/asterisk.jar worker.jar
	$(CODE_GENERATOR) -x --type=java --mode=script interface.py AsteriskInterface -o $@

.cc.o: $<
	$(CXX) $(CXXFLAGS) -c -o $@ $< 


