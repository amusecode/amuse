AMUSE_DIR?=../..
-include ${AMUSE_DIR}/config.mk

LIBPROFILING = src/profiling/libibisprofiling.a

JAVAC ?= /usr/bin/javac
JAVA ?= /usr/bin/java
JAVAH ?= /usr/bin/javah
JAR ?= /usr/bin/jar

CLASSPATH="deploy/lib/*:tmp"

all: deploy ${LIBPROFILING} lib/ibis-amuse.jar

$(LIBPROFILING):src/profiling/*.c
ifeq ($(MPI_ENABLED),yes)
	make -C src/profiling all
endif

lib/ibis-amuse.jar: src/*/*/*.java
	@echo Compiling Amuse/Ibis lib
	$(RM) -r lib tmp
	mkdir tmp
	mkdir lib
	$(JAVAC) -classpath $(CLASSPATH) -d tmp src/*/*/*.java
	#$(JAVA) -classpath $(CLASSPATH) ibis.compile.Ibisc tmp
	$(JAVAH) -classpath $(CLASSPATH) -d include ibis.amuse.JNICodeInterface
	@echo Building jar file
	$(JAR) -cf lib/ibis-amuse.jar -C tmp .
	$(JAR) -uf lib/ibis-amuse.jar images/*
	$(RM) -r tmp
	@echo done compiling Amuse/Ibis lib

clean:
	make -C src/profiling clean
	$(RM) -r lib tmp
	
distclean: clean

