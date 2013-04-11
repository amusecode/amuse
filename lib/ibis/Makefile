AMUSE_DIR?=../..
-include ${AMUSE_DIR}/config.mk

VERSION = -20120709-snapshot

CODEDIR = deploy${VERSION}
DOWNLOAD_FROM_WEB = ./downloadhttp.py
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

deploy: ${CODEDIR}
	ln -fs ${CODEDIR} deploy

$(CODEDIR):
ifdef DOWNLOAD_CODES
	$(DOWNLOAD_FROM_WEB) ${VERSION}
else
	@echo ""
	@echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	@echo ""
	@echo "DOWNLOAD_CODES is not set. Ibis will not be downloaded and build."
	@echo "If you do want Ibis, set DOWNLOAD_CODES to 1."
	@echo "bash> export DOWNLOAD_CODES=1"
	@echo "csh> setenv DOWNLOAD_CODES 1"
	@echo ""
	@echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	@echo ""
	@make -s --no-print-directory -C . raise_error
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
	$(RM) -r deploy
	$(RM) -r deploy-*
	

