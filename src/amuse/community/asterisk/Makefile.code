AMUSE_DIR?=../../../../..
-include ${AMUSE_DIR}/config.mk

all: dist/asterisk.jar

RM ?= rm
JAVAC ?= javac
JAR ?= jar

dist/asterisk.jar: $(shell find src -name "*.java")
	@echo Compiling Asterisk
	$(RM) -r tmp dist
	mkdir tmp
	mkdir dist
	$(JAVAC) -g -classpath "lib/*:lib/jogl/*"  -d tmp $(shell find src -name "*.java")

	@echo Building jar file
	$(JAR) -cf dist/asterisk.jar -C tmp .
	$(JAR) -uf dist/asterisk.jar images shaders
	$(RM) -r tmp
	@echo done compiling Asterisk

clean:
	$(RM) -r tmp dist


