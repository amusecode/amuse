AMUSE_DIR?=../../../../..
-include ${AMUSE_DIR}/config.mk

all: dist/distributed-amuse.jar

RM ?= rm
JAVAC ?= javac
JAR ?= jar

dist/distributed-amuse.jar: $(shell find src -name "*.java")
	@echo Compiling Distributed Code
	$(RM) -r tmp dist
	mkdir tmp
	mkdir dist
	$(JAVAC) -g -classpath "lib/ipl/*:lib/jetty/*:lib/xenon/*"  -d tmp $(shell find src -name "*.java")

	@echo Building jar file
	$(JAR) -cf dist/distributed-amuse.jar -C tmp .
	cp lib/ipl/*.jar lib/jetty/*.jar lib/xenon/*.jar dist
	$(RM) -r tmp
	@echo done compiling Distributed code

clean:
	$(RM) -r tmp dist


