all: worker.jar java_worker


# TestInterface
%.class: %.java
	echo ${CLASSPATH}
	$(JAVAC) $<

Code.class: CodeInterface.class

CodeInterface.java: interface.py
	amusifier --mode=interface --type=java -o $@ $< ForTestingInterface

Worker.java: interface.py
	amusifier --mode=class --type=java -o $@ $< ForTestingInterface

# Bit of a hack with the dollar sign here, sorry
worker.jar: Code.class CodeInterface.class Worker.class
	$(JAR) -cf $@ $^ Worker\$$AmuseMessage.class

java_worker: interface.py
	amusifier --mode=script --type=java -o $@ $< ForTestingInterface
	chmod +x $@


.PHONY: clean
clean:
	rm -f CodeInterface.java Worker.java *.class *.jar *_worker

.PHONY: distclean
distclean: clean

