include ../../Makefile.inc

SRC = $(shell ls *.C)
OBJ = $(SRC:.C=.o)
LIB = ../libsstar.a

lib:	$(LIB)($(OBJ))

$(LIB)($(OBJ)):	$(OBJ)
	ar -cr $(LIB) $(OBJ)

clean:
	$(RM) $(OBJ) *~
