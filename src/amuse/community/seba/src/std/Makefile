include ../Makefile.inc

SRC = $(shell ls *.C)
OBJ = $(SRC:.C=.o)
LIB = libstd.a

lib:	$(LIB)

$(LIB):	$(OBJ)
	ar -cr $(LIB) $(OBJ)

clean:
	$(RM) $(OBJ) $(LIB) *~
