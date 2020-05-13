CC=g++
CFLAGS=-Wall -g
DEPS=pso.h

SOURCE_FILES=$(shell find . -name '*.c')
OBJ_FILES=$(SOURCE_FILES:.c=.o) 
LIB=-lm 

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)
	
demopso: $(OBJ_FILES)
	$(CC) $^ -o $@ $(LIB)


.PHONY: clean
clean:
	rm -f demopso $(OBJ_FILES)
