CC=g++
CFLAGS=-Wall -g

default: all
all: clean demo
demo: demo.c pso.c
clean:
	rm demo
