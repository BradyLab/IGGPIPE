# Makefile for findMers

#local
CC=g++

CFLAGS= -O2
LFLAGS= -lm -w -g

INCLUDEDIRS= -I ../stdhdr -I../Strings -I../InputOutput -I .

OBJS= findMers.o ../Strings/Strings.o ../InputOutput/InputOutput.o

findMers:	$(OBJS)
	$(CC) -o findMers $(OBJS) $(LFLAGS)

%.o: %.cpp
	$(CC) $(INCLUDEDIRS) $(CFLAGS) -c $< -o $@

clean:
	find ".." -name "*.o" -print -delete
	rm findMers
