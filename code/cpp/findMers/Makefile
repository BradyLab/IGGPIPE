# Makefile for findMers

#local
CC=g++

CFLAGS= -O2
LFLAGS= -lm -w -g

OBJS= findMers.o Strings.o InputOutput.o

# Automatic dependency tracking definitions.

DEPDIR := .d
$(shell mkdir -p $(DEPDIR) >/dev/null)
DEPFLAGS = -MT $@ -MMD -MP -MF $(DEPDIR)/$*.Td

COMPILE.cpp = $(CC) $(DEPFLAGS) $(CFLAGS) -c
POSTCOMPILE = mv -f $(DEPDIR)/$*.Td $(DEPDIR)/$*.d

$(DEPDIR)/%.d: ;

# Targets and recipes.

findMers : $(OBJS)
	$(CC) -o findMers $(OBJS) $(LFLAGS)

%.o : %.cpp $(DEPDIR)/%.d
	$(COMPILE.cpp) -o $@ $<
	$(POSTCOMPILE)

clean:
	find ".." -name "*.o" -print -delete
	rm findMers

# Dependency files.

include $(patsubst %,$(DEPDIR)/%.d,$(basename $(SRCS)))
