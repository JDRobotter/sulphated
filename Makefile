CC = gcc

COMPOPT = -O3 -Wall -g

CFLAGS = $(COMPOPT)

LDFLAGS = `sdl-config --libs` -lm

OBJS = main.o fluid.o loutron.o

SOURCE = .

EXEC = sulphate

all: $(OBJS)
	$(CC) $(OBJS) $(COMPOPT) $(LDFLAGS) -o $(EXEC);
	./upversion.pl > /tmp/loutre; mv /tmp/loutre ./version.h

%.o:$(SOURCE)/%.c $(SOURCE)/%.h
	$(CC) $(CFLAGS) $< -c -o $@

%.o:$(SOURCE)/%.c
	$(CC) $(CFLAGS) $< -c -o $@


clean:
	rm -f $(OBJS) $(SOURCE)/*~ $(SOURCE)/\#*\# $(EXEC)

distclean:
	rm -f $(OBJS) $(SOURCE)/*~ $(SOURCE)/\#*\#
