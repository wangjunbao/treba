CC = gcc
BINPREFIX = /usr/local/bin/
MANPREFIX = /usr/local/share/man/man1/
RM = /bin/rm -f
CFLAGS = -O3 -ffast-math
LFLAGS = -lm -lpthread

all:	treba

treba:	treba.c treba.h fastlogexp.h
	$(CC) $(CFLAGS) treba.c -o treba $(LFLAGS)

clean:
	$(RM) treba

install: treba
	-@if [ ! -d $(BINPREFIX) ]; then mkdir -p $(BINPREFIX); fi
	-@if [ ! -d $(MANPREFIX) ]; then mkdir -p $(MANPREFIX); fi
	cp treba $(BINPREFIX)
	cp ./man/treba.1 $(MANPREFIX)
