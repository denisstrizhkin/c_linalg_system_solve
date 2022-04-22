CFLAGS = -Wall -Wextra -std=c99 -pedantic
LDFLAGS = -static -s -lm

all: prog
	./prog

prog: prog.c
	$(CC) $(CFLAGS) prog.c $(LDFLAGS) -o $@

clean:
	rm -f prog
